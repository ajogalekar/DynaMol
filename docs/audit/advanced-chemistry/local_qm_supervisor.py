"""Local-storage quantum supervisor. Resource bounds are monitored, never model acceptance."""
from __future__ import annotations
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import threading
import time
from contextlib import contextmanager

class Cancelled(Exception):
    """A supervisor cancellation, recorded only after its group is cleaned up."""

@contextmanager
def _cancellation_scope():
    """Latch cancellation without raising asynchronously at any boundary.

    Ordinary code checks the flag while polling. Signals therefore cannot abort
    Popen assignment, a finally block, or entry into the cleanup context itself.
    The child does not inherit a blocked signal mask.
    """
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError('The native process supervisor must run on the main thread')
    signals = {signal.SIGINT, signal.SIGTERM}
    state = {'requested': None}
    def cancel(signum, frame):
        if state['requested'] is None:
            state['requested'] = signum
    old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    handlers = {sig: signal.signal(sig, cancel) for sig in signals}
    try:
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        yield state
    finally:
        signal.pthread_sigmask(signal.SIG_BLOCK, signals)
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)

@contextmanager
def _uninterrupted_cleanup():
    """Repeated cancellation cannot interrupt termination/reaping of our group."""
    signals = {signal.SIGINT, signal.SIGTERM}
    old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    handlers = {sig: signal.signal(sig, signal.SIG_IGN) for sig in signals}
    try:
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        yield
    finally:
        signal.pthread_sigmask(signal.SIG_BLOCK, signals)
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)

def _group_is_alive(group_id):
    # A leader may already be reaped while its descendants remain. macOS can
    # return EPERM for killpg(..., 0) during teardown; inspect live group members.
    rows = subprocess.check_output(['ps', '-axo', 'pgid=,stat='], text=True,
                                   timeout=3).splitlines()
    return any(len(parts := row.split()) == 2 and int(parts[0]) == group_id and
               not parts[1].startswith('Z') for row in rows)

def _terminate_owned_group(child, grace_seconds=1.0):
    """Terminate descendants even if their leader exited; always reap the leader.

    macOS does not expose Linux subreapers: grandchildren are reaped by their
    parent or launchd. We verify there are no live executable group members.
    """
    try:
        try:
            os.killpg(child.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        except PermissionError:
            if not _group_is_alive(child.pid):
                return
            raise
        deadline = time.monotonic()+grace_seconds
        while time.monotonic() < deadline:
            child.poll()
            if not _group_is_alive(child.pid):
                return
            time.sleep(.025)
        try:
            os.killpg(child.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except PermissionError:
            if _group_is_alive(child.pid):
                raise
        deadline = time.monotonic()+3
        while _group_is_alive(child.pid):
            if time.monotonic() >= deadline:
                raise RuntimeError('Owned process group still has live members after SIGKILL')
            child.poll()
            time.sleep(.025)
    finally:
        child.wait(timeout=3)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def verify_pins(pins):
    for path, expected in pins.items():
        if sha(path) != expected:
            raise ValueError('Pinned file changed or unavailable: '+path)


def process_rows():
    text = subprocess.check_output(['ps', '-axo', 'pid=,pgid=,rss=,stat=,command='],
                                   text=True, timeout=3)
    result = []
    for line in text.splitlines():
        parts = line.strip().split(None, 4)
        if len(parts) == 5:
            result.append(dict(pid=int(parts[0]), pgid=int(parts[1]),
                               rss_bytes=int(parts[2])*1024, stat=parts[3], command=parts[4]))
    return result


def is_qm_worker(row):
    try:
        words = shlex.split(row['command'])
    except ValueError:
        return False
    if not words or 'python' not in Path(words[0]).name.lower():
        return False
    # Classify the executed script, not a similarly named input or supervisor
    # argument. The qualified geometry adapter's child flag is --worker.
    index = 1
    while index < len(words) and words[index].startswith('-'):
        flag = words[index]
        if flag in ('-c', '-m', '-'):
            return False
        index += 2 if flag in ('-W', '-X') else 1
    if index >= len(words):
        return False
    name = Path(words[index]).name
    return name in ('qm_worker.py', 'covalent_reference_worker.py') or (
        name == 'geometric_adapter.py' and '--worker' in words[index+1:])


def outside_protected_storage(path):
    """Reject File Provider / user-source storage, including symlink targets."""
    resolved = Path(path).resolve()
    home = Path.home()
    protected = (home/'Documents', home/'Desktop', home/'Downloads',
                 home/'Library/Mobile Documents', home/'Library/CloudStorage')
    if any(resolved.is_relative_to(root.resolve()) for root in protected):
        raise ValueError('Recovery execution paths must be outside protected or cloud storage')
    return resolved


def local_storage_roots():
    local_root = outside_protected_storage(Path.home()/'.cache/dynamol-research')
    runtime_root = outside_protected_storage(Path.home()/'.cache/dynamol-runtimes')
    return local_root, runtime_root


def validate_spec(spec):
    expected = {'schema_version', 'label', 'root', 'output', 'runtime', 'worker', 'input',
                'pins', 'outer_wall_seconds', 'native_memory_mb', 'max_rss_bytes',
                'entry_disk_bytes', 'running_disk_bytes', 'threads'}
    if set(spec) != expected or spec['schema_version'] != 1:
        raise ValueError('Unknown or incomplete local job specification')
    for key in ('outer_wall_seconds', 'native_memory_mb', 'max_rss_bytes',
                'entry_disk_bytes', 'running_disk_bytes', 'threads'):
        if type(spec[key]) is not int or spec[key] <= 0:
            raise ValueError('Limits must be positive integers')
    if spec['threads'] > 2 or spec['native_memory_mb'] > 5000 or spec['max_rss_bytes'] > 8_000_000_000:
        raise ValueError('Resource admission exceeds the reviewed local recovery bounds')
    if spec['outer_wall_seconds'] > 21720 or spec['entry_disk_bytes'] < 35*1024**3 or spec['running_disk_bytes'] < 8*1024**3:
        raise ValueError('Wall or disk bounds exceed the reviewed local recovery scope')
    local_root, runtime_root = local_storage_roots()
    for key in ('root', 'output', 'runtime', 'worker', 'input'):
        if not Path(spec[key]).is_absolute():
            raise ValueError('Paths must be absolute')
    for key in ('output', 'worker', 'input'):
        if not outside_protected_storage(spec[key]).is_relative_to(local_root):
            raise ValueError('Worker, request and output must use local recovery storage')
    # Python itself is a stable uv-managed symlink outside Documents.
    runtime_path = Path(spec['runtime'])
    if not outside_protected_storage(runtime_path.parent).is_relative_to(runtime_root):
        raise ValueError('Use the verified local runtime copy')
    runtime_target = outside_protected_storage(runtime_path)
    uv_python_root = outside_protected_storage(Path.home()/'.local/share/uv/python')
    if not any(runtime_target.is_relative_to(root) for root in (runtime_root, uv_python_root)):
        raise ValueError('Runtime target must belong to local runtimes or pinned uv Python')
    if not isinstance(spec['pins'], dict) or not spec['pins']:
        raise ValueError('Pinned provenance is required')
    for key in ('worker', 'input', 'runtime'):
        if spec[key] not in spec['pins']:
            raise ValueError('Executable, worker and request must be pinned')
    for path, digest in spec['pins'].items():
        if not Path(path).is_absolute() or len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('Invalid provenance pin')
    request = json.loads(Path(spec['input']).read_text())
    if request.get('threads') != spec['threads'] or request.get('max_memory_mb') != spec['native_memory_mb']:
        raise ValueError('Request resource settings disagree with the controller')
    if type(request.get('max_wall_seconds')) is not int or request['max_wall_seconds'] > spec['outer_wall_seconds']-120:
        raise ValueError('Native deadline must leave the independent cleanup margin')
    if request.get('initial_checkpoint'):
        raise ValueError('Recovery admits geometry only; no unaccepted checkpoint reuse')
    return request


def resource_sample(child_pid, max_rss_bytes, running_disk_bytes, output):
    rows = process_rows()
    group = [row for row in rows if row['pgid'] == child_pid and not row['stat'].startswith('Z')]
    maximum = max((row['rss_bytes'] for row in group), default=0)
    free = shutil.disk_usage(output).free
    sample = dict(sample_unix=time.time(), maximum_process_rss_bytes=maximum,
                  process_group_rss_bytes=sum(row['rss_bytes'] for row in group), free_disk_bytes=free)
    if maximum > max_rss_bytes:
        raise MemoryError('Owned quantum process exceeded sampled RSS limit: '+json.dumps(sample))
    if free < running_disk_bytes:
        raise RuntimeError('Free disk fell below running reserve: '+json.dumps(sample))
    return sample


def supervise(command, output, timeout, environment, on_start, on_sample,
              max_rss_bytes=8_000_000_000, running_disk_bytes=8*1024**3):
    child = None
    deadline = time.monotonic()+timeout
    next_sample = 0.0
    with _cancellation_scope() as cancellation:
        with (Path(output)/'controller-native.log').open('x') as log:
            try:
                child = subprocess.Popen(command, cwd=output, stdout=log, stderr=subprocess.STDOUT,
                                         start_new_session=True, env=environment)
                on_start(child.pid)
                while True:
                    if cancellation['requested'] is not None:
                        raise Cancelled('Cancellation signal '+str(cancellation['requested']))
                    now = time.monotonic()
                    if now >= deadline:
                        raise TimeoutError('Independent wall deadline exceeded')
                    if now >= next_sample:
                        on_sample(resource_sample(child.pid, max_rss_bytes, running_disk_bytes, output))
                        next_sample = now+1.0
                    try:
                        code = child.wait(timeout=min(.1, deadline-time.monotonic()))
                        break
                    except subprocess.TimeoutExpired:
                        pass
            finally:
                if child is not None:
                    with _uninterrupted_cleanup():
                        _terminate_owned_group(child)
            if cancellation['requested'] is not None:
                raise Cancelled('Cancellation during cleanup')
            return code


def run(spec_path):
    local_root, _ = local_storage_roots()
    spec_path = outside_protected_storage(spec_path)
    controller_path = outside_protected_storage(__file__)
    if not spec_path.is_relative_to(local_root) or not controller_path.is_relative_to(local_root):
        raise ValueError('Stage the controller and specification in local recovery storage before execution')
    spec_bytes = spec_path.read_bytes()
    spec = json.loads(spec_bytes)
    request = validate_spec(spec)
    output = Path(spec['output'])
    # A run is claimed by creating a new directory; failure never permits reuse.
    output.mkdir(parents=True, exist_ok=False)
    state = dict(schema_version=1, status='preflight', controller_pid=os.getpid(), child_pid=None,
                 started_unix=time.time(), spec_sha256=hashlib.sha256(spec_bytes).hexdigest(),
                 controller_sha256=sha(__file__), physical_acceptance=False,
                 simulation_ready=False, label=spec['label'], limits={
                     'native_memory_mb': spec['native_memory_mb'], 'max_rss_bytes': spec['max_rss_bytes'],
                     'rss_limit_kind': 'Sampled once per second; terminates on observed crossing, not an OS allocation cap',
                     'threads': spec['threads'], 'outer_wall_seconds': spec['outer_wall_seconds']})
    def update(**values):
        state.update(values)
        state['updated_unix'] = time.time()
        write(output/'controller-progress.json', state)
    try:
        update()
        verify_pins(spec['pins'])
        write(output/'spec.json', spec)
        (output/'controller-source.py').write_bytes(Path(__file__).read_bytes())
        captured = {'input': output/'input.json', 'worker': output/'qm_worker.py'}
        for key, destination in captured.items():
            raw = Path(spec[key]).read_bytes()
            if hashlib.sha256(raw).hexdigest() != spec['pins'][spec[key]]:
                raise ValueError('Pinned '+key+' changed while capturing execution bytes')
            destination.write_bytes(raw)
            if sha(destination) != spec['pins'][spec[key]]:
                raise ValueError('Captured '+key+' does not match its original pin')
        scratch = output/'scratch'
        scratch.mkdir()
        env = dict(os.environ)
        env.pop('PYTHONPATH', None)
        env.pop('PYTHONHOME', None)
        env.update(TMPDIR=str(scratch), PYSCF_TMPDIR=str(scratch), PYTHONDONTWRITEBYTECODE='1')
        for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                    'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            env[key] = str(spec['threads'])
        lock_path = Path(spec['root'])/'build/advanced-chemistry/QM-LAUNCH.lock'
        with lock_path.open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            active = [row for row in process_rows() if is_qm_worker(row)]
            free = shutil.disk_usage(output).free
            # Recovery is serial to avoid competing DF scratch and memory peaks.
            if active or free < spec['entry_disk_bytes']:
                raise RuntimeError('Recovery admission held: active quantum work or insufficient disk')
            verify_pins(spec['pins'])
            for key, destination in captured.items():
                if sha(destination) != spec['pins'][spec[key]]:
                    raise ValueError('Captured '+key+' changed before launch')
            if sha(spec_path) != state['spec_sha256'] or sha(__file__) != state['controller_sha256']:
                raise ValueError('Controller or specification changed before launch')
            def started(pid):
                update(status='running', child_pid=pid, qm_started_unix=time.time(),
                       outer_deadline_unix=time.time()+spec['outer_wall_seconds'])
                # Shared lock protects the spawn, not the entire compute lifetime.
                fcntl.flock(lock, fcntl.LOCK_UN)
            def sampled(sample):
                state['peak_sampled_rss_bytes'] = max(state.get('peak_sampled_rss_bytes', 0), sample['maximum_process_rss_bytes'])
                if time.time()-state.get('last_sample_saved_unix', 0) >= 15:
                    update(resource_sample=sample, last_sample_saved_unix=time.time())
            code = supervise([spec['runtime'], str(captured['worker']), str(captured['input']), str(output/'optimization')],
                             output, spec['outer_wall_seconds'], env, started, sampled,
                             spec['max_rss_bytes'], spec['running_disk_bytes'])
        verify_pins(spec['pins'])
        for key, destination in captured.items():
            if sha(destination) != spec['pins'][spec[key]]:
                raise ValueError('Captured '+key+' changed during execution')
        if sha(spec_path) != state['spec_sha256'] or sha(__file__) != state['controller_sha256']:
            raise ValueError('Controller or specification changed during execution')
        result_path = output/'optimization/result.json'
        result = json.loads(result_path.read_text()) if result_path.exists() else {}
        if code or not result.get('accepted') or not result.get('optimization', {}).get('converged'):
            raise ValueError('Native optimization did not complete successfully')
        if result.get('input_sha256') != spec['pins'][spec['input']] or result.get('worker_sha256') != spec['pins'][spec['worker']]:
            raise ValueError('Native result provenance mismatch')
        if sha(output/'optimization/arrays.npz') != result.get('arrays_sha256'):
            raise ValueError('Native result array mismatch')
        update(status='numerically_complete_pending_geometry_review', result_sha256=sha(result_path), child_pid=None)
    except Exception as error:
        update(status='failed', error_type=type(error).__name__, error=str(error), child_pid=None)
    write(output/'controller-result.json', state)
    return state


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('spec')
    args = parser.parse_args()
    result = run(args.spec)
    print(json.dumps(result, indent=2), flush=True)
    raise SystemExit(1 if result['status'] == 'failed' else 0)
