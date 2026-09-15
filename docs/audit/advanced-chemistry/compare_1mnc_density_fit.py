"""Bounded, resource-queued same-geometry 1MNC density-fitting comparison.

This audit computes an approximation error, not a physical force-field verdict.
It retains the conventional result and never starts geometry optimization.
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / 'build/advanced-chemistry/metal/qm-1mnc-small-screen-v1/calculation'
OUT = ROOT / 'build/advanced-chemistry/metal/qm-1mnc-df-comparison-v2'
PRIOR = ROOT / 'build/advanced-chemistry/metal/qm-1mnc-df-comparison-v1'
RUNTIME = ROOT / '.tools/qm-parallel/bin/python'
MIN_FREE_BYTES = 22 * 1024**3
CONCURRENT_DF_FREE_BYTES = 35 * 1024**3
RESERVE_BYTES = 8 * 1024**3
WAIT_SECONDS = 21600
NATIVE_SECONDS = 7200


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def active_workers():
    found = []
    for line in subprocess.check_output(['ps', '-axo', 'pid=,command='], text=True).splitlines():
        row = line.strip().split(None, 1)
        if len(row) != 2:
            continue
        if 'qm_worker.py' not in row[1]:
            continue
        try:
            args = shlex.split(row[1])
        except ValueError:
            raise ValueError('Cannot safely classify a running QM worker command')
        if (not args or Path(args[0]).name not in ('python', 'python3', 'python3.12')
                or len(args) != 4 or Path(args[1]).name != 'qm_worker.py'):
            continue
        # Unknown/unreadable inputs conservatively occupy the large-DF slot.
        density_fit = True
        request = Path(args[2])
        if request.is_file() and request.stat().st_size < 1024**2:
            try:
                density_fit = bool(json.loads(request.read_text()).get('density_fit', False))
            except (OSError, ValueError):
                pass
        found.append({'pid': int(row[0]), 'density_fit': density_fit, 'input': str(request)})
    return found


def terminate(child):
    if child.poll() is not None:
        return
    os.killpg(child.pid, signal.SIGTERM)
    try:
        child.wait(timeout=15)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGKILL)
        child.wait()


def run():
    OUT.mkdir(parents=True, exist_ok=False)
    state = {'status': 'preflight', 'controller_pid': os.getpid(), 'started_unix': time.time(),
        'scope': 'Same-geometry numerical approximation measurement only; no model fitting, optimization or MD.',
        'full_preparation_ready': False, 'minimum_start_free_bytes': MIN_FREE_BYTES,
        'concurrent_df_minimum_start_free_bytes': CONCURRENT_DF_FREE_BYTES,
        'maximum_large_df_jobs': 2,
        'minimum_running_reserve_bytes': RESERVE_BYTES, 'wait_seconds': WAIT_SECONDS,
        'native_wall_seconds': NATIVE_SECONDS, 'outer_wall_seconds': NATIVE_SECONDS + 120,
        'script_sha256': sha(__file__), 'source': str(SOURCE)}
    child = None
    try:
        import numpy as np
        original = json.loads((SOURCE / 'input.json').read_text())
        reference = json.loads((SOURCE / 'result.json').read_text())
        if (reference.get('accepted') is not True or reference.get('optimization', {}).get('requested')
                or reference.get('method', {}).get('density_fit') is not False
                or reference.get('input_sha256') != sha(SOURCE / 'input.json')
                or reference.get('arrays_sha256') != sha(SOURCE / 'arrays.npz')):
            raise ValueError('Conventional reference failed its pinned input/result checks')
        with np.load(SOURCE / 'arrays.npz', allow_pickle=False) as saved:
            xyz = np.array(saved['coords_angstrom'])
            gradient = np.array(saved['gradient_hartree_per_bohr'])
        if (not np.isfinite(gradient).all() or gradient.shape != (97, 3)
                or not np.allclose(xyz, np.asarray(original['coords_angstrom']), atol=1e-12, rtol=0)
                or original['atom_ids'] != reference['atom_ids']):
            raise ValueError('Reference atom order, geometry or gradient is invalid')
        request = dict(original, density_fit=True, auxbasis='def2-universal-jkfit',
            max_wall_seconds=NATIVE_SECONDS,
            initial_checkpoint={'result_path': str(SOURCE / 'result.json'),
                'result_sha256': sha(SOURCE / 'result.json'),
                'checkpoint_path': str(SOURCE / 'scf.chk'),
                'checkpoint_sha256': sha(SOURCE / 'scf.chk')})
        write(OUT / 'input.json', request)
        worker = OUT / 'qm_worker.py'
        worker.write_bytes((ROOT / 'backend/qm_worker.py').read_bytes())
        prior=json.loads((PRIOR / 'intentional-stop.json').read_text())
        if prior.get('status')!='intentionally_stopped_before_qm_launch':
            raise ValueError('Prior controller lacks an explicit pre-launch stop record')
        state.update(status='waiting_for_resources', input_sha256=sha(OUT / 'input.json'),
            worker_sha256=sha(worker), source_result_sha256=sha(SOURCE / 'result.json'),
            source_arrays_sha256=sha(SOURCE / 'arrays.npz'),
            reference_energy_hartree=reference['energy_hartree'],
            reference_max_gradient_component_hartree_per_bohr=float(np.abs(gradient).max()),
            reference_rms_gradient_component_hartree_per_bohr=float(np.sqrt(np.mean(gradient**2))),
            prior_controller_stop_sha256=sha(PRIOR / 'intentional-stop.json'),
            wait_deadline_unix=min(time.time() + WAIT_SECONDS, prior['prior_progress']['wait_deadline_unix']))
        if state['input_sha256']!=prior['prior_progress']['input_sha256']:
            raise ValueError('Resource-only restart changed the scientific request')
        write(OUT / 'progress.json', state)
        with (ROOT / 'build/advanced-chemistry/QM-LAUNCH.lock').open('a') as launch_lock:
            while child is None:
                if time.time() > state['wait_deadline_unix']:
                    raise TimeoutError('Resource wait ended; no QM calculation was started')
                fcntl.flock(launch_lock, fcntl.LOCK_EX)
                try:
                    workers = active_workers()
                    free = shutil.disk_usage(OUT).free
                    df_count=sum(p['density_fit'] for p in workers)
                    threshold=CONCURRENT_DF_FREE_BYTES if df_count else MIN_FREE_BYTES
                    state.update(free_disk_bytes=free, other_workers=workers, current_entry_reserve_bytes=threshold)
                    if len(workers) < 3 and df_count < 2 and free >= threshold:
                        log = (OUT / 'controller.log').open('w')
                        child = subprocess.Popen([str(RUNTIME), str(worker), str(OUT / 'input.json'),
                            str(OUT / 'calculation')], stdout=log, stderr=subprocess.STDOUT,
                            start_new_session=True)
                        state.update(status='running', worker_pid=child.pid, qm_started_unix=time.time())
                    write(OUT / 'progress.json', state)
                finally:
                    fcntl.flock(launch_lock, fcntl.LOCK_UN)
                if child is None:
                    time.sleep(30)
        deadline = time.time() + NATIVE_SECONDS + 120
        while child.poll() is None:
            if time.time() > deadline:
                raise TimeoutError('Independent QM wall-time limit reached')
            free = shutil.disk_usage(OUT).free
            state.update(free_disk_bytes=free)
            write(OUT / 'progress.json', state)
            if free < RESERVE_BYTES:
                raise RuntimeError('Stopped to preserve the documented minimum disk reserve')
            time.sleep(15)
        result = json.loads((OUT / 'calculation/result.json').read_text())
        if child.returncode or not result.get('accepted'):
            raise ValueError('DF numerical calculation failed: ' + str(result.get('error')))
        if result['input_sha256'] != state['input_sha256'] or result['worker_sha256'] != state['worker_sha256']:
            raise ValueError('DF calculation changed its pinned input or worker')
        if sha(OUT / 'calculation/arrays.npz') != result['arrays_sha256']:
            raise ValueError('DF output array hash differs')
        with np.load(OUT / 'calculation/arrays.npz', allow_pickle=False) as calculated:
            if (result['atom_ids'] != reference['atom_ids']
                    or not np.array_equal(calculated['coords_angstrom'], xyz)):
                raise ValueError('Comparison does not use identical ordered geometry')
            difference = calculated['gradient_hartree_per_bohr'] - gradient
        state.update(status='completed', numerical_calculation_converged=True,
            df_minus_conventional_energy_hartree=result['energy_hartree'] - reference['energy_hartree'],
            maximum_gradient_component_difference_hartree_per_bohr=float(np.abs(difference).max()),
            rms_gradient_component_difference_hartree_per_bohr=float(np.sqrt(np.mean(difference**2))),
            result_sha256=sha(OUT / 'calculation/result.json'),
            interpretation='One fixed geometry only. No relative conformer-energy or physical force-field validation.')
    except Exception as error:
        if child is not None:
            terminate(child)
        state.update(status='failed', error={'type': type(error).__name__, 'message': str(error)})
    state['elapsed_seconds'] = time.time() - state['started_unix']
    write(OUT / 'progress.json', state)
    write(OUT / 'result.json', state)
    print(json.dumps(state, indent=2))


if __name__ == '__main__':
    run()
