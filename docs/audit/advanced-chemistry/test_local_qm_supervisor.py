"""No-QM contract and real disposable-process tests for local recovery supervision."""
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest

SOURCE = Path(__file__).with_name('local_qm_supervisor.py')
spec = importlib.util.spec_from_file_location('local_qm_supervisor', SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def run_process(tmp_path, monkeypatch, script='import time; time.sleep(20)', **options):
    holder = []
    samples = []
    on_start = options.pop('on_start', holder.append)
    def started(pid):
        holder.append(pid)
        on_start(pid)
    try:
        return module.supervise([sys.executable, '-c', script], tmp_path,
                                options.pop('timeout', .2), dict(os.environ), started,
                                options.pop('on_sample', samples.append), **options)
    finally:
        for pid in set(holder):
            assert not module._group_is_alive(pid)
            with pytest.raises(ChildProcessError):
                os.waitpid(pid, os.WNOHANG)


def test_copied_lifecycle_functions_exact():
    proof = json.loads(SOURCE.with_name('local-qm-supervisor-lifecycle-provenance.json').read_text())
    functions = {node.name: node for node in ast.parse(SOURCE.read_text()).body if hasattr(node, 'name')}
    for name, expected in proof['ast_sha256'].items():
        assert hashlib.sha256(ast.dump(functions[name], include_attributes=False).encode()).hexdigest() == expected


def test_normal_exit(tmp_path, monkeypatch):
    assert run_process(tmp_path, monkeypatch, 'print("done")', timeout=5) == 0


def test_timeout_reaps_descendants(tmp_path, monkeypatch):
    script = 'import subprocess,sys,time; subprocess.Popen([sys.executable,"-c","import time; time.sleep(20)"]); time.sleep(20)'
    with pytest.raises(TimeoutError):
        run_process(tmp_path, monkeypatch, script)


def test_start_artifact_failure_cleans_up(tmp_path, monkeypatch):
    def error(pid):
        raise OSError('injected artifact failure')
    with pytest.raises(OSError, match='artifact'):
        run_process(tmp_path, monkeypatch, on_start=error)


def test_sample_artifact_failure_cleans_up(tmp_path, monkeypatch):
    def error(sample):
        raise OSError('injected sample write failure')
    with pytest.raises(OSError, match='sample'):
        run_process(tmp_path, monkeypatch, on_sample=error)


def test_sampled_memory_crossing_cleans_up(tmp_path, monkeypatch):
    with pytest.raises(MemoryError, match='RSS'):
        run_process(tmp_path, monkeypatch, timeout=5, max_rss_bytes=1)


def test_low_disk_cleans_up(tmp_path, monkeypatch):
    with pytest.raises(RuntimeError, match='disk'):
        run_process(tmp_path, monkeypatch, timeout=5, running_disk_bytes=10**20)


def test_cancellation_reaps(tmp_path, monkeypatch):
    def cancel(pid):
        os.kill(os.getpid(), signal.SIGTERM)
        os.kill(os.getpid(), signal.SIGINT)
    with pytest.raises(module.Cancelled):
        run_process(tmp_path, monkeypatch, on_start=cancel)


def test_leader_exit_descendant_cleanup(tmp_path, monkeypatch):
    script = 'import subprocess,sys; subprocess.Popen([sys.executable,"-c","import time; time.sleep(20)"])'
    assert run_process(tmp_path, monkeypatch, script, timeout=5) == 0


@pytest.mark.parametrize('command,expected', [
    ('/tmp/python /tmp/qm_worker.py /tmp/in.json /tmp/out', True),
    ('/tmp/python /tmp/covalent_reference_worker.py /tmp/in.json /tmp/out', True),
    ('/tmp/python /tmp/geometric_adapter.py --worker /tmp/request', True),
    ('/tmp/python /tmp/geometric_adapter.py --native-child /tmp/request', False),
    ('/tmp/python /tmp/test_geometric_adapter.py', False),
    ('/bin/bash -c "python qm_worker.py"', False),
    ('/tmp/python /tmp/continue_covalent_qm.py', False),
    ('/tmp/python -c "print(\'qm_worker.py\')"', False),
])
def test_worker_classification(command, expected):
    assert module.is_qm_worker({'command': command}) is expected


def valid_spec(tmp_path):
    # Validation resolves paths under the explicit local cache; fixture input is
    # provided through read_text monkeypatch for a no-mutation admission check.
    local = Path.home()/'.cache/dynamol-research/test-proposed-only'
    runtime = Path.home()/'.cache/dynamol-runtimes/test-proposed-only/bin/python'
    result = dict(schema_version=1, label='test', root='/tmp/root', output=str(local/'out'),
                  runtime=str(runtime), worker=str(local/'qm_worker.py'), input=str(local/'input.json'),
                  pins={}, outer_wall_seconds=21720, native_memory_mb=5000,
                  max_rss_bytes=8_000_000_000, entry_disk_bytes=35*1024**3,
                  running_disk_bytes=8*1024**3, threads=2)
    result['pins'] = {result[key]: 'a'*64 for key in ('runtime','worker','input')}
    return result


@pytest.mark.parametrize('key,value', [
    ('threads', 3), ('threads', True), ('max_rss_bytes', 8_000_000_001),
    ('native_memory_mb', 8000), ('outer_wall_seconds', 21721),
    ('entry_disk_bytes', 34*1024**3), ('running_disk_bytes', 7*1024**3),
    ('output', '/Users/ashujo/Documents/recovery'), ('worker', '/tmp/qm_worker.py'),
    ('input', 'relative.json'), ('runtime', '/Users/ashujo/Documents/python'),
])
def test_reject_unsafe_spec(tmp_path, monkeypatch, key, value):
    data = valid_spec(tmp_path)
    data[key] = value
    with pytest.raises(ValueError):
        module.validate_spec(data)


def test_spec_request_agreement_and_no_checkpoint(tmp_path, monkeypatch):
    data = valid_spec(tmp_path)
    request = dict(threads=2, max_memory_mb=5000, max_wall_seconds=21600)
    monkeypatch.setattr(Path, 'read_text', lambda self: json.dumps(request))
    assert module.validate_spec(data) == request
    for change in ({'threads': 1}, {'max_memory_mb': 8000}, {'max_wall_seconds': 21720},
                   {'initial_checkpoint': {'result_path': '/tmp/failed'}}):
        altered = dict(request, **change)
        monkeypatch.setattr(Path, 'read_text', lambda self: json.dumps(altered))
        with pytest.raises(ValueError):
            module.validate_spec(data)
