"""Independent local-storage/provenance checks; all workers are synthetic stubs."""
from pathlib import Path
import hashlib
import importlib.util
import json
from types import SimpleNamespace

import pytest

SOURCE = Path(__file__).with_name('local_qm_supervisor.py')
loader = importlib.util.spec_from_file_location('local_qm_admission_review', SOURCE)
module = importlib.util.module_from_spec(loader)
loader.loader.exec_module(module)


@pytest.fixture
def setup(tmp_path, monkeypatch):
    home = tmp_path/'home'
    home.mkdir()
    monkeypatch.setattr(Path, 'home', classmethod(lambda cls: home))
    local = home/'.cache/dynamol-research/example'
    local.mkdir(parents=True)
    runtime = home/'.cache/dynamol-runtimes/example/bin/python'
    runtime.parent.mkdir(parents=True)
    runtime.write_text('synthetic runtime: never executed')
    worker = local/'source-worker.py'
    worker.write_text('synthetic worker: never executed')
    request = local/'request.json'
    request.write_text(json.dumps(dict(threads=2, max_memory_mb=5000, max_wall_seconds=21600)))
    root = tmp_path/'project'
    (root/'build/advanced-chemistry').mkdir(parents=True)
    controller = local/'controller.py'
    controller.write_bytes(SOURCE.read_bytes())
    monkeypatch.setattr(module, '__file__', str(controller))
    spec = dict(schema_version=1, label='synthetic admission fixture', root=str(root),
                output=str(local/'output'), runtime=str(runtime), worker=str(worker),
                input=str(request), pins={}, outer_wall_seconds=21720,
                native_memory_mb=5000, max_rss_bytes=8_000_000_000,
                entry_disk_bytes=35*1024**3, running_disk_bytes=8*1024**3, threads=2)
    spec['pins'] = {spec[key]: module.sha(spec[key]) for key in ('runtime', 'worker', 'input')}
    spec_path = local/'spec.json'
    spec_path.write_text(json.dumps(spec))
    monkeypatch.setattr(module, 'process_rows', lambda: [])
    monkeypatch.setattr(module.shutil, 'disk_usage', lambda path: SimpleNamespace(free=100*1024**3))
    return SimpleNamespace(home=home, local=local, runtime=runtime, worker=worker,
                           request=request, controller=controller, spec=spec, spec_path=spec_path)


def test_pinned_uv_interpreter_symlink_allowed(setup):
    target = setup.home/'.local/share/uv/python/version/bin/python'
    target.parent.mkdir(parents=True)
    target.write_text('synthetic uv interpreter')
    setup.runtime.unlink()
    setup.runtime.symlink_to(target)
    setup.spec['pins'][str(setup.runtime)] = module.sha(target)
    assert module.validate_spec(setup.spec)['max_memory_mb'] == 5000


@pytest.mark.parametrize('protected', ['Documents', 'Desktop', 'Downloads',
                                      'Library/Mobile Documents', 'Library/CloudStorage'])
def test_runtime_symlink_to_protected_storage_rejected(setup, protected):
    target = setup.home/protected/'python'
    target.parent.mkdir(parents=True)
    target.write_text('synthetic protected interpreter')
    setup.runtime.unlink()
    setup.runtime.symlink_to(target)
    with pytest.raises(ValueError, match='protected or cloud'):
        module.validate_spec(setup.spec)


@pytest.mark.parametrize('key', ['worker', 'input', 'output'])
def test_execution_path_symlink_to_documents_rejected(setup, key):
    target = setup.home/'Documents'/key
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text('synthetic protected source')
    path = Path(setup.spec[key])
    if path.exists():
        path.unlink()
    path.symlink_to(target)
    with pytest.raises(ValueError, match='protected or cloud'):
        module.validate_spec(setup.spec)


@pytest.mark.parametrize('which', ['spec', 'controller'])
def test_run_entry_rejects_unstaged_paths_before_claim(setup, monkeypatch, which):
    outside = setup.home/'outside.py'
    outside.write_text('{}')
    if which == 'controller':
        monkeypatch.setattr(module, '__file__', str(outside))
    with pytest.raises(ValueError, match='Stage the controller'):
        module.run(outside if which == 'spec' else setup.spec_path)
    assert not Path(setup.spec['output']).exists()


@pytest.mark.parametrize('command,expected', [
    ('/tmp/python -u -X dev /tmp/geometric_adapter.py --worker --root /tmp/root', True),
    ('/tmp/python -W ignore /tmp/qm_worker.py /tmp/input /tmp/output', True),
    ('/tmp/Python /tmp/covalent_reference_worker.py /tmp/input /tmp/output', True),
    ('/tmp/python /tmp/not-a-worker.py /tmp/qm_worker.py', False),
    ('/tmp/python -c "qm_worker.py"', False),
    ('/tmp/python -m package /tmp/qm_worker.py', False),
    ('/tmp/python /tmp/geometric_adapter.py --native-child', False),
])
def test_actual_executed_script_classification(command, expected):
    assert module.is_qm_worker({'command': command}) is expected


def install_stub(setup, monkeypatch, mutation=None):
    commands = []
    def stub(command, output, timeout, environment, on_start, on_sample, *limits):
        commands.append(command)
        output = Path(output)
        assert command[1:3] == [str(output/'qm_worker.py'), str(output/'input.json')]
        assert Path(command[1]).read_bytes() == setup.worker.read_bytes()
        assert Path(command[2]).read_bytes() == setup.request.read_bytes()
        folder = output/'optimization'
        folder.mkdir()
        arrays = folder/'arrays.npz'
        arrays.write_bytes(b'synthetic contract fixture, not numerical arrays')
        (folder/'result.json').write_text(json.dumps(dict(
            accepted=True, optimization={'converged': True},
            input_sha256=module.sha(command[2]), worker_sha256=module.sha(command[1]),
            arrays_sha256=module.sha(arrays))))
        if mutation is not None:
            mutation(output)
        return 0
    monkeypatch.setattr(module, 'supervise', stub)
    return commands


def test_launch_uses_verified_captured_worker_and_input(setup, monkeypatch):
    commands = install_stub(setup, monkeypatch)
    result = module.run(setup.spec_path)
    assert result['status'] == 'numerically_complete_pending_geometry_review'
    assert result['physical_acceptance'] is False and result['simulation_ready'] is False
    assert len(commands) == 1


@pytest.mark.parametrize('key', ['worker', 'input'])
def test_changed_source_is_rejected_before_spawn(setup, monkeypatch, key):
    commands = install_stub(setup, monkeypatch)
    def rows():
        Path(setup.spec[key]).write_text('changed after capture')
        return []
    monkeypatch.setattr(module, 'process_rows', rows)
    result = module.run(setup.spec_path)
    assert result['status'] == 'failed' and 'Pinned file changed' in result['error']
    assert not commands


@pytest.mark.parametrize('target', ['runtime', 'captured-input', 'captured-worker', 'spec'])
def test_changed_completion_binding_rejected(setup, monkeypatch, target):
    def mutation(output):
        path = {'runtime': setup.runtime, 'captured-input': output/'input.json',
                'captured-worker': output/'qm_worker.py', 'spec': setup.spec_path}[target]
        path.write_text('changed during synthetic execution')
    install_stub(setup, monkeypatch, mutation)
    result = module.run(setup.spec_path)
    assert result['status'] == 'failed'
    assert 'changed' in result['error']

