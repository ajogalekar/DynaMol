"""Integrity/eligibility guards; native continuation is checked by check_recovery.py."""
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from backend import config, jobs, recovery
from backend.storage import atomic_json


@pytest.fixture
def checkpoint(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'JOBS_DIR', tmp_path)
    monkeypatch.setattr(recovery, 'runtime_identity', lambda engine: {'engine': engine, 'cpu_threads': 2})
    folder = tmp_path / 'audit-job'
    folder.mkdir()
    for name in ['config.json', 'input.pdb', 'input-state.json', 'prepared.pdb', 'system.xml', 'integrator.xml']:
        (folder / name).write_text(name)
    for name in ['ligands/LIG/parameters.xml', 'residue-parameters/phosaa.xml']:
        target = folder / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(name)
    manifest = recovery.create_manifest(folder, 'openmm')
    cp = folder / 'checkpoints/step-000000250-audit/checkpoint.chk'
    cp.parent.mkdir(parents=True)
    cp.write_bytes(b'Integrity fixture, not a native checkpoint')
    frame = folder / 'frames/frame-000000000.npz'
    frame.parent.mkdir()
    frame.write_bytes(b'Integrity fixture, not a trajectory')
    manifest['checkpoint'] = {'directory': cp.parent.name, 'sha256': recovery.digest(cp), 'step': 250, 'frames': [{'name': frame.name, 'sha256': recovery.digest(frame), 'step': 0}]}
    atomic_json(folder / 'recovery.json', manifest)
    return folder, manifest


def test_checkpoint_integrity_covers_prepared_ligand_and_modified_bundles(checkpoint):
    folder, manifest = checkpoint
    assert 'ligands/LIG/parameters.xml' in manifest['dependencies']
    assert 'residue-parameters/phosaa.xml' in manifest['dependencies']
    assert recovery.validate_manifest(folder, 'openmm') == manifest
    (folder / 'ligands/LIG/parameters.xml').write_text('Changed charge')
    with pytest.raises(ValueError, match='parameters changed'):
        recovery.validate_manifest(folder, 'openmm')


@pytest.mark.parametrize('target', ['system.xml', 'integrator.xml', 'input-state.json', 'residue-parameters/phosaa.xml', 'config.json'])
def test_changed_prepared_state_is_never_recreated_silently(checkpoint, target):
    folder, _ = checkpoint
    (folder / target).write_text('altered')
    with pytest.raises(ValueError, match='changed'):
        recovery.validate_manifest(folder, 'openmm')


def test_checkpoint_and_saved_prefix_corruption_rejected(checkpoint):
    folder, manifest = checkpoint
    cp = folder / 'checkpoints' / manifest['checkpoint']['directory'] / 'checkpoint.chk'
    original = cp.read_bytes()
    cp.write_bytes(b'broken')
    with pytest.raises(ValueError, match='missing or damaged'):
        recovery.validate_manifest(folder, 'openmm')
    cp.write_bytes(original)
    (folder / 'frames/frame-000000000.npz').unlink()
    with pytest.raises(ValueError, match='continuous trajectory'):
        recovery.validate_manifest(folder, 'openmm')


def test_runtime_change_blocks_binary_checkpoint(checkpoint, monkeypatch):
    folder, _ = checkpoint
    monkeypatch.setattr(recovery, 'runtime_identity', lambda engine: {'engine': engine, 'cpu_threads': 4})
    with pytest.raises(ValueError, match='original runtime'):
        recovery.validate_manifest(folder, 'openmm')


def test_missing_checkpoint_and_non_simulation_jobs_explain_why_unavailable(checkpoint):
    folder, _ = checkpoint
    (folder / 'recovery.json').unlink()
    info = recovery.recovery_info({'id': folder.name, 'engine': 'openmm', 'status': 'cancelled'})
    assert not info['available'] and 'Preparation and initial relaxation' in info['reason']
    assert not recovery.recovery_info({'engine': 'preparation'})['available']


def test_resume_is_disabled_while_engine_active(checkpoint):
    folder, _ = checkpoint
    info = recovery.recovery_info({'id': folder.name, 'engine': 'openmm', 'status': 'running'}, verify=True)
    assert not info['available'] and info['checkpoint_saved']


def test_orphaned_native_engine_is_stopped_before_resuming(checkpoint, monkeypatch):
    folder, _ = checkpoint
    job = {'id': folder.name, 'engine': 'gromacs', 'status': 'running', 'created_at': '2026-09-11T00:00:00+00:00', 'worker_pid': 12, 'native_pid': 34}
    atomic_json(folder / 'status.json', job)
    signals = []
    def kill(pid, sig):
        if pid == 12:
            raise ProcessLookupError()
        signals.append((pid, sig))
    monkeypatch.setattr(jobs.os, 'kill', kill)
    monkeypatch.setattr(jobs.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=f'gmx mdrun -deffnm {folder}/production'))
    observed = jobs.get_job(folder.name)
    assert observed['status'] == 'cancelling' and observed['orphan_stop_requested']
    assert signals == [(34, jobs.signal.SIGTERM)]
    jobs.get_job(folder.name)
    assert len(signals) == 1
    monkeypatch.setattr(jobs.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=''))
    assert jobs.get_job(folder.name)['status'] == 'interrupted'


def test_repeated_stop_does_not_interrupt_native_checkpoint_flush(monkeypatch):
    job = {'id': 'audit', 'status': 'cancelling'}
    monkeypatch.setattr(jobs, 'get_job', lambda job_id: job)
    monkeypatch.setattr(jobs.os, 'kill', lambda *args: pytest.fail('Repeated stop must not send another signal'))
    assert jobs.cancel_job('audit') == job


@pytest.mark.parametrize('engine', ['openmm', 'gromacs', 'preparation', 'solvation'])
def test_stop_during_worker_startup_finishes_cancelled(checkpoint, monkeypatch, engine):
    folder, _ = checkpoint
    job = {'id': folder.name, 'engine': engine, 'status': 'cancelling', 'created_at': '2026-09-11T00:00:00+00:00', 'worker_pid': 12, 'elapsed_seconds': 0}
    atomic_json(folder / 'status.json', job)
    (folder / 'cancel.request').touch()
    def stopped_before_handler(pid, sig):
        raise ProcessLookupError()
    monkeypatch.setattr(jobs.os, 'kill', stopped_before_handler)
    observed = jobs.get_job(folder.name)
    assert observed['status'] == 'cancelled'
    assert observed['stage'] == 'Cancelled'
    assert 'by user' in observed['error']
    assert observed['elapsed_seconds'] > 0
    assert json.loads((folder / 'status.json').read_text())['status'] == 'cancelled'


def test_unrequested_worker_loss_stays_interrupted(checkpoint, monkeypatch):
    folder, _ = checkpoint
    atomic_json(folder / 'status.json', {'id': folder.name, 'engine': 'openmm', 'status': 'running', 'created_at': '2026-09-11T00:00:00+00:00', 'worker_pid': 12})
    def unexpectedly_stopped(pid, sig):
        raise ProcessLookupError()
    monkeypatch.setattr(jobs.os, 'kill', unexpectedly_stopped)
    assert jobs.get_job(folder.name)['status'] == 'interrupted'


def test_requested_stop_waits_for_orphan_before_finishing_cancelled(checkpoint, monkeypatch):
    folder, _ = checkpoint
    atomic_json(folder / 'status.json', {'id': folder.name, 'engine': 'gromacs', 'status': 'cancelling', 'created_at': '2026-09-11T00:00:00+00:00', 'worker_pid': 12, 'native_pid': 34})
    (folder / 'cancel.request').touch()
    signals = []
    def kill(pid, sig):
        if pid == 12:
            raise ProcessLookupError()
        signals.append((pid, sig))
    monkeypatch.setattr(jobs.os, 'kill', kill)
    monkeypatch.setattr(jobs.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=f'gmx mdrun -deffnm {folder}/production'))
    assert jobs.get_job(folder.name)['status'] == 'cancelling'
    assert signals == [(34, jobs.signal.SIGTERM)]
    with pytest.raises(ValueError, match='Only interrupted'):
        jobs.resume_job(folder.name)
    assert signals == [(34, jobs.signal.SIGTERM)]
    monkeypatch.setattr(jobs.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=''))
    observed = jobs.get_job(folder.name)
    assert observed['status'] == 'cancelled' and 'native_pid' not in observed


def test_output_processing_retry_reuses_reserved_dataset_identity(checkpoint, monkeypatch, tmp_path):
    import os
    import mdtraj as md
    import backend.worker as worker_module
    folder, _ = checkpoint
    job = {'id': folder.name, 'name': 'Audit', 'engine': 'openmm', 'status': 'running', 'worker_pid': os.getpid(), 'config': {'dataset_id': 'source', 'name': 'Audit', 'engine': 'openmm', 'temperature_k': 300, 'seed': 2026}, 'total_steps': 1, 'completed_steps': 0, 'logs': []}
    (folder / 'input-state.json').write_text('{}')
    atomic_json(folder / 'status.json', job)
    trajectory = md.load(str(config.ROOT / 'data/datasets/ubiquitin-start/topology.pdb'))
    monkeypatch.setattr(worker_module.Worker, 'run_openmm', lambda self: trajectory)
    calls = []
    def interrupted_save(*args, **kwargs):
        calls.append(kwargs['dataset_id'])
        raise RuntimeError('Injected interruption during result conversion')
    monkeypatch.setattr(worker_module, 'save_dataset', interrupted_save)
    first = worker_module.Worker(folder.name)
    first.run()
    assert first.job['status'] == 'failed'
    assert first.job['output_dataset_id'] == calls[0]
    def completed_save(*args, **kwargs):
        calls.append(kwargs['dataset_id'])
        return {'id': kwargs['dataset_id']}
    output = tmp_path / 'output'
    output.mkdir()
    monkeypatch.setattr(worker_module, 'save_dataset', completed_save)
    monkeypatch.setattr(worker_module, 'dataset_dir', lambda dataset_id: output)
    second = worker_module.Worker(folder.name)
    second.run()
    assert second.job['status'] == 'completed'
    assert calls == [calls[0], calls[0]]


def test_cancel_during_native_child_registration_still_cleans_child(tmp_path, monkeypatch):
    import backend.worker as worker_module
    events = []
    class Child:
        pid = 123
        returncode = None
        def poll(self): return self.returncode
        def send_signal(self, sig): events.append(('signal', sig)); self.returncode = 0
        def wait(self, timeout=None): events.append(('wait', timeout)); return self.returncode
    child = Child()
    worker = object.__new__(worker_module.Worker)
    worker.folder = tmp_path
    worker.job = {}
    worker.provenance = {'commands': []}
    worker.update = lambda **fields: None
    monkeypatch.setattr(worker_module.subprocess, 'Popen', lambda *args, **kwargs: child)
    def fail_registration(path, payload):
        if payload.get('native_pid'):
            raise worker_module.Cancelled('Injected stop during ownership registration')
    monkeypatch.setattr(worker_module, 'atomic_json', fail_registration)
    with pytest.raises(worker_module.Cancelled):
        worker.run_command(['gmx', 'mdrun'], 'Audit')
    assert events == [('signal', jobs.signal.SIGTERM), ('wait', 30)]
    assert 'native_pid' not in worker.job


def test_completed_worker_generation_is_not_overwritten_by_stale_poll(checkpoint, monkeypatch):
    folder, _ = checkpoint
    job = {'id': folder.name, 'engine': 'openmm', 'status': 'running', 'created_at': '2026-09-11T00:00:00+00:00', 'worker_pid': 12}
    atomic_json(folder / 'status.json', job)
    def exited_after_final_write(pid, sig):
        atomic_json(folder / 'status.json', {**job, 'status': 'completed', 'dataset_id': 'result'})
        raise ProcessLookupError()
    monkeypatch.setattr(jobs.os, 'kill', exited_after_final_write)
    assert jobs.get_job(folder.name)['status'] == 'completed'
    assert json.loads((folder / 'status.json').read_text())['dataset_id'] == 'result'


def test_poll_and_resume_transitions_share_a_lock(checkpoint, monkeypatch):
    import threading
    folder, _ = checkpoint
    job = {'id': folder.name, 'engine': 'openmm', 'status': 'running', 'created_at': '2026-09-11T00:00:00+00:00', 'worker_pid': 12}
    atomic_json(folder / 'status.json', job)
    probed, release, replacement_written = threading.Event(), threading.Event(), threading.Event()
    def stale_probe(pid, sig):
        probed.set()
        assert release.wait(2)
        raise ProcessLookupError()
    monkeypatch.setattr(jobs.os, 'kill', stale_probe)
    poll = threading.Thread(target=lambda: jobs.get_job(folder.name))
    def replace_worker():
        with jobs._lock:
            atomic_json(folder / 'status.json', {**job, 'status': 'queued', 'worker_pid': 34, 'restarts': [{'request': 'new'}]})
            replacement_written.set()
    poll.start()
    assert probed.wait(2)
    replacement = threading.Thread(target=replace_worker)
    replacement.start()
    try:
        assert not replacement_written.wait(.1)
    finally:
        release.set()
    poll.join(2); replacement.join(2)
    assert not poll.is_alive() and not replacement.is_alive()
    persisted = json.loads((folder / 'status.json').read_text())
    assert persisted['worker_pid'] == 34 and persisted['status'] == 'queued'
    assert persisted['restarts'] == [{'request': 'new'}]


def test_terminal_job_with_live_native_child_cannot_resume(checkpoint, monkeypatch):
    folder, _ = checkpoint
    job = {'id': folder.name, 'engine': 'gromacs', 'status': 'failed', 'created_at': '2026-09-11T00:00:00+00:00', 'worker_pid': 12, 'native_pid': 34}
    atomic_json(folder / 'status.json', job)
    def kill(pid, sig):
        if pid == 12: raise ProcessLookupError()
    monkeypatch.setattr(jobs.os, 'kill', kill)
    monkeypatch.setattr(jobs.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=f'gmx mdrun -deffnm {folder}/production'))
    assert jobs.get_job(folder.name)['status'] == 'cancelling'
    with pytest.raises(ValueError, match='Only interrupted'):
        jobs.resume_job(folder.name)
