"""Recovery admission uses native counts/current limits and retains rejection reasons."""
import json
from types import SimpleNamespace

import pytest

from backend import config, jobs, recovery, resources
from backend.storage import atomic_json
from test_recovery import checkpoint


@pytest.fixture
def recoverable(checkpoint, monkeypatch):
    folder, manifest = checkpoint
    monkeypatch.setattr(resources.shutil, 'disk_usage', lambda _: SimpleNamespace(free=10 * 1024**3))
    return folder, manifest


@pytest.mark.parametrize('source', ['provenance', 'check-log'])
def test_stereochemistry_rejection_blocks_inspection_and_resume_without_mutation(recoverable, monkeypatch, source):
    folder, manifest = recoverable
    report = {'passed': False, 'stage': 'Production step 500', 'violations': [{'residue': 'ILE', 'center': 'CA'}]}
    if source == 'provenance':
        atomic_json(folder / 'provenance.json', {'stereochemistry_failure': report})
    else:
        atomic_json(folder / 'stereochemistry-checks.json', {'checks': [report]})
    job = {'id': folder.name, 'engine': 'openmm', 'status': 'failed', 'logs': []}
    monkeypatch.setattr(jobs, 'get_job', lambda _: dict(job))
    monkeypatch.setattr(jobs, 'list_jobs', lambda: [job])
    monkeypatch.setattr(jobs, '_launch_worker', lambda *args: pytest.fail('A chemistry-rejected checkpoint must not launch'))
    before = {str(p): p.read_bytes() for p in folder.rglob('*') if p.is_file()}
    for verify in [False, True]:
        info = recovery.recovery_info(job, verify=verify)
        assert not info['available']
        assert 'Repair and inspect' in info['reason']
    with pytest.raises(ValueError, match='cannot be resumed'):
        jobs.resume_job(folder.name)
    assert before == {str(p): p.read_bytes() for p in folder.rglob('*') if p.is_file()}


def test_early_stereo_failure_reason_precedes_missing_checkpoint(recoverable):
    folder, _ = recoverable
    (folder / 'recovery.json').unlink()
    atomic_json(folder / 'provenance.json', {'stereochemistry_failure': {'stage': 'After minimization'}})
    info = recovery.recovery_info({'id': folder.name, 'engine': 'openmm', 'status': 'failed'}, verify=True)
    assert not info['available'] and 'stereochemistry' in info['reason']


@pytest.mark.parametrize('engine', ['openmm', 'gromacs'])
def test_current_atom_profile_is_enforced_from_native_topology(recoverable, monkeypatch, engine):
    folder, manifest = recoverable
    count = manifest['resource_requirements']['atoms']
    if engine == 'gromacs':
        # Integrity-only native checkpoint fixture; no engine execution occurs.
        (folder / 'system.gro').write_text(f'Native particle-count fixture\n{count}\n')
        (folder / 'production.tpr').write_bytes(b'TPR integrity fixture')
        (folder / 'topol.top').write_text('Topology integrity fixture')
        (folder / 'production.cpt').write_bytes(b'GROMACS checkpoint integrity fixture')
        manifest = recovery.create_manifest(folder, engine)
    monkeypatch.setattr(config, 'MAX_ATOMS', count - 1)
    with pytest.raises(resources.ResourceLimitError, match='current resource limits') as rejected:
        recovery.validate_manifest(folder, engine)
    assert rejected.value.resources['input_atoms'] == count
    assert rejected.value.resources['estimated_atoms'] == count  # no second solvent estimate
    job = {'id': folder.name, 'engine': engine, 'status': 'cancelled', 'logs': []}
    assert not recovery.recovery_info(job)['available']
    monkeypatch.setattr(jobs, 'get_job', lambda _: dict(job))
    monkeypatch.setattr(jobs, 'list_jobs', lambda: [job])
    monkeypatch.setattr(jobs, '_launch_worker', lambda *args: pytest.fail('Oversized resume must not launch'))
    with pytest.raises(resources.ResourceLimitError):
        jobs.resume_job(folder.name)
    monkeypatch.setattr(config, 'MAX_ATOMS', count)
    assert recovery.validate_manifest(folder, engine) == manifest


@pytest.mark.parametrize('limit', ['coordinates', 'frames', 'disk'])
def test_independent_current_recovery_limits_cannot_be_bypassed(recoverable, monkeypatch, limit):
    folder, manifest = recoverable
    requirements = manifest['resource_requirements']
    if limit == 'coordinates':
        monkeypatch.setattr(config, 'MAX_COORD_BYTES', requirements['coordinate_bytes'] - 1)
    elif limit == 'frames':
        monkeypatch.setattr(config, 'MAX_FRAMES', requirements['saved_frames'] - 1)
    else:
        monkeypatch.setattr(resources.shutil, 'disk_usage', lambda _: SimpleNamespace(free=100))
    with pytest.raises(resources.ResourceLimitError):
        recovery.validate_manifest(folder, 'openmm')


def test_legacy_manifest_without_cached_resources_uses_same_exact_validation(recoverable, monkeypatch):
    folder, manifest = recoverable
    count = manifest.pop('resource_requirements')['atoms']
    atomic_json(folder / 'recovery.json', manifest)
    before = (folder / 'recovery.json').read_bytes()
    assert recovery.validate_manifest(folder, 'openmm') == manifest
    monkeypatch.setattr(config, 'MAX_ATOMS', count - 1)
    with pytest.raises(resources.ResourceLimitError):
        recovery.validate_manifest(folder, 'openmm')
    assert (folder / 'recovery.json').read_bytes() == before


def test_successful_old_stereo_checks_do_not_block_recovery(recoverable):
    folder, manifest = recoverable
    atomic_json(folder / 'stereochemistry-checks.json', {'checks': [{'stage': 'Before minimization', 'passed': True, 'violations': []}]})
    atomic_json(folder / 'provenance.json', {'purpose': 'Interrupted software fixture'})
    assert recovery.validate_manifest(folder, 'openmm') == manifest
