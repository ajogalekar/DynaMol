"""HTTP submission cannot bypass the resource limits displayed by readiness."""
from types import SimpleNamespace

import mdtraj as md
import pytest
from fastapi.testclient import TestClient

from backend import config, jobs, readiness, resources, storage
from backend.main import app
from backend.models import SimulationConfig


@pytest.fixture
def input_system(tmp_path, monkeypatch):
    for key, name in [('DATA_ROOT', ''), ('DATASETS_DIR', 'datasets'), ('JOBS_DIR', 'jobs')]:
        path = tmp_path / name
        path.mkdir(exist_ok=True)
        monkeypatch.setattr(config, key, path)
    metadata = storage.save_dataset(md.load(str(config.ROOT / 'data/datasets/ubiquitin-start/topology.pdb')), 'Resource audit', 'software fixture', 'Resource admission tests')
    monkeypatch.setattr(jobs, 'health', lambda: {'engines': [{'id': 'openmm', 'available': True}]})
    monkeypatch.setattr(jobs, '_launch_worker', lambda *args: pytest.fail('A rejected job must never start a worker'))
    monkeypatch.setattr(resources.shutil, 'disk_usage', lambda _: SimpleNamespace(free=10 * 1024 ** 3))
    return metadata


@pytest.mark.parametrize('limit,expected', [('atoms', 'atom limit'), ('coordinates', 'coordinate limit'), ('disk', 'disk space')])
def test_resource_blockers_match_readiness_and_actual_http_submission(input_system, monkeypatch, limit, expected):
    settings = SimulationConfig(dataset_id=input_system['id'], duration_ps=.02, report_interval=5)
    if limit == 'atoms':
        monkeypatch.setattr(config, 'MAX_ATOMS', input_system['n_atoms'] - 1)
    elif limit == 'coordinates':
        monkeypatch.setattr(config, 'MAX_COORD_BYTES', input_system['n_atoms'] * 3 * 4)
    else:
        monkeypatch.setattr(resources.shutil, 'disk_usage', lambda _: SimpleNamespace(free=100))
    before = {str(path): path.read_bytes() for path in config.DATA_ROOT.rglob('*') if path.is_file()}
    with pytest.raises(resources.ResourceLimitError) as rejected:
        jobs.validate_simulation(settings)
    assert expected in str(rejected.value)
    displayed = readiness.readiness(input_system['id'], readiness.ReadinessRequest(mode='simulation', settings=settings.model_dump()))
    assert not displayed['ready']
    assert displayed['blockers'] == rejected.value.blockers
    assert displayed['resources'] == rejected.value.resources
    with TestClient(app) as client:
        response = client.post('/api/jobs', json=settings.model_dump())
    assert response.status_code == 422 and expected in response.json()['detail']
    assert list(config.JOBS_DIR.iterdir()) == []
    assert before == {str(path): path.read_bytes() for path in config.DATA_ROOT.rglob('*') if path.is_file()}


def test_valid_resource_check_is_read_only_and_returns_shared_estimate(input_system):
    settings = SimulationConfig(dataset_id=input_system['id'], duration_ps=.02, report_interval=5)
    before = {str(path): path.read_bytes() for path in config.DATA_ROOT.rglob('*') if path.is_file()}
    validated = jobs.validate_simulation(settings)
    displayed = readiness.readiness(input_system['id'], readiness.ReadinessRequest(mode='simulation', settings=settings.model_dump()))
    assert displayed['ready']
    assert displayed['resources'] == validated['resources']
    assert displayed['resources']['saved_frames'] == 3
    assert before == {str(path): path.read_bytes() for path in config.DATA_ROOT.rglob('*') if path.is_file()}


def test_unavailable_storage_check_does_not_silently_start(input_system, monkeypatch):
    def unavailable(_):
        raise OSError('Injected unavailable volume')
    monkeypatch.setattr(resources.shutil, 'disk_usage', unavailable)
    with TestClient(app) as client:
        response = client.post('/api/jobs', json={'dataset_id': input_system['id']})
    assert response.status_code == 422 and 'Resource availability could not be checked' in response.json()['detail']
    assert list(config.JOBS_DIR.iterdir()) == []
