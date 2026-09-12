import io
import json
import zipfile
from pathlib import Path
from types import SimpleNamespace
import pytest
from fastapi.testclient import TestClient
from backend import config, diagnostics, jobs, readiness, resources, storage
from backend.main import app

@pytest.fixture(autouse=True)
def isolated(tmp_path,monkeypatch):
    for key,name in [('DATA_ROOT',''),('DATASETS_DIR','datasets'),('JOBS_DIR','jobs')]:
        path=tmp_path/name;path.mkdir(exist_ok=True);monkeypatch.setattr(config,key,path)


def job_record():
    return {'id':'diag','engine':'openmm','status':'completed','stage':'Complete','config':{'temperature_k':300},'elapsed_seconds':20,'total_steps':1000,'completed_steps':1000,'logs':[]}

def test_diagnostics_partial_rows_missing_temperature_and_invalid_data(monkeypatch):
    monkeypatch.setattr(jobs,'get_job',lambda _:job_record()); folder=config.JOBS_DIR/'diag';folder.mkdir()
    (folder/'energies.csv').write_text('step,time_ps,potential_kj_mol,kinetic_kj_mol\n0,0,-5,2\n1,0.002,nan,2\n2,0.004,-4,3\n3,0.006,-')
    r=diagnostics.read_diagnostics('diag')
    assert r['recorded_points']==2 and r['invalid_rows']==1
    assert r['latest']['step']==2 and r['latest']['temperature_k'] is None
    assert 'temperature_k' not in r['available'] and r['estimated_remaining_seconds'] is None

def test_diagnostic_downsampling_keeps_peaks_and_endpoints(monkeypatch):
    job=job_record();job.update(status='running',completed_steps=500)
    monkeypatch.setattr(jobs,'get_job',lambda _:job); folder=config.JOBS_DIR/'diag';folder.mkdir()
    with (folder/'energies.csv').open('w') as out:
        out.write('step,time_ps,potential_kj_mol,kinetic_kj_mol,temperature_k\n')
        for i in range(1000):out.write(f'{i},{i*.002},{-100 if i==417 else -1},2,{800 if i==561 else 300}\n')
    r=diagnostics.read_diagnostics('diag',20)
    assert len(r['points'])<=20 and r['recorded_points']==1000
    assert {0,417,561,999}.issubset({p['step'] for p in r['points']})
    assert r['estimated_remaining_seconds']==20

def test_diagnostics_reject_nonfinite_fractional_or_negative_steps(monkeypatch):
    monkeypatch.setattr(jobs,'get_job',lambda _:job_record())
    folder=config.JOBS_DIR/'diag';folder.mkdir()
    (folder/'energies.csv').write_text('step,time_ps,potential_kj_mol\ninf,0,-5\n1e999,0,-5\n-1,0,-5\n.5,0,-5\n1,-.2,-5\n2,0.004,-4\n')
    with TestClient(app, base_url="http://127.0.0.1") as client:
        response=client.get('/api/jobs/diag/diagnostics')
    assert response.status_code==200
    result=response.json()
    assert result['recorded_points']==1 and result['invalid_rows']==5
    assert result['latest']['step']==2

def test_readiness_uses_same_submission_guard_and_does_not_queue(monkeypatch):
    meta={'id':'sample','n_atoms':10,'n_frames':1,'atoms':[],'preparation':{'ph':7},'solvation':{}}
    monkeypatch.setattr(storage,'get_dataset',lambda _:meta)
    monkeypatch.setattr(jobs,'list_jobs',lambda:[])
    def reject(_):raise ValueError('Exact system is not ready')
    monkeypatch.setattr(jobs,'validate_simulation',reject)
    r=readiness.readiness('sample',readiness.ReadinessRequest(mode='simulation',settings={}))
    assert not r['ready'] and 'Exact system is not ready' in r['blockers']
    assert not list(config.JOBS_DIR.iterdir())

def test_readiness_blocks_resource_overflow_and_low_disk(monkeypatch):
    meta={'id':'sample','n_atoms':100000,'n_frames':1,'atoms':[]}
    monkeypatch.setattr(storage,'get_dataset',lambda _:meta);monkeypatch.setattr(jobs,'list_jobs',lambda:[])
    monkeypatch.setattr(jobs,'validate_simulation',lambda _: {})
    monkeypatch.setattr(resources.shutil,'disk_usage',lambda _:SimpleNamespace(free=100))
    r=readiness.readiness('sample',readiness.ReadinessRequest(mode='simulation',settings={'duration_ps':10,'report_interval':1}))
    assert not r['ready'] and any('coordinate limit' in b for b in r['blockers']) and any('disk space' in b for b in r['blockers'])

@pytest.mark.parametrize('settings',[{'timestep_fs':None},{'duration_ps':'bad'},{'report_interval':0}])
def test_invalid_readiness_numbers_return_blockers_without_server_error(monkeypatch,settings):
    monkeypatch.setattr(storage,'get_dataset',lambda _:{'id':'x','n_atoms':10,'n_frames':1,'atoms':[]});monkeypatch.setattr(jobs,'list_jobs',lambda:[])
    r=readiness.readiness('x',readiness.ReadinessRequest(mode='simulation',settings=settings))
    assert not r['ready']

def test_job_download_is_private_snapshot_and_cleans_tempfile(monkeypatch,tmp_path):
    folder=config.JOBS_DIR/'diag';folder.mkdir();(folder/'energies.csv').write_text('original')
    (folder/'dynamol-output.zip').write_text('stale archive')
    monkeypatch.setattr(jobs,'get_job',lambda _:job_record())
    with TestClient(app, base_url="http://127.0.0.1") as client:
        first=client.get('/api/jobs/diag/download');second=client.get('/api/jobs/diag/download')
    for response in (first,second):
        assert response.status_code==200
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            assert z.namelist()==['energies.csv'] and z.read('energies.csv')==b'original'
    assert not (folder/'checkpoint.chk').exists()
