"""Read-only numerical, native-file and result-handoff checks for fresh native runs."""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
DATA = ROOT / 'data/integration-checks/final-release-2026-09-12/engines'
os.environ['DYNAMOL_DATA_DIR'] = str(DATA)
os.environ['DYNAMOL_GMX'] = str(ROOT / '.gromacs/bin/gmx')
os.environ['PATH'] = '/usr/bin:/bin:/usr/sbin:/sbin'
sys.path.insert(0, str(ROOT))
from fastapi.testclient import TestClient
from backend.main import app
from backend import storage
import numpy as np

native = json.loads((OUT / 'native-recovery.json').read_text())
report = {'scope': 'Fresh native 10 ps explicit-water runs; result opening/download via real HTTP routes in an isolated in-process client. Browser mode transition is covered separately.', 'engines': {}}
client = TestClient(app, base_url='http://127.0.0.1')

for engine, run in native['engines'].items():
    jid = run['job_id']
    folder = DATA / 'jobs' / jid
    target = OUT / f'fresh-{engine}-validation.json'
    validation = subprocess.run([sys.executable, str(ROOT/'scripts/validate_pressure_test.py'), '--data-root', str(DATA), '--job', jid, '--output', str(target)], text=True, capture_output=True, timeout=120)
    print(engine, validation.stdout, flush=True)
    result = {'artifact_validation_returncode': validation.returncode, 'checks': {}}
    report['engines'][engine] = result

    job = client.get(f'/api/jobs/{jid}')
    assert job.status_code == 200 and job.json()['status'] == 'completed', job.text
    did = job.json()['dataset_id']
    meta = client.get(f'/api/datasets/{did}')
    assert meta.status_code == 200
    meta = meta.json()
    physical = storage.load_physical(did)
    topology = client.get(f'/api/datasets/{did}/topology')
    coordinates = client.get(f'/api/datasets/{did}/coordinates')
    result['checks']['completed_output_opens'] = topology.status_code == 200 and coordinates.status_code == 200 and len(coordinates.content) == meta['n_frames'] * meta['n_atoms'] * 12
    result['checks']['metadata_frame_times_match'] = meta['n_frames'] == physical.n_frames and meta['n_atoms'] == physical.n_atoms

    snapshot = client.get(f'/api/jobs/{jid}/measurements')
    assert snapshot.status_code == 200
    series = snapshot.json()['measurements']
    result['checks']['measurement_definitions_retained'] = {item['id'] for item in series} == {'backbone-distance','backbone-angle','backbone-dihedral'}
    for item in series:
        response = client.post(f'/api/datasets/{did}/measurements', json={'kind':item['kind'], 'atoms':item['output_atoms']})
        assert response.status_code == 200, response.text
        actual = np.asarray(response.json()['values'])
        expected = np.asarray(item['values'])
        difference = actual - expected
        if item['kind'] == 'dihedral':
            difference = (difference + 180) % 360 - 180
        result['checks'][item['id'] + '_posthoc_matches_live'] = actual.shape == expected.shape and np.all(np.isfinite(actual)) and float(np.max(np.abs(difference))) < 1e-4

    archive = client.get(f'/api/jobs/{jid}/download')
    assert archive.status_code == 200, archive.text
    with zipfile.ZipFile(io.BytesIO(archive.content)) as bundle:
        result['checks']['download_zip_integrity'] = bundle.testzip() is None
        for name in ['config.json','provenance.json','energies.csv', 'trajectory.xtc' if engine == 'openmm' else 'production.xtc']:
            result['checks']['download_' + name] = bundle.read(name) == (folder / name).read_bytes()
        result['download'] = {'bytes':len(archive.content), 'sha256':hashlib.sha256(archive.content).hexdigest(), 'files':len(bundle.namelist())}
    provenance = json.loads((folder / 'provenance.json').read_text())
    result['application'] = provenance['application']
    result['checks']['current_app_version'] = provenance['application'] == 'DynaMol 0.1.1'
    result['passed'] = validation.returncode == 0 and all(result['checks'].values())
    print(engine, 'handoff', result['passed'], len(result['checks']), 'checks', flush=True)

report['passed'] = native['passed'] and all(r['passed'] for r in report['engines'].values())
(OUT / 'fresh-handoff-validation.json').write_text(json.dumps(report, indent=2, default=lambda v: bool(v) if isinstance(v,np.bool_) else str(v)) + '\n')
raise SystemExit(0 if report['passed'] else 1)
