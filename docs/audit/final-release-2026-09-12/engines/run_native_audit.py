#!/usr/bin/env python3
"""Native checkpoint continuity checks, deliberately isolated from user projects."""
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
os.environ['DYNAMOL_DATA_DIR'] = str(ROOT / 'data/integration-checks/final-release-2026-09-12/engines')
os.environ['DYNAMOL_CPU_THREADS'] = '2'
os.environ['PATH'] = '/usr/bin:/bin:/usr/sbin:/sbin'
os.environ['DYNAMOL_GMX'] = str(ROOT / '.gromacs/bin/gmx')
from backend import config, jobs, storage
from backend.models import SimulationConfig
from backend.recovery import validate_manifest
import mdtraj as md
import numpy as np

REPORT = ROOT / 'docs/audit/final-release-2026-09-12/engines/native-recovery.json'
report = {'passed': False, 'scope': 'Fresh 10 ps explicit-water ubiquitin NVT runs in each engine, with worker interruption and cooperative stop/resume. Software, numerical health and artifact continuity audit; no equilibration, convergence, biological stability, or bitwise reproducibility claim.', 'started_at': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), 'seed': 20260912, 'engines': {}}
REPORT.parent.mkdir(parents=True, exist_ok=True)
if REPORT.exists():
    report['engines'] = json.loads(REPORT.read_text()).get('engines', {})

def save():
    REPORT.write_text(json.dumps(report, indent=2) + '\n')

def wait(job, predicate, timeout=900):
    started = time.monotonic()
    while time.monotonic() - started < timeout:
        state = jobs.get_job(job['id'])
        if predicate(state):
            return state
        if state['status'] in {'failed', 'completed', 'cancelled', 'interrupted'}:
            raise RuntimeError(f"Unexpected state: {state['status']} {state.get('error')}")
        time.sleep(.2)
    jobs.cancel_job(job['id'])
    raise TimeoutError(job['id'])

def terminal(state):
    return state['status'] in {'failed', 'completed', 'cancelled', 'interrupted'}

source = ROOT / 'data/datasets/ubiquitin-start/topology.pdb'
fixture = storage.save_dataset(md.load(str(source)), 'Final release audit ubiquitin', 'local audit fixture', 'Native restart software check')
report['input_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
atoms = fixture['atoms']
ca = [a for a in atoms if a['category'] == 'protein' and a['name'] == 'CA']
r = ca[20]
by_name = {a['name']: a['index'] for a in atoms if a['chain'] == r['chain'] and a['resid'] == r['resid']}
next_n = next(b if a == by_name['C'] else a for a,b in fixture['bonds'] if (a == by_name['C'] or b == by_name['C']) and atoms[b if a == by_name['C'] else a]['name'] == 'N')
measurements = [
    {'id':'backbone-distance','kind':'distance','atoms':[ca[10]['index'],ca[25]['index']],'label':'C-alpha separation'},
    {'id':'backbone-angle','kind':'angle','atoms':[by_name[k] for k in ('N','CA','C')],'label':'Backbone N-CA-C'},
    {'id':'backbone-dihedral','kind':'dihedral','atoms':[by_name[k] for k in ('N','CA','C')] + [next_n],'label':'Backbone psi'},
]
report['runtime'] = {'python': sys.executable, 'path': os.environ['PATH'], 'gromacs': jobs.gromacs_executable(), 'health': jobs.health()}
report['source_hashes'] = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'backend/worker.py',ROOT/'backend/jobs.py',ROOT/'backend/recovery.py',ROOT/'backend/live_measurements.py']}
try:
    for engine in sys.argv[1:] or ['openmm', 'gromacs']:
        settings = SimulationConfig(dataset_id=fixture['id'], engine=engine, name=f'Final audit {engine} 10 ps', duration_ps=10, solvent='explicit', report_interval=100, minimize=True, equilibration_steps=500, seed=20260912, measurements=measurements)
        job = jobs.submit_job(settings)
        folder = config.JOBS_DIR / job['id']
        print(engine, 'started', job['id'], flush=True)
        result = report['engines'][engine] = {'job_id': job['id'], 'config': settings.model_dump()}
        retained_native_prefix = None
        if engine == 'openmm':
            state = wait(job, lambda j: j.get('recovery', {}).get('checkpoint_saved') and j.get('recovery', {}).get('step', 0) >= 500)
            manifest = validate_manifest(folder, engine)
            result['interrupted_checkpoint_step'] = manifest['checkpoint']['step']
            original_frames = {record['name']: (folder / 'frames' / record['name']).read_bytes() for record in manifest['checkpoint']['frames']}
            result['diagnostics_available_while_running'] = (folder / 'energies.csv').is_file() and len((folder / 'energies.csv').read_text().splitlines()) > 1
            os.kill(state['worker_pid'], signal.SIGKILL)
            state = wait(job, lambda j: j['status'] == 'interrupted')
            assert state['recovery']['available']
        else:
            state = wait(job, lambda j: (folder / 'production.cpt').is_file())
            original_frames = {}
            result['interrupted_checkpoint_step'] = state['completed_steps']
            assert state.get('native_pid'), 'Native child ownership not recorded'
            result['diagnostics_available_while_running'] = (folder / 'energies.csv').is_file() and len((folder / 'energies.csv').read_text().splitlines()) > 1
            os.kill(state['worker_pid'], signal.SIGKILL)
            state = wait(job, terminal)
            assert state['status'] == 'interrupted', state
            result['orphaned_native_child_stopped_before_resume'] = True
            retained_native_prefix = md.load(str(folder / 'production.xtc'), top=str(folder / 'system.gro'))
        result['live_diagnostics_before_resume'] = (folder / 'energies.csv').is_file() and len((folder / 'energies.csv').read_text().splitlines()) > 1
        assert result['live_diagnostics_before_resume']
        protected = folder / ('system.xml' if engine == 'openmm' else 'production.tpr')
        original = protected.read_bytes()
        protected.write_bytes(original + b'\n')
        try:
            jobs.resume_job(job['id'])
        except ValueError as exc:
            assert 'changed' in str(exc)
            result['tampered_system_rejected'] = True
        else:
            raise AssertionError('Tampered system accepted')
        finally:
            protected.write_bytes(original)
        saved_runtime = (folder / 'recovery.json').read_text()
        manifest = json.loads(saved_runtime)
        manifest['runtime']['cpu_threads'] = 99
        storage.atomic_json(folder / 'recovery.json', manifest)
        try:
            jobs.resume_job(job['id'])
        except ValueError as exc:
            assert 'changed' in str(exc)
            result['different_runtime_rejected'] = True
        else:
            raise AssertionError('Changed runtime accepted')
        finally:
            (folder / 'recovery.json').write_text(saved_runtime)
        resumed = jobs.resume_job(job['id'])
        if engine == 'openmm':
            running = wait(resumed, lambda j: j['status'] == 'running' and j.get('recovery', {}).get('step', 0) >= result['interrupted_checkpoint_step'] + 500)
        else:
            running = wait(resumed, lambda j: j['status'] == 'running' and j['completed_steps'] >= result['interrupted_checkpoint_step'] + 500)
        jobs.cancel_job(job['id'])
        stopped = wait(resumed, terminal)
        assert stopped['status'] == 'cancelled', stopped
        result['cooperative_stop_checkpoint_step'] = validate_manifest(folder, engine).get('checkpoint', {}).get('step', stopped['completed_steps'])
        resumed = jobs.resume_job(job['id'])
        state = wait(resumed, terminal)
        assert state['status'] == 'completed', state.get('error')
        physical = storage.load_physical(state['dataset_id'])
        assert np.isfinite(physical.xyz).all()
        if retained_native_prefix is not None:
            count = retained_native_prefix.n_frames
            assert np.array_equal(physical.xyz[:count], retained_native_prefix.xyz)
            assert np.allclose(physical.time[:count], retained_native_prefix.time, atol=1e-6)
            result['verified_retained_prefix_frames'] = count
        assert np.all(np.diff(physical.time) > 0)
        assert len(physical.time) == len(np.unique(physical.time))
        assert abs(float(physical.time[-1]) - settings.duration_ps) < 1e-5
        with (folder / 'energies.csv').open() as stream:
            energies = list(csv.DictReader(stream))
        assert energies and all(np.isfinite([float(row[key]) for key in ('potential_kj_mol', 'kinetic_kj_mol', 'temperature_k')]).all() for row in energies)
        assert all((folder / 'frames' / name).read_bytes() == data for name, data in original_frames.items())
        provenance = json.loads((folder / 'provenance.json').read_text())
        assert len(provenance['restarts']) == 2
        if engine == 'gromacs':
            assert any('-cpi' in command and '-append' in command for command in provenance['commands'])
        result.update(status=state['status'], saved_frames=physical.n_frames, atoms=physical.n_atoms, final_time_ps=float(physical.time[-1]), unique_strictly_increasing_times=True, saved_prefix_unchanged=True, finite_energy_temperature_coordinates=True, energy_rows=len(energies), restart_provenance=provenance['restarts'], outputs_sha256={name: hashlib.sha256((folder / name).read_bytes()).hexdigest() for name in ('energies.csv','provenance.json')})
        result['elapsed_seconds'] = state['elapsed_seconds']
        result['dataset_id'] = state['dataset_id']
        result['completed_steps'] = state['completed_steps']
        print(engine, 'PASS', result['saved_frames'], 'frames', state['elapsed_seconds'], 'seconds', flush=True)
        save()
    report['passed'] = True
finally:
    save()
