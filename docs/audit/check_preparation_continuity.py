"""One-step periodic continuation from actual tiny preparation worker output."""
from pathlib import Path
import hashlib,json,shutil,sys,tempfile,time
import numpy as np
from openmm import app,unit
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from backend import config,sources,solvent,storage,jobs
from backend.models import SimulationConfig
report={'scope':'One-step OpenMM continuation with minimization disabled solely to test exact saved input geometry/box continuity, not physical stability or equilibrium.'}
evidence=ROOT/'docs/audit/preparation-worker-evidence'
try:
 prep_report=json.loads((ROOT/'docs/audit/preparation-pipeline-checks.json').read_text())
 row=next(x for x in prep_report['results'] if x['name'].startswith('Repreparing'))['evidence']
 pp=evidence/row['ph2_job']/'prepared.pdb';state=json.loads((evidence/row['ph2_job']/'preparation.json').read_text())
 with tempfile.TemporaryDirectory(prefix='dynamol-preview-continuity-') as tmp:
  config.DATA_ROOT=Path(tmp);config.DATASETS_DIR=Path(tmp)/'datasets';config.JOBS_DIR=Path(tmp)/'jobs';config.DATASETS_DIR.mkdir();config.JOBS_DIR.mkdir()
  m=sources.import_structure(pp);folder=storage.dataset_dir(m['id']);shutil.copy2(pp,folder/'prepared.pdb');m['preparation']=state;storage.atomic_json(folder/'metadata.json',m)
  preview=solvent.solvate_dataset(m['id'],padding_nm=1,seed=2026,ph=2)
  inputpath=storage.dataset_dir(preview['id'])/'prepared.pdb';inputbytes=inputpath.read_bytes();inp=app.PDBFile(str(inputpath));ix=np.asarray(inp.positions.value_in_unit(unit.nanometer),dtype=np.float32);ibox=np.asarray(inp.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer),dtype=np.float32)
  j=jobs.submit_job(SimulationConfig(dataset_id=preview['id'],engine='openmm',solvent='explicit',minimize=False,equilibration_steps=0,duration_ps=.002,report_interval=1))
  assert (config.JOBS_DIR/j['id']/'input.pdb').read_bytes()==inputbytes
  start=time.monotonic()
  while time.monotonic()-start<40:
   j=jobs.get_job(j['id'])
   if j['status'] not in {'queued','running','cancelling'}:break
   time.sleep(.1)
  if j['status']!='completed':
   jobs.cancel_job(j['id']);raise AssertionError(j)
  out=storage.load_physical(j['dataset_id']);metadata=storage.get_dataset(j['dataset_id'])
  assert out.n_atoms==preview['n_atoms'] and out.n_frames==2
  assert np.array_equal(out.xyz[0],ix)
  assert np.allclose(out.unitcell_vectors[0],ibox,rtol=0,atol=1e-6)
  p=json.loads((config.JOBS_DIR/j['id']/'provenance.json').read_text())
  assert p['input_solvation'] and any('no water or ions added' in x for x in p['preparation'])
  assert metadata.get('preparation',{}).get('ph')==2 and metadata.get('solvation'), 'Result dataset dropped prepared state/solvent markers'
  exact=storage.dataset_dir(metadata['id'])/'prepared.pdb';assert exact.is_file(), 'Result lacks exact first-frame prepared topology'
  assert metadata.get('parent_dataset_id')==preview['id']
  try:jobs.submit_job(SimulationConfig(dataset_id=metadata['id'],engine='gromacs',solvent='explicit',duration_ps=.002,equilibration_steps=0))
  except ValueError as exc:assert 'protonation' in str(exc) or 'prepared' in str(exc)
  else:raise AssertionError('Result dataset bypassed prepared GROMACS guard')
  report.update(status='pass',job_id=j['id'],atoms=out.n_atoms,frames=out.n_frames,times_ps=out.time.tolist(),input_pdb_sha256=hashlib.sha256(inputbytes).hexdigest(),first_frame_exact_saved_preview_coordinates=True,box_retained_within_1e_6_nm=True,atom_count_unchanged=True,actual_preparation_provenance_carried=True,result_dataset_state_and_exact_first_frame_retained=True,result_still_protected_from_gromacs_reset=True,preparation_log=p['preparation'])
  target=evidence/j['id'];target.mkdir()
  for path in (config.JOBS_DIR/j['id']).iterdir():
   if path.is_file() and path.suffix in {'.json','.pdb','.log','.csv'} and path.stat().st_size<2_000_000:(target/path.name).write_bytes(path.read_bytes())
except Exception as exc:report.update(status='fail',error=repr(exc))
(ROOT/'docs/audit/preparation-continuity-checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));sys.exit(report['status']!='pass')
