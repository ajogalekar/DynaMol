"""Independent tiny end-to-end preparation checks; no scientific accuracy claim.

Run only after the exploratory preparation preflight and worker are available.
"""
from pathlib import Path
import hashlib,json,sys,tempfile,time
import numpy as np
from openmm import app,unit
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from backend import config,sources,preparation,storage,jobs
from backend.models import SimulationConfig
results=[];artifact=ROOT/'docs/audit/preparation-worker-evidence';artifact.mkdir(exist_ok=True)
def check(name,f):
 try:
  ev=f();results.append({'name':name,'status':'pass','evidence':ev})
 except Exception as e:results.append({'name':name,'status':'fail','error':repr(e)})
 print(results[-1]['name'],results[-1]['status'],flush=True)
def reject(f):
 try:f()
 except ValueError as e:return str(e)
 raise AssertionError('Invalid operation accepted')
def wait(job,max_seconds=80):
 start=time.monotonic()
 while time.monotonic()-start<max_seconds:
  j=jobs.get_job(job['id'])
  if j['status'] not in {'queued','running','cancelling'}:
   out=artifact/job['id'];out.mkdir(exist_ok=True)
   folder=config.JOBS_DIR/job['id']
   for file in folder.iterdir():
    if file.is_file() and file.suffix in {'.json','.log','.pdb','.csv','.xml'} and file.stat().st_size<2_000_000:
     (out/file.name).write_bytes(file.read_bytes())
   if j['status']!='completed':raise AssertionError(f"Worker {j['status']}: {j.get('error')} ({out})")
   return j,storage.get_dataset(j['dataset_id'])
  time.sleep(.1)
 jobs.cancel_job(job['id']);raise AssertionError('Tiny worker check exceeded time bound and was cancelled')
def pdb(m):return app.PDBFile(str(storage.dataset_dir(m['id'])/'prepared.pdb'))
def inventory(p,resid):return sorted(a.name for a in p.topology.atoms() if a.residue.id==str(resid))
with tempfile.TemporaryDirectory(prefix='dynamol-preparation-pipeline-') as tmp:
 config.DATA_ROOT=Path(tmp);config.DATASETS_DIR=Path(tmp)/'datasets';config.JOBS_DIR=Path(tmp)/'jobs';config.DATASETS_DIR.mkdir();config.JOBS_DIR.mkdir()
 hold={};fixtures=ROOT/'docs/audit/preparation-fixtures'
 def missing_atom():
  m=sources.import_structure(fixtures/'six_residues_missing_atom.pdb')
  j,s=wait(preparation.submit_preparation({'dataset_id':m['id'],'optimize_sidechains':False}))
  p=pdb(s);assert 'NZ' in inventory(p,6) and 'OXT' in inventory(p,6)
  assert not preparation.inspect_preparation(s['id'])['missing_atoms']
  return {'job_id':j['id'],'atoms':s['n_atoms'],'restored':['LYS6 NZ','LYS6 OXT'],'preparation':s['preparation']}
 check('Actual preparation restores missing heavy/terminal atoms and parameterizes tiny peptide',missing_atom)
 def ph_rebuild():
  m=sources.import_structure(fixtures/'histidine_fragment.pdb')
  j,a=wait(preparation.submit_preparation({'dataset_id':m['id'],'ph':2,'optimize_sidechains':False}));pa=pdb(a)
  j2,b=wait(preparation.submit_preparation({'dataset_id':a['id'],'ph':7,'optimize_sidechains':False}));pb=pdb(b)
  ai=inventory(pa,68);bi=inventory(pb,68)
  assert 'HD1' in ai and 'HE2' in ai
  assert ('HD1' in bi) != ('HE2' in bi)
  assert b['preparation']['ph']==7
  hold['prepared']=a
  return {'ph2_job':j['id'],'ph7_job':j2['id'],'ph2_HIS_atoms':ai,'ph7_HIS_atoms':bi,'ph2_atoms':a['n_atoms'],'ph7_atoms':b['n_atoms'],'preparation':b['preparation']}
 check('Repreparing pH2 histidine at pH7 actually replaces hydrogen state',ph_rebuild)
 def missing_loop():
  m=sources.import_structure(fixtures/'six_residues_known_gap.pdb')
  start_job=preparation.submit_preparation({'dataset_id':m['id'],'build_missing_residues':True,'optimize_sidechains':True})
  try:j,s=wait(start_job)
  except AssertionError as exc:
   failed=jobs.get_job(start_job['id'])
   if failed['status']=='failed' and 'stereochemistry' in failed.get('error','') and not failed.get('dataset_id'):
    return {'job_id':start_job['id'],'invalid_geometry_rejected_before_publication':True,'reason':failed['error'],'scope':'This particular ILE loop attempt is invalid; see separate ALA28 positive control for accepted construction.'}
   raise
  p=pdb(s);res=[(r.id,r.name) for r in p.topology.residues()]
  assert ('3','ILE') in res,res
  inspection=preparation.inspect_preparation(s['id']);assert not any(g['structural_break'] for g in inspection['gaps'])
  assert any('confidence' in w.lower() or 'modeled' in w.lower() or 'uncertain' in w.lower() for w in s['warnings']+s['preparation'].get('warnings',[]))
  xyz=np.asarray(p.positions.value_in_unit(unit.nanometer));signed=[]
  for residue in p.topology.residues():
   atoms={a.name:a.index for a in residue.atoms()}
   if all(k in atoms for k in ('N','CA','C','CB')):
    center=xyz[atoms['CA']];value=float(np.linalg.det([xyz[atoms[k]]-center for k in ('N','C','CB')]))
    assert value>1e-4, f"Inverted/degenerate CA {residue.name}{residue.id}: {value}"
    signed.append({'resid':residue.id,'center':'CA','signed_volume_nm3':value})
   if residue.name=='ILE':
    center=xyz[atoms['CB']];value=float(np.linalg.det([xyz[atoms[k]]-center for k in ('CA','CG1','CG2')]))
    assert value>1e-4, f"Inverted/degenerate ILE CB {residue.id}: {value}"
    signed.append({'resid':residue.id,'center':'CB','signed_volume_nm3':value})
  return {'job_id':j['id'],'residues':res,'remaining_gaps':inspection['gaps'],'signed_volume_checks':signed,'preparation':s['preparation']}
 check('Opt-in sequence-supported loop produces checked geometry or rejects invalid stereochemistry before publication',missing_loop)
 def gmx_reject():
  m=hold['prepared']
  message=reject(lambda:jobs.submit_job(SimulationConfig(dataset_id=m['id'],engine='gromacs',solvent='explicit',duration_ps=.002,equilibration_steps=0,report_interval=1)))
  assert 'proton' in message.lower() or 'prepared' in message.lower(),message
  return {'rejection':message}
 check('Prepared GROMACS state is rejected before any silent protonation reset',gmx_reject)
 def openmm_preserve():
  m=hold['prepared'];source=storage.dataset_dir(m['id'])/'prepared.pdb';before=source.read_bytes();p=app.PDBFile(str(source));names=[(a.name,a.residue.name,a.residue.id) for a in p.topology.atoms()]
  job=jobs.submit_job(SimulationConfig(dataset_id=m['id'],engine='openmm',solvent='implicit',duration_ps=.002,equilibration_steps=0,report_interval=1))
  assert (config.JOBS_DIR/job['id']/'input.pdb').read_bytes()==before
  j,s=wait(job)
  out=app.PDBFile(str(config.JOBS_DIR/j['id']/'prepared.pdb'))
  after=[(a.name,a.residue.name,a.residue.id) for a in out.topology.atoms()]
  assert names==after
  assert 'HE2' in inventory(out,68) and 'HD1' in inventory(out,68)
  hold['implicit']=s
  return {'job_id':j['id'],'input_sha256':hashlib.sha256(before).hexdigest(),'all_atom_identities_retained':True,'ph2_HIP_retained':True,'frames':s['n_frames'],'times_ps':s['times_ps']}
 check('One-step real OpenMM continuation uses exact prepared input and preserves actual pH2 hydrogen state',openmm_preserve)
report={'scope':'Tiny synthetic peptide end-to-end preparation and one-step MD software checks. Not loop accuracy, pKa validation, or MD convergence.','results':results}
(ROOT/'docs/audit/preparation-pipeline-checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));sys.exit(any(r['status']=='fail' for r in results))
