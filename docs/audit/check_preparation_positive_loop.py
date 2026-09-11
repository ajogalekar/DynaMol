"""One bounded positive control: synthetic missing ALA28 in 1UBQ residues 26–30."""
from pathlib import Path
import hashlib,json,sys,tempfile,time
import numpy as np
from openmm import app,unit
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from backend import config,sources,preparation,storage,jobs
fixture=ROOT/'docs/audit/preparation-fixtures/known_alanine_gap.pdb'
lines=(ROOT/'data/sources/1UBQ.pdb').read_text().splitlines();atoms=[s for s in lines if s.startswith('ATOM  ') and 26<=int(s[22:26])<=30 and int(s[22:26])!=28]
fixture.write_text('SEQRES   1 A    5  VAL LYS ALA LYS ILE\n'+'\n'.join(atoms)+'\nTER\nEND\n')
r={'scope':'Single fixed-seed positive-control attempt; a synthetic 1UBQ 26–30 fragment missing known ALA28, not a real observed truncated construct or accurate-loop benchmark.','seed':2026,'optimize_sidechains':False,'fixture':str(fixture.relative_to(ROOT)),'fixture_sha256':hashlib.sha256(fixture.read_bytes()).hexdigest()}
try:
 with tempfile.TemporaryDirectory(prefix='dynamol-positive-loop-') as tmp:
  config.DATA_ROOT=Path(tmp);config.DATASETS_DIR=Path(tmp)/'datasets';config.JOBS_DIR=Path(tmp)/'jobs';config.DATASETS_DIR.mkdir();config.JOBS_DIR.mkdir()
  m=sources.import_structure(fixture);i=preparation.inspect_preparation(m['id']);assert any(x['residues']==['ALA'] for x in i['missing_residues'])
  j=preparation.submit_preparation({'dataset_id':m['id'],'build_missing_residues':True,'optimize_sidechains':False,'seed':2026})
  start=time.monotonic()
  while time.monotonic()-start<45:
   j=jobs.get_job(j['id'])
   if j['status'] not in {'queued','running','cancelling'}:break
   time.sleep(.1)
  if j['status'] in {'queued','running','cancelling'}:jobs.cancel_job(j['id'])
  out=ROOT/'docs/audit/preparation-worker-evidence'/j['id'];out.mkdir()
  for f in (config.JOBS_DIR/j['id']).iterdir():
   if f.is_file() and f.suffix in {'.json','.pdb','.log','.xml'} and f.stat().st_size<2_000_000:(out/f.name).write_bytes(f.read_bytes())
  r['job_id']=j['id'];assert j['status']=='completed',j.get('error',j)
  result=storage.get_dataset(j['dataset_id']);p=app.PDBFile(str(storage.dataset_dir(result['id'])/'prepared.pdb'));xyz=np.asarray(p.positions.value_in_unit(unit.nanometer));volumes=[]
  for res in p.topology.residues():
   named={a.name:a.index for a in res.atoms()}
   for center,names in [('CA',['N','C','CB'])]+([('CB',['CA','CG1','CG2'])] if res.name=='ILE' else []):
    v=float(np.linalg.det([xyz[named[n]]-xyz[named[center]] for n in names]));assert v>1e-4,(res.id,center,v);volumes.append({'resid':res.id,'residue':res.name,'center':center,'signed_volume_nm3':v})
  assert any(res.name=='ALA' and res.id=='28' for res in p.topology.residues())
  assert not any(x['structural_break'] for x in preparation.inspect_preparation(result['id'])['gaps'])
  assert result['preparation']['stereochemistry']['violations']==[]
  r.update(status='pass',atoms=result['n_atoms'],rebuilt_segments=result['preparation']['rebuilt_segments'],validated_center_volumes=volumes,no_long_backbone_gap=True,uncertainty_warnings=result['warnings'],elapsed_seconds=j['elapsed_seconds'])
except Exception as exc:r.update(status='fail',error=repr(exc))
(ROOT/'docs/audit/preparation-positive-loop.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));sys.exit(r['status']!='pass')
