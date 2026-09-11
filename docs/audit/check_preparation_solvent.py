"""Small direct solvent checks. Input is a modeled 7-residue software fixture."""
from pathlib import Path
import hashlib,json,shutil,sys,tempfile
import numpy as np
from openmm import app,unit
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from backend import config,sources,solvent,storage
results=[]
def check(name,f):
 try:results.append({'name':name,'status':'pass','evidence':f()})
 except Exception as e:results.append({'name':name,'status':'fail','error':repr(e)})
with tempfile.TemporaryDirectory(prefix='dynamol-solvent-review-') as tmp:
 config.DATA_ROOT=Path(tmp);config.DATASETS_DIR=Path(tmp)/'datasets';config.DATASETS_DIR.mkdir()
 inp=ROOT/'docs/audit/preparation-gromacs-probes/histidine-acid-plain/input.pdb'
 m=sources.import_structure(inp,name='Synthetic pH2 HIS fragment')
 folder=storage.dataset_dir(m['id']);shutil.copy2(inp,folder/'prepared.pdb')
 m['preparation']={'ph':2.0,'method':'Independent upstream Modeller probe fixture, not DynaMol worker','prepared_pdb':'prepared.pdb'}
 storage.atomic_json(folder/'metadata.json',m)
 before=app.PDBFile(str(inp));bx=np.asarray(before.positions.value_in_unit(unit.nanometer));n=len(bx)
 hold={}
 def preview():
  s=solvent.solvate_dataset(m['id'],padding_nm=1,seed=2026,ph=2);hold['preview']=s
  physical=storage.load_physical(s['id']);sp=app.PDBFile(str(storage.dataset_dir(s['id'])/'prepared.pdb'));sx=np.asarray(sp.positions.value_in_unit(unit.nanometer))
  assert np.array_equal(physical.xyz[0,:n],bx.astype(np.float32))
  assert np.max(np.abs(sx-physical.xyz[0]))<=5.1e-5
  assert np.array_equal(sx[:n],bx)
  assert s['has_unitcell'] and s['n_atoms']>n
  old_ids=[(a.name,a.residue.name,a.residue.id,a.element.symbol) for a in before.topology.atoms()]
  new_ids=[(a.name,a.residue.name,a.residue.id,a.element.symbol) for a in list(sp.topology.atoms())[:n]]
  assert old_ids==new_ids
  assert any(a.name=='HE2' and a.residue.id=='68' for a in sp.topology.atoms())
  return {'solute_atoms':n,'total_atoms':s['n_atoms'],'solute_prefix_exact':True,'ph':s['preparation']['ph'],'histidine_HE2_retained':True,'max_pdb_roundtrip_error_nm':float(np.max(np.abs(sx-physical.xyz[0]))),'box_vectors_nm':physical.unitcell_vectors.tolist(),'solvation':s['solvation']}
 check('Real TIP3P preview retains prepared atom order, hydrogen state, and physical solute coordinates',preview)
 def reuse():
  s=hold['preview'];again=solvent.solvate_dataset(s['id'],padding_nm=1,seed=2026,ph=2)
  assert again['id']==s['id'] and again['n_atoms']==s['n_atoms']
  return {'same_dataset_id':True,'atoms_unchanged':True}
 check('Identical solvent-preview settings return the same dataset without double solvation',reuse)
 def rebuild():
  s=hold['preview'];changed=solvent.solvate_dataset(s['id'],padding_nm=1.1,seed=2026,ph=7)
  assert changed['solvation']['parent_dataset_id']==m['id']
  assert changed['preparation']['ph']==2.0
  assert any('pH 2' in w for w in changed['warnings'])
  physical=storage.load_physical(changed['id']);assert np.array_equal(physical.xyz[0,:n],bx.astype(np.float32))
  return {'new_total_atoms':changed['n_atoms'],'rebuilt_from_original_prepared_parent':True,'original_ph_retained':True,'ph_mismatch_warned':True}
 check('Changed padding rebuilds from original prepared solute and warns if requested pH differs',rebuild)
report={'scope':'Direct app solvent helper with synthetic upstream-prepared tiny peptide. No dynamics or equilibration; temporary generated solvent outputs.','input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'results':results}
(ROOT/'docs/audit/preparation-solvent-checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));sys.exit(any(x['status']=='fail' for x in results))
