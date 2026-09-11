"""Independent boundary checks; synthetic fixtures, no scientific accuracy claim."""
from pathlib import Path
import hashlib,json,sys,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from backend import config,sources,preparation,storage
results=[]
def check(name,f):
 try:results.append({'name':name,'status':'pass','evidence':f()})
 except Exception as e:results.append({'name':name,'status':'fail','error':repr(e)})
def reject(f):
 try:f()
 except ValueError as e:return str(e)
 raise AssertionError('Invalid operation accepted')
with tempfile.TemporaryDirectory(prefix='dynamol-preparation-review-') as temp:
 config.DATA_ROOT=Path(temp);config.DATASETS_DIR=Path(temp)/'datasets';config.JOBS_DIR=Path(temp)/'jobs'
 config.DATASETS_DIR.mkdir();config.JOBS_DIR.mkdir()
 fixtures=ROOT/'docs/audit/preparation-fixtures'
 metadata={}
 def seqres():
  out={}
  for name in ['six_residues_intact','six_residues_missing_atom','six_residues_known_gap','six_residues_numbering_gap_only','histidine_fragment']:
   src=fixtures/(name+'.pdb');m=sources.import_structure(src);metadata[name]=m;i=preparation.inspect_preparation(m['id'])
   assert (storage.dataset_dir(m['id'])/m['source_file']).read_bytes()==src.read_bytes()
   if name=='six_residues_known_gap':assert i['has_sequence'] and any(e['residues']==['ILE'] for e in i['missing_residues'])
   if name=='six_residues_numbering_gap_only':assert not i['has_sequence'] and not i['missing_residues'] and any(g['structural_break'] for g in i['gaps'])
   if name=='six_residues_missing_atom':assert any('NZ' in a['atoms'] for a in i['missing_atoms'])
   out[name]={'has_sequence':i['has_sequence'],'missing_residues':i['missing_residues'],'missing_atoms':i['missing_atoms'],'gaps':i['gaps'],'source_hash':hashlib.sha256(src.read_bytes()).hexdigest()}
  return out
 check('Original SEQRES preserved; known sequence gap distinguished from numbering-only gap',seqres)
 def cif_chain_mapping():
  m=sources.import_structure(fixtures/'known_gap_multichar_chain.cif');i=preparation.inspect_preparation(m['id'])
  assert i['has_sequence'] and any(e['residues']==['ILE'] for e in i['missing_residues']),i
  return {'display_chains':sorted({a['chain'] for a in m['atoms']}),'has_sequence':i['has_sequence'],'missing_residues':i['missing_residues']}
 check('mmCIF original polymer sequence maps after multicharacter chain normalization',cif_chain_mapping)
 def cif_author_label_mapping():
  m=sources.import_structure(fixtures/'known_gap_label_auth_chain.cif');i=preparation.inspect_preparation(m['id'])
  assert i['has_sequence'] and any(e['residues']==['ILE'] for e in i['missing_residues']),i
  return {'display_chains':sorted({a['chain'] for a in m['atoms']}),'has_sequence':i['has_sequence'],'missing_residues':i['missing_residues']}
 check('mmCIF label/auth chain distinction preserves known missing sequence',cif_author_label_mapping)
 def gaps():
  original=preparation._submit;preparation._submit=lambda s,op:{'accepted':True,'operation':op}
  try:
   known=metadata['six_residues_known_gap']['id'];unknown=metadata['six_residues_numbering_gap_only']['id']
   messages=[reject(lambda:preparation.submit_preparation({'dataset_id':known})),reject(lambda:preparation.submit_preparation({'dataset_id':unknown,'build_missing_residues':True}))]
   accepted=preparation.submit_preparation({'dataset_id':known,'build_missing_residues':True})
   assert accepted['accepted']
   return {'rejections':messages,'known_short_loop_optin_accepted':True,'worker_launch':'mocked; no preparation performed by this check'}
  finally:preparation._submit=original
 check('Unresolved gaps fail before worker launch; known short-loop building requires opt-in',gaps)
 def smiles():
  smi='N[C@@H](C)C(=O)O';a=sources.create_smiles(smi,seed=2026);b=sources.create_smiles(smi,seed=2026)
  pa=storage.load_physical(a['id']);pb=storage.load_physical(b['id'])
  assert np.array_equal(pa.xyz,pb.xyz)
  c=a['chemistry'];assert '@' in c['canonical_isomeric_smiles'] and c['unassigned_stereocenters']==0
  unknown=sources.create_smiles('NC(C)C(=O)O');assert unknown['chemistry']['unassigned_stereocenters']==1
  return {'canonical_isomeric_smiles':c['canonical_isomeric_smiles'],'identical_seed_coordinates':True,'atoms':a['n_atoms'],'stereo_warning':unknown['warnings'],'salt_rejection':reject(lambda:sources.create_smiles('CC(=O)[O-].[Na+]'))}
 check('SMILES defined stereochemistry, unspecified-center warning, seeded reproducibility',smiles)
 def mol2_identity():
  # Interleaved residue atoms are legal enough to reach a parser; may be rejected
  # explicitly but must never silently scramble coordinate/atom identity.
  text='''@<TRIPOS>MOLECULE
fixture
4 2 2 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 0 0 0 C.3 1 LIG1
2 O1 10 0 0 O.3 2 LIG2
3 C2 1 0 0 C.3 1 LIG1
4 O2 11 0 0 O.3 2 LIG2
@<TRIPOS>BOND
1 1 3 1
2 2 4 1
'''
  p=Path(temp)/'interleaved.mol2';p.write_text(text)
  try:m=sources.import_structure(p)
  except ValueError as e:return {'correctly_rejected':str(e)}
  t=storage.load_physical(m['id']);expected={'C1':0,'O1':10,'C2':1,'O2':11}
  pairs=[(a.name,float(t.xyz[0,a.index,0]*10)) for a in t.topology.atoms]
  assert all(abs(x-expected[name])<1e-4 for name,x in pairs),pairs
  return {'atom_positions_angstrom':pairs}
 check('MOL2 interleaved residue atom-coordinate association preserved or rejected',mol2_identity)
 def pdb_cif_equivalence():
  pdb_path=ROOT/'data/sources/1UBQ.pdb';cif_path=ROOT/'data/datasets/28d695aa1b804e44/source.cif'
  assert cif_path.is_file(), 'Retained 1UBQ mmCIF source from source smoke test unavailable'
  a=sources.import_structure(pdb_path);b=sources.import_structure(cif_path)
  pa=storage.load_physical(a['id']);pb=storage.load_physical(b['id'])
  assert pa.n_atoms==pb.n_atoms==660 and np.array_equal(pa.xyz,pb.xyz)
  assert np.array_equal(pa.unitcell_lengths,pb.unitcell_lengths)
  assert [a.name for a in pa.topology.atoms]==[a.name for a in pb.topology.atoms]
  return {'atoms':660,'coordinate_arrays_identical':True,'unitcell_lengths_nm':pa.unitcell_lengths.tolist(),'pdb_sha256':hashlib.sha256(pdb_path.read_bytes()).hexdigest(),'cif_sha256':hashlib.sha256(cif_path.read_bytes()).hexdigest(),'source_cif':str(cif_path.relative_to(ROOT))}
 check('Retained RCSB 1UBQ PDB/mmCIF inputs agree in atoms, coordinates, and crystal box',pdb_cif_equivalence)
 def remote_reject():
  return {'pdb_url':reject(lambda:sources.fetch_structure('pdb','https://example.com/x')),'pubchem_url':reject(lambda:sources.fetch_structure('pubchem','https://example.com/x')),'arbitrary_provider':reject(lambda:sources.fetch_structure('other','1'))}
 check('Remote source identifiers cannot become arbitrary URLs',remote_reject)
report={'scope':'Independent import/inspection boundary checks; temporary data; no MD or loop reconstruction in this script.','results':results}
(ROOT/'docs/audit/preparation-source-checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
sys.exit(any(r['status']=='fail' for r in results))
