"""Read-only independent comparison of tiny already-computed artifacts; no QM imports or launches."""
from pathlib import Path
import hashlib,json,math
import numpy as np
P=Path('/Users/ashujo/.cache/dynamol-research/qm-journal-native-qualification-v2')
O=Path(__file__).parent

def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def envelope(path):
 e=read(path);assert set(e)=={'payload','payload_sha256'} and digest(e['payload'])==e['payload_sha256'];return e['payload'],e['payload_sha256']
def arrays(path):
 with np.load(path,allow_pickle=False) as f:d={k:f[k].copy() for k in f.files}
 assert set(d)=={'coords_bohr','coords_angstrom','gradient_hartree_per_bohr','atomic_numbers','effective_nuclear_charges'}
 for k,v in d.items():assert v.shape==((3,) if k in ('atomic_numbers','effective_nuclear_charges') else (3,3)) and v.dtype.kind in 'if' and np.isfinite(v).all()
 assert np.array_equal(d['atomic_numbers'],[8,1,1]) and np.array_equal(d['effective_nuclear_charges'],[8,1,1])
 return d

assert sha(P/'package-manifest.json')=='aa9e15d38e9a45f60a4258e3b153f4058797b595beec99df7ac672f9c699176e'
assert sha(P/'artifact-manifest.json')=='0b9aa5bea3ae0a088a5f4f3bbea43332418d977944b4b66f7133d896356ebce8'
outer=read(P/'artifact-manifest.json')
for e in outer['files']:
 path=P/e['file'];assert path.is_file() and not path.is_symlink() and path.stat().st_size==e['bytes'] and sha(path)==e['sha256']
for f,h in read(P/'package-manifest.json')['files'].items():assert sha(P/f)==h
original='ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1'
candidate='780b5bb968fa9d45987597be9412331d071e17ea38b2299749dc5e92be9f12b5'
helper='b3a0543b5b2394373e7a787ad01ecd77865cc507f85bec8ec99efe6efed26528'
assert sha(P/'sources/original/qm_worker.py')==original and sha(P/'sources/candidate/qm_worker.py')==candidate
assert sha(P/'sources/candidate/qm_evaluation_journal.py')==helper
plan=read(P/'plan.json');campaign=read(P/'run/result.json');summary=read(P/'SUMMARY.json')
report={'scope':'Independent read-only verification of saved three-atom numerical qualification. No native calculations.',
 'outer_manifest_sha256':sha(P/'artifact-manifest.json'),'package_manifest_sha256':sha(P/'package-manifest.json'),
 'outer_artifacts_verified':len(outer['files']),'original_candidate_helper_sources_verified':True,'profiles':{},'QM_launched':False,'physical_acceptance':False}
for name in ['rhf','rks']:
 req=read(P/'requests'/f'{name}.json');assert req['atom_ids']==['water:O','water:H1','water:H2'] and req['charge']==0 and req['spin']==0
 assert 'initial_checkpoint' not in req and req['density_fit'] is True and req['basis']=='sto-3g'
 assert req['method']=={'rhf':'RHF','rks':'RKS'}[name]
 if name=='rks':assert req['functional']=='B3LYPG' and req['grid_level']==2
 pair=[]
 for role,whash in [('original',original),('candidate',candidate)]:
  folder=P/'run'/f'{name}-{role}'/'native';r=read(folder/'result.json');a=arrays(folder/'arrays.npz')
  assert r['accepted'] is True and r['scf']['converged'] is True and r['optimization']['converged'] is True
  assert r['input_sha256']==sha(P/'requests'/f'{name}.json') and r['worker_sha256']==whash
  for f,k in [('arrays.npz','arrays_sha256'),('scf.chk','checkpoint_sha256'),('resolved-method.json','resolved_method_file_sha256')]:assert sha(folder/f)==r[k]
  assert (folder/'constraints.txt').read_bytes()==b'$freeze\nxyz 1\n'
  assert r['actual_pyscf_threads']==1 and r['nao']==7 and r['explicit_electron_count']==10
  assert r['method']['method']==req['method'] and r['method']['density_fit'] is True and r['method']['basis_requested']=='sto-3g'
  for k,v in req['optimization'].items():assert r['optimization']['settings'][k]==v
  assert np.linalg.norm(a['coords_bohr'][0])<1e-5
  pair.append((r,a))
 lr,la=pair[0];rr,ra=pair[1]
 assert lr['energy_hartree']==rr['energy_hartree'] and lr['method']==rr['method']
 assert all(np.array_equal(la[k],ra[k]) for k in la)
 assert len(lr['optimization']['steps'])==len(rr['optimization']['steps'])==5
 for l,r in zip(lr['optimization']['steps'],rr['optimization']['steps']):
  assert l['cycle']==r['cycle'] and l['energy_hartree']==r['energy_hartree'] and l['gradient_norm_hartree_per_bohr']==r['gradient_norm_hartree_per_bohr']
 folder=P/'run'/f'{name}-candidate'/'native';manifest,mhash=envelope(folder/'evaluation-journal/manifest.json')
 assert manifest['source_input_file_sha256']==sha(P/'requests'/f'{name}.json') and manifest['worker_source_file_sha256']==candidate
 assert manifest['callback_contract']['helper_source_sha256']==helper
 assert manifest['requested_method_fingerprint_sha256']==digest(manifest['requested_method_and_settings'])
 assert manifest['callback_contract']['native_initial_method_sha256']==digest(manifest['callback_contract']['native_initial_method'])
 identity=manifest['identity'];assert identity=={'atom_ids':req['atom_ids'],'elements':req['elements'],'charge':0,'spin':0}
 rows=[];geometries=[]
 for step in rr['optimization']['steps']:
  cycle=step['cycle'];record_path=folder/f'evaluation-journal/evaluation-{cycle:06d}.json';row,_=envelope(record_path)
  assert row['identity']==identity and row['manifest_sha256']==mhash and row['scf_converged'] is True
  assert row['requested_method_fingerprint_sha256']==manifest['requested_method_fingerprint_sha256']
  for k in ['accepted','optimization_converged','physical_acceptance','checkpoint_reuse_authorized']:assert row[k] is False and manifest[k] is False
  assert row['status']==manifest['status']=='UNCONVERGED' and row['scf_checkpoint_same_geometry_validity']==manifest['scf_checkpoint_same_geometry_validity']=='NOT_ESTABLISHED'
  assert row['geometry_sha256']==digest({'identity':identity,'unit':'bohr','coords':row['coords_bohr']})
  assert row['gradient_sha256']==digest({'unit':'hartree/bohr','gradient':row['gradient_hartree_per_bohr']})
  assert row['units']=={'coordinates':'bohr','energy':'hartree','gradient':'hartree/bohr'}
  assert row['energy_hartree']==step['energy_hartree'] and np.linalg.norm(row['gradient_hartree_per_bohr'])==step['gradient_norm_hartree_per_bohr']
  assert step['evaluation_record']['file_sha256']==sha(record_path) and step['evaluation_record']['accepted'] is False
  replay=P/'run'/f'{name}-replay-{cycle:06d}'/'native';request=P/'run'/f'{name}-replay-{cycle:06d}-input.json'
  r=read(replay/'result.json');a=arrays(replay/'arrays.npz');q=read(request)
  assert q['coords_bohr']==row['coords_bohr'] and q['optimization']=={'enabled':False} and 'initial_checkpoint' not in q
  assert r['input_sha256']==sha(request) and r['worker_sha256']==original and sha(replay/'arrays.npz')==r['arrays_sha256']
  assert r['accepted'] is True and r['scf']['converged'] is True and r['method']==rr['method']
  assert np.array_equal(a['coords_bohr'],row['coords_bohr'])
  de=abs(r['energy_hartree']-row['energy_hartree']);dg=float(np.max(abs(a['gradient_hartree_per_bohr']-row['gradient_hartree_per_bohr'])))
  assert de<=1e-8 and dg<=1e-7
  rows.append({'cycle':cycle,'energy_difference_hartree':de,'gradient_max_difference_hartree_per_bohr':dg,'exact_geometry':True})
  geometries.append(row['geometry_sha256'])
 assert len(set(geometries))==5
 report['profiles'][name]={'pair_all_arrays_and_final_energy_exact':True,'moving_evaluations':5,'replays':rows}
f=P/'run/induced-persistence-failure/native';r=read(f/'result.json');inject=read(f/'failure-injection.json')
assert r['accepted'] is False and r['status']=='failed' and r['error']['type']=='OSError' and 'Injected qualification failure' in r['error']['message']
assert r['optimization']['steps']==[] and not (f/'arrays.npz').exists() and len(inject['injected'])==1 and inject['injected'][0]['errno']==28
assert not list((f/'evaluation-journal').glob('evaluation-*.json'))
tmp=list((f/'evaluation-journal').glob('.evaluation-000001.json.*.tmp'));assert len(tmp)==1
payload,_=envelope(tmp[0]);assert payload['scf_converged'] is True and payload['accepted'] is False
report['ENOSPC_before_commit_rejects_acceptance_and_retains_temporary']=True
children=campaign['children'];assert len(children)==15 and all(c['finished_unix']-c['started_unix']<1 for c in children)
report['memory_evidence']='All children <1s; once-per-second samples observe startup, not true peak RSS. No large-system memory qualification.'
report['maximum_native_child_elapsed_seconds']=max(c['finished_unix']-c['started_unix'] for c in children)
report['status']='no_blocker_for_scoped_three_atom_journal_candidate'
assert sha(P/'artifact-manifest.json')==report['outer_manifest_sha256']
(O/'review.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='profiles'},indent=2))
