"""Analytic geometry checks for the discrete sidechain heuristic, not pose quality."""
from pathlib import Path
import json,sys
import numpy as np
from openmm import app,unit
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from backend.preparation_worker import adjust_sidechains,torsion_degrees,proper_overlay_points,stereochemistry_report,require_valid_stereochemistry
results=[]
def check(name,f):
 try:results.append({'name':name,'status':'pass','evidence':f()})
 except Exception as exc:results.append({'name':name,'status':'fail','error':repr(exc)})
def topology(name='SER',link=False):
 t=app.Topology();c=t.addChain('A');r=t.addResidue(name,c,'1');atoms=[]
 defs=[('N',app.element.nitrogen),('CA',app.element.carbon),('C',app.element.carbon),('O',app.element.oxygen),('CB',app.element.carbon),('OG' if name=='SER' else 'SG',app.element.oxygen if name=='SER' else app.element.sulfur)]
 for n,e in defs:atoms.append(t.addAtom(n,e,r))
 for a,b in [(0,1),(1,2),(2,3),(1,4),(4,5)]:t.addBond(atoms[a],atoms[b])
 r2=t.addResidue('LIG',c,'2');other=t.addAtom('C1',app.element.carbon,r2)
 if link:t.addBond(atoms[5],other)
 x=np.asarray([[-.10,0,0],[0,0,0],[.10,-.07,0],[.20,-.10,0],[.08,.08,.08],[.18,.14,.10],[.1805,.1405,.1005]])
 return t,x

def lengths(t,x):return np.asarray([np.linalg.norm(x[a.index]-x[b.index]) for a,b in t.bonds()])
def chirality(x):return float(np.linalg.det(np.asarray([x[0]-x[1],x[2]-x[1],x[4]-x[1]])))
def rigid_stereo():
 t,x=topology();after,report=adjust_sidechains(t,x*unit.nanometer);y=np.asarray(after.value_in_unit(unit.nanometer));r=report['adjustments']
 assert r and all(v['clash_score_after']<v['clash_score_before'] for v in r)
 assert np.max(np.abs(lengths(t,x)-lengths(t,y)))<1e-12
 assert np.array_equal(x[:4],y[:4]) and np.allclose(x[4],y[4],atol=1e-12,rtol=0)
 assert chirality(x)==chirality(y) and chirality(x)!=0
 assert np.array_equal(x[6],y[6])
 for record in r:
  assert abs(((torsion_degrees(y[[0,1,4,5]])-record['selected_degrees']+180)%360)-180)<1e-8
 return {'adjustment_count':len(r),'score_before':r[0]['clash_score_before'],'score_after':r[0]['clash_score_after'],'max_bond_length_change_nm':float(np.max(np.abs(lengths(t,x)-lengths(t,y)))),'backbone_exactly_unchanged':True,'beta_carbon_unchanged_within_1e-12_nm':True,'CA_signed_volume_preserved':chirality(y),'target_chi_reached':True,'note':'Synthetic steric collision demonstrates numerical invariants, not biological geometry accuracy.'}
check('Clash-lowering rotation preserves covalent lengths, fixed backbone, and CA handedness',rigid_stereo)
def crosslink():
 t,x=topology('CYS',link=True);after,report=adjust_sidechains(t,x*unit.nanometer);y=np.asarray(after.value_in_unit(unit.nanometer))
 assert np.array_equal(x,y) and report['adjusted_chi_count']==0
 assert any('crosslink' in s for s in report['skipped'])
 return {'all_coordinates_unchanged':True,'skipped':report['skipped'],'note':'Generic covalent crosslink fixture, not a real cysteine molecular model.'}
check('Distal covalent crosslink prevents rotating a sidechain into another residue',crosslink)
def kabsch_guard():
 moving=np.asarray([[0.,0,0],[1.,0,0],[0,2.,0],[0,0,3.]])
 R=np.asarray([[0.,-1,0],[1,0,0],[0,0,1.]])
 reference=moving@R.T+np.asarray([3.,4.,5.])
 shift,rotation,center=proper_overlay_points(reference,moving)
 assert np.linalg.det(rotation)>0.999999 and np.allclose((moving+shift)@rotation.T+center,reference,atol=1e-12,rtol=0)
 mirrored=moving*np.asarray([-1.,1.,1.]);shift,rotation,center=proper_overlay_points(mirrored,moving)
 assert np.linalg.det(rotation)>0.999999
 return {'known_proper_rotation_and_translation_recovered':True,'mirror_target_still_gets_det_plus_one_rotation':True}
check('Template alignment recovers rigid motion while forbidding reflections',kabsch_guard)
def reject_original_bad_loop():
 from pdbfixer import PDBFixer
 p=app.PDBFile(str(ROOT/'docs/audit/preparation-worker-evidence/dfafc6cf4b0649aa/prepared.pdb'))
 fixer=PDBFixer(filename=str(ROOT/'docs/audit/preparation-fixtures/six_residues_known_gap.pdb'))
 r=stereochemistry_report(p.topology,p.positions,fixer.templates)
 assert {(x['resid'],x['center']) for x in r['violations']}=={('3','CA'),('3','CB')}
 try:require_valid_stereochemistry(r,'Regression fixture')
 except ValueError as e:return {'violations':r['violations'],'rejection':str(e)}
 raise AssertionError('The historical inverted loop was accepted')
check('Historical inverted ILE alpha and beta centers are both detected and rejected',reject_original_bad_loop)
report={'scope':'Pure geometry/graph fixtures. No empirical rotamer or native-pose validation.','results':results};(ROOT/'docs/audit/preparation-geometry-checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));sys.exit(any(x['status']=='fail' for x in results))
