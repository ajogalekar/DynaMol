"""Independent geometric witnesses for new-branch stereo repair, no dynamics."""
from pathlib import Path
import copy,hashlib,json,sys,time
import numpy as np
from openmm import app,unit
from pdbfixer import PDBFixer
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from backend.added_atom_stereo import correct_added_branches
from backend.preparation_worker import atom_key,stereochemistry_report
OUT=Path(__file__).resolve().parent
fixture=ROOT/'docs/audit/preparation-fixtures/six_residues_intact.pdb'
templates=PDBFixer(filename=str(fixture)).templates
results=[]
def reflect(xyz,moving,center,left,right):
 n=np.cross(xyz[left]-xyz[center],xyz[right]-xyz[center]);n/=np.linalg.norm(n)
 result=xyz.copy();d=result[moving]-xyz[center];result[moving]-=2*np.outer(d@n,n);return result
for name in ['THR','ILE']:
 template=templates[name];top=copy.deepcopy(template.topology);xyz=np.asarray(template.positions.value_in_unit(unit.nanometer)).copy();ids={a.name:a.index for a in top.atoms()}
 members=[a.index for a in top.atoms() if a.name not in {'N','CA','C','O','OXT'}]
 wrong=reflect(xyz,members,ids['CA'],ids['N'],ids['C'])
 observed={atom_key(a) for a in top.atoms() if a.index not in members}
 fixed,report=correct_added_branches(top,wrong*unit.nanometer,templates,observed);actual=fixed.value_in_unit(unit.nanometer)
 assert not report['stereochemistry_after']['violations']
 for a,b in top.bonds():assert abs(np.linalg.norm(wrong[a.index]-wrong[b.index])-np.linalg.norm(actual[a.index]-actual[b.index]))<1e-12
 fixedids=[a.index for a in top.atoms() if atom_key(a) in observed];assert np.array_equal(actual[fixedids],wrong[fixedids])
 results.append({'case':name+' full new sidechain CA inversion','passed':True,'corrections':report['corrections'],'all_bond_lengths_unchanged':True,'observed_atoms_unchanged':True})
# An independently placed observed environment atom can occupy the corrected
# branch position. Stereo validity cannot be used as a steric safety check.
template=templates['THR'];top=copy.deepcopy(template.topology);xyz=np.asarray(template.positions.value_in_unit(unit.nanometer)).copy();ids={a.name:a.index for a in top.atoms()}
wrong=reflect(xyz,[ids['CG2']],ids['CB'],ids['CA'],ids['OG1'])
other=top.addResidue('LIG',top.addChain('Z'),'99');external=top.addAtom('C1',app.element.carbon,other)
wrong=np.vstack([wrong,xyz[ids['CG2']]])
observed={atom_key(a) for a in top.atoms() if a.index!=ids['CG2']}
fixed,report=correct_added_branches(top,wrong*unit.nanometer,templates,observed);actual=fixed.value_in_unit(unit.nanometer)
before=float(np.linalg.norm(wrong[ids['CG2']]-wrong[external.index])*10);after=float(np.linalg.norm(actual[ids['CG2']]-actual[external.index])*10)
assert before>1.3 and after<1e-10 and not report['stereochemistry_after']['violations']
assert np.array_equal(actual[external.index],wrong[external.index])
results.append({'case':'observed environment atom occupies corrected methyl position','witness_reproduced':True,'before_contact_angstrom':before,'after_contact_angstrom':after,'stereo_checks_pass':True,'observed_atom_unchanged':True,'conclusion':'A final complete-environment gross-contact gate for corrected atoms is needed; reflection and stereochemistry alone do not ensure acceptable contacts.'})
report={'scope':'Independent synthetic geometry review of new-branch correction. No native-loop accuracy or production stability claim; no dynamics.','helper_sha256':hashlib.sha256((ROOT/'backend/added_atom_stereo.py').read_bytes()).hexdigest(),'fixture_sha256':hashlib.sha256(fixture.read_bytes()).hexdigest(),'cases':results}
(OUT/'benchmark-added-branch-review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
