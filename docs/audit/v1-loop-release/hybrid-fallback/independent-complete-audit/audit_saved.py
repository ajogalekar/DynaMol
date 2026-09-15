"""Read-only independent saved-complex audit. No minimization or dynamics."""
import hashlib,json,sys,time,resource
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np
import openmm as mm
from openmm import app,unit
from pdbfixer import PDBFixer
ROOT=Path('/Users/ashujo/Documents/Science/DynaMol');sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'docs/audit/advanced-chemistry/covalent-v1-focused'))
from backend.loop_geometry import loop_geometry_report
from backend.preparation_worker import stereochemistry_report,require_valid_stereochemistry
from backend.prepared_system import load_prepared_forcefield
from backend.residue_identity import residue_key
from rama_reference import rama_report,dihedral64
RUN=Path('/Users/ashujo/.cache/dynamol-research/v1-loop-release/outer-torsion-anchors-v1')
OUT=RUN.parent/'outer-anchors-independent-v1';OUT.mkdir(exist_ok=True)
assert not (OUT/'audit-summary.json').exists(), 'Preserve completed audit output'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def key(a):return (*residue_key(a.residue),a.name)
def bonds(t):return {tuple(sorted([key(a),key(b)])) for a,b in t.bonds()}
def xyz(s):return np.asarray(s.positions.value_in_unit(unit.nanometer))
def parameter_comparison(first_xml,second_xml):
 # Keep exact indexed parameters/multiplicity; ignore only bond-term serialization order.
 first,second=ET.fromstring(first_xml),ET.fromstring(second_xml)
 forces1,forces2=first.find('Forces'),second.find('Forces')
 hb1=[f for f in forces1 if f.attrib.get('type')=='HarmonicBondForce']
 hb2=[f for f in forces2 if f.attrib.get('type')=='HarmonicBondForce']
 assert len(hb1)==len(hb2)==1
 terms=lambda f:Counter(tuple(sorted(term.attrib.items())) for term in f.find('Bonds'))
 termsequal=terms(hb1[0])==terms(hb2[0])
 assert hb1[0].attrib==hb2[0].attrib
 count=len(hb1[0].find('Bonds'))
 forces1.remove(hb1[0]);forces2.remove(hb2[0])
 canonical=lambda e:ET.canonicalize(ET.tostring(e,encoding='unicode'),strip_text=True)
 restequal=canonical(first)==canonical(second)
 return {'harmonic_bond_parameter_multisets_equal':termsequal,'bond_terms':count,
  'other_system_fields_and_terms_exactly_equal':restequal,
  'difference_is_only_harmonic_bond_term_order':first_xml!=second_xml and termsequal and restequal,
  'all_atom_indexed_parameters_equal':termsequal and restequal}
start=time.monotonic();plan=load(RUN/'plan.json');dataset=Path(plan['dataset']);parameters=load(dataset/'loop-refinement-parameters.json')
hashchecks={p:h for k in ['source_sha256','implementation_sha256'] for p,h in plan[k].items()}
for a in plan['attempts']:hashchecks[a['candidate']]=a['candidate_sha256']
assert all(sha(p)==h for p,h in hashchecks.items())
arch=app.PDBxFile(str(dataset/'loop-refinement-topology.cif'));top=arch.topology;initial=np.load(dataset/'loop-refinement-input.npz')['xyz_nm'];atoms=list(top.atoms());keys=[key(a) for a in atoms];indices={k:i for i,k in enumerate(keys)}
assert len(keys)==len(set(keys))==len(initial)
modeled=set(map(tuple,parameters['modeled_residue_keys']));entries={a['name']:a for a in plan['attempts']};summary=[]
templates=PDBFixer(filename=str(dataset/'prepared.pdb')).templates
ff,files=load_prepared_forcefield(dataset,parameters['preparation'],solvent=parameters['solvent'])
def fresh(t):return ff.createSystem(t,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
physical=fresh(top);fresh_xml=mm.XmlSerializer.serialize(physical);fresh_sha=hashlib.sha256(fresh_xml.encode()).hexdigest()
masses=[physical.getParticleMass(i).value_in_unit(unit.dalton) for i in range(physical.getNumParticles())]
physical_report={'fresh_native_system_sha256':fresh_sha,'particle_count':physical.getNumParticles(),'zero_mass_particles':[keys[i] for i,m in enumerate(masses) if m==0],
 'force_classes':[type(f).__name__ for f in physical.getForces()],'constraint_count':physical.getNumConstraints(),'parameter_files':files,
 'no_construction_force_classes':not any(type(f).__name__ in ['CustomExternalForce','CustomCompoundBondForce','CustomTorsionForce'] for f in physical.getForces())}
assert not physical_report['zero_mass_particles'] and physical_report['no_construction_force_classes']
assert physical.getNumParticles()==len(atoms)
for name in ['04-conditioned','08-shifted']:
 folder=RUN/name/'validation';result=load(folder/'result.json');recorded=load(folder/'reference.json');savedrecord=load(folder/'saved-reference.json')
 assert result['accepted_by_recorded_checks'] and all(result['checks'].values())
 assert set(result['checks'])=={'identity','retained_environment','geometry','stereochemistry','backbone_reference','observed_displacement','serialized_output','parameter_integrity'}
 verify={**result['source_sha256'],**result['artifacts_sha256']};assert all(sha(p)==h for p,h in verify.items());hashchecks.update(verify)
 hashchecks[str(folder/'result.json')]=sha(folder/'result.json')
 final=np.load(folder/'refined.npz')['xyz_nm'];saved=app.PDBFile(str(folder/'refined.pdb'));savedxyz=xyz(saved);sa=list(saved.topology.atoms())
 assert [(key(a),a.element.symbol) for a in sa]==[(key(a),a.element.symbol) for a in atoms]
 assert bonds(top)==bonds(saved.topology)
 assert final.shape==initial.shape==savedxyz.shape
 assert sha(folder/'native-system.xml')==fresh_sha==result['native_system_sha256']
 saved_fresh_xml=mm.XmlSerializer.serialize(fresh(saved.topology))
 saved_fresh_sha=hashlib.sha256(saved_fresh_xml.encode()).hexdigest()
 saved_parameters=parameter_comparison(fresh_xml,saved_fresh_xml)
 assert saved_parameters['all_atom_indexed_parameters_equal']
 source=app.PDBFile(result['scaffold_binding']['source']);sourcexyz=xyz(source);source_atoms={key(a):a for a in source.topology.atoms()}
 flanks=set(map(tuple,result['observed_context_residue_keys']));selected=modeled|flanks
 outside=[i for i,a in enumerate(atoms) if a.element!=app.element.hydrogen and residue_key(a.residue) not in selected]
 assert np.array_equal(initial[outside],final[outside])
 # Independent full-topology selection from each actual defining atom; do not call context_integrity.
 fullref=rama_report(top,final*unit.nanometer);expected=set(modeled)
 atom_by_res={r:{a.name:a.index for a in r.atoms()} for r in top.residues()}
 for row in fullref['rows']:
  before=atoms[row['phi_atoms'][0]].residue;defines=sorted(set(row['phi_atoms']+row['psi_atoms']+[atom_by_res[before]['CA']]))
  if any(residue_key(atoms[i].residue) in modeled for i in defines) or not np.array_equal(initial[defines],final[defines]):expected.add(tuple(row['residue']))
 assert expected=={tuple(r['residue']) for r in recorded['rows']}
 reference=rama_report(top,final*unit.nanometer,expected);serialized_reference=rama_report(saved.topology,saved.positions,expected)
 assert reference['all_selected_scored_without_outliers'] and serialized_reference['all_selected_scored_without_outliers']
 for freshref,oldref in [(reference,recorded),(serialized_reference,savedrecord)]:
  assert [(tuple(r['residue']),r['classification']) for r in freshref['rows']]==[(tuple(r['residue']),r['classification']) for r in oldref['rows']]
 # ARG57 exact phi, psi and preceding omega support across source/archive/refined/reloaded file.
 definitions={'phi':[('A','56','','ASN','C'),('A','57','','ARG','N'),('A','57','','ARG','CA'),('A','57','','ARG','C')],
  'psi':[('A','57','','ARG','N'),('A','57','','ARG','CA'),('A','57','','ARG','C'),('A','58','','THR','N')],
  'preceding_omega':[('A','56','','ASN','CA'),('A','56','','ASN','C'),('A','57','','ARG','N'),('A','57','','ARG','CA')]}
 boundary={}
 for label,ks in definitions.items():
  ix=[indices[k] for k in ks];src=np.array([sourcexyz[source_atoms[k].index] for k in ks])
  assert np.array_equal(src,initial[ix]) and np.array_equal(initial[ix],final[ix]) and np.array_equal(initial[ix],savedxyz[ix])
  boundary[label]={'atoms':ks,'all_four_coordinate_sets_exactly_equal':True,'degrees':dihedral64(src),'source_xyz_nm':src.tolist(),'maximum_change_A':0.0}
 arg=('A','57','','ARG');assert arg not in expected
 argref=rama_report(top,final*unit.nanometer,{arg});assert argref['outliers']
 fixed=[indices[tuple(k)] for k in result['additional_fixed_observed_atoms']]
 assert np.array_equal(initial[fixed],final[fixed]) and np.array_equal(initial[fixed],savedxyz[fixed])
 retained=[]
 for r in top.residues():
  if r.name not in {'ATP','TPO'}:continue
  rk=residue_key(r);ra=list(r.atoms());ridx=[a.index for a in ra];heavy=[a.index for a in ra if a.element!=app.element.hydrogen]
  archiveheavy={key(atoms[i]) for i in heavy};sourceheavy={k for k,a in source_atoms.items() if k[:4]==rk and a.element!=app.element.hydrogen}
  assert archiveheavy==sourceheavy
  assert np.array_equal(initial[heavy],final[heavy])
  sourcepoints=np.array([sourcexyz[source_atoms[key(atoms[i])].index] for i in heavy])
  # CIF and PDB parsers can differ by floating point roundoff after unit conversion.
  # Original PDB and saved PDB are bit-identical; archive and refined NPZ are too.
  assert np.array_equal(sourcepoints,savedxyz[heavy])
  crossformat_delta=float(np.max(np.linalg.norm(sourcepoints-initial[heavy],axis=1))*10)
  assert crossformat_delta<1e-12
  incident={b for b in bonds(top) if any(k[:4]==rk for k in b)};savedincident={b for b in bonds(saved.topology) if any(k[:4]==rk for k in b)}
  assert incident==savedincident
  hydrogen_parent_pairs=[]
  for a,b in top.bonds():
   for h,parent in [(a,b),(b,a)]:
    if h.residue is r and h.element==app.element.hydrogen:hydrogen_parent_pairs.append([key(h),key(parent)])
  retained.append({'identity':rk,'atoms':len(ra),'heavy_atoms':len(heavy),'hydrogens':len(ra)-len(heavy),
    'source_and_archive_heavy_identities_exact':True,'archive_refined_heavy_coordinates_exact':True,
    'source_pdb_saved_pdb_heavy_coordinates_exact':True,
    'source_pdb_archive_cif_maximum_numeric_roundoff_A':crossformat_delta,'incident_bonds':len(incident),
    'saved_bonds_equal_archive':True,'hydrogen_parent_pairs':hydrogen_parent_pairs,'archive_saved_atom_keys':[key(a) for a in ra],
    'maximum_hydrogen_coordinate_change_A':float(max([np.linalg.norm(final[i]-initial[i])*10 for i in ridx if atoms[i].element==app.element.hydrogen] or [0]))})
 assert {tuple(r['identity']) for r in retained}=={('A','170','','TPO'),('E','381','','ATP')}
 geometry=loop_geometry_report(top,final*unit.nanometer,selected);sg=loop_geometry_report(saved.topology,saved.positions,selected)
 stereo=stereochemistry_report(top,final*unit.nanometer,templates);ss=stereochemistry_report(saved.topology,saved.positions,templates)
 require_valid_stereochemistry(stereo,'Independent saved audit');require_valid_stereochemistry(ss,'Independent serialized audit')
 assert geometry['accepted'] and sg['accepted'] and not geometry['gross_collisions'] and not sg['gross_collisions']
 moving=[i for i,a in enumerate(atoms) if residue_key(a.residue) in flanks and a.element!=app.element.hydrogen]
 maxmove=float(np.max(np.linalg.norm(final[moving]-initial[moving],axis=1))*10);savedmove=float(np.max(np.linalg.norm(savedxyz[moving]-initial[moving],axis=1))*10)
 assert maxmove<=1 and savedmove<=1
 details={'geometry':geometry,'saved_geometry':sg,'stereo':stereo,'saved_stereo':ss,'reference':reference,'saved_reference':serialized_reference}
 sub=OUT/name;sub.mkdir()
 for label,data in details.items():save(sub/(label+'.json'),data)
 summary.append({'name':name,'result_sha256':sha(folder/'result.json'),'verified_artifact_sha256':result['artifacts_sha256'],
 'atom_count':len(atoms),'bond_count':len(bonds(top)),'all_atom_identity_elements_bonds_preserved':True,
 'outside_selected_observed_heavy_atoms_exactly_preserved':True,'selected_reference_keys':sorted(expected),'selected_reference_count':len(expected),
 'recomputed_geometry_stereo_reference_pass':True,'maximum_saved_rounding_difference_A':float(np.max(np.linalg.norm(savedxyz-final,axis=1))*10),
 'observed_heavy_maximum_displacement_A':maxmove,'saved_observed_heavy_maximum_displacement_A':savedmove,
 'ARG57_definitions':boundary,'ARG57_current_reference':argref['rows'][0],
 'ARG57_excluded_only_because_all_defining_coordinates_unchanged':True,'retained_components':retained,
 'fresh_system_matches_archive_and_saved_topology_parameters':True,
 'saved_topology_fresh_system_sha256':saved_fresh_sha,'saved_topology_parameter_comparison':saved_parameters,
 'source_hashes_verified':True,'artifact_files_unchanged':True,
 'scope':'Independent private static saved-artifact verification, not preparation rerun or MD/release qualification.'})
assert all(sha(p)==h for p,h in hashchecks.items())
report={'plan_sha256':sha(RUN/'plan.json'),'audit_implementation_sha256':sha(__file__),'verified_hashes':hashchecks,'physical_system':physical_report,
 'results':summary,'native_library_sha256':reference['native_library_sha256'],'elapsed_seconds':time.monotonic()-start,
 'max_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'all_requested_checks_pass':True,
 'app_ready':False,'physical_model_validated':False,'new_md_or_qm':False,
 'parameter_limit':'Same current frozen native force field with archived ligand/TPO manifests; historical protein/water distribution equivalence was not established.'}
save(OUT/'audit-summary.json',report);save(Path(__file__).with_name('audit-summary.json'),report)
print(json.dumps({'results':[{k:r[k] for k in ['name','atom_count','bond_count','selected_reference_count','observed_heavy_maximum_displacement_A','maximum_saved_rounding_difference_A']} for r in summary],'physical_system':physical_report,'elapsed_seconds':report['elapsed_seconds'],'max_rss_bytes':report['max_rss_bytes']},indent=2))
