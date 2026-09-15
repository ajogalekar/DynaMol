"""Freeze source-evidence panels; never selects entries by preparation outcome."""
from pathlib import Path
import datetime, hashlib, json
from inventory_sources import METALS

OUT = Path(__file__).resolve().parent
COVALENT = [
 ('6OIM','KRAS–sotorasib','MOV','Cys thioether; acrylamide Michael product','GDP/Mg and its coordinating waters are retained; one three-residue internal loop.'),
 ('6UT0','KRAS–adagrasib','M1X','Cys thioether; fluoroacrylamide Michael product','Assembly1 selects deposited monomer A with its GDP/Mg; other independent monomers remain in the original input.'),
 ('5P9J','BTK–ibrutinib','8E8','Cys thioether; acrylamide Michael product','Three-residue loop and alternate conformers must be handled explicitly.'),
 ('6J6M','BTK–zanubrutinib','BA0','Cys thioether; acrylamide Michael product','Complete observed kinase backbone; retain IMD crystallization component.'),
 ('8A4V','Cathepsin L–E64','E64','Cys thioether; epoxide-opening product','Deposited E64 CCD already represents the open-chain product; keep its guanidinium/carboxyl state separate from pH-model assumptions.'),
 ('5TDI','Cathepsin K–odanacatib','7AS','Cys thioimidate; nitrile-addition product','CCD specifies iminomethyl bound form and stereochemistry; complete monomer.'),
 ('6QBS','Cathepsin K–alkyne inhibitor','HVW','Cys vinyl thioether; alkyne-addition product','Retain the declared dimer, Ca and Cl; the second chain uses author Cys240, not Cys25.'),
 ('3LJ7','Humanized-rat FAAH–URB597','OHO','Ser carbamoyl adduct','Carbamoyl fragment is deposited; leaving group is absent. Author-defined monomer assembly1; source separately proposes a software dimer.'),
 ('3LJ6','Humanized-rat FAAH–PF-3845','PIX','Ser carbamoyl adduct from urea inhibitor','Retain only the deposited reacted fragment, with chemically declared leaving atoms; no free inhibitor reconstruction.'),
 ('6AX1','Monoacylglycerol lipase–carbamate','C0S','Ser carbamoyl adduct','C0S CCD contains the free leaving group and labels its atoms as leaving; observed product omits them. Assembly1 selects the complete-backbone monomer.'),
 ('2H5I','Caspase-3–Ac-DEVD-CHO','ASJ','Cys thiohemiacetal; peptide-aldehyde product','Inhibitor is a polymeric peptide chain with terminal ASJ and ACE. Generate both assembly operators and retain both caspase subunits.'),
 ('3BJM','DPP4–saxagliptin','BJM','Ser imidate; nitrile-addition product','Full deposited dimer and all NAG retained. Source annotation of glycan connections requires validation, not silent glycan removal.'),
 ('5O4A','FGFR1–m-FSBA','9K8','Lys sulfonamide; sulfonyl-fluoride substitution','CCD F is explicitly a leaving atom; deposited adduct lacks it. Assembly1 has an eleven-residue internal loop.'),
 ('5LF3','Human 20S proteasome–bortezomib','BO2','N-terminal Thr boronate adduct','Combined-chemistry challenge: full 28-chain assembly, six inhibitor links, other modified Cys residues, ions and PEG; no reduction to an isolated active-site fragment.'),
 ('5F19','COX2–aspirin acetylation','OAS','Ser O-acetyl ester encoded as one modified residue','Combined-chemistry challenge: full dimer, glycans, cobalt protoporphyrin, OAS, and other deposited components retained. This tests product-as-modified-residue representation.'),
]
METAL = [
 ('1MNC','MMP8–hydroxamate','Catalytic Zn with three His and bidentate inhibitor; separate structural Zn and Ca.'),
 ('1HFC','MMP1–hydroxamate','Two Zn sites plus Ca; PLH hydroxamate bound-state protonation must be declared.'),
 ('1GKC','MMP9–hydroxamate','Catalytic and structural Zn plus Ca sites; all assembly1 metal sites retained.'),
 ('1XUC','MMP13–non-zinc-binding inhibitor','Separates protein-bound Zn/water sites from an inhibitor that does not chelate catalytic Zn.'),
 ('1BKC','ADAM17–hydroxamate','Zn metalloprotease from a different catalytic family; bound hydroxamate state requires explicit choice.'),
 ('1CA2','Carbonic anhydrase II','Zn with three His and deposited solvent oxygen. Water versus hydroxide cannot be concluded from X-ray oxygen positions alone.'),
 ('1T64','HDAC8–trichostatin A','Catalytic Zn plus deposited Ca/Na; retain each site with its own model assignment.'),
 ('1CLL','Calmodulin','Four Ca sites with carboxylate/backbone/solvent donors.'),
 ('1HET','Horse liver alcohol dehydrogenase','Catalytic and structural Zn with sulfur donors; retain NAD and MRD. Source describes a hydroxide adduct to NADH, requiring cofactor-state review.'),
 ('1A6Q','PPM1A phosphatase','Two Mn sites and phosphate; four-residue internal gap.'),
 ('1HL5','Cu/Zn superoxide dismutase','Assembly1 is the deposited A/H dimer with Cu/Zn sites. CCD declares Cu(II), but redox and bridging-His protonation are separate model decisions.'),
 ('1HCK','CDK2–ATP/Mg','ATP-associated Mg and coordinating waters; four-residue missing loop. Nucleotide formal charge must be modeled explicitly.'),
 ('1A6M','Oxy-myoglobin','Heme Fe and bound O2: redox, spin, Fe–O2 bonding and oxygen occupancy require a heme-specific model. CCD does not declare the Fe formal charge.'),
 ('1LUV','Human MnSOD H30V','Two deposited Mn sites expanded to the declared tetramer; CCD Mn(II) is recorded without claiming the experiment determines oxidation/spin state.'),
 ('1OKL','Carbonic anhydrase II–dansylamide','Official MCPB reference structure with Zn and a distinct Hg site; Hg is retained as a separate explicit model requirement, never silently removed.'),
]

def read_json(path):
 return json.loads(path.read_text())

def source_ccd(code):
 data=read_json(OUT/'inventories'/'ccd'/f'{code}.json')
 return {'id':code,'source':data['input'],'name':data['chem_comp'][0]['name'],
         'component_formal_charge':data['chem_comp'][0].get('pdbx_formal_charge'),
         'metal_atom_formal_charges':{x['atom_id']:x.get('charge') for x in data['atoms'] if x['type_symbol'].upper() in METALS},
         'leaving_atom_names':[x['atom_id'] for x in data['atoms'] if x.get('pdbx_leaving_atom_flag')=='Y']}

def belongs(endpoint, selected):
 return endpoint['label_chain'] in selected

def metal_match(atom, endpoint):
 return (atom['label_chain'],atom['author_chain'],atom['author_residue_id'],atom['component'],atom['atom']) == (
     endpoint['label_chain'],endpoint['author_chain'],endpoint['author_residue_id'],endpoint['component'],endpoint['atom'])

def base_case(code,label,notes):
 d=read_json(OUT/'inventories'/f'{code}.json')
 generators=[x for x in d['assembly_generators'] if x['assembly_id']=='1']
 selected=set(','.join(x['asym_id_list'] for x in generators).split(','))
 retained=[x for x in d['nonstandard_residues'] if x['label_chain'] in selected]
 ccds={x['component']:source_ccd(x['component']) for x in retained}
 deficits=[]
 for residue in retained:
  ccd=read_json(OUT/'inventories'/'ccd'/f"{residue['component']}.json")
  heavy={a['atom_id'] for a in ccd['atoms'] if a['type_symbol'] not in ('H','D')}
  leaving={a['atom_id'] for a in ccd['atoms'] if a.get('pdbx_leaving_atom_flag')=='Y'}
  observed=set(residue['observed_heavy_atom_names'])
  missing=heavy-observed
  if missing or observed-heavy:
   deficits.append({k:v for k,v in residue.items() if k!='observed_heavy_atom_names'} | {
       'missing_ccd_heavy_atom_names':sorted(missing),'missing_declared_leaving_atoms':sorted(missing & leaving),
       'missing_other_ccd_heavy_atoms':sorted(missing-leaving),'extra_heavy_atom_names':sorted(observed-heavy),
       'interpretation':'Inventory only. Polymer terminal leaving atoms and declared reaction leaving groups are distinct from unresolved coordinate omissions.'})
 connections=[x for x in d['connections'] if belongs(x['partner1'],selected) and belongs(x['partner2'],selected)]
 sites=[]
 for atom in d['metal_atoms']:
  if atom['label_chain'] not in selected:continue
  contacts=[]
  for row in connections:
   if row['type']!='metalc':continue
   for own,other in [('partner1','partner2'),('partner2','partner1')]:
    if metal_match(atom,row[own]):
     contacts.append({'connection_id':row['id'],'donor':row[other],
         'metal_symmetry':row[own]['symmetry'],'distance_angstrom':row['deposited_distance_angstrom'],
         'same_explicit_image':row[own]['symmetry'] is not None and row[own]['symmetry']==row[other]['symmetry']})
  ccd=read_json(OUT/'inventories'/'ccd'/f"{atom['component']}.json")
  intra=[x for x in ccd['bonds'] if atom['atom'] in (x['atom_id_1'],x['atom_id_2'])]
  sites.append(atom | {'ccd':ccds[atom['component']], 'deposited_contacts':contacts,
      'deposited_coordinating_waters':[x for x in contacts if x['donor']['component']=='HOH'],
      'intra_component_ccd_metal_bonds':intra,
      'coordination_basis':'Original struct_conn plus explicit same-component CCD bonds; no distance-derived donor assignment.',
      'chosen_oxidation_state':None,'chosen_spin_state':None,'chosen_model_id':None,
      'model_status':'Not assigned by curation; source formal charges are evidence, not an automatically validated bound-state model.'})
 return {'id':code,'label':label,'status':'frozen_source_case_not_prepared','notes':notes,
     'source':d['input'],'rcsb_entry_url':f'https://www.rcsb.org/structure/{code}',
     'pdb_doi':f'https://doi.org/10.2210/pdb{code}/pdb', 'title':d['title'],
     'organisms':d['organisms'],'resolution_angstrom':float(d['resolution_angstrom']),
     'primary_citations':d['primary_citations'],
     'assembly_choice':{'id':'1','declaration':next(x for x in d['assemblies'] if x['id']=='1'),
         'generators':generators,'operators':d['assembly_operators'],
         'label_chains':sorted(selected),'protein_chains':[x for x in d['chains'] if x['label_chain'] in selected],
         'policy':'Use every component and operator in deposited assembly1. Other independent assemblies are excluded by this predeclared assembly choice; original input is preserved.'},
     'retained_nonstandard_residues':retained,'ccd_sources':list(ccds.values()),
     'ccd_coordinate_differences':deficits,'deposited_connections':connections,'metal_sites':sites,
     'excluded_source_label_chains':sorted({x['label_chain'] for x in d['chains']+d['nonstandard_residues']}-selected),
     'future_preparation_policy':{'remove_heterogens':False,'remove_coordinating_waters':False,
         'build_missing_internal_residues':'Only with existing supported sequence and geometry safeguards; retained gaps listed above.',
         'chemical_mutation':'No undeclared protein/ligand mutation or ligand deletion; any accepted bound-state transformation must be recorded and validated.'}}

def main():
 targets=[OUT/'covalent-panel.json',OUT/'metal-panel.json',OUT/'panel-freeze.json']
 if any(p.exists() for p in targets):raise RuntimeError('Frozen panel exists; create a versioned amendment rather than silently overwriting it.')
 stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
 panels={}
 for panel,config in [('covalent',COVALENT),('metal',METAL)]:
  cases=[]
  for values in config:
   code,label=values[:2]
   if panel=='covalent':
    component,adduct,notes=values[2:];case=base_case(code,label,notes)
    target=[x for x in case['deposited_connections'] if x['type']=='covale' and component in (x['partner1']['component'],x['partner2']['component'])]
    if code=='5F19':
     target=[]
     for r in case['retained_nonstandard_residues']:
      if r['component']=='OAS':target.append({'representation':'intra_component_CCD_bond','residue':r,'atoms':['OG','C1A'],'order':'SING','source_ccd':'OAS'})
    assert target,code
    case['covalent_target']={'component':component,'adduct_class':adduct,'connections':target,
        'state_evidence':'Exact deposited bond endpoints, deposited bound atoms, CCD bond graph/leaving flags, and linked primary structure publication. Hydrogens, net charge and missing coordinate decisions remain explicit model requirements.',
        'chosen_product_model':None}
   else:case=base_case(code,label,values[2])
   cases.append(case)
  panels[panel]={'schema_version':1,'frozen_at_utc':stamp,'panel':panel,'case_count':15,
    'selection_basis':'Purposeful chemistry-diversity challenge set from actual deposited structures. Selected without consulting full-case preparation outcomes; early prototype feasibility work is separate. This is not a random PDB survey or an estimated general success rate.',
    'retention_policy':'Preserve original source bytes, declared assembly components and relevant waters; do not omit chemistry to produce a passing result.',
    'cases':cases}
  (OUT/f'{panel}-panel.json').write_text(json.dumps(panels[panel],indent=2)+'\n')
 selected={x[0] for x in COVALENT+METAL}
 reasons={'4IFG':'Mistaken remembered accession: organism is Toxoplasma gondii, outside benign-host scope; excluded after metadata inventory and never modeled.',
    '4ZAU':'No deposited covalent link and substantial multiple sequence gaps; a title mentioning osimertinib does not establish a supported reacted graph.',
    '3QJ9':'No deposited covalent link; alternate FAAH structures with explicit reacted products selected.',
    '3S3P':'Multiple missing loops including a 13-residue gap introduce unrelated repair barriers; preserved as a screened candidate.',
    '1ONE':'Deposited equilibrium substrate/product mixture and numerous mutually overlapping cross-component covale records need dedicated source disambiguation, so not an ordinary metal-model reference.',
    '1CGL':'Same MMP1 family represented by higher-resolution 1HFC with a complete observed backbone.',
    '3JWE':'MAGL represented by higher-resolution 6AX1 assembly1 with a complete observed backbone.',
    '3W2T':'Additional DPP4 example would duplicate the Ser nitrile product class and introduce larger glycan inventory; 3BJM already retains DPP4 glycans.'}
 excluded=[]
 for path in sorted((OUT/'inventories').glob('*.json')):
  d=read_json(path)
  if d['id'] not in selected:excluded.append({'id':d['id'],'source':d.get('input'),'reason':reasons.get(d['id'],'Screened candidate outside fixed thirty-case diversity allocation.'),'inventory':str(path.relative_to(OUT))})
 (OUT/'excluded-candidates.json').write_text(json.dumps({'recorded_at_utc':stamp,'candidates':excluded},indent=2)+'\n')
 hashes={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in targets[:2]}
 inputs={case['source']['path']:case['source']['sha256'] for p in panels.values() for case in p['cases']}
 for p in panels.values():
  for case in p['cases']:
   for ccd in case['ccd_sources']:inputs[ccd['source']['path']]=ccd['source']['sha256']
 targets[2].write_text(json.dumps({'frozen_at_utc':stamp,'panels':hashes,'original_inputs':inputs,
    'excluded_candidates_sha256':hashlib.sha256((OUT/'excluded-candidates.json').read_bytes()).hexdigest(),
    'curation_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'native_outcomes_used_for_selection':False},indent=2)+'\n')
 print(json.dumps(hashes,indent=2))

if __name__=='__main__':main()
