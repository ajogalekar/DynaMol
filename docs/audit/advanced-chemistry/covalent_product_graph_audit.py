"""Graph-only pressure test over the frozen deposited covalent product panel.

No parameters, quantum calculations, MD, or full-complex preparation are done.
CCD charges define explicit candidate states; pH-dependent protonation and
chemical accuracy remain unvalidated. Missing non-leaving atoms are never built.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'docs/audit/advanced-chemistry'
ORDER={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
SIDECHAINS={'CYS':[('CA','CB','C'),('CB','SG','S')],
            'SER':[('CA','CB','C'),('CB','OG','O')],
            'LYS':[('CA','CB','C'),('CB','CG','C'),('CG','CD','C'),('CD','CE','C'),('CE','NZ','N')]}


def value(x):return '' if x in (False,None,'.','?') else str(x)


def source_rows(block):
    columns=block.get_mmcif_category('_atom_site.');result=[]
    for i in range(len(columns['id'])):
        row={key:value(items[i]) for key,items in columns.items()}
        if row.get('pdbx_PDB_model_num','1')!='1' or row['type_symbol']=='H':continue
        row['xyz']=np.array([float(row[k]) for k in ['Cartn_x','Cartn_y','Cartn_z']]);result.append(row)
    return result


def select_residue(rows,identity):
    matches=[r for r in rows if r['label_asym_id']==identity['label_chain'] and r['auth_seq_id']==value(identity['author_residue_id'])
             and r['label_comp_id']==identity['component'] and r.get('pdbx_PDB_ins_code','')==value(identity.get('insertion_code'))]
    if not matches:raise ValueError('No observed heavy atoms for declared identity '+json.dumps(identity))
    alts=sorted({r['label_alt_id'] for r in matches if r['label_alt_id']})
    declared=value(identity.get('altloc'))
    chosen=declared or (max(alts,key=lambda alt:(sum(float(r.get('occupancy') or 1) for r in matches if r['label_alt_id']==alt),-alts.index(alt))) if alts else '')
    selected=[r for r in matches if not r['label_alt_id'] or r['label_alt_id']==chosen]
    atoms={r['label_atom_id']:r for r in selected}
    if len(atoms)!=len(selected):raise ValueError('Duplicate heavy atom identity after coherent altloc selection.')
    return atoms,chosen


def neighbor_identity(rows,protein,delta):
    seq=int(protein['label_sequence_id'])+delta
    matches=[r for r in rows if r['label_asym_id']==protein['label_chain'] and r['label_seq_id']==str(seq)]
    keys={(r['label_asym_id'],r['auth_seq_id'],r['label_comp_id'],r.get('pdbx_PDB_ins_code','')) for r in matches}
    if len(keys)!=1:raise ValueError('Peptide cap requires an unambiguous observed immediate neighbor at label sequence '+str(seq))
    chain,resid,name,insertion=next(iter(keys))
    return {'label_chain':chain,'label_sequence_id':str(seq),'author_residue_id':resid,'component':name,'insertion_code':insertion}


def build(case,connection,rows,ccd):
    integrated=connection.get('representation')=='intra_component_CCD_bond'
    if not integrated and any(p['component']==case['covalent_target']['component'] and value(p.get('label_sequence_id')) for p in (connection['partner1'],connection['partner2'])):
        module_spec=importlib.util.spec_from_file_location('covalent_polymer_graph',BASE/'covalent_polymer_graph.py')
        module=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(module)
        return module.expand(case,connection,rows,build,select_residue)
    if integrated:protein=connection['residue'];ligand=protein;protein_end,ligand_end=connection['atoms']
    else:
        partners=[connection['partner1'],connection['partner2']]
        if any(value(p.get('symmetry')) not in ('','1_555') for p in partners):raise ValueError('Symmetry-transformed product requires explicit assembly coordinates.')
        proteins=[p for p in partners if p['component'] in SIDECHAINS or p['component']=='THR']
        if len(proteins)!=1:raise ValueError('This product needs a connected polymer/multi-residue model, not a one-residue cap.')
        protein=proteins[0];ligand=partners[1-partners.index(protein)];protein_end=protein['atom'];ligand_end=ligand['atom']
    if protein['component']=='THR':raise ValueError('Tetrahedral boronate and N-terminal Thr charge state require an explicitly reviewed joint model; neutral trigonal boron is not reused.')
    if not integrated and value(ligand.get('label_sequence_id')):raise ValueError('Polymeric ligand retains additional peptide endpoints; a single ligand-component cap would truncate its chemistry.')
    observed_protein,altp=select_residue(rows,protein);observed_ligand,altl=select_residue(rows,ligand)
    left=neighbor_identity(rows,protein,-1);right=neighbor_identity(rows,protein,1)
    observed_left,altleft=select_residue(rows,left);observed_right,altright=select_residue(rows,right)
    atoms=ccd.get_mmcif_category('_chem_comp_atom.');bonds=ccd.get_mmcif_category('_chem_comp_bond.')
    heavy={name:i for i,name in enumerate(atoms['atom_id']) if atoms['type_symbol'][i]!='H'}
    missing=set(heavy)-set(observed_ligand);extra=set(observed_ligand)-set(heavy)
    allowed={name for name in missing if atoms['pdbx_leaving_atom_flag'][heavy[name]]=='Y'}
    if missing-allowed or extra:raise ValueError(f'Unresolved ligand heavy inventory: missing non-leaving {sorted(missing-allowed)}, extra {sorted(extra)}')
    mol=Chem.RWMol();coordinates=[];records=[];lookup={}
    def add(name,symbol,xyz,role,source,charge=0):
        atom=Chem.Atom(symbol.title());atom.SetFormalCharge(int(charge));index=mol.AddAtom(atom);lookup[name]=index
        coordinates.append(xyz);records.append({'index':index,'name':name,'role':role,'source':source});return index
    def connect(a,b,order=Chem.BondType.SINGLE):mol.AddBond(lookup[a],lookup[b],order)
    for name,symbol,residue,atom in [('ACE_CH3','C',observed_left,'CA'),('ACE_C','C',observed_left,'C'),('ACE_O','O',observed_left,'O'),
                                    ('NME_N','N',observed_right,'N'),('NME_CH3','C',observed_right,'CA')]:
        if atom not in residue:raise ValueError('Observed cap anchor atom is missing: '+atom)
        add(name,symbol,residue[atom]['xyz'],'temporary_cap',{'label_chain':residue[atom]['label_asym_id'],'author_resid':residue[atom]['auth_seq_id'],'atom':atom})
    connect('ACE_CH3','ACE_C');connect('ACE_C','ACE_O',Chem.BondType.DOUBLE);connect('NME_N','NME_CH3')
    if not integrated:
        required={'N':'N','CA':'C','C':'C','O':'O',**{b:symbol for _,b,symbol in SIDECHAINS[protein['component']]}}
        if set(observed_protein)!=set(required):raise ValueError('Protein residue requires complete standard heavy inventory before graph construction.')
        for name,symbol in required.items():add('P_'+name,symbol,observed_protein[name]['xyz'],'observed_product',{'identity':protein,'atom':name})
        protein_bonds=[('N','CA'),('CA','C')]+[(a,b) for a,b,_ in SIDECHAINS[protein['component']]]
        for a,b in protein_bonds:connect('P_'+a,'P_'+b)
        connect('P_C','P_O',Chem.BondType.DOUBLE)
    for name,index in heavy.items():
        if name in observed_ligand:add('L_'+name,atoms['type_symbol'][index],observed_ligand[name]['xyz'],'observed_product',{'identity':ligand,'atom':name},atoms['charge'][index])
    for a,b,order in zip(bonds['atom_id_1'],bonds['atom_id_2'],bonds['value_order']):
        if 'L_'+a in lookup and 'L_'+b in lookup:connect('L_'+a,'L_'+b,ORDER[order])
    prefix='L_' if integrated else 'P_'
    connect('ACE_C',prefix+'N');connect(prefix+'C','NME_N')
    if not integrated:connect('P_'+protein_end,'L_'+ligand_end)
    endpoint_names=[prefix+protein_end,'L_'+ligand_end]
    result=mol.GetMol();Chem.SanitizeMol(result)
    conformer=Chem.Conformer(len(coordinates))
    for i,xyz in enumerate(coordinates):conformer.SetAtomPosition(i,xyz)
    result.AddConformer(conformer);Chem.AssignStereochemistryFrom3D(result);Chem.AssignStereochemistry(result,cleanIt=True,force=True)
    alpha=result.GetAtomWithIdx(lookup[prefix+'CA']);expected_alpha='R' if protein['component']=='CYS' else 'S'
    if not alpha.HasProp('_CIPCode') or alpha.GetProp('_CIPCode')!=expected_alpha:raise ValueError('Observed alpha-carbon stereochemistry differs from the declared L residue.')
    stereo=[]
    for atom in result.GetAtoms():
        if atom.HasProp('_CIPCode'):
            record=records[atom.GetIdx()];ccd_name=record['name'][2:] if record['name'].startswith('L_') else None
            declared=atoms['pdbx_stereo_config'][heavy[ccd_name]] if ccd_name in heavy else None
            stereo.append({'name':record['name'],'product_cip_from_observed_coordinates':atom.GetProp('_CIPCode'),'ccd_stereo_label':declared})
    result=Chem.AddHs(result,addCoords=True)
    heavy_xyz=np.asarray(result.GetConformer().GetPositions())[:len(coordinates)]
    if not np.array_equal(heavy_xyz,np.asarray(coordinates)):raise ValueError('Hydrogen addition changed observed product or cap coordinates.')
    endpoints=[]
    for name in endpoint_names:
        atom=result.GetAtomWithIdx(lookup[name]);endpoints.append({'name':name,'element':atom.GetSymbol(),'formal_charge':atom.GetFormalCharge(),
            'bond_order_sum':float(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())),
            'attached_hydrogens':sum(a.GetAtomicNum()==1 for a in atom.GetNeighbors()),
            'heavy_neighbors':[records[a.GetIdx()]['name'] for a in atom.GetNeighbors() if a.GetAtomicNum()!=1]})
    cap_distances=[float(np.linalg.norm(np.asarray(coordinates[lookup['ACE_C']])-coordinates[lookup[prefix+'N']])),
                   float(np.linalg.norm(np.asarray(coordinates[lookup[prefix+'C']])-coordinates[lookup['NME_N']]))]
    if any(not 1.1<distance<1.8 for distance in cap_distances):raise ValueError('Observed peptide cap boundary is broken: '+str(cap_distances))
    report={'status':'graph_candidate_constructed','chemical_state_accepted':False,'full_preparation_accepted':False,
        'scope':'Graph/valence/identity candidate under deposited CCD formal charges, with regenerated hydrogens; pH protonation, parameterization and full retained-complex chemistry are not validated.',
        'formal_charge_from_ccd_candidate':Chem.GetFormalCharge(result),'atom_count_with_hydrogens':result.GetNumAtoms(),'heavy_atom_count':len(coordinates),
        'explicit_hydrogen_count':sum(a.GetAtomicNum()==1 for a in result.GetAtoms()),'original_product_heavy_coordinates_preserved_exactly':True,
        'omitted_absent_ccd_leaving_atoms':sorted(allowed),'heavy_bond_order_changes':[],
        'observed_connection_distance_angstrom':float(np.linalg.norm(np.asarray(coordinates[lookup[endpoint_names[0]]])-coordinates[lookup[endpoint_names[1]]])),
        'endpoint_valence_and_protons':endpoints,'observed_stereochemistry':stereo,'alpha_stereo':expected_alpha,
        'double_bond_stereochemistry':[{'atoms':[records[b.GetBeginAtomIdx()]['name'],records[b.GetEndAtomIdx()]['name']],
            'assigned_from_observed_coordinates':str(b.GetStereo())} for b in result.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE],
        'caps':{'type':'ACE/NME','source_left':left,'source_right':right,'peptide_join_lengths_angstrom':cap_distances,'never_export_to_protein':True},
        'chosen_altlocs':{'protein':altp,'ligand':altl,'left_cap':altleft,'right_cap':altright},'heavy_atom_map':records,
        'mapped_product_smiles':Chem.MolToSmiles(Chem.RemoveHs(result),isomericSmiles=True)}
    return result,report


def run(output):
    output=output.resolve();output.mkdir(parents=True,exist_ok=False)
    panel_path=BASE/'covalent-panel.json';panel=json.loads(panel_path.read_text());results=[]
    for case in panel['cases']:
        source=BASE/case['source']['path']
        if hashlib.sha256(source.read_bytes()).hexdigest()!=case['source']['sha256']:raise ValueError('Frozen source hash changed.')
        component=case['covalent_target']['component'];ccd_info=next(x for x in case['ccd_sources'] if x['id']==component);ccd_path=BASE/ccd_info['source']['path']
        if hashlib.sha256(ccd_path.read_bytes()).hexdigest()!=ccd_info['source']['sha256']:raise ValueError('Frozen CCD hash changed.')
        rows=source_rows(gemmi.cif.read_file(str(source)).sole_block());ccd=gemmi.cif.read_file(str(ccd_path)).sole_block();folder=output/case['id'];folder.mkdir()
        result={'id':case['id'],'label':case['label'],'adduct_class':case['covalent_target']['adduct_class'],
            'source_sha256':case['source']['sha256'],'ccd_sha256':ccd_info['source']['sha256'],'full_preparation_accepted':False,'connections':[]}
        for index,connection in enumerate(case['covalent_target']['connections']):
            name=connection.get('id',f'integrated-{index+1}')
            try:
                molecule,report=build(case,connection,rows,ccd)
                (folder/(name+'.sdf')).write_text(Chem.MolToMolBlock(molecule)+'\n$$$$\n')
            except Exception as error:report={'status':'blocked','reason':str(error),'error_type':type(error).__name__,'full_preparation_accepted':False}
            report['connection']=connection;(folder/(name+'.json')).write_text(json.dumps(report,indent=2)+'\n');result['connections'].append({'name':name,**{k:v for k,v in report.items() if k in ['status','reason','heavy_atom_count','atom_count_with_hydrogens']}})
        result['all_target_graph_candidates_constructed']=all(x['status']=='graph_candidate_constructed' for x in result['connections'])
        (folder/'result.json').write_text(json.dumps(result,indent=2)+'\n');results.append(result)
    summary={'scope':'Frozen15-case deposited-product graph audit only; no force-field or full-preparation success claims.',
        'panel_sha256':hashlib.sha256(panel_path.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'polymer_expansion_source_sha256':hashlib.sha256((BASE/'covalent_polymer_graph.py').read_bytes()).hexdigest(),
        'additional_standard_ccd_manifest_sha256':hashlib.sha256((BASE/'polymer-ccd/source-manifest.json').read_bytes()).hexdigest(),
        'case_count':len(results),'all_target_graphs_constructed':sum(x['all_target_graph_candidates_constructed'] for x in results),
        'blocked_cases':[x['id'] for x in results if not x['all_target_graph_candidates_constructed']],'cases':results}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);args=parser.parse_args()
    print(json.dumps({k:v for k,v in run(args.output).items() if k!='cases'},indent=2))
