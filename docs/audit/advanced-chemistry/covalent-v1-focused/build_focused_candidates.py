"""Explicit human Cys–acrylamide adduct candidates for the focused v1 audit.

Original CIF/CCD bytes remain immutable. Deposited linkage plus a narrowly
matched acrylamide graph distinguishes a free CCD from an already-reacted CCD.
Protonation choices are explicit test hypotheses, not bound-state pKa predictions.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from prepare_adduct import native, digest, mol2_sections, prepare


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def declared_connection(block):
    cols = block.get_mmcif_category('_struct_conn.')
    found = []
    for i, kind in enumerate(cols['conn_type_id']):
        if kind != 'covale':
            continue
        row = {k:v[i] for k,v in cols.items()}
        if row['ptnr1_label_comp_id'] != 'CYS' or row['ptnr1_label_atom_id'] != 'SG':
            continue
        partners = []
        for n in [1,2]:
            prefix = f'ptnr{n}_'
            partners.append({'label_chain':row[prefix+'label_asym_id'],
                'label_sequence_id':row[prefix+'label_seq_id'],
                'component':row[prefix+'label_comp_id'], 'atom':row[prefix+'label_atom_id'],
                'author_chain':row[prefix+'auth_asym_id'],
                'author_residue_id':row[prefix+'auth_seq_id'],
                'insertion_code':row.get(f'pdbx_ptnr{n}_PDB_ins_code'),
                'altloc':row.get(f'pdbx_ptnr{n}_label_alt_id'),
                'symmetry':row[prefix+'symmetry']})
        found.append({'id':row['id'], 'partner1':partners[0], 'partner2':partners[1]})
    if len(found) != 1:
        raise ValueError('Focused case needs exactly one explicit cysteine–ligand connection')
    return found[0]


def reacted_ccd(original, beta, protonated_atom=None):
    # Work on a fresh in-memory copy; the downloaded CCD is never changed.
    doc = gemmi.cif.read_string(original.as_string()); ccd = doc.sole_block()
    atoms = ccd.get_mmcif_category('_chem_comp_atom.')
    bonds = ccd.get_mmcif_category('_chem_comp_bond.')
    elements = dict(zip(atoms['atom_id'], atoms['type_symbol']))
    edges = list(zip(bonds['atom_id_1'], bonds['atom_id_2'], bonds['value_order']))
    def neighbors(name):
        return [(b if a == name else a, order, i) for i,(a,b,order) in enumerate(edges)
                if (a == name or b == name) and elements[b if a == name else a] != 'H']
    adjacent = neighbors(beta)
    if elements[beta] != 'C' or len(adjacent) != 1 or elements[adjacent[0][0]] != 'C':
        raise ValueError('Not a terminal carbon in the supported acrylamide linkage class')
    alpha, order, bond_index = adjacent[0]
    carbonyls = [n for n,o,_ in neighbors(alpha) if n != beta and elements[n] == 'C'
                 and any(elements[x]=='O' and bo=='DOUB' for x,bo,_ in neighbors(n))
                 and any(elements[x]=='N' and bo=='SING' for x,bo,_ in neighbors(n))]
    if len(carbonyls) != 1 or order not in ['SING','DOUB']:
        raise ValueError('The observed linkage is outside the Cys–acrylamide test scope')
    edits = []
    if order == 'DOUB':
        bonds['value_order'][bond_index] = 'SING'
        edits.append({'atoms':[alpha,beta], 'before':'DOUB', 'after':'SING',
            'reason':'Deposited cysteine linkage and terminal acrylamide graph specify the saturated Michael-addition product.'})
    ccd.set_mmcif_category('_chem_comp_bond.', bonds)
    if protonated_atom:
        i = atoms['atom_id'].index(protonated_atom)
        ns = neighbors(protonated_atom)
        if elements[protonated_atom] != 'N' or len(ns) != 3 or any(elements[n]!='C' or o!='SING' for n,o,_ in ns):
            raise ValueError('Declared protonation site is not the expected tertiary amine')
        atoms['charge'][i] = '1'
        ccd.set_mmcif_category('_chem_comp_atom.', atoms)
    return ccd, edits


def generate(root, code, amber):
    graph = module(root/'sources/covalent_product_graph_audit.py', 'focused_graph')
    coverage = module(root/'sources/covalent_native_coverage.py', 'focused_coverage')
    original = root/'inputs'/f'{code}.cif'
    block = gemmi.cif.read_file(str(original)).sole_block()
    connection = declared_connection(block)
    component = connection['partner2']['component']
    source_ccd = root/'inputs'/f'{component}.cif'
    # Free-drug experimental data identify this aliphatic amine as basic.
    # The bound-adduct state is still an explicitly declared +1 hypothesis.
    protonated = 'N2' if component == 'YY3' else None
    ccd, edits = reacted_ccd(gemmi.cif.read_file(str(source_ccd)).sole_block(),
                            connection['partner2']['atom'], protonated)
    case = {'covalent_target':{'component':component}}
    mol, info = graph.build(case, connection, graph.source_rows(block), ccd)
    folder = root/(code.lower()+'-baseline-v1')
    folder.mkdir(parents=True, exist_ok=False)
    (folder/'capped-native').mkdir(); (folder/'resp').mkdir(); (folder/'interface-native-scope-v4').mkdir()
    heavy = info['heavy_atom_map']
    # Relieve only generated hydrogen contacts before semiempirical charging.
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    if props is not None:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)
        for row in heavy:
            ff.AddFixedPoint(row['index'])
        ff.Initialize(); ff.Minimize(maxIts=400)
    sdf = folder/'capped-native/input.sdf'
    sdf.write_text(Chem.MolToMolBlock(mol)+'\n$$$$\n')
    info.update({'source_sha256':digest(original), 'ccd_source_sha256':digest(source_ccd),
        'connection':connection, 'explicit_product_edits':edits,
        'protonated_ccd_atom':protonated,
        'state_basis':'6JXT: singly protonated distal aliphatic amine candidate; 5P9J: neutral deposited-adduct candidate. These are not computed bound-state pKas.'})
    (folder/'graph.json').write_text(json.dumps(info,indent=2)+'\n')
    commands = [native(amber, folder/'capped-native', ['antechamber','-i','input.sdf','-fi','sdf',
        '-o','charged.mol2','-fo','mol2','-at','gaff2','-c','bcc','-nc',str(Chem.GetFormalCharge(mol)),
        '-s','2','-seq','n','-pf','n'], seconds=900)]
    typed = mol2_sections(folder/'capped-native/charged.mol2')
    if len(typed['ATOM']) != mol.GetNumAtoms():
        raise ValueError('Native charging changed atom inventory')
    input_edges = {tuple(sorted((b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1))) for b in mol.GetBonds()}
    if {tuple(sorted((int(b[1]),int(b[2])))) for b in typed['BOND']} != input_edges:
        raise ValueError('Native charging changed product connectivity')
    old = json.loads((root/'baseline-6oim/resp/constraints.json').read_text())
    template = {(a['template'],a['atom']):a for a in old['frozen_atoms']}
    atom_map, frozen, caps = [], [], []
    role_template = {'ACE_CH3':('ACE','CH3'),'ACE_C':('ACE','C'),'ACE_O':('ACE','O'),
                     'NME_N':('NME','N'),'NME_CH3':('NME','C')}
    for i,a in enumerate(mol.GetAtoms()):
        if i < len(heavy):
            row = heavy[i]; name = row['name']; origin = row['source']
            iscap = row['role'] == 'temporary_cap'
            protein = name.startswith('P_')
            record = {'index':i, 'name':name, 'role':'temporary-cap' if iscap else 'observed-adduct',
                'resname':'CYS' if protein else component, 'atom':name[2:] if protein or name.startswith('L_') else name,
                'source':origin}
            atom_map.append(record)
            if iscap:
                caps.append(i); t,n = role_template[name]; frozen.append(dict(template[t,n],index=i))
            elif protein and name[2:] in ['N','CA','C','O']:
                frozen.append(dict(template['CYS',name[2:]],index=i))
        else:
            parent = a.GetNeighbors()[0].GetIdx(); pname = heavy[parent]['name']; iscap = parent in caps
            atom_map.append({'index':i,'name':f'H{i+1}','role':'generated-hydrogen','parent_index':parent,'element':'H'})
            if iscap:
                caps.append(i); t,n = ('ACE','H1') if pname.startswith('ACE') else ('NME','H' if pname=='NME_N' else 'H1')
                frozen.append(dict(template[t,n],index=i))
            elif pname in ['P_N','P_CA']:
                frozen.append(dict(template['CYS','H' if pname=='P_N' else 'HA'],index=i))
    constraints = {'atom_count':mol.GetNumAtoms(),'formal_charge':Chem.GetFormalCharge(mol),
                   'frozen_atoms':frozen,'cap_indices':sorted(caps),
                   'scope':'Canonical type/boundary inventory only; candidate charges follow native AM1-BCC/PREPGEN, not RESP.'}
    (folder/'capped-adduct.json').write_text(json.dumps({'atom_map':atom_map},indent=2)+'\n')
    (folder/'resp/constraints.json').write_text(json.dumps(constraints,indent=2)+'\n')
    override = {r['index']:r['type'] for r in frozen}
    lines = (folder/'capped-native/charged.mol2').read_text().splitlines(); inside=False
    for j,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):
            inside=line.endswith('ATOM'); continue
        if inside and line.strip():
            row=line.split(); i=int(row[0])-1
            if i in override:
                row[5]=override[i]; lines[j]=' '.join(row)
    interface = folder/'interface-native-scope-v4'
    (interface/'hybrid.mol2').write_text('\n'.join(lines)+'\n')
    commands.append(native(amber, interface, ['parmchk2','-i','hybrid.mol2','-f','mol2','-o',
        'gaff-all.frcmod','-s','gaff2','-a','Y']))
    parsed = coverage.parse_frcmod(interface/'gaff-all.frcmod')
    by_line = {x['line']:x for x in parsed['entries']}
    canon = set(override.values())|{'X'}; selected=[]
    for number,line in enumerate((interface/'gaff-all.frcmod').read_text().splitlines(),1):
        row=by_line.get(number)
        if row and set(row['atom_types']) <= canon:
            continue
        if row and (row['attention'] or row['invalid']):
            raise ValueError('Unassigned or invalid adduct parameter: '+line)
        selected.append(line)
    (interface/'gaff-adduct.frcmod').write_text('\n'.join(selected)+'\n')
    (folder/'parameter-assignments.json').write_text(json.dumps(parsed,indent=2)+'\n')
    (folder/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    return prepare(folder, root/(code.lower()+'-prepgen-v1'), amber)


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--amber',type=Path,required=True);p.add_argument('--case',choices=['5P9J','6JXT','6JX0'],required=True)
    a=p.parse_args(); print(json.dumps(generate(a.root,a.case,a.amber),indent=2))
