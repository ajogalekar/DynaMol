"""Source/sequence-supported complete peptide-ligand graph expansion.

Combines every residue of the declared peptide entity before attachment to
the protein cap. Virtual names are internal only; source identities are
restored in the output atom map. No PDB-specific residues or mutations.
"""
from __future__ import annotations
import copy
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np

BASE=Path(__file__).resolve().parent


def value(x):return '' if x in (False,None,'.','?') else str(x)


def expand(case,requested,rows,build_single,select_residue):
    component=case['covalent_target']['component']
    sites=[p for p in [requested['partner1'],requested['partner2']] if p['component']==component and value(p.get('label_sequence_id'))]
    if len(sites)!=1:raise ValueError('Polymeric target residue identity is ambiguous.')
    ligand_site=sites[0];chain=ligand_site['label_chain'];links=[]
    for connection in case['deposited_connections']:
        if connection.get('type')!='covale':continue
        sides=[connection['partner1'],connection['partner2']]
        if sum(p['label_chain']==chain for p in sides)==1:links.append(connection)
    if len(links)!=1:raise ValueError('Complete peptide ligand has multiple/ambiguous external covalent links; a multi-residue site model is required.')
    primary=copy.deepcopy(links[0]);protein=next(p for p in [primary['partner1'],primary['partner2']] if p['label_chain']!=chain)
    if protein['component'] not in ('CYS','SER','LYS'):raise ValueError('External peptide-adduct protein residue requires another reviewed cap chemistry.')
    bound=next(p for p in [primary['partner1'],primary['partner2']] if p['label_chain']==chain)
    path=BASE/case['source']['path']
    if hashlib.sha256(path.read_bytes()).hexdigest()!=case['source']['sha256']:raise ValueError('Frozen source changed before polymer expansion.')
    block=gemmi.cif.read_file(str(path)).sole_block();polyrows=[r for r in rows if r['label_asym_id']==chain]
    entities={r['label_entity_id'] for r in polyrows}
    if len(entities)!=1:raise ValueError('Peptide ligand chain does not map to one declared entity.')
    entity=next(iter(entities));polymers=block.get_mmcif_category('_entity_poly.')
    kinds=[kind for eid,kind in zip(polymers['entity_id'],polymers['type']) if eid==entity]
    if kinds not in (['polypeptide(L)'],['polypeptide(D)']):raise ValueError('Polymeric ligand is not a conventional peptide entity.')
    seq=block.get_mmcif_category('_entity_poly_seq.');sequence=[]
    for eid,number,monomer in zip(seq['entity_id'],seq['num'],seq['mon_id']):
        if eid==entity:sequence.append((int(number),monomer))
    sequence.sort()
    if len({number for number,_ in sequence})!=len(sequence) or [n for n,_ in sequence]!=list(range(1,len(sequence)+1)):
        raise ValueError('Peptide ligand sequence is ambiguous or noncontiguous.')
    if {int(r['label_seq_id']) for r in polyrows}!={number for number,_ in sequence}:
        raise ValueError('Peptide ligand has unobserved residues; none will be trimmed or guessed.')
    extra_manifest=json.loads((BASE/'polymer-ccd/source-manifest.json').read_text())
    ccd_sources={x['id']:(BASE/x['source']['path'],x['source']['sha256']) for x in case['ccd_sources']}
    ccd_sources.update({x['id']:(BASE.parents[2]/x['path'],x['sha256']) for x in extra_manifest['sources']})
    atom_columns={key:[] for key in ['atom_id','type_symbol','charge','pdbx_leaving_atom_flag','pdbx_stereo_config']}
    bond_columns={key:[] for key in ['atom_id_1','atom_id_2','value_order']}
    virtual_rows=[];source_map={};identities={};sequence_edges=[];source_records=[];observed_by_seq={}
    for number,name in sequence:
        matches=[r for r in polyrows if r['label_seq_id']==str(number)]
        identity_keys={(r['auth_seq_id'],r['label_comp_id'],r.get('pdbx_PDB_ins_code','')) for r in matches}
        if len(identity_keys)!=1:raise ValueError('Ambiguous observed peptide residue identity.')
        resid,observed_name,insertion=next(iter(identity_keys))
        if observed_name!=name:raise ValueError('Observed peptide residue differs from declared sequence.')
        identity={'label_chain':chain,'label_sequence_id':str(number),'author_chain':matches[0]['auth_asym_id'],
                  'author_residue_id':resid,'component':name,'insertion_code':insertion}
        selected,alt=select_residue(rows,identity);observed_by_seq[number]=selected;identities[number]=identity
        if name not in ccd_sources:raise ValueError('A standard/modified peptide CCD source must be explicitly recorded before graph expansion: '+name)
        source,digest=ccd_sources[name]
        if hashlib.sha256(source.read_bytes()).hexdigest()!=digest:raise ValueError('Peptide CCD source hash changed.')
        ccd=gemmi.cif.read_file(str(source)).sole_block();atoms=ccd.get_mmcif_category('_chem_comp_atom.');bonds=ccd.get_mmcif_category('_chem_comp_bond.')
        heavy={atom for atom,symbol in zip(atoms['atom_id'],atoms['type_symbol']) if symbol!='H'}
        for i,atom in enumerate(atoms['atom_id']):
            if atom not in heavy:continue
            for key in atom_columns:atom_columns[key].append(f'R{number}_{atom}' if key=='atom_id' else atoms[key][i])
        for a,b,order in zip(bonds['atom_id_1'],bonds['atom_id_2'],bonds['value_order']):
            if a in heavy and b in heavy:
                for key,val in [('atom_id_1',f'R{number}_{a}'),('atom_id_2',f'R{number}_{b}'),('value_order',order)]:bond_columns[key].append(val)
        for atom,row in selected.items():
            mapped=f'R{number}_{atom}';source_map[mapped]={'identity':identity,'atom':atom,'selected_altloc':alt}
            virtual=copy.deepcopy(row);virtual.update(label_comp_id='POLY',auth_seq_id=bound['author_residue_id'],label_seq_id='',label_atom_id=mapped,label_alt_id='',pdbx_PDB_ins_code='');virtual_rows.append(virtual)
        source_records.append({'component':name,'source':str(source),'sha256':digest,'identity':identity,'selected_altloc':alt})
    for (left,_),(right,_) in zip(sequence,sequence[1:]):
        if 'C' not in observed_by_seq[left] or 'N' not in observed_by_seq[right]:raise ValueError('Observed peptide connector atom is missing.')
        distance=float(np.linalg.norm(observed_by_seq[left]['C']['xyz']-observed_by_seq[right]['N']['xyz']))
        if not 1.1<distance<1.8:raise ValueError('Declared peptide sequence has a broken observed connector: '+str(distance))
        a,b=f'R{left}_C',f'R{right}_N'
        for key,val in [('atom_id_1',a),('atom_id_2',b),('value_order','SING')]:bond_columns[key].append(val)
        sequence_edges.append({'atoms':[source_map[a],source_map[b]],'distance_angstrom':distance,'basis':'Consecutive positions in deposited peptide entity sequence plus observed C/N connector atoms.'})
    virtual=gemmi.cif.Block('POLY');virtual.set_mmcif_category('_chem_comp_atom.',atom_columns);virtual.set_mmcif_category('_chem_comp_bond.',bond_columns)
    virtual_identity={**bound,'component':'POLY','label_sequence_id':None,'insertion_code':'','atom':f"R{bound['label_sequence_id']}_{bound['atom']}"}
    child_connection=copy.deepcopy(primary)
    for key in ('partner1','partner2'):
        if child_connection[key]['label_chain']==chain:child_connection[key]=virtual_identity
    input_rows=[r for r in rows if r['label_asym_id']!=chain]+virtual_rows
    molecule,report=build_single(case,child_connection,input_rows,virtual)
    for record in report['heavy_atom_map']:
        if record['name'].startswith('L_'):
            record['source']=source_map[record['name'][2:]]
    report['polymeric_ligand']={'expanded_entire_declared_chain':True,'label_chain':chain,'entity_id':entity,
        'sequence':sequence,'source_ccds':source_records,'peptide_edges':sequence_edges,
        'primary_external_connection':primary,'requested_connection':requested,
        'observed_ligand_heavy_atoms_retained':len(virtual_rows),'virtual_names_are_internal_only':True}
    # Every additional deposited peptide target connection must be represented
    # by the complete graph, including the terminal modified-residue linkage.
    indices={record['name']:record['index'] for record in report['heavy_atom_map']}
    for connection in case['covalent_target']['connections']:
        endpoints=[]
        for endpoint in (connection['partner1'],connection['partner2']):
            if endpoint['label_chain']==chain:endpoints.append(f"L_R{endpoint['label_sequence_id']}_{endpoint['atom']}")
            elif endpoint['label_chain']==protein['label_chain'] and endpoint['author_residue_id']==protein['author_residue_id']:endpoints.append('P_'+endpoint['atom'])
            else:raise ValueError('A deposited target connection leaves the expanded covalent component.')
        if any(name not in indices for name in endpoints) or molecule.GetBondBetweenAtoms(*[indices[name] for name in endpoints]) is None:
            raise ValueError('Expanded graph is missing a deposited covalent connection.')
    return molecule,report
