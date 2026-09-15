"""Graph contracts on the frozen public human/rat covalent product cases."""
import copy
import importlib.util
import json
from pathlib import Path

import gemmi
import numpy as np
import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('covalent_graph',HERE/'covalent_product_graph_audit.py')
graph=importlib.util.module_from_spec(spec);spec.loader.exec_module(graph)
PANEL=json.loads((HERE/'covalent-panel.json').read_text())


def fixture(name):
    case=next(case for case in PANEL['cases'] if case['id']==name)
    rows=graph.source_rows(gemmi.cif.read_file(str(HERE/case['source']['path'])).sole_block())
    source=next(x for x in case['ccd_sources'] if x['id']==case['covalent_target']['component'])
    ccd=gemmi.cif.read_file(str(HERE/source['source']['path'])).sole_block()
    return case,copy.deepcopy(case['covalent_target']['connections'][0]),rows,ccd


@pytest.mark.parametrize('name',['6OIM','6UT0','5P9J','6J6M'])
def test_deposited_michael_graph_has_reacted_protons_and_preserves_coordinates(name):
    mol,report=graph.build(*fixture(name))
    assert report['status']=='graph_candidate_constructed'
    assert report['chemical_state_accepted'] is False
    assert report['full_preparation_accepted'] is False
    assert report['heavy_bond_order_changes']==[]
    assert report['original_product_heavy_coordinates_preserved_exactly']
    sulfur,carbon=report['endpoint_valence_and_protons']
    assert sulfur['element']=='S' and sulfur['attached_hydrogens']==0 and sulfur['bond_order_sum']==2
    assert carbon['element']=='C' and carbon['attached_hydrogens']==2 and carbon['bond_order_sum']==4
    assert mol.GetNumAtoms()==report['atom_count_with_hydrogens']


def test_missing_declared_leaving_group_is_distinct_from_missing_retained_atom():
    args=fixture('6AX1');_,report=graph.build(*args)
    assert set(report['omitted_absent_ccd_leaving_atoms'])=={'O3','C13','C14','C15','F1','F2','F3','F4','F5','F6'}
    case,connection,rows,ccd=args
    ligand=connection['partner2'];rows=[row for row in rows if not(row['label_asym_id']==ligand['label_chain'] and row['label_atom_id']=='C12')]
    with pytest.raises(ValueError,match='missing non-leaving'):
        graph.build(case,connection,rows,ccd)


def test_explicit_label_identity_prevents_wrong_chain_attachment():
    case,connection,rows,ccd=fixture('6OIM');connection['partner1']['label_chain']='UNDECLARED'
    with pytest.raises(ValueError,match='No observed heavy atoms'):
        graph.build(case,connection,rows,ccd)


def test_wrong_carbonyl_endpoint_is_rejected_by_product_valence():
    case,connection,rows,ccd=fixture('6OIM');connection['partner2']['atom']='C23'
    with pytest.raises(Exception,match='valence'):
        graph.build(case,connection,rows,ccd)


def test_changed_observed_alpha_stereo_is_rejected():
    case,connection,rows,ccd=fixture('6OIM');target=connection['partner1']
    selected,_=graph.select_residue(rows,target)
    ca,n,c,cb=[selected[name]['xyz'] for name in ['CA','N','C','CB']]
    normal=np.cross(n-ca,c-ca);normal/=np.linalg.norm(normal)
    selected['CB']['xyz']=cb-2*np.dot(cb-ca,normal)*normal
    with pytest.raises(ValueError,match='alpha-carbon stereochemistry'):
        graph.build(case,connection,rows,ccd)


def test_boronate_remains_an_explicit_electronic_state_requirement():
    with pytest.raises(ValueError,match='boronate'):graph.build(*fixture('5LF3'))


def test_complete_peptide_ligand_retains_sequence_and_both_deposited_connections():
    mol,report=graph.build(*fixture('2H5I'))
    polymer=report['polymeric_ligand']
    assert polymer['sequence']==[(1,'ACE'),(2,'ASP'),(3,'GLU'),(4,'VAL'),(5,'ASJ')]
    assert polymer['observed_ligand_heavy_atoms_retained']==35
    assert len(polymer['peptide_edges'])==4
    assert report['heavy_atom_count']==46
    indices={row['name']:row['index'] for row in report['heavy_atom_map']}
    assert mol.GetBondBetweenAtoms(indices['P_SG'],indices['L_R5_C']) is not None
    assert mol.GetBondBetweenAtoms(indices['L_R4_C'],indices['L_R5_N']) is not None
    actual=[row['source']['identity']['component'] for row in report['heavy_atom_map'] if row['name'].startswith('L_')]
    assert set(actual)=={'ACE','ASP','GLU','VAL','ASJ'}
    assert 'POLY' not in actual
    assert report['original_product_heavy_coordinates_preserved_exactly']
    assert report['full_preparation_accepted'] is False


def test_missing_peptide_residue_is_never_truncated():
    case,connection,rows,ccd=fixture('2H5I')
    rows=[row for row in rows if not(row['label_asym_id']=='C' and row['label_seq_id']=='3')]
    with pytest.raises(ValueError,match='unobserved residues'):graph.build(case,connection,rows,ccd)


def test_broken_observed_peptide_connector_is_rejected():
    case,connection,rows,ccd=fixture('2H5I')
    row=next(row for row in rows if row['label_asym_id']=='C' and row['label_seq_id']=='3' and row['label_atom_id']=='C')
    row['xyz']=row['xyz']+np.array([10.,0,0])
    with pytest.raises(ValueError,match='broken observed connector'):graph.build(case,connection,rows,ccd)


def test_extra_external_peptide_crosslink_needs_a_multi_residue_model():
    case,connection,rows,ccd=fixture('2H5I');case=copy.deepcopy(case)
    case['deposited_connections'].append(copy.deepcopy(connection))
    with pytest.raises(ValueError,match='multiple/ambiguous external covalent links'):graph.build(case,connection,rows,ccd)


def test_integrated_modified_residue_keeps_internal_product_bond():
    mol,report=graph.build(*fixture('5F19'))
    oxygen,carbon=report['endpoint_valence_and_protons']
    assert oxygen['name']=='L_OG' and oxygen['attached_hydrogens']==0
    assert carbon['name']=='L_C1A' and carbon['bond_order_sum']==4
    assert report['omitted_absent_ccd_leaving_atoms']==['OXT']
    assert report['alpha_stereo']=='S'
