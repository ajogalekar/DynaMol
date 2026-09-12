"""Protein identity is retained independently of residue parameter availability."""
import hashlib

import mdtraj as md
import numpy as np
import pytest
from openmm import app, unit

from backend import config, ligands, preparation, storage
from backend.residue_identity import protein_residue_keys, residue_key


def peptide(names=('ALA', 'ZZZ', 'GLY'), connected=True):
    topology = app.Topology()
    chain = topology.addChain('A')
    positions = []
    previous = None
    for i, name in enumerate(names):
        residue = topology.addResidue(name, chain, str(i + 1))
        atoms = {}
        for atom_name, element, offset in [('N', app.element.nitrogen, 0), ('CA', app.element.carbon, .08), ('C', app.element.carbon, .17), ('O', app.element.oxygen, .21)]:
            atoms[atom_name] = topology.addAtom(atom_name, element, residue)
            positions.append([i * .3 + offset, .05 if atom_name == 'O' else 0, 0])
        topology.addBond(atoms['N'], atoms['CA'])
        topology.addBond(atoms['CA'], atoms['C'])
        topology.addBond(atoms['C'], atoms['O'])
        if previous and connected:
            topology.addBond(previous, atoms['N'])
        previous = atoms['C']
    return topology, np.array(positions) * unit.nanometer


@pytest.fixture(autouse=True)
def isolated_data(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'DATA_ROOT', tmp_path)
    monkeypatch.setattr(config, 'DATASETS_DIR', tmp_path / 'datasets')
    monkeypatch.setattr(config, 'JOBS_DIR', tmp_path / 'jobs')
    config.DATASETS_DIR.mkdir()
    config.JOBS_DIR.mkdir()


def test_unknown_amino_acid_is_polymer_when_linked_to_protein():
    top, positions = peptide()
    unknown = list(top.residues())[1]
    assert residue_key(unknown) in protein_residue_keys(None, top)
    top_without_bonds, positions = peptide(connected=False)
    unknown = list(top_without_bonds.residues())[1]
    assert residue_key(unknown) not in protein_residue_keys(None, top_without_bonds)
    assert residue_key(unknown) in protein_residue_keys(None, top_without_bonds, positions)


def test_isolated_unknown_amino_acid_shaped_ligand_is_not_swept_into_protein():
    top, positions = peptide(('ZZZ',))
    assert protein_residue_keys(None, top, positions) == set()


@pytest.mark.parametrize('as_mdtraj', [False, True])
def test_monatomic_calcium_alias_is_not_a_protein_residue(as_mdtraj):
    top = app.Topology()
    residue = top.addResidue('CAL', top.addChain('I'), '1')
    top.addAtom('CA', app.element.calcium, residue)
    if as_mdtraj:
        top = md.Topology.from_openmm(top)
    assert protein_residue_keys(None, top) == set()


def test_multiatom_cal_amino_acid_identity_is_not_removed_as_an_ion():
    top, positions = peptide(('ALA', 'CAL', 'GLY'))
    residue = list(top.residues())[1]
    assert residue_key(residue) in protein_residue_keys(None, top, positions)


def test_backbone_inspection_does_not_skip_phosphorylated_residue():
    top, positions = peptide(('ALA', 'TPO', 'GLY'))
    assert preparation.backbone_gaps(top, positions) == []
    moved = np.array(positions.value_in_unit(unit.nanometer))
    moved[8:] += 2
    gaps = preparation.backbone_gaps(top, moved * unit.nanometer)
    assert len(gaps) == 1
    assert (gaps[0]['after'], gaps[0]['before']) == ('2', '3')


@pytest.mark.parametrize('sequence,type_,expected', [('1', 'L-peptide linking', True), ('.', 'L-peptide linking', False), ('1', 'NON-POLYMER', False)])
def test_mmcif_polymer_evidence_recognizes_incomplete_uaa_without_claiming_all_heterogens(tmp_path, monkeypatch, sequence, type_, expected):
    top = app.Topology()
    residue = top.addResidue('ZZZ', top.addChain('A'), '1')
    top.addAtom('CA', app.element.carbon, residue)
    source = tmp_path / 'source.cif'
    source.write_text(f'''data_test
loop_
_chem_comp.id
_chem_comp.type
ZZZ '{type_}'
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.label_seq_id
1 C CA CA ZZZ A A 1 {sequence}
''')
    monkeypatch.setattr(ligands, '_source_cif', lambda _: source)
    assert (residue_key(residue) in protein_residue_keys('fixture', top)) is expected


@pytest.mark.parametrize('missing', ['.', '?'])
@pytest.mark.parametrize('polymer_type,expected', [("'polypeptide(L)'", True), ('.', False), ('?', False)])
def test_nullable_cif_types_preserve_available_polymer_evidence(tmp_path, monkeypatch, missing, polymer_type, expected):
    top = app.Topology()
    residue = top.addResidue('ZZZ', top.addChain('A'), '1')
    top.addAtom('CA', app.element.carbon, residue)
    source = tmp_path / 'source.cif'
    source.write_text(f'''data_nullable
loop_
_chem_comp.id
_chem_comp.type
ZZZ {missing}
loop_
_entity_poly.entity_id
_entity_poly.type
1 {polymer_type}
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.auth_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.label_seq_id
_atom_site.label_entity_id
1 C CA ZZZ A A 1 1 1
''')
    monkeypatch.setattr(ligands, '_source_cif', lambda _: source)
    assert (residue_key(residue) in protein_residue_keys('fixture', top)) is expected


@pytest.mark.parametrize('name', ['ZZZ', 'MSE', 'ALY', 'CSO', 'FME', 'DLY'])
def test_unsupported_uaa_is_reported_per_residue_and_never_sent_to_ligand_params(name):
    top, positions = peptide(('ALA', name, 'GLY'))
    traj = md.Trajectory(np.array(positions.value_in_unit(unit.nanometer))[None], md.Topology.from_openmm(top))
    meta = storage.save_dataset(traj, 'Unknown peptide', 'Analytical fixture', 'Preservation check')
    assert all(atom['category'] == 'protein' for atom in meta['atoms'])
    path = storage.dataset_dir(meta['id']) / 'topology.pdb'
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    inspection = preparation.inspect_preparation(meta['id'])
    assert inspection['ligands'] == []
    record = next(row for row in inspection['modified_residues'] if row['residue'] == name)
    assert record['supported'] is False
    assert f'{name} A:2' in record['error']
    assert 'covalent amino-acid template' in record['error']
    for remove in [False, True]:
        with pytest.raises(ValueError, match='identity and atoms are preserved'):
            preparation.submit_preparation({'dataset_id': meta['id'], 'remove_heterogens': remove})
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
    assert list(config.JOBS_DIR.iterdir()) == []
