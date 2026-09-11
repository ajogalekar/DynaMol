"""Chemical identity and native force-field checks for the ligand prep path."""
import json
from pathlib import Path

import numpy as np
import pytest
from openmm import app, unit
from rdkit import Chem
from rdkit.Chem import AllChem

from backend import config, sources, storage
from backend.ligands import (LigandModel, _add_hydrogens, _graph_from_override,
                             _mapped_variant, _normalize_native_charge, _parameterize,
                             _select_state, inspect_ligands, ligand_runtime_status,
                             prepare_ligands)


@pytest.fixture
def isolated_data(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'DATA_ROOT', tmp_path)
    monkeypatch.setattr(config, 'DATASETS_DIR', tmp_path / 'datasets')
    monkeypatch.setattr(config, 'JOBS_DIR', tmp_path / 'jobs')
    return tmp_path


def embedded(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=2026) == 0
    return mol


def test_smiles_source_charge_and_bound_pose_preserved(isolated_data):
    meta = sources.create_smiles('CC(=O)O', name='acetate reference')
    inspected = inspect_ligands(meta['id'], ph=7)
    assert len(inspected) == 1
    assert inspected[0]['error'] is None
    assert inspected[0]['original_formal_charge'] == 0
    assert inspected[0]['formal_charge'] == -1
    assert inspected[0]['reference']['provider'].startswith('Retained authoritative')


def test_hydrogen_only_relaxation_preserves_bound_heavy_atoms():
    mol = Chem.RemoveHs(embedded('C[C@H](O)C(=O)O'))
    names = [f'{atom.GetSymbol()}{i + 1}' for i, atom in enumerate(mol.GetAtoms())]
    model = LigandModel(mol, None, list(range(mol.GetNumAtoms())), names, {})
    prepared, method = _add_hydrogens(model)
    np.testing.assert_allclose(prepared.GetConformer().GetPositions()[:mol.GetNumAtoms()], mol.GetConformer().GetPositions(), atol=1e-9)
    assert prepared.GetNumAtoms() > mol.GetNumAtoms()
    assert 'hydrogen-only' in method


def test_ph_heuristic_preserves_amide_charge_and_bound_stereotags():
    mol = Chem.RemoveHs(embedded('CC(=O)N[C@@H](C)CN(C)C'))
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    names = [f'{a.GetSymbol()}{i + 1}' for i, a in enumerate(mol.GetAtoms())]
    selected, report = _select_state(mol, names, 'UNKNOWN', 7, None)
    amide_n = mol.GetSubstructMatch(Chem.MolFromSmarts('[N]-[C](=O)'))[0]
    assert selected.GetAtomWithIdx(amide_n).GetFormalCharge() == 0
    assert Chem.GetFormalCharge(selected) == 1
    assert [a.GetChiralTag() for a in selected.GetAtoms()] == [a.GetChiralTag() for a in mol.GetAtoms()]
    np.testing.assert_array_equal(selected.GetConformer().GetPositions(), mol.GetConformer().GetPositions())
    assert 'no microscopic pKa' in report['protonation_method']


def test_opposite_stereo_override_is_rejected():
    mol = Chem.RemoveHs(embedded('C[C@H](O)CO'))
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    with pytest.raises(ValueError, match='stereochemistry'):
        _mapped_variant(mol, 'C[C@@H](O)CO')


def test_unknown_pdb_override_requires_unique_or_explicit_mapping():
    mol = Chem.RemoveHs(embedded('CCO'))
    top = app.Topology(); residue = top.addResidue('LIG', top.addChain('A'), '1')
    for i, atom in enumerate(mol.GetAtoms()):
        top.addAtom(f'{atom.GetSymbol()}{i + 1}', app.element.Element.getBySymbol(atom.GetSymbol()), residue)
    prepared, atoms, names, stereo, method = _graph_from_override(residue, mol.GetConformer().GetPositions(), 'CCO')
    assert Chem.MolToSmiles(prepared) == 'CCO'
    np.testing.assert_array_equal(prepared.GetConformer().GetPositions(), mol.GetConformer().GetPositions())
    assert 'Unique' in method
    benzene = Chem.RemoveHs(embedded('c1ccccc1'))
    top = app.Topology(); residue = top.addResidue('LIG', top.addChain('A'), '1')
    for i in range(6):
        top.addAtom(f'C{i + 1}', app.element.carbon, residue)
    with pytest.raises(ValueError, match='uniquely'):
        _graph_from_override(residue, benzene.GetConformer().GetPositions(), 'c1ccccc1')


def test_charge_normalization_rejects_unexplained_residual(tmp_path):
    (tmp_path / 'sqm.out').write_text('geometry converged\n Atom Element Mulliken Charge\n 1 C 0.000\n Total Mulliken Charge = 0.000\n')
    (tmp_path / 'charged.mol2').write_text('@<TRIPOS>ATOM\n1 C1 0 0 0 c3 1 LIG 0.004\n@<TRIPOS>BOND\n')
    with pytest.raises(ValueError, match='cannot be attributed'):
        _normalize_native_charge(tmp_path, 0)


@pytest.mark.skipif(not ligand_runtime_status()['available'], reason='private AmberTools runtime not installed')
def test_native_gaff2_xml_equivalence_and_combined_templates(tmp_path):
    import openmm as mm
    all_models = []
    xmls = []
    for label, smiles in [('ethanol', 'CCO'), ('acetate', 'CC(=O)[O-]')]:
        mol = embedded(smiles)
        path = tmp_path / label
        xml, names, charges = _parameterize(mol, path, 'DML_' + label, None, None)
        assert abs(sum(charges) - Chem.GetFormalCharge(mol)) < 1e-6
        assert np.ptp(charges) > .01
        amber = app.AmberPrmtopFile(str(path / 'ligand.prmtop'))
        positions = app.AmberInpcrdFile(str(path / 'ligand.inpcrd')).positions
        native = amber.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
        converted = app.ForceField(str(xml)).createSystem(amber.topology, nonbondedMethod=app.NoCutoff, constraints=None)
        energies = []
        for system in [native, converted]:
            integrator = mm.VerletIntegrator(.001)
            context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
            context.setPositions(positions)
            energies.append(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            del context, integrator
        assert np.isfinite(energies).all()
        assert abs(energies[0] - energies[1]) < 1e-4
        xmls.append(str(xml)); all_models.append((amber.topology, positions))
        assert (path / 'sqm.out').is_file()
        assert (path / 'charged-native.mol2').is_file()
    combined = app.Modeller(*all_models[0])
    combined.add(*all_models[1])
    ff = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', *xmls)
    system = ff.createSystem(combined.topology, nonbondedMethod=app.NoCutoff)
    assert system.getNumParticles() == combined.topology.getNumAtoms()


@pytest.mark.skipif(not ligand_runtime_status()['available'], reason='private AmberTools runtime not installed')
def test_native_prepared_atoms_survive_openmm_hydrogen_assignment(isolated_data):
    meta = sources.create_smiles('CCO', name='native ethanol')
    original = storage.load_physical(meta['id'])
    prepared = prepare_ligands(meta['id'], 7, 2026, isolated_data / 'artifacts')[0]
    modeller = app.Modeller(prepared.topology, prepared.positions)
    ff = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', str(prepared.ffxml))
    before = np.asarray(modeller.positions.value_in_unit(unit.angstrom)).copy()
    modeller.addHydrogens(ff, pH=7)
    np.testing.assert_allclose(np.asarray(modeller.positions.value_in_unit(unit.angstrom)), before, atol=1e-6, rtol=0)
    for atom_map in prepared.atom_map:
        np.testing.assert_allclose(before[atom_map['ligand_index']], original.xyz[0, atom_map['original_index']] * 10, atol=.001)
    assert prepared.original_residue_key == ('A', '1', '', 'LIG')


def test_ez_override_mismatch_rejected():
    mol = Chem.RemoveHs(embedded('C/C=C/Cl'))
    with pytest.raises(ValueError, match='E/Z'):
        _mapped_variant(mol, 'C/C=C\\Cl')


def test_ambiguous_state_override_rejected():
    mol = Chem.RemoveHs(embedded('c1ccccc1'))
    with pytest.raises(ValueError, match='uniquely'):
        _mapped_variant(mol, 'c1ccccc1')


def test_list_valued_legacy_ancestry_inspects(isolated_data):
    meta = sources.create_smiles('CCO')
    folder = storage.dataset_dir(meta['id'])
    saved = json.loads((folder / 'metadata.json').read_text())
    saved['preparation'] = ['legacy description']
    storage.atomic_json(folder / 'metadata.json', saved)
    assert inspect_ligands(meta['id'])[0]['error'] is None


def test_original_pdb_link_blocks_covalent_ligand(isolated_data):
    meta = sources.create_smiles('CCO')
    folder = storage.dataset_dir(meta['id'])
    link = list(' ' * 80)
    link[:6] = 'LINK  '
    link[17:20] = 'LIG'
    link[47:50] = 'CYS'
    (folder / 'source.pdb').write_text(''.join(link) + '\n')
    record = inspect_ligands(meta['id'])[0]
    assert 'covalent' in record['error']


@pytest.mark.skipif(not ligand_runtime_status()['available'], reason='private AmberTools runtime not installed')
def test_identical_copies_keep_parameter_source_and_individual_provenance(isolated_data):
    import hashlib
    import mdtraj as md
    meta = sources.create_smiles('CCO', name='ethanol source')
    source = storage.load_physical(meta['id'])
    chemical = json.loads((storage.dataset_dir(meta['id']) / 'chemistry.json').read_text())
    top = source.topology.copy()
    second = top.add_residue('LIG', top.add_chain('B'), resSeq=2)
    count = source.n_atoms
    copies = [top.add_atom(a.name, a.element, second) for a in source.topology.atoms]
    for bond in source.topology.bonds:
        top.add_bond(copies[bond.atom1.index], copies[bond.atom2.index])
    trajectory = md.Trajectory(np.concatenate([source.xyz, source.xyz + np.array([2., 0., 0.])], axis=1), top)
    paired = storage.save_dataset(trajectory, 'Two ethanol copies', 'test', 'Provenance reuse regression')
    paired_chemistry = {'schema_version': 1, 'authoritative_bonds': True,
                        'atoms': chemical['atoms'] + [{**a, 'index': a['index'] + count} for a in chemical['atoms']],
                        'bonds': chemical['bonds'] + [{**b, 'atoms': [i + count for i in b['atoms']]} for b in chemical['bonds']]}
    storage.atomic_json(storage.dataset_dir(paired['id']) / 'chemistry.json', paired_chemistry)
    results = prepare_ligands(paired['id'], 7, 2026, isolated_data / 'paired-parameters')
    assert len(results) == 2 and results[0].ffxml == results[1].ffxml
    folder = results[0].ffxml.parent
    source_annotation = json.loads((folder / 'provenance.json').read_text())
    assert source_annotation['key'] == results[0].provenance['key'] == 'A:1::LIG'
    assert source_annotation['parameterization_source_residue_key'] == ['A', '1', '', 'LIG']
    for result in results:
        instance = hashlib.sha256(result.provenance['key'].encode()).hexdigest()[:12]
        annotation = json.loads((folder / f'residue-{instance}.json').read_text())
        assert annotation['key'] == result.provenance['key']
        assert annotation['parameterization_source_residue_key'] == ['A', '1', '', 'LIG']
    native = app.AmberInpcrdFile(str(folder / 'ligand.inpcrd')).positions.value_in_unit(unit.angstrom)
    first = np.asarray(results[0].positions.value_in_unit(unit.angstrom))
    second = np.asarray(results[1].positions.value_in_unit(unit.angstrom))
    np.testing.assert_allclose(native, first, atol=.002, rtol=0)
    assert np.max(np.abs(np.asarray(native) - second)) > 10
    commands = json.loads((folder / 'commands.json').read_text())
    assert sum(Path(command['argv'][0]).name == 'antechamber' for command in commands) == 1
