"""Sulfate retains a declared source state; neutralization is never inferred."""
import numpy as np
import pytest
from rdkit import Chem

from backend.ligands import LigandModel, _add_hydrogens, _select_state, _unused_gaff14_defaults


def sulfate():
    molecule = Chem.MolFromSmiles('[O-]S(=O)(=O)[O-]')
    conformer = Chem.Conformer(5)
    for index, position in enumerate(((1.4, 0, 0), (0, 0, 0), (-.5, 1.3, 0),
                                       (-.5, -.6, 1.2), (-.5, -.6, -1.2))):
        conformer.SetAtomPosition(index, position)
    molecule.AddConformer(conformer)
    return molecule


@pytest.mark.parametrize('ph', [6, 7, 8])
def test_sulfate_reference_state_keeps_identity_charge_bonds_and_coordinates(monkeypatch, ph):
    import dimorphite_dl
    monkeypatch.setattr(dimorphite_dl, 'protonate_smiles',
                        lambda *a, **k: pytest.fail('Do not neutralize the named sulfate reference.'))
    original = sulfate()
    before = Chem.MolToMolBlock(original)
    selected, report = _select_state(original, ['O1', 'S', 'O2', 'O3', 'O4'], 'SO4', ph, None)
    assert selected is not original
    assert Chem.GetFormalCharge(selected) == -2
    assert [a.GetFormalCharge() for a in selected.GetAtoms()] == [-1, 0, 0, 0, -1]
    assert sum(a.GetTotalNumHs() for a in selected.GetAtoms()) == 0
    assert Chem.MolToMolBlock(selected) == before == Chem.MolToMolBlock(original)
    np.testing.assert_array_equal(selected.GetConformer().GetPositions(), original.GetConformer().GetPositions())
    assert report['reference_state']['original_graph_retained']
    assert report['reference_state']['ph_heuristic_applied'] is False
    assert 'fixed reference-state' in report['protonation_method']
    assert report['candidate_smiles'] == [report['original_smiles']]


@pytest.mark.parametrize('ph', [5.99, 8.01])
def test_sulfate_outside_declared_range_needs_an_explicit_state(ph):
    with pytest.raises(ValueError, match='explicit atom-mapped'):
        _select_state(sulfate(), ['O1', 'S', 'O2', 'O3', 'O4'], 'SO4', ph, None)


@pytest.mark.parametrize('smiles', ['OS(=O)(=O)O', 'OS(=O)(=O)[O-]',
                                     'COS(=O)(=O)[O-]', '[O-]P(=O)([O-])[O-]'])
def test_component_label_does_not_authorize_changing_a_different_graph(smiles):
    original = Chem.MolFromSmiles(smiles)
    before = Chem.MolToSmiles(original)
    with pytest.raises(ValueError, match='does not match'):
        _select_state(original, [str(i) for i in range(original.GetNumAtoms())], 'SO4', 7, None)
    assert Chem.MolToSmiles(original) == before


def test_explicit_mapped_state_takes_precedence_over_reference_ph_policy():
    original = sulfate()
    selected, report = _select_state(original, ['O1', 'S', 'O2', 'O3', 'O4'], 'SO4', 2,
                                    '[O-:1][S:2](=[O:3])(=[O:4])[O-:5]')
    assert Chem.GetFormalCharge(selected) == -2
    assert report['protonation_method'].startswith('Explicit user')
    assert 'reference_state' not in report


def test_unregistered_component_does_not_receive_a_sulfate_fallback(monkeypatch):
    import dimorphite_dl
    def fail(*args, **kwargs):
        raise RuntimeError('ordinary heuristic path exercised')
    monkeypatch.setattr(dimorphite_dl, 'protonate_smiles', fail)
    with pytest.raises(RuntimeError, match='ordinary heuristic'):
        _select_state(sulfate(), ['O1', 'S', 'O2', 'O3', 'O4'], 'UNKNOWN', 7, None)


@pytest.mark.parametrize('smiles', ['[O-]S(=O)(=O)[O-]', 'ClC(Cl)(Cl)Cl'])
def test_no_added_hydrogens_skips_all_fixed_minimization(monkeypatch, smiles):
    from rdkit.Chem import AllChem
    original = Chem.MolFromSmiles(smiles)
    original.AddConformer(sulfate().GetConformer())
    names = [f'{atom.GetSymbol()}{i + 1}' for i, atom in enumerate(original.GetAtoms())]
    before = Chem.MolToMolBlock(original)
    def fail(*args, **kwargs):
        pytest.fail('No optimizer is needed when there are no new hydrogens.')
    monkeypatch.setattr(AllChem, 'MMFFHasAllMoleculeParams', fail)
    monkeypatch.setattr(AllChem, 'UFFHasAllMoleculeParams', fail)
    model = LigandModel(original, None, list(range(5)), names, {})
    selected, method = _add_hydrogens(model)
    assert selected.GetNumAtoms() == 5
    assert Chem.MolToMolBlock(selected) == before == Chem.MolToMolBlock(original)
    np.testing.assert_array_equal(selected.GetConformer().GetPositions(), original.GetConformer().GetPositions())
    assert [a.GetProp('_TriposAtomName') for a in selected.GetAtoms()] == names
    assert 'No hydrogens added' in method


def graph_structure(edges, count):
    from types import SimpleNamespace
    atoms = [SimpleNamespace(idx=index) for index in range(count)]
    return SimpleNamespace(atoms=atoms, bonds=[SimpleNamespace(atom1=atoms[a], atom2=atoms[b]) for a, b in edges])


@pytest.mark.parametrize('edges,count', [([(0, 1), (0, 2), (0, 3), (0, 4)], 5),
                                        ([(0, 1), (1, 2), (2, 3), (3, 0)], 4)])
def test_unused_global_14_defaults_are_compatible_without_changing_pairs(edges, count):
    from parmed.parameters import ParameterSet
    params = ParameterSet()
    report = _unused_gaff14_defaults(graph_structure(edges, count), params)
    assert params.default_scee == 1.2 and params.default_scnb == 2
    assert report['shortest_distance_three_pairs'] == 0
    assert report['native_energy_force_validation_required']
    assert report['original_generic_defaults'] == {'scee': 1, 'scnb': 1}


def test_graph_14_pairs_prevent_normalization_even_without_torsion_terms():
    from parmed.parameters import ParameterSet
    params = ParameterSet()
    assert _unused_gaff14_defaults(graph_structure([(0, 1), (1, 2), (2, 3)], 4), params) is None
    assert params.default_scee == params.default_scnb == 1


def test_existing_torsion_parameters_are_never_rescaled():
    from parmed.parameters import ParameterSet
    params = ParameterSet()
    params.dihedral_types[('a', 'b', 'c', 'd')] = object()
    assert _unused_gaff14_defaults(graph_structure([(0, 1)], 2), params) is None
    assert params.default_scee == params.default_scnb == 1
