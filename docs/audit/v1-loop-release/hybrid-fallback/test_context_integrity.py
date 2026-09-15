import numpy as np
import pytest
from openmm import app

from context_integrity import (residue_key, outer_torsion_anchor_indices, affected_reference_keys,
                               outside_observed_torsions, AnchoredConstructionForceField)


@pytest.fixture
def protein():
    top = app.Topology()
    chain = top.addChain('Z')
    residues, named, previous = [], {}, None
    for i in range(1, 9):
        residue = top.addResidue('PRO' if i == 7 else 'ALA', chain, str(i + 100))
        residues.append(residue)
        names = {}
        for name in ['N', 'CA', 'C', 'O']:
            element = app.element.nitrogen if name == 'N' else app.element.oxygen if name == 'O' else app.element.carbon
            names[name] = top.addAtom(name, element, residue)
        for a, b in [('N', 'CA'), ('CA', 'C'), ('C', 'O')]:
            top.addBond(names[a], names[b])
        if previous:
            top.addBond(previous, names['N'])
        previous = names['C']
        named[i] = {name: atom.index for name, atom in names.items()}
    initial = np.arange(top.getNumAtoms() * 3, dtype=float).reshape(-1, 3)
    modeled = {residue_key(residues[i - 1]) for i in [4, 5]}
    flanks = {residue_key(residues[i - 1]) for i in [3, 6]}
    def row(i):
        return {'residue': list(residue_key(residues[i - 1])),
                'phi_atoms': [named[i - 1]['C'], named[i]['N'], named[i]['CA'], named[i]['C']],
                'psi_atoms': [named[i]['N'], named[i]['CA'], named[i]['C'], named[i + 1]['N']],
                'classification': 'outlier'}
    return top, residues, named, initial, modeled, flanks, row


def test_anchor_support_follows_connectivity_and_preserves_all_external_definitions(protein):
    top, residues, names, initial, modeled, flanks, row = protein
    anchors = outer_torsion_anchor_indices(top, modeled, flanks)
    assert anchors == {names[3]['N'], names[3]['CA'], names[6]['CA'], names[6]['C']}
    final = initial.copy()
    for atom in top.atoms():
        if residue_key(atom.residue) in modeled | flanks and atom.index not in anchors:
            final[atom.index] += [.2, -.1, .3]
    report = outside_observed_torsions(top, initial, final, modeled, flanks)
    assert report['definitions_checked'] > 0
    assert report['all_defining_coordinates_exactly_preserved']
    final[names[6]['C'], 0] += 1e-9
    report = outside_observed_torsions(top, initial, final, modeled, flanks)
    assert not report['all_defining_coordinates_exactly_preserved']
    assert {r['torsion'] for r in report['changed_definitions']} == {'phi', 'omega'}


def test_unchanged_external_outlier_is_excluded_but_new_boundary_is_always_checked(protein):
    top, residues, names, initial, modeled, flanks, row = protein
    report = {'rows': [row(i) for i in [3, 4, 5, 6, 7]], 'unavailable': []}
    selected = affected_reference_keys(top, initial, initial.copy(), modeled, report)
    assert selected == {residue_key(residues[i - 1]) for i in [3, 4, 5, 6]}


@pytest.mark.parametrize('atom_name', ['C', 'CA'])
def test_changed_external_phi_or_omega_cannot_hide_behind_preexisting_outlier(protein, atom_name):
    top, residues, names, initial, modeled, flanks, row = protein
    final = initial.copy()
    final[names[6][atom_name], 0] += 1e-9
    report = {'rows': [row(7)], 'unavailable': []}
    assert residue_key(residues[6]) in affected_reference_keys(top, initial, final, modeled, report)


def test_unavailable_reference_is_never_silently_removed(protein):
    top, residues, names, initial, modeled, flanks, row = protein
    key = residue_key(residues[6])
    report = {'rows': [], 'unavailable': [{'residue': key, 'reason': 'missing definition'}]}
    assert key in affected_reference_keys(top, initial, initial.copy(), modeled, report)


def test_fixed_masses_are_confined_to_fresh_construction_system():
    class System:
        def __init__(self):
            self.masses = [12., 14., 16.]
        def setParticleMass(self, index, value):
            self.masses[index] = value
    class ForceField:
        def createSystem(self, *args, **kwargs):
            return System()
    physical = ForceField()
    construction = AnchoredConstructionForceField(physical, {1})
    assert construction.createSystem().masses == [12., 0., 16.]
    assert physical.createSystem().masses == [12., 14., 16.]


def test_ambiguous_peptide_branch_is_rejected(protein):
    top, residues, names, initial, modeled, flanks, row = protein
    atoms = list(top.atoms())
    top.addBond(atoms[names[6]['C']], atoms[names[8]['N']])
    with pytest.raises(ValueError, match='Ambiguous peptide'):
        outer_torsion_anchor_indices(top, modeled, flanks)
