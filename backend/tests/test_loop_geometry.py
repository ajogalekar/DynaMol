"""Analytical controls for all-loop geometry, independent of stereo correction."""
import json

import numpy as np
import pytest
from openmm import app, unit

from backend import config
from backend.loop_geometry import loop_geometry_report
from backend.residue_identity import residue_key


@pytest.fixture
def peptide():
    pdb = app.PDBFile(str(config.ROOT / 'docs/audit/preparation-fixtures/six_residues_intact.pdb'))
    return pdb.topology, np.asarray(pdb.positions.value_in_unit(unit.nanometer)).copy()


def selected(topology):
    return {residue_key(list(topology.residues())[2])}


def named(topology, resid=3):
    return {a.name: a.index for a in topology.atoms() if a.residue.id == str(resid)}


def copy_without(topology, xyz, *, atom_to_remove=None, bond_to_remove=None):
    result = app.Topology()
    mapping = {}
    positions = []
    for chain in topology.chains():
        added_chain = result.addChain(chain.id)
        for residue in chain.residues():
            added_residue = result.addResidue(residue.name, added_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                if atom.index == atom_to_remove:
                    continue
                mapping[atom] = result.addAtom(atom.name, atom.element, added_residue)
                positions.append(xyz[atom.index])
    for a, b in topology.bonds():
        if a in mapping and b in mapping and {a.index, b.index} != bond_to_remove:
            result.addBond(mapping[a], mapping[b])
    return result, np.asarray(positions)


def test_intact_experimental_fragment_checks_both_anchors_and_preserves_inputs(peptide):
    topology, xyz = peptide
    original = xyz.copy()
    identity = [(a.name, residue_key(a.residue)) for a in topology.atoms()]
    report = loop_geometry_report(topology, xyz * unit.nanometer, selected(topology))
    assert report['accepted'], report['errors']
    assert report['modeled_heavy_atom_count'] == 8
    assert len(report['omega_checks']) == 2
    assert sum(row['kind'] == 'peptide_C-N' for row in report['bond_checks']) == 2
    assert len(report['angle_checks']) == 8
    assert len(report['phi_psi']) == 1
    assert all(row['state'] == 'trans' for row in report['omega_checks'])
    np.testing.assert_array_equal(xyz, original)
    assert identity == [(a.name, residue_key(a.residue)) for a in topology.atoms()]
    json.dumps(report, allow_nan=False)


def test_distorted_sidechain_rejected_even_without_a_stereo_correction(peptide):
    topology, xyz = peptide
    atoms = named(topology)
    xyz[atoms['CD1']] += [1, 0, 0]
    report = loop_geometry_report(topology, xyz, selected(topology))
    assert not report['accepted']
    assert any(row['kind'] == 'other_heavy' and not row['accepted'] for row in report['bond_checks'])


def test_backbone_bond_length_is_checked_independently_of_collision(peptide):
    topology, xyz = peptide
    atoms = named(topology)
    xyz[atoms['O']] = xyz[atoms['C']] + 2 * (xyz[atoms['O']] - xyz[atoms['C']])
    report = loop_geometry_report(topology, xyz, selected(topology))
    assert not report['accepted']
    assert any(set(row['atoms']) == {atoms['C'], atoms['O']} and not row['accepted'] for row in report['bond_checks'])


def test_missing_loop_backbone_atom_rejected(peptide):
    topology, xyz = peptide
    keys = selected(topology)
    topology, xyz = copy_without(topology, xyz, atom_to_remove=named(topology)['O'])
    report = loop_geometry_report(topology, xyz, keys)
    assert not report['accepted']
    assert any('lacks required backbone atoms: O' in message for message in report['errors'])


def test_geometrically_closed_but_unbonded_peptide_is_rejected(peptide):
    topology, xyz = peptide
    keys = selected(topology)
    pair = {named(topology, 2)['C'], named(topology, 3)['N']}
    topology, xyz = copy_without(topology, xyz, bond_to_remove=pair)
    report = loop_geometry_report(topology, xyz, keys)
    assert not report['accepted']
    assert any('Missing expected peptide bond' in message for message in report['errors'])
    assert all(row['accepted'] for row in report['omega_checks'])


def test_disconnected_sidechain_is_rejected_even_at_original_coordinates(peptide):
    topology, xyz = peptide
    keys = selected(topology)
    atoms = named(topology)
    topology, xyz = copy_without(topology, xyz, bond_to_remove={atoms['CG1'], atoms['CD1']})
    report = loop_geometry_report(topology, xyz, keys)
    assert not report['accepted']
    assert any('disconnected from its backbone' in message for message in report['errors'])


def test_backbone_angle_rejection_does_not_require_bond_stretch(peptide):
    topology, xyz = peptide
    atoms = named(topology)
    ca = xyz[atoms['CA']].copy()
    n_direction = xyz[atoms['N']] - ca
    xyz[atoms['C']] = ca - n_direction / np.linalg.norm(n_direction) * np.linalg.norm(xyz[atoms['C']] - ca)
    report = loop_geometry_report(topology, xyz, selected(topology))
    assert not report['accepted']
    assert any(row['kind'] == 'N-CA-C' and not row['accepted'] for row in report['angle_checks'])


def test_nonplanar_peptide_detected_with_unchanged_bond_lengths_and_angles(peptide):
    topology, xyz = peptide
    before = named(topology, 2)
    after = named(topology, 3)
    pivot = xyz[after['N']].copy()
    axis = pivot - xyz[before['C']]
    axis /= np.linalg.norm(axis)
    indices = [a.index for a in topology.atoms() if a.residue.index >= 2]
    vectors = xyz[indices] - pivot
    # A 90-degree rigid rotation of the downstream chain around its peptide bond.
    xyz[indices] = pivot + np.cross(axis, vectors) + np.outer(vectors @ axis, axis)
    report = loop_geometry_report(topology, xyz, selected(topology))
    assert not report['accepted']
    assert all(row['accepted'] for row in report['bond_checks'])
    assert all(row['accepted'] for row in report['angle_checks'])
    assert not report['omega_checks'][0]['accepted']


def test_planar_nonproline_cis_is_recorded_as_uncertain(peptide):
    topology, xyz = peptide
    before, after = named(topology, 2), named(topology, 3)
    pivot = xyz[after['N']].copy()
    axis = pivot - xyz[before['C']]
    axis /= np.linalg.norm(axis)
    indices = [a.index for a in topology.atoms() if a.residue.index >= 2]
    vectors = xyz[indices] - pivot
    xyz[indices] = pivot - vectors + 2 * np.outer(vectors @ axis, axis)
    report = loop_geometry_report(topology, xyz, selected(topology))
    assert report['omega_checks'][0]['accepted']
    assert report['omega_checks'][0]['state'] == 'cis'
    assert any('non-proline cis' in message for message in report['warnings'])


def test_external_ligand_collision_uses_every_modeled_heavy_atom(peptide):
    topology, xyz = peptide
    atoms = named(topology)
    # Backbone was never a target of sidechain correction; it still needs checking.
    report = loop_geometry_report(topology, xyz, selected(topology), [(xyz[atoms['CA']].copy(), .17)])
    assert not report['accepted']
    assert any(row['distance_nm'] == 0 for row in report['gross_collisions'])
    assert report['external_environment_atoms'] == 1


def test_soft_overlap_is_reported_without_becoming_a_gross_collision(peptide):
    topology, xyz = peptide
    atoms = named(topology)
    direction = xyz[atoms['CD1']] - xyz[atoms['CG1']]
    direction /= np.linalg.norm(direction)
    point = xyz[atoms['CD1']] + .20 * direction
    report = loop_geometry_report(topology, xyz, selected(topology), [(point, .17)])
    assert report['accepted'], report['errors']
    assert report['requires_minimization']
    assert any(row['external_environment_index'] == 0 for row in report['soft_overlaps'])
    assert report['gross_collisions'] == []


def test_direct_bonds_and_one_three_pairs_are_excluded_from_contacts(peptide):
    topology, xyz = peptide
    report = loop_geometry_report(topology, xyz, selected(topology))
    excluded = set()
    adjacency = [set() for _ in topology.atoms()]
    for a, b in topology.bonds():
        adjacency[a.index].add(b.index)
        adjacency[b.index].add(a.index)
    for i, neighbors in enumerate(adjacency):
        for j in neighbors | set().union(*(adjacency[n] for n in neighbors)):
            excluded.add(tuple(sorted((i, j))))
    assert all(tuple(row['atoms']) not in excluded for row in report['soft_overlaps'])


@pytest.mark.parametrize('value', [np.nan, np.inf, -np.inf])
def test_nonfinite_coordinates_return_json_safe_rejection(peptide, value):
    topology, xyz = peptide
    xyz[0, 0] = value
    report = loop_geometry_report(topology, xyz, selected(topology))
    assert not report['accepted']
    assert 'finite' in report['errors'][0]
    json.dumps(report, allow_nan=False)


def test_insertion_codes_are_required_for_exact_residue_selection(peptide):
    topology, xyz = peptide
    residue = list(topology.residues())[2]
    residue.insertionCode = 'A'
    wrong = loop_geometry_report(topology, xyz, {('A', '3', '', 'ILE')})
    correct = loop_geometry_report(topology, xyz, {('A', '3', 'A', 'ILE')})
    assert not wrong['accepted']
    assert correct['accepted'], correct['errors']


def test_empty_modeled_residue_returns_rejection_instead_of_tree_error():
    topology = app.Topology()
    residue = topology.addResidue('GLY', topology.addChain('A'), '3')
    report = loop_geometry_report(topology, np.empty((0, 3)), {residue_key(residue)})
    assert not report['accepted']
    assert 'no heavy atoms' in report['errors'][0]


@pytest.mark.parametrize('external', [[([0, 0, np.nan], .17)], [([0, 0, 0], 0)], [([0, 0, 0], np.inf)],
                                    [([0, 0, 0, 0, 0, 0], .17)], [([0, 0, 0], [.17, .17])]])
def test_invalid_external_environment_is_not_silently_ignored(peptide, external):
    topology, xyz = peptide
    report = loop_geometry_report(topology, xyz, selected(topology), external)
    assert not report['accepted']
    json.dumps(report, allow_nan=False)
