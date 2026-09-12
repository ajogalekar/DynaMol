"""Analytical controls for narrow, low-confidence new-loop stereo correction."""
import numpy as np
import pytest
from openmm import app, unit
from pdbfixer import PDBFixer

from backend import config
from backend.preparation_worker import atom_key, restore_modeled_loop_sidechains, stereochemistry_report


@pytest.fixture(scope='module')
def templates():
    return PDBFixer(filename=str(config.ROOT / 'docs/audit/preparation-fixtures/six_residues_intact.pdb')).templates


def inverted_template(templates, name):
    template = templates[name]
    topology = template.topology
    xyz = np.asarray(template.positions.value_in_unit(unit.nanometer)).copy()
    named = {a.name: a.index for a in topology.atoms()}
    center = xyz[named['CA']]
    normal = np.cross(xyz[named['N']] - center, xyz[named['C']] - center)
    normal /= np.linalg.norm(normal)
    for atom in topology.atoms():
        if atom.name not in {'N', 'CA', 'C', 'O', 'OXT'}:
            xyz[atom.index] -= 2 * np.dot(xyz[atom.index] - center, normal) * normal
    return topology, xyz


@pytest.mark.parametrize('name', ['ILE', 'THR', 'ARG'])
def test_wholly_modeled_sidechain_corrects_all_centers_without_moving_backbone(templates, name):
    topology, xyz = inverted_template(templates, name)
    original = xyz.copy()
    before = stereochemistry_report(topology, xyz * unit.nanometer, templates)
    assert before['violations']
    corrected, report = restore_modeled_loop_sidechains(topology, xyz * unit.nanometer, templates, set())
    actual = np.asarray(corrected.value_in_unit(unit.nanometer))
    assert report['accepted'], report['errors']
    assert not report['stereochemistry_after']['violations']
    assert report['corrections'][0]['rotation_determinant'] == pytest.approx(1)
    backbone = [a.index for a in topology.atoms() if a.name in {'N', 'CA', 'C', 'O', 'OXT'}]
    assert np.array_equal(actual[backbone], original[backbone])
    assert np.array_equal(xyz, original)
    for bond in report['bond_checks']:
        assert bond['distance_nm'] == pytest.approx(bond['template_distance_nm'], abs=1e-8)
    if name in {'ILE', 'THR'}:
        assert {v['center'] for v in before['violations']} == {'CA', 'CB'}


def test_even_one_observed_atom_protects_entire_residue(templates):
    topology, xyz = inverted_template(templates, 'ILE')
    observed = {atom_key(next(topology.atoms()))}
    _, report = restore_modeled_loop_sidechains(topology, xyz * unit.nanometer, templates, observed)
    assert not report['accepted']
    assert 'observed residue' in report['errors'][0]
    assert report['corrections'] == []


def test_proline_is_not_rebuilt_by_a_sidechain_shortcut(templates):
    topology, xyz = inverted_template(templates, 'PRO')
    _, report = restore_modeled_loop_sidechains(topology, xyz * unit.nanometer, templates, set())
    assert not report['accepted']
    assert 'non-Pro' in report['errors'][0]


def test_degenerate_backbone_frame_rejected(templates):
    topology, xyz = inverted_template(templates, 'ARG')
    names = {a.name: a.index for a in topology.atoms()}
    xyz[names['N']] = xyz[names['CA']] + .5 * (xyz[names['C']] - xyz[names['CA']])
    _, report = restore_modeled_loop_sidechains(topology, xyz * unit.nanometer, templates, set())
    assert not report['accepted']
    assert 'collinear' in report['errors'][0]


def test_retained_external_molecule_collision_rejects_candidate(templates):
    topology, xyz = inverted_template(templates, 'ALA')
    reference = np.asarray(templates['ALA'].positions.value_in_unit(unit.nanometer))
    cb = next(a.index for a in topology.atoms() if a.name == 'CB')
    _, report = restore_modeled_loop_sidechains(topology, xyz * unit.nanometer, templates, set(), [(reference[cb], .17)])
    assert not report['accepted']
    assert report['steric_screen']['gross_collisions']
    assert report['steric_screen']['external_environment_atoms'] == 1


def test_soft_overlap_is_reported_as_requiring_minimization(templates):
    topology, xyz = inverted_template(templates, 'ALA')
    reference = np.asarray(templates['ALA'].positions.value_in_unit(unit.nanometer))
    cb = next(a.index for a in topology.atoms() if a.name == 'CB')
    _, report = restore_modeled_loop_sidechains(topology, xyz * unit.nanometer, templates, set(), [(reference[cb] + [0, 0, .2], .17)])
    assert report['accepted'], report['errors']
    assert report['steric_screen']['requires_forcefield_minimization']
    assert report['steric_screen']['soft_overlap_score'] > 0


def test_insertion_code_keeps_distinct_same_name_residues_separate(templates):
    source = templates['ALA']; topology = app.Topology(); chain = topology.addChain('A')
    blocks = []; observed = set()
    _, inverted = inverted_template(templates, 'ALA')
    for code, coordinates in [('A', np.asarray(source.positions.value_in_unit(unit.nanometer))), ('B', inverted + [3, 0, 0])]:
        residue = topology.addResidue('ALA', chain, '1', insertionCode=code)
        mapping = {a: topology.addAtom(a.name, a.element, residue) for a in source.topology.atoms()}
        for a, b in source.topology.bonds(): topology.addBond(mapping[a], mapping[b])
        if code == 'A': observed.update(atom_key(a) for a in residue.atoms())
        blocks.append(coordinates)
    xyz = np.concatenate(blocks)
    corrected, report = restore_modeled_loop_sidechains(topology, xyz * unit.nanometer, templates, observed)
    # Deliberately separated same-chain residues also exercise the unchanged
    # backbone gap guard; identity selection must still target only code B.
    assert report['errors'] == ['A long backbone connection remains after loop correction.']
    assert [r['insertion_code'] for r in report['corrections']] == ['B']
    assert np.array_equal(np.asarray(corrected.value_in_unit(unit.nanometer))[:len(blocks[0])], blocks[0])
