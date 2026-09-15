"""Added branches may recover known standard stereo; observed coordinates may not."""
import numpy as np
import pytest
from openmm import app, unit
from pdbfixer import PDBFixer

from backend import config
from backend.added_atom_stereo import correct_added_branches, corrected_branch_contacts
from backend.preparation_worker import atom_key, stereochemistry_report


@pytest.fixture
def templates():
    return PDBFixer(filename=str(config.ROOT / 'docs/audit/preparation-fixtures/six_residues_intact.pdb')).templates


def invert(template, branch, plane):
    top = template.topology
    xyz = np.asarray(template.positions.value_in_unit(unit.nanometer)).copy()
    named = {a.name: a.index for a in top.atoms()}
    center, left, right = [named[n] for n in plane]
    normal = np.cross(xyz[left] - xyz[center], xyz[right] - xyz[center])
    normal /= np.linalg.norm(normal)
    index = named[branch]
    xyz[index] -= 2 * np.dot(xyz[index] - xyz[center], normal) * normal
    return top, xyz


@pytest.mark.parametrize('residue,other', [('THR', 'OG1'), ('ILE', 'CG1')])
def test_new_methyl_branch_corrected_without_moving_any_observed_atom(templates, residue, other):
    top, xyz = invert(templates[residue], 'CG2', ('CB', 'CA', other))
    assert stereochemistry_report(top, xyz * unit.nanometer, templates)['violations']
    observed = {atom_key(a) for a in top.atoms() if a.name != 'CG2'}
    corrected, report = correct_added_branches(top, xyz * unit.nanometer, templates, observed)
    actual = corrected.value_in_unit(unit.nanometer)
    assert not report['stereochemistry_after']['violations']
    assert len(report['corrections']) == 1
    fixed = [a.index for a in top.atoms() if atom_key(a) in observed]
    np.testing.assert_array_equal(actual[fixed], xyz[fixed])
    for a, b in top.bonds():
        assert np.linalg.norm(actual[a.index] - actual[b.index]) == pytest.approx(np.linalg.norm(xyz[a.index] - xyz[b.index]), abs=1e-12)


def test_observed_wrong_stereochemistry_is_never_silently_changed(templates):
    top, xyz = invert(templates['THR'], 'CG2', ('CB', 'CA', 'OG1'))
    observed = {atom_key(a) for a in top.atoms()}
    corrected, report = correct_added_branches(top, xyz * unit.nanometer, templates, observed)
    np.testing.assert_array_equal(corrected.value_in_unit(unit.nanometer), xyz)
    assert not report['corrections'] and report['stereochemistry_after']['violations']


def test_crosslinked_added_branch_is_not_reflected(templates):
    top, xyz = invert(templates['THR'], 'CG2', ('CB', 'CA', 'OG1'))
    named = {a.name: a for a in top.atoms()}
    top.addBond(named['CG2'], named['N'])
    observed = {atom_key(a) for a in top.atoms() if a.name != 'CG2'}
    corrected, report = correct_added_branches(top, xyz * unit.nanometer, templates, observed)
    np.testing.assert_array_equal(corrected.value_in_unit(unit.nanometer), xyz)
    assert not report['corrections'] and report['stereochemistry_after']['violations']


def test_correct_stereo_can_still_collide_with_observed_environment_and_is_rejected(templates):
    top, xyz = invert(templates['THR'], 'CG2', ('CB', 'CA', 'OG1'))
    named = {a.name: a.index for a in top.atoms()}
    reference = np.asarray(templates['THR'].positions.value_in_unit(unit.nanometer))
    external = top.addAtom('C1', app.element.carbon, top.addResidue('LIG', top.addChain('Z'), '99'))
    xyz = np.vstack([xyz, reference[named['CG2']]])
    observed = {atom_key(a) for a in top.atoms() if a.index != named['CG2']}
    corrected, report = correct_added_branches(top, xyz * unit.nanometer, templates, observed)
    assert not report['stereochemistry_after']['violations']
    assert np.linalg.norm(xyz[named['CG2']] - xyz[external.index]) > .13
    # Atom indices are deliberately stale: the final check must resolve names.
    report['corrections'][0]['moved_atom_indices'] = [999999]
    saved = corrected.value_in_unit(unit.nanometer).copy()
    contacts = corrected_branch_contacts(top, corrected, report['corrections'])
    assert not contacts['accepted']
    assert len(contacts['gross_collisions']) == 1
    assert contacts['gross_collisions'][0]['distance_angstrom'] < 1e-10
    np.testing.assert_array_equal(corrected.value_in_unit(unit.nanometer), saved)
    np.testing.assert_array_equal(saved[external.index], xyz[external.index])


def test_branch_contact_gate_does_not_reject_unrelated_observed_contacts(templates):
    top, xyz = invert(templates['THR'], 'CG2', ('CB', 'CA', 'OG1'))
    named = {a.name: a.index for a in top.atoms()}
    external = top.addResidue('LIG', top.addChain('Z'), '99')
    top.addAtom('C1', app.element.carbon, external)
    top.addAtom('C2', app.element.carbon, external)
    xyz = np.vstack([xyz, [5,5,5], [5,5,5]])
    observed = {atom_key(a) for a in top.atoms() if a.index != named['CG2']}
    corrected, report = correct_added_branches(top, xyz * unit.nanometer, templates, observed)
    contacts = corrected_branch_contacts(top, corrected, report['corrections'])
    assert contacts['accepted'] and not contacts['gross_collisions']
    assert contacts['checked_corrected_heavy_atoms'] == 1
    assert np.array_equal(corrected.value_in_unit(unit.nanometer)[-2:], xyz[-2:])


def test_branch_contact_gate_excludes_bonded_and_one_three_pairs(templates):
    top = templates['THR'].topology
    xyz = np.asarray(templates['THR'].positions.value_in_unit(unit.nanometer)).copy()
    named = {a.name: a.index for a in top.atoms()}
    # OG1 and CG2 share CB (a 1–3 pair); do not reinterpret bonded geometry as
    # an environment collision. Bond/angle validation has its own contract.
    xyz[named['OG1']] = xyz[named['CG2']]
    correction = {'identity': list(atom_key(next(a for a in top.atoms() if a.name=='CG2'))[:4]),
                  'moved_atom_names': ['CG2'], 'moved_atom_indices': []}
    assert corrected_branch_contacts(top, xyz * unit.nanometer, [correction])['accepted']
