from copy import deepcopy
from xml.etree import ElementTree as ET

import numpy as np
import openmm as mm
import pytest

from backend.loop_backbone_reference import dihedral64
from backend.loop_candidate_validation import physical_parameter_comparison


def system_xml():
    system = mm.System()
    for _ in range(3):
        system.addParticle(12)
    force = mm.HarmonicBondForce()
    force.addBond(0, 1, .14, 100)
    force.addBond(1, 2, .15, 200)
    system.addForce(force)
    return mm.XmlSerializer.serialize(system)


def test_only_bond_order_is_ignored():
    first = system_xml()
    tree = ET.fromstring(first)
    bonds = tree.find('Forces/Force/Bonds')
    terms = list(bonds)
    bonds[:] = list(reversed(terms))
    second = ET.tostring(tree, encoding='unicode')
    assert physical_parameter_comparison(first, second)['all_atom_indexed_parameters_equal']
    bonds[0].set('k', '201')
    assert not physical_parameter_comparison(first, ET.tostring(tree, encoding='unicode'))['all_atom_indexed_parameters_equal']


def test_undirected_bond_endpoint_order_is_equivalent_but_multiplicity_is_not():
    first = system_xml()
    tree = ET.fromstring(first)
    bonds = tree.find('Forces/Force/Bonds')
    for bond in bonds:
        left, right = bond.get('p1'), bond.get('p2')
        bond.set('p1', right); bond.set('p2', left)
    assert physical_parameter_comparison(first, ET.tostring(tree, encoding='unicode'))['all_atom_indexed_parameters_equal']
    bonds.append(deepcopy(bonds[0]))
    assert not physical_parameter_comparison(first, ET.tostring(tree, encoding='unicode'))['all_atom_indexed_parameters_equal']


@pytest.mark.parametrize('change', ['duplicate', 'atom', 'mass', 'other_child', 'container_attribute'])
def test_order_normalization_does_not_hide_physical_changes(change):
    first = system_xml()
    tree = ET.fromstring(first)
    bonds = tree.find('Forces/Force/Bonds')
    if change == 'duplicate':
        bonds.append(deepcopy(bonds[0]))
    elif change == 'atom':
        bonds[0].set('p2', '2')
    elif change == 'mass':
        tree.find('Particles/Particle').set('mass', '0')
    elif change == 'other_child':
        ET.SubElement(tree.find('Forces/Force'), 'Extra').set('changed', 'true')
    else:
        bonds.set('changed', 'true')
    assert not physical_parameter_comparison(first, ET.tostring(tree, encoding='unicode'))['all_atom_indexed_parameters_equal']


def test_malformed_force_inventory_is_not_a_parameter_match():
    with pytest.raises(ValueError):
        physical_parameter_comparison('<System/>', '<System/>')
    tree = ET.fromstring(system_xml())
    tree.find('Forces').append(deepcopy(tree.find('Forces/Force')))
    with pytest.raises(ValueError):
        physical_parameter_comparison(system_xml(), ET.tostring(tree, encoding='unicode'))


def test_backbone_angle_orientation_and_degenerate_geometry():
    points = np.array([[0., 1., 0.], [0., 0., 0.], [1., 0., 0.], [1., 0., 1.]])
    assert dihedral64(points) == pytest.approx(90)
    points[3, 2] = -1
    assert dihedral64(points) == pytest.approx(-90)
    with pytest.raises(ValueError):
        dihedral64(np.zeros((4, 3)))
