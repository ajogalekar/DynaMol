"""Observed metal environments must survive protein-only repair operations."""
import numpy as np
from openmm import app, unit

from backend.complex_topology import metal_environment, residue_key, subset
from backend.preparation_worker import adjust_sidechains


def test_near_water_and_donor_are_preserved_without_making_metal_bonds():
    topology = app.Topology()
    chain = topology.addChain("A")
    serine = topology.addResidue("SER", chain, "17")
    carbon = topology.addAtom("CB", app.element.carbon, serine)
    oxygen = topology.addAtom("OG", app.element.oxygen, serine)
    topology.addBond(carbon, oxygen)
    magnesium = topology.addResidue("MG", topology.addChain("M"), "202")
    metal = topology.addAtom("MG", app.element.magnesium, magnesium)
    close_water = topology.addResidue("HOH", topology.addChain("W"), "301")
    near = topology.addAtom("O", app.element.oxygen, close_water)
    far_water = topology.addResidue("HOH", close_water.chain, "302")
    topology.addAtom("O", app.element.oxygen, far_water)
    topology.addBond(metal, oxygen)
    topology.addBond(metal, near)
    positions = np.array([[.36, .03, 0], [.21, 0, 0], [0, 0, 0], [-.21, 0, 0], [2, 0, 0]]) * unit.nanometer
    environment = metal_environment(topology, positions)
    assert environment["water_keys"] == {residue_key(close_water)}
    assert environment["protected_keys"] == {residue_key(serine)}
    selected = {residue_key(serine), residue_key(magnesium), *environment["water_keys"]}
    result = subset(topology, positions, selected)
    assert result.topology.getNumAtoms() == 4
    assert [(a.name, b.name) for a, b in result.topology.bonds()] == [("CB", "OG")]
    np.testing.assert_array_equal(result.positions.value_in_unit(unit.nanometer), positions[:4].value_in_unit(unit.nanometer))
    moved, report = adjust_sidechains(result.topology, result.positions, protected_residues=environment["protected_keys"])
    np.testing.assert_array_equal(moved.value_in_unit(unit.nanometer), result.positions.value_in_unit(unit.nanometer))
    assert report["adjusted_chi_count"] == 0
    assert any("observed metal coordination" in row for row in report["skipped"])
