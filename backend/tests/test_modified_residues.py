from io import StringIO
from pathlib import Path
import hashlib
import json

import numpy as np
import pytest
from openmm import app, Context, Platform, VerletIntegrator, XmlSerializer, unit
from pdbfixer import PDBFixer

from backend.modified_residues import (
    DATA, SUPPORTED_MODIFIED, inspect_modified, modified_forcefield_files,
    modified_forcefield_provenance, modified_stereochemistry_report,
    register_fixer_templates, register_modified_forcefield, register_topology_definitions,
)

NATIVE = Path(__file__).parents[2] / "docs/audit/modified-residues/native"


def ff():
    return register_modified_forcefield(app.ForceField("amber14/protein.ff14SB.xml", *modified_forcefield_files()))


def model(name):
    register_topology_definitions()
    return app.PDBFile(str(NATIVE / f"{name}.pdb"))


def state(system, positions):
    integrator = VerletIntegrator(.001)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    result = context.getState(getEnergy=True, getForces=True)
    return result.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), result.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)


def reorder(topology, positions):
    replacement, mapped, order = app.Topology(), {}, []
    for chain in topology.chains():
        new_chain = replacement.addChain(chain.id)
        for residue in chain.residues():
            new_residue = replacement.addResidue(residue.name, new_chain, residue.id, residue.insertionCode)
            for atom in reversed(list(residue.atoms())):
                mapped[atom.index] = replacement.addAtom(atom.name, atom.element, new_residue)
                order.append(atom.index)
    for first, second in reversed(list(topology.bonds())):
        replacement.addBond(mapped[first.index], mapped[second.index])
    return replacement, np.asarray(positions)[order], order


@pytest.mark.parametrize("name", sorted(SUPPORTED_MODIFIED))
def test_native_amber_parameters_energy_forces_and_atom_order(name):
    pdb = model(name)
    native = app.AmberPrmtopFile(str(NATIVE / f"{name}.prmtop"))
    xyz = np.asarray(app.AmberInpcrdFile(str(NATIVE / f"{name}.inpcrd")).positions.value_in_unit(unit.nanometer))
    options = {"nonbondedMethod": app.NoCutoff, "constraints": None, "rigidWater": False, "removeCMMotion": False}
    native_system = native.createSystem(**options)
    xml_system = ff().createSystem(pdb.topology, **options)
    shuffled, shuffled_xyz, order = reorder(pdb.topology, xyz)
    shuffled_system = ff().createSystem(shuffled, **options)
    rng = np.random.default_rng(2026)
    for pose in range(3):
        coordinates = xyz if not pose else xyz + rng.normal(0, .002, xyz.shape)
        native_energy, native_forces = state(native_system, coordinates)
        energy, forces = state(xml_system, coordinates)
        reordered_energy, reordered_forces = state(shuffled_system, coordinates[order])
        # TLEAP serializes angular constants using its finite-precision pi;
        # the official XML uses mathematical pi. These bounds cover measured
        # serialization differences, not chemistry uncertainty or MD accuracy.
        assert abs(native_energy - energy) < .002
        assert np.max(np.abs(native_forces - forces)) < .02
        assert abs(reordered_energy - energy) < 1e-7
        assert np.max(np.abs(reordered_forces - forces[order])) < 1e-6
    charge = lambda system: np.array([force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for force in system.getForces() if type(force).__name__ == "NonbondedForce" for i in range(system.getNumParticles())])
    assert np.max(np.abs(charge(native_system) - charge(xml_system))) < 1e-7
    assert abs(charge(xml_system).sum() - (0 if name == "HYP" else -2)) < 1e-5
    # Corrected terms survive serialization independently of Python callbacks.
    restored = XmlSerializer.deserialize(XmlSerializer.serialize(xml_system))
    assert abs(state(restored, xyz)[0] - state(xml_system, xyz)[0]) < 1e-8


@pytest.mark.parametrize("name,missing", [("SEP", "O3P"), ("TPO", "O3P"), ("PTR", "O3P"), ("HYP", "OD1")])
def test_repairs_modified_heavy_atom_and_adds_fixed_hydrogens(name, missing):
    pdb = model(name)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.delete([atom for atom in modeller.topology.atoms() if atom.element == app.element.hydrogen or atom.residue.name == name and atom.name == missing])
    output = StringIO()
    app.PDBFile.writeFile(modeller.topology, modeller.positions, output)
    fixer = PDBFixer(pdbfile=StringIO(output.getvalue()))
    register_fixer_templates(fixer)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    assert [(r.name, [a.name for a in atoms]) for r, atoms in fixer.missingAtoms.items()] == [(name, [missing])]
    fixer.addMissingAtoms(seed=2026)
    modeller = app.Modeller(fixer.topology, fixer.positions)
    forcefield = ff()
    modeller.addHydrogens(forcefield, platform=Platform.getPlatformByName("Reference"))
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff)
    assert np.isfinite(state(system, modeller.positions)[0])
    assert not modified_stereochemistry_report(modeller.topology, modeller.positions)["violations"]
    inspection = inspect_modified(modeller.topology, 7.0)
    assert not inspection["blockers"]
    assert inspection["residues"][0]["residue"] == name
    assert inspection["residues"][0]["formal_charge"] == (0 if name == "HYP" else -2)


def test_inspection_is_explicit_about_fixed_states_terminal_support_and_insertion_codes():
    pdb = model("TPO")
    residue = next(r for r in pdb.topology.residues() if r.name == "TPO")
    residue.insertionCode = "A"
    result = inspect_modified(pdb.topology, 5.0)
    assert result["residues"][0]["formal_charge"] == -2
    assert result["residues"][0]["insertion_code"] == "A"
    assert "does not calculate or change" in result["residues"][0]["protonation_note"]
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.delete([r for r in modeller.topology.residues() if r.name == "ACE"])
    result = inspect_modified(modeller.topology)
    assert any("unsupported terminal" in text for text in result["blockers"])
    assert result["residues"][0]["supported"] is False
    assert modified_forcefield_files({"modified_residues": []}) == []
    assert modified_forcefield_files({"modified_residues": [{"residue": "HYP"}]}) == []


def test_tpo_and_hyp_secondary_stereocenters_are_checked():
    for name, atom_name in [("TPO", "CG2"), ("HYP", "OD1")]:
        pdb = model(name)
        report = modified_stereochemistry_report(pdb.topology, pdb.positions)
        assert report["checked_centers"] == 2 and not report["violations"]
        xyz = np.array(pdb.positions.value_in_unit(unit.nanometer))
        residue = next(r for r in pdb.topology.residues() if r.name == name)
        atoms = {a.name: a.index for a in residue.atoms()}
        center = atoms["CB" if name == "TPO" else "CG"]
        xyz[atoms[atom_name]] = 2 * xyz[center] - xyz[atoms[atom_name]]
        assert modified_stereochemistry_report(pdb.topology, xyz * unit.nanometer)["violations"]


def test_vendor_sources_and_generated_xml_are_reproducible():
    from scripts.generate_modified_residue_xml import generate
    manifest = json.loads((DATA / "manifest.json").read_text())
    for name, metadata in manifest["files"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == metadata["sha256"]
    assert generate() == (DATA / "phosaa14SB.xml").read_text()
    provenance = modified_forcefield_provenance()
    assert "/Users/" not in json.dumps(provenance)
    assert provenance["implementation_sha256"]
    assert "System.xml" in provenance["ptr_improper_correction"]
