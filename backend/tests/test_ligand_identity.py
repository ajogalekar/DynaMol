"""Native 6DBK regression: atom identity must survive parameter transfer.

The fixture is already parameterized; these tests never run charge fitting.
Agreement tests parameter transfer, not the physical accuracy of GAFF2.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import shutil
from xml.etree import ElementTree

import numpy as np
import openmm as mm
from openmm import app, unit
import pytest

from backend.forcefield_identity import IdentityForceField
from backend.ligands import _validate_parameter_conversion
from backend.prepared_system import load_prepared_forcefield


FIXTURE = Path(__file__).parent / "fixtures" / "ligand-6dbk"
NAMESPACE = "DML_1f22f93cdebf"
ENERGY_TOLERANCE = 1e-4
FORCE_TOLERANCE = 1e-3


def native_fixture():
    amber = app.AmberPrmtopFile(str(FIXTURE / "ligand.prmtop"))
    xyz = np.asarray(app.AmberInpcrdFile(str(FIXTURE / "ligand.inpcrd")).positions.value_in_unit(unit.nanometer))
    system = amber.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
    return amber.topology, xyz, system


def atom_binding(topology, *, native_order=None):
    xml_atoms = ElementTree.parse(FIXTURE / "ligand.xml").findall("./Residues/Residue/Atom")
    atoms = list(topology.atoms())
    residue = atoms[0].residue
    native_order = list(range(len(atoms))) if native_order is None else native_order
    return {
        "residue_key": [residue.chain.id, residue.id, residue.insertionCode, residue.name],
        "template_name": NAMESPACE,
        "atoms": [{"prepared_name": atom.name, "template_atom_name": xml_atoms[index].get("name"),
                   "native_index": index} for atom, index in zip(atoms, native_order)],
    }


def clone_topology(source, *, order=None, names=None, chain="A", residue_id="1", reverse_bonds=False):
    """Return a new topology and the native indices in its actual atom order."""
    original = list(source.atoms())
    order = list(range(len(original))) if order is None else list(order)
    result = app.Topology()
    residue = result.addResidue("LIG", result.addChain(chain), residue_id)
    copied = {}
    for index in order:
        atom = original[index]
        copied[index] = result.addAtom(atom.name if names is None else names[index], atom.element, residue)
    bonds = list(source.bonds())
    for a, b in reversed(bonds) if reverse_bonds else bonds:
        result.addBond(copied[b.index] if reverse_bonds else copied[a.index],
                       copied[a.index] if reverse_bonds else copied[b.index])
    return result, order


def identity_system(topology, bindings):
    forcefield = IdentityForceField(str(FIXTURE / "ligand.xml"), ligand_atom_maps=bindings)
    return forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None)


def state_values(system, xyz):
    integrator = mm.VerletIntegrator(.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    try:
        context.setPositions(xyz)
        state = context.getState(getEnergy=True, getForces=True)
        return (state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)))
    finally:
        del context, integrator


def assert_equivalent(native, converted, xyz, native_order=None):
    native_order = np.arange(len(xyz)) if native_order is None else np.asarray(native_order)
    for pose in range(3):
        displacement = 0 if pose == 0 else np.random.default_rng(2026 + pose).normal(0, .0003, xyz.shape)
        expected_energy, expected_forces = state_values(native, xyz + displacement)
        actual_energy, actual_forces = state_values(converted, (xyz + displacement)[native_order])
        assert np.isfinite([expected_energy, actual_energy]).all()
        assert np.isfinite(expected_forces).all() and np.isfinite(actual_forces).all()
        assert abs(actual_energy - expected_energy) < ENERGY_TOLERANCE
        np.testing.assert_allclose(actual_forces, expected_forces[native_order], rtol=0, atol=FORCE_TOLERANCE)


def parameter_terms(system, native_order=None):
    """Canonical nonzero force terms keyed by native atom identity.

    XML legitimately omits zero-amplitude torsions. A quartet may be reversed
    in Amber's sign encoding, but other improper permutations are not merged.
    Six decimal places absorb source-file precision in mixed LJ exceptions;
    the separate three-pose energy/force check retains production tolerances.
    """
    indices = list(range(system.getNumParticles())) if native_order is None else native_order
    terms = Counter()

    def numeric(value):
        return round(value.value_in_unit_system(unit.md_unit_system) if unit.is_quantity(value) else float(value), 6)

    def add(kind, atoms, values, reversible=True):
        atoms = tuple(indices[index] for index in atoms)
        if reversible:
            atoms = min(atoms, atoms[::-1])
        terms[(kind, atoms, tuple(numeric(value) for value in values))] += 1

    for index in range(system.getNumParticles()):
        add("mass", [index], [system.getParticleMass(index)])
    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            for index in range(force.getNumBonds()):
                values = force.getBondParameters(index)
                add("bond", values[:2], values[2:])
        elif isinstance(force, mm.HarmonicAngleForce):
            for index in range(force.getNumAngles()):
                values = force.getAngleParameters(index)
                add("angle", values[:3], values[3:])
        elif isinstance(force, mm.PeriodicTorsionForce):
            for index in range(force.getNumTorsions()):
                values = force.getTorsionParameters(index)
                if values[-1].value_in_unit(unit.kilojoule_per_mole) != 0:
                    add("torsion", values[:4], values[4:])
        elif isinstance(force, mm.NonbondedForce):
            for index in range(force.getNumParticles()):
                add("nonbonded", [index], force.getParticleParameters(index))
            for index in range(force.getNumExceptions()):
                values = list(force.getExceptionParameters(index))
                # Sigma has no effect for an excluded (zero-epsilon) pair;
                # native Amber and XML use different placeholder lengths.
                if values[-1].value_in_unit(unit.kilojoule_per_mole) == 0:
                    values[-2] = 0 * unit.nanometer
                add("exception", values[:2], values[2:])
    return terms


def prepared_state(folder, topology, binding, *, legacy=False):
    target = folder / "ligands" / "6dbk" / "ligand.xml"
    target.parent.mkdir(parents=True)
    shutil.copyfile(FIXTURE / "ligand.xml", target)
    ligand = {"key": ":".join(binding["residue_key"]), "template_name": NAMESPACE,
              "ffxml_sha256": hashlib.sha256(target.read_bytes()).hexdigest()}
    if legacy:
        ligand["atom_map"] = [{"ligand_index": entry["native_index"], "prepared_name": entry["prepared_name"],
                               "original_name": entry["prepared_name"], "original_index": entry["native_index"]}
                              for atom, entry in zip(topology.atoms(), binding["atoms"])
                              if atom.element != app.element.hydrogen]
        ligand["hydrogens_added"] = sum(atom.element == app.element.hydrogen for atom in topology.atoms())
    else:
        ligand["parameter_atom_map"] = binding["atoms"]
    return {"simulation_ready": True, "ligand_parameters": {
        "forcefield": "GAFF2", "charge_method": "AM1-BCC", "requires_explicit_solvent": True,
        "ligands": [ligand], "files": [{"path": "ligands/6dbk/ligand.xml", "sha256": ligand["ffxml_sha256"]}],
    }}


def test_retained_fixture_hashes_match_source_provenance():
    provenance = json.loads((FIXTURE / "provenance.json").read_text())
    for name, record in provenance["files"].items():
        data = (FIXTURE / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == record["sha256"]
        assert len(data) == record["bytes"]


def test_graph_only_baseline_reproduces_6dbk_failure_and_four_wrong_impropers():
    topology, xyz, native = native_fixture()
    baseline = app.ForceField(str(FIXTURE / "ligand.xml")).createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None)
    native_energy, native_forces = state_values(native, xyz)
    baseline_energy, baseline_forces = state_values(baseline, xyz)
    assert abs(baseline_energy - native_energy) > ENERGY_TOLERANCE
    assert np.max(np.abs(baseline_forces - native_forces)) > 1
    missing = parameter_terms(native) - parameter_terms(baseline)
    extra = parameter_terms(baseline) - parameter_terms(native)
    assert sum(missing.values()) == sum(extra.values()) == 4
    assert {term[0] for term in missing | extra} == {"torsion"}


def test_verified_identity_matches_all_native_terms_and_three_pose_energies_forces():
    topology, xyz, native = native_fixture()
    converted = identity_system(topology, [atom_binding(topology)])
    assert parameter_terms(converted) == parameter_terms(native)
    assert_equivalent(native, converted, xyz)


def test_production_conversion_check_preserves_strict_tolerances(tmp_path):
    for name in ("ligand.prmtop", "ligand.inpcrd", "ligand.xml"):
        shutil.copyfile(FIXTURE / name, tmp_path / name)
    report = _validate_parameter_conversion(tmp_path, tmp_path / "ligand.xml")
    assert report["passed"] is True
    assert report["energy_tolerance_kj_mol"] == ENERGY_TOLERANCE
    assert report["force_tolerance_kj_mol_nm"] == FORCE_TOLERANCE
    assert len(report["conformations"]) == 3
    assert json.loads((tmp_path / "conversion-validation.json").read_text()) == report


def test_corrupted_numeric_parameters_still_fail_native_conversion_check(tmp_path):
    for name in ("ligand.prmtop", "ligand.inpcrd", "ligand.xml"):
        shutil.copyfile(FIXTURE / name, tmp_path / name)
    path = tmp_path / "ligand.xml"
    document = ElementTree.parse(path)
    atom = document.find("./Residues/Residue/Atom")
    atom.set("charge", str(float(atom.get("charge")) + .01))
    document.write(path)
    with pytest.raises(ValueError, match="energies/forces"):
        _validate_parameter_conversion(tmp_path, path)
    report = json.loads((tmp_path / "conversion-validation.json").read_text())
    assert report["passed"] is False
    assert report["energy_tolerance_kj_mol"] == ENERGY_TOLERANCE
    assert report["force_tolerance_kj_mol_nm"] == FORCE_TOLERANCE


def test_reordered_atoms_and_bonds_keep_native_parameter_identity():
    topology, xyz, native = native_fixture()
    order = np.random.default_rng(914).permutation(topology.getNumAtoms()).tolist()
    reordered, order = clone_topology(topology, order=order, reverse_bonds=True)
    converted = identity_system(reordered, [atom_binding(reordered, native_order=order)])
    assert parameter_terms(converted, order) == parameter_terms(native)
    assert_equivalent(native, converted, xyz, order)


@pytest.mark.parametrize("mismatch", ["name", "duplicate-name", "element", "missing-bond", "external-bond"])
def test_mismatched_runtime_identity_or_connectivity_fails_closed(mismatch):
    topology, _, _ = native_fixture()
    topology, _ = clone_topology(topology)
    binding = atom_binding(topology)
    atoms = list(topology.atoms())
    if mismatch == "name":
        atoms[0].name = "UNMAPPED"
    elif mismatch == "duplicate-name":
        atoms[0].name = atoms[1].name
    elif mismatch == "element":
        atoms[0].element = app.element.carbon
    elif mismatch == "missing-bond":
        topology._bonds.pop()
    else:
        other = topology.addResidue("UNK", topology.addChain("Z"), "9")
        external = topology.addAtom("Z1", app.element.carbon, other)
        topology.addBond(atoms[0], external)
    with pytest.raises(ValueError):
        identity_system(topology, [binding])


@pytest.mark.parametrize("mismatch", ["duplicate-residue", "duplicate-template-atom", "duplicate-native-index",
                                      "swapped-native-indices", "swapped-template-names", "boolean-index",
                                      "missing-field", "missing-atom"])
def test_incomplete_or_ambiguous_provenance_fails_closed(mismatch):
    topology, _, _ = native_fixture()
    bindings = [atom_binding(topology)]
    if mismatch == "duplicate-residue":
        bindings.append(deepcopy(bindings[0]))
    elif mismatch == "duplicate-template-atom":
        bindings[0]["atoms"][1]["template_atom_name"] = bindings[0]["atoms"][0]["template_atom_name"]
    elif mismatch == "duplicate-native-index":
        bindings[0]["atoms"][1]["native_index"] = bindings[0]["atoms"][0]["native_index"]
    elif mismatch in {"swapped-native-indices", "swapped-template-names"}:
        field = "native_index" if mismatch == "swapped-native-indices" else "template_atom_name"
        first, second = bindings[0]["atoms"][:2]
        first[field], second[field] = second[field], first[field]
    elif mismatch == "boolean-index":
        bindings[0]["atoms"][0]["native_index"] = False
    elif mismatch == "missing-field":
        bindings[0]["atoms"][0].pop("prepared_name")
    else:
        bindings[0]["atoms"].pop()
    with pytest.raises(ValueError):
        identity_system(topology, bindings)


@pytest.mark.parametrize("explicit", [False, True], ids=["automatic-selection", "explicit-selection"])
def test_template_selection_cannot_bypass_missing_identity_binding(explicit):
    topology, _, _ = native_fixture()
    residue = next(topology.residues())
    with pytest.raises(ValueError):
        forcefield = IdentityForceField(str(FIXTURE / "ligand.xml"), ligand_atom_maps=[])
        forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff,
                                residueTemplates={residue: NAMESPACE} if explicit else {})


def test_explicit_template_selection_preserves_the_recorded_atom_map():
    topology, xyz, native = native_fixture()
    residue = next(topology.residues())
    forcefield = IdentityForceField(str(FIXTURE / "ligand.xml"), ligand_atom_maps=[atom_binding(topology)])
    converted = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None,
                                       residueTemplates={residue: NAMESPACE})
    assert parameter_terms(converted) == parameter_terms(native)
    assert_equivalent(native, converted, xyz)


def test_modeller_keeps_existing_ligand_hydrogens_and_parameter_identity():
    topology, xyz, native = native_fixture()
    forcefield = IdentityForceField(str(FIXTURE / "ligand.xml"), ligand_atom_maps=[atom_binding(topology)])
    modeller = app.Modeller(topology, xyz * unit.nanometer)
    names = [atom.name for atom in topology.atoms()]
    modeller.addHydrogens(forcefield, pH=7, platform=mm.Platform.getPlatformByName("Reference"))
    assert [atom.name for atom in modeller.topology.atoms()] == names
    np.testing.assert_array_equal(modeller.positions.value_in_unit(unit.nanometer), xyz)
    converted = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    assert parameter_terms(converted) == parameter_terms(native)


@pytest.mark.parametrize("legacy", [False, True], ids=["all-atom-map", "verified-legacy-map"])
def test_factory_and_pdb_reload_preserve_renamed_prepared_atom_identity(tmp_path, legacy):
    topology, xyz, native = native_fixture()
    atoms = list(topology.atoms())
    names = [f"R{index + 1}" if atom.element != app.element.hydrogen else f"H{index - 29}"
             for index, atom in enumerate(atoms)]
    prepared, _ = clone_topology(topology, names=names, chain="C", residue_id="401")
    state = prepared_state(tmp_path, prepared, atom_binding(prepared), legacy=legacy)
    output = io.StringIO()
    app.PDBFile.writeFile(prepared, xyz * unit.nanometer, output, keepIds=True)
    reloaded = app.PDBFile(io.StringIO(output.getvalue()))
    assert [atom.name for atom in reloaded.topology.atoms()] == names
    forcefield, files = load_prepared_forcefield(tmp_path, state)
    converted = forcefield.createSystem(reloaded.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    assert files[-1] == "ligands/6dbk/ligand.xml"
    assert parameter_terms(converted) == parameter_terms(native)
    assert_equivalent(native, converted, np.asarray(reloaded.positions.value_in_unit(unit.nanometer)))


def test_two_ligand_copies_keep_separate_residue_identity_maps():
    topology, xyz, native = native_fixture()
    first, _ = clone_topology(topology, chain="A", residue_id="101")
    order = np.random.default_rng(915).permutation(len(xyz)).tolist()
    second, _ = clone_topology(topology, order=order, chain="B", residue_id="202", reverse_bonds=True)
    modeller = app.Modeller(first, xyz * unit.nanometer)
    modeller.add(second, (xyz[order] + np.array([4., 0., 0.])) * unit.nanometer)
    bindings = [atom_binding(first), atom_binding(second, native_order=order)]
    converted = identity_system(modeller.topology, bindings)
    actual = parameter_terms(converted, list(range(len(xyz))) + order)
    expected = Counter({term: count * 2 for term, count in parameter_terms(native).items()})
    assert actual == expected
    assert np.isfinite(state_values(converted, modeller.positions)[0])


def test_factory_rejects_missing_all_atom_and_unprovable_legacy_binding(tmp_path):
    topology, _, _ = native_fixture()
    state = prepared_state(tmp_path, topology, atom_binding(topology))
    state["ligand_parameters"]["ligands"][0].pop("parameter_atom_map")
    with pytest.raises(ValueError):
        forcefield, _ = load_prepared_forcefield(tmp_path, state)
        forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff)
