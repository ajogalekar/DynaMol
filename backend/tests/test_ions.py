"""Named ion identity and actual Amber TIP3P parameter continuity; no MD."""
from io import StringIO
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pytest
from openmm import app, unit, NonbondedForce, Platform

from backend.ions import ION_STATES, inspect_ions, validate_ion_system
from backend.prepared_system import load_prepared_forcefield


def topology(name, element, charge=None, *, atom_count=1):
    top = app.Topology()
    residue = top.addResidue(name, top.addChain("I"), "3003")
    for index in range(atom_count):
        top.addAtom(name + str(index), app.element.Element.getBySymbol(element), residue, formalCharge=charge)
    return top


def forcefield():
    return app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")


@pytest.mark.parametrize("name", sorted(ION_STATES))
def test_named_states_have_matching_template_and_real_nonbonded_particle(name):
    element, charge = ION_STATES[name]
    top = topology(name, element, charge)
    before = [(atom.name, atom.element.symbol, atom.formalCharge) for atom in top.atoms()]
    record = inspect_ions(top)[0]
    assert record["supported"] and record["formal_charge"] == charge
    params = record["parameters"]
    assert len(params["forcefield_sha256"]) == 64
    assert params["forcefield_file"] == "amber14/tip3p.xml"
    assert params["charge_e"] == charge and params["sigma_nm"] > 0 and params["epsilon_kj_mol"] > 0
    assert params["identity_source"].startswith("https://www.rcsb.org/ligand/")
    ff = forcefield()
    assert ff.getMatchingTemplates(top)[0].name == params["template"]
    result = validate_ion_system(top, ff.createSystem(top), [record])
    assert result[0]["parameters"]["prepared_system_verified"]
    assert [(atom.name, atom.element.symbol, atom.formalCharge) for atom in top.atoms()] == before


def test_zinc_is_the_installed_named_divalent_model_not_a_generic_metal_acceptance():
    row = inspect_ions(topology("ZN", "Zn"))[0]
    assert row["supported"]
    assert row["parameters"]["template"] == "ZN"
    assert row["parameters"]["atom_type"] == "tip3p_standard-Zn2+"
    assert row["parameters"]["charge_e"] == 2
    assert row["parameters"]["sigma_nm"] == pytest.approx(0.22646645415127425)
    assert row["parameters"]["epsilon_kj_mol"] == pytest.approx(0.01381916624)
    assert "not a bonded coordination" in row["model"]
    assert any("Li, P." in source for source in row["parameters"]["references"])


@pytest.mark.parametrize("name,element,charge,atom_count", [
    ("ZN", "Zn", 1, 1), ("ZN", "Zn", 0, 1), ("ZN", "Zn", -2, 1),
    ("ZN", "Zn", 2.5, 1), ("ZN", "Fe", 2, 1), ("ZN", "Zn", 2, 2),
    ("NA", "Mg", 1, 1), ("MG", "Mg", 1, 1), ("CA", "Ca", 0, 1),
    ("METAL", "Zn", 2, 1), ("FE", "Fe", None, 1), ("FE2", "Fe", 2, 1),
    ("AG", "Ag", 1, 1), ("CU1", "Cu", 1, 1),
])
def test_conflicting_or_unregistered_states_are_blocked_without_rewriting(name, element, charge, atom_count):
    top = topology(name, element, charge, atom_count=atom_count)
    row = inspect_ions(top)[0]
    assert not row["supported"] and row["error"]
    assert row["formal_charge"] is None and row["parameters"] is None
    assert [atom.formalCharge for atom in top.atoms()] == [charge] * atom_count
    with pytest.raises(ValueError):
        validate_ion_system(top, forcefield().createSystem(topology("ZN", "Zn")))


@pytest.mark.parametrize("field,value", [("charge", "1"), ("sigma", "0.3"), ("epsilon", "0")])
def test_inspection_rejects_a_selected_forcefield_with_conflicting_particle_parameters(field, value):
    path = Path(app.__file__).parent / "data/amber14/tip3p.xml"
    tree = ElementTree.parse(path)
    if field == "charge":
        target = tree.find("./Residues/Residue[@name='ZN']/Atom")
    else:
        target = tree.find("./NonbondedForce/Atom[@type='tip3p_standard-Zn2+']")
    target.set(field, value)
    changed = app.ForceField(StringIO(ElementTree.tostring(tree.getroot(), encoding="unicode")))
    row = inspect_ions(topology("ZN", "Zn"), forcefield=changed)[0]
    assert not row["supported"] and "verification failed" in row["error"]


def test_no_ion_template_or_ambiguous_residue_identifiers_is_not_ready():
    row = inspect_ions(topology("ZN", "Zn"), forcefield=app.ForceField("amber14/protein.ff14SB.xml"))[0]
    assert not row["supported"]
    top = topology("ZN", "Zn")
    residue = top.addResidue("ZN", next(top.chains()), "3003")
    top.addAtom("ZN", app.element.zinc, residue)
    assert all(not row["supported"] and "not unique" in row["error"] for row in inspect_ions(top))


@pytest.mark.parametrize("parameter", [0, 1, 2])
def test_final_system_validation_rejects_changed_charge_or_lj(parameter):
    top = topology("ZN", "Zn", 2)
    expected = inspect_ions(top)
    system = forcefield().createSystem(top)
    nonbonded = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    values = list(nonbonded.getParticleParameters(0))
    values[parameter] *= 1.1
    nonbonded.setParticleParameters(0, *values)
    with pytest.raises(ValueError, match="actual ion"):
        validate_ion_system(top, system, expected)


def test_final_system_validation_rejects_a_missing_or_substituted_retained_ion():
    expected = inspect_ions(topology("ZN", "Zn", 2))
    top = topology("MG", "Mg", 2)
    with pytest.raises(ValueError, match="inventory changed"):
        validate_ion_system(top, forcefield().createSystem(top), expected)


def test_prepared_protein_retains_zinc_coordinates_identity_and_parameters_through_pdb(tmp_path):
    from pdbfixer import PDBFixer
    from backend import config
    path = config.ROOT / "docs/audit/preparation-fixtures/six_residues_intact.pdb"
    fixer = PDBFixer(filename=str(path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=42)
    model = app.Modeller(fixer.topology, fixer.positions)
    zinc = topology("ZN", "Zn", 2)
    location = np.asarray(model.positions.value_in_unit(unit.nanometer)).max(axis=0) + .8
    model.add(zinc, [location] * unit.nanometer)
    expected = inspect_ions(model.topology)
    ff, _ = load_prepared_forcefield(tmp_path, {"ions": expected}, solvent="explicit")
    model.addHydrogens(ff, pH=7, platform=Platform.getPlatformByName("Reference"))
    system = ff.createSystem(model.topology, nonbondedMethod=app.NoCutoff)
    result = validate_ion_system(model.topology, system, expected)
    index = result[0]["parameters"]["prepared_particle_index"]
    np.testing.assert_allclose(model.positions[index].value_in_unit(unit.nanometer), location, atol=1e-12)
    assert next(residue for residue in model.topology.residues() if residue.name == "ZN").chain.id == "I"
    with (tmp_path / "prepared.pdb").open("w") as handle:
        app.PDBFile.writeFile(model.topology, model.positions, handle, keepIds=True)
    loaded = app.PDBFile(str(tmp_path / "prepared.pdb"))
    reloaded = validate_ion_system(loaded.topology, ff.createSystem(loaded.topology), expected)
    assert reloaded[0]["formal_charge"] == 2
    np.testing.assert_allclose(loaded.positions[reloaded[0]["parameters"]["prepared_particle_index"]].value_in_unit(unit.nanometer), location, atol=0.000051)
