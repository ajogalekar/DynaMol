"""Preparation invariants and validation; real job evidence lives in docs/audit."""
import shutil
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest
from openmm import app, unit
from pdbfixer import PDBFixer

from backend import config, jobs, preparation, storage
from backend.models import SimulationConfig
from backend.preparation_worker import (
    proper_overlay_points, protonation_inventory, require_valid_stereochemistry,
    rotate_about_axis, stereochemistry_report, torsion_degrees,
)

FIXTURES = config.ROOT / "docs" / "audit" / "preparation-fixtures"


@pytest.fixture(autouse=True)
def preparation_data(tmp_path, monkeypatch):
    # Eligibility controls isolate installation; native runtime has separate checks.
    monkeypatch.setattr("backend.loop_modeling.loop_runtime_status", lambda: {"available": True})
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(config, "DATASETS_DIR", tmp_path / "datasets")
    monkeypatch.setattr(config, "JOBS_DIR", tmp_path / "jobs")
    config.DATASETS_DIR.mkdir()
    config.JOBS_DIR.mkdir()


def imported_fixture(filename):
    source = FIXTURES / filename
    trajectory = md.load(str(source))
    result = storage.save_dataset(trajectory, filename, "Analytical preparation fixture", "Public 1UBQ-derived test fixture")
    shutil.copy2(source, storage.dataset_dir(result["id"]) / "source.pdb")
    return result


def test_inspector_reports_missing_heavy_atoms_and_sequence():
    dataset = imported_fixture("six_residues_missing_atom.pdb")
    result = preparation.inspect_preparation(dataset["id"])
    assert result["can_prepare"]
    assert result["has_sequence"]
    assert any(row["residue"] == "LYS" and "NZ" in row["atoms"] for row in result["missing_atoms"])
    assert result["sequence_source_sha256"]


def test_known_internal_gap_requires_explicit_building(monkeypatch):
    dataset = imported_fixture("six_residues_known_gap.pdb")
    inspection = preparation.inspect_preparation(dataset["id"])
    assert inspection["missing_residues"][0]["residues"] == ["ILE"]
    assert inspection["missing_residues"][0]["buildable"]
    with pytest.raises(ValueError, match="Unresolved internal sequence gaps"):
        preparation.submit_preparation({"dataset_id": dataset["id"]})
    calls = []
    monkeypatch.setattr(preparation, "_submit", lambda settings, operation: calls.append((settings, operation)) or {"status": "queued"})
    result = preparation.submit_preparation({"dataset_id": dataset["id"], "build_missing_residues": True})
    assert result["status"] == "queued"
    assert calls[0][0]["build_missing_residues"] is True
    assert calls[0][1] == "prepare"


def test_numbering_does_not_supply_missing_residue_identities():
    dataset = imported_fixture("six_residues_numbering_gap_only.pdb")
    inspection = preparation.inspect_preparation(dataset["id"])
    assert not inspection["has_sequence"]
    assert inspection["missing_residues"] == []
    assert any(gap["structural_break"] for gap in inspection["gaps"])
    with pytest.raises(ValueError, match="will not be guessed"):
        preparation.submit_preparation({"dataset_id": dataset["id"], "build_missing_residues": True})


def test_raw_simulation_cannot_bypass_structural_gap_guard():
    dataset = imported_fixture("six_residues_numbering_gap_only.pdb")
    with pytest.raises(ValueError, match="artificial stretched peptide bond"):
        jobs.submit_job(SimulationConfig(dataset_id=dataset["id"]))
    assert jobs.list_jobs() == []


def test_prepared_gromacs_requires_validated_state_conversion():
    dataset = imported_fixture("six_residues_intact.pdb")
    dataset["preparation"] = {"ph": 2, "simulation_ready": True}
    storage.atomic_json(storage.dataset_dir(dataset["id"]) / "metadata.json", dataset)
    with pytest.raises(ValueError, match="cannot yet preserve"):
        jobs.submit_job(SimulationConfig(dataset_id=dataset["id"], engine="gromacs", solvent="explicit"))


@pytest.mark.parametrize("target", [-60, 60, 180])
def test_discrete_chi_rotation_reaches_target_without_stretching_bonds(target):
    original = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]], dtype=float)
    changed = original.copy()
    changed[3:] = rotate_about_axis(original[3:], original[1], original[2] - original[1], target - torsion_degrees(original))
    wrapped_difference = (torsion_degrees(changed) - target + 180) % 360 - 180
    assert abs(wrapped_difference) < 1e-8
    assert np.allclose(changed[:3], original[:3])
    assert np.linalg.norm(changed[3] - changed[2]) == pytest.approx(1)


def test_template_alignment_never_uses_an_improper_reflection():
    source = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    reflected = source * [-1, 1, 1]
    translate, rotation, center = proper_overlay_points(reflected, source)
    aligned = (rotation @ (source + translate).T).T + center
    assert np.linalg.det(rotation) == pytest.approx(1)
    assert np.linalg.det(np.stack([aligned[i] - aligned[0] for i in (1, 2, 3)])) > 0
    # Exact overlap of opposite-handed tetrahedra would require an invalid reflection.
    assert not np.allclose(aligned, reflected)


def test_inverted_standard_residue_is_rejected_before_marking_ready():
    fixer = PDBFixer(filename=str(FIXTURES / "six_residues_intact.pdb"))
    baseline = stereochemistry_report(fixer.topology, fixer.positions, fixer.templates)
    assert baseline["violations"] == []
    xyz = np.asarray(fixer.positions.value_in_unit(unit.nanometer)).copy()
    residue = next(fixer.topology.residues())
    indices = {atom.name: atom.index for atom in residue.atoms()}
    center = xyz[indices["CA"]]
    normal = np.cross(xyz[indices["N"]] - center, xyz[indices["C"]] - center)
    normal /= np.linalg.norm(normal)
    xyz[indices["CB"]] -= 2 * np.dot(xyz[indices["CB"]] - center, normal) * normal
    inverted = stereochemistry_report(fixer.topology, xyz * unit.nanometer, fixer.templates)
    assert inverted["violations"][0]["center"] == "CA"
    with pytest.raises(ValueError, match="inverted or near-planar"):
        require_valid_stereochemistry(inverted, "Test structure")


def test_actual_hydrogen_inventory_captures_states_when_openmm_variant_is_none():
    top = app.Topology()
    chain = top.addChain("A")
    definitions = [("ASP", "OD2", ["HD2"]), ("GLU", "OE2", ["HE2"]), ("LYS", "NZ", ["HZ2", "HZ3"])]
    for index, (name, parent_name, hydrogens) in enumerate(definitions):
        residue = top.addResidue(name, chain, str(index + 1))
        parent = top.addAtom(parent_name, app.element.nitrogen if parent_name == "NZ" else app.element.oxygen, residue)
        for hydrogen_name in hydrogens:
            hydrogen = top.addAtom(hydrogen_name, app.element.hydrogen, residue)
            top.addBond(parent, hydrogen)
    states = protonation_inventory(top, [None, None, None])
    assert [record["state"] for record in states] == ["ASH", "GLH", "LYN"]
    assert all(record["openmm_returned_variant"] is None for record in states)
    assert states[2]["bonded_hydrogens"]["NZ"] == ["HZ2", "HZ3"]


def test_prepared_explicit_simulation_requires_visible_saved_preview():
    dataset = imported_fixture("six_residues_intact.pdb")
    dataset["preparation"] = {"ph": 7, "simulation_ready": True}
    storage.atomic_json(storage.dataset_dir(dataset["id"]) / "metadata.json", dataset)
    with pytest.raises(ValueError, match="explicit-water preview first"):
        jobs.submit_job(SimulationConfig(dataset_id=dataset["id"], engine="openmm", solvent="explicit"))
    assert jobs.list_jobs() == []


@pytest.mark.parametrize("requires_minimization,minimize,ready", [(True, False, False), (True, True, True), (False, False, True)])
def test_repaired_loop_minimization_gate_matches_readiness_and_submission(requires_minimization, minimize, ready):
    from backend import readiness
    dataset = imported_fixture("six_residues_intact.pdb")
    folder = storage.dataset_dir(dataset["id"])
    shutil.copy2(folder / "source.pdb", folder / "prepared.pdb")
    dataset["preparation"] = {"ph": 7, "simulation_ready": True, "requires_minimization": requires_minimization}
    storage.atomic_json(folder / "metadata.json", dataset)
    settings = SimulationConfig(dataset_id=dataset["id"], minimize=minimize)
    displayed = readiness.readiness(dataset["id"], readiness.ReadinessRequest(mode="simulation", settings=settings.model_dump()))
    assert displayed["ready"] is ready
    if ready:
        assert jobs.validate_simulation(settings)["metadata"]["id"] == dataset["id"]
    else:
        assert any("Enable energy minimization" in item for item in displayed["blockers"])
        with pytest.raises(ValueError, match="Enable energy minimization"):
            jobs.submit_job(settings)
    assert jobs.list_jobs() == []


def test_recommends_explicit_solvent_by_size_for_protein_only(monkeypatch):
    from backend.preparation_worker import recommends_explicit_solvent
    monkeypatch.setattr(config, "RECOMMEND_EXPLICIT_ATOMS", 4000)
    # Large protein-only system: prefer explicit (implicit GBn2 is O(N^2) on CPU).
    assert recommends_explicit_solvent(4000, requires_explicit=False) is True
    assert recommends_explicit_solvent(50000, requires_explicit=False) is True
    # Small protein-only system: keep the implicit default.
    assert recommends_explicit_solvent(3999, requires_explicit=False) is False
    # Never override a hard explicit requirement (ligand / ion / modified residue).
    assert recommends_explicit_solvent(50000, requires_explicit=True) is False
    # Threshold is the env-configurable constant.
    monkeypatch.setattr(config, "RECOMMEND_EXPLICIT_ATOMS", 100000)
    assert recommends_explicit_solvent(5000, requires_explicit=False) is False
