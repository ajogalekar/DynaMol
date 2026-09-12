import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest
from fastapi.testclient import TestClient

from backend import config, storage
from backend.analysis import measure
from backend.main import app
from backend.models import MeasurementRequest, SimulationConfig


@pytest.fixture(autouse=True)
def isolated_data(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(config, "DATASETS_DIR", tmp_path / "datasets")
    monkeypatch.setattr(config, "JOBS_DIR", tmp_path / "jobs")
    config.DATASETS_DIR.mkdir()
    config.JOBS_DIR.mkdir()


def molecule(xyz, cell=None):
    top = md.Topology()
    chain = top.add_chain("A")
    residue = top.add_residue("TST", chain, resSeq=1)
    names = ["N", "H", "O", "C"]
    elements = [md.element.nitrogen, md.element.hydrogen, md.element.oxygen, md.element.carbon]
    for name, element in zip(names, elements):
        top.add_atom(name, element, residue)
    top.add_bond(top.atom(0), top.atom(1))
    xyz = np.asarray(xyz, dtype=np.float32)
    traj = md.Trajectory(xyz, top, time=np.arange(len(xyz)) * 0.2)
    if cell:
        traj.unitcell_lengths = np.tile(cell, (len(xyz), 1))
        traj.unitcell_angles = np.tile([90, 90, 90], (len(xyz), 1))
    return traj


def save(traj):
    return storage.save_dataset(traj, "Fixture", "Test", "Analytical geometry fixture")


def test_physical_units_and_periodic_distance():
    traj = molecule([[[0.05, 0, 0], [0.95, 0, 0], [0, 0.2, 0], [0, 0, 0.3]]], [1, 1, 1])
    meta = save(traj)
    value = measure(meta["id"], MeasurementRequest(kind="distance", atoms=[0, 1]))
    assert value["values"] == pytest.approx([1.0], abs=1e-5)
    assert value["unit"] == "Å"
    loaded = storage.load_physical(meta["id"])
    assert np.allclose(loaded.xyz, traj.xyz)
    assert meta["has_unitcell"]
    raw = np.fromfile(storage.dataset_dir(meta["id"]) / "coordinates.bin", dtype="<f4")
    assert raw.size == traj.n_frames * traj.n_atoms * 3


def test_angle_and_signed_dihedral():
    meta = save(molecule([[[0, 0, 0], [0.1, 0, 0], [0.1, 0.1, 0], [0.1, 0.1, 0.1]]]))
    angle = measure(meta["id"], MeasurementRequest(kind="angle", atoms=[0, 1, 2]))
    torsion = measure(meta["id"], MeasurementRequest(kind="dihedral", atoms=[0, 1, 2, 3]))
    assert angle["values"] == pytest.approx([90], abs=1e-5)
    assert torsion["values"] == pytest.approx([90], abs=1e-5)


def test_hbond_explicit_connectivity_and_occupancy():
    meta = save(molecule([[[0, 0, 0], [0.1, 0, 0], [0.3, 0, 0], [0, 0.4, 0]], [[0, 0, 0], [0.1, 0, 0], [0.4, 0, 0], [0, 0.4, 0]]]))
    result = measure(meta["id"], MeasurementRequest(kind="hbond", atoms=[0, 1, 2]))
    assert result["occupancy"] == 0.5
    assert result["angle_values"] == pytest.approx([180, 180])
    assert result["values"] == pytest.approx([3, 4])
    with pytest.raises(ValueError, match="bonded to the donor"):
        measure(meta["id"], MeasurementRequest(kind="hbond", atoms=[2, 1, 0]))


@pytest.mark.parametrize("atoms", [[0, 0], [-1, 1], [0, 99], [0, 1, 2]])
def test_invalid_atom_selection(atoms):
    meta = save(molecule([[[0, 0, 0], [0.1, 0, 0], [0.2, 0.1, 0], [0, 0, 0.1]]]))
    with pytest.raises(ValueError):
        measure(meta["id"], MeasurementRequest(kind="distance", atoms=atoms))


def test_undefined_torsion_rejected():
    meta = save(molecule([[[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [0.3, 0, 0]]]))
    with pytest.raises(ValueError, match="collinear"):
        measure(meta["id"], MeasurementRequest(kind="dihedral", atoms=[0, 1, 2, 3]))


def test_dcd_timestamps_are_not_invented(tmp_path):
    traj = molecule([[[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [0.3, 0, 0]]] * 3)
    top, frames = tmp_path / "top.pdb", tmp_path / "traj.dcd"
    traj[0].save_pdb(str(top))
    traj.save_dcd(str(frames))
    loaded, warnings = storage.load_uploaded(top, frames, 1, None)
    assert any("Physical timestamps are unavailable" in warning for warning in warnings)
    meta = storage.save_dataset(loaded, "DCD", "test", "test", warnings=warnings)
    assert meta["time_unit"] == "frame"
    loaded, warnings = storage.load_uploaded(top, frames, 2, 0.2)
    assert loaded.time == pytest.approx([0, 0.4])


def test_same_size_reordered_self_describing_file_rejected(tmp_path):
    traj = molecule([[[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [0.3, 0, 0]]])
    top, frames = tmp_path / "top.pdb", tmp_path / "swapped.pdb"
    traj.save_pdb(str(top))
    traj.topology.atom(0).name = "NX"
    traj.save_pdb(str(frames))
    with pytest.raises(ValueError, match="identities or ordering differ"):
        storage.load_uploaded(top, frames, 1, None)


def test_api_upload_persists_originals_and_binary_matches_metadata(tmp_path):
    traj = molecule([[[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [0.3, 0, 0]]])
    source = tmp_path / "original.pdb"
    traj.save_pdb(str(source))
    with TestClient(app, base_url="http://127.0.0.1") as client:
        response = client.post("/api/datasets/upload", files={"topology": ("original.pdb", source.read_bytes(), "chemical/x-pdb")})
        assert response.status_code == 200, response.text
        meta = response.json()
        folder = storage.dataset_dir(meta["id"])
        assert (folder / "originals" / "topology.pdb").read_bytes() == source.read_bytes()
        assert json.loads((folder / "provenance.json").read_text())["sha256"]["topology"]
        assert len(client.get(meta["coordinates_url"]).content) == meta["n_atoms"] * meta["n_frames"] * 12
        assert client.get(meta["topology_url"]).status_code == 200
        assert client.get("/api/datasets").json()[0]["n_atoms"] == 4


def test_local_origin_and_errors():
    with TestClient(app, base_url="http://127.0.0.1") as client:
        assert client.post("/api/jobs", headers={"Origin": "https://untrusted.example"}, json={}).status_code == 403
        result = client.post("/api/jobs", json={"duration_ps": -1})
        assert result.status_code == 422
        assert isinstance(result.json()["detail"], str)
        assert client.get("/api/datasets/missing").status_code == 404


def test_simulation_caps_and_supported_modes():
    with pytest.raises(ValueError, match="explicit TIP3P"):
        SimulationConfig(dataset_id="test", engine="gromacs", solvent="implicit")
    with pytest.raises(ValueError):
        SimulationConfig(dataset_id="test", timestep_fs=4)
    with pytest.raises(ValueError, match="10,000 frames"):
        SimulationConfig(dataset_id="test", duration_ps=100, report_interval=1)
    with pytest.raises(ValueError):
        SimulationConfig(dataset_id="../etc")


def test_unsupported_chemistry_is_explicit():
    meta = save(molecule([[[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [0.3, 0, 0]]]))
    with TestClient(app, base_url="http://127.0.0.1") as client:
        result = client.post("/api/jobs", json={"dataset_id": meta["id"], "engine": "openmm"})
        assert result.status_code == 422
        assert "no atoms were removed" in result.json()["detail"].lower()


@pytest.mark.parametrize("residue_name", ["A", "DA"])
@pytest.mark.parametrize("include_protein", [False, True])
def test_nucleic_acids_are_viewable_but_rejected_by_protein_simulation(residue_name, include_protein):
    top = md.Topology()
    chain = top.add_chain("A")
    residue = top.add_residue(residue_name, chain, resSeq=1)
    top.add_atom("P", md.element.phosphorus, residue)
    if include_protein:
        protein = top.add_residue("ALA", top.add_chain("B"), resSeq=1)
        top.add_atom("CA", md.element.carbon, protein)
    coordinates = np.zeros((1, top.n_atoms, 3), dtype=np.float32)
    coordinates[0, :, 0] = np.arange(top.n_atoms, dtype=np.float32)
    meta = save(md.Trajectory(coordinates, top))
    assert meta["atoms"][0]["category"] == "nucleic"
    if include_protein:
        assert meta["atoms"][1]["category"] == "protein"
    with TestClient(app, base_url="http://127.0.0.1") as client:
        # The same imported coordinates remain available to the real viewer API.
        assert client.get(meta["topology_url"]).status_code == 200
        assert len(client.get(meta["coordinates_url"]).content) == top.n_atoms * 12
        response = client.post("/api/jobs", json={"dataset_id": meta["id"], "engine": "openmm"})
        assert response.status_code == 422
        assert "standard proteins only" in response.json()["detail"]
        assert "No atoms were removed" in response.json()["detail"]
        assert client.get("/api/jobs").json() == []


def test_health_reuses_subprocess_probe_for_thirty_seconds(monkeypatch):
    from types import SimpleNamespace
    from backend import jobs
    clock = [100.0]
    calls = []
    monkeypatch.setattr(jobs, "_health_cache", None)
    monkeypatch.setattr(jobs.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(jobs, "gromacs_executable", lambda: "/test/gmx")

    def probe(arguments, **kwargs):
        calls.append(arguments)
        return SimpleNamespace(returncode=0, stdout="GROMACS version: 2025.4\n", stderr="")

    monkeypatch.setattr(jobs.subprocess, "run", probe)
    first = jobs.health()
    first["engines"][0]["message"] = "Caller mutation must not alter the cache"
    clock[0] += 29.9
    second = jobs.health()
    assert len(calls) == 1
    assert second["engines"][0]["message"] != first["engines"][0]["message"]
    assert second["engines"][1]["version"] == "2025.4"
    clock[0] += 0.2
    jobs.health()
    assert len(calls) == 2
