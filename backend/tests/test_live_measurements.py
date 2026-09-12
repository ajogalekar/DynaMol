"""Saved-frame tracking, atom identity, incremental XTC and restart controls."""
import json

import mdtraj as md
import numpy as np
import pytest
from fastapi.testclient import TestClient
from mdtraj.formats import XTCTrajectoryFile
from pydantic import ValidationError

from backend import config, live_measurements as live, storage, workspaces
from backend.analysis import _measure
from backend.main import app
from backend.models import LiveMeasurement, MeasurementRequest, RemapMeasurementsRequest, SimulationConfig
from backend.worker import Worker


@pytest.fixture(autouse=True)
def isolated(tmp_path, monkeypatch):
    for name, suffix in [("DATA_ROOT", ""), ("DATASETS_DIR", "datasets"), ("JOBS_DIR", "jobs")]:
        path = tmp_path / suffix
        path.mkdir(exist_ok=True)
        monkeypatch.setattr(config, name, path)


def molecule(*, extra=False, periodic=False, frames=3):
    top = md.Topology()
    residue = top.add_residue("GLY", top.add_chain("A"), resSeq=1)
    for name, element in [("N", md.element.nitrogen), ("H", md.element.hydrogen), ("O", md.element.oxygen), ("C", md.element.carbon)]:
        top.add_atom(name, element, residue)
    top.add_bond(top.atom(0), top.atom(1))
    xyz = [[0, 0, 0], [.1, 0, 0], [.1, .1, 0], [.1, .1, .1]]
    if extra:
        for i in range(12):
            top.add_atom(f"Q{i}", md.element.carbon, residue)
            xyz.append([.01 * i, .02 * i, .03 * i])
    array = np.asarray([xyz] * frames, dtype=np.float32)
    for i in range(frames):
        array[i, 2, 0] += .01 * i
    traj = md.Trajectory(array, top, time=np.arange(frames) * .02)
    if periodic:
        traj.unitcell_vectors = np.tile(np.eye(3, dtype=np.float32), (frames, 1, 1))
    return traj


def definition(kind="distance", atoms=None, id=None):
    return {"id": id or kind, "kind": kind, "atoms": atoms or list(range({"distance": 2, "angle": 3, "dihedral": 4, "hbond": 3}[kind])), "label": f"Test {kind}", "color": "#74dfc5"}


def monitor(tmp_path, traj, definitions=None):
    definitions = definitions or [definition()]
    folder = tmp_path / "monitor"
    folder.mkdir(exist_ok=True)
    traj[0].save_pdb(str(folder / "input.pdb"))
    job = {"id": "tracking-test", "engine": "openmm", "status": "running", "config": {"dataset_id": "source", "measurements": definitions}}
    state = {"measurement_source_atoms": {str(i): i for item in definitions for i in item["atoms"]}}
    return live.LiveMeasurements(folder, job, state, traj.topology), job, state


@pytest.mark.parametrize("change", [
    {"atoms": [0, 0]}, {"atoms": [-1, 1]}, {"atoms": [0.2, 1]}, {"atoms": [True, 1]},
    {"atoms": [0, 1, 2]}, {"kind": "unknown"}, {"id": "../escape"}, {"color": "red"},
])
def test_request_rejects_invalid_selection(change):
    with pytest.raises(ValidationError):
        LiveMeasurement(**{**definition(), **change})


def test_request_bounds_and_unique_ids():
    with pytest.raises(ValidationError, match="unique"):
        SimulationConfig(dataset_id="source", measurements=[definition(), definition()])
    with pytest.raises(ValidationError):
        SimulationConfig(dataset_id="source", measurements=[definition(id=f"m{i}") for i in range(13)])
    assert SimulationConfig(dataset_id="source").measurements == []


def reordered(top, order, *, add_hydrogen=False):
    target = md.Topology()
    residue = target.add_residue("GLY", target.add_chain("A"), resSeq=1)
    if add_hydrogen:
        target.add_atom("H2", md.element.hydrogen, residue)
    mapping = {}
    for index in order:
        atom = top.atom(index)
        mapping[index] = target.add_atom(atom.name, atom.element, residue)
    for a, b in top.bonds:
        if a.index in mapping and b.index in mapping:
            target.add_bond(mapping[a.index], mapping[b.index])
    return target


def test_mapping_survives_added_hydrogens_and_reordered_atoms():
    source = molecule().topology
    target = reordered(source, [3, 1, 0, 2], add_hydrogen=True)
    assert live.map_atoms(source, target, [0, 2, 3]) == {0: 3, 2: 4, 3: 1}


def test_mapping_refuses_lost_or_duplicate_atom_names():
    source = molecule().topology
    target = reordered(source, [0, 1, 3])
    with pytest.raises(ValueError, match="missing or ambiguous"):
        live.map_atoms(source, target, [2])
    target = source.copy()
    target.add_atom("N", md.element.nitrogen, target.residue(0))
    with pytest.raises(ValueError, match="missing or ambiguous"):
        live.map_atoms(source, target, [0])


def test_identical_chain_permutation_is_not_a_valid_ordinal_map():
    def chains(names):
        top = md.Topology()
        for name in names:
            residue = top.add_residue("GLY", top.add_chain(name), resSeq=1)
            top.add_atom("CA", md.element.carbon, residue)
        return top
    with pytest.raises(ValueError, match="chain order"):
        live.map_atoms(chains(["A", "B"]), chains(["B", "A"]), [0])


def test_native_coordinate_validation_rejects_permutation_without_chain_ids(tmp_path):
    traj = molecule()
    traj[0].save_pdb(str(tmp_path / "input.pdb"))
    changed = traj[0]
    for chain in changed.topology.chains:
        chain.chain_id = None
    changed.xyz[0, [0, 2]] = changed.xyz[0, [2, 0]]
    with pytest.raises(ValueError, match="identity or order"):
        live.verify_native_identity(tmp_path, {"measurement_source_atoms": {"0": 0}}, changed)


def test_native_coordinate_validation_respects_gro_rounding(tmp_path):
    traj = molecule()
    traj.xyz[0, 0] += .00049
    traj[0].save_pdb(str(tmp_path / "input.pdb"))
    rounded = traj[0]
    rounded.xyz = np.round(rounded.xyz, 3)
    state = {"measurement_source_atoms": {"0": 0}}
    with pytest.raises(ValueError, match="identity or order"):
        live.verify_native_identity(tmp_path, state, rounded)
    assert live.verify_native_identity(tmp_path, state, rounded, coordinate_tolerance=6e-4) == {"0": 0}


@pytest.mark.parametrize("kind", ["distance", "angle", "dihedral", "hbond"])
def test_actual_saved_frames_match_offline_geometry(kind, tmp_path):
    traj = molecule(periodic=True)
    tracker, job, _ = monitor(tmp_path, traj, [definition(kind)])
    for i in range(traj.n_frames):
        tracker.add_frame(traj.xyz[i], float(traj.time[i]), i * 10, traj.unitcell_vectors[i])
    tracker.finish(traj)
    saved = live.read_snapshot(tracker.folder, job)["measurements"][0]
    expected = _measure(traj, MeasurementRequest(kind=kind, atoms=definition(kind)["atoms"]), strict=False)
    assert saved["values"] == pytest.approx(expected["values"])
    assert saved["times_ps"] == pytest.approx(expected["times_ps"])
    assert saved["output_atoms"] == definition(kind)["atoms"]
    if kind == "hbond":
        assert saved["angle_values"] == pytest.approx(expected["angle_values"])
        assert saved["occupancy"] == sum(expected["geometry_passes"]) / len(expected["geometry_passes"])


def test_periodic_geometry_and_undefined_samples_keep_exact_frame_slots(tmp_path):
    traj = molecule(periodic=True)
    traj.xyz[0, 0] = [.05, 0, 0]
    traj.xyz[0, 1] = [.95, 0, 0]
    traj.xyz[1, 1] = traj.xyz[1, 0]
    tracker, job, _ = monitor(tmp_path, traj)
    for i in range(traj.n_frames):
        tracker.add_frame(traj.xyz[i], float(traj.time[i]), i, traj.unitcell_vectors[i])
    tracker.finish(traj)
    saved = live.read_snapshot(tracker.folder, job)["measurements"][0]
    assert saved["values"][0] == pytest.approx(1, abs=1e-5)
    assert saved["values"][1] is None
    assert "coincident" in saved["frame_errors"][1]
    assert len(saved["values"]) == traj.n_frames
    json.dumps(saved, allow_nan=False)


def test_checkpoint_reconstruction_drops_uncommitted_snapshot_tail(tmp_path):
    traj = molecule()
    tracker, job, state = monitor(tmp_path, traj)
    for i in range(3):
        tracker.add_frame(traj.xyz[i], float(traj.time[i]), i)
    tracker.flush(force=True)
    # The native checkpoint commits only frames 0 and 1. Restart rebuilds those
    # frames, even if an old snapshot contained a later complete frame.
    restored = live.LiveMeasurements(tracker.folder, job, state, traj.topology)
    for i in range(2):
        restored.add_frame(traj.xyz[i], float(traj.time[i]), i)
    restored.add_frame(traj.xyz[1], float(traj.time[1]), 1)
    restored.add_frame(traj.xyz[2], float(traj.time[2]), 2)
    restored.finish(traj)
    assert restored.steps == [0, 1, 2]
    assert restored.snapshot["measurements"][0]["times_ps"] == pytest.approx(traj.time)


def write_xtc(path, traj):
    with XTCTrajectoryFile(str(path), "w") as stream:
        stream.write(traj.xyz, time=traj.time.astype(np.float32), step=np.arange(traj.n_frames, dtype=np.int32) * 10, box=traj.unitcell_vectors)


def test_xtc_retries_incomplete_tail_without_duplicates_or_recomputing_prefix(tmp_path, monkeypatch):
    traj = molecule(extra=True, periodic=True, frames=4)
    tracker, job, _ = monitor(tmp_path, traj, [definition("distance", [0, 2])])
    full = tmp_path / "complete.xtc"
    write_xtc(full, traj)
    with XTCTrajectoryFile(str(full)) as reader:
        offsets = reader.offsets
    raw = full.read_bytes()
    path = tracker.folder / "production.xtc"
    path.write_bytes(raw[:offsets[2] + 24])
    calls = []
    original = tracker.add_frame
    monkeypatch.setattr(tracker, "add_frame", lambda xyz, time, step, box: (calls.append(step), original(xyz, time, step, box))[1])
    tracker.read_xtc()
    assert tracker.xtc_frames == 2
    tracker.read_xtc()
    assert calls == [0, 10]
    path.write_bytes(raw)
    tracker.read_xtc(final=True)
    tracker.finish(md.load(str(full), top=traj.topology))
    assert calls == [0, 10, 20, 30]
    assert tracker.snapshot["measurements"][0]["times_ps"] == pytest.approx(traj.time)


def test_xtc_native_truncate_and_regrow_rebuilds_replaced_tail(tmp_path):
    traj = molecule(extra=True, periodic=True, frames=4)
    tracker, _, _ = monitor(tmp_path, traj, [definition("distance", [0, 2])])
    path = tracker.folder / "production.xtc"
    write_xtc(path, traj)
    tracker.read_xtc(final=True)
    # A native restart can roll back trailing complete frames, and can regrow
    # beyond the previous size before the next observer poll.
    changed = traj[:]
    changed.xyz[2:, 2, 0] += .05
    write_xtc(path, changed)
    tracker.read_xtc(final=True)
    final = md.load(str(path), top=traj.topology)
    tracker.finish(final)
    result = tracker.snapshot["measurements"][0]
    expected = _measure(final, MeasurementRequest(kind="distance", atoms=[0, 2]), strict=True)
    assert tracker.steps == [0, 10, 20, 30]
    assert result["values"] == pytest.approx(expected["values"])


def test_source_validation_rejects_rebuilt_gromacs_hydrogens_and_bad_geometry(tmp_path):
    traj = molecule()
    dataset = storage.save_dataset(traj[0], "Source", "Test", "Analytical fixture")
    path = storage.dataset_dir(dataset["id"]) / "topology.pdb"
    settings = SimulationConfig(dataset_id=dataset["id"], measurements=[definition("distance", [0, 2])])
    assert live.validate_source(settings, path) == {"0": 0, "2": 2}
    settings = SimulationConfig(dataset_id=dataset["id"], engine="gromacs", solvent="explicit", measurements=[definition("hbond")])
    with pytest.raises(ValueError, match="rebuilds input hydrogens"):
        live.validate_source(settings, path)
    bad = SimulationConfig(dataset_id=dataset["id"], measurements=[definition("distance", [0, 999])])
    with pytest.raises(ValueError, match="outside"):
        live.validate_source(bad, path)


def test_poll_endpoint_only_reads_snapshot_and_overlays_current_status(tmp_path, monkeypatch):
    traj = molecule()
    tracker, job, _ = monitor(tmp_path, traj)
    folder = config.JOBS_DIR / job["id"]
    tracker.folder.rename(folder)
    job.update(status="completed", dataset_id="finished")
    storage.atomic_json(folder / "status.json", job)
    monkeypatch.setattr(live, "load_physical", lambda *args: pytest.fail("Polling cannot read trajectories"))
    monkeypatch.setattr(live, "_measure", lambda *args: pytest.fail("Polling cannot calculate geometry"))
    with TestClient(app, base_url="http://127.0.0.1") as client:
        response = client.get(f"/api/jobs/{job['id']}/measurements")
    assert response.status_code == 200
    result = response.json()
    assert result["status"] == "completed" and result["output_dataset_id"] == "finished"
    assert result["source_dataset_id"] == "source"


def test_measurement_failure_does_not_cancel_worker_dynamics(tmp_path):
    tracker, job, _ = monitor(tmp_path, molecule())
    worker = Worker.__new__(Worker)
    worker.job = job
    worker.live_measurements = tracker
    messages = []
    worker.update = lambda **values: messages.append(values)
    worker.report_measurements("add_frame", np.zeros((2, 3)), 0, 0)
    assert job["status"] == "running"
    assert not tracker.active
    assert "invalid coordinates" in live.read_snapshot(tracker.folder, job)["errors"][0]
    assert "Dynamics can continue" in messages[0]["message"]


def test_resume_preserves_failed_native_identity_proof(tmp_path):
    traj = molecule()
    tracker, job, state = monitor(tmp_path, traj)
    worker = Worker.__new__(Worker)
    worker.job, worker.folder, worker.settings, worker.input_state = job, tracker.folder, job["config"], state
    worker.update = lambda **values: None
    changed = traj[0]
    changed.xyz[0, 0] += .1
    worker.verify_measurement_identity(changed)
    assert not json.loads((tracker.folder / "measurement-identity.json").read_text())["verified"]
    resumed = Worker.__new__(Worker)
    resumed.job, resumed.folder, resumed.settings, resumed.input_state = job, tracker.folder, job["config"], state
    resumed.update = lambda **values: None
    resumed.start_measurements(traj.topology)
    assert resumed.live_measurements is None
    assert "identity or order" in live.read_snapshot(tracker.folder, job)["errors"][0]
    assert job["status"] == "running"


def test_prep_remap_keeps_heavy_atoms_but_drops_replaced_hydrogen(tmp_path):
    traj = molecule()
    source = storage.save_dataset(traj[0], "Source", "Test", "Analytical fixture")
    target_top = reordered(traj.topology, [3, 1, 0, 2], add_hydrogen=True)
    xyz = np.asarray([[.01, .02, .03], *[traj.xyz[0, i] for i in [3, 1, 0, 2]]], dtype=np.float32)
    target = storage.save_dataset(md.Trajectory(xyz[None], target_top), "Prepared", "Test", "Analytical fixture")
    target["preparation"] = {"ph": 7}
    storage.atomic_json(storage.dataset_dir(target["id"]) / "metadata.json", target)
    request = RemapMeasurementsRequest(source_dataset_id=source["id"], measurements=[definition("distance", [0, 2]), definition("hbond")])
    result = live.remap_definitions(target["id"], request)
    assert result["measurements"][0]["atoms"] == [3, 4]
    assert len(result["measurements"]) == 1
    assert "replace hydrogen identities" in result["errors"][0]


def test_workspace_retains_tracking_opt_in_and_null_error_slots():
    traj = molecule()
    dataset = storage.save_dataset(traj, "Source", "Test", "Analytical fixture")
    item = {**definition(), "values": [1, None, 1], "times_ps": dataset["times_ps"], "unit": "Å", "frame_errors": [None, "Undefined geometry", None], "trackDuringRun": True}
    state = {"version": 1, "dataset_id": dataset["id"], "measurements": [item]}
    clean = workspaces._validate_state(state)
    assert clean["measurements"][0]["trackDuringRun"] is True
    assert clean["measurements"][0]["frame_errors"] == item["frame_errors"]
    with pytest.raises(ValueError, match="frame errors"):
        workspaces._validate_state({**state, "measurements": [{**item, "frame_errors": [None]}]})
