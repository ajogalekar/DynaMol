"""Project portability and recovery checks in isolated temporary data roots."""
import hashlib
import io
import json
import shutil
import stat
import zipfile
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest
from fastapi import FastAPI, HTTPException
from fastapi.responses import JSONResponse
from fastapi.testclient import TestClient

from backend import config, storage, workspaces


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    def use_root(path):
        path.mkdir(exist_ok=True)
        monkeypatch.setattr(config, "DATA_ROOT", path)
        monkeypatch.setattr(config, "DATASETS_DIR", path / "datasets")
        monkeypatch.setattr(config, "JOBS_DIR", path / "jobs")
        config.DATASETS_DIR.mkdir(exist_ok=True)
        config.JOBS_DIR.mkdir(exist_ok=True)
    use_root(tmp_path / "original")
    top = md.Topology()
    chain = top.add_chain("A")
    residue = top.add_residue("ALA", chain, resSeq=1)
    atoms = [top.add_atom(name, element, residue) for name, element in zip(("N", "CA", "C", "O"), (md.element.nitrogen, md.element.carbon, md.element.carbon, md.element.oxygen))]
    for a, b in ((0, 1), (1, 2), (2, 3)):
        top.add_bond(atoms[a], atoms[b])
    xyz = np.array([[[0, 0, 0], [.14, 0, 0], [.2, .14, 0], [.3, .15, 0]] for _ in range(3)], dtype=np.float32)
    xyz[1, 3, 1] += .03
    xyz[2, 3, 1] += .05
    dataset = storage.save_dataset(md.Trajectory(xyz, top, time=[0, 1, 2]), "Test alanine", "Isolated test", "Synthetic coordinate fixture", dataset_id="project-fixture")
    state = {"version": 1, "dataset_id": dataset["id"], "frame": 2, "camera": np.eye(4).reshape(-1).tolist(), "visibility": {"protein": True, "water": True, "ligands": False, "ions": True, "hydrogens": "polar"}, "representation": "ball+stick", "named_selections": [{"id": "site", "name": "Backbone", "atoms": [0, 1, 2]}], "selected_atoms": [0, 1], "measurements": [{"id": "distance", "kind": "distance", "atoms": [0, 1], "label": "N–CA", "unit": "Å", "values": [1.4, 1.4, 1.4], "times_ps": [0, 1, 2], "color": "#74dfc5", "visible": False}], "active_measurement": "distance"}
    return dataset, state, use_root, tmp_path


def rewrite_archive(source, target, transform):
    with zipfile.ZipFile(source) as archive:
        files = {name: archive.read(name) for name in archive.namelist()}
    manifest = json.loads(files.pop("dynamol-project.json"))
    transform(manifest, files)
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("dynamol-project.json", json.dumps(manifest))
        for name, data in files.items():
            archive.writestr(name, data)


def test_autosave_restores_complete_scene_and_detects_changed_atom_order(workspace):
    dataset, state, _, _ = workspace
    workspaces.save_workspace({"state": state})
    reopened = workspaces.get_workspace()["state"]
    assert reopened["frame"] == 2 and reopened["camera"] == state["camera"]
    assert reopened["named_selections"] == state["named_selections"]
    assert reopened["visibility"]["hydrogens"] == "polar"
    assert reopened["measurements"][0]["visible"] is False
    assert reopened["active_measurement"] == "distance"
    dataset["atoms"][0]["name"] = "NX"
    storage.atomic_json(storage.dataset_dir(dataset["id"]) / "metadata.json", dataset)
    result = workspaces.get_workspace()
    assert result["state"] is None
    assert "atom order" in result["warning"]


def test_fresh_workspace_defaults_to_elements_and_preserves_explicit_color(workspace):
    dataset, state, _, _ = workspace
    workspaces.save_workspace({"state": {**state, "color_scheme": "chain"}})
    assert workspaces.get_workspace()["state"]["color_scheme"] == "chain"
    # A fresh scene should not inherit the previous scene's colors or selections.
    workspaces.save_workspace({"state": {"version": 1, "dataset_id": dataset["id"]}})
    fresh = workspaces.get_workspace()["state"]
    assert fresh["color_scheme"] == "element"
    assert fresh["selected_atoms"] == []
    assert fresh["measurements"] == []


@pytest.mark.parametrize("change", [
    {"selected_atoms": [100]}, {"frame": 3}, {"camera": [0] * 16}, {"camera": [float("nan")] * 16},
    {"visibility": {"hydrogens": "polar"}}, {"named_selections": [{"id": "../escape", "name": "X", "atoms": [0]}]},
    {"measurements": [{"kind": "distance", "atoms": [0, 1], "values": [1], "times_ps": [0]}]},
])
def test_rejects_invalid_or_mismatched_saved_state(workspace, change):
    _, state, _, _ = workspace
    with pytest.raises(ValueError):
        workspaces.save_workspace({"state": {**state, **change}})
    assert workspaces.get_workspace()["state"] is None


def test_named_project_is_a_snapshot_not_changed_by_autosave(workspace):
    _, state, _, _ = workspace
    project = workspaces.save_project({"name": "Backbone study", "state": state})
    workspaces.save_workspace({"state": {**state, "frame": 0, "selected_atoms": []}})
    assert workspaces.get_project(project["id"])["state"]["frame"] == 2
    assert workspaces.get_workspace()["state"]["frame"] == 0
    assert workspaces.list_projects()[0]["name"] == "Backbone study"
    updated = workspaces.save_project({"id": project["id"], "name": "Renamed project", "state": {**state, "frame": 1}})
    assert updated["created_at"] == project["created_at"]
    assert workspaces.get_project(project["id"])["state"]["frame"] == 1


def test_autosave_rejects_changed_physical_coordinates_even_with_identical_atoms(workspace):
    dataset, state, _, _ = workspace
    workspaces.save_workspace({"state": state})
    saved = workspaces.get_workspace()["state"]
    physical = storage.dataset_dir(dataset["id"]) / "physical.npz"
    with np.load(physical, allow_pickle=False) as source:
        data = {key: source[key].copy() for key in source.files}
    data["xyz"][0, 0, 0] += .02
    np.savez_compressed(physical, **data)
    with pytest.raises(ValueError, match="different trajectory coordinates"):
        workspaces.save_workspace({"state": saved})
    assert "different trajectory coordinates" in workspaces.get_workspace()["warning"]


def test_structural_analysis_settings_and_layout_are_saved_without_large_result_arrays(workspace):
    _, state, _, _ = workspace
    settings = {"kind": "rmsf", "selection": "named:site", "alignment": "backbone", "reference": 1, "start": 1, "end": 3, "stride": 1, "periodic": "whole", "residue_average": True}
    workspaces.save_workspace({"state": {**state, "analysis_settings": settings, "analysis_expanded": True}})
    restored = workspaces.get_workspace()["state"]
    assert restored["analysis_settings"] == settings
    assert restored["analysis_expanded"] is True
    workspaces.save_workspace({"state": {**restored, "named_selections": []}})
    assert workspaces.get_workspace()["state"]["analysis_settings"]["selection"] == "ca"


def test_corrupt_current_or_project_record_does_not_hide_the_remaining_library(workspace):
    _, state, _, _ = workspace
    project = workspaces.save_project({"name": "Intact project", "state": state})
    (config.DATA_ROOT / "workspaces" / "projects" / "damaged.json").write_text("partial json")
    (config.DATA_ROOT / "workspaces" / "current.json").write_text("partial json")
    current = workspaces.get_workspace()
    assert current["state"] is None and "could not be read" in current["warning"]
    assert {p["id"] for p in workspaces.list_projects()} == {project["id"], "damaged"}
    assert workspaces.get_project(project["id"])["name"] == "Intact project"


def test_real_tpo_parameter_bundle_remains_loadable_after_portable_import(workspace):
    from openmm import XmlSerializer, app
    from backend.modified_residues import register_topology_definitions, inspect_modified, modified_forcefield_provenance
    from backend.prepared_system import snapshot_modified_parameters, load_prepared_forcefield
    _, _, use_root, tmp_path = workspace
    source = Path(__file__).resolve().parents[2] / "docs" / "audit" / "modified-residues" / "native" / "TPO.pdb"
    register_topology_definitions()
    pdb = app.PDBFile(str(source))
    prepared = {"modified_residues": inspect_modified(pdb.topology, 7)["residues"], "modified_residue_parameters": modified_forcefield_provenance(), "requires_explicit_solvent": True}
    assert prepared["modified_residues"][0]["residue"] == "TPO"
    data = storage.save_dataset(md.load(str(source)), "Capped TPO reference", "AmberTools native fixture", "A reference geometry, not a prepared user protein.", dataset_id="tpo-portable")
    data["preparation"] = prepared
    folder = storage.dataset_dir(data["id"])
    storage.atomic_json(folder / "metadata.json", data)
    shutil.copy2(source, folder / "prepared.pdb")
    snapshot_modified_parameters(folder, prepared)
    ff, _ = load_prepared_forcefield(folder, prepared)
    before = XmlSerializer.serialize(ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff))
    project = workspaces.save_project({"name": "TPO portable parameters", "state": {"version": 1, "dataset_id": data["id"]}})
    backup = tmp_path / "tpo.zip"
    shutil.copy2(workspaces.export_project(project["id"]), backup)
    use_root(tmp_path / "tpo-destination")
    result = workspaces.import_project(backup)
    imported = storage.get_dataset(result["project"]["dataset_id"])
    ff, _ = load_prepared_forcefield(storage.dataset_dir(imported["id"]), imported["preparation"])
    after = XmlSerializer.serialize(ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff))
    assert after == before
    assert imported["preparation"] == prepared


def test_archive_export_import_preserves_sources_parameters_and_run_provenance(workspace):
    dataset, state, use_root, tmp_path = workspace
    folder = storage.dataset_dir(dataset["id"])
    original = folder / "originals"
    original.mkdir()
    (original / "source.pdb").write_bytes((folder / "topology.pdb").read_bytes())
    (folder / "provenance.json").write_text(json.dumps({"source_sha256": "recorded-source", "preparation": "fixture only"}))
    job = config.JOBS_DIR / "complete-fixture"
    job.mkdir()
    (job / "status.json").write_text(json.dumps({"id": "complete-fixture", "name": "Test job", "status": "completed", "config": {"dataset_id": dataset["id"]}}))
    (job / "energy.csv").write_text("step,energy\n1,-2.0\n")
    project = workspaces.save_project({"name": "Portable study", "state": state})
    archive = workspaces.export_project(project["id"])
    backup = tmp_path / "backup.zip"
    shutil.copy2(archive, backup)
    hashes = {str(path.relative_to(folder)): workspaces._digest(path) for path in folder.rglob("*") if path.is_file()}
    use_root(tmp_path / "new-machine")
    result = workspaces.import_project(backup)
    assert result["imported_datasets"] == [dataset["id"]]
    reopened = workspaces.get_project(result["project"]["id"])
    assert reopened["state"]["named_selections"] == state["named_selections"]
    assert reopened["state"]["measurements"][0]["values"] == [1.4, 1.4, 1.4]
    imported = storage.dataset_dir(dataset["id"])
    assert hashes == {str(path.relative_to(imported)): workspaces._digest(path) for path in imported.rglob("*") if path.is_file()}
    assert not list(config.JOBS_DIR.iterdir())
    records = config.DATA_ROOT / "workspaces" / "project-records" / result["project"]["id"]
    assert (records / "complete-fixture" / "energy.csv").is_file()
    repeated = workspaces.import_project(backup)
    assert repeated["reused_datasets"] == [dataset["id"]]
    reexported = workspaces.export_project(result["project"]["id"])
    with zipfile.ZipFile(reexported) as exported:
        assert exported.read("run-records/complete-fixture/energy.csv") == b"step,energy\n1,-2.0\n"


def test_backup_requires_ancestry_and_refuses_conflicting_existing_files(workspace):
    dataset, state, _, tmp_path = workspace
    project = workspaces.save_project({"name": "Study", "state": state})
    archive = workspaces.export_project(project["id"])
    metadata = storage.dataset_dir(dataset["id"]) / "metadata.json"
    data = json.loads(metadata.read_text())
    data["name"] = "Changed locally"
    storage.atomic_json(metadata, data)
    with pytest.raises(HTTPException, match="409"):
        workspaces.import_project(archive)
    assert storage.get_dataset(dataset["id"])["name"] == "Changed locally"
    data["parent_dataset_id"] = "missing-source"
    storage.atomic_json(metadata, data)
    with pytest.raises(FileNotFoundError):
        workspaces.export_project(project["id"])


def test_backup_refuses_active_related_job(workspace):
    dataset, state, _, _ = workspace
    folder = config.JOBS_DIR / "active"
    folder.mkdir()
    storage.atomic_json(folder / "status.json", {"status": "running", "config": {"dataset_id": dataset["id"]}})
    project = workspaces.save_project({"name": "Study", "state": state})
    with pytest.raises(HTTPException, match="409"):
        workspaces.export_project(project["id"])


@pytest.mark.parametrize("attack", ["traversal", "checksum", "coordinate-count", "atom-order", "remote-url", "executable-xml"])
def test_import_rejects_corrupt_and_unsafe_archives_without_partial_writes(workspace, attack):
    dataset, state, use_root, tmp_path = workspace
    project = workspaces.save_project({"name": "Study", "state": state})
    archive = workspaces.export_project(project["id"])
    changed = tmp_path / "attacked.zip"
    prefix = f"datasets/{dataset['id']}/"
    def mutate(manifest, files):
        if attack == "traversal":
            path, data = "datasets/../../escaped.txt", b"bad"
        elif attack == "checksum":
            files[prefix + "coordinates.bin"] += b"bad"
            return
        elif attack == "coordinate-count":
            path, data = prefix + "coordinates.bin", b"too short"
        elif attack == "executable-xml":
            path, data = prefix + "untrusted.xml", b"<ForceField><Script>raise RuntimeError('not executed')</Script></ForceField>"
        else:
            path = prefix + "metadata.json"
            meta = json.loads(files[path])
            if attack == "atom-order":
                meta["atoms"][0]["name"] = "WRONG"
                manifest["datasets"][dataset["id"]]["atom_signature"] = workspaces.atom_signature(meta)
                manifest["project"]["state"]["atom_signature"] = workspaces.atom_signature(meta)
            else:
                meta["coordinates_url"] = "https://example.invalid/collect"
            data = json.dumps(meta).encode()
        files[path] = data
        manifest["files"][path] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    rewrite_archive(archive, changed, mutate)
    use_root(tmp_path / "destination")
    with pytest.raises(ValueError):
        workspaces.import_project(changed)
    assert not list(config.DATASETS_DIR.iterdir())
    assert not workspaces.list_projects()
    assert not (config.DATA_ROOT / "escaped.txt").exists()


def test_import_rejects_symlink(workspace):
    _, state, _, tmp_path = workspace
    project = workspaces.save_project({"name": "Study", "state": state})
    source = workspaces.export_project(project["id"])
    target = tmp_path / "symlink.zip"
    with zipfile.ZipFile(source) as origin, zipfile.ZipFile(target, "w") as output:
        for info in origin.infolist():
            if info.filename.endswith("coordinates.bin"):
                info.external_attr = (stat.S_IFLNK | 0o777) << 16
            output.writestr(info, origin.read(info.filename))
    with pytest.raises(ValueError, match="symlinks"):
        workspaces.import_project(target)


def test_import_rejects_case_collisions_and_archive_expansion(workspace, monkeypatch):
    _, state, _, tmp_path = workspace
    project = workspaces.save_project({"name": "Study", "state": state})
    source = workspaces.export_project(project["id"])
    target = tmp_path / "case-collision.zip"
    with zipfile.ZipFile(source) as origin, zipfile.ZipFile(target, "w") as output:
        for info in origin.infolist():
            output.writestr(info, origin.read(info.filename))
        output.writestr("DYNAMOL-PROJECT.JSON", "duplicate")
    with pytest.raises(ValueError, match="duplicate"):
        workspaces.import_project(target)
    monkeypatch.setattr(workspaces, "MAX_EXPANDED_BYTES", 1)
    with pytest.raises(ValueError, match="expands"):
        workspaces.import_project(source)


def test_import_rejects_nested_npz_bomb(workspace):
    dataset, state, use_root, tmp_path = workspace
    project = workspaces.save_project({"name": "Study", "state": state})
    source = workspaces.export_project(project["id"])
    nested = tmp_path / "physical.npz"
    with zipfile.ZipFile(nested, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("unknown.npy", b"invalid")
    changed = tmp_path / "nested.zip"
    def mutate(manifest, files):
        name = f"datasets/{dataset['id']}/physical.npz"
        data = nested.read_bytes()
        files[name] = data
        manifest["files"][name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    rewrite_archive(source, changed, mutate)
    use_root(tmp_path / "nested-destination")
    with pytest.raises(ValueError, match="coordinate archive exceeds"):
        workspaces.import_project(changed)
    assert not list(config.DATASETS_DIR.iterdir())


def test_import_checks_npy_shapes_before_numpy_can_allocate(workspace, monkeypatch):
    dataset, state, use_root, tmp_path = workspace
    project = workspaces.save_project({"name": "Study", "state": state})
    source = workspaces.export_project(project["id"])
    header = io.BytesIO()
    np.lib.format.write_array_header_1_0(header, {"descr": "<f4", "fortran_order": False, "shape": (10**12, 4, 3)})
    nested = io.BytesIO()
    with zipfile.ZipFile(storage.dataset_dir(dataset["id"]) / "physical.npz") as original, zipfile.ZipFile(nested, "w", zipfile.ZIP_DEFLATED) as archive:
        for entry in original.infolist():
            archive.writestr(entry.filename, header.getvalue() if entry.filename == "xyz.npy" else original.read(entry.filename))
    changed = tmp_path / "huge-header.zip"
    def mutate(manifest, files):
        name = f"datasets/{dataset['id']}/physical.npz"
        data = nested.getvalue()
        files[name] = data
        manifest["files"][name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    rewrite_archive(source, changed, mutate)
    use_root(tmp_path / "header-destination")
    def forbidden_load(*args, **kwargs):
        raise AssertionError("NumPy allocation must not precede header validation")
    monkeypatch.setattr(workspaces.np, "load", forbidden_load)
    with pytest.raises(ValueError, match="array headers"):
        workspaces.import_project(changed)
    assert not list(config.DATASETS_DIR.iterdir())


def test_library_rename_archive_trash_restore_and_protection(workspace):
    dataset, state, _, _ = workspace
    identifier = dataset["id"]
    original = (storage.dataset_dir(identifier) / "physical.npz").read_bytes()
    workspaces.rename_dataset(identifier, {"name": "Renamed alanine"})
    workspaces.archive_dataset(identifier, {"archived": True})
    index = workspaces.library()
    assert index["datasets"][0]["name"] == "Renamed alanine"
    assert index["datasets"][0]["archived"] is True
    assert index["disk"]["datasets_bytes"] > len(original)
    project = workspaces.save_project({"name": "Keep this", "state": state})
    with pytest.raises(HTTPException, match="409"):
        workspaces.trash_dataset(identifier, {"confirm_name": "Renamed alanine"})
    workspaces.delete_project(project["id"], {"confirm": True})
    with pytest.raises(ValueError, match="Type"):
        workspaces.trash_dataset(identifier, {"confirm_name": "wrong"})
    trashed = workspaces.trash_dataset(identifier, {"confirm_name": "Renamed alanine"})
    assert not storage.dataset_dir(identifier).exists()
    assert workspaces.library()["disk"]["trash_bytes"] > len(original)
    workspaces.restore_dataset(trashed["id"])
    assert storage.get_dataset(identifier)["name"] == "Renamed alanine"
    assert (storage.dataset_dir(identifier) / "physical.npz").read_bytes() == original
    assert not workspaces.library()["trash"]


def test_current_workspace_and_job_sources_cannot_be_trashed(workspace):
    dataset, state, _, _ = workspace
    workspaces.save_workspace({"state": state})
    with pytest.raises(HTTPException, match="current workspace"):
        workspaces.trash_dataset(dataset["id"], {"confirm_name": dataset["name"]})
    (config.DATA_ROOT / "workspaces" / "current.json").unlink()
    job = config.JOBS_DIR / "done"
    job.mkdir()
    storage.atomic_json(job / "status.json", {"status": "completed", "config": {"dataset_id": dataset["id"]}})
    with pytest.raises(HTTPException, match="retains this source"):
        workspaces.trash_dataset(dataset["id"], {"confirm_name": dataset["name"]})


def test_router_json_download_and_multipart_import_roundtrip(workspace):
    dataset, state, _, _ = workspace
    app = FastAPI()
    app.include_router(workspaces.router)
    @app.exception_handler(ValueError)
    async def invalid(request, exc):
        return JSONResponse({"detail": str(exc)}, status_code=422)
    with TestClient(app, base_url="http://127.0.0.1") as client:
        saved = client.post("/api/workspace", json={"state": state})
        assert saved.status_code == 200 and len(saved.json()["trajectory_signature"]) == 64
        assert client.get("/api/workspace").json()["state"]["frame"] == 2
        response = client.post("/api/projects", json={"name": "API roundtrip", "state": state})
        assert response.status_code == 200
        project_id = response.json()["id"]
        archive = client.get(f"/api/projects/{project_id}/export")
        assert archive.headers["content-type"] == "application/zip"
        assert not list((config.DATA_ROOT / "workspaces").glob("export-*.zip"))
        imported = client.post("/api/projects/import", files={"file": ("project.zip", archive.content, "application/zip")})
        assert imported.status_code == 200, imported.text
        assert imported.json()["reused_datasets"] == [dataset["id"]]
        assert len(client.get("/api/projects").json()) == 2
        bad = client.post("/api/projects/import", files={"file": ("invalid.zip", b"not a zip", "application/zip")})
        assert bad.status_code == 422
        assert "malformed or unreadable" in bad.json()["detail"]
