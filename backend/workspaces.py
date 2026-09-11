"""Durable viewer projects and bounded, lossless local backups.

Archives are data, never installations: no executable scripts, external XML
includes, symlinks, or job resumption are imported. Dataset IDs and source bytes
are preserved; a conflicting existing ID is an error rather than an overwrite.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import stat
import tempfile
import threading
import unicodedata
import uuid
import zipfile
from functools import lru_cache
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree

import mdtraj as md
import numpy as np
from fastapi import APIRouter, File, HTTPException, UploadFile
from .file_responses import TemporaryFileResponse
from starlette.concurrency import run_in_threadpool

from . import config, storage

router = APIRouter()
_LOCK = threading.RLock()
MAX_STATE_BYTES = 8 * 1024 * 1024
MAX_ARCHIVE_BYTES = 500 * 1024 * 1024
MAX_EXPANDED_BYTES = 2 * 1024 * 1024 * 1024
MAX_ARCHIVE_FILES = 10_000
TERMINAL = {"completed", "failed", "cancelled", "interrupted"}


def _root() -> Path:
    path = config.DATA_ROOT / "workspaces"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _now():
    return datetime.now(timezone.utc).isoformat()


def _json(path, default=None):
    return json.loads(path.read_text()) if path.is_file() else default


def _digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def atom_signature(dataset: dict) -> str:
    identities = [[a.get(key) for key in ("index", "name", "element", "residue", "resid", "chain")] for a in dataset["atoms"]]
    return hashlib.sha256(json.dumps(identities, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


@lru_cache(maxsize=128)
def _cached_digest(path: str, modified_ns: int, size: int, inode: int) -> str:
    return _digest(Path(path))


def _trajectory_signature(path: Path) -> str:
    info = path.stat()
    return _cached_digest(str(path), info.st_mtime_ns, info.st_size, info.st_ino)


def _name(value) -> str:
    if not isinstance(value, str) or not value.strip() or len(value.strip()) > 120:
        raise ValueError("Choose a name between 1 and 120 characters.")
    if any(ord(c) < 32 for c in value):
        raise ValueError("Names cannot contain control characters.")
    return value.strip()


def _number(value, low, high, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not low <= value <= high:
        raise ValueError(f"Invalid {label} in the saved workspace.")
    return value


def _indices(values, count, *, maximum=100_000):
    if not isinstance(values, list) or len(values) > maximum or any(type(a) is not int or not 0 <= a < count for a in values) or len(set(values)) != len(values):
        raise ValueError("Saved atom selections do not match the dataset atom order.")
    return values


def _validate_state(state: dict, dataset: dict | None = None, physical_path: Path | None = None) -> dict:
    try:
        encoded = json.dumps(state, allow_nan=False)
    except (TypeError, ValueError, RecursionError) as exc:
        raise ValueError("Workspace state must be finite JSON data.") from exc
    if len(encoded.encode()) > MAX_STATE_BYTES or not isinstance(state, dict) or state.get("version") != 1:
        raise ValueError("Unsupported or oversized workspace state (8 MiB maximum).")
    dataset = dataset or storage.get_dataset(storage.safe_id(state.get("dataset_id", "")))
    signature = atom_signature(dataset)
    if state.get("atom_signature") not in (None, signature):
        raise ValueError("Saved selections belong to a different topology or atom order. Open the structure afresh.")
    if state["dataset_id"] != dataset["id"]:
        raise ValueError("Workspace and dataset identifiers differ.")
    count, frames = dataset["n_atoms"], dataset["n_frames"]
    trajectory_signature = _trajectory_signature(physical_path or storage.dataset_dir(dataset["id"]) / "physical.npz")
    if state.get("trajectory_signature") not in (None, trajectory_signature):
        raise ValueError("Saved plots belong to different trajectory coordinates. Open the structure afresh and calculate the plots again.")
    result = {"version": 1, "dataset_id": dataset["id"], "atom_signature": signature, "trajectory_signature": trajectory_signature}
    result["frame"] = _number(state.get("frame", 0), 0, frames - 1, "frame")
    camera = state.get("camera")
    if camera is not None:
        if not isinstance(camera, list) or len(camera) != 16:
            raise ValueError("Saved camera requires a 16-number orientation matrix.")
        for value in camera:
            _number(value, -1e6, 1e6, "camera")
        matrix = np.asarray(camera).reshape(4, 4)
        if abs(float(np.linalg.det(matrix[:3, :3]))) < 1e-12:
            raise ValueError("The saved camera orientation is singular.")
    result["camera"] = camera
    for field, allowed, default in (
        ("representation", {"cartoon", "ball+stick", "licorice", "surface"}, "cartoon"),
        ("color_scheme", {"chain", "residue", "element"}, "residue"),
        ("measure_kind", {"distance", "angle", "dihedral", "hbond"}, "distance"),
    ):
        value = state.get(field, default)
        if value not in allowed:
            raise ValueError(f"Invalid saved {field}.")
        result[field] = value
    visibility = state.get("visibility", {"protein": True, "water": False, "ligands": True, "ions": True, "hydrogens": "none"})
    if not isinstance(visibility, dict) or any(type(visibility.get(k)) is not bool for k in ("protein", "water", "ligands", "ions")) or visibility.get("hydrogens") not in {"all", "polar", "none"}:
        raise ValueError("Invalid saved visibility controls.")
    result["visibility"] = {k: visibility[k] for k in ("protein", "water", "ligands", "ions", "hydrogens")}
    result["selected_atoms"] = _indices(state.get("selected_atoms", []), count)
    result["speed"] = _number(state.get("speed", 1), 0.25, 4, "playback speed")
    for field, default in (("loop", True), ("inspector", True), ("analysis_expanded", False)):
        if type(state.get(field, default)) is not bool:
            raise ValueError(f"Invalid saved {field}.")
        result[field] = state.get(field, default)
    selections = state.get("named_selections", [])
    if not isinstance(selections, list) or len(selections) > 100:
        raise ValueError("A project can contain at most 100 named selections.")
    result["named_selections"] = []
    for selection in selections:
        result["named_selections"].append({"id": storage.safe_id(selection["id"]), "name": _name(selection["name"]), "atoms": _indices(selection["atoms"], count)})
    if len({s["id"] for s in result["named_selections"]}) != len(selections):
        raise ValueError("Named selection identifiers must be unique.")
    measurements = state.get("measurements", [])
    if not isinstance(measurements, list) or len(measurements) > 100:
        raise ValueError("A project can contain at most 100 measurement plots.")
    result["measurements"] = []
    for measurement in measurements:
        kind = measurement.get("kind")
        if kind not in {"distance", "angle", "dihedral", "hbond"}:
            raise ValueError("Invalid saved measurement kind.")
        indices = _indices(measurement.get("atoms"), count, maximum=4)
        if len(indices) != {"distance": 2, "angle": 3, "dihedral": 4, "hbond": 3}[kind]:
            raise ValueError("Saved measurement atom count is incorrect.")
        values = measurement.get("values")
        times = measurement.get("times_ps")
        if not isinstance(values, list) or len(values) != frames or times != dataset["times_ps"]:
            raise ValueError("Saved measurement frames do not match the trajectory.")
        for value in values:
            if value is not None:
                _number(value, -1e12, 1e12, "measurement")
        clean = {k: measurement[k] for k in ("id", "kind", "atoms", "label", "unit", "values", "times_ps", "color")}
        clean["id"] = storage.safe_id(clean["id"])
        clean["label"] = str(clean["label"])[:500]
        if clean["unit"] not in {"Å", "°", "angstrom", "degrees"} or not isinstance(clean["color"], str) or not __import__("re").fullmatch(r"#[0-9a-fA-F]{6}", clean["color"]):
            raise ValueError("Invalid saved measurement appearance.")
        clean["visible"] = measurement.get("visible") is not False
        if "trackDuringRun" in measurement:
            if type(measurement["trackDuringRun"]) is not bool:
                raise ValueError("Saved live measurement selection must be a boolean.")
            clean["trackDuringRun"] = measurement["trackDuringRun"]
        if "frame_errors" in measurement:
            errors = measurement["frame_errors"]
            if not isinstance(errors, list) or len(errors) != frames or any(value is not None and (not isinstance(value, str) or len(value) > 2000) for value in errors):
                raise ValueError("Saved measurement frame errors do not match the trajectory.")
            clean["frame_errors"] = errors
        if "angle_values" in measurement:
            if not isinstance(measurement["angle_values"], list) or len(measurement["angle_values"]) != frames:
                raise ValueError("Saved hydrogen-bond angles have incorrect dimensions.")
            clean["angle_values"] = [_number(v, -360, 360, "angle") if v is not None else None for v in measurement["angle_values"]]
        if "occupancy" in measurement:
            clean["occupancy"] = _number(measurement["occupancy"], 0, 1, "occupancy") if measurement["occupancy"] is not None else None
        clean["warnings"] = [str(v)[:2000] for v in measurement.get("warnings", [])[:100]]
        result["measurements"].append(clean)
    if len({m["id"] for m in result["measurements"]}) != len(measurements):
        raise ValueError("Measurement identifiers must be unique.")
    active = state.get("active_measurement")
    result["active_measurement"] = active if active in {m["id"] for m in result["measurements"]} else None
    analysis = state.get("analysis_settings")
    result["analysis_settings"] = None
    if analysis is not None:
        if not isinstance(analysis, dict) or analysis.get("kind") not in {"rmsd", "rmsf"} or analysis.get("alignment") not in {"ca", "backbone", "selection", "none"} or analysis.get("periodic") not in {"whole", "cartesian"} or type(analysis.get("residue_average")) is not bool:
            raise ValueError("Invalid saved structural-analysis settings.")
        selection = analysis.get("selection")
        allowed_selections = {"ca", "backbone", "protein-heavy", "solute-heavy", "current"} | {"named:" + s["id"] for s in result["named_selections"]}
        # A removed named selection resets the analysis selector, leaving the
        # rest of a workspace savable. The analysis API validates a draft's
        # actual frame range before calculating anything.
        if selection not in allowed_selections:
            selection = "ca"
        clean_analysis = {key: analysis[key] for key in ("kind", "alignment", "periodic", "residue_average")}
        clean_analysis["selection"] = selection
        for key in ("reference", "start", "end", "stride"):
            clean_analysis[key] = _number(analysis.get(key), 0, 1_000_000_000, "analysis " + key)
        result["analysis_settings"] = clean_analysis
    return result


def validate_state(state: dict, dataset: dict | None = None, physical_path: Path | None = None) -> dict:
    try:
        return _validate_state(state, dataset, physical_path)
    except (KeyError, TypeError, AttributeError, OverflowError, RecursionError) as exc:
        raise ValueError("The saved workspace has missing or malformed fields.") from exc


@router.get("/api/workspace")
def get_workspace():
    try:
        saved = _json(_root() / "current.json")
    except (ValueError, OSError) as exc:
        source = _root() / "current.json"
        if source.is_file():
            recovery = _root() / "recovered-records"
            recovery.mkdir(exist_ok=True)
            target = recovery / f"current-{_digest(source)[:16]}.json"
            if not target.exists():
                shutil.copy2(source, target)
        return {"state": None, "warning": f"The previous workspace record could not be read: {exc}. A recovery copy has been retained."}
    if saved:
        try:
            saved["state"] = validate_state(saved["state"])
        except (ValueError, FileNotFoundError, KeyError, TypeError) as exc:
            return {"state": None, "warning": f"The previous workspace could not be restored: {exc}"}
    return saved or {"state": None}


@router.post("/api/workspace")
def save_workspace(payload: dict):
    with _LOCK:
        state = validate_state(payload.get("state"))
        record = {"state": state, "updated_at": _now()}
        storage.atomic_json(_root() / "current.json", record)
        return {"updated_at": record["updated_at"], "atom_signature": state["atom_signature"], "trajectory_signature": state["trajectory_signature"]}


def _project_path(project_id):
    return _root() / "projects" / f"{storage.safe_id(project_id)}.json"


@router.get("/api/projects")
def list_projects():
    result = []
    for path in sorted((_root() / "projects").glob("*.json"), key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            record = _json(path)
            result.append({k: record[k] for k in ("id", "name", "created_at", "updated_at", "dataset_id")})
        except (ValueError, OSError, KeyError, TypeError):
            result.append({"id": path.stem, "name": "Unreadable project · " + path.stem, "created_at": _now(), "updated_at": _now(), "dataset_id": None, "unavailable": True})
    return result


@router.post("/api/projects")
def save_project(payload: dict):
    with _LOCK:
        state = validate_state(payload.get("state"))
        project_id = storage.safe_id(payload.get("id") or uuid.uuid4().hex[:16])
        path = _project_path(project_id)
        previous = _json(path, {})
        record = {"id": project_id, "name": _name(payload.get("name")), "created_at": previous.get("created_at", _now()), "updated_at": _now(), "dataset_id": state["dataset_id"], "state": state}
        path.parent.mkdir(exist_ok=True)
        storage.atomic_json(path, record)
        return record


@router.get("/api/projects/{project_id}")
def get_project(project_id: str):
    record = _json(_project_path(project_id))
    if not record:
        raise FileNotFoundError("This project was not found.")
    record["state"] = validate_state(record["state"])
    return record


@router.post("/api/projects/{project_id}/delete")
def delete_project(project_id: str, payload: dict):
    with _LOCK:
        path = _project_path(project_id)
        if payload.get("confirm") is not True:
            raise ValueError("Confirm removal of this named project. Its molecular datasets remain in the library.")
        if not path.is_file():
            raise FileNotFoundError("This project was not found.")
        trash = _root() / "removed-projects"
        trash.mkdir(exist_ok=True)
        path.rename(trash / f"{path.stem}-{uuid.uuid4().hex[:8]}.json")
        return {"removed": project_id}


def _library_flags():
    return _json(_root() / "library.json", {})


def _files(folder: Path):
    for path in sorted(folder.rglob("*")):
        if path.is_symlink() or not path.resolve().is_relative_to(folder.resolve()):
            raise ValueError("A molecular archive contains a symlink or escaped path.")
        if path.is_file():
            if not stat.S_ISREG(path.stat().st_mode):
                raise ValueError("Only regular files may be archived.")
            yield path
        elif not path.is_dir():
            raise ValueError("Only regular files and folders may be archived.")


def _bytes(folder):
    return sum(p.stat().st_size for p in _files(folder)) if folder.exists() else 0


@router.get("/api/library")
def library():
    flags = _library_flags()
    datasets = []
    for meta in storage.list_datasets():
        datasets.append({**meta, "archived": flags.get(meta["id"], {}).get("archived", False), "bytes": _bytes(storage.dataset_dir(meta["id"]))})
    trash = []
    for path in (_root() / "trash").glob("*/entry.json"):
        entry = _json(path)
        trash.append({**entry, "bytes": _bytes(path.parent / "dataset")})
    total_datasets = sum(d["bytes"] for d in datasets)
    return {"datasets": datasets, "trash": trash, "disk": {"datasets_bytes": total_datasets, "jobs_bytes": _bytes(config.JOBS_DIR), "trash_bytes": sum(t["bytes"] for t in trash), "free_bytes": shutil.disk_usage(config.DATA_ROOT).free}}


@router.post("/api/library/{dataset_id}/rename")
def rename_dataset(dataset_id: str, payload: dict):
    with _LOCK:
        meta = storage.get_dataset(dataset_id)
        meta["name"] = _name(payload.get("name"))
        storage.atomic_json(storage.dataset_dir(dataset_id) / "metadata.json", meta)
        return meta


@router.post("/api/library/{dataset_id}/archive")
def archive_dataset(dataset_id: str, payload: dict):
    with _LOCK:
        storage.get_dataset(dataset_id)
        if type(payload.get("archived")) is not bool:
            raise ValueError("Choose archive or restore.")
        flags = _library_flags()
        flags[dataset_id] = {"archived": payload["archived"]}
        storage.atomic_json(_root() / "library.json", flags)
        return flags[dataset_id]


def _references(value, target):
    if isinstance(value, dict):
        return any((key in {"dataset_id", "parent_dataset_id", "source_dataset_id"} and child == target) or _references(child, target) for key, child in value.items())
    return isinstance(value, list) and any(_references(child, target) for child in value)


def _deletion_blockers(dataset_id):
    blockers = []
    if dataset_id in {"demo", "ubiquitin-start"}:
        blockers.append("This bundled example is retained by DynaMol.")
    current = _json(_root() / "current.json", {})
    if _references(current, dataset_id):
        blockers.append("This structure is open in the current workspace. Open a different structure first.")
    for project in list_projects():
        if project["dataset_id"] == dataset_id:
            blockers.append(f"Named project “{project['name']}” uses this structure. Remove the project or archive the structure instead.")
    for path in config.DATASETS_DIR.glob("*/metadata.json"):
        if path.parent.name == dataset_id:
            continue
        meta = _json(path, {})
        if _references(meta, dataset_id) or _references(_json(path.parent / "provenance.json", {}), dataset_id):
            blockers.append(f"Derived structure “{meta.get('name', path.parent.name)}” needs this source.")
    for path in config.JOBS_DIR.glob("*/status.json"):
        record = _json(path, {})
        if _references(record, dataset_id):
            blockers.append(f"Simulation/preparation “{record.get('name', path.parent.name)}” retains this source or output. Archive it instead.")
    return blockers


@router.post("/api/library/{dataset_id}/trash")
def trash_dataset(dataset_id: str, payload: dict):
    from . import jobs
    with jobs._lock, _LOCK:
        meta = storage.get_dataset(dataset_id)
        if payload.get("confirm_name") != meta["name"]:
            raise ValueError("Type the dataset name to confirm moving it to recoverable trash.")
        blockers = _deletion_blockers(dataset_id)
        if blockers:
            raise HTTPException(409, " ".join(blockers))
        trash_id = uuid.uuid4().hex[:16]
        target = _root() / "trash" / trash_id
        target.mkdir(parents=True)
        entry = {"id": trash_id, "dataset_id": dataset_id, "name": meta["name"], "trashed_at": _now()}
        storage.atomic_json(target / "entry.json", entry)
        storage.dataset_dir(dataset_id).rename(target / "dataset")
        return entry


@router.post("/api/library/trash/{trash_id}/restore")
def restore_dataset(trash_id: str):
    with _LOCK:
        folder = _root() / "trash" / storage.safe_id(trash_id)
        entry = _json(folder / "entry.json")
        if not entry:
            raise FileNotFoundError("This trash entry was not found.")
        destination = storage.dataset_dir(entry["dataset_id"])
        if destination.exists():
            raise HTTPException(409, "A dataset with this identifier already exists; no files were overwritten.")
        (folder / "dataset").rename(destination)
        shutil.rmtree(folder)
        return storage.get_dataset(entry["dataset_id"])


def _ancestry(dataset_id):
    pending, found = [dataset_id], {}
    while pending:
        current = pending.pop()
        if current in found:
            continue
        meta = storage.get_dataset(current)
        found[current] = meta
        if len(found) > 100:
            raise ValueError("A project backup can include at most 100 related datasets.")
        def collect(value):
            if isinstance(value, dict):
                for key, child in value.items():
                    if key in {"parent_dataset_id", "source_dataset_id"} and isinstance(child, str) and child not in found:
                        pending.append(storage.safe_id(child))
                    else:
                        collect(child)
            elif isinstance(value, list):
                for child in value:
                    collect(child)
        collect(meta)
        collect(_json(storage.dataset_dir(current) / "provenance.json", {}))
    return found


def export_project(project_id: str) -> Path:
    with _LOCK:
        project = get_project(project_id)
        datasets = _ancestry(project["dataset_id"])
        entries = []
        for dataset_id in datasets:
            folder = storage.dataset_dir(dataset_id)
            entries.extend((f"datasets/{dataset_id}/{path.relative_to(folder).as_posix()}", path) for path in _files(folder) if not path.name.endswith(".tmp"))
        # Run records are retained as read-only provenance. Importing a backup
        # never inserts a process or resumable job into the live job queue.
        for status in config.JOBS_DIR.glob("*/status.json"):
            job = _json(status, {})
            if not any(_references(job, identifier) for identifier in datasets):
                continue
            if job.get("status") not in TERMINAL:
                raise HTTPException(409, "Wait for related simulations or preparations to stop before exporting a consistent backup.")
            entries.extend((f"run-records/{status.parent.name}/{path.relative_to(status.parent).as_posix()}", path) for path in _files(status.parent) if path.suffix not in {".zip", ".tmp"})
        previous_records = _root() / "project-records" / storage.safe_id(project_id)
        if previous_records.exists():
            existing_names = {name for name, _ in entries}
            entries.extend((f"run-records/{path.relative_to(previous_records).as_posix()}", path) for path in _files(previous_records) if f"run-records/{path.relative_to(previous_records).as_posix()}" not in existing_names)
        if len(entries) > MAX_ARCHIVE_FILES or sum(path.stat().st_size for _, path in entries) > MAX_EXPANDED_BYTES:
            raise ValueError("Project backup exceeds the 2 GiB / 10,000-file archive limit. Use a smaller project.")
        manifest = {"format": "DynaMol project", "version": 1, "created_at": _now(), "project": project, "datasets": {key: {"atom_signature": atom_signature(meta), "n_atoms": meta["n_atoms"], "n_frames": meta["n_frames"]} for key, meta in datasets.items()}, "files": {name: {"bytes": path.stat().st_size, "sha256": _digest(path)} for name, path in entries}}
        target = _root() / f"export-{uuid.uuid4().hex}.zip"
        try:
            with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
                archive.writestr("dynamol-project.json", json.dumps(manifest, allow_nan=False))
                for name, path in entries:
                    archive.write(path, name)
            # A worker writes outside this module's lock. Do not deliver a ZIP
            # if any source changed while its bytes were being collected.
            with zipfile.ZipFile(target) as archive:
                for name, expected in manifest["files"].items():
                    digest = hashlib.sha256()
                    with archive.open(name) as handle:
                        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                            digest.update(chunk)
                    if digest.hexdigest() != expected["sha256"]:
                        raise ValueError("A source file changed during export. Wait for related work to finish and export again.")
            for status in config.JOBS_DIR.glob("*/status.json"):
                record = _json(status, {})
                if record.get("status") not in TERMINAL and any(_references(record, identifier) for identifier in datasets):
                    raise HTTPException(409, "A related run started during export. Wait for it to stop and export again.")
            if target.stat().st_size > MAX_ARCHIVE_BYTES:
                raise ValueError("This backup exceeds the 500 MiB compressed import limit. Use a smaller project.")
            return target
        except Exception:
            target.unlink(missing_ok=True)
            raise


@router.get("/api/projects/{project_id}/export")
def project_download(project_id: str):
    path = export_project(project_id)
    return TemporaryFileResponse(path, media_type="application/zip", filename=f"DynaMol-project-{project_id}.zip")


def _archive_path(value):
    if not isinstance(value, str) or "\\" in value or ":" in value or "\x00" in value:
        raise ValueError("The project archive contains an unsafe path.")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) != value or len(path.parts) < 3 or path.parts[0] not in {"datasets", "run-records"}:
        raise ValueError("The project archive contains an unsafe or unrecognized path.")
    storage.safe_id(path.parts[1])
    return path


def _validate_molecular_folder(folder, expected):
    if (folder / "metadata.json").stat().st_size > 64 * 1024 * 1024:
        raise ValueError("Imported topology metadata exceeds the size limit.")
    meta = _json(folder / "metadata.json")
    if not isinstance(meta, dict) or meta.get("id") != folder.name or atom_signature(meta) != expected.get("atom_signature"):
        raise ValueError("Imported topology metadata has an incorrect atom-order fingerprint.")
    count, frames = meta.get("n_atoms"), meta.get("n_frames")
    if type(count) is not int or type(frames) is not int or not 0 < count <= config.MAX_ATOMS or not 0 < frames <= config.MAX_FRAMES or count * frames * 12 > config.MAX_COORD_BYTES:
        raise ValueError("Imported molecular coordinates exceed local viewer limits.")
    if expected.get("n_atoms") != count or expected.get("n_frames") != frames or len(meta["atoms"]) != count:
        raise ValueError("Imported molecular dimensions disagree with their manifest.")
    coordinates = folder / "coordinates.bin"
    if coordinates.stat().st_size != count * frames * 12 or not np.isfinite(np.memmap(coordinates, dtype="<f4", mode="r")).all():
        raise ValueError("Imported display coordinates have invalid dimensions or values.")
    # Bound nested NPZ decompression before NumPy allocates any coordinate array.
    with zipfile.ZipFile(folder / "physical.npz") as archive:
        allowed = {"xyz.npy", "time.npy", "lengths.npy", "angles.npy"}
        if len(archive.namelist()) != len(set(archive.namelist())) or any(entry.filename not in allowed or entry.compress_type not in {zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED} or entry.file_size > config.MAX_COORD_BYTES + 1024 * 1024 for entry in archive.infolist()) or sum(entry.file_size for entry in archive.infolist()) > config.MAX_COORD_BYTES + 2 * 1024 * 1024:
            raise ValueError("Imported physical coordinate archive exceeds its bounds.")
        if not {"xyz.npy", "time.npy"}.issubset(archive.namelist()):
            raise ValueError("Imported physical coordinates or timestamps are missing.")
        for entry in archive.infolist():
            with archive.open(entry) as raw:
                version = np.lib.format.read_magic(raw)
                if version == (1, 0):
                    shape, _, dtype = np.lib.format.read_array_header_1_0(raw, max_header_size=10_000)
                elif version == (2, 0):
                    shape, _, dtype = np.lib.format.read_array_header_2_0(raw, max_header_size=10_000)
                else:
                    raise ValueError("Unsupported physical coordinate array encoding.")
                expected_shape = (frames, count, 3) if entry.filename == "xyz.npy" else (frames,) if entry.filename == "time.npy" else (frames, 3)
                if shape != expected_shape or dtype.kind not in "fi" or dtype.itemsize not in (4, 8) or dtype.hasobject or math.prod(shape) * dtype.itemsize + raw.tell() != entry.file_size:
                    raise ValueError("Imported physical array headers disagree with bounded molecular dimensions.")
                if entry.filename == "xyz.npy" and (dtype.kind != "f" or math.prod(shape) * dtype.itemsize > config.MAX_COORD_BYTES):
                    raise ValueError("Imported physical coordinates exceed their memory limit.")
    with np.load(folder / "physical.npz", allow_pickle=False) as physical:
        if physical["xyz"].dtype.kind != "f" or physical["time"].dtype.kind not in "fi" or physical["xyz"].shape != (frames, count, 3) or physical["time"].shape != (frames,) or not np.isfinite(physical["xyz"]).all() or not np.isfinite(physical["time"]).all() or physical["time"].tolist() != meta["times_ps"]:
            raise ValueError("Imported physical trajectory does not match its saved frames.")
        if "lengths" in physical or "angles" in physical:
            for key in ("lengths", "angles"):
                if physical[key].shape != (frames, 3) or not np.isfinite(physical[key]).all() or not (physical[key] > 0).all():
                    raise ValueError("Imported periodic cell is invalid.")
    if (folder / "topology.pdb").stat().st_size > 32 * 1024 * 1024:
        raise ValueError("Imported canonical topology exceeds the size limit.")
    topology = md.load_topology(str(folder / "topology.pdb"))
    if topology.n_atoms != count:
        raise ValueError("Imported topology atom count differs from its trajectory.")
    for index, atom in enumerate(topology.atoms):
        stored = meta["atoms"][index]
        if stored["index"] != index or atom.name.upper() != stored["name"].upper() or atom.residue.name.upper() != stored["residue"][:3].upper() or (stored["element"] != "X" and (atom.element is None or atom.element.symbol.upper() != stored["element"].upper())):
            raise ValueError("Imported topology atom identities or order differ from its metadata.")
    for bond in meta.get("bonds", []):
        if len(_indices(bond, count, maximum=2)) != 2:
            raise ValueError("Invalid imported topology bond.")
    chemical_path = folder / "chemistry.json"
    if chemical_path.is_file():
        if chemical_path.stat().st_size > 64 * 1024 * 1024:
            raise ValueError("Imported chemical graph exceeds the size limit.")
        chemical = _json(chemical_path)
        if chemical.get("authoritative_bonds"):
            if len(chemical.get("atoms", [])) != count or any(atom.get("index") != index for index, atom in enumerate(chemical["atoms"])):
                raise ValueError("Imported chemical graph has inconsistent atom ordering.")
            chemical_bonds = []
            for bond in chemical.get("bonds", []):
                indices = _indices(bond.get("atoms"), count, maximum=2)
                if len(indices) != 2:
                    raise ValueError("Imported chemical graph has an invalid bond.")
                if bond.get("order") is not None:
                    _number(bond["order"], 0, 6, "chemical bond order")
                chemical_bonds.append(indices)
            if chemical_bonds != meta.get("bonds"):
                raise ValueError("Imported chemical graph and displayed bonds disagree.")
    # Existing validators check checksums, paths, source version and executable
    # XML restrictions without loading or running an OpenMM force field.
    from .prepared_system import ligand_parameter_files, _validate_modified_snapshot
    preparation = meta.get("preparation")
    ligand_parameter_files(folder, preparation)
    _validate_modified_snapshot(folder, preparation)
    # Canonical URL paths prevent archives from making the UI fetch remote data.
    meta["topology_url"] = f"/api/datasets/{meta['id']}/topology"
    meta["coordinates_url"] = f"/api/datasets/{meta['id']}/coordinates"
    return meta


def import_project(source: Path) -> dict:
    with _LOCK, tempfile.TemporaryDirectory(prefix="project-import-", dir=_root()) as temporary:
        staging = Path(temporary)
        try:
            with zipfile.ZipFile(source) as archive:
                infos = archive.infolist()
                names = [entry.filename for entry in infos]
                portable_names = {unicodedata.normalize("NFC", name).casefold() for name in names}
                if len(infos) > MAX_ARCHIVE_FILES + 1 or len(names) != len(portable_names) or "dynamol-project.json" not in names:
                    raise ValueError("Invalid project archive manifest or duplicate filenames.")
                if sum(entry.file_size for entry in infos) > MAX_EXPANDED_BYTES + MAX_STATE_BYTES or any(entry.file_size > MAX_EXPANDED_BYTES for entry in infos):
                    raise ValueError("Project archive expands beyond the 2 GiB limit.")
                if archive.getinfo("dynamol-project.json").file_size > MAX_STATE_BYTES + 4 * 1024 * 1024:
                    raise ValueError("Project manifest is too large.")
                manifest = json.loads(archive.read("dynamol-project.json"))
                if manifest.get("format") != "DynaMol project" or manifest.get("version") != 1 or not isinstance(manifest.get("files"), dict) or set(names) != {"dynamol-project.json", *manifest["files"]}:
                    raise ValueError("Unsupported project format or incomplete file manifest.")
                for entry in infos:
                    mode = entry.external_attr >> 16
                    if entry.is_dir() or stat.S_ISLNK(mode) or (stat.S_IFMT(mode) not in (0, stat.S_IFREG)) or entry.flag_bits & 1 or entry.compress_type not in {zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED}:
                        raise ValueError("Project backups cannot contain symlinks, special files, directories or encrypted entries.")
                    if entry.filename == "dynamol-project.json":
                        continue
                    relative = _archive_path(entry.filename)
                    expected = manifest["files"][entry.filename]
                    if entry.file_size != expected["bytes"]:
                        raise ValueError("Project file size does not match its manifest.")
                    destination = staging.joinpath(*relative.parts)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    digest, size = hashlib.sha256(), 0
                    with archive.open(entry) as origin, destination.open("wb") as output:
                        for chunk in iter(lambda: origin.read(1024 * 1024), b""):
                            size += len(chunk)
                            if size > entry.file_size:
                                raise ValueError("Project file expands beyond its declared size.")
                            digest.update(chunk)
                            output.write(chunk)
                    if digest.hexdigest() != expected["sha256"]:
                        raise ValueError("Project file checksum mismatch; no molecular files were imported.")
                    if destination.suffix.lower() == ".xml":
                        if size > 64 * 1024 * 1024:
                            raise ValueError("An imported XML file exceeds the size limit.")
                        data = destination.read_bytes()
                        if b"<!DOCTYPE" in data.upper() or b"<!ENTITY" in data.upper():
                            raise ValueError("Imported XML cannot contain document types or entities.")
                        document = ElementTree.fromstring(data)
                        if any(node.tag.rsplit("}", 1)[-1] in {"Script", "InitializationScript", "Include"} for node in document.iter()):
                            raise ValueError("Imported force-field XML cannot contain scripts or external includes.")
                dataset_manifest = manifest.get("datasets", {})
                folders = {path.name: path for path in (staging / "datasets").iterdir()}
                if not isinstance(dataset_manifest, dict) or not dataset_manifest or len(dataset_manifest) > 100 or set(folders) != set(dataset_manifest):
                    raise ValueError("The project dataset manifest is incomplete.")
                verified = {key: _validate_molecular_folder(folder, dataset_manifest[key]) for key, folder in folders.items()}
                project = manifest["project"]
                if project["dataset_id"] not in verified:
                    raise ValueError("The project's active dataset is missing from its backup.")
                state = validate_state(project["state"], verified[project["dataset_id"]], folders[project["dataset_id"]] / "physical.npz")
                for meta in verified.values():
                    if meta["topology_url"] != _json(folders[meta["id"]] / "metadata.json")["topology_url"] or meta["coordinates_url"] != _json(folders[meta["id"]] / "metadata.json")["coordinates_url"]:
                        raise ValueError("Imported dataset URLs must point to their canonical local files.")
                    # All ancestry must travel together, even when this machine
                    # already happens to contain the referenced parent.
                    def check_parents(value):
                        if isinstance(value, dict):
                            for key, child in value.items():
                                if key in {"parent_dataset_id", "source_dataset_id"} and isinstance(child, str) and child not in verified:
                                    raise ValueError("A source dataset required by this project is missing from the backup.")
                                check_parents(child)
                        elif isinstance(value, list):
                            for child in value:
                                check_parents(child)
                    check_parents(meta)
                    check_parents(_json(folders[meta["id"]] / "provenance.json", {}))
                reuse = []
                for key, folder in folders.items():
                    target = storage.dataset_dir(key)
                    if target.exists():
                        incoming = {str(path.relative_to(folder)): _digest(path) for path in _files(folder)}
                        existing = {str(path.relative_to(target)): _digest(path) for path in _files(target)}
                        if incoming != existing:
                            raise HTTPException(409, f"Dataset '{key}' already exists with different files. No files were overwritten. Import this backup into a separate DynaMol workspace.")
                        reuse.append(key)
                # Stage all project records before committing dataset folders.
                project_id = uuid.uuid4().hex[:16]
                record = {"id": project_id, "name": _name(project["name"]), "created_at": project.get("created_at", _now()), "updated_at": _now(), "dataset_id": state["dataset_id"], "state": state, "imported_from": {"project_id": project.get("id"), "archive_sha256": _digest(source), "imported_at": _now()}}
                project_path = _project_path(project_id)
                project_path.parent.mkdir(exist_ok=True)
                moved = []
                imported_records = _root() / "project-records" / project_id
                try:
                    for key, folder in folders.items():
                        if key not in reuse:
                            target = storage.dataset_dir(key)
                            folder.rename(target)
                            moved.append(target)
                    if (staging / "run-records").exists():
                        imported_records.parent.mkdir(exist_ok=True)
                        (staging / "run-records").rename(imported_records)
                    storage.atomic_json(project_path, record)
                except Exception:
                    for target in moved:
                        target.rename(staging / "datasets" / target.name)
                    shutil.rmtree(imported_records, ignore_errors=True)
                    project_path.unlink(missing_ok=True)
                    raise
                return {"project": record, "imported_datasets": [key for key in folders if key not in reuse], "reused_datasets": reuse, "run_records": "Historical run files were preserved for provenance; jobs were not started or added to the live queue."}
        except (zipfile.BadZipFile, KeyError, TypeError, OSError, ElementTree.ParseError, RecursionError) as exc:
            raise ValueError(f"The project backup is malformed or unreadable: {exc}") from exc


@router.post("/api/projects/import")
async def project_upload(file: UploadFile = File(...)):
    with tempfile.NamedTemporaryFile(prefix="project-upload-", suffix=".zip", dir=_root()) as temporary:
        size = 0
        while chunk := await file.read(1024 * 1024):
            size += len(chunk)
            if size > MAX_ARCHIVE_BYTES:
                raise HTTPException(413, "Project backups must be at most 500 MiB compressed.")
            temporary.write(chunk)
        temporary.flush()
        return await run_in_threadpool(import_project, Path(temporary.name))
