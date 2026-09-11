"""Read-only geometry on saved production frames, with explicit atom identities.

Workers publish atomic, bounded JSON snapshots. Browser polling never opens a
trajectory or recomputes an old frame. Restarting a worker reconstructs the
series from the engine's retained frame prefix, not from a stale snapshot.
"""
from __future__ import annotations

import json
import hashlib
import time
from collections import defaultdict
from pathlib import Path

import mdtraj as md
import numpy as np

from .analysis import _measure
from .models import MeasurementRequest, SimulationConfig, RemapMeasurementsRequest
from .storage import atomic_json, load_physical, get_dataset


_RESIDUE_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "WAT": "HOH", "SOL": "HOH",
}


def _residue_name(residue):
    return _RESIDUE_ALIASES.get(residue.name, residue.name)


def map_atoms(source: md.Topology, target: md.Topology, selected: list[int]) -> dict[int, int]:
    """Map named atoms within a verified residue sequence, never by atom offset.

    Engine adapters preserve input residue order and append solvent. GRO drops
    chain IDs, so residue ordinal is the invariant; names and element distinguish
    atoms after hydrogen insertion/reordering. Ambiguities are refused, including
    duplicate names even when one happens to occupy the previous numeric index.
    """
    source_residues, target_residues = list(source.residues), list(target.residues)
    if len(target_residues) < len(source_residues) or any(
        _residue_name(a) != _residue_name(b) for a, b in zip(source_residues, target_residues)
    ):
        raise ValueError("Tracked atom identities cannot be preserved: the engine changed the input residue sequence.")
    if any(a.chain.chain_id and b.chain.chain_id and a.chain.chain_id != b.chain.chain_id
           for a, b in zip(source_residues, target_residues)):
        raise ValueError("Tracked atom identities cannot be preserved: the engine changed the input chain order.")
    def key(atom):
        return atom.residue.index, atom.name, atom.element.symbol if atom.element else None
    source_keys, target_keys = defaultdict(list), defaultdict(list)
    for atom in source.atoms:
        source_keys[key(atom)].append(atom.index)
    for atom in target.atoms:
        target_keys[key(atom)].append(atom.index)
    result = {}
    for index in selected:
        if index < 0 or index >= source.n_atoms:
            raise ValueError("A tracked atom index is outside the source topology.")
        atom = source.atom(index)
        identity = key(atom)
        if len(source_keys[identity]) != 1 or len(target_keys[identity]) != 1:
            raise ValueError(f"Tracked atom {atom.residue.name} {atom.residue.resSeq} {atom.name} is missing or ambiguous after engine preparation. Choose an unambiguous retained atom; no substitute was selected.")
        result[index] = target_keys[identity][0]
    return result


def verify_native_identity(folder: Path, input_state: dict, target: md.Trajectory, *, coordinate_tolerance=2e-4):
    """Verify selected retained atoms before translation or native integration."""
    indices = list(input_state.get("measurement_source_atoms", {}).values())
    if not indices:
        return {}
    source = md.load(str(folder / "input.pdb"))
    mapping = map_atoms(source.topology, target.topology, indices)
    for index, output in mapping.items():
        if source.topology.atom(index).element != md.element.hydrogen and not np.allclose(
            source.xyz[0, index], target.xyz[0, output], atol=coordinate_tolerance, rtol=0
        ):
            raise ValueError("Engine preparation changed the identity or order of a tracked heavy atom. Its live measurement was not accepted.")
    return {str(index): output for index, output in mapping.items()}


def remap_definitions(target_id: str, request: RemapMeasurementsRequest) -> dict:
    source_metadata, target_metadata = get_dataset(request.source_dataset_id), get_dataset(target_id)
    source, target = load_physical(request.source_dataset_id)[0], load_physical(target_id)[0]
    def key(atom):
        residue = atom.residue
        return residue.chain.chain_id, residue.resSeq, _residue_name(residue), atom.name, atom.element.symbol if atom.element else None
    source_keys, target_keys = defaultdict(list), defaultdict(list)
    for atom in source.topology.atoms:
        source_keys[key(atom)].append(atom.index)
    for atom in target.topology.atoms:
        target_keys[key(atom)].append(atom.index)
    result = {"measurements": [], "warnings": [], "errors": []}
    for definition in request.measurements:
        try:
            if any(i >= source.n_atoms for i in definition.atoms):
                raise ValueError("selected atoms are outside the source structure")
            atoms = [source.topology.atom(i) for i in definition.atoms]
            if source_metadata.get("preparation") != target_metadata.get("preparation") and any(atom.element == md.element.hydrogen for atom in atoms):
                raise ValueError("protein preparation can replace hydrogen identities; select the hydrogens again")
            mapped = []
            for atom in atoms:
                identity = key(atom)
                if len(source_keys[identity]) != 1 or len(target_keys[identity]) != 1:
                    raise ValueError(f"{atom.residue.name} {atom.residue.resSeq} {atom.name} is missing or ambiguous in the new structure")
                mapped.append(target_keys[identity][0])
            _measure(target, MeasurementRequest(kind=definition.kind, atoms=mapped), strict=True)
            result["measurements"].append({**definition.model_dump(exclude_none=True), "atoms": mapped})
        except ValueError as exc:
            result["errors"].append(f"{definition.label}: {exc}.")
    return result


def validate_source(settings: SimulationConfig, input_path: Path) -> dict[str, int]:
    """Validate selections before submission and pin their exact input identity."""
    if not settings.measurements:
        return {}
    source = load_physical(settings.dataset_id)[0]
    selected = sorted({i for definition in settings.measurements for i in definition.atoms})
    if any(i >= source.n_atoms for i in selected):
        raise ValueError("A tracked atom index is outside the source topology.")
    if settings.engine == "gromacs" and any(
        source.topology.atom(i).element == md.element.hydrogen for i in selected
    ):
        raise ValueError("GROMACS rebuilds input hydrogens in this setup. Track heavy-atom distances, angles or dihedrals, or choose OpenMM to retain explicit hydrogen and hydrogen-bond selections.")
    exact = md.load(str(input_path))
    mapping = map_atoms(source.topology, exact.topology, selected)
    # The canonical viewer and exact prepared PDB may differ by PDB rounding,
    # but a changed coordinate/identity must not silently become another atom.
    for index, target in mapping.items():
        if not np.allclose(source.xyz[0, index], exact.xyz[0, target], atol=2e-4, rtol=0):
            raise ValueError("The tracked source atom does not match the exact simulation input. Reload the prepared structure and select the measurement again.")
    for definition in settings.measurements:
        _measure(source, MeasurementRequest(kind=definition.kind, atoms=definition.atoms), strict=True)
    return {str(index): target for index, target in mapping.items()}


def empty_snapshot(job: dict) -> dict:
    return {
        "job_id": job["id"], "status": job["status"],
        "source_dataset_id": job.get("config", {}).get("dataset_id"),
        "output_dataset_id": job.get("dataset_id"),
        "measurements": [
            {**definition, "output_atoms": None, "values": [], "times_ps": [],
             "unit": "°" if definition["kind"] in {"angle", "dihedral"} else "Å",
             "frame_errors": [], "warnings": []}
            for definition in job.get("config", {}).get("measurements", [])
        ],
        "warnings": [], "errors": [], "last_time_ps": None,
    }


def read_snapshot(folder: Path, job: dict) -> dict:
    """Bounded JSON read only; status comes from current job supervision."""
    result = empty_snapshot(job)
    path = folder / "measurements.json"
    if path.exists():
        # Twelve measurements × at most 10,000 saved frames, including geometry
        # errors and H-bond angles. Never trust an unbounded local snapshot.
        if path.stat().st_size > 32 * 1024 * 1024:
            result["errors"] = ["The saved measurement snapshot exceeds the local size limit."]
        else:
            try:
                saved = json.loads(path.read_text())
                if saved.get("job_id") != job["id"]:
                    raise ValueError("Measurement snapshot belongs to another job.")
                result.update(saved)
            except (ValueError, OSError) as exc:
                result["errors"] = [f"Could not read the saved measurement snapshot: {exc}"]
    result.update(status=job["status"], output_dataset_id=job.get("dataset_id"))
    return result


class LiveMeasurements:
    def __init__(self, folder: Path, job: dict, input_state: dict, topology: md.Topology):
        self.folder, self.topology = folder, topology
        self.snapshot = empty_snapshot(job)
        self.steps: list[int] = []
        self.last_flush = 0.0
        self.xtc_frames = 0
        self.xtc_size = -1
        self.xtc_generation = None
        self.xtc_signature = None
        self.active = bool(self.snapshot["measurements"])
        if not self.active:
            return
        source_mapping = input_state.get("measurement_source_atoms", {})
        selected = sorted({i for item in self.snapshot["measurements"] for i in item["atoms"]})
        if any(str(index) not in source_mapping for index in selected):
            raise ValueError("The saved live measurement atom identities are missing; start a new simulation from the source structure.")
        input_topology = md.load_topology(str(folder / "input.pdb"))
        input_indices = [source_mapping[str(index)] for index in selected]
        output_mapping = map_atoms(input_topology, topology, input_indices)
        verified_output = input_state.get("measurement_verified_output_atoms")
        if verified_output is not None and any(verified_output.get(str(i)) != output_mapping[i] for i in input_indices):
            raise ValueError("Tracked atom ordering differs from the verified native preparation; no replacement atom was selected.")
        for item in self.snapshot["measurements"]:
            item["output_atoms"] = [output_mapping[source_mapping[str(index)]] for index in item["atoms"]]
        # Persist the verified correspondence with the job outputs for inspection.
        atomic_json(folder / "measurement-atoms.json", {
            "source_dataset_id": job["config"]["dataset_id"],
            "source_to_input": source_mapping,
            "source_to_output": {str(i): output_mapping[source_mapping[str(i)]] for i in selected},
            "method": "Unique atom name and element within verified input residue order; no coordinate-nearest or numeric-index substitution.",
        })
        self.flush(force=True)

    def stop(self, error: Exception):
        """An unavailable monitor never cancels otherwise valid MD."""
        self.snapshot["errors"] = [f"Live measurements stopped: {error}"]
        self.active = False
        self.flush(force=True)

    def add_frame(self, xyz, time_ps: float, step: int, box=None):
        if not self.active:
            return
        if self.steps and step <= self.steps[-1]:
            if step == self.steps[-1]:
                return
            raise ValueError("Saved production frames are not in increasing step order.")
        if len(self.steps) >= 10000:
            raise ValueError("Live measurements exceed the saved-frame limit.")
        coordinates = np.asarray(xyz, dtype=np.float32)
        if coordinates.shape != (self.topology.n_atoms, 3) or not np.isfinite(coordinates).all() or not np.isfinite(time_ps):
            raise ValueError("A saved production frame has invalid coordinates or time.")
        if self.steps and time_ps <= self.snapshot["last_time_ps"]:
            raise ValueError("Saved production frames are not in increasing time order.")
        frame = md.Trajectory(coordinates[None], self.topology, time=[time_ps])
        if box is not None:
            cell = np.asarray(box, dtype=np.float32)
            if cell.shape != (3, 3) or not np.isfinite(cell).all() or np.linalg.det(cell) <= 0:
                raise ValueError("A saved production frame has an invalid periodic box.")
            frame.unitcell_vectors = cell[None]
        # Compute every result before appending any, keeping all series aligned.
        results = [_measure(frame, MeasurementRequest(kind=item["kind"], atoms=item["output_atoms"]), strict=False)
                   for item in self.snapshot["measurements"]]
        for item, result in zip(self.snapshot["measurements"], results):
            for key in ("values", "times_ps", "frame_errors", "angle_values", "geometry_passes"):
                if key in result:
                    item.setdefault(key, []).extend(result[key])
            item["warnings"] = result["warnings"]
            if "geometry_passes" in item:
                valid = [value for value in item["geometry_passes"] if value is not None]
                item["occupancy"] = sum(valid) / len(valid) if valid else None
        self.steps.append(int(step))
        self.snapshot["last_time_ps"] = float(time_ps)

    def flush(self, *, force=False):
        if self.snapshot["measurements"] and (force or (self.active and time.monotonic() - self.last_flush >= 1)):
            atomic_json(self.folder / "measurements.json", self.snapshot)
            self.last_flush = time.monotonic()

    def read_xtc(self, *, final=False):
        """Read only new complete frames; native frame seeking skips old data.

        A growing XTC may stop within a frame. A fresh handle on the next file
        growth retries that same frame. Only fully decoded frames advance the
        cursor; no fabricated interpolation or partially decoded coordinates.
        """
        if not self.active:
            return
        path = self.folder / "production.xtc"
        if not path.is_file():
            return
        info = path.stat()
        size = info.st_size
        generation = (size, info.st_mtime_ns, info.st_ino)
        if not size:
            if self.xtc_frames:
                self._reset_xtc()
                self.flush(force=True)
            return
        if not final and generation == self.xtc_generation:
            return
        from mdtraj.formats import XTCTrajectoryFile
        reached_limit = False
        def signature(xyz, time_ps, step, box):
            digest = hashlib.sha256()
            for value in (xyz, np.asarray([time_ps], dtype=np.float64), np.asarray([step], dtype=np.int64), box):
                digest.update(np.asarray(value).tobytes())
            return digest.hexdigest()
        try:
            with XTCTrajectoryFile(str(path)) as reader:
                if self.xtc_frames:
                    # Native resume can truncate and regenerate the trailing
                    # frame prefix. Verify its boundary before advancing, even
                    # when truncation and regrowth happened between polls.
                    try:
                        reader.seek(self.xtc_frames - 1)
                        xyz, times, steps, boxes = reader.read(1)
                        same = len(xyz) == 1 and signature(xyz[0], times[0], steps[0], boxes[0]) == self.xtc_signature
                    except (OSError, RuntimeError):
                        same = False
                    if not same:
                        self._reset_xtc()
                        reader.seek(0)
                # A bounded batch during dynamics; finalization drains the tail.
                for _ in range(10000 if final else 256):
                    try:
                        xyz, times, steps, boxes = reader.read(1)
                    except (RuntimeError, OSError):
                        if final:
                            raise
                        break
                    if not len(xyz):
                        break
                    self.add_frame(xyz[0], float(times[0]), int(steps[0]), boxes[0])
                    self.xtc_frames += 1
                    self.xtc_signature = signature(xyz[0], times[0], steps[0], boxes[0])
                else:
                    reached_limit = True
        except (RuntimeError, OSError):
            if final:
                raise ValueError("The final production trajectory contains an unreadable frame; live measurements were not completed.")
        self.xtc_size = -1 if reached_limit else size
        self.xtc_generation = None if reached_limit else generation
        self.flush(force=final)

    def _reset_xtc(self):
        self.steps = []
        self.xtc_frames = 0
        self.xtc_size = -1
        self.xtc_generation = None
        self.xtc_signature = None
        self.snapshot["last_time_ps"] = None
        for item in self.snapshot["measurements"]:
            for key in ("values", "times_ps", "frame_errors", "angle_values", "geometry_passes"):
                if key in item:
                    item[key] = []
            item.pop("occupancy", None)

    def finish(self, trajectory: md.Trajectory):
        if not self.active:
            return
        if len(self.steps) != trajectory.n_frames or not np.allclose(
            self.snapshot["measurements"][0]["times_ps"], trajectory.time, atol=1e-7, rtol=1e-6
        ):
            raise ValueError("Live measurements do not match the final saved trajectory frames.")
        self.flush(force=True)
