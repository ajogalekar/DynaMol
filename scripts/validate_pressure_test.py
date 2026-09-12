#!/usr/bin/env python3
"""Audit saved DynaMol simulation artifacts; never starts or resumes dynamics.

Run in a fresh process so DYNAMOL_DATA_DIR is set before backend imports.
Checks establish software/artifact consistency and screen gross numerical
problems. They do not establish equilibration, native poses, or convergence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np


def series_comparison(actual, expected, *, circular=False, tolerance=1e-4):
    """Compare nullable geometry; torsions use shortest signed angular distance."""
    if len(actual) != len(expected):
        return {"passed": False, "reason": "frame count mismatch", "actual": len(actual), "expected": len(expected)}
    mask_a = np.asarray([v is None for v in actual])
    mask_b = np.asarray([v is None for v in expected])
    if not np.array_equal(mask_a, mask_b):
        return {"passed": False, "reason": "undefined frame masks differ"}
    a = np.asarray([v for v in actual if v is not None], dtype=float)
    b = np.asarray([v for v in expected if v is not None], dtype=float)
    difference = (a - b + 180) % 360 - 180 if circular else a - b
    maximum = float(np.max(np.abs(difference))) if len(difference) else 0.0
    return {"passed": bool(np.isfinite(a).all() and np.isfinite(b).all() and maximum <= tolerance),
            "max_absolute_difference": maximum, "tolerance": tolerance,
            "valid_frames": int(len(a)), "undefined_frames": int(mask_a.sum()), "circular": circular}


def protein_bond_screen(trajectory, protein_indices, *, maximum_nm=.30):
    import mdtraj as md
    protein = set(protein_indices)
    pairs = [(a.index, b.index) for a, b in trajectory.topology.bonds
             if a.index in protein and b.index in protein and a.element is not None and b.element is not None
             and a.element.symbol not in {"H", "D"} and b.element.symbol not in {"H", "D"}]
    if not pairs:
        return {"passed": False, "reason": "No protein heavy-atom bonds available to inspect"}
    distances = md.compute_distances(trajectory, pairs, periodic=trajectory.unitcell_vectors is not None)
    frame, bond = np.unravel_index(np.argmax(distances), distances.shape)
    return {"passed": bool(np.isfinite(distances).all() and np.all(distances > .03) and distances.max() < maximum_nm),
            "bonds": len(pairs), "frames": trajectory.n_frames, "maximum_angstrom": float(distances[frame, bond] * 10),
            "maximum_frame": int(frame), "maximum_pair": list(pairs[bond]),
            "minimum_angstrom": float(distances.min() * 10), "threshold_angstrom": maximum_nm * 10,
            "note": "Gross bond-length screen with minimum-image distances; not a force-field or structural-quality validation."}


def validate(data_root: Path, job_id: str):
    data_root = data_root.resolve()
    if not job_id or Path(job_id).name != job_id or job_id in {".", ".."}:
        raise ValueError("Job ID must be one directory name")
    os.environ["DYNAMOL_DATA_DIR"] = str(data_root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from backend import config
    from backend.storage import load_physical
    from backend.analysis import _measure
    from backend.models import MeasurementRequest
    import mdtraj as md
    if config.DATA_ROOT != data_root:
        raise ValueError("Run validator in a fresh Python process; backend was already configured for another data root")
    folder = data_root / "jobs" / job_id
    report = {"job_id": job_id, "data_root": str(data_root), "checks": {}, "evidence": {},
              "scope": "Saved-artifact consistency and gross numerical health for exploratory NVT software pressure testing; no equilibration, convergence, biological stability, or pose-accuracy claim.",
              "limits": ["Temperature bands and bond limits are gross alarms, not statistical convergence tests.",
                         "NVT temperature fluctuations and total-energy changes are expected; NVE energy conservation is not tested.",
                         "Live/posthoc comparison uses the app's geometry implementation; its independent analytic controls are tested separately."]}

    def check(name, passed, **details):
        report["checks"][name] = {"passed": bool(passed), **details}

    def digest(path):
        with path.open("rb") as handle:
            return hashlib.file_digest(handle, "sha256").hexdigest()

    def evidence(path):
        report["evidence"][str(path.relative_to(data_root))] = {"sha256": digest(path), "bytes": path.stat().st_size}

    def read_json(path):
        value = json.loads(path.read_text())
        evidence(path)
        return value

    try:
        status = read_json(folder / "status.json")
        settings = read_json(folder / "config.json")
        provenance = read_json(folder / "provenance.json")
        steps = round(settings["duration_ps"] * 1000 / settings["timestep_fs"])
        expected_frames = math.ceil(steps / settings["report_interval"]) + 1
        final_time = steps * settings["timestep_fs"] / 1000
        report.update(engine=settings["engine"], config=settings, elapsed_seconds=status.get("elapsed_seconds"),
                      versions=provenance.get("versions"), cpu_threads=provenance.get("cpu_threads"),
                      expected={"production_steps": steps, "frames": expected_frames, "final_time_ps": final_time})
        check("completed", status.get("status") == "completed" and status.get("completed_steps") == steps and status.get("total_steps") == steps,
              status=status.get("status"), completed_steps=status.get("completed_steps"), error=status.get("error"))
        dataset_id = status.get("dataset_id")
        if not dataset_id or Path(dataset_id).name != dataset_id:
            raise ValueError("Completed output dataset ID is absent or invalid")
        output = data_root / "datasets" / dataset_id
        metadata = read_json(output / "metadata.json")
        physical = load_physical(dataset_id)
        evidence(output / "physical.npz")
        evidence(output / "topology.pdb")
        check("physical_trajectory", physical.n_frames == expected_frames and np.isfinite(physical.xyz).all()
              and np.isfinite(physical.time).all() and np.all(np.diff(physical.time) > 0)
              and abs(float(physical.time[0])) < 1e-6 and abs(float(physical.time[-1]) - final_time) < 1e-4
              and physical.n_atoms == metadata["n_atoms"] and metadata["n_frames"] == physical.n_frames,
              dataset_id=dataset_id, atoms=physical.n_atoms, frames=physical.n_frames, final_time_ps=float(physical.time[-1]))
        moved = np.linalg.norm(physical.xyz[-1] - physical.xyz[0], axis=1)
        check("coordinates_evolve", bool(np.any(moved > 1e-5)), maximum_displacement_nm=float(moved.max()))
        explicit = settings["solvent"] == "explicit"
        boxes = physical.unitcell_vectors
        valid_box = boxes is not None and np.isfinite(boxes).all() and np.all(np.linalg.det(boxes) > 0)
        check("solvent_box", bool(valid_box and np.allclose(boxes, boxes[0], atol=1e-5, rtol=1e-5)) if explicit else boxes is None,
              explicit=explicit, has_box=boxes is not None, fixed_volume=True)
        # Display coordinates may be imaged/aligned. They must never substitute
        # for physical.npz in any numerical check below.
        check("display_payload", (output / "coordinates.bin").stat().st_size == physical.n_frames * physical.n_atoms * 12)
        report["checks"]["protein_heavy_bonds"] = protein_bond_screen(physical, [a["index"] for a in metadata["atoms"] if a["category"] == "protein"])
        from backend.stereo_monitor import StereoMonitor
        stereo = StereoMonitor(physical.topology.to_openmm(), physical.xyz[0], folder / 'input.pdb')
        stereo_failures = []
        for index, coordinates in enumerate(physical.xyz):
            result = stereo.check(coordinates, boxes[index] if boxes is not None else None, f'Frame {index}')
            if not result['passed']:
                stereo_failures.append(result)
        check('residue_stereochemistry', len(stereo.centers) > 0 and not stereo_failures,
              centers_per_frame=len(stereo.centers), frames=physical.n_frames, failures=stereo_failures)

        native_path = folder / ("trajectory.xtc" if settings["engine"] == "openmm" else "production.xtc")
        native_topology = folder / ("prepared.pdb" if settings["engine"] == "openmm" else "system.gro")
        native = md.load(str(native_path), top=str(native_topology))
        evidence(native_path)
        evidence(native_topology)
        compatible_shape = native.xyz.shape == physical.xyz.shape
        delta = float(np.max(np.abs(native.xyz - physical.xyz))) if compatible_shape else None
        tolerance = 6e-4 if settings["engine"] == "openmm" else 1e-7  # XTC writes at 0.001 nm precision.
        native_names = [(a.name, a.element.symbol if a.element else None) for a in native.topology.atoms]
        physical_names = [(a.name, a.element.symbol if a.element else None) for a in physical.topology.atoms]
        check("native_trajectory", compatible_shape and native_names == physical_names and np.isfinite(native.xyz).all()
              and delta <= tolerance and np.allclose(native.time, physical.time, atol=1e-4, rtol=1e-6),
              max_coordinate_difference_nm=delta, tolerance_nm=tolerance, frames=native.n_frames, atoms=native.n_atoms)
        if explicit:
            check("native_box", native.unitcell_vectors is not None and valid_box and np.allclose(native.unitcell_vectors, boxes, atol=1e-5, rtol=1e-5))

        with (folder / "energies.csv").open() as handle:
            energies = list(csv.DictReader(handle))
        evidence(folder / "energies.csv")
        columns = ["step", "time_ps", "potential_kj_mol", "kinetic_kj_mol", "temperature_k"]
        numeric = np.asarray([[float(row[key]) for key in columns] for row in energies])
        if numeric.ndim != 2 or not len(numeric):
            raise ValueError("Native energy diagnostics are empty")
        aligned = len(numeric) == physical.n_frames and np.allclose(numeric[:, 1], physical.time, atol=1e-4, rtol=1e-6)
        temperatures = numeric[:, 4]
        target = settings["temperature_k"]
        check("energy_diagnostics", np.isfinite(numeric).all() and aligned and np.all(numeric[:, 3] > 0)
              and np.all(temperatures > 0) and np.all(np.diff(numeric[:, 0]) > 0)
              and np.allclose(numeric[:, 0], np.round(numeric[:, 1] * 1000 / settings["timestep_fs"]), atol=1e-4),
              rows=len(numeric), min_potential_kj_mol=float(numeric[:, 2].min()), max_potential_kj_mol=float(numeric[:, 2].max()))
        check("gross_temperature_screen", .5 * target <= temperatures.mean() <= 1.5 * target and temperatures.max() < 10 * target,
              target_k=target, mean_k=float(temperatures.mean()), min_k=float(temperatures.min()), max_k=float(temperatures.max()),
              note="Broad numerical alarm only; no claim of equilibration or expected distribution.")

        if settings["engine"] == "openmm":
            import openmm as mm
            from openmm import unit
            system = mm.XmlSerializer.deserialize((folder / "system.xml").read_text())
            integrator = mm.XmlSerializer.deserialize((folder / "integrator.xml").read_text())
            evidence(folder / "system.xml")
            evidence(folder / "integrator.xml")
            check("native_system", system.getNumParticles() == physical.n_atoms
                  and abs(integrator.getStepSize().value_in_unit(unit.femtosecond) - settings["timestep_fs"]) < 1e-9
                  and abs(integrator.getTemperature().value_in_unit(unit.kelvin) - target) < 1e-9
                  and integrator.getRandomNumberSeed() == settings["seed"], particles=system.getNumParticles())
            dof = sum(3 for i in range(system.getNumParticles()) if system.getParticleMass(i) > 0 * unit.dalton)
            dof -= sum(system.getParticleMass(a) > 0 * unit.dalton and system.getParticleMass(b) > 0 * unit.dalton for a, b, _ in (system.getConstraintParameters(i) for i in range(system.getNumConstraints())))
            dof -= 3 if any(isinstance(f, mm.CMMotionRemover) for f in system.getForces()) else 0
            calculated = 2 * numeric[:, 3] / (.00831446261815324 * dof)
            check("temperature_from_kinetic_energy", np.allclose(calculated, temperatures, atol=1e-5, rtol=1e-7), degrees_of_freedom=dof)
            recovery = read_json(folder / "recovery.json")
            checkpoint = recovery["checkpoint"]
            committed = checkpoint["frames"]
            check("committed_checkpoint", checkpoint["step"] == steps and len(committed) == expected_frames
                  and digest(folder / "checkpoints" / checkpoint["directory"] / "checkpoint.chk") == checkpoint["sha256"])
            prefix_ok = len(committed) == physical.n_frames
            for i, record in enumerate(committed):
                path = folder / "frames" / record["name"]
                with np.load(path, allow_pickle=False) as frame:
                    prefix_ok = prefix_ok and digest(path) == record["sha256"] and np.array_equal(frame["xyz"], physical.xyz[i])
                    prefix_ok = prefix_ok and abs(float(frame["time"]) - physical.time[i]) < 1e-7 and int(frame["step"]) == numeric[i, 0]
                    prefix_ok = prefix_ok and np.allclose([frame["potential"], frame["kinetic"], frame["temperature"]], numeric[i, 2:], atol=1e-7, rtol=0)
            check("committed_physical_frames", prefix_ok)
            dcd = md.load(str(folder / "trajectory.dcd"), top=str(native_topology))
            times = np.loadtxt(folder / "frame_times_ps.csv", skiprows=1, ndmin=1)
            check("dcd_with_time_sidecar", dcd.xyz.shape == physical.xyz.shape and np.allclose(dcd.xyz, physical.xyz, atol=2e-6, rtol=1e-6)
                  and np.allclose(times, physical.time, atol=1e-7, rtol=1e-7))

        definitions = settings.get("measurements", [])
        if definitions:
            snapshot = read_json(folder / "measurements.json")
            identity = read_json(folder / "measurement-identity.json")
            input_state = read_json(folder / "input-state.json")
            items = {item["id"]: item for item in snapshot.get("measurements", [])}
            check("live_measurement_coverage", not snapshot.get("errors") and identity.get("verified")
                  and set(items) == {item["id"] for item in definitions}, errors=snapshot.get("errors"))
            for definition in definitions:
                item = items.get(definition["id"])
                if not item or not item.get("output_atoms"):
                    check("measurement:" + definition["id"], False, reason="Missing verified output atom identities")
                    continue
                recomputed = _measure(physical, MeasurementRequest(kind=item["kind"], atoms=item["output_atoms"]), strict=False)
                comparison = series_comparison(item.get("values", []), recomputed["values"], circular=item["kind"] == "dihedral")
                verified_atoms = [identity["input_to_native"][str(input_state["measurement_source_atoms"][str(index)])] for index in definition["atoms"]]
                comparison["identities_match"] = item["kind"] == definition["kind"] and item.get("atoms") == definition["atoms"] and item["output_atoms"] == verified_atoms
                comparison["times_match"] = bool(len(item.get("times_ps", [])) == physical.n_frames and np.allclose(item["times_ps"], physical.time, atol=1e-4, rtol=1e-6))
                comparison["undefined_reasons_match"] = item.get("frame_errors") == recomputed.get("frame_errors")
                comparison["passed"] &= comparison["times_match"] and comparison["undefined_reasons_match"] and comparison["identities_match"] and comparison.get("valid_frames", 0) > 0
                if item["kind"] == "hbond":
                    angle = series_comparison(item.get("angle_values", []), recomputed["angle_values"])
                    valid = [v for v in recomputed["geometry_passes"] if v is not None]
                    occupancy = sum(valid) / len(valid) if valid else None
                    same_occupancy = item.get("occupancy") is None if occupancy is None else item.get("occupancy") is not None and abs(item["occupancy"] - occupancy) < 1e-8
                    comparison.update(angle_comparison=angle, occupancy=occupancy, occupancy_matches=same_occupancy)
                    comparison["passed"] &= angle["passed"] and same_occupancy and item.get("geometry_passes") == recomputed["geometry_passes"]
                report["checks"]["measurement:" + definition["id"]] = comparison
        else:
            report["limits"].append("No live measurements were configured; this run does not exercise live tracking.")
        hashes = provenance.get("outputs", {})
        bad = []
        for name, expected_hash in hashes.items():
            path = (folder / name).resolve()
            if not path.is_relative_to(folder.resolve()) or not path.is_file() or digest(path) != expected_hash:
                bad.append(name)
        check("recorded_artifact_hashes", bool(hashes) and not bad, checked=len(hashes), mismatches=bad)
    except Exception as exc:
        check("artifact_validation_exception", False, error=f"{type(exc).__name__}: {exc}")
    report["passed"] = bool(report["checks"]) and all(value["passed"] for value in report["checks"].values())
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--job", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = validate(args.data_root, args.job)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def json_safe(value):
        if isinstance(value, dict):
            return {key: json_safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [json_safe(item) for item in value]
        if isinstance(value, (float, np.floating)) and not math.isfinite(value):
            return None
        if isinstance(value, np.generic):
            return value.item()
        return value
    args.output.write_text(json.dumps(json_safe(report), indent=2, allow_nan=False) + "\n")
    print(json.dumps({"job": args.job, "passed": report["passed"], "checks": len(report["checks"]),
                      "failed": [key for key, value in report["checks"].items() if not value["passed"]], "report": str(args.output)}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
