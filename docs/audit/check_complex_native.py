"""Bounded live API regression, using an explicitly selected prepared cluster.

Run from the project root with .venv/bin/python docs/audit/check_complex_native.py DATASET_ID.
This checks software/state continuity, not equilibration or scientific validity.
"""
import csv
import hashlib
import io
import json
import math
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
from openmm import XmlSerializer, app, unit, NonbondedForce

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from backend import config, storage
from backend.prepared_system import load_prepared_forcefield, ligand_parameter_files

BASE = "http://127.0.0.1:8765/api/"
REPORT = ROOT / "docs" / "audit" / "complex-native-run.json"
report = {"purpose": "Software continuity and finite short CPU dynamics; no equilibrium, affinity, or metal-model accuracy claim.", "checks": {}}


def request(path, body=None):
    req = urllib.request.Request(BASE + path, data=json.dumps(body).encode() if body is not None else None,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.load(response)


def record():
    REPORT.write_text(json.dumps(report, indent=2) + "\n")


def wait_job(job):
    started, previous = time.monotonic(), None
    while time.monotonic() - started < 900:
        current = request("jobs/" + job["id"])
        status = (current["status"], current["stage"])
        if status != previous:
            print(job["id"], *status, flush=True)
            previous = status
        if current["status"] not in {"queued", "running", "cancelling"}:
            assert current["status"] == "completed", current.get("error", current)
            return current
        time.sleep(.5)
    raise TimeoutError("Native validation job exceeded 15 minutes; inspect its live status before proceeding.")


def verify_bundle(dataset):
    paths = ligand_parameter_files(storage.dataset_dir(dataset["id"]), dataset["preparation"], required=True)
    hashes = sorted(hashlib.sha256(path.read_bytes()).hexdigest() for path in paths)
    assert hashes == report["parameter_hashes"]
    return hashes


def actual_ligand_charges(dataset, system=None):
    folder = storage.dataset_dir(dataset["id"])
    pdb = app.PDBFile(str(folder / "prepared.pdb"))
    if system is None:
        ff, _ = load_prepared_forcefield(folder, dataset["preparation"])
        system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    force = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    results = []
    for ligand in dataset["preparation"]["ligand_parameters"]["ligands"]:
        residue = next(r for r in pdb.topology.residues() if ":".join((r.chain.id, r.id, r.insertionCode.strip(), r.name)) == ligand["key"])
        actual = [force.getParticleParameters(atom.index)[0].value_in_unit(unit.elementary_charge) for atom in residue.atoms()]
        np.testing.assert_allclose(actual, ligand["charges_e"], atol=1e-6, rtol=0)
        assert abs(sum(actual) - ligand["formal_charge"]) < 1e-4
        results.append({"key": ligand["key"], "atoms": len(actual), "charge_e": sum(actual)})
    return results


def main(dataset_id):
    parent = request("datasets/" + dataset_id)
    report["fixture_dataset_id"] = parent["id"]
    report["parameter_hashes"] = sorted(entry["sha256"] for entry in parent["preparation"]["ligand_parameters"]["files"])
    report["prepared_charges"] = actual_ligand_charges(parent)
    job = request("datasets/" + parent["id"] + "/solvate", {"padding_nm": 1, "seed": 2026, "ph": 7})
    report["solvation_job_id"] = job["id"]; record()
    job = wait_job(job)
    solvated = request("datasets/" + job["dataset_id"])
    report["solvated_dataset_id"] = solvated["id"]
    original = app.PDBFile(str(storage.dataset_dir(parent["id"]) / "prepared.pdb"))
    box = app.PDBFile(str(storage.dataset_dir(solvated["id"]) / "prepared.pdb"))
    delta = np.asarray(original.positions.value_in_unit(unit.angstrom)) - np.asarray(box.positions.value_in_unit(unit.angstrom))[:original.topology.getNumAtoms()]
    report["checks"]["solvation_solute_max_displacement_angstrom"] = float(np.max(np.linalg.norm(delta, axis=1)))
    assert report["checks"]["solvation_solute_max_displacement_angstrom"] < 1e-6
    verify_bundle(solvated)
    report["solvated_charges"] = actual_ligand_charges(solvated)
    report["solvated_atoms"] = solvated["n_atoms"]
    assert solvated["has_unitcell"] and solvated["n_atoms"] > parent["n_atoms"]
    record()
    settings = {"dataset_id": solvated["id"], "name": "e2e-complex-validation · native OpenMM", "engine": "openmm", "solvent": "explicit", "duration_ps": .02,
                "temperature_k": 300, "timestep_fs": 2, "report_interval": 2, "seed": 2026, "minimize": True, "equilibration_steps": 10}
    job = request("jobs", settings)
    report["simulation_job_id"] = job["id"]; report["simulation_settings"] = settings; record()
    job = wait_job(job)
    result = request("datasets/" + job["dataset_id"])
    verify_bundle(result)
    report["trajectory_dataset_id"] = result["id"]
    native = config.JOBS_DIR / job["id"]
    system = XmlSerializer.deserialize((native / "system.xml").read_text())
    report["native_system_charges"] = actual_ligand_charges(result, system)
    with np.load(storage.dataset_dir(result["id"]) / "physical.npz") as physical:
        assert all(np.isfinite(physical[key]).all() for key in physical.files)
    values = list(csv.reader(io.StringIO((native / "energies.csv").read_text())))
    assert len(values) > 1 and all(math.isfinite(float(value)) for row in values[1:] for value in row)
    assert result["n_frames"] == 6 and result["n_atoms"] == solvated["n_atoms"]
    report["checks"].update(parameter_hashes_preserved=True, charges_preserved=True, trajectory_finite=True, energy_records_finite=True,
                            saved_frames=result["n_frames"], production_steps=job["completed_steps"], production_ps=.02)
    report["native_files_sha256"] = {name: hashlib.sha256((native / name).read_bytes()).hexdigest() for name in ("system.xml", "prepared.pdb", "energies.csv", "config.json")}
    report["passed"] = True
    record()
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    try:
        main(sys.argv[1])
    except Exception as exc:
        report.update(passed=False, error=str(exc)); record()
        raise
