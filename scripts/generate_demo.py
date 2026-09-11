#!/usr/bin/env python3
"""Download public ubiquitin coordinates and run a real 2 ps OpenMM demonstration.

Explicit protein selection removes crystallographic waters from the downloaded
PDB before creating the implicit-solvent starting dataset. User uploads never
undergo this selection silently. At most two CPU threads by default.
"""
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import mdtraj as md
from backend import config, jobs, storage
from backend.models import SimulationConfig


def main():
    existing = config.DATASETS_DIR / "demo" / "metadata.json"
    if existing.exists() and storage.get_dataset("demo")["n_frames"] > 1 and "--force" not in sys.argv:
        print("Real MD demo already exists. Use --force to regenerate.")
        return
    source_dir = config.DATA_ROOT / "sources"
    source_dir.mkdir(exist_ok=True)
    source = source_dir / "1UBQ.pdb"
    url = "https://files.rcsb.org/download/1UBQ.pdb"
    if not source.exists():
        with urllib.request.urlopen(url, timeout=30) as response:
            source.write_bytes(response.read())
    original = md.load(str(source))
    protein = original.atom_slice(original.topology.select("protein"))
    removed = original.n_atoms - protein.n_atoms
    print(f"Demo preparation: retained {protein.n_atoms} protein atoms; explicitly excluded {removed} crystallographic water/nonprotein atoms from 1UBQ.", flush=True)
    provenance = {"pdb_id": "1UBQ", "source_url": url, "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "selection": "protein", "excluded_nonprotein_atoms": removed, "purpose": "Software demonstration only; not a converged MD study."}
    initial = storage.save_dataset(protein, "Ubiquitin · starting structure", "RCSB PDB 1UBQ", "Experimental ubiquitin structure; protein-only selection for an implicit-solvent demonstration. Crystallographic waters explicitly excluded.", warnings=[f"Demo preparation explicitly selected protein and excluded {removed} crystallographic water/nonprotein atoms."], dataset_id="ubiquitin-start", provenance=provenance)
    settings = SimulationConfig(dataset_id=initial["id"], engine="openmm", name="Ubiquitin · in motion", duration_ps=2, timestep_fs=2, report_interval=10, equilibration_steps=100, temperature_k=300, friction_ps=1, solvent="implicit", seed=2026, minimize=True)
    job = jobs.submit_job(settings)
    print(f"Started real OpenMM demo job {job['id']}", flush=True)
    previous = ""
    while True:
        job = jobs.get_job(job["id"])
        summary = f"{job['status']}: {job['stage']} — {job['completed_steps']}/{job['total_steps']}"
        if summary != previous:
            print(summary, flush=True)
            previous = summary
        if job["status"] in {"completed", "failed", "cancelled", "interrupted"}:
            break
        time.sleep(1)
    if job["status"] != "completed":
        raise RuntimeError(job.get("error", "Demo job did not finish."))
    trajectory = storage.load_physical(job["dataset_id"])
    job_provenance = json.loads((config.JOBS_DIR / job["id"] / "provenance.json").read_text())
    combined = {**job_provenance, "starting_structure": provenance, "job_id": job["id"]}
    meta = storage.save_dataset(trajectory, "Ubiquitin · in motion", "Real OpenMM MD · PDB 1UBQ", "A real 2 ps ubiquitin trajectory at 300 K. Amber ff14SB + GBn2 implicit solvent, 2 fs Langevin-middle steps, seed 2026. An exploration of motion, not an equilibrated scientific study.", warnings=["Only 2 ps of production after 0.2 ps initial relaxation: this demonstration does not establish equilibration or convergence.", f"The starting model was explicitly selected from PDB 1UBQ; {removed} crystallographic water atoms were excluded for implicit solvent."], dataset_id="demo", provenance=combined)
    print(json.dumps({"demo": meta["id"], "frames": meta["n_frames"], "atoms": meta["n_atoms"], "job_id": job["id"]}, indent=2))


if __name__ == "__main__":
    main()
