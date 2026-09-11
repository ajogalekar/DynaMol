"""Subprocess job supervision, atomic disk status and restart-aware cancellation."""
import copy
import datetime as dt
import importlib.metadata
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
from pathlib import Path

from . import config
from .models import SimulationConfig
from .storage import atomic_json, dataset_dir, get_dataset, safe_id

_lock = threading.Lock()
_health_lock = threading.Lock()
_health_cache: tuple[float, dict] | None = None
HEALTH_CACHE_SECONDS = 30.0


def gromacs_executable() -> str | None:
    configured = os.environ.get("DYNAMOL_GMX")
    candidates = [configured, shutil.which("gmx"), shutil.which("gmx_mpi"), str(config.ROOT / ".tools" / "gromacs" / "bin" / "gmx"), str(config.ROOT / ".gromacs" / "bin" / "gmx")]
    return next((p for p in candidates if p and Path(p).is_file() and os.access(p, os.X_OK)), None)


def health() -> dict:
    """Share engine availability probes across browser polling for thirty seconds."""
    global _health_cache
    with _health_lock:
        now = time.monotonic()
        if _health_cache is None or now - _health_cache[0] >= HEALTH_CACHE_SECONDS:
            _health_cache = (now, _probe_health())
        return copy.deepcopy(_health_cache[1])


def _probe_health() -> dict:
    try:
        version = importlib.metadata.version("openmm")
        openmm = {"id": "openmm", "name": "OpenMM", "available": True, "version": version, "message": "CPU engine ready · Amber ff14SB + GBn2 or TIP3P · standard proteins"}
    except importlib.metadata.PackageNotFoundError:
        openmm = {"id": "openmm", "name": "OpenMM", "available": False, "version": None, "message": "Install the Python backend dependencies to enable OpenMM."}
    executable = gromacs_executable()
    gromacs = {"id": "gromacs", "name": "GROMACS", "available": False, "version": None, "message": "GROMACS is not installed. Install gmx and put it on PATH, or set DYNAMOL_GMX."}
    if executable:
        try:
            probe = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=10)
            if probe.returncode == 0:
                line = next((line.split(":", 1)[1].strip() for line in (probe.stdout + probe.stderr).splitlines() if line.startswith("GROMACS version:")), "installed")
                gromacs.update(available=True, version=line, message="CPU engine ready · Amber99SB-ILDN / TIP3P · explicit water")
            else:
                gromacs["message"] = "GROMACS executable found but its version probe failed. Check runtime libraries."
        except (OSError, subprocess.TimeoutExpired):
            gromacs["message"] = "GROMACS executable found but its version probe failed."
    return {"status": "ok", "engines": [openmm, gromacs]}


def get_job(job_id: str) -> dict:
    folder = config.JOBS_DIR / safe_id(job_id)
    path = folder / "status.json"
    if not path.exists():
        raise FileNotFoundError(f"Job '{job_id}' was not found.")
    job = json.loads(path.read_text())
    if job["status"] in {"queued", "running", "cancelling"}:
        job["elapsed_seconds"] = round((dt.datetime.now(dt.timezone.utc) - dt.datetime.fromisoformat(job["created_at"])).total_seconds(), 2)
        pid = job.get("worker_pid")
        if pid:
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                job.update(status="interrupted", stage="Interrupted", error="The worker stopped before completing. Outputs and logs have been retained.")
                atomic_json(path, job)
    return job


def list_jobs() -> list[dict]:
    result = []
    for path in sorted(config.JOBS_DIR.glob("*/status.json"), key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            result.append(get_job(path.parent.name))
        except (OSError, ValueError):
            pass
    return result


def submit_job(settings: SimulationConfig) -> dict:
    metadata = get_dataset(settings.dataset_id)
    if any(atom["category"] == "nucleic" for atom in metadata["atoms"]):
        raise ValueError("This simulation preset supports standard proteins only. RNA/DNA can be viewed and analyzed, but nucleic-acid and protein–nucleic-acid simulations require a separately parameterized workflow. No atoms were removed.")
    standard = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN"}
    unsupported = sorted({atom["residue"] for atom in metadata["atoms"] if atom["category"] == "ligands" or (atom["category"] == "protein" and atom["residue"] not in standard)})
    if unsupported:
        raise ValueError("This protein simulation workflow does not parameterize ligands or nonstandard residues: " + ", ".join(unsupported) + ". Prepare a supported protein-only structure; no atoms were removed.")
    if not any(atom["category"] == "protein" for atom in metadata["atoms"]):
        raise ValueError("Choose a structure containing a standard amino-acid protein.")
    if settings.solvent == "implicit" and any(atom["category"] in {"water", "ions"} for atom in metadata["atoms"]):
        raise ValueError("Implicit solvent requires a protein-only input. This structure has explicit waters or ions; choose explicit solvent or upload a prepared protein-only structure. No atoms were removed.")
    engine = next(engine for engine in health()["engines"] if engine["id"] == settings.engine)
    if not engine["available"]:
        raise ValueError(engine["message"])
    with _lock:
        if any(job["status"] in {"queued", "running", "cancelling"} for job in list_jobs()):
            raise ValueError("A simulation is already active. This laptop-friendly version runs one job at a time; wait for it or cancel it first.")
        job_id = uuid.uuid4().hex[:16]
        folder = config.JOBS_DIR / job_id
        folder.mkdir()
        total = round(settings.duration_ps * 1000 / settings.timestep_fs)
        now = dt.datetime.now(dt.timezone.utc).isoformat()
        job = {"id": job_id, "name": settings.name, "engine": settings.engine, "status": "queued", "stage": "Starting worker", "progress": 0, "completed_steps": 0, "total_steps": total, "elapsed_seconds": 0, "logs": [f"Queued {total:,} production steps on {config.CPU_THREADS} CPU threads."], "config": settings.model_dump(), "created_at": now}
        atomic_json(folder / "config.json", settings.model_dump())
        shutil.copy2(dataset_dir(settings.dataset_id) / "topology.pdb", folder / "input.pdb")
        atomic_json(folder / "status.json", job)
        environment = os.environ.copy()
        environment.update(OPENMM_CPU_THREADS=str(config.CPU_THREADS), OMP_NUM_THREADS=str(config.CPU_THREADS), OPENBLAS_NUM_THREADS=str(config.CPU_THREADS), MKL_NUM_THREADS=str(config.CPU_THREADS), DYNAMOL_DATA_DIR=str(config.DATA_ROOT), PYTHONUNBUFFERED="1")
        with (folder / "worker.log").open("ab", buffering=0) as output:
            process = subprocess.Popen([sys.executable, "-m", "backend.worker", job_id], cwd=config.ROOT, env=environment, stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
        job["worker_pid"] = process.pid
        # Worker waits for PID-bearing initial status so writes cannot race.
        atomic_json(folder / "status.json", job)
        return job


def cancel_job(job_id: str) -> dict:
    with _lock:
        job = get_job(job_id)
        if job["status"] not in {"queued", "running", "cancelling"}:
            return job
        folder = config.JOBS_DIR / safe_id(job_id)
        (folder / "cancel.request").touch()
        pid = job.get("worker_pid")
        if pid:
            probe = subprocess.run(["ps", "-p", str(pid), "-o", "command="], capture_output=True, text=True)
            if "backend.worker" in probe.stdout and job_id in probe.stdout:
                try:
                    os.killpg(pid, signal.SIGTERM)
                    time.sleep(0.2)
                    still_running = subprocess.run(["ps", "-p", str(pid), "-o", "command="], capture_output=True, text=True)
                    if "backend.worker" in still_running.stdout and job_id in still_running.stdout:
                        os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        # The worker may be stopped inside a native engine; parent persists state.
        job.update(status="cancelled", stage="Cancelled", error="Cancelled by user. Partial outputs are retained.")
        job["logs"] = (job["logs"] + ["Cancellation requested; worker process group terminated."])[-150:]
        atomic_json(folder / "status.json", job)
        return job
