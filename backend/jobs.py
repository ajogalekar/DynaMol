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
from .prepared_system import copy_ligand_parameters, ligand_parameter_files, load_prepared_forcefield
from .storage import atomic_json, dataset_dir, get_dataset, safe_id

_lock = threading.RLock()
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
    with _lock:
        return _get_job(job_id)


def _get_job(job_id: str) -> dict:
    folder = config.JOBS_DIR / safe_id(job_id)
    path = folder / "status.json"
    if not path.exists():
        raise FileNotFoundError(f"Job '{job_id}' was not found.")
    original_status = path.read_text()
    job = json.loads(original_status)
    # Terminal status alone is never permission to start a second native engine.
    if job["status"] in {"failed", "cancelled", "interrupted"} and job.get("native_pid"):
        native = subprocess.run(["ps", "-p", str(job["native_pid"]), "-o", "command="], capture_output=True, text=True, timeout=2)
        if job_id in native.stdout and "mdrun" in native.stdout:
            if not job.get("orphan_stop_requested"):
                try:
                    os.kill(job["native_pid"], signal.SIGTERM)
                except ProcessLookupError:
                    pass
            job.update(status="cancelling", stage="Waiting for native engine to stop", orphan_stop_requested=True)
            if path.read_text() != original_status:
                return _get_job(job_id)
            atomic_json(path, job)
            original_status = path.read_text()
    if job["status"] in {"queued", "running", "cancelling"}:
        restart = job.get("restarts", [])[-1] if job.get("restarts") else None
        active_since = restart["requested_at"] if restart else job["created_at"]
        previous_elapsed = restart.get("previous_elapsed_seconds", 0) if restart else 0
        active_elapsed = round(previous_elapsed + (dt.datetime.now(dt.timezone.utc) - dt.datetime.fromisoformat(active_since)).total_seconds(), 2)
        pid = job.get("worker_pid")
        if pid:
            try:
                os.kill(pid, 0)
                process = subprocess.run(["ps", "-p", str(pid), "-o", "command="], capture_output=True, text=True, timeout=2)
                if not (job_id in process.stdout and any(module in process.stdout for module in ("backend.worker", "backend.preparation_worker"))):
                    raise ProcessLookupError()
            except ProcessLookupError:
                native_pid = job.get("native_pid")
                native = subprocess.run(["ps", "-p", str(native_pid), "-o", "command="], capture_output=True, text=True, timeout=2) if native_pid else None
                if native and job_id in native.stdout and "mdrun" in native.stdout:
                    # A killed supervisor can leave a native child alive. Stop
                    # that exact job's child before ever permitting continuation.
                    if not job.get("orphan_stop_requested"):
                        try:
                            os.kill(native_pid, signal.SIGTERM)
                        except ProcessLookupError:
                            pass
                    job.update(status="cancelling", stage="Saving orphaned engine progress", orphan_stop_requested=True, error="The worker stopped. Waiting for its native engine to save and exit before recovery.")
                else:
                    job.pop("native_pid", None)
                    # Stop can arrive during interpreter startup, before the
                    # worker installs its signal handler or writes any output.
                    # Preserve that explicit intent after both processes exit;
                    # orphan cleanup alone is not a user cancellation request.
                    if (folder / "cancel.request").exists():
                        job.update(status="cancelled", stage="Cancelled", error="Cancelled by user. Partial outputs and logs retained.")
                    else:
                        job.update(status="interrupted", stage="Interrupted", error="The worker stopped before completing. Outputs and logs have been retained.")
                    job["elapsed_seconds"] = active_elapsed
                # A worker may have written its final status just before the
                # liveness probe. Never overwrite that newer disk generation.
                if path.read_text() != original_status:
                    return _get_job(job_id)
                atomic_json(path, job)
        if job["status"] in {"queued", "running", "cancelling"}:
            job["elapsed_seconds"] = active_elapsed
    if job.get("engine") in {"openmm", "gromacs"}:
        from .recovery import recovery_info
        job["recovery"] = recovery_info(job)
    return job


def list_jobs() -> list[dict]:
    result = []
    for path in sorted(config.JOBS_DIR.glob("*/status.json"), key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            result.append(get_job(path.parent.name))
        except (OSError, ValueError):
            pass
    return result


def validate_simulation(settings: SimulationConfig) -> dict:
    metadata = get_dataset(settings.dataset_id)
    if any(atom["category"] == "nucleic" for atom in metadata["atoms"]):
        raise ValueError("This simulation preset supports standard proteins only. RNA/DNA can be viewed and analyzed, but nucleic-acid and protein–nucleic-acid simulations require a separately parameterized workflow. No atoms were removed.")
    preparation_state = metadata.get("preparation")
    from .modified_residues import SUPPORTED_MODIFIED, register_topology_definitions
    modified = (preparation_state or {}).get("modified_residues", [])
    solvation_state = metadata.get("solvation")
    has_ligands = any(atom["category"] == "ligands" for atom in metadata["atoms"])
    ligand_files = ligand_parameter_files(dataset_dir(settings.dataset_id), preparation_state, required=has_ligands)
    if (preparation_state or solvation_state) and settings.engine == "gromacs":
        raise ValueError("The current GROMACS adapter cannot yet preserve an explicitly prepared protonation state or solvent preview. Use OpenMM for this prepared dataset; GROMACS requires a separate validated state-conversion workflow.")
    if ligand_files and settings.solvent != "explicit":
        raise ValueError("Prepared protein–ligand complexes require explicit TIP3P water. GBn2 implicit parameters are not available for these ligands.")
    if (preparation_state or {}).get("requires_explicit_solvent") and settings.solvent != "explicit":
        raise ValueError("This prepared system requires explicit TIP3P water to retain its recorded modified-residue, ligand or ion parameters.")
    if preparation_state and settings.solvent == "explicit" and not solvation_state:
        raise ValueError("Create the explicit-water preview first so the simulation uses the periodic box you inspected. No hidden solvent box will be generated for a prepared protein.")
    if preparation_state and not preparation_state.get("simulation_ready", True):
        raise ValueError("This preparation is not marked simulation-ready; repair its unresolved structural issues first.")
    standard = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN"}
    if modified:
        standard |= SUPPORTED_MODIFIED
    unsupported = sorted({atom["residue"] for atom in metadata["atoms"] if atom["category"] == "protein" and atom["residue"] not in standard})
    if unsupported:
        raise ValueError("Protein residues need preparation or compatible residue parameters: " + ", ".join(unsupported) + ". Open protein preparation to review supported modifications and any missing templates; no atoms were removed.")
    if not any(atom["category"] == "protein" for atom in metadata["atoms"]):
        raise ValueError("Choose a structure containing a standard amino-acid protein.")
    if settings.solvent == "implicit" and any(atom["category"] in {"water", "ions"} for atom in metadata["atoms"]):
        raise ValueError("Implicit solvent requires a protein-only input. This structure has explicit waters or ions; choose explicit solvent or upload a prepared protein-only structure. No atoms were removed.")
    from .preparation import exact_input_path, backbone_gaps
    from openmm import app
    input_path = exact_input_path(settings.dataset_id)
    if (preparation_state or solvation_state) and input_path.name != "prepared.pdb":
        raise ValueError("The exact prepared topology is missing. Repeat preparation to preserve the requested protonation state.")
    register_topology_definitions()
    input_structure = app.PDBFile(str(input_path))
    if ligand_files or modified:
        forcefield, _ = load_prepared_forcefield(dataset_dir(settings.dataset_id), preparation_state)
        unmatched = forcefield.getUnmatchedResidues(input_structure.topology)
        if unmatched:
            raise ValueError("The saved complex parameters do not cover these residues: " + ", ".join(sorted({residue.name for residue in unmatched})) + ". Prepare the complex again; no molecules were removed.")
    if any(gap["structural_break"] for gap in backbone_gaps(input_structure.topology, input_structure.positions)):
        raise ValueError("A long backbone C–N connection indicates an unresolved structural gap. Inspect and repair the protein before simulation; an artificial stretched peptide bond will not be simulated.")
    engine = next(engine for engine in health()["engines"] if engine["id"] == settings.engine)
    if not engine["available"]:
        raise ValueError(engine["message"])
    from .resources import validate_resources
    resources = validate_resources(metadata, settings.model_dump(), input_path=input_path)
    return {"metadata": metadata, "input_path": input_path, "preparation_state": preparation_state, "solvation_state": solvation_state, "engine": engine, "resources": resources}


def submit_job(settings: SimulationConfig) -> dict:
    with _lock:
        validated = validate_simulation(settings)
        input_path = validated["input_path"]
        preparation_state = validated["preparation_state"]
        solvation_state = validated["solvation_state"]
        if any(job["status"] in {"queued", "running", "cancelling"} for job in list_jobs()):
            raise ValueError("A simulation is already active. This laptop-friendly version runs one job at a time; wait for it or cancel it first.")
        job_id = uuid.uuid4().hex[:16]
        folder = config.JOBS_DIR / job_id
        folder.mkdir()
        total = round(settings.duration_ps * 1000 / settings.timestep_fs)
        now = dt.datetime.now(dt.timezone.utc).isoformat()
        job = {"id": job_id, "name": settings.name, "engine": settings.engine, "status": "queued", "stage": "Starting worker", "progress": 0, "completed_steps": 0, "total_steps": total, "elapsed_seconds": 0, "logs": [f"Queued {total:,} production steps on {config.CPU_THREADS} CPU threads."], "config": settings.model_dump(), "created_at": now}
        atomic_json(folder / "config.json", settings.model_dump())
        shutil.copy2(input_path, folder / "input.pdb")
        copy_ligand_parameters(dataset_dir(settings.dataset_id), folder, preparation_state)
        atomic_json(folder / "input-state.json", {"preparation": preparation_state, "solvation": solvation_state})
        atomic_json(folder / "status.json", job)
        return _launch_worker(folder, job)


def _launch_worker(folder: Path, job: dict) -> dict:
    environment = os.environ.copy()
    environment.update(OPENMM_CPU_THREADS=str(config.CPU_THREADS), OMP_NUM_THREADS=str(config.CPU_THREADS), OPENBLAS_NUM_THREADS=str(config.CPU_THREADS), MKL_NUM_THREADS=str(config.CPU_THREADS), DYNAMOL_DATA_DIR=str(config.DATA_ROOT), PYTHONUNBUFFERED="1")
    # Remove the previous PID before spawning; a resumed worker waits for its own PID.
    job.pop("worker_pid", None)
    atomic_json(folder / "status.json", job)
    with (folder / "worker.log").open("ab", buffering=0) as output:
        process = subprocess.Popen([sys.executable, "-m", "backend.worker", job["id"]], cwd=config.ROOT, env=environment, stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
    job["worker_pid"] = process.pid
    atomic_json(folder / "status.json", job)
    return job


def resume_job(job_id: str) -> dict:
    from .recovery import TERMINAL, validate_manifest, digest
    with _lock:
        job = get_job(job_id)
        if job["status"] not in TERMINAL or job["engine"] not in {"openmm", "gromacs"}:
            raise ValueError("Only interrupted, cancelled or failed simulations can be resumed.")
        if any(j["status"] in {"queued", "running", "cancelling"} for j in list_jobs()):
            raise ValueError("Another job is active. Wait for it or stop it before resuming.")
        folder = config.JOBS_DIR / safe_id(job_id)
        manifest = validate_manifest(folder, job["engine"])
        checkpoint = folder / "checkpoints" / manifest["checkpoint"]["directory"] / "checkpoint.chk" if job["engine"] == "openmm" else folder / "production.cpt"
        attempt = {"requested_at": dt.datetime.now(dt.timezone.utc).isoformat(), "previous_status": job["status"], "checkpoint_sha256": digest(checkpoint), "checkpoint_step": manifest.get("checkpoint", {}).get("step"), "previous_elapsed_seconds": job.get("elapsed_seconds", 0), "original_configuration_retained": True}
        job.setdefault("restarts", []).append(attempt)
        (folder / "cancel.request").unlink(missing_ok=True)
        job.pop("orphan_stop_requested", None)
        job.update(status="queued", stage="Resuming saved checkpoint", resume_requested=True, error=None)
        job["logs"] = (job.get("logs", []) + ["Resume requested. Validated original inputs, native runtime and saved parameters; preparation will not repeat."])[-150:]
        return _launch_worker(folder, job)


def cancel_job(job_id: str) -> dict:
    with _lock:
        job = get_job(job_id)
        if job["status"] not in {"queued", "running"}:
            return job
        folder = config.JOBS_DIR / safe_id(job_id)
        (folder / "cancel.request").touch()
        job.update(status="cancelling", stage="Stopping and saving progress")
        job["logs"] = (job["logs"] + ["Stop requested. The engine is finishing its current operation and saving recoverable progress."])[-150:]
        atomic_json(folder / "status.json", job)
        pid = job.get("worker_pid")
        if pid and job.get("engine") in {"openmm", "gromacs"}:
            probe = subprocess.run(["ps", "-p", str(pid), "-o", "command="], capture_output=True, text=True)
            if any(module in probe.stdout for module in ("backend.worker", "backend.preparation_worker")) and job_id in probe.stdout:
                try:
                    # The worker owns graceful native-engine termination. Killing
                    # its process group here would discard a GROMACS checkpoint.
                    os.kill(pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        return job
