"""Validated native checkpoint continuation, without repeating molecular preparation."""
from __future__ import annotations
import hashlib
import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path
from fastapi import APIRouter
from . import config
from .storage import atomic_json, safe_id

router = APIRouter()
TERMINAL = {"failed", "cancelled", "interrupted"}


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def runtime_identity(engine: str) -> dict:
    identity = {"engine": engine, "machine": platform.machine(), "os": platform.system(), "os_release": platform.release(), "cpu_threads": config.CPU_THREADS, "processor": platform.processor()}
    if platform.system() == "Darwin":
        probe = subprocess.run(["/usr/sbin/sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, timeout=2)
        identity["cpu_model"] = probe.stdout.strip()
    if engine == "openmm":
        import openmm
        root = Path(openmm.version.openmm_library_path)
        native = [path for path in root.rglob("*") if path.is_file() and (path.suffix in {".dylib", ".so", ".dll"})]
        identity.update(version=importlib.metadata.version("openmm"), platform="CPU", deterministic_forces=True, native_libraries={str(path.relative_to(root)): digest(path) for path in sorted(native)})
    else:
        from .jobs import gromacs_executable
        executable = gromacs_executable()
        if not executable:
            raise ValueError("The original GROMACS runtime is unavailable.")
        native_root = Path(executable).resolve().parent.parent / "lib"
        libraries = {path.name: digest(path) for path in sorted(native_root.glob("libgromacs*")) if path.is_file() and not path.is_symlink()}
        identity.update(executable_sha256=digest(Path(executable)), native_libraries=libraries)
    # Binary checkpoints are intentionally local to this runtime, rather than a
    # portable state interchange. Native checkpoint loading is the final guard.
    identity["worker_sha256"] = digest(Path(__file__).with_name("worker.py"))
    return identity


def dependency_hashes(folder: Path, engine: str) -> dict[str, str]:
    names = ["config.json", "input.pdb", "input-state.json"]
    names += ["prepared.pdb", "system.xml", "integrator.xml"] if engine == "openmm" else ["production.tpr", "system.gro", "topol.top"]
    files = [folder / name for name in names]
    for directory in ("ligands", "residue-parameters"):
        root = folder / directory
        if root.exists():
            files += [p for p in root.rglob("*") if p.is_file()]
    if engine == "gromacs":
        files += list(folder.glob("*.itp"))
    return {str(path.relative_to(folder)): digest(path) for path in sorted(set(files))}


def create_manifest(folder: Path, engine: str) -> dict:
    manifest = {"schema": 1, "engine": engine, "runtime": runtime_identity(engine), "dependencies": dependency_hashes(folder, engine), "checkpoint_interval_steps": 250 if engine == "openmm" else None, "checkpoint_interval_seconds": 6 if engine == "gromacs" else None}
    atomic_json(folder / "recovery.json", manifest)
    return manifest


def validate_manifest(folder: Path, engine: str) -> dict:
    path = folder / "recovery.json"
    if not path.is_file():
        raise ValueError("No production checkpoint has been saved yet. Preparation and initial relaxation must finish first.")
    manifest = json.loads(path.read_text())
    if manifest.get("schema") != 1 or manifest.get("engine") != engine:
        raise ValueError("The checkpoint format or engine is incompatible.")
    if manifest.get("runtime") != runtime_identity(engine):
        raise ValueError("The engine, CPU settings, operating system or DynaMol worker changed. Restore the original runtime to resume this binary checkpoint.")
    actual = dependency_hashes(folder, engine)
    if manifest.get("dependencies") != actual:
        raise ValueError("The saved input, configuration, topology or force-field parameters changed. Resume is blocked to preserve the original molecular system.")
    if engine == "openmm":
        checkpoint = manifest.get("checkpoint", {})
        generation = checkpoint.get("directory", "")
        if not generation.startswith("step-") or Path(generation).name != generation:
            raise ValueError("No complete OpenMM checkpoint is available.")
        target = folder / "checkpoints" / generation / "checkpoint.chk"
        if not target.is_file() or digest(target) != checkpoint.get("sha256"):
            raise ValueError("The native checkpoint is missing or damaged.")
        for frame in checkpoint.get("frames", []):
            name = frame.get("name", "")
            if Path(name).name != name or not name.endswith(".npz"):
                raise ValueError("The checkpoint frame manifest is invalid.")
            frame_path = folder / "frames" / name
            if not frame_path.is_file() or digest(frame_path) != frame.get("sha256"):
                raise ValueError("A saved trajectory frame is missing or damaged; a continuous trajectory cannot be reconstructed.")
    elif not (folder / "production.cpt").is_file():
        raise ValueError("No native GROMACS production checkpoint is available yet.")
    return manifest


def recovery_info(job: dict, *, verify: bool = False) -> dict:
    if job.get("engine") not in {"openmm", "gromacs"}:
        return {"available": False, "reason": "Structure preparation jobs can be started again from their input."}
    folder = config.JOBS_DIR / safe_id(job["id"])
    try:
        manifest = validate_manifest(folder, job["engine"]) if verify else json.loads((folder / "recovery.json").read_text())
        checkpoint = manifest.get("checkpoint", {})
        exists = bool(checkpoint.get("directory")) if job["engine"] == "openmm" else (folder / "production.cpt").is_file()
        if not exists:
            raise ValueError("The first production checkpoint has not been saved yet.")
        return {"available": job["status"] in TERMINAL, "checkpoint_saved": True, "step": checkpoint.get("step"), "reason": "Resume continues the original configuration from its last complete native checkpoint." if job["status"] in TERMINAL else "Production checkpoints are being saved automatically.", "attempts": len(job.get("restarts", []))}
    except (OSError, ValueError, KeyError) as exc:
        return {"available": False, "checkpoint_saved": False, "reason": str(exc) if not isinstance(exc, FileNotFoundError) else "No production checkpoint has been saved yet. Preparation and initial relaxation must finish first."}


@router.get("/api/jobs/{job_id}/recovery")
def inspect_recovery(job_id: str):
    from .jobs import get_job
    return recovery_info(get_job(job_id), verify=True)


@router.post("/api/jobs/{job_id}/resume")
def resume(job_id: str):
    from .jobs import resume_job
    return resume_job(job_id)
