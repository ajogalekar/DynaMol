"""DynaMol's localhost-only API. Run: uvicorn backend.main:app --host 127.0.0.1 --port 8765."""
import math
import hashlib
import shutil
import tempfile
import zipfile
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile, Request
from fastapi.exceptions import RequestValidationError
from starlette.concurrency import run_in_threadpool
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from starlette.middleware.base import BaseHTTPMiddleware

from . import config, jobs, storage
from .analysis import measure
from .models import MeasurementRequest, SimulationConfig, StructureFetchRequest, SmilesRequest, PreparationRequest, SolvationRequest, InspectionRequest


class LocalOriginMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        origin = request.headers.get("origin")
        allowed = {"http://localhost:5173", "http://127.0.0.1:5173", "http://localhost:8765", "http://127.0.0.1:8765", "http://localhost:4173", "http://127.0.0.1:4173"}
        if origin and origin not in allowed:
            return JSONResponse({"detail": "This local API accepts only DynaMol's local browser origin."}, status_code=403)
        if request.method in {"POST", "PUT", "DELETE", "PATCH"}:
            length = request.headers.get("content-length")
            if length and int(length) > config.MAX_UPLOAD_BYTES * 2 + 1024 * 1024:
                return JSONResponse({"detail": "Upload exceeds the 250 MiB per file limit."}, status_code=413)
        return await call_next(request)


app = FastAPI(title="DynaMol", version="0.1.0", description="Local molecular dynamics studio — real trajectories, explicit provenance.")
app.add_middleware(CORSMiddleware, allow_origins=["http://localhost:5173", "http://127.0.0.1:5173", "http://localhost:4173", "http://127.0.0.1:4173"], allow_methods=["GET", "POST"], allow_headers=["Content-Type"])
app.add_middleware(LocalOriginMiddleware)


@app.exception_handler(RequestValidationError)
async def validation_error(request, exc):
    return JSONResponse(status_code=422, content={"detail": "; ".join(str(error["msg"]) for error in exc.errors())})


@app.exception_handler(ValueError)
async def value_error(request, exc):
    return JSONResponse(status_code=422, content={"detail": str(exc)})


@app.exception_handler(FileNotFoundError)
async def not_found(request, exc):
    return JSONResponse(status_code=404, content={"detail": str(exc)})


@app.get("/api/health")
def health():
    return jobs.health()


@app.get("/api/datasets")
def datasets():
    return storage.list_datasets()


@app.get("/api/datasets/demo")
def demo():
    try:
        return storage.get_dataset("demo")
    except FileNotFoundError:
        raise HTTPException(404, "The demonstration has not been generated. Run .venv/bin/python scripts/generate_demo.py.")


@app.get("/api/datasets/{dataset_id}")
def dataset(dataset_id: str):
    return storage.get_dataset(dataset_id)


@app.get("/api/datasets/{dataset_id}/topology")
def topology(dataset_id: str):
    storage.get_dataset(dataset_id)
    return FileResponse(storage.dataset_dir(dataset_id) / "topology.pdb", media_type="chemical/x-pdb")


@app.get("/api/datasets/{dataset_id}/coordinates")
def coordinates(dataset_id: str):
    storage.get_dataset(dataset_id)
    return FileResponse(storage.dataset_dir(dataset_id) / "coordinates.bin", media_type="application/octet-stream", headers={"Cache-Control": "private, max-age=3600"})


@app.get("/api/datasets/{dataset_id}/prepared")
def prepared_structure(dataset_id: str):
    storage.get_dataset(dataset_id)
    path = storage.dataset_dir(dataset_id) / "prepared.pdb"
    if not path.is_file():
        raise FileNotFoundError("This dataset has no prepared PDB. Prepare the protein first.")
    return FileResponse(path, media_type="chemical/x-pdb", filename="DynaMol-prepared.pdb")


async def save_upload(upload: UploadFile, folder: Path, name: str) -> Path:
    suffix = Path(upload.filename or "").suffix.lower()
    destination = folder / (name + suffix)
    size = 0
    with destination.open("wb") as handle:
        while chunk := await upload.read(1024 * 1024):
            size += len(chunk)
            if size > config.MAX_UPLOAD_BYTES:
                raise HTTPException(413, "Each uploaded file must be smaller than 250 MiB. Trim it first or load a smaller trajectory.")
            handle.write(chunk)
    if size == 0:
        raise HTTPException(422, "The uploaded file is empty.")
    return destination


@app.post("/api/datasets/upload")
async def upload(topology: UploadFile = File(...), trajectory: UploadFile | None = File(None), stride: int = Form(1), frame_interval_ps: float | None = Form(None)):
    if stride < 1 or stride > 1_000_000:
        raise HTTPException(422, "Stride must be between 1 and 1,000,000.")
    if frame_interval_ps is not None and (not math.isfinite(frame_interval_ps) or frame_interval_ps <= 0):
        raise HTTPException(422, "Frame interval must be a positive number of picoseconds.")
    # UploadFile disk spooling plus chunked trajectory decoding bounds retained data.
    with tempfile.TemporaryDirectory(prefix="upload-", dir=config.DATA_ROOT) as temporary:
        folder = Path(temporary)
        top = await save_upload(topology, folder, "topology")
        traj_path = await save_upload(trajectory, folder, "trajectory") if trajectory else None
        traj, warnings = await run_in_threadpool(storage.load_uploaded, top, traj_path, stride, frame_interval_ps)
        provenance = {"source_names": {"topology": topology.filename, "trajectory": trajectory.filename if trajectory else None}, "sha256": {"topology": hashlib.sha256(top.read_bytes()).hexdigest(), "trajectory": hashlib.sha256(traj_path.read_bytes()).hexdigest() if traj_path else None}, "stride": stride, "frame_interval_ps": frame_interval_ps, "original_files_retained": True}
        result = await run_in_threadpool(storage.save_dataset, traj, Path(topology.filename or "Uploaded structure").stem, "Local upload", "User-supplied molecular coordinates. Atom order must match between topology and trajectory.", warnings=warnings, provenance=provenance)
        originals = storage.dataset_dir(result["id"]) / "originals"
        originals.mkdir()
        shutil.copy2(top, originals / top.name)
        if traj_path:
            shutil.copy2(traj_path, originals / traj_path.name)
        return result


@app.post("/api/datasets/{dataset_id}/measurements")
def measurement(dataset_id: str, request: MeasurementRequest):
    return measure(dataset_id, request)


@app.get("/api/jobs")
def list_jobs():
    return jobs.list_jobs()


@app.post("/api/jobs", status_code=201)
def submit_job(settings: SimulationConfig):
    return jobs.submit_job(settings)


@app.get("/api/jobs/{job_id}")
def job(job_id: str):
    return jobs.get_job(job_id)


@app.post("/api/jobs/{job_id}/cancel")
def cancel(job_id: str):
    return jobs.cancel_job(job_id)


@app.get("/api/jobs/{job_id}/download")
def download(job_id: str):
    job = jobs.get_job(job_id)
    if job["status"] in {"queued", "running", "cancelling"}:
        raise HTTPException(409, "Wait for the job to finish, or cancel it, before downloading a consistent output archive.")
    folder = config.JOBS_DIR / storage.safe_id(job_id)
    target = folder / "dynamol-output.zip"
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(folder.rglob("*")):
            if path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(folder.resolve()) and path != target and path.suffix != ".tmp":
                archive.write(path, arcname=str(path.relative_to(folder)))
    return FileResponse(target, media_type="application/zip", filename=f"DynaMol-{job_id}.zip")


@app.post("/api/structures/upload")
async def import_structure(file: UploadFile = File(...)):
    from . import sources
    with tempfile.TemporaryDirectory(prefix="source-", dir=config.DATA_ROOT) as temporary:
        path = await save_upload(file, Path(temporary), "structure")
        return await run_in_threadpool(sources.import_structure, path, Path(file.filename or "Molecule").stem, {"source_name": file.filename})


@app.post("/api/structures/fetch")
async def fetch_structure(settings: StructureFetchRequest):
    from . import sources
    return await run_in_threadpool(sources.fetch_structure, **settings.model_dump())


@app.post("/api/structures/smiles")
async def smiles_structure(settings: SmilesRequest):
    from . import sources
    return await run_in_threadpool(sources.create_smiles, **settings.model_dump())


@app.get("/api/datasets/{dataset_id}/inspection")
async def inspect_protein(dataset_id: str, ph: float = 7.0):
    from . import preparation
    return await run_in_threadpool(preparation.inspect_preparation, dataset_id, ph)


@app.post("/api/preparations", status_code=201)
def prepare_protein(settings: PreparationRequest):
    from . import preparation
    return preparation.submit_preparation(settings.model_dump())


@app.post("/api/datasets/{dataset_id}/inspection")
async def inspect_complex(dataset_id: str, settings: InspectionRequest):
    from . import preparation
    return await run_in_threadpool(preparation.inspect_preparation, dataset_id, settings.ph, settings.ligand_overrides)


@app.post("/api/datasets/{dataset_id}/solvate", status_code=201)
def solvate_structure(dataset_id: str, settings: SolvationRequest):
    from . import preparation
    return preparation.submit_solvation(dataset_id, settings.model_dump())


# API routes above remain authoritative. A production build can run as one
# localhost application, without a separate Node/Vite process.
frontend_build = config.ROOT / "frontend" / "dist"
if frontend_build.is_dir():
    app.mount("/", StaticFiles(directory=frontend_build, html=True), name="frontend")
