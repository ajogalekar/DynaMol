"""Readiness uses the same validation as submission; estimates never construct solvent."""
from typing import Literal

from fastapi import APIRouter
from pydantic import BaseModel, Field, ValidationError

from . import jobs, preparation, storage, resources as resources_module
from .resources import ResourceLimitError, resource_blockers
from .models import PreparationRequest, SimulationConfig

router = APIRouter(prefix="/api")


class ReadinessRequest(BaseModel):
    mode: Literal["preparation", "simulation"]
    settings: dict = Field(default_factory=dict)


def resource_estimate(metadata, settings, mode):
    return resources_module.resource_estimate(metadata, settings, mode, input_path=preparation.exact_input_path(metadata["id"]))


def readiness(dataset_id, request: ReadinessRequest):
    metadata = storage.get_dataset(dataset_id)
    settings = {**request.settings, "dataset_id": dataset_id}
    blockers, warnings = [], []
    inspection = None
    resources = None
    try:
        if request.mode == "preparation":
            parsed = PreparationRequest.model_validate(settings)
            settings, inspection = preparation.validate_preparation(parsed.model_dump())
            warnings.extend(inspection.get("warnings", []))
        else:
            parsed = SimulationConfig.model_validate(settings)
            settings = parsed.model_dump()
            validated = jobs.validate_simulation(parsed)
            resources = (validated or {}).get("resources")
    except ResourceLimitError as exc:
        resources = exc.resources
        blockers.extend(exc.blockers)
    except (ValueError, FileNotFoundError, ValidationError) as exc:
        if isinstance(exc, ValidationError):
            blockers.extend(str(error["msg"]) for error in exc.errors())
        else:
            blockers.append(str(exc))
    if any(job["status"] in {"queued", "running", "cancelling"} for job in jobs.list_jobs()):
        blockers.append("Another computation is active. Its files are retained; wait for it to finish or cancel it before starting this job.")
    try:
        if resources is None:
            resources = resource_estimate(metadata, settings, request.mode)
        blockers.extend(resource_blockers(resources, request.mode))
    except (ValueError, TypeError, OSError, ZeroDivisionError, OverflowError):
        warnings.append("Resource estimates are unavailable until the input and settings are valid.")
    prepared = metadata.get("preparation") or {}
    ligand_state = prepared.get("ligand_parameters") or {}
    chemical_state = {"ph": settings.get("ph", prepared.get("ph", 7)),
                      "modifications": prepared.get("modified_residues", []) if inspection is None else inspection.get("modified_residues", []),
                      "ligands": ligand_state.get("ligands", []) if inspection is None else inspection.get("ligands", []),
                      "ions": prepared.get("ions", []) if inspection is None else inspection.get("ions", [])}
    return {"ready": not blockers, "mode": request.mode, "dataset_id": dataset_id,
            "blockers": list(dict.fromkeys(blockers)), "warnings": list(dict.fromkeys(warnings)),
            "resources": resources, "chemical_state": chemical_state,
            "model": "Fixed-state preparation with recorded protonation assumptions" if request.mode == "preparation" else f"{settings.get('engine', 'openmm').upper()} · CPU · fixed-volume NVT · {settings.get('solvent', 'implicit')} solvent"}


@router.post("/datasets/{dataset_id}/readiness")
def check_readiness(dataset_id: str, request: ReadinessRequest):
    return readiness(dataset_id, request)
