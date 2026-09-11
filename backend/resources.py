"""Read-only resource planning and limits shared by readiness and job submission."""
import math
import shutil

import numpy as np
from . import config


def resource_estimate(metadata, settings, mode, *, input_path=None):
    atoms = metadata["n_atoms"]
    estimated_atoms = atoms
    estimate_kind = "saved topology"
    explicit = mode == "simulation" and settings.get("solvent") == "explicit"
    if explicit and not metadata.get("solvation"):
        from openmm import app, unit
        from .modified_residues import register_topology_definitions
        register_topology_definitions()
        if input_path is None:
            raise ValueError("An exact topology is required to estimate the solvent box.")
        pdb = app.PDBFile(str(input_path))
        xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
        padding = float(settings.get("padding_nm", 1))
        if len(xyz):
            if settings.get("engine") == "gromacs":
                side = float(np.ptp(xyz, axis=0).max()) + 2 * padding
            else:
                radius = np.linalg.norm(xyz - (xyz.min(axis=0) + xyz.max(axis=0)) / 2, axis=1).max()
                side = max(2 * radius + padding, 2 * padding)
            estimated_atoms = atoms + math.ceil(side ** 3 * 110)
            if not metadata.get("preparation"):
                estimated_atoms += sum(a["element"] not in {"H", "D"} and a["category"] == "protein" for a in metadata["atoms"])
            estimate_kind = "conservative solvent/solute estimate; actual count follows preparation and box construction"
    frames = 1
    if mode == "simulation":
        steps = round(float(settings.get("duration_ps", 10)) * 1000 / float(settings.get("timestep_fs", 2)))
        frames = math.ceil(steps / int(settings.get("report_interval", 50))) + 1
    coordinate_bytes = estimated_atoms * frames * 3 * 4
    # Multiple physical/display/native representations and checkpoints coexist.
    disk_estimate = coordinate_bytes * 6 + estimated_atoms * 2048 + 64 * 1024 ** 2
    free = shutil.disk_usage(config.DATA_ROOT).free
    return {"input_atoms": atoms, "estimated_atoms": estimated_atoms, "estimate_kind": estimate_kind,
            "saved_frames": frames, "coordinate_bytes": coordinate_bytes,
            "working_memory_estimate_bytes": coordinate_bytes * 4,
            "disk_estimate_bytes": disk_estimate, "free_disk_bytes": free,
            "atom_limit": config.MAX_ATOMS, "coordinate_limit_bytes": config.MAX_COORD_BYTES,
            "note": "Planning estimates include retained coordinates and working copies; native engines can use additional memory. Runtime is measured after dynamics starts."}



def resource_blockers(resources, mode):
    blockers = []
    if mode == "simulation" and resources["estimated_atoms"] > config.MAX_ATOMS:
        blockers.append(f"The estimated system exceeds the {config.MAX_ATOMS:,}-atom limit. Reduce solvent padding or choose a smaller prepared system.")
    if resources["coordinate_bytes"] > config.MAX_COORD_BYTES:
        blockers.append(f"The requested saved frames exceed the {config.MAX_COORD_BYTES // (1024 ** 2)} MiB coordinate limit. Increase the save interval or shorten the run.")
    if resources["disk_estimate_bytes"] + 256 * 1024 ** 2 > resources["free_disk_bytes"]:
        blockers.append("Available disk space is below the estimated output plus a 256 MiB reserve. Free space or choose fewer saved frames.")
    return blockers


class ResourceLimitError(ValueError):
    def __init__(self, resources, blockers):
        super().__init__(" ".join(blockers))
        self.resources = resources
        self.blockers = blockers


def validate_resources(metadata, settings, mode="simulation", *, input_path=None):
    try:
        resources = resource_estimate(metadata, settings, mode, input_path=input_path)
    except (OSError, TypeError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError("Resource availability could not be checked. Verify the dataset and storage location before starting this job.") from exc
    blockers = resource_blockers(resources, mode)
    if blockers:
        raise ResourceLimitError(resources, blockers)
    return resources
