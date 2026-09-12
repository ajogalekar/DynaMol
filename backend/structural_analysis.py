"""Explicitly defined RMSD/RMSF on saved physical coordinates and atom identities."""
from collections import defaultdict
from typing import Literal

import mdtraj as md
import numpy as np
from fastapi import APIRouter
from pydantic import BaseModel, Field

from . import config, storage

router = APIRouter(prefix="/api")


class StructuralAnalysisRequest(BaseModel):
    kind: Literal["rmsd", "rmsf"] = "rmsd"
    selection: Literal["ca", "backbone", "protein-heavy", "solute-heavy", "custom"] = "ca"
    atoms: list[int] = Field(default_factory=list, max_length=config.MAX_ATOMS)
    alignment: Literal["ca", "backbone", "selection", "none"] = "ca"
    reference_frame: int = Field(default=0, ge=0)
    start_frame: int = Field(default=0, ge=0)
    end_frame: int | None = Field(default=None, ge=0)
    stride: int = Field(default=1, ge=1, le=10000)
    periodic: Literal["whole", "cartesian"] = "whole"
    residue_average: bool = True


def selected_indices(metadata, choice, custom=()):
    if choice == "custom":
        if not custom or len(set(custom)) != len(custom) or any(type(i) is not int or i < 0 or i >= metadata["n_atoms"] for i in custom):
            raise ValueError("Choose distinct, valid atoms for the analysis selection.")
        return np.asarray(custom, dtype=int)
    result = []
    for atom in metadata["atoms"]:
        protein = atom["category"] == "protein"
        heavy = atom["element"] not in {"H", "D"}
        wanted = ((choice == "ca" and protein and atom["name"] == "CA")
                  or (choice == "backbone" and protein and atom["name"] in {"N", "CA", "C", "O"})
                  or (choice == "protein-heavy" and protein and heavy)
                  or (choice == "solute-heavy" and atom["category"] not in {"water", "ions"} and heavy))
        if wanted:
            result.append(atom["index"])
    if not result:
        raise ValueError("This structure has no atoms in that selection. Choose another group or select atoms in the viewer.")
    return np.asarray(result, dtype=int)


def _make_whole(trajectory, anchor_indices):
    try:
        molecules = trajectory.topology.find_molecules()
        if not molecules:
            raise ValueError("No connected molecules found")
        wanted = set(map(int, anchor_indices))
        anchor = max(molecules, key=lambda molecule: (len(wanted.intersection(atom.index for atom in molecule)), len(molecule)))
        trajectory.image_molecules(inplace=True, anchor_molecules=[anchor], make_whole=True)
    except Exception as exc:
        raise ValueError("Periodic structural analysis needs a usable bonded topology to make molecules whole. Supply connectivity, or explicitly choose saved Cartesian coordinates.") from exc


def structural_analysis(dataset_id, request: StructuralAnalysisRequest):
    metadata = storage.get_dataset(dataset_id)
    trajectory = storage.load_physical(dataset_id)
    if trajectory.n_atoms != metadata["n_atoms"] or trajectory.n_frames != metadata["n_frames"]:
        raise ValueError("Stored trajectory dimensions no longer match this dataset.")
    end = trajectory.n_frames - 1 if request.end_frame is None else request.end_frame
    if not (0 <= request.start_frame <= end < trajectory.n_frames) or request.reference_frame >= trajectory.n_frames:
        raise ValueError("Choose a reference and frame window inside the saved trajectory.")
    frames = np.arange(request.start_frame, end + 1, request.stride)
    if request.kind == "rmsf" and len(frames) < 2:
        raise ValueError("RMSF needs at least two saved frames in the chosen window.")
    atoms = selected_indices(metadata, request.selection, request.atoms)
    fit = atoms if request.alignment == "selection" else selected_indices(metadata, request.alignment) if request.alignment != "none" else None
    reference = trajectory[request.reference_frame]
    target = trajectory[frames]
    del trajectory
    warnings = []
    periodic = target.unitcell_vectors is not None
    if periodic and request.periodic == "whole":
        anchor = atoms if fit is None else fit
        _make_whole(target, anchor)
        _make_whole(reference, anchor)
        warnings.append("Bonded molecules are made whole and nearby periodic images are used before analysis. This is not a diffusion/unwrapping analysis.")
    elif periodic:
        warnings.append("Using saved Cartesian coordinates without periodic imaging; box crossings can create apparent structural changes.")
    if not np.isfinite(target.xyz).all() or not np.isfinite(reference.xyz).all():
        raise ValueError("Non-finite saved coordinates prevent structural analysis.")
    if fit is not None:
        if len(fit) < 3:
            raise ValueError("Alignment needs at least three non-collinear atoms. Choose another alignment group or turn alignment off.")
        reference_fit = reference.xyz[0, fit].astype(np.float64)
        if np.linalg.matrix_rank(reference_fit - reference_fit.mean(axis=0), tol=1e-6) < 2:
            raise ValueError("The alignment atoms are collinear in a saved frame; their rotational fit is undefined.")
        # Proper Kabsch rotation in float64 also handles planar three-atom fits;
        # a quaternion fit can fail to converge for this valid limiting case.
        reference_center = reference_fit.mean(axis=0)
        centered_reference = reference_fit - reference_center
        for index in range(target.n_frames):
            mobile = target.xyz[index, fit].astype(np.float64)
            mobile_center = mobile.mean(axis=0)
            if np.linalg.matrix_rank(mobile - mobile_center, tol=1e-6) < 2:
                raise ValueError("The alignment atoms are collinear in a saved frame; their rotational fit is undefined.")
            u, _, vt = np.linalg.svd((mobile - mobile_center).T @ centered_reference)
            u[:, -1] *= 1 if np.linalg.det(u @ vt) >= 0 else -1
            target.xyz[index] = (target.xyz[index].astype(np.float64) - mobile_center) @ (u @ vt) + reference_center
    else:
        warnings.append("Alignment is off: rigid translation and rotation contribute to the reported values.")
    if request.kind == "rmsd":
        reference_atoms = reference.xyz[0, atoms].astype(np.float64)
        values = np.asarray([np.sqrt(np.mean(np.sum((xyz[atoms].astype(np.float64) - reference_atoms) ** 2, axis=1))) * 10 for xyz in target.xyz])
        x = target.time.tolist()
        x_label = "Original frame index" if metadata["time_unit"] == "frame" else "Time (ps)"
        if not np.isfinite(target.time).all() or np.any(np.diff(target.time) <= 0):
            x = frames.tolist()
            x_label = "Saved frame index"
            warnings.append("Saved timestamps are not finite and strictly increasing in this window. The plot uses saved zero-based frame indices; physical elapsed time has not been inferred.")
        labels = [f"Frame {int(frame) + 1}" for frame in frames]
        groups = None
        method = "Equal-atom RMSD to the specified saved reference after the stated fit; nanometers converted to Å."
    else:
        # Welford population variance avoids full-trajectory float64 temporaries.
        # Only the original float32 trajectory and per-atom accumulators coexist.
        mean = np.zeros((len(atoms), 3), dtype=np.float64)
        squared = np.zeros_like(mean)
        for count, xyz in enumerate(target.xyz, 1):
            coordinates = xyz[atoms].astype(np.float64)
            delta = coordinates - mean
            mean += delta / count
            squared += delta * (coordinates - mean)
        variance = np.maximum(squared.sum(axis=1) / target.n_frames, 0)
        atom_values = np.sqrt(variance) * 10
        if request.residue_average:
            grouped = defaultdict(list)
            for i, atom_index in enumerate(atoms):
                grouped[target.topology.atom(int(atom_index)).residue.index].append(i)
            groups = [[int(atoms[i]) for i in indices] for indices in grouped.values()]
            values = np.asarray([np.sqrt(np.mean(variance[indices])) * 10 for indices in grouped.values()])
            labels = []
            for group in groups:
                atom = metadata["atoms"][group[0]]
                labels.append(f"{atom['chain']}:{atom['residue']}{atom['resid']}")
        else:
            values = atom_values
            groups = [[int(i)] for i in atoms]
            labels = [f"{metadata['atoms'][i]['chain']}:{metadata['atoms'][i]['residue']}{metadata['atoms'][i]['resid']}.{metadata['atoms'][i]['name']}" for i in atoms]
        x = list(range(1, len(values) + 1))
        method = "RMSF about each atom's mean position over the selected aligned frames (population mean squared displacement). Residue values are the square root of the mean atomic variance, not the mean of RMSFs."
    if not np.isfinite(values).all():
        raise ValueError("Structural analysis produced non-finite values.")
    return {"dataset_id": dataset_id, "kind": request.kind, "unit": "Å", "values": values.tolist(),
            "x": x, "labels": labels, "frame_indices": frames.tolist(), "atom_groups": groups,
            "atoms": atoms.tolist(), "alignment_atoms": [] if fit is None else fit.tolist(),
            "x_label": x_label if request.kind == "rmsd" else "Selected residue" if request.residue_average else "Selected atom",
            "request": request.model_dump(), "method": method, "warnings": warnings,
            "summary": {"min": float(values.min()), "max": float(values.max()), "mean": float(values.mean()), "frames": len(frames), "atoms": len(atoms)}}


@router.post("/datasets/{dataset_id}/structural-analysis")
def analyze_structure(dataset_id: str, request: StructuralAnalysisRequest):
    return structural_analysis(dataset_id, request)
