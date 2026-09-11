"""Framewise geometry on original coordinates, never interpolated display frames."""
import mdtraj as md
import numpy as np

from .models import MeasurementRequest
from .storage import load_physical


def measure(dataset_id: str, request: MeasurementRequest) -> dict:
    return _measure(load_physical(dataset_id), request, strict=True)


def preview(dataset_id: str, request: MeasurementRequest) -> dict:
    """Transient saved-frame readout; undefined frames do not hide valid moments."""
    return _measure(load_physical(dataset_id), request, strict=False)


def _measure(traj: md.Trajectory, request: MeasurementRequest, *, strict: bool) -> dict:
    indices = request.atoms
    expected = {"distance": 2, "angle": 3, "dihedral": 4, "hbond": 3}[request.kind]
    if len(indices) != expected or len(set(indices)) != expected:
        raise ValueError(f"Select {expected} distinct atoms for a {request.kind} measurement.")
    if any(i < 0 or i >= traj.n_atoms for i in indices):
        raise ValueError("An atom index is outside the selected topology.")
    periodic = traj.unitcell_vectors is not None
    warnings = [] if periodic else ["No periodic box found; using Cartesian geometry."]
    # Explicit displacement checks reject zero bonds and undefined torsion planes.
    pairs = [[indices[i], indices[i + 1]] for i in range(len(indices) - 1)]
    vectors = md.compute_displacements(traj, pairs, periodic=periodic, opt=True)
    lengths = np.linalg.norm(vectors, axis=-1)
    errors = np.full(traj.n_frames, None, dtype=object)
    coincident = np.any(lengths < 1e-7, axis=1)
    if strict and np.any(coincident):
        raise ValueError("Selected atoms have coincident coordinates in at least one frame; geometry is undefined.")
    errors[coincident] = "Selected atoms have coincident coordinates in this frame; geometry is undefined."
    result = {"kind": request.kind, "atoms": indices, "times_ps": traj.time.tolist(), "warnings": warnings}
    if request.kind == "distance":
        values = md.compute_distances(traj, [indices], periodic=periodic)[:, 0] * 10
        unit = "Å"
    elif request.kind == "angle":
        values = np.rad2deg(md.compute_angles(traj, [indices], periodic=periodic)[:, 0])
        unit = "°"
    elif request.kind == "dihedral":
        collinear = (np.linalg.norm(np.cross(vectors[:, 0], vectors[:, 1]), axis=-1) < 1e-9) | (np.linalg.norm(np.cross(vectors[:, 1], vectors[:, 2]), axis=-1) < 1e-9)
        if strict and np.any(collinear):
            raise ValueError("A selected dihedral contains collinear bonds in at least one frame; its torsion is undefined.")
        errors[collinear & ~coincident] = "This frame contains collinear bonds; the dihedral is undefined."
        values = np.rad2deg(md.compute_dihedrals(traj, [indices], periodic=periodic)[:, 0])
        unit = "°"
        warnings.append("Torsions are wrapped to −180°…180°; a jump across that boundary is a coordinate convention, not a sudden molecular motion.")
    else:
        donor, hydrogen, acceptor = [traj.topology.atom(i) for i in indices]
        if hydrogen.element is None or hydrogen.element.symbol != "H":
            raise ValueError("Hydrogen-bond selection order is donor, explicit hydrogen, acceptor (D,H,A). The second atom must be hydrogen.")
        if any(atom.element is None or atom.element.symbol not in {"N", "O", "S"} for atom in (donor, acceptor)):
            raise ValueError("Select nitrogen, oxygen, or sulfur as donor and acceptor; this is a geometric screen, not an acceptor chemistry assignment.")
        bonds = {frozenset((a.index, b.index)) for a, b in traj.topology.bonds}
        if frozenset(indices[:2]) not in bonds:
            raise ValueError("The selected hydrogen is not bonded to the donor in the topology.")
        if frozenset((indices[1], indices[2])) in bonds:
            raise ValueError("The selected acceptor is covalently bonded to the hydrogen.")
        values = md.compute_distances(traj, [[indices[0], indices[2]]], periodic=periodic)[:, 0] * 10
        angles = np.rad2deg(md.compute_angles(traj, [indices], periodic=periodic)[:, 0])
        occupied = (values <= 3.5 + 1e-6) & (angles >= 150.0 - 1e-5)
        if strict:
            result.update(occupancy=float(occupied.mean()), angle_values=angles.tolist())
        else:
            errors[~np.isfinite(angles)] = "This selection produces an undefined angle in this frame."
            result.update(angle_values=angles.tolist(), geometry_passes=occupied.tolist())
        warnings.append("Geometric H-bond occupancy: D–A ≤ 3.5 Å and D–H–A ≥ 150° at the hydrogen. Element/connectivity checks do not establish chemical donor/acceptor eligibility; protonation and electronic state need review.")
        unit = "Å"
    if strict:
        if not np.isfinite(values).all():
            raise ValueError("This selection produces undefined geometry in one or more frames.")
        return {**result, "values": values.tolist(), "unit": unit}
    errors[~np.isfinite(values)] = "This selection produces undefined geometry in this frame."
    valid = errors == None  # noqa: E711 — elementwise comparison for the object array
    # JSON must contain null, never NaN/Infinity or a plausible number for an undefined frame.
    if "angle_values" in result:
        result["angle_values"] = [angle if ok else None for angle, ok in zip(result["angle_values"], valid)]
        result["geometry_passes"] = [occupied if ok else None for occupied, ok in zip(result["geometry_passes"], valid)]
    return {**result, "values": [float(value) if ok else None for value, ok in zip(values, valid)], "unit": unit, "frame_errors": errors.tolist()}
