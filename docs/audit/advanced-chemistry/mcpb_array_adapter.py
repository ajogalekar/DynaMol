"""Audit prototype: explicit atomic-unit arrays -> native MCPB/RESP data.

Run with the isolated AmberTools Python. This module neither infers chemistry nor
claims that a supplied Hessian represents an adequate physical site model.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


def _real_finite(value, label, shape=None):
    value = np.asarray(value)
    if np.iscomplexobj(value) and np.max(np.abs(value.imag), initial=0) > 1e-10:
        raise ValueError(f"{label} contains material imaginary components")
    value = np.asarray(value.real, dtype=np.float64)
    if shape is not None and value.shape != shape:
        raise ValueError(f"{label} shape {value.shape} != {shape}")
    if not np.isfinite(value).all():
        raise ValueError(f"{label} is not finite")
    return value


def native_provenance():
    from pymsmt.mcpb import gene_final_frcmod_file as native
    from pymsmt.mol import constants

    path = Path(native.__file__).resolve()
    return {
        "implementation": "native pymsmt.mcpb.gene_final_frcmod_file",
        "source_file": str(path),
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "B_TO_A": constants.B_TO_A,
        "H_TO_KCAL_MOL": constants.H_TO_KCAL_MOL,
        "HB2_TO_KCAL_MOL_A2": constants.HB2_TO_KCAL_MOL_A2,
        "numpy_version": np.__version__,
        "parameter_energy_convention": "Amber k*(x-x0)^2, not one-half*k*(x-x0)^2",
    }


def seminario_parameters(*, atom_ids, coords_bohr, hessian_hartree_per_bohr2,
                        bonds, angles, frequency_scale=1.0, average_blocks=False):
    """Extract parameters for an explicitly supplied graph (zero-based indices).

    Equilibrium values are from the supplied QM coordinates. A caller must keep
    any desired crystal-reference equilibrium policy separate and explicit.
    Names/types are not inferred: each requested term needs a unique `name` and
    ordered `atoms`, with the second atom the angle center.
    """
    from pymsmt.mcpb.gene_final_frcmod_file import (
        get_ang_fc_with_sem, get_bond_fc_with_sem,
    )

    if not atom_ids or len(set(atom_ids)) != len(atom_ids):
        raise ValueError("atom_ids must be nonempty and unique")
    n = len(atom_ids)
    coords = _real_finite(coords_bohr, "coords_bohr", (n, 3))
    hessian = _real_finite(hessian_hartree_per_bohr2, "Hessian", (3*n, 3*n))
    asymmetry = float(np.max(np.abs(hessian - hessian.T)))
    if asymmetry > 1e-8:
        raise ValueError(f"Hessian is not symmetric: max residual {asymmetry}")
    if not np.isfinite(frequency_scale) or not 0 < frequency_scale <= 2:
        raise ValueError("frequency_scale must be explicit, positive, and <= 2")
    output = {"atom_ids": list(atom_ids), "bonds": [], "angles": [],
              "method": "native MCPB classic Seminario",
              "frequency_scale": float(frequency_scale),
              "force_constant_scale": float(frequency_scale**2),
              "average_blocks": bool(average_blocks),
              "equilibrium_reference": "supplied QM coordinates",
              "hessian_max_asymmetry": asymmetry,
              "native": native_provenance(),
              "scope": "parameter extraction only; chemistry/QM/model adequacy not established"}
    seen = set()
    bond_graph = set()
    for kind, count, terms, function in (
        ("bonds", 2, bonds, get_bond_fc_with_sem),
        ("angles", 3, angles, get_ang_fc_with_sem),
    ):
        for term in terms:
            name = term["name"]
            ids = term["atoms"]
            if (not isinstance(name, str) or not name or (kind, name) in seen
                    or len(ids) != count or len(set(ids)) != count
                    or any(type(i) is not int or not 0 <= i < n for i in ids)):
                raise ValueError(f"Invalid/duplicate {kind} definition: {term}")
            seen.add((kind, name))
            if count == 2:
                bond_graph.add(tuple(sorted(ids)))
                if np.linalg.norm(coords[ids[0]]-coords[ids[1]]) < 1e-8:
                    raise ValueError(f"Coincident bonded atoms in {name}")
            else:
                # All graph bonds, including ordinary donor-neighbor bonds,
                # must be explicitly supplied. This avoids invented angles.
                for edge in ((ids[0], ids[1]), (ids[1], ids[2])):
                    if tuple(sorted(edge)) not in bond_graph:
                        raise ValueError(f"Angle {name} lacks explicit graph bond {edge}")
                u, v = coords[ids[0]]-coords[ids[1]], coords[ids[2]]-coords[ids[1]]
                if np.linalg.norm(np.cross(u, v)) < 1e-8*np.linalg.norm(u)*np.linalg.norm(v):
                    raise ValueError(f"Classic Seminario undefined for linear angle {name}")
            result = _real_finite(function([], coords.ravel().tolist(), hessian,
                                          *(i+1 for i in ids), frequency_scale,
                                          int(average_blocks)), name)
            if result[1] <= 0:
                raise ValueError(f"Nonpositive native force constant for {name}; no zero replacement")
            row = {"name": name, "atoms": list(ids),
                   "atom_ids": [atom_ids[i] for i in ids],
                   "equilibrium": float(result[0]), "force_constant": float(result[1]),
                   "equilibrium_unit": "angstrom" if count == 2 else "degree",
                   "force_constant_unit": "kcal mol^-1 angstrom^-2" if count == 2 else "kcal mol^-1 rad^-2"}
            if average_blocks:
                row["block_standard_deviation"] = float(result[2])
            output[kind].append(row)
    return output


def write_resp_esp(path, *, coords_bohr, esp_points_bohr, esp_potential_hartree_per_e):
    """Write native RESP ESP data; this is not a Gaussian output/log writer."""
    atoms = _real_finite(coords_bohr, "coords_bohr")
    points = _real_finite(esp_points_bohr, "esp_points_bohr")
    if atoms.ndim != 2 or atoms.shape[1] != 3 or not len(atoms):
        raise ValueError("coords_bohr must be nonempty N x 3")
    if points.ndim != 2 or points.shape[1] != 3 or not len(points):
        raise ValueError("esp_points_bohr must be nonempty M x 3")
    potential = _real_finite(esp_potential_hartree_per_e, "ESP", (len(points),))
    if len(atoms) > 99999 or len(points) > 99999:
        raise ValueError("Native RESP fixed-width count field would overflow")
    with Path(path).open("w") as stream:
        print("%5d%5d%5d" % (len(atoms), len(points), 0), file=stream)
        for x, y, z in atoms:
            print("%16s %15.7E %15.7E %15.7E" % (" ", x, y, z), file=stream)
        for p, (x, y, z) in zip(potential, points, strict=True):
            print("%16.7E %15.7E %15.7E %15.7E" % (p, x, y, z), file=stream)
    return {"atom_count": len(atoms), "point_count": len(points),
            "coordinate_unit": "bohr", "potential_unit": "hartree per elementary charge",
            "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}


def run_native_resp(*, executable, workdir, stage1_input, stage2_input, esp_file,
                    atom_ids, total_charge, timeout_seconds=120):
    """Fit with native RESP; native MCPB generates the supplied constraint files.

    No shell and no substitution of another charge method. This intentionally
    requires a fresh output folder and one-to-one identities for the full large
    model, including caps; transplanting standard-model charges is separate.
    """
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=False)
    if type(total_charge) is not int or not atom_ids or len(set(atom_ids)) != len(atom_ids):
        raise ValueError("RESP requires an explicit integer charge and unique atom identities")
    executable = Path(executable).resolve(strict=True)
    inputs = [Path(p).resolve(strict=True) for p in (stage1_input, stage2_input, esp_file)]
    # Native RESP has short fixed-width filename buffers. Preserve originals
    # byte-for-byte, but pass short working-directory filenames to the binary.
    for original, local in zip(inputs, ("stage1.in", "stage2.in", "input.esp"), strict=True):
        shutil.copyfile(original, workdir/local)
    results = []
    for step in (1, 2):
        command = [str(executable), "-O", "-i", f"stage{step}.in",
                   "-o", f"resp{step}.out", "-p", f"resp{step}.pch",
                   "-t", f"resp{step}.chg", "-e", "input.esp",
                   "-s", f"resp{step}-calc.esp"]
        if step == 2:
            command += ["-q", "resp1.chg"]
        completed = subprocess.run(command, cwd=workdir, capture_output=True,
                                   text=True, timeout=timeout_seconds)
        (workdir / f"resp{step}-process.log").write_text(completed.stdout+completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"Native RESP stage {step} failed: exit {completed.returncode}; see process log")
        output = (workdir / f"resp{step}.out").read_text()
        convergence = re.search(r"Convergence in\s+(\d+) iterations", output)
        if not convergence:
            raise RuntimeError(f"Native RESP stage {step} did not report convergence")
        charges = _real_finite([float(s) for s in (workdir / f"resp{step}.chg").read_text().split()],
                               "RESP charges", (len(atom_ids),))
        # Native .chg output uses six decimals. Bound its worst-case rounding;
        # this is only a total-charge serialization allowance, not fit tolerance.
        rounding_bound = 0.5e-6*len(atom_ids)+1e-10
        if abs(float(charges.sum())-total_charge) > rounding_bound:
            raise RuntimeError("RESP fitted charge sum differs from declared large-model charge")
        statistics = re.search(r"ESP relative RMS \(SQRT\(chipot/ssvpot\)\)\s+([\d.E+-]+)", output)
        results.append({"stage": step, "iterations": int(convergence.group(1)),
                        "charge_sum": float(charges.sum()), "charge_rounding_bound": rounding_bound,
                        "esp_relative_rms": float(statistics.group(1)) if statistics else None,
                        "charges_sha256": hashlib.sha256((workdir/f"resp{step}.chg").read_bytes()).hexdigest()})
    return {"method": "native AmberTools two-stage RESP", "atom_ids": atom_ids,
            "charges": charges.tolist(), "total_charge": total_charge, "stages": results,
            "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            "input_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
            "scope": "numerically converged charge fit for supplied ESP and constraints; model adequacy separate"}
