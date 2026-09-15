"""Reproduce published 1OKL MCPB terms from neutral arrays, without new QM.

AMBERHOME=$PWD/.tools/ambertools .tools/ambertools/bin/python \
  docs/audit/advanced-chemistry/verify_mcpb_array_reference.py
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np

from mcpb_array_adapter import seminario_parameters, write_resp_esp
from pymsmt.mol import constants
from pymsmt.mol.gauio import get_crds_from_fchk, get_esp_from_gau, get_matrix_from_fchk


BASE = Path(__file__).resolve().parent
REF = BASE / "reference-mcpb-1okl"


def fchk_fields(path):
    """Small independent reader, using typed FCHK field headers and counts."""
    lines = path.read_text().splitlines()
    fields = {}
    row = 2
    while row < len(lines):
        match = re.fullmatch(r"(.{43})([IRC])\s+(?:N=\s*(\d+)|(.+))", lines[row])
        row += 1
        if not match:
            continue
        key, kind, count, scalar = match.groups()
        cast = int if kind == "I" else float
        if kind == "C":
            continue
        if count is None:
            fields[key.strip()] = cast(scalar.replace("D", "E"))
            continue
        tokens = []
        while len(tokens) < int(count):
            tokens.extend(lines[row].split())
            row += 1
        if len(tokens) != int(count):
            raise ValueError(f"FCHK array count overrun: {key}")
        fields[key.strip()] = np.asarray([cast(t.replace("D", "E")) for t in tokens])
    return fields


def gaussian_esp_arrays(path):
    """Reference-only independent extraction from a genuine published log."""
    content = path.read_text()
    geometry = content.rsplit("Electrostatic Properties Using The SCF Density", 1)[1]
    atom_coords, points = [], []
    for line in geometry.splitlines():
        if "Atomic Center" in line or "ESP Fit Center" in line:
            # These are the last three whitespace-separated fields in this
            # reference, not provider-specific objects or a fabricated log.
            xyz = [float(s) / constants.B_TO_A for s in line.split()[-3:]]
            (atom_coords if "Atomic Center" in line else points).append(xyz)
    values = content.rsplit("Electrostatic Properties (Atomic Units)", 1)[1]
    potentials = [float(line.split()[-1]) for line in values.splitlines()
                  if re.match(r"\s*\d+\s+Fit\s+", line)]
    if len(points) != len(potentials):
        raise ValueError(f"ESP point/potential count differs: {len(points)}/{len(potentials)}")
    return np.asarray(atom_coords), np.asarray(points), np.asarray(potentials)


def main():
    fchk = REF / "1OKL_small_opt.fchk"
    fields = fchk_fields(fchk)
    n = fields["Number of atoms"]
    coords = fields["Current cartesian coordinates"].reshape(n, 3)
    hessian = np.zeros((3*n, 3*n))
    lower = np.tril_indices(3*n)
    hessian[lower] = fields["Cartesian Force Constants"]
    hessian[(lower[1], lower[0])] = hessian[lower]
    native_coords = np.asarray(get_crds_from_fchk(str(fchk), n)).reshape(n, 3)
    native_hessian = get_matrix_from_fchk(str(fchk), 3*n)
    assert np.array_equal(coords, native_coords)
    assert np.array_equal(hessian, native_hessian)
    atom_ids = [f"published-small-model:{i+1}" for i in range(n)]
    # Explicit indices are from the original Gaussian small-model atom order.
    # Site-donor atom types are confirmed by the published standard fingerprint;
    # neighboring atom identities are visible in the small-model input.
    site = {"Y1": 12, "Y2": 27, "Y3": 38, "M1": 46, "Y4": 63,
            "CR1": 10, "CV1": 13, "CR2": 25, "CV2": 28,
            "CC3": 37, "CR3": 39, "hn4": 64, "sy4": 60}
    bond_defs = [("M1-Y4", "M1", "Y4"), ("Y1-M1", "Y1", "M1"),
                 ("Y2-M1", "Y2", "M1"), ("Y3-M1", "Y3", "M1"),
                 ("CR1-Y1", "CR1", "Y1"), ("Y1-CV1", "Y1", "CV1"),
                 ("CR2-Y2", "CR2", "Y2"), ("Y2-CV2", "Y2", "CV2"),
                 ("CC3-Y3", "CC3", "Y3"), ("Y3-CR3", "Y3", "CR3"),
                 ("Y4-hn4", "Y4", "hn4"), ("sy4-Y4", "sy4", "Y4")]
    angle_defs = [("CC-Y3-M1", "CC3", "Y3", "M1"),
                  ("CR-Y1-M1", "CR1", "Y1", "M1"),
                  ("CR-Y2-M1", "CR2", "Y2", "M1"),
                  ("M1-Y1-CV", "M1", "Y1", "CV1"),
                  ("M1-Y2-CV", "M1", "Y2", "CV2"),
                  ("M1-Y3-CR", "M1", "Y3", "CR3"),
                  ("M1-Y4-hn", "M1", "Y4", "hn4"),
                  ("M1-Y4-sy", "M1", "Y4", "sy4"),
                  ("Y1-M1-Y2", "Y1", "M1", "Y2"),
                  ("Y1-M1-Y3", "Y1", "M1", "Y3"),
                  ("Y1-M1-Y4", "Y1", "M1", "Y4"),
                  ("Y2-M1-Y3", "Y2", "M1", "Y3"),
                  ("Y2-M1-Y4", "Y2", "M1", "Y4"),
                  ("Y3-M1-Y4", "Y3", "M1", "Y4")]
    bonds = [{"name": row[0], "atoms": [site[s]-1 for s in row[1:]]} for row in bond_defs]
    angles = [{"name": row[0], "atoms": [site[s]-1 for s in row[1:]]} for row in angle_defs]
    result = seminario_parameters(atom_ids=atom_ids, coords_bohr=coords,
                                 hessian_hartree_per_bohr2=hessian,
                                 bonds=bonds, angles=angles)
    native_rows = {}
    for line in (REF / "1OKL_mcpbpy.frcmod").read_text().splitlines():
        if "Created by Seminario method" in line:
            name, force, eq = line.split()[:3]
            native_rows[name] = (float(force), float(eq))
    comparisons = []
    for kind in ("bonds", "angles"):
        for term in result[kind]:
            if term["name"] not in native_rows:
                continue  # ordinary bond QM values are not assigned to the FF
            expected_k, expected_eq = native_rows[term["name"]]
            k_dp, eq_dp = (1, 4) if kind == "bonds" else (2, 2)
            passed = round(term["force_constant"], k_dp) == expected_k and round(term["equilibrium"], eq_dp) == expected_eq
            comparisons.append({"kind": kind, "name": term["name"], "pass": passed,
                                "observed_k": term["force_constant"], "published_k": expected_k,
                                "observed_equilibrium": term["equilibrium"], "published_equilibrium": expected_eq,
                                "comparison": "exact agreement at native published decimal precision"})
    assert len(comparisons) == len(native_rows) == 18
    assert all(row["pass"] for row in comparisons), comparisons
    atoms, points, potential = gaussian_esp_arrays(REF / "1OKL_large_mk.log")
    native_esp, adapter_esp = REF / "native-extracted.esp", REF / "array-adapter.esp"
    get_esp_from_gau(str(REF / "1OKL_large_mk.log"), str(native_esp))
    esp_report = write_resp_esp(adapter_esp, coords_bohr=atoms, esp_points_bohr=points,
                               esp_potential_hartree_per_e=potential)
    assert native_esp.read_bytes() == adapter_esp.read_bytes()
    np.savez_compressed(REF / "reference-small-qm-arrays.npz", atom_ids=atom_ids,
                        atomic_numbers=fields["Atomic numbers"], coords_bohr=coords,
                        hessian_hartree_per_bohr2=hessian,
                        gradient_hartree_per_bohr=fields["Cartesian Gradient"].reshape(n, 3),
                        energy_hartree=fields["SCF Energy"])
    np.savez_compressed(REF / "reference-large-esp-arrays.npz", coords_bohr=atoms,
                        esp_points_bohr=points, esp_potential_hartree_per_e=potential)
    # Fail-closed guards are independent perturbations of the reference input.
    guards = {}
    for label, kwargs in (
        ("wrong_hessian_shape", {"hessian_hartree_per_bohr2": hessian[:-1]}),
        ("nonfinite_hessian", {"hessian_hartree_per_bohr2": np.full_like(hessian, np.nan)}),
        ("duplicate_atom_identity", {"atom_ids": [atom_ids[0]]*n}),
        ("unmapped_angle_graph", {"bonds": bonds[:4]}),
    ):
        call = dict(atom_ids=atom_ids, coords_bohr=coords, hessian_hartree_per_bohr2=hessian,
                    bonds=bonds, angles=angles)
        call.update(kwargs)
        try:
            seminario_parameters(**call)
        except ValueError as error:
            guards[label] = {"rejected": True, "reason": str(error)}
        else:
            raise AssertionError(f"Guard failed: {label}")
    report = {
        "status": "passed", "reference": "official Amber MCPB 1OKL tutorial",
        "scope": "Native parameter/ESP adapter parity only; no new QM or new-site accuracy validation",
        "reference_fchk_sha256": hashlib.sha256(fchk.read_bytes()).hexdigest(),
        "independent_fchk_coordinates_equal_native": True,
        "independent_fchk_hessian_equal_native": True,
        "published_parameter_comparisons": comparisons,
        "esp": {**esp_report, "byte_identical_to_native_gaussian_extractor": True},
        "guard_checks": guards, "native": result["native"],
        "reference_qm": {key: fields[key] for key in (
            "Number of atoms", "Charge", "Multiplicity", "Number of electrons",
            "Number of basis functions", "Pure/Cartesian d shells", "SCF Energy", "RMS Force")},
    }
    (REF / "adapter-parity-report.json").write_text(json.dumps(report, indent=2)+"\n")
    (REF / "adapter-parameters.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "passed", "published_terms": len(comparisons),
                      "esp": report["esp"], "guards": len(guards)}, indent=2))


if __name__ == "__main__":
    main()
