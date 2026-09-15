"""Read-only numeric lane of the periodic validation contract.

This does not inspect topology, validate chemistry, certify conventions supplied
by a caller, or authorize a simulation. It never fits/subtracts an energy offset
and has no default pass thresholds. Production policy must be independently
frozen and supported by the evidence described in PERIODIC-VALIDATION.md.
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any

import numpy as np


UNITS = {"coordinates": "nm", "box": "nm", "forces": "kJ/mol/nm", "energy": "kJ/mol"}


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def _finite_number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return float(value)


def _snapshot(value):
    if not isinstance(value, dict):
        raise ValueError("Each snapshot must be a structured record")
    if value.get("units") != UNITS:
        raise ValueError("Explicit coordinate, box, force and energy units must match the contract")
    ids = value.get("atom_ids")
    if (not isinstance(ids, list) or not ids or any(not isinstance(x, str) or not x for x in ids)
            or len(ids) != len(set(ids))):
        raise ValueError("Atom IDs must be nonempty, unique strings in actual saved-array order")
    arrays = {}
    for key in ("coords_nm", "forces_kj_mol_nm", "box_nm"):
        a = np.asarray(value[key], dtype=np.float64)
        shape = (3, 3) if key == "box_nm" else (len(ids), 3)
        if a.shape != shape or not np.isfinite(a).all():
            raise ValueError(f"{key} must be a finite array of shape {shape}")
        arrays[key] = a
    volume = float(np.linalg.det(arrays["box_nm"]))
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError("The declared periodic box must have finite positive volume")
    energy = _finite_number(value.get("potential_energy_kj_mol"), "potential energy")
    metadata = value.get("metadata")
    if not isinstance(metadata, dict):
        raise ValueError("Runtime and convention metadata are required")
    for key in ("runtime", "force_conventions", "energy_conventions", "evidence"):
        if not isinstance(metadata.get(key), dict) or not metadata[key]:
            raise ValueError(f"Nonempty {key} metadata is required")
    metadata_hash = canonical_hash(metadata)  # Also rejects NaN metadata.
    return ids, arrays, energy, json.loads(json.dumps(metadata, allow_nan=False)), metadata_hash


def _force_metrics(reference, candidate):
    with np.errstate(over="ignore", invalid="ignore"):
        delta = np.linalg.norm(candidate - reference, axis=1)
        ref = np.linalg.norm(reference, axis=1)
        rms = float(np.sqrt(np.mean(delta * delta)))
        ref_rms = float(np.sqrt(np.mean(ref * ref)))
    if not np.isfinite(delta).all() or not math.isfinite(rms) or not math.isfinite(ref_rms):
        raise ValueError("Force metric overflow; no finite-agreement claim can be made")
    return {"rms_vector_difference_kj_mol_nm": rms,
            "maximum_vector_difference_kj_mol_nm": float(delta.max()),
            "p99_vector_difference_kj_mol_nm": float(np.percentile(delta, 99)),
            "reference_rms_vector_force_kj_mol_nm": ref_rms,
            "relative_rms_vector_difference": rms / ref_rms if ref_rms else None}


def analyze_pair(reference, candidate, *, policy=None, focus_groups=None):
    """Compare one physical configuration; preserve raw metrics and separate lanes.

    `force_conventions` must include every force-affecting setting in the parent
    contract. `energy_conventions` additionally declares analytical tail/self/
    shift conventions. These are upstream assertions, not facts inferred here.
    Runtime precision/compiler/thread metadata may differ and remain recorded.
    A supplied policy requires both absolute RMS and maximum force limits; an
    RMS-only test can hide a local site error in a large solvent box.
    """
    result = {"schema_version": 1, "status": "invalid_input", "release_approved": False,
              "parameter_mapping_validated_by_this_analyzer": False,
              "scientific_model_validated_by_this_analyzer": False,
              "energy_offset_applied_kj_mol": 0.0, "failures": []}
    try:
        ids, ref, re, rm, rh = _snapshot(reference)
        other_ids, obs, oe, om, oh = _snapshot(candidate)
        result["metadata"] = {"reference": rm, "candidate": om}
        result["metadata_sha256"] = {"reference": rh, "candidate": oh}
        result["atom_ids_sha256"] = canonical_hash(ids)
        if ids != other_ids:
            raise ValueError("Atom identity/order mismatch; this analyzer never silently reorders arrays")
        for key in ("coords_nm", "box_nm"):
            if not np.array_equal(ref[key], obs[key]):
                raise ValueError(f"{key} differs; recompute both engines at the same saved geometry/box")
        result["identical_saved_geometry_and_box"] = True
        result["array_sha256"] = {
            label: {key: hashlib.sha256(a.astype("<f8").tobytes()).hexdigest()
                    for key, a in arrays.items()} for label, arrays in (("reference", ref), ("candidate", obs))}
        result["force_metrics"] = _force_metrics(ref["forces_kj_mol_nm"], obs["forces_kj_mol_nm"])
        delta_e = oe - re
        if not math.isfinite(delta_e):
            raise ValueError("Raw energy difference overflow")
        result["raw_energies_kj_mol"] = {"reference": re, "candidate": oe, "candidate_minus_reference": delta_e}
        result["focus_group_metrics"] = {}
        lookup = {name: i for i, name in enumerate(ids)}
        if focus_groups is not None and not isinstance(focus_groups, dict):
            raise ValueError("Focus groups must be a named mapping of atom-ID lists")
        for name, members in (focus_groups or {}).items():
            if (not isinstance(name, str) or not name or not isinstance(members, list) or not members or
                    len(members) != len(set(members)) or any(x not in lookup for x in members)):
                raise ValueError("Focus groups must contain explicit unique known atom IDs")
            indices = [lookup[x] for x in members]
            result["focus_group_metrics"][name] = _force_metrics(
                ref["forces_kj_mol_nm"][indices], obs["forces_kj_mol_nm"][indices])
        force_comparable = canonical_hash(rm["force_conventions"]) == canonical_hash(om["force_conventions"])
        energy_comparable = force_comparable and canonical_hash(rm["energy_conventions"]) == canonical_hash(om["energy_conventions"])
        result["comparability"] = {"atomic_forces": force_comparable, "raw_absolute_energy": energy_comparable,
                                   "basis": "Caller-supplied complete convention metadata; upstream source verification is mandatory."}
        result["force_assessment"] = "not_assessed" if force_comparable else "not_comparable"
        result["energy_assessment"] = "not_assessed" if energy_comparable else "not_comparable"
        if policy is not None:
            if not isinstance(policy, dict) or set(policy) != {"profile_id", "basis", "limits"}:
                raise ValueError("Only named policy, calibration basis and explicit limits are accepted; no energy offsets")
            if (not isinstance(policy.get("profile_id"), str) or not policy["profile_id"] or
                    not isinstance(policy.get("basis"), str) or not policy["basis"]):
                raise ValueError("A named policy and independent calibration basis are required")
            limits = policy["limits"]
            expected = {"force_rms_kj_mol_nm", "force_max_kj_mol_nm", "energy_absolute_kj_mol"}
            if not isinstance(limits, dict) or set(limits) != expected:
                raise ValueError("A policy must explicitly declare RMS force, maximum force and raw absolute-energy limits")
            for name, value in limits.items():
                if _finite_number(value, name) < 0:
                    raise ValueError("Policy limits cannot be negative")
            result["policy"] = json.loads(json.dumps(policy, allow_nan=False))
            result["policy_sha256"] = canonical_hash(policy)
            if force_comparable:
                m = result["force_metrics"]
                result["force_assessment"] = "pass" if (
                    m["rms_vector_difference_kj_mol_nm"] <= limits["force_rms_kj_mol_nm"] and
                    m["maximum_vector_difference_kj_mol_nm"] <= limits["force_max_kj_mol_nm"]) else "fail"
            if energy_comparable:
                result["energy_assessment"] = "pass" if abs(delta_e) <= limits["energy_absolute_kj_mol"] else "fail"
        result["status"] = "analyzed"
    except (ValueError, TypeError, KeyError, OverflowError) as exc:
        result["failures"].append(str(exc))
        result["status"] = "invalid_input"
        # An invalid later field must never leave an earlier apparent pass.
        result["force_assessment"] = result["energy_assessment"] = "invalid_input"
    return result
