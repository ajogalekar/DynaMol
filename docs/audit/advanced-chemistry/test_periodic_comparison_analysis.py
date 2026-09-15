"""Adversarial numerical witnesses; synthetic policies are not release policies."""
from copy import deepcopy

import numpy as np
import pytest

from periodic_comparison_analysis import UNITS, analyze_pair


POLICY = {"profile_id": "synthetic-test-only", "basis": "Closed-form two-particle harmonic witness; not a native runtime calibration",
          "limits": {"force_rms_kj_mol_nm": 0.01, "force_max_kj_mol_nm": 0.1, "energy_absolute_kj_mol": 0.001}}


def harmonic_snapshot(n=2):
    # k=1000 kJ mol^-1 nm^-2, r0=.15 nm, r=.16 nm:
    # E=.05 kJ/mol; forces are +10 and -10 kJ/mol/nm along x.
    xyz = np.zeros((n, 3)); xyz[1, 0] = 0.16
    force = np.zeros_like(xyz); force[:2, 0] = (10, -10)
    return {"atom_ids": [f"fixture:{i}" for i in range(n)], "coords_nm": xyz,
            "forces_kj_mol_nm": force, "box_nm": np.eye(3) * 3,
            "potential_energy_kj_mol": 0.05, "units": dict(UNITS),
            "metadata": {"runtime": {"engine": "analytic-harmonic-fixture", "precision": "double"},
                         "force_conventions": {"model": "unconstrained harmonic bond", "k": 1000.0, "r0_nm": 0.15},
                         "energy_conventions": {"zero": "minimum at zero"},
                         "evidence": {"scope": "Synthetic calculation; no protein or native engine claim"}}}


def test_exact_witness_and_no_implicit_release_approval():
    a = harmonic_snapshot(); r = analyze_pair(a, deepcopy(a), policy=POLICY, focus_groups={"bond": a["atom_ids"]})
    assert r["status"] == "analyzed"
    assert r["force_assessment"] == r["energy_assessment"] == "pass"
    assert r["force_metrics"]["reference_rms_vector_force_kj_mol_nm"] == 10
    assert not r["release_approved"] and not r["parameter_mapping_validated_by_this_analyzer"]


def test_one_bad_site_cannot_hide_in_bulk_rms():
    a = harmonic_snapshot(1000); b = deepcopy(a); b["forces_kj_mol_nm"][0, 1] = 0.2
    r = analyze_pair(a, b, policy=POLICY, focus_groups={"site": ["fixture:0", "fixture:1"]})
    assert r["force_metrics"]["rms_vector_difference_kj_mol_nm"] < POLICY["limits"]["force_rms_kj_mol_nm"]
    assert r["force_assessment"] == "fail"
    assert r["focus_group_metrics"]["site"]["maximum_vector_difference_kj_mol_nm"] == 0.2


def test_constant_energy_error_is_preserved_and_fails_raw_gate():
    a = harmonic_snapshot(); b = deepcopy(a); b["potential_energy_kj_mol"] += 10
    r = analyze_pair(a, b, policy=POLICY)
    assert r["force_assessment"] == "pass" and r["energy_assessment"] == "fail"
    assert r["raw_energies_kj_mol"]["candidate_minus_reference"] == 10
    assert r["energy_offset_applied_kj_mol"] == 0


def test_tail_convention_difference_is_not_an_energy_pass():
    a = harmonic_snapshot(); b = deepcopy(a)
    b["metadata"]["energy_conventions"]["tail"] = "explicit volume-dependent native convention"
    b["potential_energy_kj_mol"] -= 0.457
    r = analyze_pair(a, b, policy=POLICY)
    assert r["force_assessment"] == "pass"
    assert r["energy_assessment"] == "not_comparable"
    assert r["raw_energies_kj_mol"]["candidate_minus_reference"] == pytest.approx(-0.457)


def test_force_convention_difference_cannot_pass_at_degenerate_pose():
    a = harmonic_snapshot(); b = deepcopy(a)
    b["metadata"]["force_conventions"]["undisclosed_exception_change"] = True
    r = analyze_pair(a, b, policy=POLICY)
    assert r["force_metrics"]["maximum_vector_difference_kj_mol_nm"] == 0
    assert r["force_assessment"] == r["energy_assessment"] == "not_comparable"


@pytest.mark.parametrize("change", ["position", "box", "order", "units", "nan", "duplicate"])
def test_input_corruption_rejected(change):
    a = harmonic_snapshot(); b = deepcopy(a)
    if change == "position": b["coords_nm"][1, 0] += 1e-12
    elif change == "box": b["box_nm"][0, 0] += 1e-12
    elif change == "order": b["atom_ids"].reverse()
    elif change == "units": b["units"]["coordinates"] = "angstrom"
    elif change == "nan": b["forces_kj_mol_nm"][0, 0] = np.nan
    elif change == "duplicate": b["atom_ids"][1] = b["atom_ids"][0]
    r = analyze_pair(a, b, policy=POLICY)
    assert r["status"] == "invalid_input" and r["failures"]
    assert r["force_assessment"] == r["energy_assessment"] == "invalid_input"


def test_no_profile_means_no_pass_threshold_is_invented():
    a = harmonic_snapshot(); r = analyze_pair(a, deepcopy(a))
    assert r["status"] == "analyzed"
    assert r["force_assessment"] == r["energy_assessment"] == "not_assessed"


@pytest.mark.parametrize("change", ["missing_max", "fitted_offset", "negative_limit"])
def test_invalid_or_offset_policy_rejected(change):
    a = harmonic_snapshot(); policy = deepcopy(POLICY)
    if change == "missing_max": policy["limits"].pop("force_max_kj_mol_nm")
    elif change == "fitted_offset": policy["energy_offset_kj_mol"] = 10
    else: policy["limits"]["energy_absolute_kj_mol"] = -1
    r = analyze_pair(a, deepcopy(a), policy=policy)
    assert r["status"] == "invalid_input" and not r["release_approved"]


def test_unknown_site_id_is_not_silently_dropped():
    a = harmonic_snapshot()
    r = analyze_pair(a, deepcopy(a), policy=POLICY, focus_groups={"zinc": ["missing_atom"]})
    assert r["status"] == "invalid_input"


@pytest.mark.parametrize("field", ["groups", "limits"])
def test_malformed_mapping_returns_invalid_instead_of_crashing(field):
    a = harmonic_snapshot(); policy = deepcopy(POLICY)
    if field == "groups":
        r = analyze_pair(a, deepcopy(a), policy=policy, focus_groups=["fixture:0"])
    else:
        policy["limits"] = list(policy["limits"])
        r = analyze_pair(a, deepcopy(a), policy=policy)
    assert r["status"] == "invalid_input"


def test_metadata_report_is_not_rewritten_by_later_input_changes():
    a = harmonic_snapshot(); b = deepcopy(a); r = analyze_pair(a, b, policy=POLICY)
    b["metadata"]["energy_conventions"]["zero"] = "changed after audit"
    assert r["metadata"]["candidate"]["energy_conventions"]["zero"] == "minimum at zero"
