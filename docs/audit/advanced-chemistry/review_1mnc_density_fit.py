"""Read-only post-analysis of the saved 1MNC fixed-geometry DF comparison.

This script invokes no native chemistry executable and submits no calculation.
The output directory must be new. Source artifacts are hash checked before and
after analysis, including protection against incomplete File Provider reads.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import re
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "docs/audit/advanced-chemistry"
REFERENCE = ROOT / "build/advanced-chemistry/metal/qm-1mnc-small-screen-v1/calculation"
CANDIDATE = ROOT / "build/advanced-chemistry/metal/qm-1mnc-df-comparison-v2"
CLUSTER = BASE / "prototype-1mnc-inputs"
OUT = BASE / "prototype-1mnc-df-independent-review"


def main():
    evidence = {}

    def read(path):
        path = Path(path).resolve()
        contents = path.read_bytes()
        if not contents or len(contents) != path.stat().st_size:
            raise RuntimeError(f"Incomplete or empty evidence read: {path}")
        record = {"bytes": len(contents), "sha256": hashlib.sha256(contents).hexdigest()}
        if str(path) in evidence and evidence[str(path)] != record:
            raise ValueError(f"Evidence changed during review: {path}")
        evidence[str(path)] = record
        return contents

    def js(path):
        return json.loads(read(path))

    def sha(path):
        read(path)
        return evidence[str(Path(path).resolve())]["sha256"]

    def require(ok, label):
        if not ok:
            raise ValueError(label)

    controller = js(CANDIDATE / "result.json")
    ref = js(REFERENCE / "result.json")
    calc = CANDIDATE / "calculation"
    df = js(calc / "result.json")
    ri, di = js(REFERENCE / "input.json"), js(calc / "input.json")
    rm, dm = js(REFERENCE / "resolved-method.json"), js(calc / "resolved-method.json")
    amap = js(CLUSTER / "small-atom-map.json")
    original = js(CLUSTER / "small-screening-qm-request.json")
    cluster_report = js(CLUSTER / "cluster-input-report.json")
    source = ROOT / ".tools/ambertools/lib/python3.14/site-packages/pymsmt/mcpb/gene_model_files.py"
    native_source = read(source).decode()
    require("def write_sc(" in native_source and "CA atom --> CH3 atom" in native_source,
            "Native cap construction source lacks expected scope")
    require(controller["status"] == "completed", "Controller did not complete")
    require(controller["script_sha256"] == sha(BASE / "compare_1mnc_density_fit.py"), "Controller source hash differs")
    for folder, result in [(REFERENCE, ref), (calc, df)]:
        require(result["accepted"] is True and result["scf"]["converged"] is True,
                "Saved SCF did not converge")
        require(not result.get("optimization", {}).get("requested", False), "Geometry optimization is outside this comparison")
        for artifact, field in [("input.json", "input_sha256"), ("arrays.npz", "arrays_sha256"),
                                ("resolved-method.json", "resolved_method_file_sha256")]:
            require(sha(folder / artifact) == result[field], f"{folder.name} {artifact} hash differs")
    require(sha(calc / "result.json") == controller["result_sha256"], "DF result hash differs")
    require(sha(REFERENCE / "result.json") == controller["source_result_sha256"], "Reference result hash differs")
    require(sha(REFERENCE / "arrays.npz") == controller["source_arrays_sha256"], "Reference array hash differs")
    require(sha(CANDIDATE / "qm_worker.py") == df["worker_sha256"] == controller["worker_sha256"], "DF worker hash differs")
    require(sha(CANDIDATE / "input.json") == df["input_sha256"] == controller["input_sha256"], "DF request hash differs")
    require(ri == original, "Reference request differs from original native-cluster request")
    changes = sorted(k for k in set(ri) | set(di) if ri.get(k) != di.get(k))
    require(changes == ["auxbasis", "density_fit", "initial_checkpoint"], "Unexpected scientific request changes")
    require(di["auxbasis"] == "def2-universal-jkfit" and not ri["density_fit"] and di["density_fit"], "Incorrect DF change")
    require(di["operations"] == ["gradient"] and di["charge"] == 1 and di["spin"] == 0, "Operation or state differs")
    for field in ["basis_resolved", "basis_sha256", "ecp_resolved", "ecp_sha256", "functional_resolved",
                  "functional_resolved_sha256", "grid_level", "spherical_basis", "symmetry"]:
        require(rm[field] == dm[field], f"Resolved method differs in {field}")
    require(dm["spherical_basis"] is True and dm["grid_level"] == 4 and dm["functional_requested"] == "B3LYPG",
            "Declared spherical B3LYPG/grid model differs")
    require(ref["atom_ids"] == df["atom_ids"] == di["atom_ids"] == [a["stable_id"] for a in amap], "Atom-order mismatch")
    require(ref["elements"] == df["elements"] == ri["elements"] == di["elements"], "Element mismatch")
    require(len(amap) == len(set(df["atom_ids"])) == 97, "Incorrect or duplicate atom identities")
    checkpoint = df["initial_checkpoint"]
    require(checkpoint["accepted_as_initial_guess"] is True and checkpoint["maximum_coordinate_difference_bohr"] == 0,
            "Initial checkpoint was not the same geometry")
    for original_file, copied_file, hash_key in [
            (REFERENCE / "result.json", calc / "initial-checkpoint/result.json", "source_result_sha256"),
            (REFERENCE / "arrays.npz", calc / "initial-checkpoint/arrays.npz", "source_arrays_sha256"),
            (REFERENCE / "scf.chk", calc / "initial-checkpoint/scf.chk", "source_checkpoint_sha256"),
            (REFERENCE / "resolved-method.json", calc / "initial-checkpoint/resolved-method.json", "source_resolved_method_sha256")]:
        require(sha(original_file) == sha(copied_file) == checkpoint[hash_key], "Checkpoint copy/hash mismatch")
    require(sha(calc / "scf.chk") == df["checkpoint_sha256"], "Final checkpoint hash differs")
    require(sha(calc / "runtime-manifest.json") == df["runtime_manifest"]["sha256"], "Runtime declaration hash differs")
    with np.load(io.BytesIO(read(REFERENCE / "arrays.npz")), allow_pickle=False) as f:
        ra = {k: f[k].copy() for k in f.files}
    with np.load(io.BytesIO(read(calc / "arrays.npz")), allow_pickle=False) as f:
        da = {k: f[k].copy() for k in f.files}
    require(set(ra) == set(da) == {"coords_bohr", "coords_angstrom", "atomic_numbers", "effective_nuclear_charges", "gradient_hartree_per_bohr"}, "Unexpected array scope")
    for name in ra:
        require(np.isfinite(ra[name]).all() and np.isfinite(da[name]).all(), "Non-finite numerical artifact")
        if name != "gradient_hartree_per_bohr":
            require(np.array_equal(ra[name], da[name]), f"Same-geometry comparison changed {name}")
    require(da["gradient_hartree_per_bohr"].shape == ra["gradient_hartree_per_bohr"].shape == (97, 3), "Gradient shape mismatch")
    input_xyz = np.asarray(ri["coords_angstrom"], dtype=float)
    bohr_constant = df["units"]["bohr_to_angstrom"]
    require(np.array_equal(input_xyz, [a["xyz"] for a in amap]), "Native input/map coordinate mismatch")
    # PySCF internally multiplies Angstrom input by 1/BOHR and converts its
    # saved coordinates back. Verify those actual operations exactly rather
    # than imposing decimal-input identity or allowing coordinate relaxation.
    require(np.array_equal(da["coords_bohr"], input_xyz * (1.0 / bohr_constant)) and
            np.array_equal(da["coords_angstrom"], da["coords_bohr"] * bohr_constant),
            "Input/actual coordinates differ beyond the exact unit roundtrip")
    for result in [ref, df]:
        require(result["units"]["energy"] == "hartree" and result["units"]["gradient_hartree_per_bohr"] == "hartree/bohr", "Unknown numerical units")
        require(result["charge"] == 1 and result["spin_2S"] == 0 and result["nao"] == 778 and result["explicit_electron_count"] == 372, "State/basis dimension mismatch")
        require(np.isfinite(result["energy_hartree"]), "Non-finite energy")
    native_input = read(CLUSTER / "site_small_opt.com").decode()
    native_rows = [line.split() for line in native_input.splitlines()
                   if re.fullmatch(r"[A-Z][a-z]?", (line.split() or [""])[0]) and len(line.split()) == 4]
    require([r[0] for r in native_rows] == di["elements"] and
            np.array_equal([[float(x) for x in r[1:]] for r in native_rows], input_xyz), "Native QM atom order/geometry differs")
    require(" Opt" in native_input and not any(token in native_input.lower() for token in ["modredundant", "freeze", "opt=readfreeze"]), "Native optimization constraints need separate review")
    source_pdb = read(CLUSTER / "site-input.pdb").decode()
    source_atoms = {int(line[6:11]): {"atom_name": line[12:16].strip(), "xyz": [float(line[i:i+8]) for i in (30, 38, 46)]}
                    for line in source_pdb.splitlines() if line.startswith(("ATOM  ", "HETATM"))}
    caps = [a for a in amap if a["role"] == "native MCPB cap"]
    require(len(caps) == 12 and len(amap) - len(caps) == 85, "Cap inventory differs")
    cap_carbons = []
    for a in caps:
        if a["name"] == "CH3":
            source_atom = source_atoms[a["serial"]]
            require(source_atom["atom_name"] == "CA" and np.array_equal(source_atom["xyz"], a["xyz"]), "Cap carbon is not the unchanged native CA boundary")
            cap_carbons.append({"atom_id": a["stable_id"], "source_atom_id": a["source_or_placed_id"], "coords_angstrom": a["xyz"]})
    require(len(cap_carbons) == 3, "Three known histidine cap carbons were not found")
    gradient = da["gradient_hartree_per_bohr"]
    difference = gradient - ra["gradient_hartree_per_bohr"]
    energy_difference = df["energy_hartree"] - ref["energy_hartree"]
    metrics = {"energy_difference_hartree": energy_difference,
        "gradient_max_component_difference_hartree_per_bohr": float(np.abs(difference).max()),
        "gradient_rms_component_difference_hartree_per_bohr": float(np.sqrt(np.mean(difference**2))),
        "gradient_max_vector_difference_hartree_per_bohr": float(np.linalg.norm(difference, axis=1).max()),
        "gradient_rms_vector_difference_hartree_per_bohr": float(np.sqrt(np.mean(np.sum(difference**2, axis=1)))),
        "largest_vector_difference_atom": df["atom_ids"][np.argmax(np.linalg.norm(difference, axis=1))],
        "df_absolute_gradient_max_component_hartree_per_bohr": float(np.abs(gradient).max()),
        "df_absolute_gradient_rms_component_hartree_per_bohr": float(np.sqrt(np.mean(gradient**2))),
        "gradient_difference_norm_over_reference_gradient_norm": float(np.linalg.norm(difference) / np.linalg.norm(ra["gradient_hartree_per_bohr"]))}
    require(metrics["energy_difference_hartree"] == controller["df_minus_conventional_energy_hartree"] and
            metrics["gradient_max_component_difference_hartree_per_bohr"] == controller["maximum_gradient_component_difference_hartree_per_bohr"] and
            metrics["gradient_rms_component_difference_hartree_per_bohr"] == controller["rms_gradient_component_difference_hartree_per_bohr"], "Independent metrics differ from controller")
    readiness = {"status": "model_boundary_choice_required_before_constrained_request", "launch_ready": False,
        "native_request": "Original MCPB small-model Gaussian/GAMESS optimization is unconstrained; no cap-freeze list is inherited.",
        "verified_cap_carbons_at_original_CA": cap_carbons,
        "remaining_native_cap_atoms": [a["stable_id"] for a in caps if a["name"] != "CH3"],
        "candidate_only_not_selected": "Freezing only the three cap carbons would preserve their source CA positions while allowing methyl H relaxation, but is a new constrained-cluster approximation. Freezing all twelve cap atoms adds a different constraint and would preserve unrelaxed native H guesses. Neither policy is silently selected.",
        "specified_model": cluster_report["declared_state"],
        "required_decisions": ["Select and document unconstrained versus exact cap-freeze policy; do not infer it from native input.",
            "Define geometry/coordination/stereochemistry and cap-displacement diagnostics before optimization; convergence alone is not model acceptance.",
            "Specify accepted treatment of missing nonlocal protein/electrostatic environment and full-PLH relaxation for this isolated model.",
            "Choose and freeze numerical convergence/resource limits and subsequent exact-gradient verification; the measured initial-geometry DF error is not an endpoint guarantee."],
        "later_unfinished_stages": ["Converged, checked optimized geometry", "Hessian and parameter extraction/diagnostics", "Larger-cluster ESP/RESP", "Other structural Zn and Ca models", "Full-complex parameter mapping and native energy/force validation", "Unrestrained motion evaluation"],
        "optimization_or_hessian_launched_by_review": False}
    report = {"reviewed_unix": time.time(), "status": "saved_numerical_comparison_independently_reproduced",
        "full_preparation_ready": False, "model_accuracy_validated": False, "native_reruns": 0,
        "checks": {"hash_bound_saved_inputs_and_outputs": True, "same_97_ordered_atom_coordinates_exact": True,
            "same_resolved_orbital_basis_functional_grid_charge_spin": True, "all_numerical_arrays_finite": True,
            "native_cluster_and_checkpoint_mappings_verified": True, "independent_metrics_match_controller_exactly": True},
        "metrics": metrics,
        "input_coordinate_unit_roundtrip": {"actual_both_runs_bohr_and_angstrom_exactly_equal": True,
            "input_to_internal_bohr_exact_formula": "input_angstrom * (1 / pyscf.lib.param.BOHR)",
            "output_angstrom_exact_formula": "internal_bohr * pyscf.lib.param.BOHR",
            "maximum_roundtrip_component_angstrom": float(np.abs(da["coords_angstrom"] - input_xyz).max()),
            "coordinate_relaxation_allowed": False},
        "reference_native_seconds": ref["elapsed_seconds"], "df_native_seconds": df["elapsed_seconds"],
        "reference_scf_cycles": ref["scf"]["cycles"], "df_scf_cycles": df["scf"]["cycles"],
        "method": {"orbital_basis": "PySCF spherical 6-31g*", "functional": "B3LYPG", "grid_level": 4,
            "DF_auxbasis": di["auxbasis"], "charge": 1, "spin_2S": 0, "atom_count": 97, "electron_count": 372, "nao": 778},
        "limitations": ["One fixed geometry with appreciable gradients, not an optimized structure, force-field fit, or physical validation.",
            "DF calculation began from the converged conventional density; wall times are not a controlled speedup comparison.",
            "Worker/runtime builds differ. Both declare the same PySCF/libxc versions and resolved scientific method; loaded binary provenance is not independently re-attested here.",
            "No relative conformer-energy, Hessian, convergence-path, protonation-state or full-complex equivalence claim.",
            "The explicitly selected Zn(II)/HID/hydroxamate state remains a model hypothesis; structural Zn, Ca, protein environment and solvent are excluded.",
            "Spherical basis is a distinct documented method, not a reproduction of Gaussian mixed 6D/7F results."],
        "next_request": readiness, "evidence": evidence, "review_script_sha256": sha(__file__)}
    for path in list(evidence):
        read(path)
    OUT.mkdir(exist_ok=False)
    (OUT / "review.json").write_text(json.dumps(report, indent=2) + "\n")
    (OUT / "optimization-readiness.json").write_text(json.dumps(readiness, indent=2) + "\n")
    (OUT / "native-cap-construction-source.py").write_bytes(read(source))
    print(json.dumps({"status": report["status"], "metrics": metrics, "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
