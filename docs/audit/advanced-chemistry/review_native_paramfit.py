"""Corrected, read-only analysis of existing native Paramfit qualification files.

No native executable is invoked. Original runner failures and all native inputs
and outputs are checked against their saved hashes and are never rewritten.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
AUDIT = ROOT / "build/advanced-chemistry/paramfit-native-audit-v1"
CACHE = Path.home() / ".cache/dynamol-native-audits/paramfit-v1"
OUT = ROOT / "docs/audit/advanced-chemistry/paramfit-native-corrected-review"


def read(path):
    path = Path(path)
    data = path.read_bytes()
    if len(data) != path.stat().st_size or not data:
        raise ValueError(f"Unavailable or empty evidence file: {path}")
    return data


def sha(path):
    return hashlib.sha256(read(path)).hexdigest()


def matrix(path, columns):
    rows = []
    for line_no, line in enumerate(read(path).decode().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != columns:
            raise ValueError(f"{path}:{line_no}: expected {columns} numeric columns")
        rows.append([float(x) for x in fields])
    result = np.array(rows, dtype=float)
    if result.ndim != 2 or result.shape[1] != columns or not np.isfinite(result).all():
        raise ValueError(f"Malformed/nonfinite numeric table: {path}")
    return result


def saved_column(path):
    values = []
    for line in read(path).decode().splitlines():
        if not line.strip() or line.strip() == "Amber+K":
            continue
        fields = line.split()
        if len(fields) != 1:
            raise ValueError("The upstream energy.out reference must contain one selected column")
        values.append(float(fields[0]))
    values = np.array(values)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError("Invalid upstream selected-column reference")
    return values


def parameter_section(text, heading):
    section = text.split(heading, 1)[1]
    records = {}
    pattern = re.compile(r"(?P<field>Kr|r_eq|Kt|th_eq|Kp|Np|Phase)\s*=\s*(?P<value>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)")
    for line in section.splitlines():
        if line.lstrip().startswith("(* means parameter"):
            continue
        match = re.match(r"\s*(IMP\s+)?\(([^)]+)\)(.*)", line)
        if not match:
            if records and ("---" in line or "Called the fitness" in line):
                break
            continue
        atoms = tuple(x.strip() for x in match.group(2).split("-"))
        values = {m.group("field"): float(m.group("value")) for m in pattern.finditer(match.group(3))}
        kind = "bond" if "Kr" in values else "angle" if "Kt" in values else "improper" if match.group(1) else "proper"
        key = "|".join([kind, *atoms, str(values.get("Np", ""))])
        expected = {"Kr", "r_eq"} if kind == "bond" else {"Kt", "th_eq"} if kind == "angle" else {"Kp", "Np", "Phase"}
        if set(values) != expected or key in records or not all(math.isfinite(x) for x in values.values()):
            raise ValueError("Incomplete or duplicate printed parameter identity")
        records[key] = values
    if len(records) != 30:
        raise ValueError("Expected 7 bonds, 10 angles and 13 printed torsion terms")
    return records


def parameter_difference(a, b):
    if set(a) != set(b):
        raise ValueError("Printed parameter inventories differ")
    changes = []
    for key in a:
        for field in a[key]:
            raw = b[key][field] - a[key][field]
            difference = math.remainder(raw, 360.0) if field == "Phase" else raw
            if abs(difference) > 1e-10:
                changes.append({"parameter": key, "field": field, "reference": a[key][field],
                                "actual": b[key][field], "signed_difference": difference,
                                "raw_signed_difference": raw})
    return changes


def main():
    OUT.mkdir(exist_ok=False)
    old_path = AUDIT / "result.json"
    old_hash = sha(old_path)
    original = json.loads(read(old_path))
    assert original["running"] is False and len(original["cases"]) == 2
    report = {"status": "corrected_postanalysis_only", "native_reruns": False,
              "chemical_model_validated": False, "original_failed_analysis": str(old_path),
              "original_failed_analysis_sha256": old_hash, "cases": [], "source_evidence": {}}
    counts = {"dihedral_least_squares_NMA": 500, "simplex_perfect_fit_NMA": 1000}
    for old in original["cases"]:
        name = old["name"]; path = CACHE / name; count = counts[name]
        assert old["returncode"] == 0 and old["input_unchanged"] is True
        assert old["status"] == "failed" and old["error"]["message"] == "Unexpected native energy-table format"
        hashes = {}
        for filename, expected in {**old["input_sha256"], **old["output_sha256"]}.items():
            observed = sha(path / filename)
            if observed != expected:
                raise ValueError(f"Original native evidence changed: {name}/{filename}")
            hashes[filename] = observed
        a = matrix(path / "energy.dat", 4)
        if a.shape != (count, 4) or not np.array_equal(a[:, 0], np.arange(count)):
            raise ValueError("Native row IDs/count do not match the frozen fixture")
        control = read(path / "Job_Control.in").decode()
        assert re.search(rf"^NSTRUCTURES={count}$", control, re.MULTILINE)
        assert "QM_ENERGY_UNITS=KCALMOL" in control
        target_path = path / ("amber_energy" if count == 500 else "mdcrd_creation/amber_energy.dat")
        targets = matrix(target_path, 1)[:, 0]
        if not np.array_equal(np.sort(targets[:count]), a[:, 2]):
            raise ValueError("Saved targets disagree with the first requested input structures after declared sorting")
        residual = a[:, 1] - a[:, 2]
        log = read(path / "native.log").decode()
        assert "Program Execution Completed" in log
        initial = parameter_section(log, "INITIAL PARAMETERS")
        final = parameter_section(log, "FINAL PARAMETERS")
        historical_log = read(path / "saved_output/prog_out.saved").decode()
        historical_parameters = parameter_section(historical_log, "FINAL PARAMETERS")
        selected = saved_column(path / "saved_output/energy.out.saved")
        if len(selected) != count:
            raise ValueError("Selected upstream reference column has a different row count")
        case = {"name": name, "structures": count, "columns": 4,
                "native_column_header": "Num, Amber+K, Quantum, Initial Amber+K",
                "all_four_columns_finite": True, "row_ids_exact_zero_through_n_minus_one": True,
                "original_native_returncode": old["returncode"], "original_analysis_error_preserved": old["error"],
                "input_targets_count": len(targets), "sorted_first_requested_targets_exact": True,
                "energy_units": "kcal/mol", "energy_print_decimal_places": 5,
                "fitted_minus_target_rmse_kcal_mol": float(np.sqrt(np.mean(residual**2))),
                "fitted_minus_target_max_abs_kcal_mol": float(np.abs(residual).max()),
                "fitted_minus_target_mean_kcal_mol": float(residual.mean()),
                "fitted_minus_target_sum_squared_kcal2_mol2": float(np.sum(residual**2)),
                "target_range_kcal_mol": [float(a[:, 2].min()), float(a[:, 2].max())],
                "reported_fourth_column_range": [float(a[:, 3].min()), float(a[:, 3].max())],
                "fourth_column_used_as_fit_or_initial_energy": False,
                "printed_parameter_count": len(final),
                "printed_initial_to_final_parameter_changes": parameter_difference(initial, final),
                "historical_log_parameter_comparison": {"comparable_regression": False,
                    "reason": "The preserved historical log uses a different structure count and/or fitted-parameter selection.",
                    "current_displayed_fit_dimensions": int(re.search(r"Total dimensions of fit =\s*(\d+)", log).group(1)),
                    "historical_fit_dimensions": int(re.search(r"Total dimensions of fit =\s*(\d+)", historical_log).group(1)),
                    "diagnostic_parameter_differences_only": parameter_difference(historical_parameters, final)},
                "warnings": [{"line": i, "text": line.strip()} for i, line in enumerate(log.splitlines(), 1)
                             if re.search(r"\bwarning\b", line, re.IGNORECASE)],
                "original_evidence_hashes_verified": hashes}
        if count == 500:
            legacy = matrix(path / "saved_output/energy.dat.saved", 3)
            assert len(legacy) == 100
            case["legacy_energy_dat_reference"] = {"rows": len(legacy), "columns": 3,
                "compared_as_current_fitted_reference": False, "reason": "Older 100-structure artifact is incompatible with this 500-structure run."}
            case["upstream_current_regression"] = {"energy_out_column": "Quantum target (tab field3 in upstream run_test.x)",
                "target_column_exact": bool(np.array_equal(selected, a[:, 2])),
                "frcmod_exact_bytes": read(path / "frcmod") == read(path / "saved_output/frcmod.saved"),
                "caveat": "The selected energy column checks target/sorting, not the fitted energy; the independently saved frcmod supplies fitted-parameter regression evidence."}
            assert case["upstream_current_regression"]["target_column_exact"] and case["upstream_current_regression"]["frcmod_exact_bytes"]
        else:
            difference = np.abs(a[:, 1] - selected)
            # Exact native ndiff.awk rule for these finite, nonzero values.
            relative = difference / np.minimum(np.abs(a[:, 1]), np.abs(selected))
            case["upstream_current_regression"] = {"energy_out_column": "Amber+K fitted energy (tab field2)",
                "maximum_fitted_energy_difference_kcal_mol": float(difference.max()),
                "maximum_native_ndiff_relative_error": float(relative.max()),
                "upstream_relative_limit": 5e-6, "upstream_relative_criterion_met": bool(np.all(relative <= 5e-6))}
            assert case["upstream_current_regression"]["upstream_relative_criterion_met"]
            case["frcmod_status"] = "Not requested by current Job_Control.in; native final parameters are preserved in native.log. No missing-output failure is inferred."
            assert "WRITE_FRCMOD" not in control and not (path / "frcmod").exists()
            source_parameter_file = AUDIT / "source/amber24_src/AmberTools/test/paramfit" / name / "prms.in"
            case["parameter_selection_input"] = {"path": str(path / "prms.in"), "sha256": sha(path / "prms.in"),
                "current_copy_matches_pinned_source": read(path / "prms.in") == read(source_parameter_file),
                "historical_pre_run_hash_available": False,
                "note": "Original runner did not hash prms.in; current source identity and native six-dimension log corroborate selection, but do not create a missing historical hash."}
            case["simplex_native_convergence"] = {"cycles": 150, "function_evaluations_at_convergence": 269,
                "reported_fitness_calls_total": 270, "criterion": 0.1, "ratio": 2.0978e-6,
                "initial_logged_sum_squares_kcal2_mol2": 3.4784e-6}
        extras = ["saved_output/energy.out.saved", "saved_output/prog_out.saved", "run_test.x", "info.txt"]
        if count == 500:
            extras += ["saved_output/energy.dat.saved", "saved_output/frcmod.saved"]
        else:
            extras += ["prms.in"]
        case["additional_current_source_hashes"] = {name: sha(path / name) for name in extras}
        (OUT / f"{name}-printed-parameters.json").write_text(json.dumps({"initial": initial, "final": final}, indent=2) + "\n")
        report["cases"].append(case)
    sources = [AUDIT / "source-manifest.json", ROOT / "docs/audit/advanced-chemistry/check_native_paramfit.py",
               ROOT / ".tools/ambertools/test/ndiff.awk", ROOT / ".tools/ambertools/test/dacdif"]
    sources += [AUDIT / "source/amber24_src/AmberTools/src/paramfit" / name for name in
                ["write_input.c", "eval_amber_std.c", "fitting_control.c", "dihedral_fitting.c"]]
    report["source_evidence"] = {str(path): sha(path) for path in sources}
    report["reviewer_script_sha256"] = sha(__file__)
    report["fourth_column_source_interpretation"] = {
        "writer": "write_input.c::write_energy prints init_energy+K under Initial Amber+K.",
        "sum_squares_path": "eval_amber_std.c starts individual_sum at K-target, adds Amber terms, then stores that residual in init_energy. It is not the initial Amber energy.",
        "least_squares_path": "The native LS log reports zero fitness calls and the saved fourth column is all zero. No initial-energy curve is derived from it.",
        "action": "Use Amber+K minus Quantum for actual fit residuals. Retain the misleading original header; do not silently reinterpret or overwrite native output."}
    assert sha(old_path) == old_hash
    assert sha(original["executable"]) == original["executable_sha256"]
    report["original_failed_analysis_and_executable_unchanged"] = True
    (OUT / "review.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": report["status"], "native_reruns": False, "cases": [
        {k: c[k] for k in ["name", "structures", "fitted_minus_target_rmse_kcal_mol", "fitted_minus_target_max_abs_kcal_mol"]}
        for c in report["cases"]]}, indent=2))


if __name__ == "__main__":
    main()
