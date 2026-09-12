#!/usr/bin/env python3
"""Summarize saved DynaMol runs as JSON and a temperature/C-alpha RMSD figure.

Manifest: {"cases": [{"dataRoot": "...", "jobId": "...", "pdbId": "...",
"name": "...", "class": "..."}]}. dataRoot is absolute or relative to the
working directory. Each analysis uses a fresh backend process; no MD is run.
Matplotlib can live in a separate --plot-python environment.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

# Apply before importing numerical libraries in either analysis or plotting.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "DYNAMOL_CPU_THREADS"):
    os.environ[variable] = "1"

ROOT = Path(__file__).resolve().parents[1]
COLORS = ["#087f8c", "#5465ca", "#c77828", "#9659a5", "#c95862"]
METHOD = ("Equal-weight protein Cα RMSD to the first saved production frame. Bonded molecules are made whole with the app's periodic imaging, "
          "using the connected molecule containing the most selected Cα atoms as the anchor and nearby images of other molecules. "
          "One proper Kabsch fit uses all protein Cα atoms together. For oligomers this includes relative subunit motion; "
          "there are no separate chain fits, subunit permutations, or symmetry matching. This is not a diffusion/unwrapping analysis.")
LIMITS = ("These short fixed-volume runs test software and gross numerical behavior. Temperature and RMSD do not establish equilibration, "
          "convergence, correct ligand poses, or biological stability. Different proteins and engines are not ranked by RMSD. "
          "Cross-engine energy comparison is intentionally omitted.")


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def collect_case(case):
    """Called in an isolated child so backend storage resolves one data root."""
    import numpy as np
    root = Path(case["dataRoot"]).resolve()
    job_id = case["jobId"]
    if Path(job_id).name != job_id or not job_id:
        raise ValueError("Invalid jobId")
    folder = root / "jobs" / job_id
    status = json.loads((folder / "status.json").read_text())
    result = {**case, "dataRoot": str(root), "status": status["status"], "engine": status.get("engine"),
              "elapsedSeconds": status.get("elapsed_seconds")}
    if status["status"] != "completed":
        return {**result, "included": False, "reason": "Only completed runs are summarized", "stage": status.get("stage")}
    os.environ["DYNAMOL_DATA_DIR"] = str(root)
    sys.path.insert(0, str(ROOT))
    from backend import storage
    from backend.structural_analysis import structural_analysis, StructuralAnalysisRequest
    settings = json.loads((folder / "config.json").read_text())
    dataset_id = status["dataset_id"]
    metadata = storage.get_dataset(dataset_id)
    analysis = structural_analysis(dataset_id, StructuralAnalysisRequest(kind="rmsd", selection="ca", alignment="ca", periodic="whole", reference_frame=0))
    times = np.asarray(analysis["x"], dtype=float)
    rmsd = np.asarray(analysis["values"], dtype=float)
    with (folder / "energies.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    energy_times = np.asarray([float(row["time_ps"]) for row in rows])
    temperatures = np.asarray([float(row["temperature_k"]) for row in rows])
    expected_steps = round(settings["duration_ps"] * 1000 / settings["timestep_fs"])
    expected_frames = (expected_steps + settings["report_interval"] - 1) // settings["report_interval"] + 1
    if status.get("completed_steps") != expected_steps or len(times) != expected_frames:
        raise ValueError("Completed status or saved frame count disagrees with run configuration")
    if analysis["x_label"] != "Time (ps)" or not len(times) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("Physical production timestamps are unavailable or invalid")
    if abs(times[0]) > 1e-6 or abs(times[-1] - expected_steps * settings["timestep_fs"] / 1000) > 1e-4:
        raise ValueError("Saved trajectory does not cover the configured production duration")
    if len(energy_times) != len(times) or not np.allclose(energy_times, times, atol=1e-4, rtol=1e-6):
        raise ValueError("Native temperature diagnostics are not aligned with saved trajectory frames")
    if not np.isfinite(temperatures).all() or np.any(temperatures <= 0) or not np.isfinite(rmsd).all():
        raise ValueError("Non-finite or nonpositive diagnostics cannot be plotted")
    if settings["solvent"] == "explicit" and not metadata.get("has_unitcell"):
        raise ValueError("An explicit-solvent run lacks the saved periodic cell required for whole-molecule RMSD")
    ca = [metadata["atoms"][i] for i in analysis["atoms"]]
    chains = sorted({atom["chain"] for atom in ca})
    topology = storage.load_physical(dataset_id).topology
    selected = set(analysis["atoms"])
    component_ca_counts = sorted([sum(atom.index in selected for atom in molecule) for molecule in topology.find_molecules()], reverse=True)
    component_ca_counts = [count for count in component_ca_counts if count]
    def summary(values):
        return {"mean": float(values.mean()), "stdPopulation": float(values.std()), "minimum": float(values.min()),
                "maximum": float(values.max()), "final": float(values[-1])}
    output = root / "datasets" / dataset_id
    evidence = {str(path.relative_to(root)): sha256(path) for path in [folder / "config.json", folder / "energies.csv", folder / "provenance.json", output / "physical.npz", output / "topology.pdb"]}
    return {**result, "included": True, "datasetId": dataset_id, "durationPs": float(times[-1]), "frames": len(times),
            "atoms": metadata["n_atoms"], "proteinCaAtoms": len(ca), "savedProteinChainLabels": chains,
            "connectedProteinComponents": len(component_ca_counts), "caAtomsPerConnectedComponent": component_ca_counts,
            "targetTemperatureK": settings["temperature_k"], "temperatureK": summary(temperatures), "caRmsdAngstrom": summary(rmsd),
            "series": {"timePs": np.round(times, 6).tolist(), "temperatureK": np.round(temperatures, 5).tolist(), "caRmsdAngstrom": np.round(rmsd, 6).tolist()},
            "structuralAnalysisMethod": analysis["method"], "warnings": analysis["warnings"],
            "evidenceSha256": evidence, "analysisImplementationSha256": sha256(ROOT / "backend/structural_analysis.py")}


def render(report, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.labelcolor": "#394755",
                         "text.color": "#152a3b", "xtick.color": "#667684", "ytick.color": "#667684", "svg.fonttype": "none"})
    background = "#f3f6f8"
    fig, axes = plt.subplots(1, 2, figsize=(13, 7.2), facecolor=background)
    fig.subplots_adjust(left=.075, right=.968, top=.725, bottom=.25, wspace=.24)
    cases = [case for case in report["cases"] if case.get("included")]
    durations = sorted({round(case["durationPs"], 3) for case in cases})
    duration_label = f"{durations[0]:g} ps" if len(durations) == 1 else "Saved trajectory"
    fig.text(.068, .944, f"DynaMol  /  {duration_label} checks", fontsize=24, fontweight="bold")
    fig.text(.068, .897, "Native temperatures and protein motion from saved production frames", fontsize=12, color="#607482")
    handles = []
    for case in cases:
        color = case["color"]
        engine_label = {"openmm": "OpenMM", "gromacs": "GROMACS"}.get(case['engine'], case['engine'])
        label = f"{case.get('pdbId', '')} · {case.get('name', case['jobId'])} · {engine_label}"
        handles.append(Line2D([0], [0], color=color, lw=2.5, label=label))
        series = case["series"]
        axes[0].plot(series["timePs"], series["temperatureK"], color=color, lw=1.7, alpha=.86)
        axes[1].plot(series["timePs"], series["caRmsdAngstrom"], color=color, lw=2.1, alpha=.95)
    if handles:
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.062, .87), ncol=2 if len(handles) > 1 else 1,
                   frameon=False, fontsize=9.2, handlelength=2.7, columnspacing=2.8)
    for index, (axis, title, ylabel) in enumerate(zip(axes, ["Temperature", "Cα RMSD after alignment"], ["Temperature (K)", "RMSD to first frame (Å)"])):
        axis.set_facecolor("white")
        axis.set_title(title, loc="left", pad=15, fontweight="bold", fontsize=13)
        axis.set_xlabel("Production time (ps)", labelpad=10)
        axis.set_ylabel(ylabel, labelpad=10)
        axis.grid(axis="y", color="#e2e8ed", linewidth=.7)
        axis.set_axisbelow(True)
        axis.tick_params(length=0, pad=7)
        for spine in axis.spines.values():
            spine.set_visible(False)
        if durations:
            axis.set_xlim(0, max(durations))
        if index == 1:
            axis.set_ylim(bottom=0)
        if not cases:
            axis.text(.5, .5, "No completed runs yet", ha="center", va="center", transform=axis.transAxes, color="#667684")
    targets = {case["targetTemperatureK"] for case in cases}
    if len(targets) == 1:
        target = next(iter(targets))
        axes[0].axhline(target, color="#7a8994", lw=.8, ls=(0, (3, 4)), zorder=0)
        axes[0].text(.98, .98, f"Target {target:g} K", transform=axes[0].transAxes, ha="right", va="top", fontsize=9, color="#687a87")
    fig.text(.075, .145, "RMSD method: molecules made whole; one fit over all protein Cα atoms. Relative subunit motion remains included.", fontsize=10, color="#516573")
    fig.text(.075, .095, "Short-run diagnostics do not establish equilibration, convergence, or biological stability.", fontsize=10, color="#516573")
    count = len(cases)
    omitted = len(report["cases"]) - count
    footer = f"{count}/{len(report['cases'])} runs included · No display-coordinate analysis · No cross-engine energy comparison"
    if omitted:
        footer += f" · {omitted} pending or unavailable; see JSON"
    fig.text(.075, .04, footer, fontsize=9, color="#7c8b96")
    for extension in ("png", "svg"):
        fig.savefig(output_dir / f"pressure-test-summary.{extension}", dpi=180, facecolor=background)
    plt.close(fig)
    return {"matplotlib": matplotlib.__version__}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--backend-python", default=str(ROOT / ".venv/bin/python"))
    parser.add_argument("--plot-python", default=str(ROOT / ".tools/pressure-plot/bin/python"))
    parser.add_argument("--collect-case", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--render-json", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.collect_case:
        try:
            result = collect_case(json.load(sys.stdin))
        except Exception as exc:
            result = {"included": False, "status": "analysis-failed", "reason": f"{type(exc).__name__}: {exc}"}
        print(json.dumps(result, allow_nan=False))
        return 0
    if args.render_json:
        print(json.dumps(render(json.loads(args.render_json.read_text()), args.output_dir)))
        return 0
    if args.manifest is None or args.output_dir is None:
        parser.error("--manifest and --output-dir are required")
    manifest = json.loads(args.manifest.read_text())
    cases = manifest["cases"] if isinstance(manifest, dict) else manifest
    if not isinstance(cases, list) or not cases:
        parser.error("Manifest must contain a nonempty cases list")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {"createdAt": dt.datetime.now(dt.timezone.utc).isoformat(), "method": METHOD, "limitations": LIMITS,
              "cpuThreads": 1, "manifestSha256": sha256(args.manifest), "scriptSha256": sha256(Path(__file__)), "cases": []}
    for index, case in enumerate(cases):
        try:
            result = subprocess.run([args.backend_python, str(Path(__file__).resolve()), "--collect-case"], input=json.dumps(case),
                                    text=True, capture_output=True, env=os.environ.copy(), timeout=180)
            record = json.loads(result.stdout)
        except (OSError, subprocess.TimeoutExpired) as exc:
            record = {"included": False, "status": "analysis-failed", "reason": str(exc)}
        except (ValueError, TypeError):
            record = {"included": False, "status": "analysis-failed", "reason": result.stderr[-2000:] or "Analysis child returned invalid JSON"}
        report["cases"].append({**case, **record, "color": COLORS[index % len(COLORS)]})
    report["includedRuns"] = sum(case.get("included", False) for case in report["cases"])
    report["complete"] = report["includedRuns"] == len(cases)
    target = output_dir / "pressure-test-summary.json"
    target.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    plotted = subprocess.run([args.plot_python, str(Path(__file__).resolve()), "--render-json", str(target), "--output-dir", str(output_dir)],
                             text=True, capture_output=True, env=os.environ.copy(), timeout=120)
    if plotted.returncode:
        raise RuntimeError("Plotting failed; install matplotlib in --plot-python environment. " + plotted.stderr[-3000:])
    report["plotRuntime"] = json.loads(plotted.stdout)
    target.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"included": report["includedRuns"], "total": len(cases), "complete": report["complete"], "outputDir": str(output_dir)}))
    return 0 if report["includedRuns"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
