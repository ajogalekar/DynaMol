#!/usr/bin/env python3
"""Validate one PDB entry through the real DynaMol workflow in an isolated data dir.

Fetch -> (monomer if multi-chain/large) -> prepare -> short MD (implicit for a
protein, explicit for a retained-ligand complex). Classifies the outcome as
passed / expected-block (covalent, metal-site model, oversize, >12-residue gap,
nucleic acid, unsupported chemistry -- the V2/out-of-scope set) / FAILED (an
unexpected crash, stereochemistry inversion, NaN, geometry rejection, etc.).

Run as its own process so DYNAMOL_DATA_DIR isolates it; prep/solvation/dynamics
spawn their own workers that inherit the same data dir.
"""
from __future__ import annotations
import argparse, json, os, re, sys, time, traceback
from pathlib import Path

# Patterns that mean "correctly refused, out of v1 scope" -- not a workflow bug.
EXPECTED = [
    r"covalent", r"specialized force-field model", r"CCD .*mismatch", r"heavy-atom identity",
    r"unsupported element", r"radical", r"exceeds? .*12-residue", r"up to 12 per gap",
    r"gap of \d+ residues", r"100,?000-atom", r"atom limit", r"solvent box may exceed",
    r"nucleic", r"no supported covalent amino-acid template", r"modified or unnatural",
    r"incomplete ligand", r"inconsistent stereochemical", r"terminal", r"supply a complete model",
    r"needs a specialized", r"coordinate limit", r"disk space", r"across the periodic boundary",
    r"no protein chain", r"requires another preparation workflow",
    # Deposited-coordinate defect: a genuinely inverted/near-planar standard-residue
    # centre in the INPUT (e.g. the 1974 entry 1LYZ, ARG14 CA). Distinct from the
    # MD-time guard message ("inverted or flattened ... Dynamics stopped").
    r"inverted or near-planar standard residue", r"The model is rejected",
]


def classify(msg: str) -> str:
    m = (msg or "").lower()
    return "expected_block" if any(re.search(p, m) for p in EXPECTED) else "FAILED"


def poll(jobs, job_id, timeout, label, result):
    start = time.monotonic()
    last = None
    while True:
        state = jobs.get_job(job_id)
        if state.get("stage") != last:
            last = state.get("stage")
            print(f"  [{label}] {state.get('status')} :: {last}", flush=True)
        if state["status"] in {"completed", "failed", "cancelled", "interrupted"}:
            return state
        if time.monotonic() - start > timeout:
            try:
                jobs.cancel_job(job_id)
            except Exception:
                pass
            state = jobs.get_job(job_id)
            state["error"] = f"TIMEOUT after {timeout}s during {label}: " + (state.get("error") or "")
            state["status"] = "failed"
            return state
        time.sleep(2)


def run(pdb_id, workdir, md_ps, timeout):
    data = Path(workdir) / "data"
    data.mkdir(parents=True, exist_ok=True)
    os.environ.update(DYNAMOL_DATA_DIR=str(data), DYNAMOL_CPU_THREADS="2",
                      OPENMM_CPU_THREADS="2", OPENBLAS_NUM_THREADS="2", OMP_NUM_THREADS="2",
                      PYTHONDONTWRITEBYTECODE="1")
    from backend import sources, monomers, preparation, jobs, storage
    from backend.monomers import MonomerRequest
    from backend.models import SimulationConfig

    result = {"pdb_id": pdb_id, "phase": "fetch", "outcome": "FAILED", "error": "",
              "stages": {}, "started_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    t0 = time.monotonic()
    try:
        meta = sources.fetch_structure("pdb", pdb_id)
        dataset = meta["id"]
        result["raw_atoms"] = meta.get("n_atoms")

        insp = preparation.inspect_preparation(dataset, 7.0)
        chains = monomers.list_monomers(dataset)["chains"]
        result["n_chains"] = len(chains)
        result["ligands"] = [l.get("key") for l in insp.get("ligands", [])]
        result["blockers_at_inspection"] = insp.get("blockers", [])
        protein_atoms = insp.get("protein_atoms", meta.get("n_atoms", 0))

        # Select a monomer when multi-chain or when a solvated whole structure
        # would obviously blow the atom budget. Keep associated molecules so a
        # ligand travels with its chain; fall back to protein-only on a PBC clash.
        result["phase"] = "monomer"
        if len(chains) > 1 or protein_atoms > 12000:
            idx = monomers.list_monomers(dataset).get("recommended_chain_index") or 0
            for keep in (True, False):
                try:
                    m = monomers.create_monomer(dataset, MonomerRequest(chain_index=idx, keep_associated_molecules=keep))
                    dataset = m["id"]
                    result["monomer"] = {"chain_index": idx, "kept_ligands": keep}
                    break
                except Exception as exc:  # PBC contact etc.
                    if keep and re.search(r"periodic boundary", str(exc)):
                        continue
                    raise

        result["phase"] = "prepare"
        prep = preparation.submit_preparation({
            **{k: v for k, v in {
                "dataset_id": dataset, "name": f"validate {pdb_id}", "ph": 7.0,
                "add_missing_atoms": True, "build_missing_residues": True,
                "optimize_sidechains": True, "remove_waters": True, "remove_heterogens": False,
            }.items()}})
        pstate = poll(jobs, prep["id"], timeout, "prepare", result)
        result["stages"]["prepare"] = {"status": pstate["status"], "seconds": round(time.monotonic() - t0, 1)}
        if pstate["status"] != "completed":
            result["error"] = pstate.get("error", "preparation did not complete")
            result["outcome"] = classify(result["error"])
            return result
        prepared = pstate["dataset_id"]

        # Short MD. Try implicit first (fast, protein-only); if the prepared
        # system needs explicit water (retained ligand/ions/modified residues),
        # solvate and rerun explicit.
        result["phase"] = "dynamics"
        report_interval = max(200, int(md_ps * 1000 / 2 / 6))  # ~6 saved frames

        def submit_md(ds, solvent):
            cfg = SimulationConfig(dataset_id=ds, engine="openmm", name=f"MD {pdb_id}",
                                   duration_ps=md_ps, timestep_fs=2, report_interval=report_interval,
                                   minimize=True, equilibration_steps=100, solvent=solvent,
                                   temperature_k=300, padding_nm=1)
            return jobs.submit_job(cfg)

        def solvate_explicit():
            sstate = poll(jobs, preparation.submit_solvation(prepared, {"padding_nm": 1, "ph": 7})["id"], timeout, "solvate", result)
            result["stages"]["solvate"] = {"status": sstate["status"]}
            if sstate["status"] != "completed":
                result["error"] = sstate.get("error", "solvation failed")
                result["outcome"] = classify(result["error"]); result["phase"] = "solvate"
                return None
            return sstate["dataset_id"]

        # Prefer explicit TIP3P/PME for large protein-only systems (implicit GBn2
        # is NoCutoff/O(N^2) and slow on CPU at that size), matching the app's
        # preparation recommendation; otherwise implicit first, with an explicit
        # fallback for retained-ligand/ion/modified complexes.
        recommend_explicit = bool(storage.get_dataset(prepared).get("preparation", {}).get("recommend_explicit_solvent"))
        if recommend_explicit:
            md_ds = solvate_explicit()
            if md_ds is None:
                return result
            solvent = "explicit"
            job = submit_md(md_ds, "explicit")
        else:
            md_ds, solvent = prepared, "implicit"
            try:
                job = submit_md(prepared, "implicit")
            except Exception as exc:
                if re.search(r"explicit|protein-only|preview", str(exc), re.I):
                    md_ds = solvate_explicit()
                    if md_ds is None:
                        return result
                    solvent = "explicit"
                    job = submit_md(md_ds, "explicit")
                else:
                    raise
        result["solvent"] = solvent
        mstate = poll(jobs, job["id"], timeout, "dynamics", result)
        result["stages"]["dynamics"] = {"status": mstate["status"], "seconds": round(time.monotonic() - t0, 1)}
        if mstate["status"] == "completed":
            result["outcome"] = "passed"; result["phase"] = "complete"
        else:
            result["error"] = mstate.get("error", "dynamics did not complete")
            result["outcome"] = classify(result["error"])
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()[-2000:]
        result["outcome"] = classify(result["error"])
    finally:
        result["elapsed_seconds"] = round(time.monotonic() - t0, 1)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb_id")
    ap.add_argument("workdir")
    ap.add_argument("--md-ps", type=float, default=10.0)
    ap.add_argument("--timeout", type=int, default=1200)
    args = ap.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    result = run(args.pdb_id, args.workdir, args.md_ps, args.timeout)
    out = Path(args.workdir) / "result.json"
    out.write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(f"RESULT {args.pdb_id} {result['outcome']} phase={result['phase']} {result.get('error','')[:160]}", flush=True)
    sys.exit(0 if result["outcome"] in {"passed", "expected_block"} else 1)


if __name__ == "__main__":
    main()
