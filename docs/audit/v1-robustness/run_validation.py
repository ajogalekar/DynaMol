#!/usr/bin/env python3
"""Drive validate_structure.py across a diverse in-scope PDB panel in parallel.

Each entry runs in its own isolated data dir through the real fetch -> prepare
-> short-MD workflow. Passed / expected-block workspaces are pruned to save
disk; FAILED ones are kept for inspection. Writes summary.json + REPORT.md.
"""
from __future__ import annotations
import argparse, concurrent.futures as cf, json, shutil, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PY = REPO / ".venv/bin/python"

# (pdb_id, target_class, ligand/chemistry note). Diverse classes and ligand
# chemistry; excludes the V2/out-of-scope set (covalent inhibitors, heme,
# bonded metalloenzymes, membrane proteins, nucleic acids, >12-residue gaps).
PANEL = [
    # Tyrosine kinases (heteroaromatic ATP-site inhibitors)
    ("2HYY", "tyrosine kinase", "imatinib"), ("1M17", "tyrosine kinase", "erlotinib"),
    ("1IEP", "tyrosine kinase", "STI-571"), ("3LCK", "tyrosine kinase", "apo Lck"),
    # Ser/Thr kinases (+ nucleotide / staurosporine)
    ("1AQ1", "ser/thr kinase", "staurosporine"), ("1FIN", "ser/thr kinase", "CDK2/cyclin"),
    ("1ATP", "ser/thr kinase", "ATP + Mg"), ("1STC", "ser/thr kinase", "staurosporine"),
    ("3HEG", "ser/thr kinase", "sorafenib"), ("1P38", "ser/thr kinase", "apo p38"),
    # Aspartic proteases (peptidomimetics)
    ("1HVR", "aspartic protease", "cyclic urea XK263"), ("1HSG", "aspartic protease", "indinavir"),
    ("1HPX", "aspartic protease", "KNI-272"), ("3PSG", "aspartic protease", "pepsinogen"),
    # Serine proteases (amidines / small inhibitors)
    ("3PTB", "serine protease", "benzamidine"), ("1DWD", "serine protease", "thrombin inhibitor"),
    ("1ELA", "serine protease", "elastase"), ("1FJS", "serine protease", "factor Xa inhibitor"),
    ("1BRA", "serine protease", "apo trypsin"),
    # Nuclear receptors (steroids / TZDs)
    ("3ERT", "nuclear receptor", "4-OH-tamoxifen"), ("1E3G", "nuclear receptor", "androgen R"),
    ("1FM6", "nuclear receptor", "rosiglitazone"), ("2AM9", "nuclear receptor", "androgen R"),
    # Immunophilins (macrolides)
    ("1FKF", "immunophilin", "FK506"), ("2CPL", "immunophilin", "cyclophilin A"),
    ("1FKB", "immunophilin", "rapamycin"),
    # Chaperones (HSP90 inhibitors)
    ("1YET", "chaperone", "geldanamycin"), ("1UYD", "chaperone", "PU3"), ("1BYQ", "chaperone", "ADP"),
    # Small / model folds (apo, varied topology)
    ("1PGB", "model fold", "protein G B1"), ("5PTI", "model fold", "BPTI"),
    ("1CRN", "model fold", "crambin"), ("1SHG", "model fold", "SH3"),
    ("2CI2", "model fold", "CI2"), ("1VII", "model fold", "villin"),
    ("1UBQ", "model fold", "ubiquitin"), ("1CSP", "model fold", "cold-shock"),
    ("1L2Y", "model fold", "trp-cage"),
    # Hydrolases
    ("2LZM", "hydrolase", "T4 lysozyme"), ("1LYZ", "hydrolase", "HEW lysozyme"),
    ("7RSA", "hydrolase", "RNase A"), ("1AKI", "hydrolase", "lysozyme"),
    # Oxidoreductases / isomerases (non-heme)
    ("4DFR", "oxidoreductase", "methotrexate"), ("1RX2", "oxidoreductase", "folate/NADP"),
    ("1TIM", "isomerase", "TIM"), ("2YPI", "isomerase", "TIM + PGA"),
    # Binding proteins (sugars / vitamins / biotin)
    ("1STP", "binding protein", "biotin"), ("1ANF", "binding protein", "maltose"),
    ("2DRI", "binding protein", "ribose"), ("1CBS", "binding protein", "retinoic acid"),
    ("1RBP", "binding protein", "retinol"),
    # Bromodomains
    ("3MXF", "bromodomain", "BRD4 inhibitor"), ("4LYW", "bromodomain", "BRD4 + JQ1"),
    ("2OSS", "bromodomain", "CBP"),
    # Phosphatase
    ("2HNP", "phosphatase", "PTP1B"), ("1PTY", "phosphatase", "PTP1B + inhibitor"),
    # Others (diverse chemistry)
    ("1EVE", "hydrolase", "acetylcholinesterase + tacrine"),
    ("2CNA", "lectin", "concanavalin A (Ca/Mn ions)"),
    ("1CBR", "binding protein", "retinoic acid"),
    ("1DHF", "oxidoreductase", "human DHFR"),
]


def run_one(entry, root, md_ps, timeout):
    pdb, klass, note = entry
    workdir = root / pdb
    workdir.mkdir(parents=True, exist_ok=True)
    log = workdir / "driver.log"
    with log.open("w") as stream:
        proc = subprocess.run([str(PY), "-B", str(HERE / "validate_structure.py"), pdb, str(workdir),
                               "--md-ps", str(md_ps), "--timeout", str(timeout)],
                              cwd=str(REPO), stdout=stream, stderr=subprocess.STDOUT)
    result_path = workdir / "result.json"
    if result_path.is_file():
        result = json.loads(result_path.read_text())
    else:
        result = {"pdb_id": pdb, "outcome": "FAILED", "phase": "harness",
                  "error": f"no result.json (exit {proc.returncode}); see driver.log"}
    result["class"], result["note"] = klass, note
    # Prune bulk data for non-failures; keep failures for inspection.
    if result.get("outcome") != "FAILED":
        shutil.rmtree(workdir / "data", ignore_errors=True)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path.home() / ".cache/dynamol-research/v1-robustness")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--md-ps", type=float, default=10.0)
    ap.add_argument("--timeout", type=int, default=1200)
    ap.add_argument("--ids", default="", help="comma-separated subset for a pilot")
    args = ap.parse_args()
    panel = [e for e in PANEL if not args.ids or e[0] in set(args.ids.split(","))]
    root = args.out
    root.mkdir(parents=True, exist_ok=True)
    started = time.time()
    results = []
    print(f"Validating {len(panel)} structures, {args.workers} parallel, {args.md_ps:g} ps MD each", flush=True)
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(run_one, e, root, args.md_ps, args.timeout): e for e in panel}
        for fut in cf.as_completed(futures):
            r = fut.result()
            results.append(r)
            done = len(results)
            print(f"[{done}/{len(panel)}] {r['pdb_id']} ({r.get('class')}) -> {r['outcome']} "
                  f"phase={r.get('phase')} {('| ' + r.get('error','')[:120]) if r.get('error') else ''}", flush=True)
            _write_summary(root, results, started, len(panel))
    _write_summary(root, results, started, len(panel), final=True)
    fails = [r for r in results if r["outcome"] == "FAILED"]
    print(f"\nDONE: {sum(r['outcome']=='passed' for r in results)} passed, "
          f"{sum(r['outcome']=='expected_block' for r in results)} expected-block, {len(fails)} FAILED "
          f"in {(time.time()-started)/60:.1f} min", flush=True)


def _write_summary(root, results, started, total, final=False):
    by = {}
    for r in results:
        by.setdefault(r["outcome"], []).append(r)
    summary = {"total_panel": total, "completed": len(results),
               "passed": len(by.get("passed", [])), "expected_block": len(by.get("expected_block", [])),
               "failed": len(by.get("FAILED", [])), "elapsed_min": round((time.time() - started) / 60, 1),
               "results": sorted(results, key=lambda r: (r["outcome"] != "FAILED", r.get("class", ""), r["pdb_id"]))}
    (root / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    lines = [f"# v1 robustness validation", "",
             f"Panel {total}, completed {len(results)}: **{summary['passed']} passed**, "
             f"{summary['expected_block']} expected-block (out of v1 scope), **{summary['failed']} FAILED**. "
             f"{summary['elapsed_min']} min.", ""]
    if by.get("FAILED"):
        lines += ["## Failures (workflow bugs to fix)", "", "| PDB | class | phase | error |", "|---|---|---|---|"]
        for r in sorted(by["FAILED"], key=lambda r: r["pdb_id"]):
            lines.append(f"| {r['pdb_id']} | {r.get('class','')} | {r.get('phase','')} | {(r.get('error','') or '')[:180].replace('|','/')} |")
        lines.append("")
    lines += ["## Expected blocks (correctly refused; out of v1 scope)", "", "| PDB | class | phase | reason |", "|---|---|---|---|"]
    for r in sorted(by.get("expected_block", []), key=lambda r: r["pdb_id"]):
        lines.append(f"| {r['pdb_id']} | {r.get('class','')} | {r.get('phase','')} | {(r.get('error','') or '')[:160].replace('|','/')} |")
    lines += ["", "## Passed", "", "| PDB | class | ligand/note | solvent | s |", "|---|---|---|---|---|"]
    for r in sorted(by.get("passed", []), key=lambda r: (r.get("class",""), r["pdb_id"])):
        lines.append(f"| {r['pdb_id']} | {r.get('class','')} | {r.get('note','')} | {r.get('solvent','')} | {r.get('elapsed_seconds','')} |")
    (root / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
