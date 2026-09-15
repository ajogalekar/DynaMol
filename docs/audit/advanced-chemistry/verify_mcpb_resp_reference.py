"""Compare native RESP using array-written ESP against published 1OKL charges.

Requires preserved native step-1/3b input artifacts in native-regenerated; output
directory is new on each invocation and every failed attempt is retained.
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np

from mcpb_array_adapter import run_native_resp


def published_mol2_charges(path):
    block = path.read_text().split("@<TRIPOS>ATOM", 1)[1].split("@<TRIPOS>", 1)[0]
    rows = [line.split() for line in block.splitlines() if line.strip()]
    return {row[1]: float(row[8]) for row in rows}


def main():
    base = Path(__file__).resolve().parent
    ref = base / "reference-mcpb-1okl"
    native = ref / "native-regenerated"
    atom_ids = (native / "1OKL_large.fingerprint").read_text().splitlines()
    root = base.parents[2]
    # Native MCPB's own input generation with ff14SB and 3b constraints is used;
    # ordinary charge/protonation decisions are not reverse-engineered here.
    output = Path(tempfile.mkdtemp(prefix="resp-parity-", dir=ref)) / "fit"
    result = run_native_resp(executable=root / ".tools/ambertools/bin/resp", workdir=output,
                             stage1_input=native / "resp1.in", stage2_input=native / "resp2.in",
                             esp_file=ref / "array-adapter.esp", atom_ids=atom_ids, total_charge=1)
    charge_map = dict(zip(atom_ids, result["charges"], strict=True))
    residue_files = {"90-HID": "HD1.mol2", "92-HID": "HD2.mol2", "115-HIE": "HE1.mol2",
                     "256-ZN": "ZN1.mol2", "257-MNS": "MS1.mol2"}
    comparisons = []
    for residue, filename in residue_files.items():
        for atom, published in published_mol2_charges(ref / filename).items():
            key = f"{residue}-{atom}"
            actual = charge_map[key]
            comparisons.append({"atom_id": key, "published": published, "actual": actual,
                                "absolute_difference": abs(actual-published),
                                "pass": actual == published})
    independent_direct = np.asarray([float(s) for s in (native / "direct-resp2.chg").read_text().split()])
    direct_equal = np.array_equal(independent_direct, np.asarray(result["charges"]))
    report = {"status": "passed" if all(row["pass"] for row in comparisons) and direct_equal else "mismatch",
              "scope": "Published-reference RESP parity only; not new-site model validation",
              "method": "Native MCPB 3b input constraints, ff14SB; direct native RESP over array-written ESP",
              "published_atom_count": len(comparisons), "published_comparisons": comparisons,
              "same_charges_as_direct_native_original_esp": direct_equal,
              "max_published_absolute_difference": max(row["absolute_difference"] for row in comparisons),
              "fit": result, "output_directory": str(output),
              "preserved_failure": "native-regenerated/step3b.log: transient dyld shared-cache failure during shell-invoked second RESP; direct subprocess stages succeeded"}
    (ref / "resp-parity-report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k: report[k] for k in ("status", "published_atom_count",
                     "max_published_absolute_difference", "same_charges_as_direct_native_original_esp")}, indent=2))
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
