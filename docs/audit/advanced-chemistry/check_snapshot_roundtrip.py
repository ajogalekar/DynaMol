"""Freeze and reopen one research model; does not authorize production MD."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import traceback

import numpy as np
from openmm import XmlSerializer, unit

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from backend.advanced_parameters import create_snapshot, copy_snapshot, load_native_snapshot
from backend.bonded_site_validation import check_finite_forces


def run(request, output):
    request, output = Path(request).resolve(), Path(output).resolve()
    data = json.loads(request.read_text())
    output.mkdir(parents=True, exist_ok=False)
    report = {"scope": "Research artifact continuity and finite forces only; no model-accuracy, stability or app-readiness claim.",
              "request_sha256": hashlib.sha256(request.read_bytes()).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "accepted": False}
    (output / "request.json").write_bytes(request.read_bytes())
    try:
        record = create_snapshot(data["files"], output / "original", model_id=data["model_id"],
                                 kind=data["kind"], provenance=data["provenance"])
        first = load_native_snapshot(output / "original", expected_manifest_sha256=record["manifest_sha256"])
        copy_snapshot(output / "original", output / "reopened", expected_manifest_sha256=record["manifest_sha256"])
        second = load_native_snapshot(output / "reopened", expected_manifest_sha256=record["manifest_sha256"])
        original_xml = XmlSerializer.serialize(first.system)
        reopened_xml = XmlSerializer.serialize(second.system)
        xyz = second.positions.value_in_unit(unit.nanometer)
        report.update(snapshot=record, atom_count=first.system.getNumParticles(),
                      system_xml_sha256=hashlib.sha256(reopened_xml.encode()).hexdigest(),
                      system_xml_identical=original_xml == reopened_xml,
                      coordinates_identical=np.array_equal(first.positions.value_in_unit(unit.nanometer), xyz),
                      finite_force_check=check_finite_forces(second.system, xyz),
                      mechanics_accepted=second.mechanics_report["accepted"],
                      simulation_ready=second.simulation_ready)
        report["accepted"] = (report["system_xml_identical"] and report["coordinates_identical"]
                              and report["finite_force_check"]["accepted"] and report["mechanics_accepted"])
        (output / "mechanics.json").write_text(json.dumps(second.mechanics_report, indent=2) + "\n")
    except Exception as exc:
        report["error"] = f"{type(exc).__name__}: {exc}"
        (output / "traceback.txt").write_text(traceback.format_exc())
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = run(args.request, args.output)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["accepted"] else 1)
