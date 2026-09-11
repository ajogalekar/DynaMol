"""Preserve a compact, provenance-linked record of an actual Playwright audit run.

Run from the repository root after `npx playwright test` in frontend. This reads
results only; it does not execute tests or invent pass/fail observations.
"""
from __future__ import annotations

import base64
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "frontend/test-results/e2e-results.json"
DESTINATION = ROOT / "docs/audit"
SELECTED_IMAGES = {
    "desktop-ribbons": "ui-audit-desktop-ribbons.png",
    "polar-hydrogens-with-context": "ui-audit-polar-hydrogens.png",
    "real-explicit-water-visible": "ui-audit-explicit-water.png",
    "preparation-complete-narrow-window": "ui-audit-preparation-complete.png",
    "phone-studio.png": "ui-audit-phone-studio.png",
    "Surface": "ui-audit-surface.png",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def specs(suites):
    for suite in suites:
        yield from suite.get("specs", [])
        yield from specs(suite.get("suites", []))


def main():
    raw = REPORT.read_bytes()
    report = json.loads(raw)
    tests = []
    artifacts = []
    for spec in specs(report["suites"]):
        for test in spec["tests"]:
            result_record = {
                "file": "frontend/tests/" + spec["file"],
                "line": spec["line"],
                "title": spec["title"],
                "outcome": test["status"],
                "expected_status": test["expectedStatus"],
                "runs": [],
            }
            for result in test["results"]:
                attempt = {key: result.get(key) for key in ("status", "duration", "startTime", "retry")}
                attempt["attachments"] = []
                if result.get("errors"):
                    attempt["errors"] = result["errors"]
                for attachment in result.get("attachments", []):
                    data = (
                        base64.b64decode(attachment["body"])
                        if "body" in attachment
                        else Path(attachment["path"]).read_bytes()
                    )
                    record = {
                        "test": spec["title"],
                        "name": attachment["name"],
                        "content_type": attachment["contentType"],
                        "bytes": len(data),
                        "sha256": sha256(data),
                    }
                    if attachment["contentType"] == "application/json":
                        record["data"] = json.loads(data)
                    if attachment["name"] in SELECTED_IMAGES:
                        path = DESTINATION / SELECTED_IMAGES[attachment["name"]]
                        path.write_bytes(data)
                        record["file"] = str(path.relative_to(ROOT))
                    attempt["attachments"].append({key: record[key] for key in ("name", "sha256")})
                    artifacts.append(record)
                result_record["runs"].append(attempt)
            tests.append(result_record)
    sources = sorted(
        set(ROOT.glob("frontend/src/**/*.ts"))
        | set(ROOT.glob("frontend/src/**/*.tsx"))
        | set(ROOT.glob("frontend/src/**/*.css"))
        | set(ROOT.glob("frontend/tests/*.spec.ts"))
        | set(ROOT.glob("backend/*.py"))
    )
    evidence = {
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "command": "DYNAMOL_BASE_URL=http://127.0.0.1:8765 DYNAMOL_AUDIT_COMPLEX_DATASET=4cf5bcedf38847c1 npx playwright test",
        "command_working_directory": "frontend",
        "platform": platform.platform(),
        "playwright_version": report["config"]["version"],
        "browser_channel": "chrome",
        "workers": report["config"]["workers"],
        "complex_dataset_id": "4cf5bcedf38847c1",
        "stats": report["stats"],
        "test_file_count": len({test["file"] for test in tests}),
        "new_audit_test_count": sum("/ui-audit-" in test["file"] for test in tests),
        "existing_regression_count": sum("/ui-audit-" not in test["file"] for test in tests),
        "report_sha256": sha256(raw),
        "full_local_report": "frontend/test-results/e2e-results.json",
        "git_base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "working_source_sha256": {str(path.relative_to(ROOT)): sha256(path.read_bytes()) for path in sources},
        "tests": tests,
        "artifacts": artifacts,
    }
    path = DESTINATION / "ui-audit-final-results.json"
    path.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps({"file": str(path), "stats": report["stats"], "tests": len(tests), "artifacts": len(artifacts)}))


if __name__ == "__main__":
    main()
