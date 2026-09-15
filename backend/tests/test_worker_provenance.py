"""Simulation exports identify the shipped release and preserve prior run provenance."""
import json
import os
import tomllib

from backend import config
from backend.worker import Worker


def test_worker_exports_current_release_and_preserves_original_on_resume(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "JOBS_DIR", tmp_path)
    folder = tmp_path / "provenance-test"
    folder.mkdir()
    status = {"worker_pid": os.getpid(), "config": {"dataset_id": "input"}}
    (folder / "status.json").write_text(json.dumps(status))
    (folder / "input.pdb").write_text("REMARK provenance-only fixture\n")

    worker = Worker(folder.name)
    version = tomllib.loads((config.ROOT / "pyproject.toml").read_text())["project"]["version"]
    assert worker.provenance["application"] == f"DynaMol {version}"
    assert worker.provenance["worker_source_sha256"]

    original = {**worker.provenance, "application": "DynaMol original release"}
    (folder / "provenance.json").write_text(json.dumps(original))
    (folder / "status.json").write_text(json.dumps({**status, "resume_requested": True}))
    assert Worker(folder.name).provenance["application"] == original["application"]
