"""Regressions found by exercising the import dialog with real molecular files."""
import hashlib
import json
from pathlib import Path

import mdtraj as md
import numpy as np
from fastapi.testclient import TestClient

from backend import config, storage
from backend.main import app


def test_mmcif_upload_canonicalizes_serials_without_changing_atoms_or_source(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(config, "DATASETS_DIR", tmp_path / "datasets")
    monkeypatch.setattr(config, "JOBS_DIR", tmp_path / "jobs")
    config.DATASETS_DIR.mkdir()
    config.JOBS_DIR.mkdir()
    fixture = Path(__file__).resolve().parents[2] / "docs/audit/ui-audit-fixtures/short.cif"
    original = fixture.read_bytes()
    expected = md.load(str(fixture))
    assert isinstance(expected.topology.atom(0).serial, str)
    with TestClient(app) as client:
        response = client.post(
            "/api/datasets/upload",
            files={"topology": ("audit.cif", original, "chemical/x-mmcif")},
            data={"stride": "2", "frame_interval_ps": "0.25"},
        )
    assert response.status_code == 200, response.text
    result = response.json()
    assert result["n_atoms"] == expected.n_atoms
    folder = storage.dataset_dir(result["id"])
    assert (folder / "originals/topology.cif").read_bytes() == original
    provenance = json.loads((folder / "provenance.json").read_text())
    assert provenance["sha256"]["topology"] == hashlib.sha256(original).hexdigest()
    actual = storage.load_physical(result["id"])
    assert np.array_equal(actual.xyz, expected.xyz)
    identity = lambda top: [(a.name, a.residue.name, a.element.symbol) for a in top.atoms]
    assert identity(actual.topology) == identity(expected.topology)
    assert [a.serial for a in actual.topology.atoms] == list(range(1, expected.n_atoms + 1))
