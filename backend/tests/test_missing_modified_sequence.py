"""Absent modifications must not become standard amino acids during loop building."""
from io import StringIO
import shutil

import mdtraj as md
import pytest
from pdbfixer import PDBFixer
from pdbfixer.pdbfixer import substitutions

from backend import config, preparation, storage
from backend.preparation_worker import PreparationWorker

FIXTURE = config.ROOT / "docs/audit/preparation-fixtures/six_residues_known_gap.pdb"


@pytest.fixture
def isolated_data(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(config, "DATASETS_DIR", tmp_path / "datasets")
    monkeypatch.setattr(config, "JOBS_DIR", tmp_path / "jobs")
    config.DATASETS_DIR.mkdir()
    config.JOBS_DIR.mkdir()
    return tmp_path


def dataset_with_missing(name, temporary):
    contents = FIXTURE.read_text().replace("MET GLN ILE PHE VAL LYS", f"MET GLN {name} PHE VAL LYS")
    source = temporary / f"missing-{name}.pdb"
    source.write_text(contents)
    dataset = storage.save_dataset(md.load(str(source)), f"Missing {name}", "sequence identity regression", "Known sequence with one unobserved residue")
    shutil.copy2(source, storage.dataset_dir(dataset["id"]) / "source.pdb")
    return dataset


@pytest.mark.parametrize("name", ["SEP", "TPO", "PTR", "HYP", "MSE", "ALY", "ZZZ"])
def test_missing_modified_identity_is_retained_and_cannot_be_built(name, isolated_data):
    dataset = dataset_with_missing(name, isolated_data)
    before = dict(substitutions)
    inspection = preparation.inspect_preparation(dataset["id"])
    assert substitutions == before
    gap = inspection["missing_residues"][0]
    assert gap["residues"] == [name]
    assert gap["buildable"] is False
    assert name in gap["reason"]
    assert "no parent-residue substitution" in gap["reason"]
    assert not inspection["can_prepare"]
    for remove_ligands in (True, False):
        with pytest.raises(ValueError, match="original identities are retained"):
            preparation.submit_preparation({"dataset_id": dataset["id"], "build_missing_residues": True, "remove_heterogens": remove_ligands})
    assert not list(config.JOBS_DIR.iterdir())


def test_worker_rejects_missing_tpo_even_if_submission_validation_is_bypassed(isolated_data):
    dataset = dataset_with_missing("TPO", isolated_data)
    worker = object.__new__(PreparationWorker)
    worker.folder = isolated_data / "worker"
    worker.folder.mkdir()
    shutil.copy2(preparation.exact_input_path(dataset["id"]), worker.folder / "input.pdb")
    worker.settings = preparation._validated({"dataset_id": dataset["id"], "build_missing_residues": True})
    worker.update = lambda **kwargs: None
    with pytest.raises(ValueError, match="TPO.*original identities are retained"):
        worker.prepare()
    assert not (worker.folder / "prepared.pdb").exists()
    assert not (worker.folder / "residue-parameters").exists()


def test_standard_missing_loop_stays_buildable(isolated_data, monkeypatch):
    dataset = dataset_with_missing("ILE", isolated_data)
    inspection = preparation.inspect_preparation(dataset["id"])
    assert inspection["missing_residues"][0]["residues"] == ["ILE"]
    assert inspection["missing_residues"][0]["buildable"]
    called = []
    monkeypatch.setattr(preparation, "_submit", lambda settings, operation: called.append(settings) or {"status": "queued"})
    assert preparation.submit_preparation({"dataset_id": dataset["id"], "build_missing_residues": True})["status"] == "queued"
    assert called[0]["build_missing_residues"] is True


def test_identity_preserving_finder_leaves_upstream_substitutions_unchanged():
    contents = FIXTURE.read_text().replace("MET GLN ILE PHE VAL LYS", "MET GLN TPO PHE VAL LYS")
    native = PDBFixer(pdbfile=StringIO(contents))
    native.findMissingResidues()
    assert native.missingResidues == {(0, 2): ["THR"]}
    preserving = PDBFixer(pdbfile=StringIO(contents))
    preparation.find_missing_residues_preserving_identity(preserving)
    assert preserving.missingResidues == {(0, 2): ["TPO"]}
    native.findMissingResidues()
    assert native.missingResidues == {(0, 2): ["THR"]}
