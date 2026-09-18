"""Loop eligibility is per region; repeated monomers share only a work budget."""
import shutil

import mdtraj as md
import pytest

from backend import config, jobs, preparation, storage
from backend.models import SimulationConfig
from backend.readiness import ReadinessRequest, readiness
from backend.preparation_worker import PreparationWorker


@pytest.fixture
def data(tmp_path, monkeypatch):
    # Eligibility controls isolate installation; native runtime has separate checks.
    monkeypatch.setattr("backend.loop_modeling.loop_runtime_status", lambda: {"available": True})
    for name, value in (("DATA_ROOT", tmp_path), ("DATASETS_DIR", tmp_path / "datasets"),
                        ("JOBS_DIR", tmp_path / "jobs")):
        monkeypatch.setattr(config, name, value)
        value.mkdir(exist_ok=True)
    return tmp_path


def imported_gap(data, count, chains="A", spatially_closed=False):
    filename = "six_residues_intact.pdb" if spatially_closed else "six_residues_known_gap.pdb"
    original = (config.ROOT / "docs/audit/preparation-fixtures" / filename).read_text()
    sequence = ["MET", "GLN"] + ["ALA"] * count + (["ILE"] if spatially_closed else []) + ["PHE", "VAL", "LYS"]
    lines = []
    for chain in chains:
        for offset in range(0, len(sequence), 13):
            lines.append(f"SEQRES {offset // 13 + 1:3d} {chain} {len(sequence):4d}  " + " ".join(sequence[offset:offset + 13]))
    for chain in chains:
        for line in original.splitlines():
            if line.startswith("ATOM"):
                resid = int(line[22:26])
                if resid >= (3 if spatially_closed else 4):
                    resid += count if spatially_closed else count - 1
                lines.append(line[:21] + chain + f"{resid:4d}" + line[26:])
        lines.append("TER")
    source = data / f"gap-{count}-{chains}.pdb"
    source.write_text("\n".join(lines) + "\nEND\n")
    dataset = storage.save_dataset(md.load(str(source)), "Declared loop", "Eligibility fixture",
                                   "Synthetic sequence gap on observed 1UBQ-derived anchors; no modeled accuracy claim")
    shutil.copy2(source, storage.dataset_dir(dataset["id"]) / "source.pdb")
    return dataset


@pytest.mark.parametrize("count,eligible,modeling", [(6, True, "short"), (7, True, "extended"),
                                                     (12, True, "extended"), (13, False, "unsupported")])
def test_inspection_and_submission_agree_at_per_loop_boundary(data, count, eligible, modeling):
    dataset = imported_gap(data, count)
    inspection = preparation.inspect_preparation(dataset["id"])
    gap, = inspection["missing_residues"]
    assert gap["count"] == count
    assert gap["buildable"] == eligible
    assert gap["modeling"] == modeling
    assert inspection["loop_policy"]["max_gap_residues"] == 12
    settings = {"dataset_id": dataset["id"], "build_missing_residues": True}
    if eligible:
        preparation.validate_preparation(settings)
    else:
        with pytest.raises(ValueError, match="can't be repaired here"):
            preparation.validate_preparation(settings)
    assert not list(config.JOBS_DIR.iterdir())


def test_four_monomers_with_twelve_residue_loops_can_prepare(data):
    dataset = imported_gap(data, 12, "ABCD")
    _, inspection = preparation.validate_preparation({"dataset_id": dataset["id"], "build_missing_residues": True})
    assert [gap["count"] for gap in inspection["missing_residues"]] == [12] * 4
    with pytest.raises(ValueError, match="Build supported missing loops"):
        preparation.validate_preparation({"dataset_id": dataset["id"]})


def test_missing_loop_runtime_is_reported_before_queuing(data, monkeypatch):
    dataset = imported_gap(data, 12)
    monkeypatch.setattr("backend.loop_modeling.loop_runtime_status",
                        lambda: {"available": False, "error": "Loop runtime is missing."})
    result = readiness(dataset["id"], ReadinessRequest(mode="preparation", settings={"build_missing_residues": True}))
    assert not result["ready"]
    assert any("Loop runtime is missing" in error for error in result["blockers"])
    with pytest.raises(ValueError, match="Loop runtime is missing"):
        preparation.submit_preparation({"dataset_id": dataset["id"], "build_missing_residues": True})
    assert not list(config.JOBS_DIR.iterdir())


def test_total_work_budget_is_distinct_from_gap_length():
    loops = [{"chain": str(i), "residues": ["ALA"] * 12} for i in range(9)]
    with pytest.raises(ValueError, match="loop-building budget"):
        preparation.validate_loop_selection(loops, True)


def test_raw_simulation_and_readiness_reject_sequence_gaps_even_with_close_endpoints(data):
    dataset = imported_gap(data, 12, spatially_closed=True)
    fixer, _ = preparation.current_fixer(dataset["id"])
    assert not any(gap["structural_break"] for gap in preparation.backbone_gaps(fixer.topology, fixer.positions))
    result = readiness(dataset["id"], ReadinessRequest(mode="simulation", settings={"solvent": "implicit"}))
    assert not result["ready"]
    assert any("source sequence contains unresolved internal residues" in error for error in result["blockers"])
    with pytest.raises(ValueError, match="source sequence contains unresolved internal residues"):
        jobs.submit_job(SimulationConfig(dataset_id=dataset["id"], solvent="implicit"))
    assert not list(config.JOBS_DIR.iterdir())


def test_worker_rechecks_length_when_submission_is_bypassed(data):
    dataset = imported_gap(data, 13)
    worker = object.__new__(PreparationWorker)
    worker.folder = data / "bypass-worker"
    worker.folder.mkdir()
    shutil.copy2(preparation.exact_input_path(dataset["id"]), worker.folder / "input.pdb")
    worker.settings = preparation._validated({"dataset_id": dataset["id"], "build_missing_residues": True})
    worker.update = lambda **kwargs: None
    with pytest.raises(ValueError, match="can't be repaired here"):
        worker.prepare()
    assert not (worker.folder / "prepared.pdb").exists()
