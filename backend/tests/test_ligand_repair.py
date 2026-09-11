"""General CCD fragment repair and explicit residue decisions, without native MD."""
import numpy as np
import mdtraj as md
import pytest
from openmm import app, unit

from backend import config, ligands, storage
from backend.ligand_repair import inspect_repair, repair_ligand


class CCD:
    def __init__(self, names=None, elements=None, points=None, bonds=None):
        names = names or ["C1", "C2", "C3", "O4"]
        elements = elements or ["C", "C", "C", "O"]
        points = points or [[0, 0, 0], [1.5, 0, 0], [2, 1.4, 0], [3.4, 1.5, 0]]
        bonds = bonds or [(0, 1, "SING"), (1, 2, "SING"), (2, 3, "SING")]
        self.atoms = {"atom_id": names, "type_symbol": elements, "charge": ["0"] * len(names), "pdbx_stereo_config": ["N"] * len(names),
                      **{f"pdbx_model_Cartn_{axis}_ideal": [str(point[index]) for point in points] for index, axis in enumerate("xyz")}}
        self.bonds = {"atom_id_1": [names[a] for a, b, kind in bonds], "atom_id_2": [names[b] for a, b, kind in bonds],
                      "value_order": [kind for a, b, kind in bonds], "pdbx_stereo_config": ["N"] * len(bonds)}

    def get_mmcif_category(self, name):
        return self.atoms if name == "_chem_comp_atom." else self.bonds


def observed(block=None, present=(0, 1, 2), *, chain="A", resid="1"):
    block = block or CCD()
    top = app.Topology()
    residue = top.addResidue("LIG", top.addChain(chain), resid)
    xyz = []
    for index in present:
        top.addAtom(block.atoms["atom_id"][index], app.element.Element.getBySymbol(block.atoms["type_symbol"][index]), residue)
        xyz.append([float(block.atoms[f"pdbx_model_Cartn_{axis}_ideal"][index]) for axis in "xyz"])
    return residue, np.asarray(xyz, dtype=float), block


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    for field, suffix in (("DATA_ROOT", ""), ("DATASETS_DIR", "datasets"), ("JOBS_DIR", "jobs")):
        folder = tmp_path / suffix
        folder.mkdir(exist_ok=True)
        monkeypatch.setattr(config, field, folder)
    return tmp_path


def save_fixture(residue, xyz):
    top = md.Topology.from_openmm(residue.chain.topology)
    return storage.save_dataset(md.Trajectory(xyz[None] / 10, top), "Partial reference ligand", "Analytical fixture", "Known complete CCD graph with one unobserved terminal atom")


def test_eligible_missing_terminal_atom_is_reported_without_modeling(monkeypatch):
    residue, xyz, block = observed()
    monkeypatch.setattr(ligands.AllChem, "EmbedMolecule", lambda *args: pytest.fail("Inspection must not model coordinates"))
    result = inspect_repair(residue, xyz, block)
    assert result["can_repair"] and result["missing_heavy_atoms"] == ["O4"]
    assert "low confidence" in result["repair_reason"]


def test_constrained_completion_preserves_every_observed_coordinate_and_records_uncertainty():
    residue, xyz, block = observed()
    before = xyz.copy()
    mol, names, mapping, stereo, report = repair_ligand(residue, xyz, block, seed=2026)
    assert names == ["C1", "C2", "C3", "O4"] and mapping == [0, 1, 2, None]
    assert mol.GetNumAtoms() == 4 and mol.GetNumBonds() == 3
    np.testing.assert_array_equal(mol.GetConformer().GetPositions()[:3], before)
    np.testing.assert_array_equal(xyz, before)
    assert report["modeled_heavy_atoms"] == ["O4"] and report["confidence"] == "low"
    assert report["observed_heavy_max_displacement_angstrom"] == 0
    assert any(row.get("optimization_converged") and row["status"] == "accepted" for row in report["attempts"])


@pytest.mark.parametrize("condition", ["extra", "collinear", "disconnected", "too_few"])
def test_unsupported_or_ambiguous_repair_remains_actionable(condition):
    residue, xyz, block = observed(present=(0, 2, 3) if condition == "disconnected" else (0, 1) if condition == "too_few" else (0, 1, 2))
    if condition == "extra":
        list(residue.atoms())[0].name = "UNKNOWN"
    if condition == "collinear":
        xyz[:] = [[0, 0, 0], [1.5, 0, 0], [3, 0, 0]]
    result = inspect_repair(residue, xyz, block)
    assert not result["can_repair"] and result["repair_reason"]
    with pytest.raises(ValueError):
        repair_ligand(residue, xyz, block)


def test_missing_ring_atom_is_not_presented_as_supported_local_repair():
    points = [[1.4 * np.cos(i * np.pi / 3), 1.4 * np.sin(i * np.pi / 3), 0] for i in range(6)]
    block = CCD(names=[f"C{i}" for i in range(6)], elements=["C"] * 6, points=points, bonds=[(i, (i + 1) % 6, "AROM") for i in range(6)])
    residue, xyz, _ = observed(block, present=(0, 1, 2, 3, 4))
    result = inspect_repair(residue, xyz, block)
    assert not result["can_repair"] and "ring" in result["repair_reason"]


def test_complete_reference_is_unchanged_and_repair_is_not_required():
    residue, xyz, block = observed(present=(0, 1, 2, 3))
    assert inspect_repair(residue, xyz, block) == {"missing_heavy_atoms": [], "extra_heavy_atoms": [], "can_repair": False, "repair_reason": None}
    mol, atoms, names, _ = ligands._graph(residue, xyz, block=block)
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(), xyz)


def test_explicit_repair_and_remove_are_dry_inspection_choices(isolated, monkeypatch):
    residue, xyz, block = observed()
    dataset = save_fixture(residue, xyz)
    monkeypatch.setattr(ligands, "_ccd", lambda _: (block, {"provider": "Synthetic CCD control", "component_id": "LIG"}))
    monkeypatch.setattr(ligands.AllChem, "EmbedMolecule", lambda *args: pytest.fail("Dry inspection must not embed"))
    initial = ligands.inspect_ligands(dataset["id"])[0]
    assert initial["error"] and initial["can_repair"] and initial["can_remove"]
    repaired = ligands.inspect_ligands(dataset["id"], actions={initial["key"]: "repair"})[0]
    assert repaired["error"] is None and repaired["repair_pending"]
    removed = ligands.inspect_ligands(dataset["id"], actions={initial["key"]: "remove"})[0]
    assert removed["removed"] and removed["error"] is None
    monkeypatch.setattr(ligands, "_parameterize", lambda *args: pytest.fail("Explicitly removed residue must not be parameterized"))
    assert ligands.prepare_ligands(dataset["id"], 7, 2026, isolated / "parameters", actions={initial["key"]: "remove"}) == []
    with pytest.raises(ValueError, match="does not match"):
        ligands.inspect_ligands(dataset["id"], actions={"wrong:key": "remove"})


def test_covalent_ligand_cannot_be_removed_or_repaired_via_actions(isolated, monkeypatch):
    residue, xyz, block = observed()
    dataset = save_fixture(residue, xyz)
    monkeypatch.setattr(ligands, "_original_connection_blocks", lambda *args: {"LIG"})
    for action in ("remove", "repair"):
        result = ligands.inspect_ligands(dataset["id"], actions={"A:1::LIG": action})[0]
        assert result["error"] and not result["can_remove"] and not result["can_repair"] and not result["removed"]


def test_repair_environment_ignores_explicitly_removed_residue_and_old_hydrogens(isolated, monkeypatch):
    residue, xyz, block = observed()
    top = residue.chain.topology
    top.addAtom("H1", app.element.hydrogen, residue)
    other = top.addResidue("DRG", top.addChain("B"), "2")
    top.addAtom("C1", app.element.carbon, other)
    xyz = np.concatenate([xyz, [[.1, .1, .1], [3.4, 1.5, 0]]])
    dataset = save_fixture(residue, xyz)
    monkeypatch.setattr(ligands, "_ccd", lambda _: (block, {"provider": "Synthetic CCD control", "component_id": "LIG"}))
    import backend.ligand_repair as repair
    captured = []
    def inspect_environment(*args, **kwargs):
        captured.append(kwargs["environment_indices"])
        raise ValueError("Stop after checking retained environment")
    monkeypatch.setattr(repair, "repair_ligand", inspect_environment)
    ligands._models(dataset["id"], actions={"A:1::LIG": "repair", "B:2::DRG": "remove"}, execute_repair=True)
    assert captured == [[0, 1, 2]]
