"""Ion display categories share named identity while preserving chemical compounds."""
import mdtraj as md
import numpy as np
import pytest

from backend import config, storage
from backend.ions import ION_STATES, SUPPORTED_IONS


@pytest.fixture(autouse=True)
def isolated_storage(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATASETS_DIR", tmp_path / "datasets")


def save_residue(name, elements):
    top = md.Topology()
    residue = top.add_residue(name, top.add_chain("I"), resSeq=301)
    for index, element in enumerate(elements):
        top.add_atom(f"{element}{index + 1}", md.element.get_by_symbol(element), residue)
    coordinates = np.arange(len(elements) * 3, dtype=np.float32).reshape(1, -1, 3) * .1
    return storage.save_dataset(md.Trajectory(coordinates, top), "Ion classification", "Test", "No MD")


def test_storage_ion_names_follow_the_shared_preparation_registry():
    assert storage.IONS == SUPPORTED_IONS == frozenset(ION_STATES)


@pytest.mark.parametrize("name", sorted(ION_STATES))
def test_named_monatomic_ion_is_visible_in_ion_category_after_storage_roundtrip(name):
    symbol, _charge = ION_STATES[name]
    dataset = save_residue(name, [symbol])
    assert [(atom["category"], atom["element"]) for atom in dataset["atoms"]] == [("ions", symbol)]
    assert storage.get_dataset(dataset["id"])["atoms"] == dataset["atoms"]
    physical = storage.load_physical(dataset["id"])
    assert physical.n_atoms == 1 and physical.topology.atom(0).element.symbol == symbol


@pytest.mark.parametrize("name,elements", [
    ("ZN", ["Zn", "C"]), ("CU", ["Cu", "N", "C"]),
    ("NA", ["C"]), ("ZN", ["Fe"]), ("MG", ["Ca"]),
    ("AG", ["Ag"]), ("FE", ["Fe"]), ("FE2", ["Fe"]),
    ("LIG", ["Zn"]),
])
def test_named_compounds_mismatched_elements_and_unregistered_states_remain_ligands(name, elements):
    dataset = save_residue(name, elements)
    assert [atom["category"] for atom in dataset["atoms"]] == ["ligands"] * len(elements)
