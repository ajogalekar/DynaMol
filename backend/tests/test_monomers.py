"""Bounded chain extraction preserves physical/source identity and covalent links."""
import json

import mdtraj as md
import numpy as np
import pytest

from backend import config, monomers, preparation, storage
from backend.ligands import _original_connection_blocks


@pytest.fixture(autouse=True)
def isolated(tmp_path, monkeypatch):
    for name, suffix in (("DATA_ROOT", ""), ("DATASETS_DIR", "datasets"), ("JOBS_DIR", "jobs")):
        path = tmp_path / suffix
        path.mkdir(exist_ok=True)
        monkeypatch.setattr(config, name, path)


def source(*, linked=False, ligand=False):
    top, xyz = md.Topology(), []
    for chain_id, shift in (("A", 0.), ("B", 3.)):
        chain = top.add_chain(chain_id)
        residue = top.add_residue("CYS", chain, resSeq=1)
        atoms = []
        for name, element, point in (("N", md.element.nitrogen, [0, 0, 0]), ("CA", md.element.carbon, [.14, 0, 0]), ("C", md.element.carbon, [.24, .12, 0]), ("O", md.element.oxygen, [.35, .12, 0])):
            atoms.append(top.add_atom(name, element, residue))
            xyz.append([point[0] + shift, point[1], point[2]])
        for first, second in ((0, 1), (1, 2), (2, 3)):
            top.add_bond(atoms[first], atoms[second])
        if ligand:
            residue = top.add_residue("LIG", chain, resSeq=500)
            first = top.add_atom("C1", md.element.carbon, residue)
            second = top.add_atom("C2", md.element.carbon, residue)
            xyz += [[.45 + shift, .12, 0], [.85 + shift, .12, 0]]
            top.add_bond(first, second)
    if linked:
        top.add_bond(top.atom(0), top.atom(4))
    traj = md.Trajectory(np.asarray(xyz, dtype=np.float32)[None], top, time=[4.])
    dataset = storage.save_dataset(traj, "Two chains", "Analytical fixture", "No native simulation")
    folder = storage.dataset_dir(dataset["id"])
    if linked:
        # PDB serialization may omit nonstandard bonds between standard residue
        # atoms. A supplied authoritative graph must still block this cut.
        storage.atomic_json(folder / "chemistry.json", {"authoritative_bonds": True,
                            "atoms": [{"index": atom.index} for atom in top.atoms],
                            "bonds": [{"atoms": [a.index, b.index], "order": 1} for a, b in top.bonds]})
    (folder / "source.pdb").write_text((folder / "topology.pdb").read_text())
    return dataset


def add_record(dataset, record):
    path = storage.dataset_dir(dataset["id"]) / "source.pdb"
    path.write_text(record + "\n" + path.read_text())


def extract(dataset, **options):
    return monomers.create_monomer(dataset["id"], monomers.MonomerRequest(chain_index=0, **options))


def test_extraction_preserves_source_coordinates_whole_ligand_and_resets_state():
    dataset = source(ligand=True)
    dataset.update(preparation={"ph": 7}, solvation={"padding_nm": 1})
    folder = storage.dataset_dir(dataset["id"])
    storage.atomic_json(folder / "metadata.json", dataset)
    before = {path.name: path.read_bytes() for path in folder.iterdir() if path.is_file()}
    listing = monomers.list_monomers(dataset["id"])
    assert [item["chain_id"] for item in listing["chains"]] == ["A", "B"]
    assert listing["recommended_chain_index"] == 0
    assert listing["chains"][0]["sequence_group"] == listing["chains"][1]["sequence_group"]
    result = extract(dataset)
    assert result["id"] != dataset["id"]
    assert "preparation" not in result and "solvation" not in result
    assert result["n_atoms"] == 6
    selection = result["monomer_selection"]
    assert selection["retained_atom_indices"] == list(range(6))
    assert selection["source_frame"] == 0 and selection["source_time"] == 4
    np.testing.assert_array_equal(storage.load_physical(result["id"]).xyz[0], storage.load_physical(dataset["id"]).xyz[0, :6])
    assert before == {path.name: path.read_bytes() for path in folder.iterdir() if path.is_file()}
    assert (storage.dataset_dir(result["id"]) / "atom-map.json").is_file()


def test_opt_out_excludes_unbound_associated_molecule():
    result = extract(source(ligand=True), keep_associated_molecules=False)
    assert result["n_atoms"] == 4
    assert result["monomer_selection"]["retained_associated_molecules"] == []


def test_topology_cross_chain_covalent_bond_blocks_without_publishing():
    dataset = source(linked=True)
    before = set(config.DATASETS_DIR.iterdir())
    with pytest.raises(ValueError, match="covalently connected"):
        extract(dataset)
    assert set(config.DATASETS_DIR.iterdir()) == before


def test_ssbond_blocks_cross_chain_cut_even_when_sg_atoms_are_missing():
    dataset = source()
    add_record(dataset, "SSBOND   1 CYS A    1    CYS B    1                          1555   1555  2.03")
    with pytest.raises(ValueError, match="covalently connected"):
        extract(dataset)


def test_symmetry_mate_ssbond_is_recorded_without_linking_displayed_chains():
    dataset = source()
    add_record(dataset, "SSBOND   1 CYS A    1    CYS B    1                          1555   3555  2.03")
    result = extract(dataset)
    assert result["n_atoms"] == 4
    records = result["monomer_selection"]["excluded_symmetry_connections"]
    assert len(records) == 1 and records[0]["relation"] == "different_asymmetric_units"
    assert records[0]["partners"][1]["chain"] == "B"
    assert records[0]["source_sha256"]
    assert any("symmetry mates" in message for message in result["warnings"])
    # Re-extraction must retain the explanation without restoring the ancestor
    # connection to an absent symmetry mate.
    again = extract(result)
    assert again["monomer_selection"]["excluded_symmetry_connections"] == records


def test_same_image_ssbond_remains_a_required_connection():
    dataset = source()
    add_record(dataset, "SSBOND   1 CYS A    1    CYS B    1                          3555   3555  2.03")
    with pytest.raises(ValueError, match="covalently connected"):
        extract(dataset)


def test_symmetry_record_does_not_override_an_actual_supplied_topology_bond():
    dataset = source(linked=True)
    add_record(dataset, "SSBOND   1 CYS A    1    CYS B    1                          1555   3555  2.03")
    with pytest.raises(ValueError, match="covalently connected"):
        extract(dataset)


def test_dangling_conect_on_selected_chain_is_refused():
    dataset = source()
    add_record(dataset, "CONECT    1 9999")
    before = set(config.DATASETS_DIR.iterdir())
    with pytest.raises(ValueError, match="missing residue endpoint"):
        extract(dataset)
    assert set(config.DATASETS_DIR.iterdir()) == before


def test_dangling_conect_on_excluded_chain_does_not_block():
    dataset = source()
    # PDB TER consumes a serial, so chain B's first atom is serial 6.
    text = (storage.dataset_dir(dataset["id"]) / "source.pdb").read_text()
    serial = next(int(line[6:11]) for line in text.splitlines() if line.startswith("ATOM") and line[21] == "B")
    add_record(dataset, f"CONECT{serial:5d}{9999:5d}")
    assert extract(dataset)["n_atoms"] == 4


def test_selected_sequence_gap_is_preserved_without_excluded_chain_sequence():
    fixture = config.ROOT / "docs/audit/preparation-fixtures/six_residues_known_gap.pdb"
    data = storage.save_dataset(md.load(str(fixture)), "Known gap", "Fixture", "Sequence-supported missing residue")
    original = fixture.read_bytes()
    (storage.dataset_dir(data["id"]) / "source.pdb").write_bytes(original)
    result = extract(data)
    sequence = preparation.find_sequence_source(result["id"])
    assert sequence.read_bytes() == original
    fixer, _ = preparation.current_fixer(result["id"])
    preparation.find_missing_residues_preserving_identity(fixer)
    assert fixer.missingResidues == {(0, 2): ["ILE"]}


def test_excluded_same_name_ligand_link_does_not_block_retained_ligand():
    dataset = source(ligand=True)
    add_record(dataset, "LINK         N   CYS B   1                 C1  LIG B 500                        ")
    result = extract(dataset)
    assert result["n_atoms"] == 6
    assert _original_connection_blocks(result["id"], None, {}) == set()


def test_periodic_association_is_refused_instead_of_dropping_wrapped_ligand():
    dataset = source(ligand=True)
    traj = storage.load_physical(dataset["id"])
    traj.xyz[0, 4:6, 0] = [3.95, 3.85]
    traj.unitcell_vectors = np.eye(3, dtype=np.float32)[None] * 4
    dataset = storage.save_dataset(traj, "Periodic association", "Fixture", "Boundary control")
    with pytest.raises(ValueError, match="periodic boundary"):
        extract(dataset)
