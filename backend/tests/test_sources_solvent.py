import hashlib
import json
import shutil
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest

from backend import config, sources, storage
from backend.solvent import solvate_dataset


@pytest.fixture(autouse=True)
def isolated_sources(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(config, "DATASETS_DIR", tmp_path / "datasets")
    monkeypatch.setattr(config, "JOBS_DIR", tmp_path / "jobs")
    config.DATASETS_DIR.mkdir()
    config.JOBS_DIR.mkdir()


ETHANOL_MOL2 = """@<TRIPOS>MOLECULE
ethanol
3 2 1 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.000 0.000 0.000 C.3 1 LIG1 -0.1
2 C2 1.500 0.000 0.000 C.3 1 LIG1 0.3
3 O1 2.000 1.200 0.000 O.3 1 LIG1 -0.2
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
@<TRIPOS>SUBSTRUCTURE
1 LIG1 1 GROUP 0 A
"""


def test_smiles_reproducible_coordinates_chemistry_and_source():
    first = sources.create_smiles("C[C@H](O)C(=O)[O-]", seed=47)
    second = sources.create_smiles("C[C@H](O)C(=O)[O-]", seed=47)
    traj = storage.load_physical(first["id"])
    assert np.array_equal(traj.xyz, storage.load_physical(second["id"]).xyz)
    chemistry = json.loads((storage.dataset_dir(first["id"]) / "chemistry.json").read_text())
    assert chemistry["total_formal_charge"] == -1
    assert "@" in chemistry["canonical_isomeric_smiles"]
    assert chemistry["unassigned_stereocenters"] == 0
    assert chemistry["conformer"]["seed"] == 47
    assert any(bond["order"] == 2 for bond in chemistry["bonds"])
    assert traj.n_atoms == first["n_atoms"]
    assert sum(atom["nonpolar_hydrogen"] for atom in first["atoms"]) == 4
    assert all(atom["category"] == "ligands" for atom in first["atoms"])
    source = storage.dataset_dir(first["id"]) / "source.smiles"
    assert source.read_text() == "C[C@H](O)C(=O)[O-]"
    provenance = json.loads((source.parent / "provenance.json").read_text())
    assert provenance["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()


def test_aromatic_graph_survives_pdb_roundtrip_in_measurement_topology():
    result = sources.create_smiles("c1ccccc1", seed=47)
    loaded = storage.load_physical(result["id"])
    graph = result["bond_orders"]
    assert sorted(tuple(sorted(b["atoms"])) for b in graph) == sorted(tuple(sorted((a.index, b.index))) for a, b in loaded.topology.bonds)
    assert sum(bond["aromatic"] for bond in graph) == 6
    assert sum(bond.type == md.core.topology.Aromatic for bond in loaded.topology.bonds) == 6


@pytest.mark.parametrize("smiles", ["not a molecule", "", "C.C", "CCO\nCCC"])
def test_invalid_smiles_has_no_partial_dataset(smiles):
    with pytest.raises(ValueError):
        sources.create_smiles(smiles)
    assert list(config.DATASETS_DIR.iterdir()) == []


def test_mol2_preserves_graph_input_charge_residue_and_units(tmp_path):
    path = tmp_path / "ethanol.mol2"
    path.write_text(ETHANOL_MOL2)
    result = sources.import_structure(path)
    assert result["bonds"] == [[0, 1], [1, 2]]
    assert [a["partial_charge"] for a in result["atoms"]] == [-0.1, 0.3, -0.2]
    assert result["atoms"][0]["input_residue"] == "LIG1"
    assert result["atoms"][0]["category"] == "ligands"
    assert storage.load_physical(result["id"]).xyz[0, 1, 0] == pytest.approx(0.15)
    assert (storage.dataset_dir(result["id"]) / "source.mol2").read_text() == ETHANOL_MOL2
    assert sorted((a.index, b.index) for a, b in storage.load_physical(result["id"]).topology.bonds) == [(0, 1), (1, 2)]


def test_mol2_does_not_silently_call_a_malformed_protein_a_ligand(tmp_path):
    path = tmp_path / "protein.mol2"
    path.write_text(ETHANOL_MOL2.replace("LIG1", "ALA1"))
    with pytest.raises(ValueError, match="backbone atom identities"):
        sources.import_structure(path)
    assert not list(config.DATASETS_DIR.iterdir())


def test_unsupported_mol2_bond_type_rejected(tmp_path):
    path = tmp_path / "unknown.mol2"
    path.write_text(ETHANOL_MOL2.replace("2 2 3 1", "2 2 3 un"))
    with pytest.raises(ValueError, match="bond type"):
        sources.import_structure(path)


@pytest.mark.parametrize("data", [ETHANOL_MOL2.replace("2 C2", "2 C1"), ETHANOL_MOL2.replace("C.3 1 LIG1 0.3", "C.3 2 LIG2 0.3")])
def test_mol2_ambiguous_atom_identity_leaves_no_dataset(tmp_path, data):
    path = tmp_path / "ambiguous.mol2"
    path.write_text(data)
    with pytest.raises(ValueError, match="unique|contiguous"):
        sources.import_structure(path)
    assert not list(config.DATASETS_DIR.iterdir())


def test_pdb_retains_seqres_original_bytes(tmp_path):
    fixture = config.ROOT / "examples/ubiquitin-start/topology.pdb"
    data = b"SEQRES   1 A   76  MET GLN ILE PHE VAL LYS THR LEU THR GLY LYS THR ILE\n" + fixture.read_bytes()
    source = tmp_path / "protein.pdb"
    source.write_bytes(data)
    result = sources.import_structure(source)
    folder = storage.dataset_dir(result["id"])
    assert (folder / "source.pdb").read_bytes() == data
    assert "SEQRES" not in (folder / "topology.pdb").read_text()
    assert result["n_atoms"] == 602


def test_fetch_restricts_identifiers_before_network(monkeypatch):
    def unexpected(*args):
        raise AssertionError("network should not be invoked")
    monkeypatch.setattr(sources, "_download", unexpected)
    for provider, identifier in [("pdb", "../../1ubq"), ("pdb", "http://localhost/"), ("pubchem", "http://127.0.0.1"), ("pubchem", "a\nb"), ("unknown", "1UBQ")]:
        with pytest.raises(ValueError):
            sources.fetch_structure(provider, identifier)


def test_download_rejects_untrusted_hosts_and_redirects():
    for url in ["http://files.rcsb.org/download/1ubq.cif", "https://files.rcsb.org.evil.test/x", "https://127.0.0.1/x", "https://files.rcsb.org@evil.test/x"]:
        with pytest.raises(ValueError, match="fixed"):
            sources._download(url)
    with pytest.raises(ValueError, match="redirect"):
        sources._NoRedirect().redirect_request(None, None, 302, "", {}, "http://127.0.0.1")


def test_pubchem_ambiguous_name_requires_cid(monkeypatch):
    monkeypatch.setattr(sources, "_download", lambda _: b'{"IdentifierList":{"CID":[1,2]}}')
    with pytest.raises(ValueError, match="multiple"):
        sources.fetch_structure("pubchem", "ambiguous compound")


def test_solvation_requires_preparation_even_if_crystal_has_a_cell():
    result = sources.import_structure(config.ROOT / "examples/ubiquitin-start/topology.pdb")
    assert result["has_unitcell"]
    with pytest.raises(ValueError, match="Prepare the protein first"):
        solvate_dataset(result["id"])


def test_real_tip3p_preview_reuses_exact_parent_and_rejects_oversize(monkeypatch):
    from openmm import app, unit
    source = config.ROOT / "examples/demo/topology.pdb"
    parent = sources.import_structure(source, name="Prepared ubiquitin fixture")
    parent["preparation"] = {"ph": 7, "method": "Previously simulated demonstration topology", "summary": []}
    folder = storage.dataset_dir(parent["id"])
    shutil.copyfile(source, folder / "prepared.pdb")
    storage.atomic_json(folder / "metadata.json", parent)
    saved_cap = config.MAX_ATOMS
    monkeypatch.setattr(config, "MAX_ATOMS", 1300)
    with pytest.raises(ValueError, match="exceed"):
        solvate_dataset(parent["id"])
    monkeypatch.setattr(config, "MAX_ATOMS", saved_cap)
    preview = solvate_dataset(parent["id"], padding_nm=1, seed=17)
    assert preview["has_unitcell"]
    assert preview["solvation"]["water_atoms"] > 3000
    assert preview["solvation"]["water_atoms"] % 3 == 0
    assert preview["solvation"]["parent_dataset_id"] == parent["id"]
    assert preview["preparation"] == parent["preparation"]
    assert preview["n_atoms"] > parent["n_atoms"]
    assert solvate_dataset(preview["id"], padding_nm=1, seed=17)["id"] == preview["id"]
    coordinates = storage.load_physical(preview["id"])
    assert np.allclose(coordinates.xyz[0, :parent["n_atoms"]], storage.load_physical(parent["id"]).xyz[0], atol=1e-6)
    exact = app.PDBFile(str(storage.dataset_dir(preview["id"]) / "prepared.pdb"))
    assert exact.topology.getNumAtoms() == preview["n_atoms"]
    assert np.allclose(exact.positions.value_in_unit(unit.nanometer), coordinates.xyz[0], atol=0.000051)
    assert exact.topology.getPeriodicBoxVectors() is not None
    water_indices = [a["index"] for a in preview["atoms"] if a["category"] == "water"]
    assert len(water_indices) == preview["solvation"]["water_atoms"]
    assert sum(a["category"] == "water" for a in storage.get_dataset(parent["id"])["atoms"]) == 0
