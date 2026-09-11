"""Build a documented, exact prepared-complex slice for bounded native validation.

Run from the project root with .venv/bin/python docs/audit/complex-native-fixture.py.
This is an audit fixture builder, not a user-facing selection or preparation path.
It does not run solvation, charge assignment, minimization, or MD.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import mdtraj as md
import numpy as np
import openmm
from openmm import app, unit

from backend import storage
from backend.complex_topology import metal_environment, residue_key, subset
from backend.prepared_system import copy_ligand_parameters, load_prepared_forcefield

PARENT = "4cf5bcedf38847c1"
FIXTURE = "9ax6_cluster_native"
SELECTED = {"A", "C", "E", "F", "I", "K"}
CHAIN_MAP = [
    {"canonical_chain": "A", "source_label_asym_id": "A", "source_auth_asym_id": "A", "role": "KRAS protein, all prepared residues"},
    {"canonical_chain": "C", "source_label_asym_id": "C", "source_auth_asym_id": "D", "role": "PPIA protein, all prepared residues"},
    {"canonical_chain": "E", "source_label_asym_id": "E", "source_auth_asym_id": "A", "resid": "201", "component_id": "GNP", "role": "bound GNP, fixed charge -4"},
    {"canonical_chain": "F", "source_label_asym_id": "F", "source_auth_asym_id": "A", "resid": "202", "component_id": "MG", "role": "bound magnesium, nonbonded +2"},
    {"canonical_chain": "I", "source_label_asym_id": "I", "source_auth_asym_id": "D", "resid": "201", "canonical_residue_name": "A1A", "component_id": "A1AHB", "role": "bound RMC-6236, chosen nominal-pH state +1"},
    {"canonical_chain": "K", "source_label_asym_id": "K", "source_auth_asym_id": "A", "resids": ["307", "348"], "component_id": "HOH", "role": "both observed Mg-coordinating waters, including prepared hydrogens"},
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parent = storage.dataset_dir(PARENT)
    destination = storage.dataset_dir(FIXTURE)
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite {destination}")
    pdb = app.PDBFile(str(parent / "prepared.pdb"))
    source_preparation = json.loads((parent / "preparation.json").read_text())
    keys = {residue_key(r) for r in pdb.topology.residues() if r.chain.id in SELECTED}
    selected_indices = [a.index for a in pdb.topology.atoms() if residue_key(a.residue) in keys]
    index_map = {old: new for new, old in enumerate(selected_indices)}
    original_residues = list(pdb.topology.residues())
    modeller = subset(pdb.topology, pdb.positions, keys)
    coordinates = np.asarray(modeller.positions.value_in_unit(unit.nanometer))
    parent_xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    assert np.array_equal(coordinates, parent_xyz[selected_indices])
    radius = float(np.linalg.norm(coordinates - (coordinates.min(0) + coordinates.max(0)) / 2, axis=1).max())
    width = max(2 * radius + 1.0, 2.0)
    estimated_atoms = len(coordinates) + math.ceil(width**3 * 110)
    assert estimated_atoms < 100_000

    # Keep the parent's detailed procedure intact as a separate artifact. The
    # fixture manifest contains only correctly scoped results and explicit ancestry.
    preparation = {key: copy.deepcopy(source_preparation[key]) for key in [
        "ph", "seed", "forcefield_files", "simulation_ready", "exact_topology_file",
        "template_placement", "relaxation", "ligand_parameters", "warnings",
    ]}
    preparation.update(
        parent_dataset_id=PARENT,
        method="Exact whole-chain slice of an already prepared 9AX6 complex; no new chemistry or geometry optimization",
        summary=["Retained complete prepared KRAS A and PPIA C with GNP E, inhibitor I, Mg F, and its two coordinating waters K307/K348.",
                 "All 5,327 atom coordinates, protein/ligand hydrogens and force-field parameters are inherited unchanged from the prepared parent."],
        original_to_prepared_atom_map_file="atom-map.json",
        fixture={"purpose": "Bounded native software validation only", "parent_dataset_id": PARENT,
                 "source_dataset_id": "9fb26d98962c453b", "parent_preparation_file": "parent-preparation.json",
                 "chain_mapping": CHAIN_MAP, "parent_atom_indices_file": "fixture-atom-map.json"},
    )
    preparation["protonation_states"] = [copy.deepcopy(x) for x in source_preparation["protonation_states"] if x["chain"] in SELECTED]
    preparation["selected_variants"] = [copy.deepcopy(x) for r, x in zip(original_residues, source_preparation["selected_variants"], strict=True) if residue_key(r) in keys]
    preparation["ligand_parameters"]["ligands"] = [copy.deepcopy(x) for x in source_preparation["ligand_parameters"]["ligands"] if x["chain"] in SELECTED]
    for ligand in preparation["ligand_parameters"]["ligands"]:
        ligand["atom_index_reference"] = "Original unprepared dataset 9fb26d98962c453b for atom_indices/original_index; ligand_index remains local to the ligand. Fixture indices are in fixture-atom-map.json."
    preparation["metal_environment"] = metal_environment(modeller.topology, modeller.positions)["report"]
    stereo = copy.deepcopy(source_preparation["stereochemistry"])
    stereo["centers"] = [x for x in stereo["centers"] if x["chain"] in SELECTED]
    stereo["violations"] = [x for x in stereo["violations"] if x["chain"] in SELECTED]
    stereo["checked_centers"] = len(stereo["centers"])
    preparation["stereochemistry"] = stereo

    forcefield, _ = load_prepared_forcefield(parent, preparation)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    charges = [nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())]
    preparation["net_forcefield_charge_e"] = sum(charges)
    ligand_charges = []
    for ligand in preparation["ligand_parameters"]["ligands"]:
        residue = next(r for r in modeller.topology.residues() if r.chain.id == ligand["chain"] and r.id == ligand["resid"])
        q = sum(charges[a.index] for a in residue.atoms())
        assert abs(q - ligand["formal_charge"]) < 1e-5
        ligand_charges.append({"key": ligand["key"], "component_id": ligand["component_id"], "charge_e": q, "atoms": len(list(residue.atoms()))})

    warnings = list(dict.fromkeys(preparation["warnings"])) + [
        "Validation fixture: one observed noncovalent complex selected from the two copies in the prepared 9AX6 asymmetric unit; the full user dataset remains intact.",
        "The selected complex loses crystal contacts with the second copy. This is a computational continuity test, not a biological-assembly or convergence claim.",
    ]
    traj = md.Trajectory(coordinates[None].astype(np.float32), md.Topology.from_openmm(modeller.topology), time=[0])
    metadata = storage.save_dataset(traj, "9AX6 · single-complex native validation", "prepared complex validation fixture", preparation["summary"][0], warnings, dataset_id=FIXTURE)
    with (destination / "prepared.pdb").open("w") as stream:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, stream, keepIds=True)
    # Re-read the exact downstream topology and verify the selection survives PDB.
    exact = app.PDBFile(str(destination / "prepared.pdb"))
    assert np.array_equal(np.asarray(exact.positions.value_in_unit(unit.nanometer)), coordinates)
    assert [(residue_key(a.residue), a.name, a.element.symbol) for a in exact.topology.atoms()] == [(residue_key(a.residue), a.name, a.element.symbol) for a in modeller.topology.atoms()]
    copy_ligand_parameters(parent, destination, preparation)
    shutil.copy2(parent / "sequence-source.cif", destination / "sequence-source.cif")
    shutil.copy2(parent / "preparation.json", destination / "parent-preparation.json")
    parent_map = json.loads((parent / "atom-map.json").read_text())
    original_to_fixture = [{**item, "output_index": index_map[item["output_index"]]} for item in parent_map["input_to_output"] if item["output_index"] in index_map]
    storage.atomic_json(destination / "atom-map.json", {"input_dataset_id": "9fb26d98962c453b", "input_to_output": original_to_fixture,
                                                       "note": "Original source to fixture map composed with the prepared-parent map. All atoms, including added H, are mapped to the prepared parent in fixture-atom-map.json."})
    atoms = list(pdb.topology.atoms())
    atom_map = [{"parent_index": old, "fixture_index": new, "identity": [*residue_key(atoms[old].residue), atoms[old].name], "element": atoms[old].element.symbol} for old, new in index_map.items()]
    storage.atomic_json(destination / "fixture-atom-map.json", {"parent_dataset_id": PARENT, "atoms": atom_map})
    storage.atomic_json(destination / "preparation.json", preparation)
    metadata.update(preparation=preparation, parent_dataset_id=PARENT)
    storage.atomic_json(destination / "metadata.json", metadata)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "dataset_id": FIXTURE, "parent_dataset_id": PARENT,
        "source_dataset_id": "9fb26d98962c453b", "builder_sha256": digest(Path(__file__)),
        "selection": "Canonical chains A,C,E,F,I,K; whole proteins and all prepared waters belonging to the selected Mg environment.",
        "contact_evidence": {"source": "Original canonical physical coordinates, all heavy-atom pairs <4 Å; contact counts do not imply bond order.",
                             "inhibitor_I_to_protein_A_pairs": 56, "inhibitor_I_to_protein_C_pairs": 73,
                             "inhibitor_I_to_protein_B_pairs": 4, "inhibitor_I_to_protein_D_pairs": 0,
                             "GNP_E_to_protein_A_pairs": 154, "GNP_E_to_other_proteins_pairs": 0},
        "chain_mapping": CHAIN_MAP, "atoms": len(coordinates), "ligand_heavy_atoms": 90, "magnesium_atoms": 1,
        "coordinating_water_molecules": 2, "residues": modeller.topology.getNumResidues(),
        "coordinate_max_displacement_from_prepared_parent_angstrom": float(np.linalg.norm(coordinates - parent_xyz[selected_indices], axis=1).max()) * 10,
        "solute_extent_nm": np.ptp(coordinates, axis=0).tolist(), "bounding_radius_nm": radius,
        "cube_padding_nm": 1.0, "cube_width_nm": width, "estimate_formula": "solute_atoms + ceil(110 * max(2*radius+padding, 2*padding)^3)",
        "estimated_atoms_upper_bound": estimated_atoms, "atom_limit_unchanged": 100_000,
        "verified_system_particles": system.getNumParticles(), "net_forcefield_charge_e": sum(charges), "ligand_charges": ligand_charges,
        "metal_environment": preparation["metal_environment"],
        "source_cif_note": "The complete original mmCIF is preserved as sequence/CCD provenance. Its unselected coordinate records are not part of this fixture. Source label and author identifiers are listed explicitly; A1A is only the canonical PDB truncation of A1AHB.",
        "limitations": ["No MD or solvation was performed by this builder.", "The complete user structure remains unchanged.", "Subset omits crystal contacts with the other complex copy.", "A short downstream native trajectory verifies software execution only; no convergence, binding affinity, metal energetics or protonation correctness claim."],
        "files": {str(path.relative_to(destination)): digest(path) for path in sorted(destination.rglob("*")) if path.is_file()},
        "parent_files": {name: digest(parent / name) for name in ["prepared.pdb", "preparation.json", "sequence-source.cif"]},
    }
    storage.atomic_json(destination / "provenance.json", {"operation": "exact_prepared_complex_fixture_selection", "parent_dataset_id": PARENT, "fixture": manifest})
    storage.atomic_json(ROOT / "docs/audit/complex-native-fixture.json", manifest)
    print(json.dumps({key: manifest[key] for key in ["dataset_id", "atoms", "estimated_atoms_upper_bound", "verified_system_particles", "net_forcefield_charge_e"]}, indent=2))


if __name__ == "__main__":
    main()
