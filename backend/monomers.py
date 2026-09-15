"""Derive one observed protein chain without changing its source or chemistry."""
from __future__ import annotations

import copy
import hashlib
import io
import json
import shutil
import uuid
from collections import defaultdict
from pathlib import Path

import mdtraj as md
import numpy as np
from fastapi import APIRouter
from pydantic import BaseModel, Field
from scipy.spatial import cKDTree

from . import storage
from .preparation import exact_input_path, find_sequence_source

router = APIRouter()
CUTOFF_NM = .5
METALS = {"Na", "K", "Mg", "Ca", "Zn", "Fe", "Mn", "Cu", "Co", "Ni", "Cd"}


class MonomerRequest(BaseModel):
    chain_index: int = Field(ge=0, strict=True)
    keep_associated_molecules: bool = True


def _load(dataset_id):
    metadata = storage.get_dataset(dataset_id)
    trajectory = storage.load_physical(dataset_id)[0]
    protein = {atom["index"] for atom in metadata["atoms"] if atom["category"] == "protein"}
    return metadata, trajectory, protein


@router.get("/api/datasets/{dataset_id}/monomers")
def list_monomers(dataset_id: str):
    _, trajectory, protein = _load(dataset_id)
    chains, groups = [], {}
    for chain in trajectory.topology.chains:
        residues = [residue for residue in chain.residues if any(atom.index in protein for atom in residue.atoms)]
        if not residues:
            continue
        sequence = tuple(residue.name for residue in residues)
        groups.setdefault(sequence, f"sequence-{len(groups) + 1}")
        chain_id = chain.chain_id or str(chain.index + 1)
        chains.append({"index": chain.index, "chain_id": chain_id, "label": f"Chain {chain_id}",
                       "n_residues": len(residues), "n_atoms": sum(atom.index in protein for residue in residues for atom in residue.atoms),
                       "sequence_group": groups[sequence], "recommended": False})
    recommended = max(chains, key=lambda item: (item["n_residues"], item["n_atoms"], -item["index"])) if chains else None
    if recommended:
        recommended["recommended"] = True
    return {"chains": chains, "recommended_chain_index": recommended["index"] if recommended else None,
            "warnings": ["Selecting one observed protein chain does not establish that it is the biological monomer. The largest observed chain is suggested; sequence groups compare observed residue names only."]}


def _ancestry(dataset_id):
    pending, visited = [dataset_id], set()
    while pending and len(visited) < 20:
        current = pending.pop(0)
        if current in visited:
            continue
        visited.add(current)
        folder = storage.dataset_dir(current)
        records = []
        for name in ("metadata.json", "provenance.json"):
            path = folder / name
            if path.is_file():
                record = json.loads(path.read_text())
                records.append(record)
                for state in (record, record.get("preparation", {}), record.get("solvation", {})):
                    if isinstance(state, dict):
                        for key in ("parent_dataset_id", "source_dataset_id"):
                            if isinstance(state.get(key), str):
                                pending.append(state[key])
        yield folder, records


def _residue_record(residue):
    return {"chain": residue.chain.id, "resid": residue.id, "insertion_code": (residue.insertionCode or "").strip(), "name": residue.name}


def _source_links(dataset_id, native, protein_residues):
    """Read explicit source links at residue level, even if an atom is missing."""
    import gemmi
    native_residues = list(native.topology.residues())
    aliases, sources, scoped = {}, [], None
    for folder, records in _ancestry(dataset_id):
        for record in records:
            for mapping in record.get("chain_id_mapping", []):
                aliases.setdefault(mapping["original"], mapping["canonical"])
            if scoped is None and "monomer_selection" in record:
                scoped = record["monomer_selection"]
        for pattern in ("source.*", "sequence-source.*", "originals/topology.*"):
            sources.extend(path for path in folder.glob(pattern) if path.suffix.lower() in {".pdb", ".cif", ".mmcif", ".pdbx"})
    source_residue_chains = {}
    def resolve(chain, resid, insertion, name):
        chain_ids = source_residue_chains.get((chain, str(resid), insertion, name), {aliases.get(chain, chain)})
        candidates = [r.index for r in native_residues if r.chain.id in chain_ids and r.id == str(resid) and (r.insertionCode or "").strip() == insertion]
        if len(candidates) > 1:
            raise ValueError("An explicit source connection has ambiguous residue identities; choose a structure with unique chain/residue identifiers before extracting one chain.")
        return candidates[0] if candidates else None
    links, seen, unresolved, symmetry_links = [], set(), set(), []
    if scoped is not None:
        # A prior extraction already scoped the original records; do not allow
        # excluded ancestor connections to re-enter subsequent derivations.
        for link in scoped.get("retained_covalent_connections", []):
            endpoints = [resolve(r["chain"], r["resid"], r.get("insertion_code", ""), r["name"]) for r in link["residues"]]
            if all(value is not None for value in endpoints):
                links.append(tuple(endpoints))
        return links, sources, unresolved, copy.deepcopy(scoped.get("excluded_symmetry_connections", []))
    for path in sources:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest in seen:
            continue
        seen.add(digest)
        source_residue_chains = {}
        cif_different_images = None
        if path.suffix.lower() in {".cif", ".mmcif", ".pdbx"}:
            from openmm import app
            from .sequence_evidence import source_cif_chain_ids, different_image_connection_ids, _text
            original = app.PDBxFile(str(path))
            label_ids = source_cif_chain_ids(path, original.topology)
            block = gemmi.cif.read_file(str(path)).sole_block()
            table = block.get_mmcif_category("_atom_site.")
            cif_different_images = different_image_connection_ids(block)
            for i, label in enumerate(table["label_asym_id"]):
                if label not in label_ids:
                    continue
                author = _text(table.get("auth_asym_id", table["label_asym_id"])[i])
                resid = _text(table.get("auth_seq_id", table.get("label_seq_id"))[i])
                code = _text(table.get("pdbx_PDB_ins_code", [None] * len(table["id"]))[i])
                name = _text(table.get("auth_comp_id", table["label_comp_id"])[i])
                original_id = label_ids[label]
                source_residue_chains.setdefault((author, resid, code, name), set()).add(aliases.get(original_id, original_id))
        structure = gemmi.read_structure(str(path))
        serial_residues = {}
        if len(structure):
            for chain in structure[0]:
                for residue in chain:
                    endpoint = resolve(chain.name, residue.seqid.num, residue.seqid.icode.strip(), residue.name)
                    for atom in residue:
                        serial_residues[atom.serial] = endpoint
            for serial, bonded in structure.conect_map.items():
                first = serial_residues.get(serial)
                for other in bonded:
                    second = serial_residues.get(other)
                    if first is not None and second is not None:
                        partners = [native_residues[first], native_residues[second]]
                        if not any(len(list(r.atoms())) == 1 and next(r.atoms()).element.symbol in METALS for r in partners):
                            links.append((first, second))
                    else:
                        unresolved.update(value for value in (first, second) if value is not None)
        for link in structure.connections:
            if link.type in {gemmi.ConnectionType.MetalC, gemmi.ConnectionType.Hydrog}:
                continue
            if (link.name in cif_different_images if cif_different_images is not None else link.asu == gemmi.Asu.Different):
                # A crystallographic symmetry mate is not the atom carrying the
                # same chain/residue name in the loaded asymmetric unit. Actual
                # bonds in the supplied topology remain independently enforced.
                symmetry_links.append({"source_sha256": digest, "source_file": path.name,
                    "connection_id": link.name, "type": link.type.name,
                    "relation": "different_asymmetric_units",
                    "partners": [{"chain": p.chain_name, "resid": str(p.res_id.seqid.num),
                        "insertion_code": p.res_id.seqid.icode.strip(), "name": p.res_id.name,
                        "atom": p.atom_name} for p in (link.partner1, link.partner2)]})
                continue
            partners = [link.partner1, link.partner2]
            endpoints = [resolve(p.chain_name, p.res_id.seqid.num, p.res_id.seqid.icode.strip(), p.res_id.name) for p in partners]
            # Monatomic metal LINK records denote coordination, not a covalent
            # requirement to drag every contacting protein chain into a subset.
            resolved = [native_residues[i] for i in endpoints if i is not None]
            if any(len(list(r.atoms())) == 1 and next(r.atoms()).element.symbol in METALS for r in resolved):
                continue
            if all(value is not None for value in endpoints):
                links.append(tuple(endpoints))
            else:
                unresolved.update(value for value in endpoints if value is not None)
    return links, sources, unresolved, symmetry_links


def _native_subset(native, trajectory, indices):
    from openmm import app, unit
    topology, atoms = app.Topology(), {}
    selected = set(indices)
    for chain in native.topology.chains():
        kept = [residue for residue in chain.residues() if any(atom.index in selected for atom in residue.atoms())]
        if not kept:
            continue
        new_chain = topology.addChain(chain.id)
        for residue in kept:
            new_residue = topology.addResidue(residue.name, new_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                if atom.index in selected:
                    atoms[atom.index] = topology.addAtom(atom.name, atom.element, new_residue, str(len(atoms) + 1))
    for a, b in trajectory.topology.bonds:
        if a.index in atoms and b.index in atoms:
            topology.addBond(atoms[a.index], atoms[b.index])
    return topology, trajectory.xyz[0, indices] * unit.nanometer


@router.post("/api/datasets/{dataset_id}/monomer", status_code=201)
def create_monomer(dataset_id: str, request: MonomerRequest):
    from openmm import app, unit
    from .modified_residues import register_topology_definitions
    metadata, trajectory, protein = _load(dataset_id)
    chains = list(trajectory.topology.chains)
    if request.chain_index >= len(chains):
        raise ValueError("The selected protein chain is not present in this structure.")
    chosen = chains[request.chain_index]
    base_atoms = {atom.index for atom in chosen.atoms if atom.index in protein}
    if not base_atoms:
        raise ValueError("Choose a chain containing protein residues.")
    register_topology_definitions()
    native = app.PDBFile(str(exact_input_path(dataset_id)))
    native_atoms, native_residues = list(native.topology.atoms()), list(native.topology.residues())
    if len(native_atoms) != trajectory.n_atoms or not np.allclose(np.asarray(native.positions.value_in_unit(unit.nanometer)), trajectory.xyz[0], atol=2e-4, rtol=0):
        raise ValueError("The exact starting topology does not match the source coordinates. Reload the structure before extracting a chain.")
    protein_residues = {native_atoms[index].residue.index for index in protein}
    selected_residues = {native_atoms[index].residue.index for index in base_atoms}
    roots = list(range(len(native_residues)))
    def root(value):
        while roots[value] != value:
            roots[value] = roots[roots[value]]
            value = roots[value]
        return value
    def connect(first, second):
        roots[root(second)] = root(first)
    links, source_paths, unresolved, symmetry_links = _source_links(dataset_id, native, protein_residues)
    residue_edges = set(tuple(sorted(pair)) for pair in links)
    for a, b in trajectory.topology.bonds:
        first, second = native_atoms[a.index], native_atoms[b.index]
        if any(len(list(atom.residue.atoms())) == 1 and atom.element.symbol in METALS for atom in (first, second)):
            continue
        residue_edges.add(tuple(sorted((first.residue.index, second.residue.index))))
    for first, second in residue_edges:
        connect(first, second)
    components = defaultdict(set)
    for residue in native_residues:
        components[root(residue.index)].add(residue.index)
    chosen_roots = {root(index) for index in selected_residues}
    other_polymer = {native_atoms[atom["index"]].residue.index for atom in metadata["atoms"] if atom["category"] in {"protein", "nucleic"}} - selected_residues
    for value in chosen_roots:
        if components[value] & other_polymer:
            raise ValueError("The selected protein chain is covalently connected to another protein or nucleic-acid chain. One-chain extraction would sever that connection; retain the linked assembly or resolve it explicitly.")
        if not request.keep_associated_molecules and components[value] - selected_residues:
            raise ValueError("A molecule is covalently attached to the selected chain. Keep associated molecules enabled; extraction will not sever that bond.")
    retained = set(selected_residues)
    if request.keep_associated_molecules:
        retained.update(index for value in chosen_roots for index in components[value])
    heavy = [index for index in base_atoms if metadata["atoms"][index]["element"] not in {"H", "D"}]
    tree = cKDTree(trajectory.xyz[0, heavy])
    periodic = trajectory.unitcell_vectors is not None
    near_atoms = set(md.compute_neighbors(trajectory, CUTOFF_NM, heavy,
                                         haystack_indices=[atom.index for atom in native_atoms if atom.element and atom.element.symbol not in {"H", "D"}],
                                         periodic=periodic)[0].tolist())
    periodic_near = {native_atoms[index].residue.index for index in near_atoms}
    distances = {}
    for residue in native_residues:
        indices = [atom.index for atom in residue.atoms() if atom.element and atom.element.symbol not in {"H", "D"}]
        distances[residue.index] = float(tree.query(trajectory.xyz[0, indices])[0].min()) if indices else float("inf")
    if request.keep_associated_molecules:
        for value, residues in components.items():
            if residues & (protein_residues | other_polymer):
                continue
            nonwater = {index for index in residues if native_residues[index].name.upper() not in storage.WATERS}
            if nonwater & periodic_near:
                retained.update(residues)
        # Keep observed first-shell water around retained monatomic metals;
        # ordinary bulk solvent is left behind before a new box is created.
        for index in list(retained):
            atoms = list(native_residues[index].atoms())
            if len(atoms) != 1 or not atoms[0].element or atoms[0].element.symbol not in METALS:
                continue
            metal = atoms[0]
            cutoff = .32 if metal.element.symbol == "Ca" else .28
            nearby = set(md.compute_neighbors(trajectory, cutoff, [metal.index], periodic=periodic)[0].tolist())
            for residue in native_residues:
                if residue.name.upper() not in storage.WATERS:
                    continue
                if any(atom.element and atom.element.symbol in {"N", "O", "S"} and atom.index in nearby for atom in residue.atoms()):
                    if periodic and not any(.10 < np.linalg.norm(trajectory.xyz[0, atom.index] - trajectory.xyz[0, metal.index]) <= cutoff for atom in residue.atoms()):
                        raise ValueError("A retained metal coordinates water across the periodic boundary. Load an unwrapped starting structure before selecting one chain; coordinates and associated molecules were not changed.")
                    retained.update(components[root(residue.index)])
    if retained & other_polymer:
        raise ValueError("An associated molecule is covalently linked to another protein or nucleic-acid chain. Extracting one chain would cut that shared molecule; retain the linked assembly.")
    if retained & unresolved:
        raise ValueError("A retained molecule has an explicit covalent source connection with a missing residue endpoint. Repair or review that connection before one-chain extraction; it will not be silently cut.")
    if periodic and any(index in periodic_near and distances[index] > CUTOFF_NM + 1e-6 for index in retained - selected_residues):
        raise ValueError("An associated molecule contacts the selected chain across the periodic boundary. Load an unwrapped starting structure before selecting one chain; coordinates and associated molecules were not changed.")
    indices = [atom.index for atom in native_atoms if atom.residue.index in retained]
    retained_atom_set = set(indices)
    bonds = [[a.index, b.index] for a, b in trajectory.topology.bonds if a.index in retained_atom_set and b.index in retained_atom_set]
    if periodic and bonds:
        pairs = np.asarray(bonds)
        direct = np.linalg.norm(trajectory.xyz[0, pairs[:, 0]] - trajectory.xyz[0, pairs[:, 1]], axis=1)
        minimum_image = md.compute_distances(trajectory, bonds, periodic=True)[0]
        if np.any(np.abs(direct - minimum_image) > 1e-4):
            raise ValueError("A retained molecule crosses the periodic boundary. Load an unwrapped starting structure before removing the old box; one-chain extraction will not split its coordinates.")
    index_map = {index: new for new, index in enumerate(indices)}
    derived = trajectory.atom_slice(indices)
    derived.unitcell_vectors = None
    derived.time = np.asarray([0.])
    chain_id = chosen.chain_id or str(chosen.index + 1)
    notes = ["One observed protein chain was selected; this is not a biological-assembly or monomer assignment.",
             "Associated molecules are retained whole when covalently attached or when any heavy atom is within 5 Å of the selected chain in the first source frame. This proximity rule does not establish biological binding.",
             "Ordinary solvent and the previous periodic box are excluded. Observed coordinating waters of retained monatomic metals are retained. Prepare the derived structure and create its solvent box again."]
    if symmetry_links:
        notes.append(f"The source records {len(symmetry_links)} covalent connection(s) to crystallographic symmetry mates outside the loaded asymmetric unit. These records were not treated as bonds between the displayed chains; review the crystal assembly if those contacts matter to your model.")
    excluded_chains = [item for item in list_monomers(dataset_id)["chains"] if item["index"] != chosen.index]
    associated = [_residue_record(native_residues[index]) for index in sorted(retained - selected_residues)]
    excluded = [_residue_record(residue) for residue in native_residues if residue.index not in retained and residue.index not in protein_residues]
    other_heavy = [atom["index"] for atom in metadata["atoms"] if atom["category"] == "protein" and atom["index"] not in base_atoms and atom["element"] not in {"H", "D"}]
    if other_heavy and associated:
        shared_atoms = set(md.compute_neighbors(trajectory, CUTOFF_NM, other_heavy, periodic=periodic)[0].tolist())
        shared = [index for index in retained - selected_residues if any(atom.index in shared_atoms for atom in native_residues[index].atoms())]
        if shared:
            notes.append("Some retained associated molecules also contact an excluded protein chain. Review those interface molecules before preparation.")
    selection = {"parent_dataset_id": dataset_id, "chain_index": chosen.index, "chain_id": chain_id, "chain_label": f"Chain {chain_id}",
                 "retained_atom_indices": indices, "excluded_protein_chains": excluded_chains, "retained_associated_molecules": associated,
                 "excluded_molecules": excluded, "association_cutoff_angstrom": 5, "keep_associated_molecules": request.keep_associated_molecules,
                 "source_frame": 0, "source_time": float(trajectory.time[0]), "source_time_unit": metadata.get("time_unit", "ps"),
                 "requires_preparation": True, "notes": notes,
                 "excluded_symmetry_connections": symmetry_links,
                 "retained_covalent_connections": [{"residues": [_residue_record(native_residues[first]), _residue_record(native_residues[second])]} for first, second in sorted(residue_edges) if first != second and first in retained and second in retained and not ({first, second} <= protein_residues)]}
    new_id = uuid.uuid4().hex[:16]
    folder = storage.dataset_dir(new_id)
    folder.mkdir()
    provenance = {"operation": "one observed protein chain", "source_dataset_id": dataset_id, "monomer_selection": selection,
                  "input_sha256": hashlib.sha256(exact_input_path(dataset_id).read_bytes()).hexdigest()}
    try:
        storage.atomic_json(folder / "provenance.json", provenance)
        sequence = find_sequence_source(dataset_id)
        if sequence:
            shutil.copy2(sequence, folder / ("sequence-source" + sequence.suffix.lower()))
        # Retain original source evidence, never source force-field/charge files.
        originals = folder / "original-source"
        originals.mkdir()
        seen = set()
        for path in source_paths:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest not in seen:
                shutil.copy2(path, originals / (digest[:12] + path.suffix.lower()))
                seen.add(digest)
        native_topology, native_positions = _native_subset(native, trajectory, indices)
        exact = io.StringIO()
        app.PDBFile.writeFile(native_topology, native_positions, exact, keepIds=True)
        result = storage.save_dataset(derived, f"{metadata['name']} · chain {chain_id}"[:100], "Selected protein chain", f"Chain {chain_id} derived from {metadata['name']}; source frame 0 retained. Prepare again before simulation.", warnings=notes, dataset_id=new_id, provenance=provenance, exact_pdb=exact.getvalue())
        result.update(parent_dataset_id=dataset_id, monomer_selection=selection)
        result["chain_id_mapping"] = copy.deepcopy(metadata.get("chain_id_mapping", []))
        # Category/chemical identities are inherited atom-by-atom, independently
        # of PDB residue-name heuristics or whether an MD template exists.
        for target, index in zip(result["atoms"], indices):
            target["category"] = metadata["atoms"][index]["category"]
        chemistry_path = storage.dataset_dir(dataset_id) / "chemistry.json"
        if chemistry_path.is_file():
            original = json.loads(chemistry_path.read_text())
            if original.get("authoritative_bonds"):
                chemical = {"schema_version": 1, "authoritative_bonds": True, "input_format": original.get("input_format"),
                            "atoms": [{**original["atoms"][index], "index": new} for new, index in enumerate(indices)],
                            "bonds": [{**bond, "atoms": [index_map[index] for index in bond["atoms"]]} for bond in original["bonds"] if all(index in index_map for index in bond["atoms"])],
                            "source_dataset_id": dataset_id, "note": "Subset of original chemical graph; no source partial charges are assigned as MD parameters."}
                storage.atomic_json(folder / "chemistry.json", chemical)
                result["bonds"] = [bond["atoms"] for bond in chemical["bonds"]]
                result["bond_orders"] = chemical["bonds"]
                result["chemistry"] = {key: value for key, value in chemical.items() if key not in {"atoms", "bonds"}}
        storage.atomic_json(folder / "metadata.json", result)
        storage.atomic_json(folder / "atom-map.json", {"parent_dataset_id": dataset_id, "source_frame": 0, "atoms": [{"source_index": index, "output_index": new} for new, index in enumerate(indices)]})
        return result
    except Exception:
        shutil.rmtree(folder, ignore_errors=True)
        raise
