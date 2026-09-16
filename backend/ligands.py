"""Chemistry-preserving preparation of noncovalent organic ligands.

Coordinates come from the current exact input, chemistry from a retained chemical
sidecar or the full CCD component ID in the original mmCIF.  Neither guessed PDB
connectivity nor a truncated component name is an acceptable chemical reference.
"""
from __future__ import annotations

from dataclasses import dataclass
import copy
import datetime as dt
import uuid
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import urllib.request

import numpy as np
from openmm import app, unit
from rdkit import Chem
from rdkit.Chem import AllChem, rdCIPLabeler

from . import config, storage
from .sources import STANDARD_RESIDUES
from .ions import SUPPORTED_IONS, inspect_ions

ORGANIC_ELEMENTS = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
MAX_LIGAND_HEAVY_ATOMS = 200


@dataclass
class PreparedLigand:
    topology: app.Topology
    positions: object
    ffxml: Path
    provenance: dict
    original_residue_key: tuple
    atom_map: list[dict]


@dataclass
class LigandModel:
    mol: object
    residue: object
    atom_indices: list[int | None]
    names: list[str]
    description: dict


def residue_key(residue):
    return (residue.chain.id, residue.id, (residue.insertionCode or "").strip(), residue.name)


def _key(residue):
    return ":".join(residue_key(residue))


def ambertools_home():
    explicit = os.environ.get("DYNAMOL_AMBERTOOLS")
    candidates = [Path(explicit)] if explicit else [config.ROOT / ".tools" / "ambertools"]
    for folder in candidates:
        if all((folder / "bin" / name).is_file() and os.access(folder / "bin" / name, os.X_OK)
               for name in ("antechamber", "parmchk2", "tleap", "sqm")):
            return folder.resolve()
    raise ValueError("Ligand parameters need AmberTools. Run scripts/install_ligand_tools.sh, or set DYNAMOL_AMBERTOOLS to a complete AmberTools installation.")


def ligand_runtime_status():
    try:
        folder = ambertools_home()
        for package in ("parmed", "gemmi", "dimorphite-dl"):
            importlib.metadata.version(package)
        return {"available": True, "path": str(folder), "method": "GAFF2 / AM1-BCC"}
    except (ValueError, importlib.metadata.PackageNotFoundError) as exc:
        return {"available": False, "message": str(exc), "method": "GAFF2 / AM1-BCC"}


def _source_cif(dataset_id):
    """Find preserved source records without silently switching to a remote model."""
    pending, visited = [dataset_id], set()
    while pending and len(visited) < 20:
        current = pending.pop(0)
        if current in visited:
            continue
        visited.add(current)
        folder = storage.dataset_dir(current)
        for path in [folder / "source.cif", folder / "source.mmcif", folder / "source.pdbx",
                     folder / "sequence-source.cif", folder / "sequence-source.mmcif",
                     *sorted((folder / "originals").glob("*.cif"))]:
            if path.is_file():
                return path
        for filename in ("metadata.json", "provenance.json"):
            path = folder / filename
            if not path.is_file():
                continue
            record = json.loads(path.read_text())
            if not isinstance(record, dict):
                continue
            for state in (record, record.get("preparation", {}), record.get("solvation", {}),
                          record.get("preparation_state", {})):
                if not isinstance(state, dict):
                    continue
                parent = state.get("parent_dataset_id") or state.get("source_dataset_id")
                if isinstance(parent, str):
                    pending.append(parent)
    return None


def _component_ids(dataset_id, topology):
    """Resolve full 5-character CCD identifiers lost by the PDB viewer format."""
    path = _source_cif(dataset_id)
    if not path:
        return {}, None
    import gemmi
    block = gemmi.cif.read_file(str(path)).sole_block()
    table = block.get_mmcif_category("_atom_site.")
    if not table:
        return {}, path
    groups = {}
    count = len(table.get("id", []))
    for i in range(count):
        if table.get("pdbx_PDB_model_num", ["1"] * count)[i] != "1":
            continue
        if table["type_symbol"][i].upper() in {"H", "D"}:
            continue
        label = table["label_asym_id"][i]
        author = (table.get("auth_asym_id") or table["label_asym_id"])[i]
        resid = (table.get("auth_seq_id") or table["label_seq_id"])[i]
        insertion = table.get("pdbx_PDB_ins_code", [None] * count)[i] or ""
        if insertion in {".", "?"}:
            insertion = ""
        group = groups.setdefault((label, author, resid, insertion), {"id": table["label_comp_id"][i], "names": set()})
        group["names"].add((table.get("auth_atom_id") or table["label_atom_id"])[i])
    mapping = {}
    for residue in topology.residues():
        if residue.name in STANDARD_RESIDUES or residue.name in storage.WATERS or residue.name in SUPPORTED_IONS:
            continue
        names = {atom.name for atom in residue.atoms() if atom.element != app.element.hydrogen}
        candidates = [(key, group) for key, group in groups.items()
                      if residue.id == key[2] and ((residue.insertionCode or "").strip()) == key[3]
                      and names == group["names"]]
        preferred = [group for key, group in candidates if key[0] == residue.chain.id]
        if not preferred:
            preferred = [group for key, group in candidates if key[1] == residue.chain.id]
        if not preferred:
            preferred = [group for _, group in candidates]
        ids = {group["id"] for group in preferred}
        if len(ids) == 1:
            mapping[_key(residue)] = ids.pop()
        elif len(ids) > 1:
            raise ValueError(f"{_key(residue)} has ambiguous full component identity in the original mmCIF.")
    return mapping, path


def _ccd(component_id):
    if not re.fullmatch(r"[A-Za-z0-9]{1,5}", component_id):
        raise ValueError("CCD component IDs must contain 1–5 letters or digits.")
    component_id = component_id.upper()
    folder = config.DATA_ROOT / "chemistry" / "ccd"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{component_id}.cif"
    url = f"https://files.rcsb.org/ligands/download/{component_id}.cif"
    if path.is_file() and path.stat().st_size > 2 * 1024 * 1024:
        raise ValueError('Cached CCD record exceeds the local 2 MiB bound.')
    if not path.is_file():
        request = urllib.request.Request(url, headers={"User-Agent": "DynaMol/0.1 (CCD chemistry retrieval)"})
        try:
            class NoRedirect(urllib.request.HTTPRedirectHandler):
                def redirect_request(self, req, fp, code, msg, headers, newurl):
                    raise ValueError('CCD retrieval redirects are not permitted.')
            with urllib.request.build_opener(NoRedirect).open(request, timeout=15) as response:
                data = response.read(2 * 1024 * 1024 + 1)
            if len(data) > 2 * 1024 * 1024:
                raise ValueError("CCD record exceeds the local 2 MiB bound.")
            temporary = path.with_suffix('.' + uuid.uuid4().hex + '.tmp')
            temporary.write_bytes(data)
            temporary.replace(path)
            storage.atomic_json(path.with_suffix('.metadata.json'), {'retrieved_at': dt.datetime.now(dt.timezone.utc).isoformat(), 'url': url, 'sha256': hashlib.sha256(data).hexdigest()})
        except Exception as exc:
            raise ValueError(f"Could not retrieve chemistry for CCD {component_id}. Supply a chemically complete SDF/MOL2 or a ligand SMILES override. {exc}") from exc
    import gemmi
    block = gemmi.cif.read_file(str(path)).sole_block()
    if block.find_value("_chem_comp.id").strip("'\"").upper() != component_id:
        raise ValueError("Cached CCD record has the wrong chemical component ID.")
    metadata_path = path.with_suffix('.metadata.json')
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    return block, {"provider": "wwPDB CCD via RCSB", "component_id": component_id, "url": url,
                   "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), 'retrieved_at': metadata.get('retrieved_at'),
                   'retrieval_note': 'Retrieval time unavailable for legacy cache.' if not metadata.get('retrieved_at') else None}


def _set_conformer(mol, positions):
    mol.RemoveAllConformers()
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in enumerate(positions):
        conformer.SetAtomPosition(i, tuple(map(float, xyz)))
    conformer.Set3D(True)
    mol.AddConformer(conformer)


def _assign_cip_labels(mol):
    """Assign accurate CIP (R/S) descriptors from already-perceived stereo tags.

    The caller has just run ``AssignStereochemistryFrom3D``, which sets both the
    atom chiral tags and bond E/Z from geometry. We then apply ``rdCIPLabeler``,
    RDKit's full CIP-algorithm implementation, to read accurate ``_CIPCode``
    values off those tags. We deliberately do NOT call legacy
    ``AssignStereochemistry`` here: its simplified priority rules mislabel
    fused-ring steroid centres (e.g. it reports R for the correctly deposited S
    at C13 of an estrane such as R1881) and its ``cleanIt`` pass strips valid
    chiral tags on other ring-junction centres (e.g. steroid C8), leaving them
    "unresolved". Either behaviour spuriously rejects a valid bound pose against
    the CCD reference; ``rdCIPLabeler`` on the geometry-derived tags agrees with
    the CCD.
    """
    rdCIPLabeler.AssignCIPLabels(mol)


def _graph(residue, coordinates, chemistry=None, block=None):
    atoms = [atom for atom in residue.atoms() if atom.element != app.element.hydrogen]
    names = [atom.name for atom in atoms]
    if len(names) != len(set(names)):
        raise ValueError("Ligand heavy atom names must be unique within each residue.")
    if not atoms or len(atoms) > MAX_LIGAND_HEAVY_ATOMS:
        raise ValueError(f"Ligands are limited to 1–{MAX_LIGAND_HEAVY_ATOMS} heavy atoms.")
    if any(atom.element.symbol not in ORGANIC_ELEMENTS for atom in atoms):
        raise ValueError("This residue contains elements outside the supported organic GAFF2 chemistry. Coordinated metals and organometallic ligands need a specialized model.")
    if chemistry and chemistry.get("canonical_isomeric_smiles"):
        reference = Chem.MolFromSmiles(chemistry["canonical_isomeric_smiles"])
        if reference is not None and any(a.GetNumRadicalElectrons() for a in reference.GetAtoms()):
            raise ValueError("Each ligand must be one connected, closed-shell chemical graph. Radical states require a specialized model; the source state was retained.")
    rw = Chem.RWMol()
    expected_stereo, expected_bonds = {}, {}
    if block is not None:
        atom_table = block.get_mmcif_category("_chem_comp_atom.")
        records = {name: i for i, name in enumerate(atom_table["atom_id"]) if atom_table["type_symbol"][i] != "H"}
        if set(names) != set(records):
            missing = sorted(set(records) - set(names)); extra = sorted(set(names) - set(records))
            raise ValueError(f"CCD heavy-atom identity mismatch (missing {missing}, extra {extra}). Missing ligand heavy atoms will not be guessed.")
        for atom in atoms:
            i = records[atom.name]
            if atom.element.symbol.upper() != atom_table["type_symbol"][i].upper():
                raise ValueError(f"CCD element mismatch for atom {atom.name}.")
            chemical = Chem.Atom(atom.element.symbol)
            chemical.SetFormalCharge(int(atom_table["charge"][i]))
            rw.AddAtom(chemical)
            stereo = atom_table.get("pdbx_stereo_config", ["N"] * len(records))[i]
            if stereo in {"R", "S"} and atom.element.symbol == "C":
                expected_stereo[atom.name] = stereo
        bond_table = block.get_mmcif_category("_chem_comp_bond.")
        for i, a in enumerate(bond_table.get("atom_id_1", [])):
            b = bond_table["atom_id_2"][i]
            if a not in names or b not in names:
                continue
            order = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE,
                     "AROM": Chem.BondType.AROMATIC}.get(bond_table["value_order"][i])
            if order is None:
                raise ValueError("The CCD contains an unsupported bond type.")
            rw.AddBond(names.index(a), names.index(b), order)
            stereo = bond_table.get('pdbx_stereo_config', ['N'] * len(bond_table['atom_id_1']))[i]
            if stereo in {'E', 'Z'}:
                expected_bonds[(names.index(a), names.index(b))] = stereo
    else:
        indices = {atom.index: i for i, atom in enumerate(atoms)}
        for atom in atoms:
            entry = chemistry["atoms"][atom.index]
            chemical = Chem.Atom(atom.element.symbol)
            # MOL2 partial charge is NOT formal charge. Recover conventional
            # charged Tripos atom types where unambiguous, otherwise sanitize.
            charge = entry.get("formal_charge")
            if charge is None:
                kind = entry.get("input_atom_type", "")
                charge = 1 if kind == "N.4" else -1 if kind == "O.co2" else 0
                if kind == "O.co2":
                    # Only singly bonded carboxylate O carries formal −1.
                    if any(atom.index in b["atoms"] and b["order"] == 2 for b in chemistry["bonds"]):
                        charge = 0
            chemical.SetFormalCharge(int(charge))
            rw.AddAtom(chemical)
        for bond in chemistry["bonds"]:
            a, b = bond["atoms"]
            if a in indices and b in indices:
                order = {1.: Chem.BondType.SINGLE, 2.: Chem.BondType.DOUBLE, 3.: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}.get(float(bond["order"]))
                if order is None:
                    raise ValueError("The source has an unsupported ligand bond order.")
                rw.AddBond(indices[a], indices[b], order)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise ValueError("Ligand valence or charge could not be established from the chemical source. Supply an explicit ligand SMILES override.") from exc
    if len(Chem.GetMolFrags(mol)) != 1 or any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        raise ValueError("Each ligand must be one connected, closed-shell chemical graph.")
    _set_conformer(mol, coordinates[[atom.index for atom in atoms]])
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    _assign_cip_labels(mol)
    for name, expected in expected_stereo.items():
        atom = mol.GetAtomWithIdx(names.index(name))
        actual = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        if actual != expected:
            raise ValueError(f"Bound ligand stereochemistry disagrees with CCD at {name}: expected {expected}, found {actual or 'unresolved'}. No coordinates were changed.")
    unverified_bond_stereo = []
    for (a, b), expected in expected_bonds.items():
        bond = mol.GetBondBetweenAtoms(a, b)
        # Compare CIP E/Z descriptors. rdCIPLabeler (run in _assign_cip_labels)
        # sets an accurate CIP "E"/"Z" _CIPCode on each stereo bond. bond.GetStereo()
        # is referenced to whichever stereo atoms RDKit picked (it returns
        # STEREOCIS/STEREOTRANS, not CIP STEREOE/STEREOZ), so it cannot be
        # compared against the CCD's CIP descriptor -- doing so falsely rejects
        # correct geometry (e.g. geldanamycin C2-C3, rapamycin C17-C18).
        actual = bond.GetProp('_CIPCode') if bond.HasProp('_CIPCode') else None
        if actual != expected:
            # CCD can describe C=NH orientation using a hydrogen/lone-pair
            # convention absent from this heavy-atom graph. RDKit requires
            # two heavy neighbors at each end for double-bond stereo:
            # https://www.rdkit.org/docs/RDKit_Book.html#stereogenic-bonds
            # This narrowly records an unverified descriptor; it does not
            # forgive an opposite or unresolved observable stereobond.
            ends = [mol.GetAtomWithIdx(a), mol.GetAtomWithIdx(b)]
            terminal_imine = (bond.GetBondType() == Chem.BondType.DOUBLE and not bond.IsInRing()
                              and sorted(atom.GetSymbol() for atom in ends) == ['C', 'N']
                              and any(atom.GetSymbol() == 'N' and atom.GetDegree() == 1
                                      and atom.GetFormalCharge() == 0 and atom.GetTotalNumHs() == 1 for atom in ends))
            if actual is None and terminal_imine:
                unverified_bond_stereo.append({'atoms': [names[a], names[b]], 'ccd_descriptor': expected,
                                              'status': 'unverified',
                                              'reason': 'Terminal C=NH has only one heavy-atom neighbor at nitrogen; its CCD E/Z descriptor cannot be verified from the retained heavy-atom geometry.'})
                continue
            raise ValueError(f'Bound ligand E/Z stereochemistry disagrees with CCD at {names[a]}–{names[b]}.')
    if unverified_bond_stereo:
        mol.SetProp('_DynaMolUnverifiedCCDBondStereo', json.dumps(unverified_bond_stereo))
    if chemistry and chemistry.get('canonical_isomeric_smiles'):
        reference = Chem.MolFromSmiles(chemistry['canonical_isomeric_smiles'])
        if reference is not None and reference.GetNumAtoms() == mol.GetNumAtoms() and not mol.HasSubstructMatch(reference, useChirality=True):
            raise ValueError('Bound ligand stereochemistry disagrees with the retained original chemical reference.')
    return mol, atoms, names, expected_stereo


def _mapped_variant(mol, smiles):
    selected = Chem.MolFromSmiles(smiles)
    if selected is None:
        raise ValueError("The ligand SMILES override is invalid.")
    selected = Chem.RemoveHs(selected)
    if len(Chem.GetMolFrags(selected)) != 1 or any(a.GetNumRadicalElectrons() for a in selected.GetAtoms()):
        raise ValueError("Ligand overrides must describe one connected, closed-shell state; radicals require a specialized model.")
    if selected.GetNumAtoms() != mol.GetNumAtoms():
        raise ValueError("Selected protonation/SMILES must retain every ligand heavy atom.")
    maps = [atom.GetAtomMapNum() for atom in selected.GetAtoms()]
    if set(maps) == set(range(1, mol.GetNumAtoms() + 1)):
        order = [maps.index(i + 1) for i in range(mol.GetNumAtoms())]
    else:
        # User SMILES can alter charge/bond order, but must have exactly the
        # same element-labeled heavy-atom connectivity as the supplied pose.
        parameters = Chem.AdjustQueryParameters.NoAdjustments()
        parameters.makeBondsGeneric = True
        query = Chem.AdjustQueryProperties(mol, parameters)
        matches = selected.GetSubstructMatches(query, useChirality=False, uniquify=False, maxMatches=128)
        if len(matches) != 1:
            raise ValueError("Ligand SMILES does not map uniquely to the bound heavy-atom graph. Supply an explicit atom-mapped SMILES using bound heavy-atom indices.")
        order = list(matches[0])
    selected = Chem.RenumberAtoms(selected, order)
    for atom in selected.GetAtoms():
        atom.SetAtomMapNum(0)
    original_edges = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))) for bond in mol.GetBonds()}
    selected_edges = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))) for bond in selected.GetBonds()}
    if original_edges != selected_edges or [a.GetSymbol() for a in mol.GetAtoms()] != [a.GetSymbol() for a in selected.GetAtoms()]:
        raise ValueError("Selected ligand state changes heavy-atom connectivity or identity.")
    Chem.AssignStereochemistry(selected, cleanIt=True, force=True)
    wanted_bonds = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetStereo() for b in selected.GetBonds() if b.GetStereo() in {Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ}}
    wanted = {a.GetIdx(): a.GetProp("_CIPCode") for a in selected.GetAtoms() if a.GetSymbol() == "C" and a.HasProp("_CIPCode")}
    _set_conformer(selected, mol.GetConformer().GetPositions())
    Chem.AssignStereochemistryFrom3D(selected, replaceExistingTags=True)
    Chem.AssignStereochemistry(selected, cleanIt=True, force=True)
    for index, desired in wanted.items():
        actual = selected.GetAtomWithIdx(index)
        if not actual.HasProp("_CIPCode") or actual.GetProp("_CIPCode") != desired:
            raise ValueError("Ligand SMILES stereochemistry disagrees with the supplied bound pose.")
    for (a, b), desired in wanted_bonds.items():
        if selected.GetBondBetweenAtoms(a, b).GetStereo() != desired:
            raise ValueError('Ligand SMILES E/Z stereochemistry disagrees with the supplied bound pose.')
    return selected


def _select_state(mol, names, component_id, ph, override):
    original = Chem.MolToSmiles(mol, isomericSmiles=True)
    warning = "Ligand protonation is a fixed model choice; it is not a binding-site pKa calculation or a tautomer ranking."
    reference_state = None
    if override:
        selected = _mapped_variant(mol, override)
        method = "Explicit user SMILES protonation/tautomer state"
        candidates = [Chem.MolToSmiles(selected, isomericSmiles=True)]
    elif component_id == "SO4":
        # Dimorphite's normalization turns isolated sulfate into sulfuric acid,
        # removing two oxygen maps without restoring the dianion at pH 7.
        # Preserve the exact named CCD state instead of guessing those maps or
        # transferring the inappropriate neutralization result to bound atoms.
        if not 6 <= ph <= 8:
            raise ValueError("Automatic sulfate preparation retains the CCD SO4 dianion only at pH 6–8. Supply an explicit atom-mapped sulfate/bisulfate state at this pH; no protonation state was changed.")
        sulfur = [a for a in mol.GetAtoms() if a.GetSymbol() == 'S']
        oxygens = [a for a in mol.GetAtoms() if a.GetSymbol() == 'O']
        valid = (mol.GetNumAtoms() == 5 and mol.GetNumBonds() == 4
                 and len(sulfur) == 1 and len(oxygens) == 4
                 and sulfur[0].GetDegree() == 4 and sulfur[0].GetFormalCharge() == 0
                 and Chem.GetFormalCharge(mol) == -2
                 and all(a.GetTotalNumHs() == 0 and a.GetNumRadicalElectrons() == 0 for a in mol.GetAtoms()))
        if valid:
            valid = all(a.GetDegree() == 1 and a.GetNeighbors()[0].GetIdx() == sulfur[0].GetIdx()
                        and (a.GetFormalCharge(), next(iter(a.GetBonds())).GetBondType())
                        in {(-1, Chem.BondType.SINGLE), (0, Chem.BondType.DOUBLE)} for a in oxygens)
        if not valid:
            raise ValueError("The supplied SO4 graph does not match the CCD sulfate dianion reference state. Supply an explicit chemically defined ligand state; no source charges or atom identities were changed.")
        selected = Chem.Mol(mol)
        method = "Retained wwPDB CCD SO4 sulfate dianion (−2), fixed reference-state policy at pH 6–8; no pH heuristic or atom-map recovery applied"
        candidates = [original]
        reference_state = {"component_id": "SO4", "url": "https://www.rcsb.org/ligand/SO4",
                           "formal_charge": -2, "automatic_ph_range": [6, 8],
                           "original_graph_retained": True, "ph_heuristic_applied": False}
        warning += " The existing sulfate dianion is retained as a named near-neutral reference state; local sulfate/bisulfate populations are not predicted."
    elif component_id == "GNP" and 6 <= ph <= 8:
        selected = Chem.Mol(mol)
        terminal_oxygens = []
        for atom in selected.GetAtoms():
            if atom.GetSymbol() == "O" and atom.GetDegree() == 1:
                bond = next(iter(atom.GetBonds()))
                if bond.GetBondType() == Chem.BondType.SINGLE and atom.GetNeighbors()[0].GetSymbol() == "P":
                    atom.SetFormalCharge(-1); atom.SetNumExplicitHs(0); atom.SetNoImplicit(True)
                    terminal_oxygens.append(names[atom.GetIdx()])
        Chem.SanitizeMol(selected)
        if len(terminal_oxygens) != 4 or Chem.GetFormalCharge(selected) != -4:
            raise ValueError("GNP −4 reference-state rule did not match four terminal phosphate hydroxyls.")
        method = "GppNHp/GNP −4 reference state at pH 6–8; imido N–H retained; model assumption supported by doi:10.1074/jbc.RA117.001110"
        candidates = [Chem.MolToSmiles(selected, isomericSmiles=True)]
    else:
        from dimorphite_dl import protonate_smiles
        mapped = Chem.Mol(mol)
        for i, atom in enumerate(mapped.GetAtoms()):
            atom.SetAtomMapNum(i + 1)
        raw = protonate_smiles(Chem.MolToSmiles(mapped, isomericSmiles=True), ph_min=ph, ph_max=ph,
                               precision=0.0, max_variants=16)
        # Preserve the input chemical graph and stereotags. Dimorphite may
        # reorder and alter stereotags while enumerating mapped molecules.
        # Transfer only mapped formal charges and H counts at supported sites.
        amide = Chem.MolFromSmarts('[N]-[C](=[O,S])')
        protected_n = {match[0] for match in mol.GetSubstructMatches(amide)}
        protected_n.update(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'N'
                           and any(n.GetSymbol() == 'N' for n in a.GetNeighbors()))
        variants, filtered, restored_maps = [], set(), False
        for smiles in raw:
            proposed = Chem.MolFromSmiles(smiles)
            if proposed is None:
                raise ValueError('The pH heuristic produced invalid ligand chemistry.')
            if {a.GetAtomMapNum() for a in proposed.GetAtoms()} != set(range(1, mol.GetNumAtoms() + 1)):
                # Dimorphite's normalization can omit an existing phosphate-O
                # map (observed for CCD NAD). Recover it only from a unique
                # element/connectivity mapping constrained by every surviving
                # map, never from atom order or geometric proximity.
                proposed = _restore_protonation_maps(mapped, proposed)
                restored_maps = True
            by_map = {a.GetAtomMapNum(): a for a in proposed.GetAtoms()}
            if set(by_map) != set(range(1, mol.GetNumAtoms() + 1)):
                raise ValueError('The pH heuristic lost the ligand atom identity map.')
            variant = Chem.Mol(mol)
            for index, atom in enumerate(variant.GetAtoms()):
                proposal = by_map[index + 1]
                if proposal.GetFormalCharge() != atom.GetFormalCharge():
                    if index in protected_n:
                        filtered.add(names[index]); continue
                    atom.SetFormalCharge(proposal.GetFormalCharge())
                    atom.SetNumExplicitHs(proposal.GetTotalNumHs())
                    atom.SetNoImplicit(True)
            Chem.SanitizeMol(variant)
            variants.append(variant)
        candidates = sorted({Chem.MolToSmiles(variant, isomericSmiles=True) for variant in variants})
        if len(candidates) != 1:
            raise ValueError(f"Ligand protonation produced {len(candidates)} unranked states. Supply an explicit ligand SMILES state; candidates: {'; '.join(candidates)}")
        selected = variants[0]
        if filtered:
            warning += " Unsupported amide/N–N protonation proposals at " + ', '.join(sorted(filtered)) + " were filtered; source states retained. Review these sites or supply an explicit ligand SMILES state."
        if restored_maps:
            warning += " A missing protonation-tool atom map was restored by a unique element/connectivity match constrained by all retained maps."
        method = "Dimorphite-DL 2 empirical site rules at requested pH, pKa precision=0; amide/N–N sites retained from source; original graph/stereotags preserved; no microscopic pKa or tautomer ranking"
    return selected, {"original_smiles": original, "original_formal_charge": Chem.GetFormalCharge(mol),
                      "selected_smiles": Chem.MolToSmiles(selected, isomericSmiles=True),
                      "formal_charge": Chem.GetFormalCharge(selected), "candidate_smiles": candidates,
                      "protonation_method": method, "ph": ph, "warnings": [warning],
                      **({"reference_state": reference_state} if reference_state else {})}


def _restore_protonation_maps(original, proposed):
    """Recover omitted maps only when graph identity has one possible solution."""
    if original.GetNumAtoms() != proposed.GetNumAtoms():
        raise ValueError("The pH heuristic changed ligand atom count.")
    known = [a.GetAtomMapNum() for a in proposed.GetAtoms() if a.GetAtomMapNum()]
    if len(known) != len(set(known)) or any(i < 1 or i > original.GetNumAtoms() for i in known):
        raise ValueError("The pH heuristic produced conflicting ligand atom identity maps.")
    query = Chem.RWMol()
    for atom in original.GetAtoms():
        query.AddAtom(Chem.AtomFromSmarts(f"[#{atom.GetAtomicNum()}]"))
    for bond in original.GetBonds():
        query.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.UNSPECIFIED)
        query.ReplaceBond(query.GetNumBonds()-1, Chem.BondFromSmarts('~'))
    matches = []
    candidates = proposed.GetSubstructMatches(query.GetMol(), uniquify=False, maxMatches=4096)
    if len(candidates) >= 4096:
        raise ValueError("Protonation atom-map recovery exceeded its bounded graph search. Supply an explicit atom-mapped state.")
    for match in candidates:
        if all(proposed.GetAtomWithIdx(j).GetAtomMapNum() in (0, i+1) for i,j in enumerate(match)):
            matches.append(match)
            if len(matches) > 1:
                break
    if len(matches) != 1:
        raise ValueError("The pH heuristic lost the ligand atom identity map and it cannot be restored uniquely. Supply an explicit atom-mapped state.")
    result = Chem.Mol(proposed)
    for i,j in enumerate(matches[0]):
        result.GetAtomWithIdx(j).SetAtomMapNum(i+1)
    edges = {tuple(sorted((result.GetAtomWithIdx(b.GetBeginAtomIdx()).GetAtomMapNum(), result.GetAtomWithIdx(b.GetEndAtomIdx()).GetAtomMapNum()))) for b in result.GetBonds()}
    expected = {tuple(sorted((b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1))) for b in original.GetBonds()}
    if edges != expected:
        raise ValueError("The pH heuristic changed ligand heavy-atom connectivity.")
    return result


def _graph_from_override(residue, coordinates, smiles):
    """Map a user-supplied chemical graph to a pose, refusing ambiguity."""
    from rdkit.Chem import rdDetermineBonds
    atoms = [a for a in residue.atoms() if a.element != app.element.hydrogen]
    names = [a.name for a in atoms]
    selected = Chem.MolFromSmiles(smiles)
    if selected is None:
        raise ValueError('The ligand SMILES override is invalid.')
    selected = Chem.RemoveHs(selected)
    if len(Chem.GetMolFrags(selected)) != 1 or any(a.GetNumRadicalElectrons() for a in selected.GetAtoms()):
        raise ValueError("Ligand overrides must describe one connected, closed-shell state; radicals require a specialized model.")
    if selected.GetNumAtoms() != len(atoms) or len(set(names)) != len(names):
        raise ValueError('The ligand SMILES must match every uniquely named bound heavy atom.')
    if len(atoms) > MAX_LIGAND_HEAVY_ATOMS or any(a.element.symbol not in ORGANIC_ELEMENTS for a in atoms):
        raise ValueError('The ligand exceeds supported GAFF2 organic chemistry/size limits.')
    mapped = [a.GetAtomMapNum() for a in selected.GetAtoms()]
    coordinates = coordinates[[a.index for a in atoms]]
    if set(mapped) == set(range(1, len(atoms) + 1)):
        selected = Chem.RenumberAtoms(selected, [mapped.index(i + 1) for i in range(len(atoms))])
        mapping_method = 'Explicit atom-map numbers refer to the bound heavy-atom order (1-based).'
    else:
        geometry = Chem.RWMol()
        for atom in atoms:
            geometry.AddAtom(Chem.Atom(atom.element.symbol))
        geometry = geometry.GetMol()
        _set_conformer(geometry, coordinates)
        rdDetermineBonds.DetermineConnectivity(geometry, useHueckel=False, covFactor=1.3)
        parameters = Chem.AdjustQueryParameters.NoAdjustments()
        parameters.makeBondsGeneric = True
        query = Chem.AdjustQueryProperties(geometry, parameters)
        matches = selected.GetSubstructMatches(query, useChirality=False, uniquify=False, maxMatches=2)
        if len(matches) != 1:
            raise ValueError('The explicit ligand SMILES cannot be mapped uniquely to the bound pose. Supply an atom-mapped SMILES using the displayed heavy-atom order; no mapping was guessed.')
        selected = Chem.RenumberAtoms(selected, list(matches[0]))
        mapping_method = 'Unique element/connectivity match between explicit SMILES and bound-distance connectivity; chemical bond orders/charge come exclusively from SMILES.'
    if [a.GetSymbol() for a in selected.GetAtoms()] != [a.element.symbol for a in atoms]:
        raise ValueError('The mapped ligand SMILES element identities disagree with the bound atom order.')
    for bond in selected.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        distance = np.linalg.norm(coordinates[a] - coordinates[b])
        if not .8 < distance < 2.5:
            raise ValueError('The mapped ligand SMILES contains a bond inconsistent with the bound pose.')
    for atom in selected.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(selected, cleanIt=True, force=True)
    expected_bonds = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetStereo() for b in selected.GetBonds() if b.GetStereo() in {Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ}}
    expected = {a.GetIdx(): a.GetProp('_CIPCode') for a in selected.GetAtoms() if a.GetSymbol() == 'C' and a.HasProp('_CIPCode')}
    _set_conformer(selected, coordinates)
    Chem.AssignStereochemistryFrom3D(selected, replaceExistingTags=True)
    Chem.AssignStereochemistry(selected, cleanIt=True, force=True)
    for i, configuration in expected.items():
        atom = selected.GetAtomWithIdx(i)
        if not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode') != configuration:
            raise ValueError('The explicit ligand SMILES stereochemistry disagrees with the bound pose.')
    for (a, b), desired in expected_bonds.items():
        if selected.GetBondBetweenAtoms(a, b).GetStereo() != desired:
            raise ValueError('The explicit ligand SMILES E/Z stereochemistry disagrees with the bound pose.')
    return selected, atoms, names, {names[i]: value for i, value in expected.items()}, mapping_method


def _original_connection_blocks(dataset_id, source_cif, component_ids):
    """Do not lose covalent ligand links during canonical PDB normalization."""
    blocked = set()
    excluded = set(STANDARD_RESIDUES) | storage.WATERS | SUPPORTED_IONS
    # A derived one-chain structure retains original source records as evidence.
    # Scope link screening to its retained residue connections: a covalent LIG
    # on an excluded chain must not block an independent LIG of the same name.
    from .monomers import _ancestry
    for _, records in _ancestry(dataset_id):
        scoped = next((record["monomer_selection"] for record in records if "monomer_selection" in record), None)
        if scoped is not None:
            for connection in scoped.get("retained_covalent_connections", []):
                for residue in connection["residues"]:
                    key = ":".join([residue["chain"], residue["resid"], residue.get("insertion_code", ""), residue["name"]])
                    name = component_ids.get(key, residue["name"])
                    if name not in excluded:
                        blocked.add(name)
            return blocked
    if source_cif:
        import gemmi
        table = gemmi.cif.read_file(str(source_cif)).sole_block().get_mmcif_category('_struct_conn.')
        if table:
            for i, kind in enumerate(table.get('conn_type_id', [])):
                if kind.lower() in {'metalc', 'hydrog', 'saltbr'}:
                    continue
                names = [table.get(f'ptnr{partner}_label_comp_id', [''] * len(table['conn_type_id']))[i] for partner in [1, 2]]
                # Modified polymer chemistry must not be treated as a removable
                # independent ligand either; its component remains blocked.
                blocked.update(name for name in names if name not in excluded)
    pending, visited = [dataset_id], set()
    while pending and len(visited) < 20:
        current = pending.pop(0)
        if current in visited:
            continue
        visited.add(current)
        folder = storage.dataset_dir(current)
        for path in [folder / 'source.pdb', folder / 'sequence-source.pdb', *sorted((folder / 'originals').glob('*.pdb'))]:
            if not path.is_file():
                continue
            for line in path.read_text(errors='replace').splitlines():
                if line.startswith('LINK  '):
                    names = [line[17:20].strip(), line[47:50].strip()]
                    if any(name in SUPPORTED_IONS for name in names):
                        continue
                    blocked.update(name for name in names if name not in excluded)
        for filename in ['metadata.json', 'provenance.json']:
            path = folder / filename
            record = json.loads(path.read_text()) if path.is_file() else {}
            if not isinstance(record, dict):
                continue
            for state in [record, record.get('preparation'), record.get('solvation'), record.get('preparation_state')]:
                if isinstance(state, dict):
                    parent = state.get('parent_dataset_id') or state.get('source_dataset_id')
                    if isinstance(parent, str):
                        pending.append(parent)
    return blocked


def _models(dataset_id, ph=7.0, overrides=None, actions=None, *, execute_repair=False, seed=2026, on_progress=None, check_cancel=None, environment_indices=None):
    from .preparation import exact_input_path
    from .modified_residues import register_topology_definitions
    from .residue_identity import protein_residue_keys
    register_topology_definitions()
    loaded = app.PDBFile(str(exact_input_path(dataset_id)))
    from .preparation import canonicalize_protonation_aliases
    canonicalize_protonation_aliases(loaded.topology)
    protein_keys = protein_residue_keys(dataset_id, loaded.topology, loaded.positions)
    coordinates = np.asarray(loaded.positions.value_in_unit(unit.angstrom))
    components, source_cif = _component_ids(dataset_id, loaded.topology)
    blocked_components = _original_connection_blocks(dataset_id, source_cif, components)
    path = storage.dataset_dir(dataset_id) / "chemistry.json"
    chemistry = json.loads(path.read_text()) if path.is_file() else None
    if chemistry and (not chemistry.get("authoritative_bonds") or len(chemistry["atoms"]) != loaded.topology.getNumAtoms()):
        chemistry = None
    models = []
    actions = actions or {}
    if not isinstance(actions, dict) or any(value not in {"repair", "remove"} for value in actions.values()):
        raise ValueError("Ligand actions must explicitly select repair or remove for a residue key.")
    visited_actions = set()
    if environment_indices is None:
        environment_indices = [atom.index for atom in loaded.topology.atoms() if atom.element != app.element.hydrogen and actions.get(_key(atom.residue)) != "remove"]
    ion_records = {record["key"]: record for record in inspect_ions(loaded.topology)}
    for residue in loaded.topology.residues():
        if residue_key(residue) in protein_keys or residue.name in STANDARD_RESIDUES or residue.name.upper() in storage.WATERS:
            continue
        if residue.name.upper() in SUPPORTED_IONS:
            if not ion_records[_key(residue)]["error"]:
                continue
        key = _key(residue)
        action = actions.get(key)
        if action:
            visited_actions.add(key)
        description = {"key": key, "chain": residue.chain.id, "resid": residue.id, "insertion_code": residue.insertionCode,
                       "residue": residue.name, "component_id": components.get(key, residue.name), "error": None,
                       "atom_indices": [atom.index for atom in residue.atoms()], "selected_action": action,
                       "can_remove": False, "can_repair": False, "removed": False}
        try:
            if key in ion_records and ion_records[key]["error"]:
                raise ValueError(ion_records[key]["error"])
            for bond in loaded.topology.bonds():
                a, b = bond
                if (a.residue == residue) != (b.residue == residue):
                    other = b if a.residue == residue else a
                    if other.residue.name.upper() not in SUPPORTED_IONS:
                        raise ValueError('Covalently linked ligand residues need a specialized force-field model; no molecule was removed.')
            component_id = components.get(key, residue.name)
            if component_id in blocked_components or residue.name in blocked_components:
                raise ValueError('Original structure records a covalent ligand/polymer link. This ligand requires a specialized covalent force-field model; no molecule was removed.')
            description["can_remove"] = True
            if action == "remove":
                description.update(removed=True, removal_reason="Explicit per-residue user choice; the complete noncovalent residue is excluded.")
                models.append(LigandModel(None, residue, [], [], description))
                continue
            override = (overrides or {}).get(key) or (overrides or {}).get(component_id)
            if isinstance(override, dict):
                override = override.get("smiles")
            explicit_model = False
            atom_indices = None
            incomplete = False
            try:
                block, reference = (None, {"provider": "Retained authoritative source chemical graph"}) if chemistry else _ccd(component_id)
                if block is not None:
                    from .ligand_repair import inspect_repair, repair_ligand
                    repair = inspect_repair(residue, coordinates, block)
                    description.update(repair)
                    incomplete = bool(repair["missing_heavy_atoms"] or repair["extra_heavy_atoms"])
                    if incomplete:
                        if action != "repair" or not repair["can_repair"]:
                            raise ValueError(f"CCD heavy-atom identity mismatch (missing {repair['missing_heavy_atoms']}, extra {repair['extra_heavy_atoms']}). " + (repair["repair_reason"] or "Review ligand identity.") + (" Select modeled repair or explicit removal for this residue." if repair["can_repair"] else ""))
                        if not execute_repair:
                            description.update(reference=reference, heavy_atoms=len(repair["missing_heavy_atoms"]) + len([atom for atom in residue.atoms() if atom.element != app.element.hydrogen]),
                                               repair_pending=True, warnings=["Missing ligand atoms will be modeled only when preparation runs; observed heavy atoms remain fixed. Modeled positions require review."])
                            models.append(LigandModel(None, residue, [], [], description))
                            continue
                        mol, names, atom_indices, stereocenters, repair_report = repair_ligand(residue, coordinates, block, seed=seed, on_progress=on_progress, check_cancel=check_cancel, environment_indices=environment_indices)
                        atoms = list(mol.GetAtoms())
                        description.update(repair=repair_report, repair_pending=False)
                if not incomplete:
                    mol, atoms, names, stereocenters = _graph(residue, coordinates, chemistry, block)
            except ValueError:
                if not override or incomplete:
                    raise
                mol, atoms, names, stereocenters, mapping_method = _graph_from_override(residue, coordinates, override)
                reference = {'provider': 'Explicit user SMILES', 'mapping_method': mapping_method}
                explicit_model = True
            if explicit_model:
                selected = mol
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                protonation = {'original_smiles': None, 'original_formal_charge': None,
                               'selected_smiles': smiles, 'formal_charge': Chem.GetFormalCharge(mol),
                               'candidate_smiles': [smiles], 'ph': ph,
                               'protonation_method': 'Explicit user SMILES defines chemistry and fixed protonation state.',
                               'warnings': ['Original input ligand chemistry was unavailable or did not match the bound pose; the explicit user chemical reference was used. It is not a binding-site pKa assignment.']}
            else:
                selected, protonation = _select_state(mol, names, component_id, ph, override)
            description.update(**protonation, heavy_atoms=len(atoms), reference=reference,
                               stereo_checked=stereocenters, source_cif=str(source_cif) if source_cif else None,
                               stereo_limitations='Carbon tetrahedral CCD/input-SMILES configurations are checked. Nitrogen inversion, atropisomerism, and phosphorus stereodescriptors after ionization are not validated; the original bound heavy-atom pose is retained.')
            if mol.HasProp('_DynaMolUnverifiedCCDBondStereo'):
                description['stereo_unverified'] = json.loads(mol.GetProp('_DynaMolUnverifiedCCDBondStereo'))
                for item in description['stereo_unverified']:
                    description['warnings'].append('Unverified CCD ' + item['ccd_descriptor'] + ' descriptor at ' + '–'.join(item['atoms']) + ': ' + item['reason'] + ' No E/Z assignment or heavy-atom coordinate change was made; review the selected ligand state.')
            if description.get("repair"):
                description["warnings"] = description.get("warnings", []) + description["repair"]["warnings"]
            models.append(LigandModel(selected, residue, atom_indices if atom_indices is not None else [atom.index for atom in atoms], names, description))
        except Exception as exc:
            description["error"] = str(exc)
            models.append(LigandModel(None, residue, [], [], description))
    if set(actions) != visited_actions:
        raise ValueError("A selected ligand action does not match a current nonprotein ligand residue: " + ", ".join(sorted(set(actions) - visited_actions)))
    return models


def inspect_ligands(dataset_id, ph=7.0, overrides=None, actions=None):
    return [model.description for model in _models(dataset_id, ph, overrides, actions)]


def _add_hydrogens(model):
    mol = Chem.AddHs(model.mol, addCoords=True)
    coordinates = mol.GetConformer().GetPositions().copy()
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp('_TriposAtomName', model.names[i] if i < len(model.names) else f"H{i - len(model.names) + 1}")
    # With no newly placed H there are no movable coordinates. RDKit's BFGS
    # optimizer can fail on an all-fixed force field (for example sulfate).
    if mol.GetNumAtoms() == model.mol.GetNumAtoms():
        return mol, "No hydrogens added for the selected ligand state; original coordinates retained; relaxation unnecessary"
    # Relax only newly placed H. Bound ligand heavy coordinates are fixed.
    if AllChem.MMFFHasAllMoleculeParams(mol):
        parameters = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        force = AllChem.MMFFGetMoleculeForceField(mol, parameters)
        method = "MMFF94s hydrogen-only relaxation"
    elif AllChem.UFFHasAllMoleculeParams(mol):
        force = AllChem.UFFGetMoleculeForceField(mol)
        method = "UFF hydrogen-only relaxation"
    else:
        force, method = None, "RDKit valence-based hydrogen placement; no H relaxation parameters"
    if force:
        for i in range(len(model.names)):
            force.AddFixedPoint(i)
        force.Initialize(); force.Minimize(maxIts=300)
    if not np.allclose(mol.GetConformer().GetPositions()[:len(model.names)], coordinates[:len(model.names)], atol=1e-9, rtol=0):
        raise ValueError("Hydrogen placement moved bound ligand heavy atoms.")
    return mol, method


def _run(command, folder, environment, on_progress, check_cancel, timeout=600):
    log = folder / (Path(command[0]).name + ".log")
    if on_progress:
        on_progress(f"{Path(command[0]).name}: {folder.name}")
    started = time.monotonic()
    commands_file = folder / 'commands.json'
    commands = json.loads(commands_file.read_text()) if commands_file.is_file() else []
    commands.append({'argv': command, 'started_at': dt.datetime.now(dt.timezone.utc).isoformat(), 'timeout_seconds': timeout})
    storage.atomic_json(commands_file, commands)
    with log.open('w') as handle:
        process = subprocess.Popen(command, cwd=folder, env=environment, stdout=handle, stderr=subprocess.STDOUT)
        try:
            while process.poll() is None:
                if check_cancel:
                    check_cancel()
                if time.monotonic() - started > timeout:
                    raise ValueError(f"{Path(command[0]).name} exceeded its {timeout}-second local limit. Review {log.name}.")
                time.sleep(.15)
        except BaseException:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill(); process.wait()
            raise
    if process.returncode:
        tail = log.read_text(errors="replace")[-1600:]
        raise ValueError(f"{Path(command[0]).name} failed while assigning ligand parameters. {tail}")


def _normalize_native_charge(folder, target):
    """Correct only proven three-decimal SQM-output rounding, retaining evidence."""
    sqm = (folder / 'sqm.out').read_text(errors='replace')
    if 'geometry converged' not in sqm or 'Total Mulliken Charge' not in sqm:
        raise ValueError('AM1 charge calculation did not report geometry convergence and a total Mulliken charge.')
    declared = re.search(r'Total Mulliken Charge\s*=\s*([-+0-9.]+)', sqm)
    section = sqm.split('Mulliken Charge', 1)[1].split('Total Mulliken', 1)[0]
    printed = [float(m.group(1)) for m in re.finditer(r'^\s*\d+\s+[A-Za-z]+\s+([-+0-9.]+)\s*$', section, re.M)]
    lines = (folder / 'charged.mol2').read_text().splitlines()
    rows, active = [], False
    for i, line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):
            active = line == '@<TRIPOS>ATOM'; continue
        if active and line.strip():
            fields = line.split()
            rows.append((i, fields, float(fields[8])))
    total = sum(row[2] for row in rows)
    if not declared or abs(float(declared.group(1)) - target) > .0001 or len(printed) != len(rows):
        raise ValueError('Native SQM charge evidence does not match the ligand formal charge or atom count.')
    if abs(sum(printed) - total) > 5e-5 or abs(target - total) > .005:
        raise ValueError('Native AM1-BCC charge residual cannot be attributed to bounded SQM text rounding.')
    correction = (target - total) / len(rows)
    shutil.copy2(folder / 'charged.mol2', folder / 'charged-native.mol2')
    for i, fields, charge in rows:
        fields[8] = f'{charge + correction:.10f}'
        lines[i] = ' '.join(fields)
    (folder / 'charged.mol2').write_text('\n'.join(lines) + '\n')
    report = {'native_charge_sum_e': total, 'sqm_declared_charge_e': float(declared.group(1)),
              'sqm_printed_atom_charge_sum_e': sum(printed), 'target_charge_e': target,
              'correction_per_atom_e': correction, 'max_total_correction_e': .005,
              'method': 'Uniform numerical correction of demonstrated three-decimal SQM Mulliken-output rounding; native AM1-BCC charge differences retained.',
              'native_file': 'charged-native.mol2', 'normalized_file': 'charged.mol2'}
    storage.atomic_json(folder / 'charge-rounding.json', report)


def _unused_gaff14_defaults(structure, parameters):
    """Use Amber's compatible XML defaults only when no 1–4 pair can exist."""
    if parameters.dihedral_types:
        return None
    neighbors = [set() for _ in structure.atoms]
    for bond in structure.bonds:
        neighbors[bond.atom1.idx].add(bond.atom2.idx)
        neighbors[bond.atom2.idx].add(bond.atom1.idx)
    for start in range(len(neighbors)):
        seen, frontier = {start}, {start}
        for distance in range(1, 4):
            frontier = {neighbor for index in frontier for neighbor in neighbors[index]} - seen
            if distance == 3 and frontier:
                return None
            seen.update(frontier)
    # Generic ParameterSet.from_structure defaults to 1/1 if there are no
    # torsions from which to recover Amber scaling. These global XML fields
    # still must agree with the protein even when unused by this ligand.
    from parmed.amber import AmberParameterSet
    amber = AmberParameterSet()
    original = {"scee": parameters.default_scee, "scnb": parameters.default_scnb}
    parameters.default_scee, parameters.default_scnb = amber.default_scee, amber.default_scnb
    return {"method": "Amber global XML 1–4 defaults for a ligand with no proper torsion types and no atom pairs at shortest bond-graph distance three; no interaction uses the changed defaults.",
            "shortest_distance_three_pairs": 0, "proper_torsion_types": 0,
            "original_generic_defaults": original,
            "amber_defaults": {"scee": parameters.default_scee, "scnb": parameters.default_scnb},
            "native_energy_force_validation_required": True}


def _parameterize(mol, folder, namespace, on_progress, check_cancel):
    import parmed
    from parmed.parameters import ParameterSet
    from parmed.modeller import ResidueTemplate
    from parmed.openmm import OpenMMParameterSet
    amber = ambertools_home()
    folder.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, AMBERHOME=str(amber), OMP_NUM_THREADS=str(config.CPU_THREADS),
                       OPENBLAS_NUM_THREADS=str(config.CPU_THREADS), VECLIB_MAXIMUM_THREADS=str(config.CPU_THREADS))
    environment['PATH'] = str(amber / 'bin') + os.pathsep + environment.get('PATH', '')
    # Sequential synthetic names avoid apostrophes/LEaP truncation; output names
    # are restored by the checked atom-order map before publication.
    names = [f"{atom.GetSymbol()}{i + 1}" for i, atom in enumerate(mol.GetAtoms())]
    input_mol = Chem.Mol(mol)
    input_mol.SetProp('_Name', 'LIG')
    (folder / "input.sdf").write_text(Chem.MolToMolBlock(input_mol) + '\n$$$$\n')
    run = lambda args: _run([str(amber / 'bin' / args[0]), *args[1:]], folder, environment, on_progress, check_cancel)
    charge = Chem.GetFormalCharge(mol)
    run(['antechamber', '-i', 'input.sdf', '-fi', 'sdf', '-o', 'charged.mol2', '-fo', 'mol2', '-c', 'bcc',
         '-nc', str(charge), '-at', 'gaff2', '-rn', 'LIG', '-s', '2', '-pf', 'n', '-dr', 'yes', '-seq', 'n'])
    if not (folder / 'charged.mol2').is_file():
        raise ValueError("Antechamber did not produce an AM1-BCC ligand file.")
    _normalize_native_charge(folder, charge)
    run(['parmchk2', '-i', 'charged.mol2', '-f', 'mol2', '-o', 'ligand.frcmod', '-s', 'gaff2'])
    if 'ATTN' in (folder / 'ligand.frcmod').read_text():
        raise ValueError("GAFF2 could not assign all ligand parameters; parmchk2 marked missing terms. No generic replacement was used.")
    (folder / 'leap.in').write_text('source leaprc.gaff2\nloadamberparams ligand.frcmod\nlig = loadmol2 charged.mol2\ncheck lig\nsaveamberparm lig ligand.prmtop ligand.inpcrd\nquit\n')
    run(['tleap', '-f', 'leap.in'])
    if not (folder / 'ligand.prmtop').is_file():
        raise ValueError("LEaP did not produce a complete ligand force field.")
    structure = parmed.load_file(str(folder / 'ligand.prmtop'), xyz=str(folder / 'ligand.inpcrd'))
    if len(structure.atoms) != mol.GetNumAtoms():
        raise ValueError("AmberTools changed the ligand atom count.")
    for i, atom in enumerate(structure.atoms):
        if atom.atomic_number != mol.GetAtomWithIdx(i).GetAtomicNum():
            raise ValueError("AmberTools changed the ligand atom order/elements.")
    original_bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    amber_bonds = {tuple(sorted((b.atom1.idx, b.atom2.idx))) for b in structure.bonds}
    if original_bonds != amber_bonds:
        raise ValueError("AmberTools changed ligand connectivity.")
    # SDF atom order is preserved by antechamber; verify actual coordinates too.
    if not np.allclose(structure.coordinates, mol.GetConformer().GetPositions(), atol=.002, rtol=0):
        raise ValueError("AmberTools changed the ligand coordinates or atom mapping.")
    charges = [float(atom.charge) for atom in structure.atoms]
    if not np.isfinite(charges).all() or abs(sum(charges) - charge) > 1e-4:
        raise ValueError("AM1-BCC charges are non-finite or do not sum to the selected formal charge.")
    for i, atom in enumerate(structure.atoms):
        atom.name = names[i]
        # Per-atom classes keep native Amber torsion assignments unambiguous,
        # including chemically similar atoms around improper centers.
        atom.type = namespace + '_' + str(i)
        atom.atom_type = copy.copy(atom.atom_type)
        atom.atom_type.name = atom.type
    structure.residues[0].name = namespace
    params = ParameterSet.from_structure(structure, allow_unequal_duplicates=False)
    # ParmEd from_structure inserts every outer-atom permutation for an
    # improper. Retain only the exact native Amber quartet: those permutations
    # have different energy/forces away from perfect planarity.
    params.improper_periodic_types.clear()
    for dihedral in structure.dihedrals:
        if dihedral.improper:
            atoms = [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]
            # Amber can reverse a quartet to avoid a negative zero atom index
            # in prmtop's sign-coded improper records. Reversing all four atoms
            # preserves the signed torsion and its energy/forces exactly.
            if all(atom in atoms[1].bond_partners for atom in (atoms[0], atoms[2], atoms[3])):
                atoms.reverse()
            if not all(atom in atoms[2].bond_partners for atom in [atoms[0], atoms[1], atoms[3]]):
                raise ValueError('An Amber improper does not use the expected central third atom.')
            key = tuple(atom.type for atom in atoms)
            if key in params.improper_periodic_types and params.improper_periodic_types[key] != dihedral.type:
                raise ValueError('Multiple unequal terms for one Amber improper require a specialized converter.')
            params.improper_periodic_types[key] = copy.copy(dihedral.type)
    template = ResidueTemplate.from_residue(structure.residues[0])
    params.residues[namespace] = template
    scale_compatibility = _unused_gaff14_defaults(structure, params)
    if scale_compatibility:
        storage.atomic_json(folder / 'nonbonded-scale-compatibility.json', scale_compatibility)
    converted = OpenMMParameterSet.from_parameterset(params, remediate_residues=False, unique_atom_types=True)
    xml_path = folder / 'ligand.xml'
    converted.write(str(xml_path), improper_dihedrals_ordering='amber')
    _validate_parameter_conversion(folder, xml_path)
    return xml_path, names, charges


def _validate_parameter_conversion(folder, xml_path):
    """Validate each new ligand against its native Amber system before reuse."""
    import openmm as mm
    from xml.etree import ElementTree
    from .forcefield_identity import IdentityForceField
    native = app.AmberPrmtopFile(str(folder / 'ligand.prmtop'))
    positions = app.AmberInpcrdFile(str(folder / 'ligand.inpcrd')).positions
    templates = ElementTree.parse(xml_path).findall('Residues/Residue')
    residues = list(native.topology.residues())
    if len(templates) != 1 or len(residues) != 1:
        raise ValueError('Ligand conversion validation requires one native residue and one XML template.')
    # These names are assigned from the checked native atom order above. Graph
    # matching alone can exchange symmetric atoms and their improper terms.
    atom_map = [{'prepared_name': atom.name,
                 'template_atom_name': f'{atom.element.symbol}{atom.index + 1}',
                 'native_index': atom.index} for atom in native.topology.atoms()]
    converted = IdentityForceField(str(xml_path), ligand_atom_maps=[{
        'residue_key': residue_key(residues[0]), 'template_name': templates[0].get('name'),
        'atoms': atom_map}])
    systems = [native.createSystem(nonbondedMethod=app.NoCutoff, constraints=None),
               converted.createSystem(native.topology, nonbondedMethod=app.NoCutoff, constraints=None)]
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    report = {'method': 'Native Amber versus converted XML OpenMM Reference energy/force comparison at bound pose plus two deterministic 0.003 Å perturbations.',
              'energy_tolerance_kj_mol': 1e-4, 'force_tolerance_kj_mol_nm': 1e-3,
              'atom_mapping': {'method': 'Verified native atom names, elements and bonds; symmetric graph permutations are not used.',
                               'atoms': atom_map}, 'conformations': []}
    compatibility = folder / 'nonbonded-scale-compatibility.json'
    if compatibility.is_file():
        report['unused_global_14_scale_compatibility'] = json.loads(compatibility.read_text())
    for index in range(3):
        displacement = 0 if index == 0 else np.random.default_rng(2026 + index).normal(0, .0003, xyz.shape)
        values = []
        for system in systems:
            integrator = mm.VerletIntegrator(.001)
            context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
            context.setPositions(xyz + displacement)
            state = context.getState(getEnergy=True, getForces=True)
            values.append((state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                           np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))))
            del context, integrator
        entry = {'pose': 'bound' if index == 0 else f'perturbed-{index}', 'native_energy_kj_mol': values[0][0],
                 'xml_energy_kj_mol': values[1][0], 'energy_delta_kj_mol': values[1][0] - values[0][0],
                 'maximum_force_delta_kj_mol_nm': float(np.max(np.abs(values[1][1] - values[0][1])))}
        report['conformations'].append(entry)
        if not np.isfinite([*values[0][1].ravel(), *values[1][1].ravel(), values[0][0], values[1][0]]).all() or abs(entry['energy_delta_kj_mol']) > 1e-4 or entry['maximum_force_delta_kj_mol_nm'] > 1e-3:
            storage.atomic_json(folder / 'conversion-validation.json', dict(report, passed=False))
            raise ValueError('Ligand force-field conversion did not reproduce native Amber energies/forces. Review conversion-validation.json; no approximate replacement was accepted.')
    report['passed'] = True
    storage.atomic_json(folder / 'conversion-validation.json', report)
    return report


def _amber_versions():
    folder = ambertools_home()
    records = list((folder / 'conda-meta').glob('ambertools-*.json'))
    record = json.loads(records[0].read_text()) if len(records) == 1 else {}
    gaff = folder / 'dat' / 'leap' / 'parm' / 'gaff2.dat'
    return {'ambertools_version': record.get('version', 'unknown external installation'),
            'ambertools_build': record.get('build'), 'gaff2_title': gaff.read_text().splitlines()[0],
            'gaff2_sha256': hashlib.sha256(gaff.read_bytes()).hexdigest(),
            'bcc_parameters_sha256': hashlib.sha256((folder / 'dat' / 'antechamber' / 'BCCPARM.DAT').read_bytes()).hexdigest()}


def _topology(model, mol, names):
    top = app.Topology()
    chain = top.addChain(model.residue.chain.id)
    residue = top.addResidue(model.residue.name, chain, model.residue.id, model.residue.insertionCode)
    atoms = [top.addAtom(name, app.element.Element.getBySymbol(atom.GetSymbol()), residue) for atom, name in zip(mol.GetAtoms(), names)]
    for bond in mol.GetBonds():
        top.addBond(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()])
    return top, np.asarray(mol.GetConformer().GetPositions()) * unit.angstrom


def prepare_ligands(dataset_id, ph, seed, folder, overrides=None, on_progress=None, check_cancel=None, actions=None, environment_indices=None):
    """Prepare all noncovalent organic residues or fail without dropping any."""
    models = _models(dataset_id, ph, overrides, actions, execute_repair=True, seed=seed, on_progress=on_progress, check_cancel=check_cancel, environment_indices=environment_indices)
    errors = [f"{model.description['key']}: {model.description['error']}" for model in models if model.description['error']]
    if errors:
        raise ValueError(' '.join(errors))
    cache, cache_sources, results = {}, {}, []
    for model in models:
        if model.description.get("removed"):
            continue
        if check_cancel:
            check_cancel()
        mol, hydrogen_method = _add_hydrogens(model)
        ordered_graph = {"smiles": model.description['selected_smiles'],
                         "atoms": [(a.GetSymbol(), a.GetFormalCharge(), int(a.GetChiralTag())) for a in mol.GetAtoms()],
                         "bonds": [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble(), int(b.GetStereo())) for b in mol.GetBonds()]}
        fingerprint = hashlib.sha256(json.dumps(ordered_graph, sort_keys=True).encode()).hexdigest()[:12]
        namespace = 'DML_' + fingerprint
        parameter_folder = Path(folder) / ('ligand-' + fingerprint)
        if fingerprint not in cache:
            cache[fingerprint] = _parameterize(mol, parameter_folder, namespace, on_progress, check_cancel)
            cache_sources[fingerprint] = residue_key(model.residue)
        xml_path, generated_names, charges = cache[fingerprint]
        # Identical chemical graphs can be differently ordered across residues.
        # Only reuse parameterization when the atom order also agrees.
        signature = [(a.GetSymbol(), a.GetFormalCharge()) for a in mol.GetAtoms()]
        identity_file = parameter_folder / 'identity.json'
        bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()) for b in mol.GetBonds()]
        identity = {'atoms': signature, 'bonds': bonds}
        if identity_file.is_file() and json.loads(identity_file.read_text()) != json.loads(json.dumps(identity)):
            raise ValueError("Copies of one ligand use different atom orderings. Reimport with consistent ligand atom names/order before preparing this complex.")
        storage.atomic_json(identity_file, identity)
        top, positions = _topology(model, mol, generated_names)
        # Restore heavy names for the viewer and retain an explicit complete
        # map to the native parameter atoms, including the added hydrogens.
        output_names = model.names + [f"H{i + 1}" for i in range(mol.GetNumAtoms() - len(model.names))]
        if len(set(output_names)) != len(output_names):
            raise ValueError('Prepared ligand atom names are not unique; an exact parameter mapping cannot be retained.')
        for atom, name in zip(top.atoms(), output_names):
            atom.name = name
        atom_map = [{'original_index': index, 'original_name': name, 'ligand_index': i, 'prepared_name': name}
                    for i, (index, name) in enumerate(zip(model.atom_indices, model.names)) if index is not None]
        parameter_atom_map = [{'native_index': i, 'prepared_name': name, 'template_atom_name': template_name}
                              for i, (name, template_name) in enumerate(zip(output_names, generated_names))]
        provenance = {**model.description, 'hydrogen_placement': hydrogen_method,
                      'hydrogens_added': mol.GetNumAtoms() - len(model.names), 'forcefield': 'GAFF2',
                      'charge_method': 'AmberTools antechamber AM1-BCC', 'charge_sum_e': sum(charges),
                      'charges_e': charges, 'charge_rounding': json.loads((parameter_folder / 'charge-rounding.json').read_text()), 'template_name': namespace, 'ffxml_sha256': hashlib.sha256(xml_path.read_bytes()).hexdigest(),
                      'heavy_coordinate_max_displacement_angstrom': 0.0, 'seed': seed,
                      'versions': {key: importlib.metadata.version(key) for key in ('rdkit', 'parmed', 'dimorphite-dl', 'gemmi')},
                      'artifacts_directory': parameter_folder.name, 'atom_map': atom_map,
                      'parameter_atom_map': parameter_atom_map,
                      'parameterization_source_residue_key': list(cache_sources[fingerprint]),
                      'parameter_reuse_note': 'Identical ordered chemical states share parameters from the first bound copy; each copy retains its own heavy-atom coordinates.',
                      'native_runtime': _amber_versions(), 'implementation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'conversion_validation': json.loads((parameter_folder / 'conversion-validation.json').read_text())}
        if residue_key(model.residue) == cache_sources[fingerprint]:
            storage.atomic_json(parameter_folder / 'provenance.json', provenance)
        instance_id = hashlib.sha256(model.description['key'].encode()).hexdigest()[:12]
        storage.atomic_json(parameter_folder / f'residue-{instance_id}.json', provenance)
        results.append(PreparedLigand(top, positions, xml_path, provenance, residue_key(model.residue), atom_map))
    return results
