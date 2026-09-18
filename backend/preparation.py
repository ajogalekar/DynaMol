"""Inspection and background submission for transparent standard-protein preparation."""
from __future__ import annotations

import copy
import datetime as dt
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from . import config, jobs, storage
from .residue_identity import STANDARD_PROTEINS, protein_residue_keys, residue_key
from .modified_residues import register_topology_definitions, register_fixer_templates, inspect_modified, SUPPORTED_MODIFIED

MAX_LOOP_LENGTH = 12
# Resource budget across independent loops, including repeated crystal monomers.
# These limits bound local work; they do not establish conformational accuracy.
MAX_REBUILT_RESIDUES = 96
SHORT_LOOP_LENGTH = 6
MAX_PROTEIN_ATOMS = 20_000
SUPPORTED_CAPS = frozenset({"ACE", "NME"})
PREPARATION_DEFAULTS = {"name": "Protein preparation", "ph": 7.0, "add_missing_atoms": True, "build_missing_residues": False, "optimize_sidechains": True, "remove_waters": True, "remove_heterogens": False, "ligand_overrides": {}, "seed": 2026}


def loop_policy() -> dict:
    return {"max_gap_residues": MAX_LOOP_LENGTH, "max_total_residues": MAX_REBUILT_RESIDUES,
            "short_gap_residues": SHORT_LOOP_LENGTH}


def validate_loop_selection(internal: list[dict], enabled: bool) -> None:
    """Shared submission/worker limit; only explicit, sequence-supported repair.

    Definitive "cannot repair here" verdicts come first, so a user is never told to
    enable loop building for a gap that building could never fix. Those gaps must be
    resolved outside DynaMol (a complete experimental or predicted model).
    """
    for entry in internal:
        names = entry["residues"]
        if any(name not in STANDARD_PROTEINS for name in names):
            raise ValueError(unsupported_missing_residue_message(entry["chain"], names))
        if len(names) > MAX_LOOP_LENGTH:
            raise ValueError(f"Chain {entry['chain']} has a {len(names)}-residue internal gap, beyond the {MAX_LOOP_LENGTH}-residue loop DynaMol can model. This gap can't be repaired here: supply a complete structure (experimental or predicted) with the loop resolved, or rebuild it externally and reload. If another chain is intact, 'Use one monomer' can prepare that one instead.")
    total = sum(len(entry["residues"]) for entry in internal)
    if total > MAX_REBUILT_RESIDUES:
        raise ValueError(f"This structure needs {total} modeled internal residues, beyond DynaMol's {MAX_REBUILT_RESIDUES}-residue loop-building budget. Prepare fewer chains with 'Use one monomer', supply a complete structure, or repair the loops externally, then reload.")
    if internal and not enabled:
        raise ValueError(f"This structure has {len(internal)} internal missing loop(s). Turn on 'Build supported missing loops / residues' to model them (provisional starting coordinates), or use 'Use one monomer' to prepare an intact chain. An internal gap left unbuilt would form an artificial stretched peptide bond, so it can't be simulated as-is.")


def _json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def exact_input_path(dataset_id: str) -> Path:
    """Prefer the exact OpenMM output over display topology normalization."""
    folder = storage.dataset_dir(dataset_id)
    prepared = folder / "prepared.pdb"
    return prepared if prepared.is_file() else folder / "topology.pdb"


def find_sequence_source(dataset_id: str) -> Path | None:
    """Follow retained local source ancestry; never fetch or infer sequence."""
    visited = set()
    pending = [dataset_id]
    while pending and len(visited) < 20:
        current = pending.pop(0)
        if current in visited:
            continue
        visited.add(current)
        try:
            folder = storage.dataset_dir(current)
        except ValueError:
            continue
        meta = _json(folder / "metadata.json")
        provenance = _json(folder / "provenance.json")
        candidates = [folder / "sequence-source.pdb", folder / "sequence-source.cif", *sorted(folder.glob("source.*")), *sorted((folder / "originals").glob("topology.*"))]
        for record in (meta, provenance, provenance.get("starting_structure", {})):
            if not isinstance(record, dict):
                continue
            value = record.get("source_file")
            if value:
                path = Path(value)
                if not path.is_absolute():
                    path = folder / path
                if path.resolve().is_relative_to(config.DATA_ROOT.resolve()):
                    candidates.append(path)
            pdb_id = record.get("pdb_id")
            if isinstance(pdb_id, str) and re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id):
                candidates.append(config.DATA_ROOT / "sources" / f"{pdb_id.upper()}.pdb")
            for key in ("parent_dataset_id", "source_dataset_id"):
                if isinstance(record.get(key), str):
                    pending.append(record[key])
        for state in (meta.get("preparation", {}), meta.get("solvation", {}), provenance.get("preparation_state", {})):
            if isinstance(state, dict) and isinstance(state.get("parent_dataset_id"), str):
                pending.append(state["parent_dataset_id"])
        for path in candidates:
            if path.is_file() and path.suffix.lower() in {".pdb", ".cif", ".mmcif", ".pdbx"}:
                # Keep scanning ancestors if a canonical file happens to lack sequence evidence.
                try:
                    from pdbfixer import PDBFixer
                    register_topology_definitions()
                    source = PDBFixer(filename=str(path))
                    if source.sequences:
                        return path
                except Exception:
                    continue
    return None


def current_fixer(dataset_id: str, sequence_source: Path | None = None):
    from pdbfixer import PDBFixer
    from openmm import Platform
    register_topology_definitions()
    fixer = PDBFixer(filename=str(exact_input_path(dataset_id)), platform=Platform.getPlatformByName("CPU"))
    fixer.protonation_aliases = canonicalize_protonation_aliases(fixer.topology)
    register_fixer_templates(fixer)
    source_path = sequence_source if sequence_source is not None else find_sequence_source(dataset_id)
    if source_path:
        original = PDBFixer(filename=str(source_path))
        register_fixer_templates(original)
        # Retained ligands may use the ID of an excluded protein chain. Only
        # protein-containing chains may inherit polymer sequence records.
        protein_keys = protein_residue_keys(dataset_id, fixer.topology, fixer.positions)
        chain_ids = {residue.chain.id for residue in fixer.topology.residues() if residue_key(residue) in protein_keys}
        label_to_author, label_to_original = {}, {}
        if source_path.suffix.lower() in {".cif", ".mmcif", ".pdbx"}:
            try:
                from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
                containers = []
                with source_path.open() as handle:
                    PdbxReader(handle).read(containers)
                atom_site = containers[0].getObj("atom_site")
                for row in range(atom_site.getRowCount()):
                    label = atom_site.getValue("label_asym_id", row)
                    author = atom_site.getValue("auth_asym_id", row)
                    if author not in {".", "?", ""}:
                        label_to_author[label] = author
            except (KeyError, ValueError, AttributeError, IndexError):
                pass
            from .sequence_evidence import source_cif_chain_ids
            label_to_original = source_cif_chain_ids(source_path, original.topology)
        original_to_canonical = {}
        pending, visited = [dataset_id], set()
        while pending and len(visited) < 20:
            parent = pending.pop(0)
            if parent in visited:
                continue
            visited.add(parent)
            parent_folder = storage.dataset_dir(parent)
            for record in (_json(parent_folder / "metadata.json"), _json(parent_folder / "provenance.json")):
                for mapping in record.get("chain_id_mapping", []):
                    original_to_canonical[mapping["original"]] = mapping["canonical"]
                for key in ("parent_dataset_id", "source_dataset_id"):
                    if isinstance(record.get(key), str):
                        pending.append(record[key])
                for key in ("preparation", "solvation"):
                    state = record.get(key, {})
                    if isinstance(state, dict) and isinstance(state.get("parent_dataset_id"), str):
                        pending.append(state["parent_dataset_id"])
        mapped = []
        for sequence in original.sequences:
            author = label_to_author.get(sequence.chainId, sequence.chainId)
            original_id = label_to_original.get(sequence.chainId, author)
            chain_id = original_to_canonical.get(original_id, original_id)
            if chain_id in chain_ids:
                sequence = copy.copy(sequence)
                sequence.chainId = chain_id
                mapped.append(sequence)
        fixer.sequences = mapped
        if source_path.suffix.lower() in {".cif", ".mmcif", ".pdbx"}:
            from .sequence_evidence import load_scheme, recover_observed_insertion_codes
            try:
                fixer.sequence_scheme = load_scheme(source_path, label_to_author, original_to_canonical, chain_ids, mapped, label_to_original)
                fixer.recovered_insertion_codes = recover_observed_insertion_codes(fixer)
            except (ValueError, KeyError, TypeError) as exc:
                fixer.sequence_scheme_error = str(exc)
    return fixer, source_path


def canonicalize_protonation_aliases(topology):
    """Use parent naming for two protonation aliases OpenMM does not normalize.

    LYN/LYS and CYM/CYS have the same heavy-atom graph. Preparation explicitly
    removes existing H and reassigns pH-dependent states; source files remain
    intact and this naming normalization is recorded. No modified sidechain
    chemistry (MSE, ALY, phosphorylated residues, etc.) is substituted.
    """
    records = []
    for residue in topology.residues():
        if residue.name not in {"LYN", "CYM"}:
            continue
        original = residue.name
        parent = {"LYN": "LYS", "CYM": "CYS"}[original]
        records.append({"chain": residue.chain.id, "resid": residue.id,
                        "insertion_code": (residue.insertionCode or "").strip(),
                        "input_alias": original, "template_name": parent})
        residue.name = parent
    if records:
        topology.createStandardBonds()
    return records


def find_missing_residues_preserving_identity(fixer):
    """Use pinned PDBFixer's sequence alignment without parent substitutions.

    PDBFixer 1.12 normally rewrites missing TPO to THR (and other modifications
    to parents). A private copy of the function's global namespace suppresses
    only that dictionary lookup, without mutating library globals shared by
    concurrent inspection requests. Existing atom/residue names are untouched.
    """
    if getattr(fixer, "sequence_scheme_error", None):
        raise ValueError(fixer.sequence_scheme_error)
    if getattr(fixer, "sequence_scheme", None):
        from .sequence_evidence import missing_from_scheme
        missing_from_scheme(fixer)
        return
    from types import FunctionType
    original = type(fixer).findMissingResidues
    preserving = FunctionType(original.__code__, {**original.__globals__, "substitutions": {}},
                              original.__name__, original.__defaults__, original.__closure__)
    preserving(fixer)


def unsupported_missing_residue_message(chain, names):
    unsupported = sorted({name for name in names if name not in STANDARD_PROTEINS})
    return (f"Missing sequence-supported residues in chain {chain} include modified or unnatural amino acids: "
            + ", ".join(unsupported)
            + ". Their original identities are retained; the local loop builder cannot reconstruct these residues. "
              "Supply coordinates for the exact modified residues or a validated complete model; no parent-residue substitution is performed.")


def backbone_gaps(topology, positions) -> list[dict]:
    """Numbering anomalies are observations; only geometry identifies broken links."""
    from openmm import unit
    import numpy as np
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    result = []
    protein_keys = protein_residue_keys(None, topology, positions)
    for chain in topology.chains():
        residues = [residue for residue in chain.residues() if residue_key(residue) in protein_keys]
        for previous, following in zip(residues, residues[1:]):
            atoms_before = {atom.name: atom for atom in previous.atoms()}
            atoms_after = {atom.name: atom for atom in following.atoms()}
            distance = None
            if "C" in atoms_before and "N" in atoms_after:
                distance = float(np.linalg.norm(xyz[atoms_before["C"].index] - xyz[atoms_after["N"].index]) * 10)
            try:
                numbering_gap = int(following.id) - int(previous.id) > 1
            except ValueError:
                numbering_gap = False
            if numbering_gap or (distance is not None and distance > 2.2):
                broken = distance is not None and distance > 2.2
                message = f"Backbone C–N separation is {distance:.2f} Å; this connection requires repair before force-field preparation." if broken else "Residue numbering skips values. This alone does not establish missing sequence or missing structure."
                result.append({"chain": chain.id, "after": previous.id, "before": following.id, "after_insertion_code": (previous.insertionCode or "").strip(), "before_insertion_code": (following.insertionCode or "").strip(), "after_index": previous.index, "before_index": following.index, "message": message, "backbone_distance_angstrom": distance, "structural_break": broken})
    return result


def inspect_preparation(dataset_id: str, ph: float = 7.0, ligand_overrides: dict | None = None, ligand_actions: dict | None = None) -> dict:
    if not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError("Choose a pH between 0 and 14.")
    metadata = storage.get_dataset(dataset_id)
    atoms = metadata["atoms"]
    warnings, blockers = [], []
    fixer, sequence_source = current_fixer(dataset_id)
    if getattr(fixer, "recovered_insertion_codes", None):
        recovered = fixer.recovered_insertion_codes
        warnings.append("Restored omitted insertion codes from unique source residue number/name matches and verified chain order: " + ", ".join(f"{r['chain']}:{r['resid']}{r['insertion_code']} {r['residue']}" for r in recovered) + ". Observed coordinates and original files are unchanged; the corrected identities are retained during preparation.")
    if fixer.protonation_aliases:
        warnings.append("Input LYN/CYM names identify protonation states of lysine/cysteine. Their heavy-atom graphs are retained under LYS/CYS template names; pressing Prep explicitly reassigns hydrogens at the selected pH and records the resulting state.")
    protein_keys = protein_residue_keys(dataset_id, fixer.topology, fixer.positions)
    protein_atoms = sum(residue_key(atom.residue) in protein_keys for atom in fixer.topology.atoms())
    heterogens = sorted({f"{residue.name} {residue.chain.id}:{residue.id}" for residue in fixer.topology.residues() if residue_key(residue) not in protein_keys and residue.name not in storage.WATERS})
    if not protein_atoms:
        blockers.append("Protein preparation requires an amino-acid polymer. Small molecules and nucleic acids remain viewable but are not prepared by this workflow.")
    if protein_atoms > MAX_PROTEIN_ATOMS:
        blockers.append(f"Preparation is capped at {MAX_PROTEIN_ATOMS:,} protein atoms on this local CPU workflow.")
    from .ligands import inspect_ligands, ligand_runtime_status
    from .complex_topology import metal_environment
    from .ions import inspect_ions
    ligands = inspect_ligands(dataset_id, ph, ligand_overrides, actions=ligand_actions)
    ligand_errors = [ligand["error"] if ligand["error"].startswith(ligand["key"] + ":") else f"{ligand['key']}: {ligand['error']}" for ligand in ligands if ligand.get("error")]
    coordination = metal_environment(fixer.topology, fixer.positions)["report"]
    ions = inspect_ions(fixer.topology)
    if ligands:
        warnings.append("Ligands are retained with their bound heavy-atom coordinates. Preparation assigns a recorded fixed protonation state and GAFF2 / AM1-BCC parameters for explicit-water OpenMM dynamics.")
    if coordination["contacts"]:
        warnings.append("Observed coordinating waters and protein donor sidechains are retained. Standard nonbonded ion parameters do not establish metal-coordination accuracy.")
    if any(ion["supported"] for ion in ions):
        warnings.append("Retained ions use declared monatomic charge states and require explicit TIP3P water; coordination energetics and alternate oxidation states are not modeled.")
    modified = inspect_modified(fixer.topology, ph)
    warnings.extend(modified["warnings"])
    blockers.extend(modified["blockers"])
    modified_residues = list(modified["residues"])
    for residue in fixer.topology.residues():
        if residue_key(residue) not in protein_keys or residue.name in STANDARD_PROTEINS | SUPPORTED_CAPS or residue.name in SUPPORTED_MODIFIED:
            continue
        insertion = (residue.insertionCode or '').strip()
        message = f"{residue.name} {residue.chain.id}:{residue.id}{insertion} is a modified or unnatural protein residue with no supported covalent amino-acid template. Its identity and atoms are preserved, including when ligand removal is selected. Supply a force-field template for this exact residue and state; automatic mutation or free-ligand GAFF treatment is not performed."
        blockers.append(message)
        modified_residues.append({"chain": residue.chain.id, "resid": residue.id, "insertion_code": insertion, "residue": residue.name, "supported": False, "error": message})
    for first, second in fixer.topology.bonds():
        if first.residue == second.residue or any(residue_key(atom.residue) not in protein_keys for atom in (first, second)):
            continue
        peptide = {first.name, second.name} == {"C", "N"}
        disulfide = all(atom.name == "SG" and atom.residue.name == "CYS" for atom in (first, second))
        if not peptide and not disulfide:
            blockers.append(f"Unsupported covalent protein crosslink: {first.residue.name} {first.residue.chain.id}:{first.residue.id}/{first.name}–{second.residue.name} {second.residue.chain.id}:{second.residue.id}/{second.name}. Standard peptide and cysteine disulfide links are supported; this link needs a matching specialized template. No bond or atom was removed.")
    try:
        find_missing_residues_preserving_identity(fixer)
    except (ValueError, IndexError) as exc:
        fixer.missingResidues = {}
        blockers.append(f"Sequence mapping could not be established ({exc}); missing residue identities will not be inferred. Supply an unambiguous original sequence mapping or a validated complete model.")
    missing_residues = []
    chains = list(fixer.topology.chains())
    for (chain_index, position), names in fixer.missingResidues.items():
        terminal = position in {0, len(list(chains[chain_index].residues()))}
        buildable = not terminal and len(names) <= MAX_LOOP_LENGTH and all(name in STANDARD_PROTEINS for name in names)
        modeling = "terminal" if terminal else "unsupported" if not buildable else "extended" if len(names) > SHORT_LOOP_LENGTH else "short"
        entry = {"chain": chains[chain_index].id, "position": position, "residues": names, "count": len(names), "terminal": terminal, "buildable": buildable, "chain_index": chain_index, "modeling": modeling}
        if (chain_index, position) in getattr(fixer, "missingResidueIdentities", {}):
            entry["source_identities"] = fixer.missingResidueIdentities[(chain_index, position)]
        if any(name not in STANDARD_PROTEINS for name in names):
            entry["reason"] = unsupported_missing_residue_message(chains[chain_index].id, names)
            if terminal:
                warnings.append(entry["reason"] + " This terminal sequence remains omitted from the prepared truncated chain.")
            else:
                blockers.append(entry["reason"])
        missing_residues.append(entry)
    fixer.findMissingAtoms()
    missing_atoms = []
    affected = set(fixer.missingAtoms) | set(fixer.missingTerminals)
    for residue in sorted(affected, key=lambda r: r.index):
        names = [atom.name for atom in fixer.missingAtoms.get(residue, [])] + fixer.missingTerminals.get(residue, [])
        missing_atoms.append({"chain": residue.chain.id, "resid": residue.id, "residue": residue.name, "atoms": sorted(set(names))})
    gaps = backbone_gaps(fixer.topology, fixer.positions)
    if not fixer.sequences:
        warnings.append("No retained SEQRES/mmCIF polymer sequence is available. Residue numbering gaps cannot identify which amino acids are absent, so loop building is unavailable.")
    if missing_residues:
        warnings.append("Missing sequence-supported residues are listed. Rebuilding internal gaps is opt-in and produces provisional coordinates; longer loops have greater conformational uncertainty. Geometry checks do not establish the native loop conformation.")
    if any(gap["structural_break"] for gap in gaps):
        warnings.append("A long backbone connection is present. Unresolved internal gaps must be repaired before force-field preparation or simulation; a stretched inferred peptide bond is not accepted.")
    if any(entry["terminal"] for entry in missing_residues):
        warnings.append("Missing terminal sequence is reported but not rebuilt. The observed structure is prepared as a truncated chain with terminal groups.")
    warnings.append("pH assignment uses OpenMM template/heuristic protonation and histidine tautomer rules, not a pKa calculation or constant-pH simulation.")
    return {"dataset_id": dataset_id, "protein_atoms": protein_atoms, "hydrogen_atoms": sum(atom["element"] == "H" for atom in atoms), "water_atoms": sum(atom["category"] == "water" for atom in atoms), "heterogen_residues": heterogens, "can_prepare": not blockers, "has_sequence": bool(fixer.sequences), "loop_policy": loop_policy(), "missing_atoms": missing_atoms, "missing_residues": missing_residues, "gaps": gaps, "warnings": warnings, "blockers": blockers, "ligands": ligands, "modified_residues": modified_residues, "ligand_errors": ligand_errors, "ligand_runtime": ligand_runtime_status() if ligands else None, "metal_environment": coordination, "ions": ions, "sequence_source": str(sequence_source) if sequence_source else None, "sequence_source_sha256": hashlib.sha256(sequence_source.read_bytes()).hexdigest() if sequence_source else None}


def _validated(settings: dict) -> dict:
    settings = {**PREPARATION_DEFAULTS, **dict(settings)}
    storage.safe_id(settings.get("dataset_id", ""))
    settings["name"] = str(settings["name"] or PREPARATION_DEFAULTS["name"])[:100]
    ph = float(settings["ph"])
    if not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError("Choose a pH between 0 and 14.")
    settings["ph"] = ph
    seed = int(settings["seed"])
    if not 1 <= seed <= 2_147_483_646:
        raise ValueError("Seed must be an integer from 1 to 2,147,483,646.")
    settings["seed"] = seed
    overrides = settings.get("ligand_overrides", {})
    if not isinstance(overrides, dict) or len(overrides) > 100 or any(not isinstance(k, str) or len(k) > 100 or not isinstance(v, str) or not 1 <= len(v) <= 10000 for k, v in overrides.items()):
        raise ValueError("Ligand overrides must map residue identifiers to explicit-state SMILES strings.")
    settings["ligand_overrides"] = overrides
    actions = settings.get("ligand_actions", {})
    if not isinstance(actions, dict) or len(actions) > 100 or any(not isinstance(k, str) or not 1 <= len(k) <= 100 or v not in ("repair", "remove") for k, v in actions.items()):
        raise ValueError("Choose repair or remove for each exact ligand residue identifier.")
    settings["ligand_actions"] = actions
    for key in ("add_missing_atoms", "build_missing_residues", "optimize_sidechains", "remove_waters", "remove_heterogens"):
        if not isinstance(settings[key], bool):
            raise ValueError(f"{key} must be true or false.")
    return settings


def _submit(settings: dict, operation: str) -> dict:
    metadata = storage.get_dataset(settings["dataset_id"])
    with jobs._lock:
        if any(job["status"] in {"queued", "running", "cancelling"} for job in jobs.list_jobs()):
            raise ValueError("A simulation or preparation is already active. Wait for it or cancel it before starting another CPU job.")
        job_id = uuid.uuid4().hex[:16]
        folder = config.JOBS_DIR / job_id
        folder.mkdir()
        settings = {**settings, "operation": operation}
        stages = 7 if operation == "prepare" else 2
        job = {"id": job_id, "name": settings["name"], "engine": "preparation" if operation == "prepare" else "solvation", "status": "queued", "stage": "Starting preparation worker", "progress": 0, "completed_steps": 0, "total_steps": stages, "elapsed_seconds": 0, "logs": ["Background preparation uses real completed stages for progress, with two CPU threads by default."], "config": settings, "created_at": dt.datetime.now(dt.timezone.utc).isoformat()}
        storage.atomic_json(folder / "config.json", settings)
        storage.atomic_json(folder / "input-metadata.json", metadata)
        shutil.copy2(exact_input_path(settings["dataset_id"]), folder / "input.pdb")
        source = find_sequence_source(settings["dataset_id"])
        if source:
            shutil.copy2(source, folder / ("sequence-source" + source.suffix.lower()))
        storage.atomic_json(folder / "status.json", job)
        environment = os.environ.copy()
        environment.update(OPENMM_CPU_THREADS=str(config.CPU_THREADS), OMP_NUM_THREADS=str(config.CPU_THREADS), OPENBLAS_NUM_THREADS=str(config.CPU_THREADS), MKL_NUM_THREADS=str(config.CPU_THREADS), DYNAMOL_DATA_DIR=str(config.DATA_ROOT), PYTHONUNBUFFERED="1")
        with (folder / "worker.log").open("ab", buffering=0) as output:
            process = subprocess.Popen([sys.executable, "-m", "backend.preparation_worker", job_id], cwd=config.ROOT, env=environment, stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
        job["worker_pid"] = process.pid
        storage.atomic_json(folder / "status.json", job)
        return job


def validate_preparation(settings: dict) -> tuple[dict, dict]:
    """Return normalized settings and inspection without creating a job.

    The readiness pane and submission share this exact eligibility check.
    Chemical-reference retrieval can populate its existing CCD cache, but this
    function never creates or modifies a dataset or job.
    """
    settings = _validated(settings)
    inspection = inspect_preparation(settings["dataset_id"], settings["ph"], settings["ligand_overrides"], settings["ligand_actions"])
    if inspection["blockers"]:
        raise ValueError(" ".join(inspection["blockers"]))
    if not settings["remove_heterogens"]:
        if inspection["ligand_errors"]:
            raise ValueError("Ligand chemistry needs attention: " + " ".join(inspection["ligand_errors"]))
        if any(not ligand.get("removed") for ligand in inspection["ligands"]) and not inspection["ligand_runtime"]["available"]:
            raise ValueError(inspection["ligand_runtime"]["message"])
    settings["complex"] = bool(any(not ligand.get("removed") for ligand in inspection["ligands"]) and not settings["remove_heterogens"])
    settings["modified_residues"] = inspection["modified_residues"]
    internal = [entry for entry in inspection["missing_residues"] if not entry["terminal"]]
    validate_loop_selection(internal, settings["build_missing_residues"])
    if internal:
        from .loop_modeling import loop_runtime_status
        runtime = loop_runtime_status()
        if not runtime["available"]:
            raise ValueError(config.setup_guidance("Loop modeling", "`bash scripts/install_loop_tools.sh`", runtime.get("error", "")))
    # Numbering alone is not missing-sequence evidence; long links with no mapped gap cannot be invented.
    mapped_pairs = set()
    fixer, _ = current_fixer(settings["dataset_id"])
    chains = list(fixer.topology.chains())
    for entry in internal:
        residues = list(chains[entry["chain_index"]].residues())
        mapped_pairs.add((residues[entry["position"] - 1].index, residues[entry["position"]].index))
    if any(gap["structural_break"] and (gap["after_index"], gap["before_index"]) not in mapped_pairs for gap in inspection["gaps"]):
        raise ValueError("An unresolved long backbone connection has no supported missing sequence. Supply the original PDB/mmCIF sequence records or repair the gap externally; residue identities will not be guessed.")
    if inspection["missing_atoms"] and not settings["add_missing_atoms"]:
        raise ValueError("Missing heavy/terminal atoms prevent force-field preparation. Enable missing-atom repair or supply a complete structure.")
    return settings, inspection


def submit_preparation(settings: dict) -> dict:
    settings, _ = validate_preparation(settings)
    return _submit(settings, "prepare")


def submit_solvation(dataset_id: str, settings: dict) -> dict:
    settings = {**settings, "dataset_id": storage.safe_id(dataset_id), "name": "Explicit water preview"}
    padding = float(settings.get("padding_nm", 1))
    ph = float(settings.get("ph", 7))
    seed = int(settings.get("seed", 2026))
    if not math.isfinite(padding) or not 1 <= padding <= 3 or not math.isfinite(ph) or not 0 <= ph <= 14 or not 1 <= seed <= 2_147_483_646:
        raise ValueError("Use 1–3 nm padding, pH 0–14, and a positive supported seed.")
    metadata = storage.get_dataset(dataset_id)
    if not metadata.get("preparation"):
        raise ValueError("Prepare the protein first so explicit solvent preserves its chosen protonation state.")
    settings.update(padding_nm=padding, ph=ph, seed=seed)
    return _submit(settings, "solvate")
