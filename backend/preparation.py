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

STANDARD_PROTEINS = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
MAX_LOOP_LENGTH = 6
MAX_REBUILT_RESIDUES = 12
MAX_PROTEIN_ATOMS = 20_000
PREPARATION_DEFAULTS = {"name": "Protein preparation", "ph": 7.0, "add_missing_atoms": True, "build_missing_residues": False, "optimize_sidechains": True, "remove_waters": True, "remove_heterogens": False, "seed": 2026}


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
                    source = PDBFixer(filename=str(path))
                    if source.sequences:
                        return path
                except Exception:
                    continue
    return None


def current_fixer(dataset_id: str, sequence_source: Path | None = None):
    from pdbfixer import PDBFixer
    from openmm import Platform
    fixer = PDBFixer(filename=str(exact_input_path(dataset_id)), platform=Platform.getPlatformByName("CPU"))
    source_path = sequence_source if sequence_source is not None else find_sequence_source(dataset_id)
    if source_path:
        original = PDBFixer(filename=str(source_path))
        chain_ids = {chain.id for chain in fixer.topology.chains()}
        label_to_author = {}
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
            chain_id = original_to_canonical.get(author, author)
            if chain_id in chain_ids:
                sequence = copy.copy(sequence)
                sequence.chainId = chain_id
                mapped.append(sequence)
        fixer.sequences = mapped
    return fixer, source_path


def backbone_gaps(topology, positions) -> list[dict]:
    """Numbering anomalies are observations; only geometry identifies broken links."""
    from openmm import unit
    import numpy as np
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    result = []
    for chain in topology.chains():
        residues = [residue for residue in chain.residues() if residue.name in STANDARD_PROTEINS]
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
                result.append({"chain": chain.id, "after": previous.id, "before": following.id, "message": message, "backbone_distance_angstrom": distance, "structural_break": broken})
    return result


def inspect_preparation(dataset_id: str) -> dict:
    metadata = storage.get_dataset(dataset_id)
    atoms = metadata["atoms"]
    warnings, blockers = [], []
    protein_atoms = sum(atom["category"] == "protein" for atom in atoms)
    heterogens = sorted({f"{atom['residue']} {atom['chain']}:{atom['resid']}" for atom in atoms if atom["category"] not in {"protein", "water"}})
    if not protein_atoms:
        blockers.append("Protein preparation requires a standard amino-acid protein. Small molecules and nucleic acids remain viewable but are not prepared by this workflow.")
    if protein_atoms > MAX_PROTEIN_ATOMS:
        blockers.append(f"Preparation is capped at {MAX_PROTEIN_ATOMS:,} protein atoms on this local CPU workflow.")
    if heterogens:
        warnings.append("Nonprotein molecules are present. This protein-only preparation cannot parameterize them; explicitly remove heterogens or prepare the complex externally.")
    fixer, sequence_source = current_fixer(dataset_id)
    protein_residues = {(atom["chain"], str(atom["resid"])) for atom in atoms if atom["category"] == "protein"}
    unsupported = sorted({residue.name for residue in fixer.topology.residues() if residue.name not in STANDARD_PROTEINS and (residue.chain.id, residue.id) in protein_residues})
    if unsupported:
        blockers.append("Unsupported protein residue templates: " + ", ".join(unsupported) + ". Residue mutation is not performed automatically.")
    try:
        fixer.findMissingResidues()
    except (ValueError, IndexError) as exc:
        fixer.missingResidues = {}
        warnings.append(f"Sequence mapping could not be established ({exc}); missing residue identities will not be inferred.")
    missing_residues = []
    chains = list(fixer.topology.chains())
    for (chain_index, position), names in fixer.missingResidues.items():
        terminal = position in {0, len(list(chains[chain_index].residues()))}
        buildable = not terminal and len(names) <= MAX_LOOP_LENGTH and all(name in STANDARD_PROTEINS for name in names)
        missing_residues.append({"chain": chains[chain_index].id, "position": position, "residues": names, "count": len(names), "terminal": terminal, "buildable": buildable, "chain_index": chain_index})
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
        warnings.append("Missing sequence-supported residues are listed. Rebuilding internal short gaps is opt-in and produces low-confidence coordinates, not an experimentally established loop.")
    if any(gap["structural_break"] for gap in gaps):
        warnings.append("A long backbone connection is present. Unresolved internal gaps must be repaired before force-field preparation or simulation; a stretched inferred peptide bond is not accepted.")
    if any(entry["terminal"] for entry in missing_residues):
        warnings.append("Missing terminal sequence is reported but not rebuilt. The observed structure is prepared as a truncated chain with terminal groups.")
    warnings.append("pH assignment uses OpenMM template/heuristic protonation and histidine tautomer rules, not a pKa calculation or constant-pH simulation.")
    return {"dataset_id": dataset_id, "protein_atoms": protein_atoms, "hydrogen_atoms": sum(atom["element"] == "H" for atom in atoms), "water_atoms": sum(atom["category"] == "water" for atom in atoms), "heterogen_residues": heterogens, "can_prepare": not blockers, "has_sequence": bool(fixer.sequences), "missing_atoms": missing_atoms, "missing_residues": missing_residues, "gaps": gaps, "warnings": warnings, "blockers": blockers, "sequence_source": str(sequence_source) if sequence_source else None, "sequence_source_sha256": hashlib.sha256(sequence_source.read_bytes()).hexdigest() if sequence_source else None}


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


def submit_preparation(settings: dict) -> dict:
    settings = _validated(settings)
    inspection = inspect_preparation(settings["dataset_id"])
    if inspection["blockers"]:
        raise ValueError(" ".join(inspection["blockers"]))
    if inspection["heterogen_residues"] and not settings["remove_heterogens"]:
        raise ValueError("Nonprotein molecules are present. Select explicit heterogen removal to prepare the protein alone, or prepare the complex externally; no molecules were removed.")
    internal = [entry for entry in inspection["missing_residues"] if not entry["terminal"]]
    if internal and not settings["build_missing_residues"]:
        raise ValueError("Unresolved internal sequence gaps would create artificial peptide connections. Enable short missing-loop building, or repair the structure externally before protein preparation.")
    if settings["build_missing_residues"] and (any(not entry["buildable"] for entry in internal) or sum(entry["count"] for entry in internal) > MAX_REBUILT_RESIDUES):
        raise ValueError(f"Only internal sequence-supported loops of up to {MAX_LOOP_LENGTH} residues each and {MAX_REBUILT_RESIDUES} total can be built in this local workflow. Repair larger gaps externally.")
    # Numbering alone is not missing-sequence evidence; long links with no mapped gap cannot be invented.
    mapped_pairs = set()
    fixer, _ = current_fixer(settings["dataset_id"])
    chains = list(fixer.topology.chains())
    for entry in internal:
        residues = list(chains[entry["chain_index"]].residues())
        mapped_pairs.add((entry["chain"], residues[entry["position"] - 1].id, residues[entry["position"]].id))
    if any(gap["structural_break"] and (gap["chain"], gap["after"], gap["before"]) not in mapped_pairs for gap in inspection["gaps"]):
        raise ValueError("An unresolved long backbone connection has no supported missing sequence. Supply the original PDB/mmCIF sequence records or repair the gap externally; residue identities will not be guessed.")
    if inspection["missing_atoms"] and not settings["add_missing_atoms"]:
        raise ValueError("Missing heavy/terminal atoms prevent force-field preparation. Enable missing-atom repair or supply a complete structure.")
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
