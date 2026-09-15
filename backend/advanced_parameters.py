"""Immutable research parameter snapshots for advanced chemistry development.

These helpers preserve native Amber inputs and their independent atom/term
manifest. They do not register a model with the app or certify model accuracy.
An accepted mechanics comparison means the implementation preserved the native
parameters, not that a chemical state or parameterization is physically sound.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path, PurePosixPath
import tempfile

import numpy as np


REQUIRED_ROLES = {"native_topology", "coordinates", "expected_terms", "chemical_state", "source_structure"}
MAX_FILE_BYTES = 128 * 1024 * 1024
MAX_TOTAL_BYTES = 256 * 1024 * 1024


def _digest(data):
    return hashlib.sha256(data).hexdigest()


def _require_digest(value):
    if not isinstance(value,str) or len(value)!=64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError('An explicit lowercase SHA-256 snapshot manifest checksum is required')


def _json(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _read_regular(root, relative):
    path = PurePosixPath(relative) if isinstance(relative, str) else None
    if (path is None or path.is_absolute() or str(path) != relative
            or len(path.parts) != 1 or relative in {".", ".."}
            or "\\" in relative or ":" in relative):
        raise ValueError("Snapshot files must have simple relative filenames")
    candidate = root / relative
    if root.is_symlink() or candidate.is_symlink() or not candidate.is_file():
        raise ValueError(f"Snapshot file is missing, not regular, or a symlink: {relative}")
    if candidate.stat().st_size > MAX_FILE_BYTES:
        raise ValueError("Snapshot file exceeds the size limit")
    return candidate.read_bytes()


def read_snapshot(folder, *, expected_manifest_sha256=None):
    """Read once and verify all declared files before a native parser sees them."""
    root = Path(folder).absolute()
    raw = _read_regular(root, "snapshot.json")
    if expected_manifest_sha256 is not None:
        _require_digest(expected_manifest_sha256)
        if _digest(raw) != expected_manifest_sha256:
            raise ValueError("Advanced parameter snapshot manifest checksum mismatch")
    manifest = json.loads(raw)
    if (not isinstance(manifest,dict) or type(manifest.get("schema_version")) is not int or manifest.get("schema_version") != 1
            or manifest.get("stage") != "research_prototype"
            or manifest.get("simulation_ready") is not False):
        raise ValueError("This loader accepts research prototypes only; it cannot authorize app simulations")
    if not isinstance(manifest.get("model_id"), str) or not manifest["model_id"]:
        raise ValueError("The snapshot requires a model ID")
    if manifest.get("kind") not in {"covalent", "metal", "combined"}:
        raise ValueError("The snapshot requires an explicit model kind")
    files = manifest.get("files")
    if not isinstance(files, dict) or not REQUIRED_ROLES <= set(files) or len(files) > 32:
        raise ValueError("Advanced parameter snapshot is missing required source/model roles")
    contents, seen, total = {}, {"snapshot.json"}, len(raw)
    for role, entry in files.items():
        if not isinstance(entry, dict):
            raise ValueError("Malformed snapshot file entry")
        filename, expected = entry.get("path"), entry.get("sha256")
        if not isinstance(filename, str) or filename in seen:
            raise ValueError("Snapshot filenames must be unique")
        seen.add(filename)
        data = _read_regular(root, filename)
        if _digest(data) != expected:
            raise ValueError(f"Advanced parameter snapshot checksum mismatch: {role}")
        total += len(data)
        if total > MAX_TOTAL_BYTES:
            raise ValueError("Advanced parameter snapshot exceeds the total size limit")
        contents[role] = data
    return {"manifest": manifest, "manifest_bytes": raw,
            "manifest_sha256": _digest(raw), "contents": contents}


def create_snapshot(files, destination, *, model_id, kind, provenance):
    """Freeze supplied artifacts without modifying their content or readiness.

    `files` maps role to a source file. The destination must be fresh. Native
    calculations and scientific review take place before/outside this function.
    """
    destination = Path(destination).absolute()
    if destination.exists() or destination.is_symlink():
        raise ValueError("Refusing to replace an existing parameter snapshot")
    if not isinstance(files, dict) or not REQUIRED_ROLES <= set(files):
        raise ValueError("Provide native topology, coordinates, expected terms, chemical state and source structure")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".parameter-snapshot-", dir=destination.parent) as temporary:
        staged = Path(temporary) / "snapshot"
        staged.mkdir()
        manifest = {"schema_version": 1, "model_id": model_id, "kind": kind,
                    "stage": "research_prototype", "simulation_ready": False,
                    "provenance": provenance, "files": {}}
        for number, (role, value) in enumerate(sorted(files.items())):
            source = Path(value).absolute()
            data = _read_regular(source.parent, source.name)
            filename = f"{number:02d}{source.suffix}"
            (staged / filename).write_bytes(data)
            manifest["files"][role] = {"path": filename, "sha256": _digest(data)}
        (staged / "snapshot.json").write_bytes(_json(manifest))
        result = read_snapshot(staged)
        staged.rename(destination)
    return {"manifest_sha256": result["manifest_sha256"], "simulation_ready": False}


def copy_snapshot(source, destination, *, expected_manifest_sha256):
    """Carry the identical, verified snapshot through an archive lifecycle."""
    _require_digest(expected_manifest_sha256)
    verified = read_snapshot(source, expected_manifest_sha256=expected_manifest_sha256)
    destination = Path(destination).absolute()
    if destination.exists() or destination.is_symlink():
        raise ValueError("Refusing to replace an existing parameter snapshot")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".parameter-copy-", dir=destination.parent) as temporary:
        staged = Path(temporary) / "snapshot"
        staged.mkdir()
        (staged / "snapshot.json").write_bytes(verified["manifest_bytes"])
        for role, data in verified["contents"].items():
            filename = verified["manifest"]["files"][role]["path"]
            (staged / filename).write_bytes(data)
        read_snapshot(staged, expected_manifest_sha256=expected_manifest_sha256)
        staged.rename(destination)


def _mapped_topology(native, expected_atoms, raw_atom_names):
    """Restore explicit identities while preserving native particle order/bonds."""
    from openmm import app
    native_atoms = list(native.atoms())
    if not isinstance(expected_atoms,list) or len(native_atoms) != len(expected_atoms):
        raise ValueError("Native topology and declared atom inventory have different sizes")
    if (any(not isinstance(atom,dict) or type(atom.get('native_index')) is not int for atom in expected_atoms) or
            [atom.get("native_index") for atom in expected_atoms] != list(range(len(native_atoms)))):
        raise ValueError("An explicit, complete native atom index map is required")
    topology = app.Topology()
    atoms, seen_residues, seen_chains = [], set(), set()
    chain = residue = None
    previous_residue = None
    previous_chain = None
    for actual, expected in zip(native_atoms, expected_atoms):
        identity = expected["identity"]
        # Native names may differ from restored source names, but every such
        # mapping must be explicit. No graph-isomorphism matching is performed.
        if actual.name != identity["atomname"]:
            raise ValueError("Imported native atom name does not match the declared index map")
        if raw_atom_names[actual.index] != expected.get("native_atom_name", identity["atomname"]):
            raise ValueError("Raw native atom name does not match the declared index map")
        if actual.residue.name != identity["resname"]:
            raise ValueError("Imported native residue name does not match the declared index map")
        key = tuple(identity[x] for x in ("chain", "resid", "insertion", "resname"))
        if key[0] != previous_chain:
            if key[0] in seen_chains:
                raise ValueError("Restored chains must be contiguous in native particle order")
            seen_chains.add(key[0])
            chain = topology.addChain(key[0])
            previous_chain = key[0]
        if key != previous_residue:
            if key in seen_residues:
                raise ValueError("Restored residues must be contiguous in native particle order")
            seen_residues.add(key)
            residue = topology.addResidue(key[3], chain, id=key[1], insertionCode=key[2])
            previous_residue = key
        atoms.append(topology.addAtom(identity["atomname"], actual.element, residue, id=str(actual.index + 1)))
    for bond in native.bonds():
        topology.addBond(atoms[bond[0].index], atoms[bond[1].index], type=bond.type, order=bond.order)
    topology.setPeriodicBoxVectors(native.getPeriodicBoxVectors())
    return topology


@dataclass
class ResearchModel:
    topology: object
    positions: object
    system: object
    mechanics_report: dict
    snapshot_sha256: str
    chemical_state: dict
    simulation_ready: bool = False


def load_native_snapshot(folder, *, expected_manifest_sha256):
    """Load a research specimen for testing; require independent term agreement.

    This deliberately constructs an unconstrained, nonperiodic comparison System.
    Production solvation/MD and GROMACS export require additional validation.
    """
    from openmm import app, unit
    from .bonded_site_validation import validate_bonded_site
    _require_digest(expected_manifest_sha256)
    verified = read_snapshot(folder, expected_manifest_sha256=expected_manifest_sha256)
    contents = verified["contents"]
    terms = json.loads(contents["expected_terms"])
    if (not isinstance(terms,dict) or terms.get("scope") != "full_system"
            or terms.get("native_reference", {}).get("sha256") != _digest(contents["native_topology"])):
        raise ValueError("Full native parameter manifest must be bound to the exact Amber topology")
    state = json.loads(contents["chemical_state"])
    if not isinstance(state, dict) or not state:
        raise ValueError("Explicit chemical-state metadata is required")
    # OpenMM's restart reader requires a filename. Parse only verified bytes in
    # a private temporary directory, rather than rereading an unchecked source.
    with tempfile.TemporaryDirectory(prefix="dynamol-native-model-") as temporary:
        coordinate_path = Path(temporary) / "model.inpcrd"
        coordinate_path.write_bytes(contents["coordinates"])
        parameter_path = Path(temporary) / "model.prmtop"
        parameter_path.write_bytes(contents["native_topology"])
        parameters = app.AmberPrmtopFile(str(parameter_path))
        coordinates = app.AmberInpcrdFile(str(coordinate_path))
        system = parameters.createSystem(nonbondedMethod=app.NoCutoff, constraints=None,
                                         rigidWater=False, removeCMMotion=False)
    topology = _mapped_topology(parameters.topology, terms["atoms"], parameters._prmtop.getAtomNames())
    xyz = np.asarray(coordinates.positions.value_in_unit(unit.nanometer))
    if xyz.shape != (system.getNumParticles(), 3) or not np.isfinite(xyz).all():
        raise ValueError("Native restart coordinates are nonfinite or do not match particle count")
    report = validate_bonded_site(topology, system, terms)
    if not report["accepted"]:
        raise ValueError(f"Native parameter snapshot failed mechanics validation: {report['failures'][:3]}")
    if report.get("additional_forces_not_used_as_parameter_evidence"):
        raise ValueError("Native model contains forces outside the current snapshot validation scope")
    return ResearchModel(topology, coordinates.positions, system, report,
                         verified["manifest_sha256"], state)
