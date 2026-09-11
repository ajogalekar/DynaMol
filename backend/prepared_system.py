"""Validate and carry the exact ligand parameter bundle between local stages.

The public helpers accept a dataset/job folder and its ``preparation`` record.
They deliberately do not discover arbitrary XML files or regenerate charges.
"""
from __future__ import annotations

import hashlib
import io
import re
import shutil
import tempfile
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree


MAX_XML_BYTES = 16 * 1024 * 1024
MAX_BUNDLE_BYTES = 128 * 1024 * 1024
MAX_BUNDLE_FILES = 1024


def _bundle(preparation: dict | None) -> dict | None:
    if preparation is not None and not isinstance(preparation, dict):
        raise ValueError("The saved preparation state is malformed. Prepare the structure again.")
    if not preparation or "ligand_parameters" not in preparation:
        return None
    bundle = preparation["ligand_parameters"]
    if not isinstance(bundle, dict) or bundle.get("forcefield") != "GAFF2" or bundle.get("charge_method") != "AM1-BCC":
        raise ValueError("The prepared ligand parameter bundle is invalid or uses an unsupported force field. Prepare the complex again.")
    if bundle.get("requires_explicit_solvent") is not True:
        raise ValueError("Prepared GAFF2 ligands require explicit TIP3P solvent; implicit ligand parameters have not been validated.")
    if not isinstance(bundle.get("ligands"), list) or not bundle["ligands"] or not all(isinstance(ligand, dict) for ligand in bundle["ligands"]):
        raise ValueError("The prepared ligand bundle is missing its component provenance. Prepare the complex again.")
    if not isinstance(bundle.get("files"), list) or not 1 <= len(bundle["files"]) <= 128:
        raise ValueError("The prepared ligand bundle has no valid parameter-file manifest. Prepare the complex again.")
    return bundle


def _safe_path(folder: Path, value: object) -> Path:
    if not isinstance(value, str) or "\\" in value or ":" in value:
        raise ValueError("Ligand parameter paths must be relative files inside ligands/.")
    relative = PurePosixPath(value)
    if (relative.is_absolute() or len(relative.parts) < 2 or relative.parts[0] != "ligands"
            or ".." in relative.parts or str(relative) != value or relative.suffix.lower() != ".xml"):
        raise ValueError("Ligand parameter paths must be relative XML files inside ligands/.")
    folder = Path(folder).resolve()
    path = folder.joinpath(*relative.parts)
    if not path.resolve().is_relative_to(folder / "ligands"):
        raise ValueError("A ligand parameter path escapes its bundle.")
    for parent in (path, *path.parents):
        if parent == folder:
            break
        if parent.is_symlink():
            raise ValueError("Symlinks are not allowed in a prepared ligand parameter bundle.")
    return path


def _validated_xml(folder: Path, preparation: dict | None, *, required: bool = False) -> list[tuple[Path, bytes]]:
    bundle = _bundle(preparation)
    if bundle is None:
        if required:
            raise ValueError("This complex has no verified ligand parameters. Run protein and ligand preparation before adding water or starting MD; no molecules were removed.")
        return []
    result = []
    seen = set()
    total_bytes = 0
    for entry in bundle["files"]:
        if not isinstance(entry, dict):
            raise ValueError("Malformed ligand parameter-file manifest.")
        path = _safe_path(folder, entry.get("path"))
        if path in seen:
            raise ValueError("A ligand parameter file is listed more than once.")
        seen.add(path)
        expected = entry.get("sha256")
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
            raise ValueError("A ligand parameter file is missing its SHA-256 checksum.")
        if not path.is_file():
            raise ValueError(f"Prepared ligand parameter file is missing: {entry['path']}. Prepare the complex again.")
        if path.stat().st_size > MAX_XML_BYTES:
            raise ValueError("A ligand parameter XML file exceeds the local size limit.")
        data = path.read_bytes()
        total_bytes += len(data)
        if total_bytes > MAX_BUNDLE_BYTES:
            raise ValueError("The ligand parameter XML files exceed the local total size limit.")
        if hashlib.sha256(data).hexdigest() != expected.lower():
            raise ValueError(f"Prepared ligand parameter checksum mismatch: {entry['path']}. The saved molecular state cannot be verified; prepare the complex again.")
        if b"<!DOCTYPE" in data.upper() or b"<!ENTITY" in data.upper():
            raise ValueError("Ligand parameter XML must not contain document types or entities.")
        try:
            document = ElementTree.fromstring(data)
        except ElementTree.ParseError as exc:
            raise ValueError("The saved ligand parameter XML cannot be parsed.") from exc
        if document.tag != "ForceField" or any(node.tag.rsplit("}", 1)[-1] in {"Include", "Script", "InitializationScript"} for node in document.iter()):
            raise ValueError("Ligand XML must be a self-contained ForceField without includes or executable scripts.")
        result.append((path, data))
    return result


def ligand_parameter_files(folder: Path, preparation: dict | None, *, required: bool = False) -> list[Path]:
    """Return verified XML paths; missing bundles return [] unless required."""
    return [path for path, _ in _validated_xml(Path(folder), preparation, required=required)]


def load_prepared_forcefield(folder: Path, preparation: dict | None, *, solvent: str = "explicit"):
    """Return ``(OpenMM ForceField, provenance_files)`` without changing state.

    Protein-only implicit GBn2 remains available. Ligand bundles are accepted
    only with explicit TIP3P; each declared XML is verified before parsing.
    """
    from .forcefield_identity import IdentityForceField
    from .modified_residues import modified_forcefield_files, modified_forcefield_provenance, register_modified_forcefield

    if solvent not in {"implicit", "explicit"}:
        raise ValueError("Choose implicit or explicit solvent.")
    if (preparation or {}).get("ions") and solvent != "explicit":
        raise ValueError("Prepared ions require explicit TIP3P water; the retained ion model has no validated implicit-solvent treatment.")
    bundle = _bundle(preparation)
    modified = (preparation or {}).get("modified_residues", [])
    if modified and solvent != "explicit":
        raise ValueError("Prepared modified amino acids require explicit TIP3P water in this workflow.")
    if modified and (preparation or {}).get("modified_residue_parameters") != modified_forcefield_provenance():
        raise ValueError("The saved modified-residue parameters do not match the installed version. Prepare the structure again to record its residue parameters.")
    if modified:
        _validate_modified_snapshot(Path(folder), preparation)
    if bundle is not None and solvent != "explicit":
        raise ValueError("Prepared protein–ligand complexes require explicit TIP3P water. GBn2 implicit parameters are not available for these ligands.")
    base = ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml" if solvent == "implicit" else "amber14/tip3p.xml"]
    verified = _validated_xml(Path(folder), preparation)
    # Feed the verified bytes to OpenMM, not an unchecked second path read.
    supplemental = modified_forcefield_files(preparation or {})
    if modified:
        supplemental = [str(Path(folder) / "residue-parameters" / Path(path).name) for path in supplemental]
    identities = _ligand_atom_maps(bundle, verified)
    forcefield = IdentityForceField(*base, *supplemental, *(io.StringIO(data.decode("utf-8")) for _, data in verified),
                                   ligand_atom_maps=identities)
    if modified:
        register_modified_forcefield(forcefield)
    return forcefield, base + ["residue-parameters/" + Path(path).name for path in supplemental] + [str(path.relative_to(Path(folder).resolve())) for path, _ in verified]


def _ligand_atom_maps(bundle, verified):
    """Recover only declared atom identities from the verified parameter bundle."""
    if bundle is None:
        return []
    templates = {}
    for _, data in verified:
        for residue in ElementTree.fromstring(data).findall('Residues/Residue'):
            name = residue.get('name')
            if name in templates:
                raise ValueError('A ligand template is duplicated in the parameter bundle.')
            templates[name] = residue
    identities = []
    for ligand in bundle['ligands']:
        template_name = ligand.get('template_name')
        # Non-DynaMol synthetic fixtures and legacy external parameter bundles
        # keep normal matching. Generated DML templates always require a map.
        if template_name is None:
            continue
        if not isinstance(template_name, str) or not isinstance(ligand.get('key'), str):
            raise ValueError('The prepared ligand template or residue identity is malformed.')
        template = templates.get(template_name)
        key = ligand['key'].split(':')
        if template is None or len(key) != 4:
            raise ValueError('The prepared ligand is missing its template or residue identity. Prepare the complex again.')
        atoms = ligand.get('parameter_atom_map')
        if atoms is None:
            # Earlier DynaMol bundles recorded every heavy atom and used H1,
            # H2,... for added hydrogens. Recover that exact documented naming
            # convention, never a new graph-isomorphism assignment.
            heavy = ligand.get('atom_map', [])
            xml_atoms = template.findall('Atom')
            types = {}
            for _, data in verified:
                for atom_type in ElementTree.fromstring(data).findall('AtomTypes/Type'):
                    types[atom_type.get('name')] = atom_type.get('element')
            if (not isinstance(heavy, list) or not heavy or not all(isinstance(atom, dict) for atom in heavy)
                    or type(ligand.get('hydrogens_added')) is not int
                    or [atom.get('ligand_index') for atom in heavy] != list(range(len(heavy)))
                    or len(heavy) + ligand.get('hydrogens_added', -1) != len(xml_atoms)
                    or any(types.get(atom.get('type')) == 'H' for atom in xml_atoms[:len(heavy)])
                    or any(types.get(atom.get('type')) != 'H' for atom in xml_atoms[len(heavy):])):
                raise ValueError('The legacy ligand atom identity cannot be recovered. Prepare the complex again.')
            names = [atom['prepared_name'] for atom in heavy] + [f'H{i + 1}' for i in range(len(xml_atoms) - len(heavy))]
            atoms = [{'native_index': i, 'prepared_name': name, 'template_atom_name': atom.get('name')}
                     for i, (name, atom) in enumerate(zip(names, xml_atoms))]
        if not isinstance(atoms, list) or not all(isinstance(atom, dict) for atom in atoms):
            raise ValueError('The prepared ligand atom identity map is malformed. Prepare the complex again.')
        identities.append({'residue_key': key, 'template_name': template_name, 'atoms': atoms})
    return identities


def _modified_manifest(preparation):
    if not isinstance(preparation, dict) or not preparation.get("modified_residues"):
        return None
    record = preparation.get("modified_residue_parameters")
    if not isinstance(record, dict) or not isinstance(record.get("files"), dict):
        raise ValueError("Modified residues lack a parameter-file manifest. Prepare the protein again.")
    for name, digest in record["files"].items():
        if not isinstance(name, str) or Path(name).name != name or name in {".", ".."} or not re.fullmatch(r"[0-9a-f]{64}", str(digest)):
            raise ValueError("The modified-residue parameter manifest contains an invalid file or checksum.")
    return record


def _validate_modified_snapshot(folder, preparation):
    import json
    record = _modified_manifest(preparation)
    if record is None:
        return
    root = Path(folder) / "residue-parameters"
    if root.is_symlink() or not root.is_dir():
        raise ValueError("The modified-residue parameter snapshot is missing or is a symlink.")
    total = 0
    for path in root.iterdir():
        if path.is_symlink() or not path.is_file():
            raise ValueError("Modified-residue parameter snapshots cannot contain symlinks or special files.")
        total += path.stat().st_size
    if total > MAX_BUNDLE_BYTES:
        raise ValueError("The modified-residue parameter snapshot exceeds the size limit.")
    for name, digest in record["files"].items():
        path = root / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Modified-residue parameter checksum mismatch: {name}.")
    try:
        manifest = json.loads((root / "manifest.json").read_text())
        if {name: entry["sha256"] for name, entry in manifest["files"].items()} != record["files"]:
            raise ValueError("The modified-residue source manifest does not match its preparation record.")
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise ValueError("The modified-residue source manifest is missing or invalid.") from exc


def snapshot_modified_parameters(folder: Path, preparation: dict) -> None:
    """Freeze the exact curated residue data in a newly prepared job archive."""
    from .modified_residues import DATA, modified_forcefield_provenance
    record = _modified_manifest(preparation)
    if record is None:
        return
    if record != modified_forcefield_provenance():
        raise ValueError("Modified-residue source versions changed before preparation could be saved.")
    folder = Path(folder)
    destination = folder / "residue-parameters"
    if destination.exists() or destination.is_symlink():
        raise ValueError("The preparation already contains modified-residue parameters; refusing to replace them.")
    folder.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".residue-copy-", dir=folder) as temporary:
        staged = Path(temporary) / "residue-parameters"
        staged.mkdir()
        for name in ["manifest.json", *record["files"]]:
            source = DATA / name
            if source.is_symlink() or not source.is_file():
                raise ValueError("Bundled modified-residue data must be regular files.")
            shutil.copy2(source, staged / name)
        _validate_modified_snapshot(Path(temporary), preparation)
        staged.rename(destination)


def copy_modified_parameters(source_folder: Path, destination_folder: Path, preparation: dict | None) -> None:
    if _modified_manifest(preparation) is None:
        return
    source_folder, destination_folder = Path(source_folder).resolve(), Path(destination_folder).resolve()
    _validate_modified_snapshot(source_folder, preparation)
    if source_folder == destination_folder:
        return
    destination_folder.mkdir(parents=True, exist_ok=True)
    target = destination_folder / "residue-parameters"
    if target.exists() or target.is_symlink():
        raise ValueError("The destination already contains modified-residue parameters; refusing to replace them.")
    with tempfile.TemporaryDirectory(prefix=".residue-copy-", dir=destination_folder) as temporary:
        staged = Path(temporary)
        shutil.copytree(source_folder / "residue-parameters", staged / "residue-parameters")
        _validate_modified_snapshot(staged, preparation)
        (staged / "residue-parameters").rename(target)


def copy_ligand_parameters(source_folder: Path, destination_folder: Path, preparation: dict | None) -> None:
    """Copy the complete ligands/ archive and verify the destination manifest.

    Native charge outputs, atom maps and logs travel alongside the XML files.
    No symlinks or special files are copied, and existing bundles are not
    overwritten. Dataset/job folders should be new when this helper is used.
    """
    copy_modified_parameters(source_folder, destination_folder, preparation)
    if not ligand_parameter_files(source_folder, preparation):
        return
    source_folder = Path(source_folder).resolve()
    destination_folder = Path(destination_folder).resolve()
    if source_folder == destination_folder:
        return
    source = source_folder / "ligands"
    entries = list(source.rglob("*"))
    total_bytes = 0
    files = []
    for entry in entries:
        if entry.is_symlink() or (not entry.is_file() and not entry.is_dir()):
            raise ValueError("A ligand archive contains a symlink or unsupported file type.")
        if entry.is_file():
            total_bytes += entry.stat().st_size
            files.append(entry)
    if len(files) > MAX_BUNDLE_FILES or total_bytes > MAX_BUNDLE_BYTES:
        raise ValueError("The ligand parameter archive exceeds local copy limits.")
    destination_folder.mkdir(parents=True, exist_ok=True)
    destination = destination_folder / "ligands"
    if destination.exists() or destination.is_symlink():
        raise ValueError("The destination already contains a ligand bundle; refusing to replace its prepared state.")
    with tempfile.TemporaryDirectory(prefix=".ligand-copy-", dir=destination_folder) as temporary:
        staging = Path(temporary)
        shutil.copytree(source, staging / "ligands")
        ligand_parameter_files(staging, preparation, required=True)
        (staging / "ligands").rename(destination)
