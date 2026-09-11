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
    from openmm import app

    if solvent not in {"implicit", "explicit"}:
        raise ValueError("Choose implicit or explicit solvent.")
    bundle = _bundle(preparation)
    if bundle is not None and solvent != "explicit":
        raise ValueError("Prepared protein–ligand complexes require explicit TIP3P water. GBn2 implicit parameters are not available for these ligands.")
    base = ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml" if solvent == "implicit" else "amber14/tip3p.xml"]
    verified = _validated_xml(Path(folder), preparation)
    # Feed the verified bytes to OpenMM, not an unchecked second path read.
    forcefield = app.ForceField(*base, *(io.StringIO(data.decode("utf-8")) for _, data in verified))
    return forcefield, base + [str(path.relative_to(Path(folder).resolve())) for path, _ in verified]


def copy_ligand_parameters(source_folder: Path, destination_folder: Path, preparation: dict | None) -> None:
    """Copy the complete ligands/ archive and verify the destination manifest.

    Native charge outputs, atom maps and logs travel alongside the XML files.
    No symlinks or special files are copied, and existing bundles are not
    overwritten. Dataset/job folders should be new when this helper is used.
    """
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
