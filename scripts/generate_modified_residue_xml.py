"""Rebuild namespaced phosaa14SB XML from pinned upstream sources.

Only type IDs, supported residue filtering, and copied phosphorus type/LJ
entries change. Numeric parameters and residue connectivity are preserved.
PTR improper atom ordering is corrected by register_modified_forcefield();
prepared System.xml is the portable, fully assembled force-field artifact.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from openmm import app

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "backend/data/modified_residues"


def generate():
    tree = ET.parse(DATA / "phosaa14SB.source.xml")
    ff = tree.getroot()
    base = ET.parse(Path(app.__file__).parent / "data/amber14/protein.ff14SB.xml").getroot()
    full = ET.parse(DATA / "ff14SB.source.xml").getroot()
    base_types = {node.attrib["name"].removeprefix("protein-"): node.attrib["name"] for node in base.find("AtomTypes")}
    ff.find("AtomTypes").append(deepcopy(next(t for t in full.find("AtomTypes") if t.attrib["name"] == "P")))
    ff.find("NonbondedForce").append(deepcopy(next(t for t in full.find("NonbondedForce") if t.attrib.get("class") == "P")))
    custom = {t.attrib["name"]: f"protein-{t.attrib['name']}" for t in ff.find("AtomTypes")}
    for node in ff.find("AtomTypes"):
        node.attrib["name"] = custom[node.attrib["name"]]
    for node in ff.iter():
        for key, value in list(node.attrib.items()):
            if key == "type" or key.startswith("type") and key[4:].isdigit():
                node.attrib[key] = custom[value] if value in custom else base_types[value]
    for residue in list(ff.find("Residues")):
        if residue.attrib["name"] not in {"SEP", "TPO", "PTR"}:
            ff.find("Residues").remove(residue)
    ET.indent(tree, space="  ")
    return ET.tostring(ff, encoding="unicode")


if __name__ == "__main__":
    path = DATA / "phosaa14SB.xml"
    path.write_text(generate())
    manifest = json.loads((DATA / "manifest.json").read_text())
    manifest["normalization"] = "Only atom type identifiers are namespaced for OpenMM amber14/protein.ff14SB.xml. Exact phosphorus type/LJ entries are copied from the same upstream amber/ff14SB.xml. Only SEP/TPO/PTR residue templates are exposed. Charges, numeric parameters, classes, and retained residue graphs are unchanged. PTR improper atom ordering is corrected separately by register_modified_forcefield(); use the assembled System.xml for portable exact terms."
    manifest["ptr_improper_quartets"] = [["CE2", "CG", "CD2", "HD2"], ["CZ", "CD2", "CE2", "HE2"], ["CD1", "CZ", "CE1", "HE1"], ["CE1", "CG", "CD1", "HD1"]]
    manifest["load_order"] = ["amber14/protein.ff14SB.xml", "phosaa14SB.xml"]
    manifest["files"].pop("phosaa14SB.types.xml", None)
    manifest["files"][path.name]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    (DATA / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
