"""Pinned, explicit residue support; never convert a modification to its parent.

This registry exposes only force-field templates whose chemistry and source are
known. It does not parameterize arbitrary unnatural amino acids. Phosphates use
fixed dianionic templates; the requested pH is recorded, not a computed pKa.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET

DATA = Path(__file__).parent / "data" / "modified_residues"


@dataclass(frozen=True)
class ModifiedResidue:
    name: str
    parent: str
    charge: int
    forcefield: str
    label: str


REGISTRY = {
    "SEP": ModifiedResidue("SEP", "SER", -2, "phosaa14SB", "Phosphoserine"),
    "TPO": ModifiedResidue("TPO", "THR", -2, "phosaa14SB", "Phosphothreonine"),
    "PTR": ModifiedResidue("PTR", "TYR", -2, "phosaa14SB", "Phosphotyrosine"),
    "HYP": ModifiedResidue("HYP", "PRO", 0, "ff14SB", "4-hydroxyproline"),
}
SUPPORTED_MODIFIED = frozenset(REGISTRY)
PARENT_RESIDUES = {name: entry.parent for name, entry in REGISTRY.items()}
PHOSPHORYLATED = frozenset({"SEP", "TPO", "PTR"})
STEREOCENTERS = {
    "SEP": (("CA", ("N", "C", "CB")),),
    "TPO": (("CA", ("N", "C", "CB")), ("CB", ("CA", "OG1", "CG2"))),
    "PTR": (("CA", ("N", "C", "CB")),),
    "HYP": (("CA", ("N", "C", "CB")), ("CG", ("CB", "CD", "OD1"))),
}


def _manifest():
    return json.loads((DATA / "manifest.json").read_text())


def _verified(name):
    data = (DATA / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != _manifest()["files"][name]["sha256"]:
        raise ValueError(f"Bundled modified-residue data checksum mismatch: {name}.")
    return data


@lru_cache(maxsize=1)
def _templates():
    from openmm import app
    ff = app.ForceField("amber14/protein.ff14SB.xml", io.StringIO(_verified("phosaa14SB.xml").decode()))
    return {name: ff._templates[name] for name in REGISTRY}


@lru_cache(maxsize=4)
def _heavy_coordinates(name):
    import gemmi
    block = gemmi.cif.read_string(_verified(f"{name}.cif").decode()).sole_block()
    rows = block.find("_chem_comp_atom.", ["atom_id", "type_symbol", "pdbx_model_Cartn_x_ideal", "pdbx_model_Cartn_y_ideal", "pdbx_model_Cartn_z_ideal"])
    return {row[0]: (row[1], tuple(float(value) * .1 for value in list(row)[2:])) for row in rows if row[1] != "H"}


@lru_cache(maxsize=1)
def register_topology_definitions():
    """Register exact internal bonds, peptide links, and fixed-state hydrogens.

    Call before reading a PDB or constructing PDBFixer. The public OpenMM loader
    APIs add these definitions globally; no polymer names or atoms are renamed.
    """
    from openmm import app
    bonds_root, hydrogens_root = ET.Element("Residues"), ET.Element("Residues")
    for name, template in _templates().items():
        bonds = ET.SubElement(bonds_root, "Residue", name=name)
        hydrogens = ET.SubElement(hydrogens_root, "Residue", name=name)
        ET.SubElement(bonds, "Bond", {"from": "-C", "to": "N"})
        ET.SubElement(bonds, "Bond", {"from": "C", "to": "OXT"})
        for first, second in template.bonds:
            a, b = template.atoms[first], template.atoms[second]
            ET.SubElement(bonds, "Bond", {"from": a.name, "to": b.name})
            if a.element == app.element.hydrogen or b.element == app.element.hydrogen:
                hydrogen, parent = (a, b) if a.element == app.element.hydrogen else (b, a)
                ET.SubElement(hydrogens, "H", name=hydrogen.name, parent=parent.name)
    app.Topology.loadBondDefinitions(io.StringIO(ET.tostring(bonds_root, encoding="unicode")))
    app.Modeller.loadHydrogenDefinitions(io.StringIO(ET.tostring(hydrogens_root, encoding="unicode")))


def register_fixer_templates(fixer):
    """Install CCD ideal heavy geometry, including OXT marked terminal only."""
    from openmm import app, unit, Vec3
    register_topology_definitions()
    for name, template in _templates().items():
        coordinates = _heavy_coordinates(name)
        topology = app.Topology()
        residue = topology.addResidue(name, topology.addChain())
        atoms, positions, terminal = {}, [], []
        for atom in template.atoms:
            if atom.element == app.element.hydrogen:
                continue
            symbol, position = coordinates[atom.name]
            atoms[atom.name] = topology.addAtom(atom.name, app.element.get_by_symbol(symbol), residue)
            positions.append(Vec3(*position))
            terminal.append(False)
        symbol, position = coordinates["OXT"]
        atoms["OXT"] = topology.addAtom("OXT", app.element.get_by_symbol(symbol), residue)
        positions.append(Vec3(*position))
        terminal.append(True)
        for first, second in template.bonds:
            a, b = template.atoms[first].name, template.atoms[second].name
            if a in atoms and b in atoms:
                topology.addBond(atoms[a], atoms[b])
        topology.addBond(atoms["C"], atoms["OXT"])
        fixer.registerTemplate(topology, positions * unit.nanometer, terminal)


def _selected_names(preparation):
    if preparation is None:
        return set(REGISTRY)
    values = preparation.get("modified_residues", []) if isinstance(preparation, dict) else preparation
    if isinstance(values, dict):
        values = values.get("residues", [])
    return {entry.get("residue", entry.get("name", "")) if isinstance(entry, dict) else str(entry) for entry in values}


def modified_forcefield_preload_files(preparation=None):
    """Reserved compatibility hook; this parameter set needs no prelude."""
    return []


def modified_forcefield_files(preparation=None):
    """Return verified packaged supplemental file paths for the recorded state."""
    if _selected_names(preparation) & PHOSPHORYLATED:
        _verified("phosaa14SB.xml")
        return [str(DATA / "phosaa14SB.xml")]
    return []


def modified_forcefield_provenance():
    """Portable immutable source identifiers; never encode machine paths."""
    manifest = _manifest()
    return {"registry_version": 1, "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "ptr_improper_correction": "Restore the four native Amber PTR ring improper quartets by atom name after generic template matching; numeric terms unchanged. Fully assembled System.xml records the portable result.", "forcefield": "ff14SB + phosaa14SB", "source_repository": manifest["upstream_repository"], "source_commit": manifest["upstream_commit"], "reference": manifest["reference"], "normalization": manifest["normalization"], "files": {name: metadata["sha256"] for name, metadata in manifest["files"].items()}, "protonation": "SEP/TPO/PTR use fixed -2 phosphate states; HYP is neutral. No residue pKa calculation or automatic phosphate titration."}


def inspect_modified(topology, ph=7.0):
    """Describe supported modifications and reject unsupported terminal links."""
    from openmm import app
    records, warnings, blockers = [], [], []
    bonds = list(topology.bonds())
    for residue in topology.residues():
        if residue.name not in REGISTRY:
            continue
        definition = REGISTRY[residue.name]
        atoms = {atom.name: atom for atom in residue.atoms()}
        expected = {atom.name for atom in _templates()[residue.name].atoms if atom.element != app.element.hydrogen}
        external = []
        for first, second in bonds:
            if (first.residue == residue) != (second.residue == residue):
                own, other = (first, second) if first.residue == residue else (second, first)
                external.append((own.name, other.name))
        required = {"N", "C"} if residue.name in PHOSPHORYLATED else {"N"}
        external_names = {first for first, _ in external}
        missing_links = required - external_names
        # Missing backbone atoms may be repaired; a present unlinked backbone
        # atom indicates a terminal residue rather than the internal template.
        missing_links.intersection_update(atoms)
        if missing_links:
            blockers.append(f"{residue.name} {residue.chain.id}:{residue.id} has an unsupported terminal or broken peptide connection ({', '.join(sorted(missing_links))}). The packaged {definition.forcefield} template supports an internal residue{' or C-terminal HYP' if residue.name == 'HYP' else ''}; supply a supported capped/connected structure or a validated terminal template.")
        unexpected = [pair for pair in external if pair not in {("N", "C"), ("C", "N")}]
        if unexpected:
            blockers.append(f"{residue.name} {residue.chain.id}:{residue.id} has an unsupported covalent crosslink {unexpected}; a matching covalent residue template is required.")
        note = (f"Fixed dianionic phosphate (-2 e) from the published phosaa14SB {residue.name} template. Requested pH {ph:g} does not calculate or change its protonation state; review bound-state chemistry." if residue.name in PHOSPHORYLATED else "Neutral 4-hydroxyproline from ff14SB; original modification is retained.")
        if residue.name in PHOSPHORYLATED:
            warnings.append(f"{residue.name} {residue.chain.id}:{residue.id}: {note}")
        records.append({"chain": residue.chain.id, "resid": residue.id, "insertion_code": (residue.insertionCode or "").strip(), "residue": residue.name, "name": definition.label, "parent_residue": definition.parent, "template": "CHYP" if residue.name == "HYP" and "C" not in external_names else residue.name, "formal_charge": -1 if residue.name == "HYP" and "C" not in external_names else definition.charge, "forcefield": definition.forcefield, "ph": ph, "protonation_note": note, "missing_heavy_atoms": sorted(expected - set(atoms)), "external_bonds": [list(pair) for pair in external], "supported": not missing_links and not unexpected})
    return {"residues": records, "warnings": warnings, "blockers": blockers}


def modified_stereochemistry_report(topology, positions):
    """Check CA and modification-specific centers against the CCD ideal template."""
    import numpy as np
    from openmm import unit
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
    centers, violations = [], []
    for residue in topology.residues():
        if residue.name not in STEREOCENTERS:
            continue
        names = {atom.name: atom.index for atom in residue.atoms()}
        ideal = _heavy_coordinates(residue.name)
        for center, neighbors in STEREOCENTERS[residue.name]:
            if not all(name in names for name in (center, *neighbors)):
                continue
            observed = float(np.linalg.det([xyz[names[n]] - xyz[names[center]] for n in neighbors]))
            expected = float(np.linalg.det([np.array(ideal[n][1]) - ideal[center][1] for n in neighbors]))
            record = {"chain": residue.chain.id, "resid": residue.id, "insertion_code": (residue.insertionCode or "").strip(), "residue": residue.name, "center": center, "signed_volume_nm3": observed, "template_signed_volume_nm3": expected}
            centers.append(record)
            if abs(observed) < 1e-4 or observed * expected <= 0:
                violations.append(record)
    return {"method": "Signed stereocenter volume versus vendored CCD ideal geometry; checks phosphate residue CA, TPO CB, and HYP CG without replacing the modification.", "checked_centers": len(centers), "centers": centers, "violations": violations}


class _NativePTRImpropers:
    """Bounded correction for OpenMM's ambiguous aromatic improper matching."""
    quartets = (("CE2", "CG", "CD2", "HD2"), ("CZ", "CD2", "CE2", "HE2"), ("CD1", "CZ", "CE1", "HE1"), ("CE1", "CG", "CD1", "HD1"))

    def createForce(self, system, data, nonbondedMethod, nonbondedCutoff, args):
        from openmm import PeriodicTorsionForce, unit
        residues = {atom.residue for atom in data.atoms if atom.residue.name == "PTR"}
        if not residues:
            return
        lookup = {}
        for residue in residues:
            atoms = {atom.name: atom.index for atom in residue.atoms()}
            for names in self.quartets:
                if not all(name in atoms for name in names):
                    raise ValueError(f"PTR {residue.chain.id}:{residue.id} is missing atoms required for its native improper terms.")
                indices = tuple(atoms[name] for name in names)
                lookup[(frozenset(indices), indices[2])] = indices
        found = set()
        for force in system.getForces():
            if not isinstance(force, PeriodicTorsionForce):
                continue
            for index in range(force.getNumTorsions()):
                first, second, center, fourth, periodicity, phase, k = force.getTorsionParameters(index)
                key = (frozenset((first, second, center, fourth)), center)
                if key in lookup:
                    if periodicity != 2 or abs(k.value_in_unit(unit.kilojoule_per_mole) - 4.6024) > 1e-8:
                        raise ValueError("Unexpected PTR improper parameter; the pinned force-field correction cannot be applied.")
                    force.setTorsionParameters(index, *lookup[key], periodicity, phase, k)
                    found.add(key)
        if found != set(lookup):
            raise ValueError("The assembled PTR system is missing native aromatic improper terms.")


def register_modified_forcefield(forcefield):
    """Apply native PTR atom ordering to every System created by this factory.

    No numeric parameters change. This uses ForceField's generator interface
    and public System/PeriodicTorsionForce methods. System XML serialization
    captures the corrected terms without requiring Python on the consumer.
    """
    if not any(isinstance(generator, _NativePTRImpropers) for generator in forcefield.getGenerators()):
        forcefield.registerGenerator(_NativePTRImpropers())
    return forcefield
