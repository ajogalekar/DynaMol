"""Gross geometry screening of modeled internal loops; no coordinate mutation.

These deliberately broad bounds detect invalid starting geometry, not native
conformations. Chirality and subsequent force-field/runtime checks are separate.
Coordinates are nonperiodic, in nm; prep must supply whole, unwrapped molecules.
"""
from __future__ import annotations

import numpy as np
from openmm import unit
from scipy.spatial import cKDTree

from .residue_identity import protein_residue_keys, residue_key


_RADII = {"C": .170, "N": .155, "O": .152, "S": .180, "P": .180}
_BACKBONE_BOUNDS = {frozenset(pair): bounds for pair, bounds in (
    (("N", "CA"), (1.25, 1.70)), (("CA", "C"), (1.30, 1.75)),
    (("C", "O"), (1.05, 1.45)), (("C", "OXT"), (1.05, 1.45)),
)}
_COVALENT_RADII_A = {"C": .76, "N": .71, "O": .66, "S": 1.05, "P": 1.07}


def _angle(points):
    first, last = points[0] - points[1], points[2] - points[1]
    denominator = np.linalg.norm(first) * np.linalg.norm(last)
    if denominator < 1e-12:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(first, last) / denominator, -1, 1))))


def _dihedral(points):
    first, middle, last = points[0] - points[1], points[2] - points[1], points[3] - points[2]
    length = np.linalg.norm(middle)
    if length < 1e-12:
        return None
    middle = middle / length
    first = first - np.dot(first, middle) * middle
    last = last - np.dot(last, middle) * middle
    if min(np.linalg.norm(first), np.linalg.norm(last)) < 1e-12:
        return None
    return float(np.degrees(np.arctan2(np.dot(np.cross(middle, first), last), np.dot(first, last))))


def loop_geometry_report(topology, positions, modeled_residue_keys, environment_positions=()):
    """Screen all heavy atoms of exact modeled residue identities.

    Keys are (chain, residue ID, insertion code, residue name). Additional
    retained molecules are (xyz_nm, elemental_vdw_radius_nm) pairs, and must
    not duplicate atoms already in topology. No external covalent connection
    can be inferred from that coordinate-only environment representation.
    """
    report = {
        "method": "Gross nonperiodic heavy-atom geometry screen for all modeled loop residues; not a native-loop, rotamer, Ramachandran or energy validation.",
        "accepted": False, "errors": [], "warnings": [], "soft_overlaps": [],
        "gross_collisions": [], "requires_minimization": False,
        "bond_checks": [], "angle_checks": [], "omega_checks": [], "phi_psi": [],
        "criteria": {"backbone_bond_bounds_angstrom": {"N-CA": [1.25, 1.70], "CA-C": [1.30, 1.75], "C-O": [1.05, 1.45], "peptide_C-N": [1.15, 1.60]},
                     "other_heavy_bonds": "0.70–1.25 times summed elemental covalent radii (C/N/O/S/P); generous gross bounds, not equilibrium lengths",
                     "backbone_angle_bounds_degrees": [80, 150], "omega_cis_or_trans_tolerance_degrees": 35,
                     "gross_contact_min_angstrom": 1.3, "soft_vdw_ratio": .78,
                     "contact_exclusions": "self, directly bonded and 1–3 atom pairs; hydrogen contacts excluded"},
    }
    errors = report["errors"]
    atoms, residues = list(topology.atoms()), list(topology.residues())
    try:
        keys = set()
        for key in modeled_residue_keys:
            if len(key) != 4:
                raise ValueError("Modeled residue keys require chain, residue ID, insertion code and name.")
            keys.add((str(key[0]), str(key[1]), str(key[2] or "").strip(), str(key[3])))
        xyz = np.asarray(positions.value_in_unit(unit.nanometer) if hasattr(positions, "value_in_unit") else positions, dtype=float)
        external = list(environment_positions)
        external_xyz = np.asarray([row[0] for row in external], dtype=float) if external else np.empty((0, 3))
        external_radii = np.asarray([row[1] for row in external], dtype=float)
    except (TypeError, ValueError, IndexError) as exc:
        errors.append(f"Invalid loop geometry input: {exc}")
        return report
    report.update(modeled_residue_keys=[list(key) for key in sorted(keys)], external_environment_atoms=len(external))
    if xyz.shape != (len(atoms), 3) or not np.isfinite(xyz).all():
        errors.append("Topology coordinates must be a finite N×3 array matching the atom count.")
        return report
    if external_xyz.shape != (len(external), 3) or external_radii.shape != (len(external),):
        errors.append("Each external environment atom requires exactly three nm coordinates and one radius.")
        return report
    if not np.isfinite(external_xyz).all() or not np.isfinite(external_radii).all() or np.any(external_radii <= 0):
        errors.append("External environment coordinates and positive radii must be finite.")
        return report
    selected = [residue for residue in residues if residue_key(residue) in keys]
    found = [residue_key(residue) for residue in selected]
    if set(found) != keys or len(found) != len(set(found)):
        errors.append("Every modeled residue key must match exactly one topology residue.")
        return report
    modeled = {atom.index for residue in selected for atom in residue.atoms()}
    heavy = {atom.index for atom in atoms if atom.element is not None and atom.element.symbol not in {"H", "D"}}
    if any(atoms[index].element is None for index in modeled):
        errors.append("A modeled atom has no chemical element.")
        return report
    modeled_heavy = modeled & heavy
    report.update(modeled_atom_count=len(modeled), modeled_heavy_atom_count=len(modeled_heavy))
    if not keys:
        report["accepted"] = True
        return report
    if not modeled_heavy:
        errors.append("Modeled loop residues contain no heavy atoms or required backbone.")
        return report

    named = {}
    for residue in residues:
        members = list(residue.atoms())
        named[residue] = {atom.name: atom for atom in members}
        if residue in selected and len(named[residue]) != len(members):
            errors.append(f"Modeled residue {residue_key(residue)} has duplicate atom names.")
    adjacency = [set() for _ in atoms]
    bonds = set()
    for first, second in topology.bonds():
        adjacency[first.index].add(second.index)
        adjacency[second.index].add(first.index)
        bonds.add(tuple(sorted((first.index, second.index))))

    def required(residue, names):
        missing = [name for name in names if name not in named[residue]]
        if missing:
            errors.append(f"Loop or anchor residue {residue_key(residue)} lacks required backbone atoms: {', '.join(missing)}.")
            return False
        wanted = {"N": "N", "CA": "C", "C": "C", "O": "O"}
        if any(named[residue][name].element is None or named[residue][name].element.symbol != wanted[name] for name in names):
            errors.append(f"Loop or anchor residue {residue_key(residue)} has incorrect backbone elements.")
            return False
        return True

    def require_bond(first, second, label):
        if tuple(sorted((first.index, second.index))) not in bonds:
            errors.append(f"Missing expected {label} bond: {residue_key(first.residue)}/{first.name}–{residue_key(second.residue)}/{second.name}.")

    def check_angle(members, label):
        value = _angle(xyz[[atom.index for atom in members]])
        accepted = value is not None and 80 <= value <= 150
        report["angle_checks"].append({"atoms": [a.index for a in members], "kind": label, "degrees": value, "accepted": accepted})
        if not accepted:
            errors.append(f"Modeled backbone angle {label} at atom {members[1].index} is undefined or outside 80–150°.")

    for residue in selected:
        if required(residue, ("N", "CA", "C", "O")):
            names = named[residue]
            for first, second in (("N", "CA"), ("CA", "C"), ("C", "O")):
                require_bond(names[first], names[second], "backbone")
            check_angle([names[name] for name in ("N", "CA", "C")], "N-CA-C")
            check_angle([names[name] for name in ("CA", "C", "O")], "CA-C-O")
            residue_heavy = {a.index for a in residue.atoms()} & heavy
            reached, pending = set(), [names["CA"].index]
            while pending:
                index = pending.pop()
                if index not in reached:
                    reached.add(index)
                    pending.extend((adjacency[index] & residue_heavy) - reached)
            if residue_heavy - reached:
                errors.append(f"Modeled residue {residue_key(residue)} has heavy atoms disconnected from its backbone.")

    # Use chain objects/order, never author-number arithmetic. Modified protein
    # neighbors remain anchors even when not eligible for loop construction.
    protein_keys = protein_residue_keys(None, topology, positions)
    peptide_pairs = set()
    for chain in topology.chains():
        protein = [r for r in chain.residues() if residue_key(r) in protein_keys or r in selected]
        for index, residue in enumerate(protein):
            if residue not in selected:
                continue
            if index == 0 or index == len(protein) - 1:
                errors.append(f"Modeled internal residue {residue_key(residue)} lacks a flanking protein residue.")
            if index:
                peptide_pairs.add((protein[index - 1], residue))
            if index + 1 < len(protein):
                peptide_pairs.add((residue, protein[index + 1]))
            if 0 < index < len(protein) - 1:
                previous, following = protein[index - 1], protein[index + 1]
                phi = [named[previous].get("C"), *[named[residue].get(n) for n in ("N", "CA", "C")]]
                psi = [*[named[residue].get(n) for n in ("N", "CA", "C")], named[following].get("N")]
                report["phi_psi"].append({"residue": list(residue_key(residue)),
                    "phi_degrees": _dihedral(xyz[[a.index for a in phi]]) if all(phi) else None,
                    "psi_degrees": _dihedral(xyz[[a.index for a in psi]]) if all(psi) else None})
    for first, second in sorted(peptide_pairs, key=lambda pair: (pair[0].index, pair[1].index)):
        if not (required(first, ("CA", "C", "O")) and required(second, ("N", "CA"))):
            continue
        a, b = named[first], named[second]
        require_bond(a["C"], b["N"], "peptide")
        check_angle([a["CA"], a["C"], b["N"]], "CA-C-N")
        check_angle([a["O"], a["C"], b["N"]], "O-C-N")
        check_angle([a["C"], b["N"], b["CA"]], "C-N-CA")
        members = [a["CA"], a["C"], b["N"], b["CA"]]
        omega = _dihedral(xyz[[atom.index for atom in members]])
        deviation = min(abs(omega), 180 - abs(omega)) if omega is not None else None
        accepted = deviation is not None and deviation <= 35
        state = "undefined" if omega is None else "cis" if abs(omega) <= 90 else "trans"
        report["omega_checks"].append({"atoms": [atom.index for atom in members], "residues": [list(residue_key(r)) for r in (first, second)],
            "degrees": omega, "state": state, "deviation_degrees": deviation, "accepted": accepted})
        if not accepted:
            errors.append(f"Peptide omega before {residue_key(second)} is undefined or more than 35° from cis/trans.")
        elif state == "cis" and second.name != "PRO":
            report["warnings"].append(f"Modeled non-proline cis peptide before {residue_key(second)} is uncertain; planarity alone does not validate its conformation.")

    for i, j in sorted(bonds):
        if not ({i, j} & modeled_heavy) or i not in heavy or j not in heavy:
            continue
        first, second = atoms[i], atoms[j]
        symbols = [a.element.symbol for a in (first, second)]
        peptide = (first.residue, second.residue) in peptide_pairs or (second.residue, first.residue) in peptide_pairs
        if peptide and {first.name, second.name} == {"C", "N"}:
            bounds, kind = (1.15, 1.60), "peptide_C-N"
        elif first.residue == second.residue and frozenset((first.name, second.name)) in _BACKBONE_BOUNDS:
            bounds, kind = _BACKBONE_BOUNDS[frozenset((first.name, second.name))], "backbone"
        elif all(symbol in _COVALENT_RADII_A for symbol in symbols):
            radius_sum = sum(_COVALENT_RADII_A[symbol] for symbol in symbols)
            bounds, kind = (.70 * radius_sum, 1.25 * radius_sum), "other_heavy"
        else:
            errors.append(f"No gross covalent-length bounds for modeled bond {i}–{j} ({'-'.join(symbols)}).")
            continue
        distance = float(np.linalg.norm(xyz[i] - xyz[j]) * 10)
        accepted = bounds[0] <= distance <= bounds[1]
        report["bond_checks"].append({"atoms": [i, j], "kind": kind, "angstrom": distance,
                                     "bounds_angstrom": list(bounds), "accepted": accepted})
        if not accepted:
            errors.append(f"Modeled {kind} bond {i}–{j} is outside gross {bounds[0]:.2f}–{bounds[1]:.2f} Å bounds.")

    indices = sorted(heavy)
    combined = np.concatenate((xyz[indices], external_xyz))
    radii = np.concatenate(([_RADII.get(atoms[i].element.symbol, .170) for i in indices], external_radii))
    tree = cKDTree(combined)
    max_radius = float(radii.max())
    seen = set()
    for index in sorted(modeled_heavy):
        excluded = {index} | adjacency[index]
        for neighbor in adjacency[index]:
            excluded |= adjacency[neighbor]
        radius = _RADII.get(atoms[index].element.symbol, .170)
        for nearby in tree.query_ball_point(xyz[index], max(.13, .78 * (radius + max_radius))):
            other = indices[nearby] if nearby < len(indices) else len(atoms) + nearby - len(indices)
            pair = tuple(sorted((index, other)))
            if other in excluded or pair in seen:
                continue
            seen.add(pair)
            distance = float(np.linalg.norm(xyz[index] - combined[nearby]))
            ratio = distance / float(radius + radii[nearby])
            record = {"atoms": list(pair), "distance_nm": distance, "vdw_distance_ratio": ratio,
                      "external_environment_index": other - len(atoms) if other >= len(atoms) else None}
            if ratio < .78:
                report["soft_overlaps"].append(record)
            if distance < .13:
                report["gross_collisions"].append(record)
    if report["gross_collisions"]:
        errors.append("Modeled loop heavy atoms retain a gross nonbonded collision below 1.3 Å.")
    report["soft_overlap_score"] = float(sum((1 - row["vdw_distance_ratio"] / .78) ** 2 for row in report["soft_overlaps"]))
    report["requires_minimization"] = bool(report["soft_overlaps"])
    report["accepted"] = not errors
    return report
