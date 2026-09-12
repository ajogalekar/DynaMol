"""Real PDBFixer/OpenMM preparation and solvation background worker.

Sidechain adjustment is bounded coordinate descent over discrete chi angles,
scored by a documented steric-overlap objective, followed by restrained local
force-field minimization. It is not exhaustive rotamer packing or loop prediction.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import random
import shutil
import signal
import sys
import time
import traceback
from collections import Counter
from pathlib import Path

import numpy as np
from openmm import app, unit

from . import config, storage
from .preparation import (STANDARD_PROTEINS, backbone_gaps, find_missing_residues_preserving_identity,
                          unsupported_missing_residue_message)
from .modified_residues import (SUPPORTED_MODIFIED, PARENT_RESIDUES, register_topology_definitions,
                                register_fixer_templates, inspect_modified, modified_forcefield_provenance,
                                modified_stereochemistry_report)
from .worker import Cancelled, Worker

# Standard chi definitions. Proline is deliberately excluded: ring closure is not sampled.
CHI = {
 "ARG": [("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","NE"),("CG","CD","NE","CZ")],
 "ASN": [("N","CA","CB","CG"),("CA","CB","CG","OD1")],
 "ASP": [("N","CA","CB","CG"),("CA","CB","CG","OD1")],
 "CYS": [("N","CA","CB","SG")],
 "GLN": [("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","OE1")],
 "GLU": [("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","OE1")],
 "HIS": [("N","CA","CB","CG"),("CA","CB","CG","ND1")],
 "ILE": [("N","CA","CB","CG1"),("CA","CB","CG1","CD1")],
 "LEU": [("N","CA","CB","CG"),("CA","CB","CG","CD1")],
 "LYS": [("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","CE"),("CG","CD","CE","NZ")],
 "MET": [("N","CA","CB","CG"),("CA","CB","CG","SD"),("CB","CG","SD","CE")],
 "PHE": [("N","CA","CB","CG"),("CA","CB","CG","CD1")],
 "SER": [("N","CA","CB","OG")],
 "THR": [("N","CA","CB","OG1")],
 "TRP": [("N","CA","CB","CG"),("CA","CB","CG","CD1")],
 "TYR": [("N","CA","CB","CG"),("CA","CB","CG","CD1")],
 "VAL": [("N","CA","CB","CG1")],
}


def atom_key(atom):
    residue = atom.residue
    return (residue.chain.id, residue.id, (residue.insertionCode or "").strip(), residue.name, atom.name)


def torsion_degrees(points):
    p0, p1, p2, p3 = np.asarray(points, dtype=float)
    axis = p2 - p1
    axis /= np.linalg.norm(axis)
    left = p0 - p1
    right = p3 - p2
    left -= np.dot(left, axis) * axis
    right -= np.dot(right, axis) * axis
    if min(np.linalg.norm(left), np.linalg.norm(right)) < 1e-10:
        raise ValueError("Undefined chi angle with collinear atoms.")
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis, left), right), np.dot(left, right))))


def rotate_about_axis(points, origin, axis, degrees):
    """Rigid Rodrigues rotation: internal distances and the bond axis are preserved."""
    axis = np.asarray(axis, dtype=float)
    length = np.linalg.norm(axis)
    if length < 1e-10:
        raise ValueError("Cannot rotate around a zero-length bond.")
    axis = axis / length
    shifted = np.asarray(points, dtype=float) - origin
    angle = np.deg2rad(degrees)
    rotated = shifted * np.cos(angle) + np.cross(axis, shifted) * np.sin(angle) + np.outer(shifted @ axis, axis) * (1 - np.cos(angle))
    return rotated + origin


def adjust_sidechains(topology, positions, check_cancel=lambda: None, max_residues=200, protected_residues=()):
    from scipy.spatial import cKDTree
    atoms = list(topology.atoms())
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float).copy()
    radii = np.array([{"C": .170, "N": .155, "O": .152, "S": .180, "P": .180, "H": .120}.get(atom.element.symbol if atom.element else "C", .170) for atom in atoms])
    adjacency = [set() for _ in atoms]
    for a, b in topology.bonds():
        adjacency[a.index].add(b.index); adjacency[b.index].add(a.index)
    excluded = []
    for i, bonded in enumerate(adjacency):
        neighbors = set(bonded) | {i}
        for j in bonded:
            neighbors.update(adjacency[j])
        excluded.append(neighbors)
    records, skipped, evaluated = [], [], 0
    for residue in topology.residues():
        definitions = CHI.get(residue.name, [])
        if not definitions:
            continue
        from .complex_topology import residue_key
        if residue_key(residue) in protected_residues:
            skipped.append(f"{residue.chain.id}:{residue.id}:{residue.name} (observed metal coordination)")
            continue
        if evaluated >= max_residues:
            skipped.append(f"{residue.chain.id}:{residue.id}:{residue.name} (residue budget)")
            continue
        evaluated += 1
        named = {atom.name: atom.index for atom in residue.atoms()}
        residue_atoms = {atom.index for atom in residue.atoms()}
        for chi_index, definition in enumerate(definitions, 1):
            check_cancel()
            if not all(name in named for name in definition):
                skipped.append(f"{residue.chain.id}:{residue.id}:chi{chi_index} (missing defining atoms)")
                continue
            indices = [named[name] for name in definition]
            a, b = indices[1:3]
            # Find the complete distal bonded component with the rotatable bond removed.
            moving, pending = set(), [b]
            while pending:
                current = pending.pop()
                if current in moving:
                    continue
                moving.add(current)
                pending.extend(neighbor for neighbor in adjacency[current] if {current, neighbor} != {a, b} and neighbor not in moving)
            if a in moving or not moving.issubset(residue_atoms) or any(atoms[i].name in {"N", "CA", "C", "O", "OXT"} for i in moving):
                skipped.append(f"{residue.chain.id}:{residue.id}:chi{chi_index} (ring/crosslink/backbone constraint)")
                continue
            moving_indices = np.array(sorted(moving), dtype=int)
            tree = cKDTree(xyz)

            def score(candidate):
                total = 0.0
                for local, index in enumerate(moving_indices):
                    neighbors = tree.query_ball_point(candidate[local], .40)
                    for other in neighbors:
                        if other in moving or other in excluded[index]:
                            continue
                        separation = np.linalg.norm(candidate[local] - xyz[other])
                        threshold = .78 * (radii[index] + radii[other])
                        total += max(0.0, (threshold - separation) / threshold) ** 2
                return total

            try:
                original_angle = torsion_degrees(xyz[indices])
            except ValueError:
                continue
            baseline = xyz[moving_indices].copy()
            before = score(baseline)
            best_score, best_angle, best = before, original_angle, baseline
            for target in (-60.0, 60.0, 180.0):
                candidate = rotate_about_axis(baseline, xyz[a], xyz[b] - xyz[a], target - original_angle)
                candidate_score = score(candidate)
                if candidate_score < best_score - 1e-8:
                    best_score, best_angle, best = candidate_score, target, candidate
            if best_score < before - 1e-8:
                xyz[moving_indices] = best
                records.append({"chain": residue.chain.id, "resid": residue.id, "residue": residue.name, "chi": chi_index, "before_degrees": original_angle, "selected_degrees": best_angle, "clash_score_before": before, "clash_score_after": best_score, "moved_atom_indices": moving_indices.tolist()})
    return xyz * unit.nanometer, {"method": "One bounded coordinate-descent pass over standard chi bonds; original and −60°, +60°, 180° states; retain only lower steric-overlap scores. No empirical rotamer probabilities or exhaustive packing.", "steric_score": "Sum of squared fractional overlaps below 0.78 times summed elemental van der Waals radii; bonded and 1–3 pairs excluded. Dimensionless heuristic, not an energy.", "residues_examined": evaluated, "max_residues": max_residues, "adjustments": records, "skipped": skipped, "adjusted_chi_count": len(records)}


def proper_overlay_points(points1, points2):
    """Kabsch alignment restricted to SO(3); never reflect a chiral template.

    Scoped replacement for PDBFixer 1.12's private placement helper, whose
    unconstrained SVD can return determinant −1. Input/output follow its API.
    """
    if not len(points1):
        return np.zeros(3), np.eye(3), np.zeros(3)
    reference, moving = np.asarray(points1, dtype=float), np.asarray(points2, dtype=float)
    center1, center2 = reference.mean(axis=0), moving.mean(axis=0)
    covariance = (moving - center2).T @ (reference - center1)
    left, _, right = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[2, 2] = 1 if np.linalg.det(left @ right) >= 0 else -1
    rotation = (left @ correction @ right).T
    return -center2, rotation, center1


def stereochemistry_report(topology, positions, templates):
    """Validate stereocenter signs against standard residue template geometry."""
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
    centers, violations = [], []
    for residue in topology.residues():
        if residue.name in SUPPORTED_MODIFIED:
            continue
        parent = PARENT_RESIDUES.get(residue.name, residue.name)
        if parent not in STANDARD_PROTEINS or parent == "GLY":
            continue
        definitions = [("CA", ("N", "C", "CB"))]
        if parent == "ILE":
            definitions.append(("CB", ("CA", "CG1", "CG2")))
        elif parent == "THR":
            definitions.append(("CB", ("CA", "OG1", "CG2")))
        named = {atom.name: atom.index for atom in residue.atoms()}
        template = templates[parent]
        template_names = {atom.name: atom.index for atom in template.topology.atoms()}
        template_xyz = np.asarray(template.positions.value_in_unit(unit.nanometer))
        for center, neighbors in definitions:
            if not all(name in named for name in (center, *neighbors)):
                continue
            volume = float(np.linalg.det(np.stack([xyz[named[name]] - xyz[named[center]] for name in neighbors])))
            expected = float(np.linalg.det(np.stack([template_xyz[template_names[name]] - template_xyz[template_names[center]] for name in neighbors])))
            record = {"chain": residue.chain.id, "resid": residue.id, "insertion_code": (residue.insertionCode or '').strip(), "residue": residue.name, "center": center, "signed_volume_nm3": volume, "template_signed_volume_nm3": expected}
            record["atom_indices"] = [named[name] for name in (center, *neighbors)]
            centers.append(record)
            if abs(volume) < 1e-4 or volume * expected <= 0:
                violations.append(record)
    modifications = modified_stereochemistry_report(topology, positions)
    centers.extend(modifications["centers"])
    violations.extend(modifications["violations"])
    return {"method": "Signed N/C/CB volume at nonglycine CA and CA/branch1/branch2 volume at ILE/THR CB versus PDBFixer standard templates; supported modified-residue centers (including TPO CB and HYP CG) versus their CCD ideal geometry. Absolute volume must exceed 0.0001 nm³.", "checked_centers": len(centers), "violations": violations, "centers": centers}


def require_valid_stereochemistry(report, stage):
    if report["violations"]:
        first = report["violations"][0]
        raise ValueError(f"{stage} has inverted or near-planar standard residue stereochemistry at {first['chain']}:{first['resid']} {first['residue']} {first['center']}. The model is rejected; use an externally validated repair instead.")


def restore_modeled_loop_sidechains(topology, positions, templates, observed_keys, environment_positions=()):
    """Correct only wholly new non-Pro loop sidechains in proper local frames.

    PDBFixer's coarse construction potential can invert newly placed centers.
    This single template placement preserves every backbone/observed coordinate;
    it is not loop prediction, packing, or force-field relaxation.
    """
    from scipy.spatial import cKDTree
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float).copy()
    original = xyz.copy()
    atoms = list(topology.atoms())
    before = stereochemistry_report(topology, positions, templates)
    affected = {(v['chain'], v['resid'], v.get('insertion_code', ''), v['residue']) for v in before['violations']}
    report = {"method": "Single rigid proper-rotation placement of complete template sidechains in the modeled N–CA–C frame, anchored at CA; only wholly new standard non-Pro residues with failed stereochemistry are eligible.",
              "confidence": "low", "stereochemistry_before": before, "corrections": [], "errors": [], "accepted": False}

    def reject(message):
        report["errors"].append(message)
        return xyz * unit.nanometer, report

    def frame(points, named):
        n = points[named['N']] - points[named['CA']]
        c = points[named['C']] - points[named['CA']]
        if min(np.linalg.norm(n), np.linalg.norm(c)) < 1e-8:
            raise ValueError("A zero-length modeled N–CA–C frame cannot orient a sidechain.")
        axis = n / np.linalg.norm(n)
        plane = c - axis * np.dot(axis, c)
        if np.linalg.norm(plane) / np.linalg.norm(c) < .1:
            raise ValueError("A near-collinear modeled N–CA–C frame cannot orient a sidechain.")
        plane /= np.linalg.norm(plane)
        return np.column_stack((axis, plane, np.cross(axis, plane)))

    moved = set()
    for residue in topology.residues():
        if (residue.chain.id, residue.id, (residue.insertionCode or '').strip(), residue.name) not in affected:
            continue
        members = list(residue.atoms())
        if any(atom_key(atom) in observed_keys for atom in members):
            return reject("A failed stereocenter belongs to an observed residue; its atoms cannot be repositioned by loop correction.")
        if residue.name not in STANDARD_PROTEINS or residue.name == 'PRO':
            return reject("Loop sidechain correction requires a standard non-Pro template; ring closure or modified stereochemistry will not be guessed.")
        template = templates[residue.name]
        reference = np.asarray(template.positions.value_in_unit(unit.nanometer))
        named = {a.name: a.index for a in members}
        reference_names = {a.name: a.index for a in template.topology.atoms()}
        template_heavy = {a.name for a in template.topology.atoms() if a.element != app.element.hydrogen and a.name != 'OXT'}
        if not template_heavy.issubset(named) or any(a.name not in reference_names for a in members if a.name not in {'OXT'}):
            return reject("A modeled residue does not match the complete standard sidechain template.")
        try:
            rotation = frame(xyz, named) @ frame(reference, reference_names).T
        except ValueError as exc:
            return reject(str(exc))
        selected = [a for a in members if a.name not in {'N', 'CA', 'C', 'O', 'OXT'}]
        for atom in selected:
            xyz[atom.index] = rotation @ (reference[reference_names[atom.name]] - reference[reference_names['CA']]) + xyz[named['CA']]
            moved.add(atom.index)
        report['corrections'].append({"chain": residue.chain.id, "resid": residue.id, "insertion_code": (residue.insertionCode or '').strip(),
                                      "residue": residue.name, "rotation_determinant": float(np.linalg.det(rotation)),
                                      "modeled_atoms": [{"index": a.index, "identity": list(atom_key(a))} for a in selected]})
    if not np.isfinite(xyz).all():
        return reject("Loop correction produced non-finite coordinates.")
    fixed = [a.index for a in atoms if a.index not in moved]
    report['unchanged_atoms_preserved_exactly'] = bool(np.array_equal(xyz[fixed], original[fixed]))
    adjacency = [set() for _ in atoms]
    bond_checks, peptide_checks = [], []
    for a, b in topology.bonds():
        adjacency[a.index].add(b.index); adjacency[b.index].add(a.index)
        length = float(np.linalg.norm(xyz[a.index] - xyz[b.index]))
        if a.residue != b.residue and {a.name, b.name} == {'N', 'C'} and all(x.residue.name in STANDARD_PROTEINS for x in (a, b)):
            peptide_checks.append({"atoms": [list(atom_key(a)), list(atom_key(b))], "distance_nm": length})
            if not .10 <= length <= .22:
                report['errors'].append("A peptide C–N distance falls outside the gross 1.0–2.2 Å construction bounds.")
        if a.index in moved or b.index in moved:
            if a.residue != b.residue:
                return reject("A corrected sidechain has an external covalent connection; template repositioning is unsupported.")
            template = templates[a.residue.name]
            ref = np.asarray(template.positions.value_in_unit(unit.nanometer))
            names = {x.name: x.index for x in template.topology.atoms()}
            expected = float(np.linalg.norm(ref[names[a.name]] - ref[names[b.name]]))
            bond_checks.append({"atoms": [a.index, b.index], "distance_nm": length, "template_distance_nm": expected})
            if abs(length - expected) > 1e-6:
                report['errors'].append("A corrected attachment/internal sidechain bond differs from its rigid template length.")
    report.update(bond_checks=bond_checks, peptide_bond_checks=peptide_checks)
    radii = np.array([{'C': .170, 'N': .155, 'O': .152, 'S': .180, 'P': .180}.get(a.element.symbol if a.element else 'C', .170) for a in atoms])
    # External retained molecules are supplied as (xyz_nm, elemental vdW radius_nm).
    external = list(environment_positions)
    combined = np.concatenate((xyz, np.asarray([x[0] for x in external]).reshape(-1, 3)))
    all_radii = np.concatenate((radii, np.asarray([x[1] for x in external])))
    tree = cKDTree(combined)
    overlaps, seen, gross = [], set(), []
    for index in sorted(moved):
        excluded = {index} | adjacency[index]
        for neighbor in adjacency[index]:
            excluded |= adjacency[neighbor]
        for other in tree.query_ball_point(xyz[index], .5):
            pair = tuple(sorted((index, other)))
            if other in excluded or pair in seen:
                continue
            seen.add(pair)
            distance = float(np.linalg.norm(xyz[index] - combined[other]))
            ratio = distance / float(all_radii[index] + all_radii[other])
            if ratio < .78:
                overlaps.append({"atoms": list(pair), "distance_nm": distance, "vdw_distance_ratio": ratio})
            if distance < .13:
                gross.append(list(pair))
    report['steric_screen'] = {"method": "Gross collisions below 1.3 Å reject; additionally report overlaps below 0.78 times summed elemental van der Waals radii, excluding bonded and 1–3 pairs. This is a geometric screen, not validated packing or an energy.",
                               "external_environment_atoms": len(external), "soft_overlaps": overlaps, "gross_collisions": gross,
                               "soft_overlap_score": float(sum((1-r['vdw_distance_ratio']/.78)**2 for r in overlaps)),
                               "requires_forcefield_minimization": bool(overlaps)}
    if gross:
        report['errors'].append("Corrected modeled sidechains retain a gross nonbonded collision below 1.3 Å.")
    report['stereochemistry_after'] = stereochemistry_report(topology, xyz * unit.nanometer, templates)
    if report['stereochemistry_after']['violations']:
        report['errors'].append("Corrected coordinates still fail the unchanged stereochemistry guard.")
    report['backbone_gaps'] = backbone_gaps(topology, xyz * unit.nanometer)
    if any(gap['structural_break'] for gap in report['backbone_gaps']):
        report['errors'].append("A long backbone connection remains after loop correction.")
    report['accepted'] = not report['errors'] and report['unchanged_atoms_preserved_exactly']
    return xyz * unit.nanometer, report


def protonation_inventory(topology, selected_variants, charges=None):
    atoms = list(topology.atoms())
    adjacency = [set() for _ in atoms]
    for a, b in topology.bonds():
        adjacency[a.index].add(b.index); adjacency[b.index].add(a.index)
    output = []
    for residue in topology.residues():
        if residue.name not in STANDARD_PROTEINS | SUPPORTED_MODIFIED:
            continue
        members = list(residue.atoms())
        inventory = {atom.name: sorted(atoms[j].name for j in adjacency[atom.index] if atoms[j].element == app.element.hydrogen) for atom in members if atom.element in {app.element.nitrogen, app.element.oxygen, app.element.sulfur}}
        returned = selected_variants[residue.index] if residue.index < len(selected_variants) else None
        state = residue.name
        if residue.name == "HIS":
            nd, ne = bool(inventory.get("ND1")), bool(inventory.get("NE2"))
            state = "HIP" if nd and ne else "HID" if nd else "HIE" if ne else "HIN"
        elif residue.name == "ASP":
            state = "ASH" if inventory.get("OD1") or inventory.get("OD2") else "ASP"
        elif residue.name == "GLU":
            state = "GLH" if inventory.get("OE1") or inventory.get("OE2") else "GLU"
        elif residue.name == "LYS":
            state = "LYN" if len(inventory.get("NZ", [])) == 2 else "LYS"
        elif residue.name == "CYS":
            disulfide = any(a.name == "SG" and b.name == "SG" and a.residue != b.residue
                            and residue in (a.residue, b.residue) for a, b in topology.bonds())
            state = "CYX" if disulfide else "CYS" if inventory.get("SG") else "CYM"
        record = {"chain": residue.chain.id, "resid": residue.id, "residue": residue.name, "state": state, "openmm_returned_variant": returned, "bonded_hydrogens": inventory, "hydrogen_count": sum(atom.element == app.element.hydrogen for atom in members)}
        if charges is not None:
            record["forcefield_charge_e"] = float(sum(charges[atom.index] for atom in members))
        output.append(record)
    return output


class PreparationWorker(Worker):
    def prepare(self):
        from pdbfixer import PDBFixer
        import openmm as mm
        import mdtraj as md

        options = self.settings
        random.seed(options["seed"]); np.random.seed(options["seed"])
        platform = mm.Platform.getPlatformByName("CPU")
        platform.setPropertyDefaultValue("Threads", str(config.CPU_THREADS))
        self.update(status="running", stage="Inspecting current structure", message="Inspecting original sequence evidence and current first-frame coordinates.")
        register_topology_definitions()
        fixer = PDBFixer(filename=str(self.folder / "input.pdb"), platform=platform)
        register_fixer_templates(fixer)
        from .preparation import current_fixer, canonicalize_protonation_aliases
        protonation_aliases = canonicalize_protonation_aliases(fixer.topology)
        from .complex_topology import metal_environment, residue_key, subset
        from .ligands import prepare_ligands, SUPPORTED_IONS
        from .residue_identity import protein_residue_keys
        matching, source_path = current_fixer(options["dataset_id"])
        fixer.sequences = matching.sequences
        fixer.sequence_scheme = getattr(matching, "sequence_scheme", {})
        fixer.sequence_scheme_error = getattr(matching, "sequence_scheme_error", None)
        if fixer.sequence_scheme_error:
            raise ValueError(fixer.sequence_scheme_error)
        from .sequence_evidence import recover_observed_insertion_codes
        fixer.recovered_insertion_codes = recover_observed_insertion_codes(fixer)
        input_atoms = list(fixer.topology.atoms())
        input_topology, input_positions = fixer.topology, fixer.positions
        protein_keys = protein_residue_keys(options["dataset_id"], input_topology, input_positions)
        modification_inspection = inspect_modified(input_topology, options["ph"])
        if modification_inspection["blockers"]:
            raise ValueError(" ".join(modification_inspection["blockers"]))
        modified = modification_inspection["residues"]
        environment = metal_environment(input_topology, input_positions)
        preserve_waters = environment["water_keys"] if not options["remove_heterogens"] else set()
        ion_keys = {residue_key(r) for r in input_topology.residues() if r.name.upper() in SUPPORTED_IONS and len(list(r.atoms())) == 1} if not options["remove_heterogens"] else set()
        input_keys = [atom_key(atom) for atom in input_atoms]
        original_backbone = {atom_key(atom) for atom in input_atoms if atom.name in {"N", "CA", "C", "O"} and residue_key(atom.residue) in protein_keys}
        summary, warnings = [], ["Protonation uses OpenMM template/heuristic pH rules, not computed residue pKa values or constant-pH dynamics."]
        if protonation_aliases:
            summary.append("Normalized LYN/CYM protonation aliases to LYS/CYS template names with unchanged heavy-atom identities; the selected pH determines newly assigned hydrogens.")
        warnings.extend(modification_inspection["warnings"])
        if modified:
            summary.append(f"Retained {len(modified)} modified protein residues with compatible named residue templates; covalent modifications are not removed as heterogens.")
        self.update(completed=1)
        self.update(stage="Applying explicit preparation choices")
        retained_keys = set()
        counts = Counter()
        removed_ligand_keys = {key for key, action in options.get("ligand_actions", {}).items() if action == "remove"}
        removed_ligands = [{"key": ":".join(residue_key(residue)), "atom_indices": [atom.index for atom in residue.atoms()],
                            "reason": "Explicit per-residue removal selected by user"}
                           for residue in input_topology.residues() if ":".join(residue_key(residue)) in removed_ligand_keys]
        for atom in input_atoms:
            water = atom.residue.name.upper() in storage.WATERS
            protein = residue_key(atom.residue) in protein_keys
            keep_water = water and (not options["remove_waters"] or residue_key(atom.residue) in preserve_waters)
            reason = "hydrogens" if atom.element == app.element.hydrogen else "water_atoms" if water and not keep_water else "heterogen_atoms" if not protein and not water and (options["remove_heterogens"] or ":".join(residue_key(atom.residue)) in removed_ligand_keys) else None
            if reason:
                counts[reason] += 1
            if protein or keep_water:
                retained_keys.add(residue_key(atom.residue))
        # Repair only the protein and retained water; ligand graphs are restored
        # from authoritative chemical data and parameterized separately below.
        modeller = subset(input_topology, input_positions, retained_keys, remove_hydrogens=True)
        fixer.topology, fixer.positions = modeller.topology, modeller.positions
        summary.append(f"Explicitly removed {counts['hydrogens']} existing hydrogens before pH reassignment, {counts['water_atoms']} water heavy atoms, and {counts['heterogen_atoms']} other heavy atoms.")
        self.update(completed=2, message=summary[-1])
        self.update(stage="Repairing missing heavy atoms and selected loops")
        find_missing_residues_preserving_identity(fixer)
        selected_loops = {}
        rebuilt = []
        chains = list(fixer.topology.chains())
        for (chain_index, position), names in fixer.missingResidues.items():
            terminal = position in {0, len(list(chains[chain_index].residues()))}
            if not terminal and any(name not in STANDARD_PROTEINS for name in names):
                raise ValueError(unsupported_missing_residue_message(chains[chain_index].id, names))
            if options["build_missing_residues"] and not terminal:
                selected_loops[(chain_index, position)] = names
                rebuilt.append({"chain": chains[chain_index].id, "position": position, "residues": names, "count": len(names), "confidence": "low", "method": "PDBFixer template placement and local optimization; no independent structure prediction or validation"})
        fixer.missingResidues = selected_loops
        from .sequence_evidence import restoration_plan, restore_residue_identities
        source_identity_plan = restoration_plan(fixer, selected_loops)
        require_valid_stereochemistry(stereochemistry_report(fixer.topology, fixer.positions, fixer.templates), "Input structure")
        before_heavy = fixer.topology.getNumAtoms()
        observed_coordinates = {atom_key(atom): np.array(fixer.positions[atom.index].value_in_unit(unit.nanometer)) for atom in fixer.topology.atoms()}
        if len(observed_coordinates) != before_heavy:
            raise ValueError("Duplicate observed atom identities prevent unambiguous loop-construction validation.")
        fixer.findMissingAtoms()
        missing_record = [{"chain": residue.chain.id, "resid": residue.id, "residue": residue.name, "atoms": [atom.name for atom in missing]} for residue, missing in fixer.missingAtoms.items()]
        if options["add_missing_atoms"] or selected_loops:
            import pdbfixer.pdbfixer as fixer_implementation
            original_overlay = fixer_implementation._overlayPoints
            fixer_implementation._overlayPoints = proper_overlay_points
            try:
                fixer.addMissingAtoms(seed=options["seed"])
            finally:
                fixer_implementation._overlayPoints = original_overlay
        restored_identities = restore_residue_identities(fixer.topology, source_identity_plan)
        loop_construction = None
        if selected_loops:
            built_xyz = np.asarray(fixer.positions.value_in_unit(unit.nanometer))
            current = {atom_key(a): a.index for a in fixer.topology.atoms()}
            unchanged = all(key in current and np.array_equal(built_xyz[current[key]], position) for key, position in observed_coordinates.items())
            initial_stereo = stereochemistry_report(fixer.topology, fixer.positions, fixer.templates)
            loop_construction = {"native_build_attempts": 1, "seed": options['seed'], "confidence": "low",
                                 "observed_coordinates_preserved_exactly": unchanged,
                                 "initial_stereochemistry": initial_stereo, "accepted": False,
                                 "coordinate_artifacts": {"before": "loop-construction-before.npz", "after": "loop-construction-after.npz", "topology": "loop-construction.pdb"}}
            np.savez_compressed(self.folder / 'loop-construction-before.npz', xyz_nm=built_xyz)
            with (self.folder / 'loop-construction.pdb').open('w') as handle:
                app.PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
            storage.atomic_json(self.folder / 'loop-construction.json', loop_construction)
            if not unchanged:
                raise ValueError("Loop construction changed or lost an observed atom; the model is rejected.")
            if initial_stereo['violations']:
                self.update(message="Checking one template sidechain correction for newly modeled loop stereocenters; observed atoms and backbone stay fixed.")
                external = [(np.asarray(input_positions[a.index].value_in_unit(unit.nanometer)),
                             {'C': .170, 'N': .155, 'O': .152, 'S': .180, 'P': .180}.get(a.element.symbol if a.element else 'C', .170))
                            for a in input_atoms if a.element != app.element.hydrogen and residue_key(a.residue) not in retained_keys
                            and not options['remove_heterogens'] and ':'.join(residue_key(a.residue)) not in removed_ligand_keys
                            and a.residue.name.upper() not in storage.WATERS]
                candidate, correction = restore_modeled_loop_sidechains(fixer.topology, fixer.positions, fixer.templates, set(observed_coordinates), external)
                loop_construction['sidechain_correction'] = correction
                np.savez_compressed(self.folder / 'loop-construction-after.npz', xyz_nm=np.asarray(candidate.value_in_unit(unit.nanometer)))
                storage.atomic_json(self.folder / 'loop-construction.json', loop_construction)
                if not correction['accepted']:
                    raise ValueError("The bounded modeled-loop sidechain correction was rejected: " + ' '.join(correction['errors']))
                fixer.positions = candidate
                warnings.append("Newly modeled loop sidechains required one proper-frame template correction after coarse construction inverted stereocenters. Observed atoms and modeled backbone were preserved exactly. This remains a low-confidence model, not validated loop or rotamer geometry.")
            else:
                np.savez_compressed(self.folder / 'loop-construction-after.npz', xyz_nm=built_xyz)
        require_valid_stereochemistry(stereochemistry_report(fixer.topology, fixer.positions, fixer.templates), "Repaired structure")
        remaining_gaps = [gap for gap in backbone_gaps(fixer.topology, fixer.positions) if gap["structural_break"]]
        if remaining_gaps:
            raise ValueError("A long backbone connection remains after optional repair. The local model is not suitable for force-field preparation; repair this gap externally. " + remaining_gaps[0]["message"])
        if loop_construction is not None:
            loop_construction['accepted'] = True
            loop_construction['artifacts_sha256'] = {name: hashlib.sha256((self.folder / name).read_bytes()).hexdigest() for name in loop_construction['coordinate_artifacts'].values()}
            storage.atomic_json(self.folder / 'loop-construction.json', loop_construction)
        summary.append(f"Added {fixer.topology.getNumAtoms() - before_heavy} heavy/terminal atoms across existing residues and {len(rebuilt)} rebuilt internal segment(s).")
        if rebuilt:
            warnings.append("LOW-CONFIDENCE REBUILT SEGMENTS: " + "; ".join(f"chain {r['chain']}, insertion position {r['position']}, {'-'.join(r['residues'])}" for r in rebuilt) + ". These coordinates are template-derived models, not experimentally resolved loops; validate them before scientific use.")
        self.update(completed=3, message=summary[-1])
        prepared_ligands, ligand_parameters = [], None
        if not options["remove_heterogens"]:
            self.update(stage="Preparing ligand states and force-field parameters", message="Resolving ligand chemistry while preserving bound heavy-atom coordinates.")
            prepared_ligands = prepare_ligands(options["dataset_id"], options["ph"], options["seed"], self.folder / "ligands",
                                               overrides=options.get("ligand_overrides"),
                                               on_progress=lambda message: self.update(message=message), check_cancel=self.check_cancel,
                                               actions=options.get("ligand_actions"),
                                               environment_indices=[atom.index for atom in input_atoms if atom.element != app.element.hydrogen
                                                                    and ":".join(residue_key(atom.residue)) not in removed_ligand_keys
                                                                    and (atom.residue.name.upper() not in storage.WATERS or not options["remove_waters"] or residue_key(atom.residue) in preserve_waters)])
        if removed_ligands:
            summary.append("Explicitly removed selected noncovalent ligand residues: " + ", ".join(record["key"] for record in removed_ligands) + ". Other molecules were retained.")
        if prepared_ligands:
            parameter_files = sorted({ligand.ffxml for ligand in prepared_ligands})
            ligand_parameters = {"forcefield": "GAFF2", "charge_method": "AM1-BCC", "requires_explicit_solvent": True,
                                 "ligands": [ligand.provenance for ligand in prepared_ligands],
                                 "files": [{"path": str(path.relative_to(self.folder)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in parameter_files]}
            for ligand in prepared_ligands:
                warnings.extend(ligand.provenance.get("warnings", []))
            summary.append(f"Prepared {len(prepared_ligands)} ligands with GAFF2 / AM1-BCC; bound heavy atoms retained, chosen formal charges and atom mappings recorded.")
        modeller = app.Modeller(fixer.topology, fixer.positions)
        for ligand in prepared_ligands:
            modeller.add(ligand.topology, ligand.positions)
        if ion_keys:
            ions = subset(input_topology, input_positions, ion_keys)
            modeller.add(ions.topology, ions.positions)
            summary.append(f"Retained {len(ion_keys)} ions and {len(preserve_waters)} observed metal-coordinating waters; coordinating protein sidechains are protected during chi sampling.")
            warnings.append("Ions use the TIP3P-compatible Amber nonbonded model; observed metal contacts are retained as coordinates, not covalent bonds. This model does not validate coordination energetics.")
        fixer.topology, fixer.positions = modeller.topology, modeller.positions
        self.update(stage="Adjusting sidechain chi angles")
        rotamers = {"method": "Skipped by user", "adjustments": [], "adjusted_chi_count": 0}
        if options["optimize_sidechains"]:
            fixer.positions, rotamers = adjust_sidechains(fixer.topology, fixer.positions, self.check_cancel, protected_residues=environment["protected_keys"] if ion_keys else ())
            summary.append(f"Bounded steric chi sampling accepted {rotamers['adjusted_chi_count']} lower-overlap angle changes across {rotamers['residues_examined']} examined residues.")
            warnings.append("Sidechain adjustment is one bounded steric chi-angle search; it does not establish the correct rotamer or perform exhaustive packing. Any subsequent relaxation is recorded separately.")
        self.update(completed=4, message=summary[-1] if options["optimize_sidechains"] else "Sidechain sampling skipped by request.")
        self.update(stage="Assigning pH-dependent hydrogens")
        modeller = app.Modeller(fixer.topology, fixer.positions)
        has_water = any(residue.name.upper() in storage.WATERS for residue in modeller.topology.residues())
        from .prepared_system import load_prepared_forcefield, copy_ligand_parameters, snapshot_modified_parameters
        parameter_state = {"ligand_parameters": ligand_parameters} if ligand_parameters else {}
        if modified:
            parameter_state.update(modified_residues=modified, modified_residue_parameters=modified_forcefield_provenance())
            snapshot_modified_parameters(self.folder, parameter_state)
        ff, files = load_prepared_forcefield(self.folder, parameter_state, solvent="explicit" if has_water or ion_keys or prepared_ligands or modified else "implicit")
        # Existing hydrogens were explicitly removed above. Auto variants alone never remove stale H.
        selected_variants = modeller.addHydrogens(ff, pH=options["ph"], platform=platform)
        system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        verified_ions = []
        if ion_keys:
            from .ions import inspect_ions, validate_ion_system
            expected_ions = [record for record in inspect_ions(input_topology) if record["key"] in {":".join(key) for key in ion_keys}]
            verified_ions = validate_ion_system(modeller.topology, system, expected_ions=expected_ions)
        charges = None
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                charges = [force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())]
                break
        states = protonation_inventory(modeller.topology, selected_variants, charges)
        if modified:
            by_key = {(r.chain.id, r.id, (r.insertionCode or "").strip(), r.name): r for r in modeller.topology.residues()}
            for record in modified:
                residue = by_key.get((record["chain"], record["resid"], record.get("insertion_code", ""), record["residue"]))
                if residue is None or charges is None:
                    raise ValueError("A modified residue lost its identity or charge parameters during preparation.")
                actual_charge = sum(charges[a.index] for a in residue.atoms())
                terminal_charge = -1 if any(a.name == "OXT" for a in residue.atoms()) else 0
                if abs(actual_charge - record["formal_charge"]) > 1e-5:
                    raise ValueError(f"Modified residue {record['residue']} {record['chain']}:{record['resid']} does not retain its expected template charge.")
                record.update(assembled_charge_e=float(actual_charge), terminal_charge_e=terminal_charge)
            (self.folder / "modified-residue-system.xml").write_text(mm.XmlSerializer.serialize(system))
        added_hydrogens = sum(atom.element == app.element.hydrogen for atom in modeller.topology.atoms())
        summary.append(f"Added {added_hydrogens} hydrogens at requested pH {options['ph']:g}; actual residue states and bonded-H inventories are recorded.")
        self.update(completed=5, message=summary[-1])
        self.update(stage="Locally relaxing sidechains with backbone restraints")
        relaxation = {"performed": False, "reason": "Sidechain adjustment not requested"}
        if options["optimize_sidechains"] and not has_water and not prepared_ligands and not ion_keys and not modified:
            restraint = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            restraint.addGlobalParameter("k", 1000 * unit.kilojoule_per_mole / unit.nanometer**2)
            for parameter in ("x0", "y0", "z0"):
                restraint.addPerParticleParameter(parameter)
            xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer))
            restrained = []
            for atom in modeller.topology.atoms():
                if atom_key(atom) in original_backbone:
                    restraint.addParticle(atom.index, xyz[atom.index].tolist()); restrained.append(atom.index)
            system.addForce(restraint)
            integrator = mm.VerletIntegrator(.001 * unit.picosecond)
            context = mm.Context(system, integrator, platform, {"Threads": str(config.CPU_THREADS)})
            context.setPositions(modeller.positions)
            energy_before = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            mm.LocalEnergyMinimizer.minimize(context, 10 * unit.kilojoule_per_mole / unit.nanometer, 150)
            state = context.getState(getPositions=True, getEnergy=True)
            modeller.positions = state.getPositions()
            energy_after = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            if not math.isfinite(energy_after) or not np.isfinite(modeller.positions.value_in_unit(unit.nanometer)).all():
                raise ValueError("Preparation relaxation produced non-finite coordinates or energy.")
            relaxation = {"performed": True, "max_iterations": 150, "tolerance_kj_mol_nm": 10, "backbone_restraint_k_kj_mol_nm2": 1000, "restrained_original_backbone_atoms": len(restrained), "energy_with_restraints_before_kj_mol": energy_before, "energy_with_restraints_after_kj_mol": energy_after, "forcefield_files": files}
            del context, integrator
            summary.append("Performed at most 150 local minimization iterations with restraints on retained original backbone atoms; restraint energy is not an affinity or quality score.")
        elif options["optimize_sidechains"]:
            relaxation = {"performed": False, "reason": "Retained ligands, ions or crystallographic waters require a complete explicit solvent box. Minimize the solvated complex before dynamics."}
            warnings.append("Protein chi sampling included the retained complex. Local implicit-solvent relaxation was skipped; use energy minimization after building explicit water.")
        self.update(completed=6, message="Local relaxation stage complete.")
        self.update(stage="Saving prepared structure and provenance")
        geometry = None
        if prepared_ligands or ion_keys or modified:
            output_named = {atom_key(atom): atom for atom in modeller.topology.atoms()}
            input_xyz = np.asarray(input_positions.value_in_unit(unit.angstrom))
            output_xyz = np.asarray(modeller.positions.value_in_unit(unit.angstrom))
            modified_keys = {residue_key(r) for r in input_topology.residues() if r.name in SUPPORTED_MODIFIED}
            retained = {tuple(ligand.original_residue_key) for ligand in prepared_ligands} | ion_keys | preserve_waters | modified_keys
            donor_keys = {tuple(contact["donor"]) for contact in environment["report"]["contacts"]} if ion_keys else set()
            checked, displacements = [], []
            for atom in input_atoms:
                identity = atom_key(atom)
                if atom.element == app.element.hydrogen or (residue_key(atom.residue) not in retained and identity not in donor_keys):
                    continue
                if identity not in output_named:
                    raise ValueError(f"Complex preparation lost a retained atom: {identity}.")
                displacement = float(np.linalg.norm(input_xyz[atom.index] - output_xyz[output_named[identity].index]))
                checked.append(list(identity)); displacements.append(displacement)
            if displacements and max(displacements) > 1e-6:
                raise ValueError("Complex preparation unexpectedly moved a bound ligand, ion, coordinating water or protected donor atom.")
            geometry = {"checked_atoms": len(checked), "max_displacement_angstrom": max(displacements, default=0), "identities": checked}
            for ligand in prepared_ligands:
                actual_atoms = [atom for atom in modeller.topology.atoms() if residue_key(atom.residue) == tuple(ligand.original_residue_key)]
                expected_charge = ligand.provenance["formal_charge"]
                actual_charge = sum(charges[atom.index] for atom in actual_atoms)
                if abs(actual_charge - expected_charge) > 1e-4 or len(actual_atoms) != ligand.topology.getNumAtoms():
                    raise ValueError(f"Prepared ligand {ligand.provenance['key']} changed atom count or charge during protein assembly.")
                # The final template match is checked, not merely the input XML.
                expected_atoms = list(ligand.topology.atoms())
                for expected_atom, expected in zip(expected_atoms, ligand.provenance["charges_e"]):
                    actual_atom = output_named.get(atom_key(expected_atom))
                    if actual_atom is None or abs(charges[actual_atom.index] - expected) > 1e-6:
                        raise ValueError(f"Prepared ligand {ligand.provenance['key']} did not retain its per-atom charge assignment.")
                maximum = max(float(np.linalg.norm(input_xyz[row["original_index"]] - output_xyz[output_named[(*ligand.original_residue_key, row["prepared_name"])].index])) for row in ligand.atom_map)
                ligand.provenance.update(assembled_forcefield_charge_e=float(actual_charge), assembled_atoms=len(actual_atoms), heavy_coordinate_max_displacement_angstrom=maximum)
        stereochemistry = stereochemistry_report(modeller.topology, modeller.positions, fixer.templates)
        require_valid_stereochemistry(stereochemistry, "Final prepared structure")
        with (self.folder / "prepared.pdb").open("w") as output:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, output, keepIds=True)
        output_atoms = list(modeller.topology.atoms())
        output_keys = {atom_key(atom): atom.index for atom in output_atoms}
        mapping = [{"input_index": i, "output_index": output_keys.get(key), "identity": list(key)} for i, key in enumerate(input_keys)]
        old_keys = set(input_keys)
        added = [{"output_index": atom.index, "identity": list(atom_key(atom))} for atom in output_atoms if atom_key(atom) not in old_keys]
        preparation = {"parent_dataset_id": options["dataset_id"], "ph": options["ph"], "method": "PDBFixer heavy-atom repair + optional short sequence-supported loops + bounded chi search + OpenMM hydrogen/template assignment", "summary": summary, "warnings": warnings, "seed": options["seed"], "forcefield_files": files, "simulation_ready": True, "exact_topology_file": "prepared.pdb", "protonation_states": states, "selected_variants": selected_variants, "removed_atom_counts": dict(counts), "rebuilt_segments": rebuilt, "repaired_atoms": missing_record, "sidechain_adjustment": rotamers, "relaxation": relaxation, "original_to_prepared_atom_map_file": "atom-map.json", "net_forcefield_charge_e": float(sum(charges)) if charges is not None else None, "stereochemistry": stereochemistry, "template_placement": "Scoped proper-rotation Kabsch alignment (det R=+1) replaces PDBFixer 1.12 reflection-capable _overlayPoints during heavy-atom/loop placement."}
        preparation.update(job_id=self.job["id"], metal_environment=environment["report"] if ion_keys else None, preserved_bound_geometry=geometry)
        if loop_construction is not None:
            preparation['loop_construction'] = loop_construction
            needs_minimization = loop_construction.get('sidechain_correction', {}).get('steric_screen', {}).get('requires_forcefield_minimization', False)
            preparation['requires_minimization'] = bool(needs_minimization and not relaxation.get('performed'))
            if preparation['requires_minimization']:
                warnings.append("Modeled loop sidechains retain soft steric overlaps and local complex relaxation was skipped. Force-field minimization is required before dynamics; steric packing and the loop pose are not validated.")
        if protonation_aliases:
            preparation["input_protonation_aliases"] = protonation_aliases
        if ion_keys:
            preparation.update(ions=verified_ions, requires_explicit_solvent=True)
        if modified:
            preparation.update(modified_residues=modified, modified_residue_parameters=modified_forcefield_provenance(), requires_explicit_solvent=True)
        if ligand_parameters:
            preparation["ligand_parameters"] = ligand_parameters
        preparation["ligand_actions"] = options.get("ligand_actions", {})
        if source_identity_plan is not None:
            preparation["sequence_mapping"] = {"method": "Exact mmCIF _pdbx_poly_seq_scheme positions; author numbering gaps are not interpreted as missing sequence", "restored_inserted_residues": restored_identities, "recovered_observed_insertion_codes": fixer.recovered_insertion_codes, "observed_residue_identities_preserved": True}
        preparation["removed_ligand_residues"] = removed_ligands
        storage.atomic_json(self.folder / "preparation.json", preparation)
        storage.atomic_json(self.folder / "atom-map.json", {"input_to_output": mapping, "added_atoms": added, "note": "Original hydrogens were removed and regenerated; matching H names do not imply retained hydrogen coordinates."})
        self.provenance.update(purpose="Exploratory protein–ligand preparation with recorded state assumptions; no claim of native rotamers, loops, binding affinity or convergence", operation="prepare", versions={**self.provenance["versions"], "pdbfixer": importlib.metadata.version("pdbfixer")}, preparation_state=preparation, source_sequence_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest() if source_path else None, preparation_worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), parent_dataset_id=options["dataset_id"])
        self.provenance["outputs"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in self.folder.iterdir() if path.is_file() and (path.name in {"input.pdb", "prepared.pdb", "preparation.json", "atom-map.json", "config.json", "modified-residue-system.xml"} or path.name.startswith('loop-construction'))}
        storage.atomic_json(self.folder / "provenance.json", self.provenance)
        traj = md.Trajectory(np.array([modeller.positions.value_in_unit(unit.nanometer)], dtype=np.float32), md.Topology.from_openmm(modeller.topology), time=[0])
        # Crystal lattice records are not advertised as a prepared periodic solvent box.
        metadata = storage.save_dataset(traj, options["name"], "PDBFixer + OpenMM + GAFF2 preparation" if prepared_ligands else "PDBFixer + OpenMM preparation", f"Prepared {'protein–ligand complex' if prepared_ligands else 'protein with modified residues' if modified else 'standard protein'} at requested pH {options['ph']:g} using recorded fixed states. " + " ".join(summary), warnings=warnings, provenance=self.provenance)
        folder = storage.dataset_dir(metadata["id"])
        for name in ("prepared.pdb", "preparation.json", "atom-map.json"):
            shutil.copy2(self.folder / name, folder / name)
        if loop_construction is not None:
            for name in ('loop-construction.json', *loop_construction['coordinate_artifacts'].values()):
                shutil.copy2(self.folder / name, folder / name)
        if modified:
            shutil.copy2(self.folder / "modified-residue-system.xml", folder / "modified-residue-system.xml")
        copy_ligand_parameters(self.folder, folder, preparation)
        if source_path:
            shutil.copy2(source_path, folder / ("sequence-source" + source_path.suffix.lower()))
        metadata.update(preparation=preparation, parent_dataset_id=options["dataset_id"])
        storage.atomic_json(folder / "metadata.json", metadata)
        self.update(completed=7, status="completed", stage="Prepared complex ready" if prepared_ligands else "Prepared protein ready", dataset_id=metadata["id"], message=f"Saved prepared dataset {metadata['id']}: {metadata['n_atoms']:,} atoms. Protonation and modeling caveats retained.")

    def run(self):
        try:
            if self.settings["operation"] == "solvate":
                from .solvent import solvate_dataset
                self.update(status="running", stage="Building explicit TIP3P solvent", completed=0, message="Solvating the exact prepared topology; preserving its hydrogens and chosen protonation.")
                result = solvate_dataset(self.settings["dataset_id"], padding_nm=self.settings["padding_nm"], seed=self.settings["seed"], ph=self.settings["ph"])
                self.update(completed=1, stage="Saving solvent preview")
                folder = storage.dataset_dir(result["id"])
                for name in ("prepared.pdb", "provenance.json"):
                    if (folder / name).exists():
                        shutil.copy2(folder / name, self.folder / name)
                from .prepared_system import copy_ligand_parameters
                copy_ligand_parameters(folder, self.folder, result.get("preparation"))
                self.update(completed=2, status="completed", stage="Explicit solvent ready", dataset_id=result["id"], message=f"Real solvent preview ready: {result['n_atoms']:,} atoms. The prepared periodic box will be reused by OpenMM.")
            else:
                self.prepare()
        except (Cancelled, KeyboardInterrupt):
            self.job.update(status="cancelled", stage="Cancelled", error="Cancelled by user. Partial preparation outputs retained.", elapsed_seconds=round(time.monotonic() - self.started, 2))
            storage.atomic_json(self.folder / "status.json", self.job)
        except Exception as exc:
            traceback.print_exc()
            self.job.update(status="failed", stage="Preparation failed", error=str(exc)[-6000:], elapsed_seconds=round(time.monotonic() - self.started, 2))
            self.job["logs"] = (self.job["logs"] + ["ERROR: " + str(exc)[-1500:]])[-150:]
            storage.atomic_json(self.folder / "status.json", self.job)
            self.provenance["failure"] = str(exc)
            storage.atomic_json(self.folder / "provenance.json", self.provenance)


def main():
    def terminate(signum, frame):
        raise Cancelled("Cancelled by user")
    signal.signal(signal.SIGTERM, terminate)
    PreparationWorker(sys.argv[1]).run()


if __name__ == "__main__":
    main()
