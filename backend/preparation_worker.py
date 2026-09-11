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
from .preparation import STANDARD_PROTEINS, backbone_gaps
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
    return (residue.chain.id, residue.id, residue.insertionCode, residue.name, atom.name)


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


def adjust_sidechains(topology, positions, check_cancel=lambda: None, max_residues=200):
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
        if residue.name not in STANDARD_PROTEINS or residue.name == "GLY":
            continue
        definitions = [("CA", ("N", "C", "CB"))]
        if residue.name == "ILE":
            definitions.append(("CB", ("CA", "CG1", "CG2")))
        elif residue.name == "THR":
            definitions.append(("CB", ("CA", "OG1", "CG2")))
        named = {atom.name: atom.index for atom in residue.atoms()}
        template = templates[residue.name]
        template_names = {atom.name: atom.index for atom in template.topology.atoms()}
        template_xyz = np.asarray(template.positions.value_in_unit(unit.nanometer))
        for center, neighbors in definitions:
            if not all(name in named for name in (center, *neighbors)):
                continue
            volume = float(np.linalg.det(np.stack([xyz[named[name]] - xyz[named[center]] for name in neighbors])))
            expected = float(np.linalg.det(np.stack([template_xyz[template_names[name]] - template_xyz[template_names[center]] for name in neighbors])))
            record = {"chain": residue.chain.id, "resid": residue.id, "residue": residue.name, "center": center, "signed_volume_nm3": volume, "template_signed_volume_nm3": expected}
            centers.append(record)
            if abs(volume) < 1e-4 or volume * expected <= 0:
                violations.append(record)
    return {"method": "Signed N/C/CB volume at nonglycine CA and CA/branch1/branch2 volume at ILE/THR CB versus PDBFixer standard templates; absolute volume must exceed 0.0001 nm³.", "checked_centers": len(centers), "violations": violations, "centers": centers}


def require_valid_stereochemistry(report, stage):
    if report["violations"]:
        first = report["violations"][0]
        raise ValueError(f"{stage} has inverted or near-planar standard residue stereochemistry at {first['chain']}:{first['resid']} {first['residue']} {first['center']}. The model is rejected; use an externally validated repair instead.")


def protonation_inventory(topology, selected_variants, charges=None):
    atoms = list(topology.atoms())
    adjacency = [set() for _ in atoms]
    for a, b in topology.bonds():
        adjacency[a.index].add(b.index); adjacency[b.index].add(a.index)
    output = []
    for residue in topology.residues():
        if residue.name not in STANDARD_PROTEINS:
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
            state = "CYS" if inventory.get("SG") else "CYX"
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
        fixer = PDBFixer(filename=str(self.folder / "input.pdb"), platform=platform)
        from .preparation import current_fixer
        matching, source_path = current_fixer(options["dataset_id"])
        fixer.sequences = matching.sequences
        input_atoms = list(fixer.topology.atoms())
        input_keys = [atom_key(atom) for atom in input_atoms]
        original_backbone = {atom_key(atom) for atom in input_atoms if atom.name in {"N", "CA", "C", "O"} and atom.residue.name in STANDARD_PROTEINS}
        summary, warnings = [], ["Protonation uses OpenMM template/heuristic pH rules, not computed residue pKa values or constant-pH dynamics."]
        self.update(completed=1)
        self.update(stage="Applying explicit preparation choices")
        remove = []
        counts = Counter()
        for atom in input_atoms:
            water = atom.residue.name.upper() in storage.WATERS
            protein = atom.residue.name in STANDARD_PROTEINS
            reason = "hydrogens" if atom.element == app.element.hydrogen else "water_atoms" if water and options["remove_waters"] else "heterogen_atoms" if not protein and not water and options["remove_heterogens"] else None
            if reason:
                remove.append(atom); counts[reason] += 1
        modeller = app.Modeller(fixer.topology, fixer.positions)
        modeller.delete(remove)
        fixer.topology, fixer.positions = modeller.topology, modeller.positions
        summary.append(f"Explicitly removed {counts['hydrogens']} existing hydrogens before pH reassignment, {counts['water_atoms']} water heavy atoms, and {counts['heterogen_atoms']} other heavy atoms.")
        self.update(completed=2, message=summary[-1])
        self.update(stage="Repairing missing heavy atoms and selected loops")
        fixer.findMissingResidues()
        selected_loops = {}
        rebuilt = []
        chains = list(fixer.topology.chains())
        for (chain_index, position), names in fixer.missingResidues.items():
            terminal = position in {0, len(list(chains[chain_index].residues()))}
            if options["build_missing_residues"] and not terminal:
                selected_loops[(chain_index, position)] = names
                rebuilt.append({"chain": chains[chain_index].id, "position": position, "residues": names, "count": len(names), "confidence": "low", "method": "PDBFixer template placement and local optimization; no independent structure prediction or validation"})
        fixer.missingResidues = selected_loops
        require_valid_stereochemistry(stereochemistry_report(fixer.topology, fixer.positions, fixer.templates), "Input structure")
        before_heavy = fixer.topology.getNumAtoms()
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
        require_valid_stereochemistry(stereochemistry_report(fixer.topology, fixer.positions, fixer.templates), "Repaired structure")
        remaining_gaps = [gap for gap in backbone_gaps(fixer.topology, fixer.positions) if gap["structural_break"]]
        if remaining_gaps:
            raise ValueError("A long backbone connection remains after optional repair. The local model is not suitable for force-field preparation; repair this gap externally. " + remaining_gaps[0]["message"])
        summary.append(f"Added {fixer.topology.getNumAtoms() - before_heavy} heavy/terminal atoms across existing residues and {len(rebuilt)} rebuilt internal segment(s).")
        if rebuilt:
            warnings.append("LOW-CONFIDENCE REBUILT SEGMENTS: " + "; ".join(f"chain {r['chain']}, insertion position {r['position']}, {'-'.join(r['residues'])}" for r in rebuilt) + ". These coordinates are template-derived models, not experimentally resolved loops; validate them before scientific use.")
        self.update(completed=3, message=summary[-1])
        self.update(stage="Adjusting sidechain chi angles")
        rotamers = {"method": "Skipped by user", "adjustments": [], "adjusted_chi_count": 0}
        if options["optimize_sidechains"]:
            fixer.positions, rotamers = adjust_sidechains(fixer.topology, fixer.positions, self.check_cancel)
            summary.append(f"Bounded steric chi sampling accepted {rotamers['adjusted_chi_count']} lower-overlap angle changes across {rotamers['residues_examined']} examined residues.")
            warnings.append("Sidechain adjustment is one bounded steric chi-angle search followed by local restrained relaxation; it does not establish the correct rotamer or perform exhaustive packing.")
        self.update(completed=4, message=summary[-1] if options["optimize_sidechains"] else "Sidechain sampling skipped by request.")
        self.update(stage="Assigning pH-dependent hydrogens")
        modeller = app.Modeller(fixer.topology, fixer.positions)
        has_water = any(residue.name.upper() in storage.WATERS for residue in modeller.topology.residues())
        files = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"] if has_water else ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml"]
        ff = app.ForceField(*files)
        # Existing hydrogens were explicitly removed above. Auto variants alone never remove stale H.
        selected_variants = modeller.addHydrogens(ff, pH=options["ph"], platform=platform)
        system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        charges = None
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                charges = [force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())]
                break
        states = protonation_inventory(modeller.topology, selected_variants, charges)
        added_hydrogens = sum(atom.element == app.element.hydrogen for atom in modeller.topology.atoms())
        summary.append(f"Added {added_hydrogens} hydrogens at requested pH {options['ph']:g}; actual residue states and bonded-H inventories are recorded.")
        self.update(completed=5, message=summary[-1])
        self.update(stage="Locally relaxing sidechains with backbone restraints")
        relaxation = {"performed": False, "reason": "Sidechain adjustment not requested"}
        if options["optimize_sidechains"] and not has_water:
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
        elif options["optimize_sidechains"] and has_water:
            relaxation = {"performed": False, "reason": "Retained crystallographic waters: avoid treating a partially hydrated structure as a validated solvent environment."}
            warnings.append("With retained crystallographic waters, discrete steric sidechain adjustment was performed but subsequent protein-only implicit-solvent relaxation was skipped.")
        self.update(completed=6, message="Local relaxation stage complete.")
        self.update(stage="Saving prepared structure and provenance")
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
        storage.atomic_json(self.folder / "preparation.json", preparation)
        storage.atomic_json(self.folder / "atom-map.json", {"input_to_output": mapping, "added_atoms": added, "note": "Original hydrogens were removed and regenerated; matching H names do not imply retained hydrogen coordinates."})
        self.provenance.update(purpose="Exploratory protein preparation, not a validated protonation, rotamer or loop prediction", operation="prepare", versions={**self.provenance["versions"], "pdbfixer": importlib.metadata.version("pdbfixer")}, preparation_state=preparation, source_sequence_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest() if source_path else None, preparation_worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), parent_dataset_id=options["dataset_id"])
        self.provenance["outputs"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in self.folder.iterdir() if path.is_file() and path.name in {"input.pdb", "prepared.pdb", "preparation.json", "atom-map.json", "config.json"}}
        storage.atomic_json(self.folder / "provenance.json", self.provenance)
        traj = md.Trajectory(np.array([modeller.positions.value_in_unit(unit.nanometer)], dtype=np.float32), md.Topology.from_openmm(modeller.topology), time=[0])
        # Crystal lattice records are not advertised as a prepared periodic solvent box.
        metadata = storage.save_dataset(traj, options["name"], "PDBFixer + OpenMM preparation", f"Prepared standard protein at template-assigned pH {options['ph']:g}. " + " ".join(summary), warnings=warnings, provenance=self.provenance)
        folder = storage.dataset_dir(metadata["id"])
        for name in ("prepared.pdb", "preparation.json", "atom-map.json"):
            shutil.copy2(self.folder / name, folder / name)
        if source_path:
            shutil.copy2(source_path, folder / ("sequence-source" + source_path.suffix.lower()))
        metadata.update(preparation=preparation, parent_dataset_id=options["dataset_id"])
        storage.atomic_json(folder / "metadata.json", metadata)
        self.update(completed=7, status="completed", stage="Prepared protein ready", dataset_id=metadata["id"], message=f"Saved prepared dataset {metadata['id']}: {metadata['n_atoms']:,} atoms. Protonation and modeling caveats retained.")

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
