"""Bounded OpenMM relaxation of modeled loops with observed heavy atoms fixed.

This creates starting geometry before solvation/dynamics. A finite, lower
energy is not acceptance: independent stereo and geometry gates must pass.
"""
from __future__ import annotations

import time
import numpy as np
import openmm as mm
from openmm import app, unit

from . import config
from .residue_identity import residue_key


def minimize_modeled_loops(topology, positions, forcefield, modeled_keys, *, max_iterations=1000):
    """Return (positions, report); move loop atoms and regenerated H only.

    Caller must have removed input hydrogens and regenerated them during prep.
    """
    if not 1 <= max_iterations <= 1000:
        raise ValueError("Loop refinement is limited to 1–1000 minimization iterations.")
    started = time.monotonic()
    keys = {tuple(key) for key in modeled_keys}
    atoms = list(topology.atoms())
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
    if xyz.shape != (len(atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Loop refinement needs finite coordinates matching its topology.")
    selected = [residue for residue in topology.residues() if residue_key(residue) in keys]
    found = [residue_key(residue) for residue in selected]
    if not keys or set(found) != keys or len(found) != len(set(found)):
        raise ValueError("Loop refinement requires exact modeled residue identities.")
    if any(not list(residue.atoms()) for residue in selected):
        raise ValueError("Every modeled loop residue must contain atoms.")
    # Recreate without constraints: constraints cannot couple massless observed
    # heavy atoms to mobile hydrogens. The same parameterized ForceField is used.
    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff,
                                    constraints=None, rigidWater=False)
    # Temporary construction restraints enforce peptide planarity while the
    # loop adapts to fixed crystal anchors. New non-Pro links start trans;
    # Pro links retain the candidate's cis/trans basin. These restraints are
    # confined to this temporary System, never exported to the MD system.
    from .loop_geometry import _dihedral
    from .modified_residues import PARENT_RESIDUES
    peptide_restraints = mm.CustomTorsionForce("k*(1-cos(theta-target))")
    peptide_restraints.addPerTorsionParameter("k")
    peptide_restraints.addPerTorsionParameter("target")
    peptide_targets = []
    for a, b in topology.bonds():
        if a.residue == b.residue or {a.name, b.name} != {"C", "N"}:
            continue
        left, right = (a, b) if a.name == "C" else (b, a)
        if not ({residue_key(left.residue), residue_key(right.residue)} & keys):
            continue
        before = {atom.name: atom for atom in left.residue.atoms()}
        after = {atom.name: atom for atom in right.residue.atoms()}
        if "CA" not in before or "CA" not in after:
            raise ValueError("Modeled peptide restraints need both alpha-carbon anchors.")
        indices = [before["CA"].index, left.index, right.index, after["CA"].index]
        omega = _dihedral(xyz[indices])
        if omega is None:
            raise ValueError("Cannot initialize a peptide restraint from a degenerate loop.")
        is_proline = PARENT_RESIDUES.get(right.residue.name, right.residue.name) == "PRO"
        target = 0.0 if is_proline and abs(omega) < 90 else np.pi
        peptide_restraints.addTorsion(*indices, [200.0, target])
        peptide_targets.append({"atoms": indices, "target_degrees": float(np.degrees(target)),
                                "initial_omega_degrees": omega,
                                "modeled_bond_before_residue": list(residue_key(right.residue))})
    if peptide_targets:
        system.addForce(peptide_restraints)
    fixed, mobile = [], []
    for atom in atoms:
        if residue_key(atom.residue) in keys or atom.element == app.element.hydrogen:
            mobile.append(atom.index)
        else:
            system.setParticleMass(atom.index, 0 * unit.dalton)
            fixed.append(atom.index)
    if not mobile:
        raise ValueError("Loop refinement requires mobile modeled atoms.")
    platform = mm.Platform.getPlatformByName("CPU")
    integrator = mm.VerletIntegrator(.001 * unit.picosecond)
    context = mm.Context(system, integrator, platform, {"Threads": str(config.CPU_THREADS)})
    try:
        context.setPositions(positions)
        before = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(before):
            raise ValueError("Loop refinement starts with non-finite energy.")
        mm.LocalEnergyMinimizer.minimize(context, 10 * unit.kilojoule_per_mole / unit.nanometer,
                                        max_iterations)
        state = context.getState(getPositions=True, getEnergy=True, getForces=True)
        after = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        refined = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        forces = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))
        unchanged = bool(np.array_equal(refined[fixed], xyz[fixed]))
        if not (np.isfinite(after) and np.isfinite(refined).all() and np.isfinite(forces).all()):
            raise ValueError("Loop refinement produced non-finite coordinates, forces or energy.")
        if not unchanged:
            raise ValueError("Loop refinement moved an observed heavy atom; the candidate is rejected.")
        if after > before + max(1e-3, abs(before) * 1e-8):
            raise ValueError("Bounded loop refinement increased the potential energy; the candidate is rejected.")
        report = {"method": "OpenMM LocalEnergyMinimizer with the prepared force field and temporary peptide construction restraints, NoCutoff, no constraints; all non-loop heavy atoms fixed by zero mass. Loop atoms and regenerated hydrogens mobile.",
                  "max_iterations": max_iterations, "tolerance_kj_mol_nm": 10,
                  "fixed_heavy_atoms": len(fixed), "mobile_atoms": len(mobile),
                  "fixed_heavy_coordinates_preserved_exactly": unchanged,
                  "energy_before_kj_mol": float(before), "energy_after_kj_mol": float(after),
                  "energy_includes_construction_restraints": bool(peptide_targets),
                  "mobile_force_rms_kj_mol_nm": float(np.sqrt(np.mean(forces[mobile] ** 2))),
                  "elapsed_seconds": time.monotonic() - started,
                  "peptide_construction_restraints": {"strength_kj_mol": 200, "targets": peptide_targets,
                      "exported_to_simulation": False,
                      "assumption": "New non-Pro peptide bonds initialized trans; Pro bonds retain the candidate cis/trans basin. Native isomer states of missing residues are unknown."},
                  "scope": "Local starting-geometry relaxation before periodic solvent equilibration; no convergence or native-loop claim. Stereo and geometric acceptance are checked separately."}
        return refined * unit.nanometer, report
    finally:
        del context, integrator
