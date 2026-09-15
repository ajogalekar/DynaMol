"""Bounded native full-complex minimization without positional restraints.

Research preparation experiment only. Source coordinates are preserved as
inputs; all movement, geometry checks and force-field hashes are recorded.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit
import parmed
import psutil

from check_complex_loops import LoopMonitor
from native_elements import inventory
from prepare_adduct import digest


def run(inputs, protein, output):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    process = psutil.Process()
    source_paths = [inputs/n for n in ['solvated.prmtop', 'solvated.inpcrd', 'source-mapping.json', 'result.json']]
    source_hashes = {str(p): digest(p) for p in source_paths}
    evidence = json.loads((inputs/'result.json').read_text())
    inventory(inputs/'solvated.prmtop')
    topology = app.AmberPrmtopFile(str(inputs/'solvated.prmtop'))
    coordinates = app.AmberInpcrdFile(str(inputs/'solvated.inpcrd'))
    native = parmed.load_file(str(inputs/'solvated.prmtop'))
    system = topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*unit.nanometer,
        constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=1e-5)
    allowed = {'HarmonicBondForce', 'HarmonicAngleForce', 'PeriodicTorsionForce', 'NonbondedForce', 'CMMotionRemover'}
    classes = [type(f).__name__ for f in system.getForces()]
    if set(classes)-allowed:
        raise ValueError('Unexpected force class in the native model')
    monitor = LoopMonitor(inputs, protein, topology.topology)
    original = np.asarray(coordinates.positions.value_in_unit(unit.nanometer))
    cov = next(r for r in native.residues if r.name == 'COV')
    centers = [(a.idx, sorted(b.idx for b in a.bond_partners)[:3]) for a in cov.atoms
               if a.atomic_number == 6 and len(a.bond_partners) == 4]
    volumes = lambda xyz: np.array([np.linalg.det(xyz[neighbors]-xyz[center]) for center, neighbors in centers])
    signs = volumes(original)
    if np.min(abs(signs)) < 1e-4:
        raise ValueError('Initial adduct tetrahedral carbon is nearly planar')
    integrator = mm.VerletIntegrator(.002*unit.picosecond)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName('CPU'),
                         {'Threads': '2', 'DeterministicForces': 'true'})
    context.setPeriodicBoxVectors(*coordinates.boxVectors)
    context.setPositions(coordinates.positions)
    monitor.check(coordinates.positions, output/'loop-checks', 'initial', 0)
    initial_state = context.getState(getEnergy=True)
    before = initial_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not np.isfinite(before):
        raise ValueError('Nonfinite initial energy')
    (output/'progress.json').write_text(json.dumps({'stage': 'unrestrained native minimization',
        'pid': process.pid, 'initial_energy_kj_mol': before, 'max_iterations': 1500,
        'positional_or_construction_restraints': False}, indent=2)+'\n')
    mm.LocalEnergyMinimizer.minimize(context, 10*unit.kilojoule_per_mole/unit.nanometer, 1500)
    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True, getParameters=True)
    xyz = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    forces = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
    after = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not (np.isfinite(xyz).all() and np.isfinite(forces).all() and np.isfinite(after)) or after > before:
        raise ValueError('Native minimization produced invalid coordinates, forces or energy')
    (output/'state.xml').write_text(mm.XmlSerializer.serialize(state))
    (output/'system.xml').write_text(mm.XmlSerializer.serialize(system))
    topology.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    with (output/'minimized.pdb').open('w') as handle:
        app.PDBFile.writeFile(topology.topology, state.getPositions(), handle, keepIds=True)
    edges = [(b.atom1.idx, b.atom2.idx, b.type.req*.1) for b in native.bonds
             if b.atom1.residue is cov and b.atom2.residue is cov]
    ratios = [np.linalg.norm(xyz[a]-xyz[b])/length for a, b, length in edges]
    if min(ratios) < .6 or max(ratios) > 1.5 or np.any(volumes(xyz)*signs <= 0):
        raise ValueError('Adduct geometry did not survive native minimization')
    loop_result = monitor.check(state.getPositions(), output/'loop-checks', 'minimized', 0)
    mapping = json.loads((inputs/'source-mapping.json').read_text())
    original_heavy = [row['native_index'] for row in mapping if native.atoms[row['native_index']].atomic_number > 1]
    displacements = np.linalg.norm(xyz[original_heavy]-original[original_heavy], axis=1)*10
    if source_hashes != {str(p): digest(p) for p in source_paths}:
        raise ValueError('Native source changed during minimization')
    report = {'stage': 'native_unrestrained_minimization_checked', 'case': evidence['case'],
        'app_ready': False, 'physical_model_validated': False, 'dynamics_completed': False,
        'positional_or_construction_restraints': False, 'force_classes': classes,
        'max_iterations': 1500, 'tolerance_kj_mol_nm': 10, 'native_threads': 2,
        'initial_energy_kj_mol': float(before), 'final_energy_kj_mol': float(after),
        'maximum_force_kj_mol_nm': float(np.linalg.norm(forces, axis=1).max()),
        'retained_input_heavy_atom_displacement_angstrom': {'maximum': float(displacements.max()),
            'rms': float(np.sqrt(np.mean(displacements**2)))},
        'adduct_bond_ratio_range': [float(min(ratios)), float(max(ratios))],
        'adduct_tetrahedral_carbon_signs_preserved': True, 'loop_checks': loop_result,
        'source_sha256': source_hashes, 'outputs_sha256': {n: digest(output/n) for n in ['state.xml', 'system.xml', 'minimized.pdb']},
        'elapsed_seconds': time.monotonic()-started,
        'scope': 'A complete native system relaxed before defining any warm-up reference. Source coordinates are unchanged on disk; this MD-minimization experiment is not a claim of native loop accuracy or satisfaction of the earlier local observed-context movement cap.'}
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['input', 'protein', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    try:
        run(args.input, args.protein, args.output)
    except Exception as exc:
        if args.output.is_dir():
            (args.output/'failure.json').write_text(json.dumps({'error': str(exc), 'type': type(exc).__name__}, indent=2)+'\n')
        raise
