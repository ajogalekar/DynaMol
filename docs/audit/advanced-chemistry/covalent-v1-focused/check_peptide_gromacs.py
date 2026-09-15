"""Compare actual GROMACS and Amber-reader mechanics for capped drug adducts.

This is a parameter-transport check, not physical force-field validation.
Inputs are immutable; each invocation requires a fresh output directory.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import mdtraj
import numpy as np
import openmm as mm
from openmm import app, unit
import parmed

from prepare_adduct import digest


MDP = """integrator = md
nsteps = 0
dt = 0.001
nstxout = 1
nstvout = 0
nstfout = 1
nstenergy = 1
nstlog = 1
cutoff-scheme = Verlet
nstlist = 1
verlet-buffer-tolerance = -1
rlist = 9.0
coulombtype = Cut-off
coulomb-modifier = None
rcoulomb = 9.0
vdwtype = Cut-off
vdw-modifier = None
rvdw = 9.0
DispCorr = no
pbc = xyz
constraints = none
tcoupl = no
pcoupl = no
comm-mode = none
gen-vel = no
continuation = yes
"""


def compare(source, output, executable):
    output.mkdir(exist_ok=False)
    start = time.monotonic()
    source_hashes = {name: digest(source/name) for name in ['peptide.prmtop', 'peptide.inpcrd']}
    native = parmed.load_file(str(source/'peptide.prmtop'), xyz=str(source/'peptide.inpcrd'))
    native.box = [200, 200, 200, 90, 90, 90]
    native.save(str(output/'system.top'))
    (output/'static.mdp').write_text(MDP)
    coordinate = np.asarray(native.coordinates)*.1
    coordinate += 10-coordinate.mean(axis=0)
    if np.max(np.linalg.norm(coordinate[:, None]-coordinate[None, :], axis=2)) >= 8:
        raise ValueError('Peptide too large for the isolated all-pairs comparison')
    topology = app.GromacsTopFile(str(output/'system.top'),
        periodicBoxVectors=tuple(mm.Vec3(*(np.eye(3)[i]*20)) for i in range(3))*unit.nanometer)
    identities = [(a.name, a.residue.name) for a in topology.topology.atoms()]
    if identities != [(a.name, a.residue.name) for a in native.atoms]:
        raise ValueError('GROMACS export changed peptide atom identity/order')
    systems = []
    for name, reader in [('amber', app.AmberPrmtopFile(str(source/'peptide.prmtop'))),
                         ('gromacs_export', topology)]:
        system = reader.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
        integrator = mm.VerletIntegrator(.001)
        context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
        systems.append((name, system, integrator, context))

    def run(args, cwd, stdin=None):
        env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', GMX_MAXBACKUP='-1')
        proc = subprocess.run([str(executable), *args], cwd=cwd, input=stdin, env=env,
                              capture_output=True, text=True, timeout=90)
        (cwd/(args[0]+'.log')).write_text(proc.stdout+'\n'+proc.stderr)
        if proc.returncode:
            raise RuntimeError(f'{args[0]} failed; output retained at {cwd}')

    run(['--version'], output)
    records = []
    for seed in [None, 2027, 2028]:
        probe = 'source' if seed is None else f'perturbed-{seed}'
        work = output/probe
        work.mkdir()
        xyz = coordinate.copy()
        if seed is not None:
            xyz += np.random.default_rng(seed).normal(0, .0003, xyz.shape)
        native.coordinates = xyz*10
        native.save(str(work/'input.gro'), precision=7)
        run(['grompp', '-f', '../static.mdp', '-p', '../system.top', '-c', 'input.gro', '-o', 'run.tpr'], work)
        run(['mdrun', '-s', 'run.tpr', '-deffnm', 'run', '-ntmpi', '1', '-ntomp', '2', '-nb', 'cpu'], work)
        run(['energy', '-f', 'run.edr', '-o', 'potential.xvg', '-xvg', 'none'], work, 'Potential\n0\n')
        energy = float(np.loadtxt(work/'potential.xvg').reshape(-1, 2)[0, 1])
        with mdtraj.formats.TRRTrajectoryFile(str(work/'run.trr')) as handle:
            arrays = handle._read(1, None, get_forces=True)
        saved = arrays[0][0].astype(float)
        force = arrays[-1][0].astype(float)
        coordinate_error = float(np.max(abs(saved-xyz)))
        if coordinate_error > 2e-6:
            raise ValueError('Native static evaluation changed coordinates')
        np.savez(work/'native-arrays.npz', coordinates_nm=saved, forces_kj_mol_nm=force, energy_kj_mol=energy)
        record = {'probe': probe, 'coordinate_serialization_max_error_nm': coordinate_error,
                  'gromacs_energy_kj_mol': energy, 'comparisons': []}
        reference = []
        for name, _, _, context in systems:
            context.setPositions(saved*unit.nanometer)
            state = context.getState(getEnergy=True, getForces=True)
            e = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            f = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
            if not np.isfinite([e, energy]).all() or not np.isfinite(f).all() or not np.isfinite(force).all():
                raise ValueError('Nonfinite comparison')
            delta = float(np.max(abs(f-force)))
            record['comparisons'].append({'reader': name, 'energy_difference_kj_mol': e-energy,
                'maximum_force_component_difference_kj_mol_nm': delta,
                'force_rms_relative_error': float(np.linalg.norm(f-force)/max(1, np.linalg.norm(f))),
                'passed_native_tolerances': abs(e-energy) <= .001 and delta <= .01})
            reference.append((e, f))
        delta_e = reference[1][0]-reference[0][0]
        delta_f = float(np.max(abs(reference[1][1]-reference[0][1])))
        record['double_precision_export_comparison'] = {'energy_difference_kj_mol': delta_e,
            'maximum_force_component_difference_kj_mol_nm': delta_f,
            'passed': abs(delta_e) <= .0001 and delta_f <= .001}
        records.append(record)
        (output/'progress.json').write_text(json.dumps(records, indent=2)+'\n')
        print(json.dumps(record), flush=True)
    if source_hashes != {name: digest(source/name) for name in source_hashes}:
        raise ValueError('Source changed during comparison')
    report = {'stage': 'native_cross_engine_peptide_comparison_complete', 'physical_model_validated': False,
        'app_ready': False, 'source': str(source), 'atom_count': len(native.atoms),
        'native_tolerances': {'energy_kj_mol': .001, 'force_component_kj_mol_nm': .01},
        'double_precision_export_tolerances': {'energy_kj_mol': .0001, 'force_component_kj_mol_nm': .001},
        'scope': 'Capped peptide only, actual GROMACS mixed precision versus OpenMM Reference at identical saved coordinates. No integration, constraints, modified charges, or omitted force-field terms. Not the solvated complex or a physical accuracy test.',
        'conformations': records, 'source_sha256': source_hashes,
        'all_passed': all(r['double_precision_export_comparison']['passed'] and
                          all(c['passed_native_tolerances'] for c in r['comparisons']) for r in records),
        'output_sha256': {name: digest(output/name) for name in ['system.top', 'static.mdp']},
        'elapsed_seconds': time.monotonic()-start}
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--gromacs', type=Path, required=True)
    args = parser.parse_args()
    try:
        report = compare(args.source, args.output, args.gromacs)
        print(json.dumps({'complete': True, 'all_passed': report['all_passed']}))
    except Exception as exc:
        if args.output.is_dir():
            (args.output/'failure.json').write_text(json.dumps({'type': type(exc).__name__, 'error': str(exc)}, indent=2)+'\n')
        raise
