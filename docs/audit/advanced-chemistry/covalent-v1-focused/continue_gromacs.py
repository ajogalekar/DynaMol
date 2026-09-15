"""Bounded native GROMACS continuation of an audited OpenMM research endpoint.

No chemistry is fitted here. This exercises the full-complex engine handoff,
not an independent preparation or a validation of the physical force field.
"""
import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

import mdtraj as md
import numpy as np
import openmm as mm
from openmm import app, unit
import parmed
import psutil

from prepare_adduct import digest
from precise_gromacs_export import write_precise_topology


def write_json(path, obj):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(obj, indent=2)+'\n')
    temporary.replace(path)


def run_engine(executable, arguments, output, timeout, stdin=None):
    log = output/(arguments[0]+'-launcher.log')
    started = time.monotonic()
    env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1')
    with log.open('w') as handle:
        process = subprocess.Popen([str(executable), *arguments], cwd=output, env=env,
            stdin=subprocess.PIPE if stdin else subprocess.DEVNULL,
            stdout=handle, stderr=subprocess.STDOUT, text=True, start_new_session=True)
        if stdin:
            process.stdin.write(stdin)
            process.stdin.close()
        peak = 0
        try:
            while process.poll() is None:
                rss = 0
                try:
                    parent = psutil.Process(process.pid)
                    for child in [parent, *parent.children(recursive=True)]:
                        try:
                            rss += child.memory_info().rss
                        except psutil.NoSuchProcess:
                            pass
                except psutil.NoSuchProcess:
                    pass
                peak = max(peak, rss)
                elapsed = time.monotonic()-started
                progress = {'command': arguments[0], 'pid': process.pid,
                    'elapsed_seconds': elapsed, 'rss_bytes': rss, 'peak_sampled_rss_bytes': peak}
                if arguments[0] == 'mdrun' and (output/'stage.log').exists():
                    text = (output/'stage.log').read_text(errors='replace')[-5000:]
                    steps = re.findall(r'Step\s+Time\s*\n\s*(\d+)\s+([\d.]+)', text)
                    if steps:
                        progress.update(step=int(steps[-1][0]), time_ps=float(steps[-1][1]))
                write_json(output/'progress.json', progress)
                if elapsed > timeout or rss > 4_000_000_000:
                    raise TimeoutError('Native engine exceeded its recorded time or memory bound')
                time.sleep(1)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise
        if process.returncode:
            raise RuntimeError(f'{arguments[0]} failed; inspect {log}')
    return {'command': arguments, 'returncode': process.returncode,
        'elapsed_seconds': time.monotonic()-started, 'peak_sampled_rss_bytes': peak}


def run(inputs, endpoint, output, executable, duration):
    if duration not in (10, 20, 50, 100):
        raise ValueError('Research continuation is bounded to 10, 20, 50 or 100 ps')
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    result = json.loads((endpoint/'result.json').read_text())
    if result.get('diagnostic_only') or result.get('qualified_for_handoff') is False:
        raise ValueError('A diagnostic-only trajectory is not a qualified engine-handoff endpoint')
    audit = json.loads((endpoint/'saved-run-audit.json').read_text())
    if result['total_ps'] != 100 or audit['frames'] != 100 or result['final']['position_restraint_k'] != 0:
        raise ValueError('Expected an audited, unrestrained endpoint')
    for name, expected in result['input_sha256'].items():
        if digest(inputs/name) != expected:
            raise ValueError('Source input changed: '+name)
    for name in ['final-state.xml', 'system.xml']:
        if digest(endpoint/name) != audit['output_sha256'][name]:
            raise ValueError('Audited endpoint changed: '+name)
    source_hashes = {str(path): digest(path) for path in
        [inputs/'solvated.prmtop', inputs/'source-mapping.json', endpoint/'final-state.xml']}
    write_json(output/'launch.json', {'pid': os.getpid(), 'started_unix': time.time(),
        'runner': str(Path(__file__).resolve()), 'runner_sha256': digest(__file__),
        'case': result['case'], 'inputs': str(inputs), 'endpoint': str(endpoint),
        'duration_ps': duration, 'mdrun_time_limit_seconds': 3600,
        'native_threads': 2, 'rss_limit_bytes': 4_000_000_000})
    native = parmed.load_file(str(inputs/'solvated.prmtop'))
    state = mm.XmlSerializer.deserialize((endpoint/'final-state.xml').read_text())
    coordinates = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    velocities = np.asarray(state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond))
    box = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    if np.max(abs(box-np.diag(np.diag(box)))) > 1e-10:
        raise ValueError('This research handoff currently requires a rectangular box')
    native.coordinates = coordinates*10
    native.velocities = velocities*10
    native.box = [*(np.diag(box)*10), 90, 90, 90]
    precision = write_precise_topology(native, output/'system.top')
    native.save(str(output/'input.gro'), precision=7)
    reread = parmed.load_file(str(output/'system.top'), xyz=str(output/'input.gro'))
    identity = [(a.name, a.residue.name, a.atomic_number) for a in native.atoms]
    if identity != [(a.name, a.residue.name, a.atomic_number) for a in reread.atoms]:
        raise ValueError('Export changed source atom ordering or identity')
    charges = np.asarray([a.charge for a in native.atoms])
    delta_charge = float(np.max(abs(charges-np.asarray([a.charge for a in reread.atoms]))))
    if delta_charge > 6e-9:
        raise ValueError('Export changed charges beyond its decimal precision')
    # Loading a topology with xyz= reads GRO coordinates but does not attach
    # velocities. Read the coordinate container independently for this check.
    saved_gro = parmed.load_file(str(output/'input.gro'))
    coordinate_error = float(np.max(abs(saved_gro.coordinates*.1-coordinates)))
    velocity_error = float(np.max(abs(saved_gro.velocities*.1-velocities)))
    if max(coordinate_error, velocity_error) > 6e-8:
        raise ValueError('Coordinate/velocity export differs beyond decimal precision')
    if {tuple(sorted((b.atom1.idx, b.atom2.idx))) for b in native.bonds} != {
        tuple(sorted((b.atom1.idx, b.atom2.idx))) for b in reread.bonds}:
        raise ValueError('Export changed the complete source bond graph')
    # Same Reference platform, cutoff and PME tolerance isolate topology transport.
    values = []
    for name, reader in [('amber', app.AmberPrmtopFile(str(inputs/'solvated.prmtop'))),
                         ('gromacs', app.GromacsTopFile(str(output/'system.top'),
                          periodicBoxVectors=tuple(mm.Vec3(*row) for row in box)*unit.nanometer))]:
        system = reader.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*unit.nanometer,
            constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=1e-5)
        system.setDefaultPeriodicBoxVectors(*(mm.Vec3(*row) for row in box))
        integrator = mm.VerletIntegrator(.002)
        context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
        context.setPositions(coordinates*unit.nanometer)
        sample = context.getState(getEnergy=True, getForces=True)
        energy = float(sample.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        force = np.asarray(sample.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
        if not np.isfinite(energy) or not np.isfinite(force).all():
            raise ValueError('Nonfinite full-complex parameter comparison')
        values.append((energy, force, system.getNumConstraints()))
        del context, integrator
    energy_error = values[1][0]-values[0][0]
    force_error = float(np.max(abs(values[1][1]-values[0][1])))
    handoff = {'atom_count': len(native.atoms), 'source_atom_order_and_bonds_preserved': True,
        'precision_export': precision,
        'charge_max_error_e': delta_charge, 'coordinate_max_error_nm': coordinate_error,
        'velocity_max_error_nm_ps': velocity_error, 'source_box_nm': box.tolist(),
        'source_energy_kj_mol': values[0][0], 'export_energy_kj_mol': values[1][0],
        'energy_difference_kj_mol': energy_error, 'max_force_difference_kj_mol_nm': force_error,
        'constraint_counts': [x[2] for x in values],
        'tolerances': {'energy_kj_mol': .001, 'force_component_kj_mol_nm': .01},
        'scope': 'Full solvated topology transport in OpenMM Reference PME, not native GROMACS arithmetic or physical model accuracy.'}
    write_json(output/'handoff.json', handoff)
    if abs(energy_error) > .001 or force_error > .01 or values[0][2] != values[1][2]:
        raise ValueError('Full-complex parameter transport failed its fixed tolerances')
    mdp = {'integrator': 'md', 'dt': '.002', 'nsteps': str(duration*500), 'tinit': '0',
        'nstxout': '500', 'nstvout': '500', 'nstfout': '500', 'nstenergy': '500', 'nstlog': '500',
        'cutoff-scheme': 'Verlet', 'nstlist': '20', 'rlist': '1.0',
        'coulombtype': 'PME', 'rcoulomb': '1.0', 'coulomb-modifier': 'None', 'ewald-rtol': '1e-5',
        'pme-order': '4', 'fourierspacing': '.10', 'vdwtype': 'Cut-off', 'rvdw': '1.0',
        'vdw-modifier': 'None', 'DispCorr': 'EnerPres', 'pbc': 'xyz',
        'constraints': 'h-bonds', 'constraint-algorithm': 'lincs', 'lincs-order': '6', 'lincs-iter': '2',
        'continuation': 'yes', 'gen-vel': 'no', 'tcoupl': 'v-rescale', 'tc-grps': 'System',
        'tau-t': '1.0', 'ref-t': '300', 'nsttcouple': '10', 'ld-seed': '3030',
        'pcoupl': 'C-rescale', 'pcoupltype': 'isotropic', 'tau-p': '2.0', 'ref-p': '1.0',
        'compressibility': '4.5e-5', 'nstpcouple': '10', 'comm-mode': 'Linear', 'nstcomm': '100'}
    (output/'stage.mdp').write_text('\n'.join(f'{k} = {v}' for k, v in mdp.items())+'\n')
    commands = [run_engine(executable, ['--version'], output, 30)]
    commands.append(run_engine(executable, ['grompp', '-f', 'stage.mdp', '-c', 'input.gro',
        '-p', 'system.top', '-o', 'stage.tpr'], output, 90))
    write_json(output/'prepared.json', {'case': result['case'], 'ready_for_native_continuation': True,
        'duration_ps': duration, 'restraints': False, 'atom_count': len(native.atoms), 'commands': commands})
    print('Full-complex transport passed; starting native unrestrained GROMACS continuation', flush=True)
    commands.append(run_engine(executable, ['mdrun', '-s', 'stage.tpr', '-deffnm', 'stage',
        '-ntmpi', '1', '-ntomp', '2', '-nb', 'cpu', '-pin', 'off'], output, 3600))
    if 'LINCS WARNING' in (output/'stage.log').read_text():
        raise ValueError('Constraint warning during native continuation')
    commands.append(run_engine(executable, ['energy', '-f', 'stage.edr', '-o', 'thermodynamics.xvg',
        '-xvg', 'none'], output, 30, 'Potential\nKinetic-En.\nTemperature\nPressure\nVolume\nDensity\n0\n'))
    thermo = np.loadtxt(output/'thermodynamics.xvg')
    with md.formats.TRRTrajectoryFile(str(output/'stage.trr')) as handle:
        frames = handle._read(duration+2, None, get_forces=True, get_velocities=True)
    xyz, times, steps, boxes = frames[:4]
    forces = frames[-1]
    if xyz.shape != (duration+1, len(native.atoms), 3) or abs(times[-1]-duration) > .001:
        raise ValueError('Incomplete saved native trajectory')
    if not all(np.isfinite(x).all() for x in [thermo, xyz, boxes, forces]):
        raise ValueError('Nonfinite saved native trajectory/forces/thermodynamics')
    trajectory = md.Trajectory(xyz, md.Topology.from_openmm(app.AmberPrmtopFile(str(inputs/'solvated.prmtop')).topology), time=times)
    trajectory.unitcell_vectors = boxes
    cov = next(r for r in native.residues if r.name == 'COV')
    bonds = [b for b in native.bonds if b.atom1.residue is cov and b.atom2.residue is cov]
    distances = md.compute_distances(trajectory, [(b.atom1.idx, b.atom2.idx) for b in bonds], periodic=True)*10
    ratios = distances/np.asarray([b.type.req for b in bonds])
    if ratios.min() < .6 or ratios.max() > 1.5:
        raise ValueError('Distorted covalent-adduct bond in native trajectory')
    attachment = [i for i, b in enumerate(bonds) if {b.atom1.atomic_number, b.atom2.atomic_number} == {6, 16}
        and {b.atom1.name, b.atom2.name} != {'CB', 'SG'}]
    if len(attachment) != 1:
        raise ValueError('Expected one recorded cysteine–drug attachment')
    for center, neighbors in result['stereocenter_indices']:
        edges = xyz[:, neighbors, :]-xyz[:, center, None, :]
        volume = np.linalg.det(edges)
        original_volume = np.linalg.det(coordinates[neighbors]-coordinates[center])
        if np.any(volume*original_volume <= 0):
            raise ValueError('Mapped stereocenter inverted')
    metal_sites = []
    for site in result.get('initial_metal_sites', []):
        donor_distances = md.compute_distances(trajectory,
            [(site['metal_index'], donor['index']) for donor in site['donors']], periodic=True)*10
        metal_sites.append({'metal_index': site['metal_index'], 'donors': [dict(donor,
            sampled_min_angstrom=float(donor_distances[:, i].min()),
            sampled_max_angstrom=float(donor_distances[:, i].max()),
            sampled_mean_angstrom=float(donor_distances[:, i].mean()))
            for i, donor in enumerate(site['donors'])]})
    for path, expected in source_hashes.items():
        if digest(path) != expected:
            raise ValueError('Source changed during native continuation')
    report = {'stage': 'native_gromacs_continuation_complete', 'case': result['case'], 'duration_ps': duration,
        'app_ready': False, 'physical_model_validated': False, 'independent_preparation': False,
        'initial_endpoint': str(endpoint), 'source_hashes': source_hashes, 'atom_count': len(native.atoms),
        'saved_frames': len(xyz), 'unrestrained': True, 'stereocenters_preserved': True,
        'metal_sites': metal_sites,
        'finite_saved_forces': True, 'no_lincs_warnings': True,
        'mean_temperature_k': float(np.mean(thermo[1:, 3])), 'final_density_g_ml': float(thermo[-1, 6]/1000),
        'attachment_range_angstrom': [float(distances[:, attachment[0]].min()), float(distances[:, attachment[0]].max())],
        'adduct_bond_ratio_range': [float(ratios.min()), float(ratios.max())], 'commands': commands,
        'scope': 'Native GROMACS 20/50/100 ps numerical continuation of the same research force field after an audited OpenMM endpoint. Not a bitwise restart, independent preparation, or physical-force-field validation.',
        'elapsed_seconds': time.monotonic()-started,
        'output_sha256': {name: digest(output/name) for name in ['system.top', 'stage.mdp', 'stage.tpr', 'stage.trr', 'stage.gro', 'stage.cpt']}}
    write_json(output/'result.json', report)
    print(json.dumps({k: report[k] for k in ['stage', 'case', 'duration_ps', 'mean_temperature_k', 'attachment_range_angstrom']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['inputs', 'endpoint', 'output', 'gromacs']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--duration-ps', type=int, default=20)
    args = parser.parse_args()
    try:
        run(args.inputs, args.endpoint, args.output, args.gromacs, args.duration_ps)
    except Exception as exc:
        if args.output.is_dir():
            write_json(args.output/'failure.json', {'type': type(exc).__name__, 'error': str(exc)})
        raise
