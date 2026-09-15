"""Retain the scope and actual diagnostics of the two warmed engine specimens."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import mdtraj

p = argparse.ArgumentParser()
p.add_argument('--source', type=Path, required=True)
p.add_argument('--gromacs', type=Path, required=True)
p.add_argument('--openmm', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
source, gm, om = a.source.resolve(), a.gromacs.resolve(), a.openmm.resolve()
review = json.loads((source / 'comparison-readiness.json').read_text())
g = json.loads((gm / 'stage-status.json').read_text())
gg = json.loads((gm / 'independent-site-geometry.json').read_text())
o = json.loads((om / 'comparison-report.json').read_text())
if review.get('ready_for_short_comparison') is not True or g['status'] != 'completed' or o['status'] != 'completed':
    raise ValueError('A reviewed endpoint and two completed runs are required')
end = np.load(source / 'endpoint.npz')
images = np.load(om / 'initial-image-mapping.npz')
with mdtraj.formats.TRRTrajectoryFile(str(gm / 'stage.trr')) as stream:
    gx, gt, gs, gb, *other = stream._read(1, None, get_forces=True)
box = end['box_nm'].astype(float)
delta = gx[0].astype(float) - end['coords_nm'].astype(float)
delta -= np.rint(delta @ np.linalg.inv(box)) @ box
gm_max_start = float(np.linalg.norm(delta, axis=1).max())
gm_box_error = float(np.max(np.abs(np.array(g['box_samples_nm']) - box)))
om_raw_exact = bool(np.array_equal(images['raw_coords_nm'], end['coords_nm']))
if gm_max_start > 2e-6 or gm_box_error > 2e-6 or not om_raw_exact or not np.array_equal(images['box_nm'], box):
    raise ValueError('The two specimens did not use the same physical coordinates/box')
if o['mobile_dof'] != g['mobile_dof']:
    raise ValueError('Mobile degrees of freedom differ')
runtime_thermo = np.loadtxt(om / 'thermodynamics.csv', delimiter=',', skiprows=1)
volume_error = float(np.max(np.abs(runtime_thermo[:, 5] - o['volume_nm3'])))
recomputed_temperature = 2 * runtime_thermo[:, 3] / (o['mobile_dof'] * .00831446261815324)
temperature_error = float(np.max(np.abs(recomputed_temperature - runtime_thermo[:, 4])))
if volume_error > 1e-6 or temperature_error > 1e-6:
    raise ValueError('Actual OpenMM runtime box/temperature differs from declared diagnostics')
periodic = source.parent / 'prototype-1ca2-periodic-diagnostic/periodic-comparison-matrix.json'
periodic_report = json.loads(periodic.read_text())
evidence = [(source, 'comparison-readiness.json'), (source, 'endpoint.npz'),
            (gm, 'stage-status.json'), (gm, 'independent-site-geometry.json'),
            (om, 'comparison-report.json'), (om, 'system.xml'), (om, 'initial-image-mapping.npz'),
            (om, 'thermodynamics.csv'), (periodic.parent, periodic.name)]
report = {
    'scope': 'Short unrestrained warmed 1CA2 neutral-water ZAFF6 engineering comparison; one starting conformation and one independent specimen per engine.',
    'limitations': ['Not a thermodynamic-convergence assessment, catalytic model validation, experimental accuracy claim, or general metal-family acceptance.',
                    'GROMACS v-rescale and OpenMM LangevinMiddle use different thermostats and velocities; trajectories are not expected to coincide.',
                    'Geometry limits are broad alarms. Bond survival in a bonded model is necessary implementation evidence, not independent proof of chemical accuracy.'],
    'source_review': review,
    'shared_initial_state': {'gromacs_maximum_minimum_image_coordinate_difference_nm': gm_max_start,
                             'gromacs_maximum_box_difference_nm': gm_box_error,
                             'openmm_raw_coordinates_and_box_exact': om_raw_exact,
                             'openmm_initial_constraint_projection_nm': o['maximum_initial_constraint_projection_nm'],
                             'actual_openmm_runtime_volume_maximum_error_nm3': volume_error,
                             'actual_openmm_kinetic_temperature_maximum_recalculation_error_kelvin': temperature_error,
                             'mobile_dof': o['mobile_dof'], 'volume_nm3': o['volume_nm3'], 'density_g_ml': o['density_g_ml']},
    'periodic_energy_comparison': {'status': periodic_report['status'], 'evidence': str(periodic),
                                    'meaning': 'Finite motion and coordination do not resolve the separately retained absolute periodic-energy discrepancy.'},
    'finite_unrestrained_motion_and_broad_coordination_checks_passed': all([
        g['finite_all_saved_forces'], gg['broad_geometry_alarms_clear'],
        o['finite_every_saved_force'], o['broad_geometry_alarms_clear'], not g['restraints'], not o['restraints']]),
    'gromacs': {**{k: g[k] for k in ['stage_duration_ps', 'last10ps_metrics', 'elapsed_seconds', 'frame_count']},
                **{k: gg[k] for k in ['Zn_min_angstrom', 'Zn_max_angstrom', 'Zn_mean_angstrom', 'Zn_angle_mean_degrees', 'maximum_CA_rmsd_angstrom']}},
    'openmm': {k: o[k] for k in ['duration_ps', 'last10ps_metrics', 'elapsed_seconds', 'Zn_min_angstrom',
                               'Zn_max_angstrom', 'Zn_mean_angstrom', 'Zn_angle_mean_degrees', 'maximum_CA_rmsd_angstrom']},
    'evidence': [{'path': str(folder / name), 'sha256': hashlib.sha256((folder / name).read_bytes()).hexdigest()}
                 for folder, name in evidence],
}
a.output.resolve().write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: report[k] for k in ['scope', 'finite_unrestrained_motion_and_broad_coordination_checks_passed', 'shared_initial_state']}, indent=2))
