"""Cross-check live loop reports against independently read saved dihedrals."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from prepare_adduct import digest


def sampling_windows(result, observations):
    released = result.get('warmup_restraints_released_before_heating', False)
    if type(released) is not bool:
        raise ValueError('Warm-up restraint metadata must be explicit boolean state')
    if len(observations) != 100:
        raise ValueError('Expected exactly 100 saved observations')
    for ps, row in enumerate(observations, 1):
        expected_k = 0. if released or ps > 30 else (1000. if ps <= 10 else 100.)
        if float(row['time_ps']) != ps or float(row['position_restraint_k']) != expected_k:
            raise ValueError('Saved restraint schedule disagrees with result metadata')
    return [(1, 10, 'unrestrained heating' if released else 'restrained heating'),
            (11, 30, 'unrestrained NPT equilibration' if released else 'restrained NPT'),
            (31, 100, 'unrestrained NPT')]


def summarize(run, planarity, output):
    result = json.loads((run/'result.json').read_text())
    if (result.get('stage') != 'diagnostic_trajectory_complete'
            or result.get('diagnostic_only') is not True
            or result.get('qualified_for_handoff') is not False):
        raise ValueError('This summary expects a completed, unqualified diagnostic run')
    with (run/'observations.csv').open() as handle:
        windows = sampling_windows(result, list(csv.DictReader(handle)))
    audit = json.loads((planarity/'result.json').read_text())
    if audit['source_sha256'].get(str(run/'trajectory.dcd')) != digest(run/'trajectory.dcd'):
        raise ValueError('Independent peptide audit is not bound to this trajectory')
    for name in ['angles.npz', 'bonds.json']:
        if audit['output_sha256'][name] != digest(planarity/name):
            raise ValueError('Independent angle output changed')
    bonds = json.loads((planarity/'bonds.json').read_text())
    data = np.load(planarity/'angles.npz')
    index = {tuple(row['indices']): i for i, row in enumerate(bonds)}
    angles = data['omega_degrees']
    if angles.shape[0] != result['total_ps'] or result['total_ps'] != 100:
        raise ValueError('Expected 100 frames saved at one ps intervals')
    output.mkdir(exist_ok=False)
    rows = []
    differences = []
    hashes = {str(p): digest(p) for p in [run/'result.json', run/'observations.csv', planarity/'result.json', planarity/'angles.npz', planarity/'bonds.json']}
    definitions = None
    for ps in range(1, 101):
        paths = list((run/'loop-checks').glob(f'{ps:03d}-*.json'))
        if len(paths) != 1:
            raise ValueError('Every saved ps needs exactly one native loop report')
        path = paths[0]
        sample = json.loads(path.read_text())
        if sample['time_ps'] != ps or (not sample['passed'] and not sample['diagnostic_continuation_permitted']):
            raise ValueError('Completed diagnostic contains an inconsistent timestamp or a hard-failure frame')
        hashes[str(path)] = digest(path)
        omega = sample['geometry']['omega_checks']
        keys = [tuple(row['atoms']) for row in omega]
        if definitions is None:
            definitions = keys
        if keys != definitions or any(key not in index for key in keys):
            raise ValueError('Loop peptide identities changed between frames or audits')
        native = np.array([row['degrees'] for row in omega])
        saved = angles[ps-1, [index[key] for key in keys]]
        error = abs((native-saved+180) % 360-180)
        differences.append(float(error.max()))
        if error.max() > .01:
            raise ValueError('Native loop angles disagree with saved DCD beyond coordinate precision')
        rows.append({'time_ps': ps, 'passed': sample['passed'],
            'diagnostic_continuation_permitted': sample['diagnostic_continuation_permitted'],
            'failed_omega_indices': [row['atoms'] for row in omega if not row['accepted']]})
    if [row['time_ps'] for row in rows if not row['passed']] != [row['time_ps'] for row in result['recorded_peptide_excursions']]:
        raise ValueError('Final result does not preserve every recorded failed geometry frame')
    summaries = []
    for key in definitions:
        i = index[key]
        deviation = data['deviation_degrees'][:, i]
        peptide_windows = []
        for lo, hi, label in windows:
            values = deviation[lo-1:hi]
            peptide_windows.append({'phase': label, 'ps': [lo, hi], 'mean_deviation_degrees': float(values.mean()),
                'maximum_deviation_degrees': float(values.max()), 'samples_over_35_degrees': int((values > 35).sum())})
        summaries.append({'indices': list(key), 'left': bonds[i]['left'], 'right': bonds[i]['right'], 'windows': peptide_windows})
    for path, expected in hashes.items():
        if digest(path) != expected:
            raise ValueError('Source changed during diagnostic summary')
    report = {'stage': 'completed_loop_diagnostic_crosscheck', 'diagnostic_only': True,
        'qualified_for_handoff': False, 'app_ready': False, 'physical_model_validated': False,
        'warmup_restraints_released_before_heating': result.get('warmup_restraints_released_before_heating', False),
        'restraint_schedule_verified_from_saved_observations': True,
        'frames': len(rows), 'loop_context_peptide_bonds': len(definitions),
        'failed_geometry_frames': [row['time_ps'] for row in rows if not row['passed']],
        'maximum_native_vs_dcd_angle_error_degrees': max(differences),
        'peptides': summaries, 'source_sha256': hashes,
        'scope': 'Independent DCD replay confirms all monitored peptide angles and preserves every flagged geometry frame. Phase labels are checked against saved positional restraints. Time trends alone do not isolate a restraint effect or establish native loop accuracy.'}
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['peptides', 'source_sha256']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['run', 'planarity', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    summarize(args.run, args.planarity, args.output)
