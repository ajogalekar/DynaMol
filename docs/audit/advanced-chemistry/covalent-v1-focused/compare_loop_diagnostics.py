"""Compare two completed, unqualified diagnostics without promoting either."""
import argparse
import json
from pathlib import Path

from prepare_adduct import digest


def load_checked(run, summary_path):
    result = json.loads((run/'result.json').read_text())
    summary = json.loads(summary_path.read_text())
    if (result.get('stage') != 'diagnostic_trajectory_complete'
            or summary.get('stage') != 'completed_loop_diagnostic_crosscheck'):
        raise ValueError('Both inputs must be completed diagnostic results')
    for report in [result, summary]:
        if report.get('diagnostic_only') is not True or report.get('qualified_for_handoff') is not False:
            raise ValueError('Diagnostic qualification flags must remain explicit')
    if not summary.get('restraint_schedule_verified_from_saved_observations'):
        raise ValueError('Phase labels need verification against actual saved restraints')
    if summary['source_sha256'].get(str(run/'result.json')) != digest(run/'result.json'):
        raise ValueError('Summary belongs to another run result')
    for name, expected in summary['source_sha256'].items():
        if digest(Path(name)) != expected:
            raise ValueError('Summary source changed: '+name)
    monitor = json.loads((run/'loop-monitor-input.json').read_text())
    launch = json.loads((run/'launch.json').read_text())
    snapshots = run/'implementation-snapshot'
    manifest = json.loads((snapshots/'manifest.json').read_text())
    entry = manifest[launch['runner']]
    if entry['sha256'] != launch['runner_sha256'] or digest(snapshots/entry['snapshot']) != entry['sha256']:
        raise ValueError('Runner snapshot does not match launched implementation')
    return result, summary, monitor


def compare(baseline_run, baseline_summary, comparison_run, comparison_summary, output):
    a, sa, ma = load_checked(baseline_run, baseline_summary)
    b, sb, mb = load_checked(comparison_run, comparison_summary)
    if a['input_sha256'] != b['input_sha256'] or a['case'] != b['case'] or a['seed'] != b['seed']:
        raise ValueError('Comparison must preserve complete model inputs, case and seed')
    for field in ['protein_result_sha256', 'source_to_native_residues', 'remodeled_native_atom_indices',
                  'native_force_classes', 'free_remodeled_region_during_warmup']:
        if ma[field] != mb[field]:
            raise ValueError('Mapped loop preparation or native model changed: '+field)
    if a.get('warmup_restraints_released_before_heating', False) is not False:
        raise ValueError('Baseline must use the original restrained warm-up')
    if b.get('warmup_restraints_released_before_heating') is not True:
        raise ValueError('Comparison must release restraints before heating')
    initial_reports = []
    for name in ['000-before-minimization.json', '000-minimized.json']:
        left = json.loads((baseline_run/'loop-checks'/name).read_text())
        right = json.loads((comparison_run/'loop-checks'/name).read_text())
        if not left['passed'] or not right['passed']:
            raise ValueError('Initial or minimized geometry did not pass')
        pairs = []
        for x, y in zip(left['geometry']['omega_checks'], right['geometry']['omega_checks'], strict=True):
            if x['atoms'] != y['atoms']:
                raise ValueError('Initial or minimized peptide identities differ')
            pairs.append({'indices': x['atoms'], 'baseline_omega_degrees': x['degrees'],
                          'comparison_omega_degrees': y['degrees'],
                          'circular_difference_degrees': (y['degrees']-x['degrees']+180) % 360-180})
        initial_reports.append({'stage': name, 'peptides': pairs})
    peptides = []
    for left, right in zip(sa['peptides'], sb['peptides'], strict=True):
        for key in ['indices', 'left', 'right']:
            if left[key] != right[key]:
                raise ValueError('Peptide identity changed between runs')
        windows = []
        for x, y in zip(left['windows'], right['windows'], strict=True):
            if x['ps'] != y['ps']:
                raise ValueError('Time windows must be identical')
            windows.append({'ps': x['ps'], 'baseline': x, 'comparison': y,
                            'mean_deviation_difference_degrees': y['mean_deviation_degrees']-x['mean_deviation_degrees']})
        peptides.append({key: left[key] for key in ['indices', 'left', 'right']} | {'windows': windows})
    paths = [baseline_summary, comparison_summary]
    for run in [baseline_run, comparison_run]:
        paths += [run/name for name in ['result.json', 'loop-monitor-input.json', 'launch.json',
                  'loop-checks/000-before-minimization.json', 'loop-checks/000-minimized.json']]
    report = {'stage': 'matched_loop_diagnostic_comparison', 'case': a['case'],
              'diagnostic_only': True, 'qualified_for_handoff': False,
              'app_ready': False, 'physical_model_validated': False,
              'input_models_identical_by_sha256': True,
              'source_sha256': {str(p): digest(p) for p in paths},
              'baseline_failed_geometry_frames': sa['failed_geometry_frames'],
              'comparison_failed_geometry_frames': sb['failed_geometry_frames'],
              'initial_and_minimized_peptides': initial_reports, 'peptides': peptides,
              'scope': 'Same input model, preparation, seed and saved time windows; positional-restraint release timing differs. CPU minimization differences are reported. One pair of short trajectories is not a replicated causal test, native loop validation or handoff admission.'}
    output.mkdir(exist_ok=False)
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['source_sha256', 'peptides', 'initial_and_minimized_peptides']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['baseline-run', 'baseline-summary', 'comparison-run', 'comparison-summary', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    compare(args.baseline_run, args.baseline_summary, args.comparison_run, args.comparison_summary, args.output)
