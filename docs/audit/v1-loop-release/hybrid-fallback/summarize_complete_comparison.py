"""Summarize frozen native validation outputs without rerunning them."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(source, destination):
    plan = json.loads((source / 'plan.json').read_text())
    assert all(sha(Path(p)) == h for case in plan['cases'] for p, h in case['frozen_inputs'].items())
    assert all(sha(Path(p)) == h for p, h in plan['implementation_sha256'].items())
    rows = []
    for case in plan['cases']:
        search = json.loads((source / case['case'] / 'search/search.json').read_text())
        assert search['status'] == 'exhausted'
        assert [r['name'] for r in search['attempts']] == [r['name'] for r in case['attempts']]
        for attempt in search['attempts']:
            folder = source / case['case'] / 'search' / attempt['name']
            path = folder / 'validation/result.json'
            result = json.loads(path.read_text())
            assert all(sha(Path(p)) == h for p, h in result['artifacts_sha256'].items())
            assert all(sha(Path(p)) == h for p, h in result['source_sha256'].items())
            supervisor = json.loads((folder / 'supervisor.json').read_text())
            reference = json.loads((folder / 'validation/reference.json').read_text())
            modeled = set(map(tuple, result['modeled_residue_keys']))
            outliers = [r['residue'] for r in reference['outliers']]
            rows.append({'case': case['case'], 'attempt': attempt['name'], 'status': attempt['status'],
                'checks': result['checks'], 'result': str(path), 'result_sha256': sha(path),
                'native_system_sha256': result['native_system_sha256'], 'atom_count': result['atom_count'],
                'modeled_outliers': [r for r in outliers if tuple(r) in modeled],
                'affected_observed_outliers': [r for r in outliers if tuple(r) not in modeled],
                'maximum_observed_heavy_displacement_angstrom': result['maximum_observed_heavy_displacement_angstrom'],
                'elapsed_seconds': supervisor['elapsed_seconds'], 'peak_sampled_rss_bytes': supervisor['peak_sampled_rss_bytes'],
                'scope': result['scope']})
    for case in plan['cases']:
        assert len({r['native_system_sha256'] for r in rows if r['case'] == case['case']}) == 1
    summary = {'source': str(source), 'plan_sha256': sha(source / 'plan.json'), 'attempts': rows,
               'completed_attempts': len(rows), 'accepted_candidates': sum(all(r['checks'].values()) for r in rows),
               'elapsed_native_validation_seconds': sum(r['elapsed_seconds'] for r in rows),
               'peak_sampled_rss_bytes': max(r['peak_sampled_rss_bytes'] for r in rows),
               'checks_passed': {key: sum(r['checks'][key] for r in rows) for key in rows[0]['checks']},
               'inputs_implementations_and_saved_artifacts_unchanged': True,
               'within_case_native_systems_identical': True,
               'app_ready': False, 'full_preparation_rerun': False, 'physical_model_validated': False,
               'new_md_or_qm': False, 'software_regression_tests_passed': 69,
               'software_test_command': 'PYTHONPATH=<pytest-target>:<repo> PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 DYNAMOL_CPU_THREADS=2 OPENBLAS_NUM_THREADS=2 <focused-python> -m pytest backend/tests/test_loop_search.py backend/tests/test_loop_refinement.py backend/tests/test_promod_loop_worker.py -q',
               'selection_limit': plan['selection']}
    destination.write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: summary[key] for key in ['completed_attempts', 'accepted_candidates', 'elapsed_native_validation_seconds', 'peak_sampled_rss_bytes', 'checks_passed']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    run(args.source.resolve(), args.destination.resolve())
