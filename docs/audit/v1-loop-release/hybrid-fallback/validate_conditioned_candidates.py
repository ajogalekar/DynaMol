"""Fixed paired complete-complex checks of source-conditioned proposals."""
import argparse
import json
from pathlib import Path
import sys
import time

from run_complete_comparison import ROOT, CACHE, PYTHON, sha, save, supervise


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    implementation = output / 'implementation'
    implementation.mkdir()
    for path in [Path(__file__), Path(__file__).with_name('validate_complete_candidate.py'), Path(__file__).with_name('run_complete_comparison.py')]:
        (implementation / path.name).write_bytes(path.read_bytes())
    native = CACHE / 'boundary-conditioned-native-v1'
    native_plan = json.loads((native / 'plan.json').read_text())
    candidates = sorted(native.glob('1UA2/*/native/*/candidate.json'))
    if len(candidates) != 6:
        raise ValueError('The frozen comparison requires every one of the six 1UA2 proposals')
    for p in candidates:
        c = json.loads(p.read_text())
        if c['sidechain_reconstruction_status'] != 'complete':
            raise ValueError('Expected coherent native heavy-atom candidate')
    policies = ['mobile_flanks', 'fixed_observed']
    entries = [{'name': '-'.join(p.relative_to(native).parts[1:4:2]) + '-' + policy,
                'candidate': str(p), 'candidate_sha256': sha(p), 'policy': policy}
               for p in candidates for policy in policies]
    assert len({e['name'] for e in entries}) == len(entries)
    dataset = ROOT / 'docs/audit/loop-fallback/runs/root-final/1UA2/workspace/datasets/4dcb5c25f5c14559'
    plan = {'case': '1UA2', 'dataset': str(dataset), 'attempts': entries,
            'native_plan_sha256': sha(native / 'plan.json'),
            'source_sha256': {str(dataset / p): sha(dataset / p) for p in ['loop-refinement-topology.cif', 'loop-refinement-input.npz', 'loop-refinement-parameters.json']},
            'implementation_sha256': {str(p): sha(p) for p in implementation.iterdir()},
            'selection': 'All six source-conditioned native proposals, each with unchanged mobile-flank refinement and existing all-observed-heavy-fixed refinement. Order frozen before full-complex scores; no filtering or resampling.',
            'limits': {'seconds_per_attempt': 180, 'total_seconds': 1200, 'threads': 2, 'rss_bytes': 4 * 1024**3},
            'invariants': 'Same archived physical parameters and existing 1000-step construction refinement; no stronger forces or relaxed acceptance thresholds. Fixed policy may have strained transplanted seed joins, which remain recorded and never accepted without final checks.',
            'app_ready': False, 'new_md_or_qm': False}
    save(output / 'plan.json', plan)
    started = time.monotonic()
    def checkpoint():
        if time.monotonic() - started > 1200:
            raise TimeoutError('Frozen complete-candidate comparison budget reached')
    results = []
    for entry in entries:
        checkpoint()
        source = Path(entry['candidate'])
        if sha(source) != entry['candidate_sha256']:
            raise ValueError('Frozen native proposal changed')
        folder = output / entry['name']
        folder.mkdir()
        command = [str(PYTHON), '-u', str(Path(__file__).with_name('validate_complete_candidate.py')),
                   '--dataset', str(dataset), '--candidate', str(source), '--output', str(folder / 'validation'),
                   '--observed-context-policy', entry['policy']]
        supervise(command, folder, checkpoint)
        result = json.loads((folder / 'validation/result.json').read_text())
        row = {'attempt': entry['name'], 'policy': entry['policy'], 'result': str(folder / 'validation/result.json'),
               'checks': result['checks'], 'accepted_by_private_checks': result['accepted_by_recorded_checks'],
               'seed_geometry_accepted': result['seed_geometry_accepted'],
               'observed_heavy_motion_A': result['maximum_observed_heavy_displacement_angstrom']}
        results.append(row)
        save(output / 'summary.json', results)
        print(json.dumps(row), flush=True)
    assert all(sha(Path(p)) == h for p, h in plan['source_sha256'].items())
    assert all(sha(Path(p)) == h for p, h in plan['implementation_sha256'].items())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
