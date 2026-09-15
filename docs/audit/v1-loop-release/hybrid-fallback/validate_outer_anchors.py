"""Frozen complete-complex test of general outer-torsion anchor support."""
import argparse
import json
from pathlib import Path
import time

from run_complete_comparison import ROOT, CACHE, PYTHON, sha, save, supervise


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    implementation = output / 'implementation'
    implementation.mkdir()
    for name in ['validate_outer_anchors.py', 'validate_complete_candidate.py', 'context_integrity.py', 'run_complete_comparison.py']:
        source = Path(__file__).with_name(name)
        (implementation / name).write_bytes(source.read_bytes())
    candidates = [('conditioned', p) for p in sorted((CACHE / 'boundary-conditioned-native-v1').glob('1UA2/*/native/*/candidate.json'))]
    candidates += [('shifted', p) for p in sorted((CACHE / 'shifted-anchor-native-v1').glob('1UA2/native/*/candidate.json'))]
    if len(candidates) != 11:
        raise ValueError('Frozen probe requires all6 conditioned and all5 converged shifted proposals')
    entries = []
    for i, (family, path) in enumerate(candidates, 1):
        data = json.loads(path.read_text())
        if data.get('sidechain_reconstruction_status') != 'complete':
            raise ValueError('Expected complete native heavy-atom proposal')
        entries.append({'name': f'{i:02d}-{family}', 'family': family, 'candidate': str(path), 'candidate_sha256': sha(path)})
    dataset = ROOT / 'docs/audit/loop-fallback/runs/root-final/1UA2/workspace/datasets/4dcb5c25f5c14559'
    plan = {'case': '1UA2', 'dataset': str(dataset), 'attempts': entries, 'policy': 'preserve_outer_torsions',
            'source_sha256': {str(dataset / p): sha(dataset / p) for p in ['loop-refinement-topology.cif', 'loop-refinement-input.npz', 'loop-refinement-parameters.json']},
            'implementation_sha256': {str(p): sha(p) for p in implementation.iterdir()},
            'selection': 'Every returned conditioned1UA2 proposal and every converged shifted-anchor1UA2 proposal. Mechanical native closure failures remain rejected. This full-complex plan is fixed before reading final refinement outcomes.',
            'invariants': 'Existing physical parameters and1000-iteration construction protocol; no stronger forces or relaxed gates. Additional observed fixed atoms are selected by peptide connectivity, never residue number. Every fully observed external torsion definition must stay exact.',
            'reference_scope': 'Modeled/new definitions always checked. External source rows excluded only if every phi/psi/preceding-omega defining coordinate is exactly unchanged; any changed outlier remains a failure.',
            'limits': {'seconds_per_attempt': 180, 'total_seconds': 1200, 'threads': 2, 'rss_bytes': 4 * 1024**3},
            'app_ready': False, 'new_md_or_qm': False}
    save(output / 'plan.json', plan)
    started = time.monotonic()
    def checkpoint():
        if time.monotonic() - started > 1200:
            raise TimeoutError('Frozen outer-anchor comparison budget exhausted')
    results = []
    for entry in entries:
        checkpoint()
        if sha(Path(entry['candidate'])) != entry['candidate_sha256']:
            raise ValueError('Frozen proposal changed')
        folder = output / entry['name']
        folder.mkdir()
        command = [str(PYTHON), '-u', str(Path(__file__).with_name('validate_complete_candidate.py')),
                   '--dataset', str(dataset), '--candidate', entry['candidate'], '--output', str(folder / 'validation'),
                   '--observed-context-policy', 'preserve_outer_torsions']
        supervise(command, folder, checkpoint)
        path = folder / 'validation/result.json'
        result = json.loads(path.read_text())
        reference = json.loads((folder / 'validation/reference.json').read_text())
        row = {'attempt': entry['name'], 'candidate': entry['candidate'], 'result': str(path),
               'checks': result['checks'], 'accepted_by_private_checks': result['accepted_by_recorded_checks'],
               'reference_outliers': [r['residue'] for r in reference['outliers']],
               'observed_heavy_motion_A': result['maximum_observed_heavy_displacement_angstrom']}
        results.append(row)
        save(output / 'summary.json', results)
        print(json.dumps(row), flush=True)
    assert all(sha(Path(p)) == h for p, h in plan['source_sha256'].items())
    for frozen in implementation.iterdir():
        assert sha(frozen) == sha(Path(__file__).with_name(frozen.name))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
