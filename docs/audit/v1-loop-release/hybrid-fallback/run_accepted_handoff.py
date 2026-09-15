"""Exercise real quality rejection and accepted-coordinate handoff privately.

The known outcomes make this an integration regression, not an unbiased
prediction benchmark. No archived artifact or application dataset is changed.
"""
import argparse
import json
from pathlib import Path

from run_complete_comparison import (
    ROOT, CACHE, PYTHON, sha, save, supervise, Attempt, Validation,
    ValidatedArtifact, run_candidate_search,
)


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    dataset = ROOT / 'docs/audit/loop-fallback/runs/root-final/1UA2/workspace/datasets/4dcb5c25f5c14559'
    source = CACHE / 'hybrid-native-mc-v1/1UA2/source.pdb'
    baseline = output / 'archive-fragment.json'
    save(baseline, {'archive_baseline': True, 'app_ready': False})
    proposals = [('archive-fragment', 'archived_fast_fragment', baseline)]
    for mode, generator in [('fragment-source-ccd', 'native_fragment_insertion_mc'),
                            ('torsion-source-ccd', 'native_torsion_mc')]:
        for seed in (2026, 2027):
            path = CACHE / 'boundary-conditioned-native-v1/1UA2' / mode / f'native/43-56-seed-{seed}/candidate.json'
            proposals.append((f'{mode}-{seed}', generator, path))
    files = [source, *[p for _, _, p in proposals], *[
        dataset / name for name in ('loop-refinement-topology.cif',
                                    'loop-refinement-input.npz',
                                    'loop-refinement-parameters.json')]]
    frozen = {str(path): sha(path) for path in files}
    implementation = output / 'implementation'
    implementation.mkdir()
    for path in [Path(__file__), Path(__file__).with_name('validate_complete_candidate.py'),
                 Path(__file__).with_name('context_integrity.py'),
                 Path(__file__).with_name('run_complete_comparison.py'),
                 ROOT / 'backend/loop_search.py']:
        (implementation / path.name).write_bytes(path.read_bytes())
    plan = {
        'purpose': 'Known-outcome integration regression: native full refinement, quality-rejection failover, final coordinate-artifact handoff.',
        'source_sha256': frozen[str(source)], 'frozen_inputs': frozen,
        'attempts': [{'name': n, 'generator': g, 'proposal': str(p)} for n, g, p in proposals],
        'policy': 'preserve_outer_torsions',
        'limits': {'max_attempts': 5, 'deadline_seconds': 900, 'seconds_per_child': 180,
                   'rss_bytes_per_child': 4 * 1024**3, 'threads': 2},
        'implementation_sha256': {str(p): sha(p) for p in implementation.iterdir()},
        'app_ready': False, 'full_preparation_rerun': False, 'new_md_or_qm': False,
    }
    save(output / 'plan.json', plan)
    by_name = {name: path for name, _, path in proposals}
    events = []

    def generate(attempt, private, checkpoint):
        checkpoint()
        original = by_name[attempt.name]
        if sha(original) != frozen[str(original)]:
            raise ValueError('A frozen proposal changed before handoff')
        copied = private / 'proposal.json'
        copied.write_bytes(original.read_bytes())
        return copied

    def validate(candidate, private, checkpoint):
        command = [str(PYTHON), '-u', str(Path(__file__).with_name('validate_complete_candidate.py')),
                   '--dataset', str(dataset), '--candidate', str(candidate),
                   '--output', str(private / 'validation'),
                   '--observed-context-policy', 'preserve_outer_torsions']
        supervise(command, private, checkpoint)
        report_path = private / 'validation/result.json'
        report = json.loads(report_path.read_text())
        artifact = None
        if all(value is True for value in report['checks'].values()):
            coordinates = private / 'validation/refined.npz'
            recorded = report['artifacts_sha256'][str(coordinates)]
            artifact = ValidatedArtifact(coordinates, recorded)
        return Validation(report['checks'], tuple(report['reasons']), artifact=artifact)

    def progress(event):
        events.append(event)
        save(output / 'progress.json', events)
        print(json.dumps(event), flush=True)

    accepted = run_candidate_search(
        [Attempt(name, generator) for name, generator, _ in proposals],
        output / 'search', generate, validate, source_sha256=frozen[str(source)],
        max_attempts=5, deadline_seconds=900, on_progress=progress,
    )
    search = json.loads((output / 'search/search.json').read_text())
    if accepted.path.name != 'refined.npz' or sha(accepted.path) != accepted.sha256:
        raise ValueError('Search did not return the final validated coordinates')
    if any(sha(Path(path)) != digest for path, digest in frozen.items()):
        raise ValueError('A frozen input changed during the integration test')
    for path, digest in plan['implementation_sha256'].items():
        if sha(Path(path)) != digest:
            raise ValueError('A saved implementation snapshot changed')
    summary = {
        'status': search['status'], 'accepted_candidate': search['accepted_candidate'],
        'accepted_artifact': search['accepted_artifact'],
        'earlier_rejected_candidates': [r['name'] for r in search['attempts'] if r['status'] == 'rejected'],
        'attempts_executed': len(search['attempts']), 'progress_events': len(events),
        'all_frozen_inputs_unchanged': True,
        'final_coordinates_returned_instead_of_proposal_or_report': True,
        'app_ready': False, 'full_preparation_rerun': False, 'new_md_or_qm': False,
    }
    save(output / 'summary.json', summary)
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    run(args.output.resolve())
