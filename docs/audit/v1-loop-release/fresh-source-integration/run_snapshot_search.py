"""Bounded private full-snapshot loop search; no app dataset publication."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

from backend.loop_process import supervise_loop_child
from backend.loop_search import (Attempt, CandidateUnavailable, SearchExhausted,
                                 Validation, ValidatedArtifact, run_candidate_search)
from backend.loop_snapshot import load_snapshot, verify_accepted_bundle

REPO = Path('/Users/ashujo/Documents/Science/DynaMol')
PYTHON = Path('/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python')
PROMOD = Path('/Users/ashujo/Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2')
REFERENCE = Path('/Users/ashujo/.cache/dynamol-runtimes/cctbx-rama-reference-v1')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    with Path(path).open('x') as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + '\n')


def run(snapshot_folder, snapshot_hash, output, seed=42, context_only=False, worker_source=None):
    snapshot = load_snapshot(snapshot_folder, expected_sha256=snapshot_hash)
    output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=False)
    sources = {f.name: f for f in snapshot.source_files}
    source = sources['loop-model-context.pdb']
    native_input = json.loads(sources['loop-model-input.json'].path.read_text())
    implementation = output / 'implementation'; implementation.mkdir()
    worker = implementation / 'promod_candidate_worker.py'
    for name in ('promod_candidate_worker.py', 'promod_loop_worker.py'):
        source_file = Path(worker_source) if worker_source and name == 'promod_candidate_worker.py' else REPO / 'backend' / name
        shutil.copyfile(source_file, implementation / name)
    validator = implementation / 'validate_snapshot.py'
    shutil.copyfile(Path(__file__).with_name('validate_snapshot.py'), validator)
    shutil.copyfile(__file__, implementation / 'run_snapshot_search.py')
    code = [*implementation.glob('*.py'), *sorted((REPO / 'backend').glob('*.py'))]
    pins = {str(p): sha(p) for p in code}
    environment = dict(os.environ)
    environment.update(OPENMM_CPU_THREADS='2', DYNAMOL_CPU_THREADS='2', OPENBLAS_NUM_THREADS='2',
                       OMP_NUM_THREADS='2', MKL_NUM_THREADS='2', VECLIB_MAXIMUM_THREADS='2',
                       PYTHONDONTWRITEBYTECODE='1', DYNAMOL_PROMOD3=str(PROMOD),
                       DYNAMOL_DATA_DIR=str(output / 'workspace'),
                       PYTHONPATH=os.pathsep.join([str(REPO), str(REFERENCE)]))
    # One fast fragment attempt, then bounded independently seeded fallbacks.
    plans = [('fragment', seed, 0), ('torsion_mc', seed, 0), ('torsion_mc', 2026, 0),
             ('torsion_mc', 2027, 0), ('fragment_mc', seed, 0), ('fragment_mc', 2027, 0),
             ('fragment', seed, 1), ('fragment', seed, 2)]
    if context_only:
        plans = [('fragment', seed, 1), ('fragment', seed, 2),
                 ('fragment', 2027, 1), ('fragment', 2027, 2)]
    plans = list(dict.fromkeys(plans))
    by_name = {f'{strategy.replace("_", "-")}{"-context-" + str(extension) if extension else ""}-{s}':
               (strategy, s, extension) for strategy, s, extension in plans}
    attempts = [Attempt(name, value[0]) for name, value in by_name.items()]
    write(output / 'plan.json', {'snapshot': str(snapshot.manifest_path), 'snapshot_sha256': snapshot_hash,
          'source_sha256': source.sha256, 'attempts': by_name, 'native_generation_timeout_seconds': 180,
          'validation_timeout_seconds': 300, 'search_deadline_seconds': 900, 'rss_bytes': 4 * 1024**3,
          'threads': 2, 'source_code_sha256': pins, 'app_publication': False, 'MD_or_parameter_fitting': False})

    def verify():
        for path, expected in pins.items():
            if sha(path) != expected:
                raise ValueError('Frozen search implementation changed: ' + path)

    def generate(attempt, private, checkpoint):
        verify(); checkpoint()
        strategy, current_seed, extension = by_name[attempt.name]
        request = {**native_input, 'input_pdb': str(source.path), 'source_sha256': source.sha256,
                   'strategy': strategy, 'seed': current_seed, 'max_res_extension': extension}
        write(private / 'native-input.json', request)
        proposal = private / 'proposal.json'
        code = supervise_loop_child([str(PROMOD / 'bin/python'), '-I', str(worker),
                                    str(private / 'native-input.json'), str(proposal)],
                                    private / 'generation-process', environment=environment,
                                    check_cancel=checkpoint, timeout_seconds=180)
        verify(); checkpoint()
        if not proposal.is_file():
            raise RuntimeError('Native generator produced no result record')
        result = json.loads(proposal.read_text())
        if code == 2 and result.get('status') == 'unavailable':
            raise CandidateUnavailable(result.get('error', result.get('reason', 'No native candidate was available')))
        if code != 0 or result.get('status') != 'candidate':
            raise RuntimeError('Native generation failed: ' + json.dumps(result))
        return {'path': str(proposal), 'sha256': sha(proposal)}

    def validate(candidate, private, checkpoint):
        verify(); checkpoint()
        request = {'snapshot': str(snapshot.manifest_path.parent), 'snapshot_sha256': snapshot_hash,
                   'proposal': candidate['path'], 'proposal_sha256': candidate['sha256'],
                   'output': str(private / 'complete')}
        write(private / 'validation-input.json', request)
        result_path = private / 'validation-result.json'
        code = supervise_loop_child([str(PYTHON), '-B', str(validator),
                                    str(private / 'validation-input.json'), str(result_path)],
                                    private / 'validation-process', environment=environment,
                                    check_cancel=checkpoint, timeout_seconds=300)
        verify(); checkpoint()
        if not result_path.is_file():
            raise RuntimeError('Complete validator produced no decision record')
        result = json.loads(result_path.read_text())
        if code not in (0, 2) or result.get('status') not in ('static_checks_passed', 'rejected'):
            raise RuntimeError('Complete validator failed: ' + json.dumps(result))
        checks = dict(result['all_required_checks'])
        checks['observed_displacement'] = checks['observed_displacement'] and result['original_source_flank_cap_passed']
        if result['status'] == 'rejected':
            if all(checks.values()):
                raise ValueError('Rejected decision failed to identify its incomplete check')
            return Validation(checks, tuple(name for name, passed in checks.items() if not passed))
        if code != 0 or not all(value is True for value in checks.values()):
            raise ValueError('Accepted decision does not pass all complete checks')
        bundle = verify_accepted_bundle(result['accepted_bundle'], expected_sha256=result['accepted_bundle_sha256'])
        if bundle.snapshot.sha256 != snapshot_hash:
            raise ValueError('Accepted bundle uses a different source preparation')
        coordinates = bundle.artifacts['coordinates']
        if len(coordinates) != 1 or coordinates[0].name != 'refined.npz':
            raise ValueError('The final coordinate handoff is ambiguous')
        return Validation(checks, artifact=ValidatedArtifact(coordinates[0].path, coordinates[0].sha256))

    def progress(row):
        print(json.dumps(row), flush=True)
        with (output / 'progress.jsonl').open('a') as handle:
            handle.write(json.dumps(row) + '\n')

    started = time.monotonic()
    try:
        artifact = run_candidate_search(attempts, output / 'search', generate, validate,
                        source_sha256=source.sha256, max_attempts=len(attempts),
                        deadline_seconds=900, on_progress=progress)
        result = {'status': 'static_checks_passed', 'coordinates': str(artifact.path),
                  'coordinates_sha256': artifact.sha256, 'full_preparation_published': False,
                  'app_default_enabled': False}
    except SearchExhausted as exc:
        result = {'status': 'exhausted', 'error': str(exc), 'full_preparation_published': False}
    except BaseException as exc:
        result = {'status': 'failed', 'error_type': type(exc).__name__, 'error': str(exc),
                  'full_preparation_published': False}
    finally:
        verify()
    result['elapsed_seconds'] = time.monotonic() - started
    write(output / 'result.json', result)
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('snapshot'); parser.add_argument('snapshot_sha256'); parser.add_argument('output')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--context-only', action='store_true', help='Only the bounded extra-context proposal comparison')
    parser.add_argument('--worker', help='Explicit private native worker source, copied and hashed for a bounded comparison')
    args = parser.parse_args()
    result = run(args.snapshot, args.snapshot_sha256, args.output, args.seed, args.context_only, args.worker)
    raise SystemExit(0 if result['status'] == 'static_checks_passed' else 2)
