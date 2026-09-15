"""Bounded private replay of frozen proposals through complete-complex checks."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import psutil

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from backend.loop_search import Attempt, SearchExhausted, Validation, ValidatedArtifact, run_candidate_search

CACHE = Path.home() / '.cache/dynamol-research/v1-loop-release'
PYTHON = Path.home() / '.cache/dynamol-runtimes/covalent-v1-focused/bin/python'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def supervise(command, folder, checkpoint):
    env = {k: v for k, v in os.environ.items() if not k.startswith(('PYTHON', 'DYLD_', 'OPENMM', 'PM3_'))}
    env.update(PYTHONPATH=str(Path.home() / '.cache/dynamol-runtimes/cctbx-rama-reference-v1'),
               OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2', DYNAMOL_CPU_THREADS='2', OPENMM_CPU_THREADS='2')
    start = time.monotonic()
    peak = 0
    error = None
    with (folder / 'worker.log').open('w') as log:
        proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True)
        monitor = psutil.Process(proc.pid)
        try:
            while proc.poll() is None:
                checkpoint()
                try:
                    peak = max(peak, sum(p.memory_info().rss for p in [monitor, *monitor.children(recursive=True)] if p.is_running()))
                except psutil.NoSuchProcess:
                    pass
                if time.monotonic() - start > 180:
                    raise TimeoutError('Native complete-candidate validation exceeded 180 seconds')
                if peak > 4 * 1024**3:
                    raise MemoryError('Native complete-candidate validation exceeded 4 GiB')
                time.sleep(.2)
            if proc.returncode:
                raise RuntimeError('Native validator failed; inspect worker.log. No quality result was inferred from an exception.')
        except BaseException as exc:
            error = f'{type(exc).__name__}: {exc}'
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGKILL)
            raise
        finally:
            code = proc.wait()
            save(folder / 'supervisor.json', {'returncode': code, 'error': error, 'elapsed_seconds': time.monotonic() - start,
                 'peak_sampled_rss_bytes': peak, 'limits': {'seconds': 180, 'rss_bytes': 4 * 1024**3, 'threads': 2}, 'command': command})


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    implementation = output / 'implementation'
    implementation.mkdir()
    for source in [Path(__file__), Path(__file__).with_name('validate_complete_candidate.py'), ROOT / 'backend/loop_search.py']:
        (implementation / source.name).write_bytes(source.read_bytes())
    plans = []
    for name, dataset_id in [('1UA2', '4dcb5c25f5c14559'), ('8K5R', '4469bcb4878a4134')]:
        folder = output / name
        folder.mkdir()
        baseline = folder / 'archive-baseline.json'
        save(baseline, {'archive_baseline': True, 'app_ready': False})
        entries = [{'name': 'archive-fragment', 'generator': 'archived_fragment_seed', 'candidate': str(baseline)}]
        if name == '1UA2':
            for mode in ['torsion_mc', 'fragment_mc']:
                for seed in [2026, 2027]:
                    candidate = CACHE / 'hybrid-native-mc-v1/1UA2' / mode / f'native/43-56-seed-{seed}/candidate.json'
                    entries.append({'name': f'{mode}-{seed}', 'generator': 'promod3_' + mode, 'candidate': str(candidate)})
        candidate = CACHE / f'disgro-evaluation-v2/{name}' / ('candidate-003.json' if name == '1UA2' else 'candidate-005.json')
        entries.append({'name': 'disgro-representative', 'generator': 'private_pydisgro', 'candidate': str(candidate)})
        dataset = ROOT / f'docs/audit/loop-fallback/runs/root-final/{name}/workspace/datasets/{dataset_id}'
        files = [dataset / p for p in ['loop-refinement-topology.cif', 'loop-refinement-input.npz', 'loop-refinement-parameters.json']]
        frozen = {str(p): sha(p) for p in files + [Path(e['candidate']) for e in entries]}
        source_sha256 = hashlib.sha256(json.dumps({str(p): frozen[str(p)] for p in files}, sort_keys=True).encode()).hexdigest()
        plans.append({'case': name, 'dataset': str(dataset), 'attempts': entries, 'frozen_inputs': frozen, 'source_sha256': source_sha256})
    plan = {'cases': plans, 'selection': 'Order frozen before full-complex refinement. All original-stem native 1UA2 candidates, plus one representative per DiSGro case selected from prior local screens. This is a feasibility probe, not a success-rate estimate.',
            'limits': {'maximum_attempts': 8, 'seconds_per_attempt': 180, 'seconds_per_case': 1200, 'threads': 2, 'rss_bytes': 4 * 1024**3},
            'scope': 'Replays archived native parameters and preserves complete complexes. All modeled archive loops are checked, even where the proposal replaces only one loop. No new MD, QM, full prep, or release admission.',
            'app_ready': False, 'implementation_sha256': {str(p): sha(p) for p in implementation.iterdir()}}
    save(output / 'plan.json', plan)
    reports = []
    for case in plans:
        candidates = {e['name']: Path(e['candidate']) for e in case['attempts']}
        def generate(attempt, folder, checkpoint):
            checkpoint()
            source = candidates[attempt.name]
            if sha(source) != case['frozen_inputs'][str(source)]:
                raise ValueError('Frozen proposal changed')
            path = folder / 'proposal.json'
            path.write_bytes(source.read_bytes())
            return path
        def validate(candidate, folder, checkpoint):
            command = [str(PYTHON), '-u', str(Path(__file__).with_name('validate_complete_candidate.py')),
                       '--dataset', case['dataset'], '--candidate', str(candidate), '--output', str(folder / 'validation')]
            supervise(command, folder, checkpoint)
            path = folder / 'validation/result.json'
            result = json.loads(path.read_text())
            artifact = ValidatedArtifact(path, sha(path)) if all(v is True for v in result['checks'].values()) else None
            return Validation(result['checks'], tuple(result['reasons']), artifact=artifact)
        def progress(event):
            print(json.dumps({'case': case['case'], **event}), flush=True)
        try:
            result = run_candidate_search([Attempt(e['name'], e['generator']) for e in case['attempts']],
                      output / case['case'] / 'search', generate, validate, source_sha256=case['source_sha256'],
                      max_attempts=8, deadline_seconds=1200, on_progress=progress)
            status = 'accepted_by_private_checks'
        except SearchExhausted:
            result = None
            status = 'exhausted'
        except Exception as exc:
            reports.append({'case': case['case'], 'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'})
            save(output / 'summary.json', reports)
            raise
        if any(sha(Path(p)) != h for p, h in case['frozen_inputs'].items()):
            raise ValueError('A frozen comparison input changed')
        reports.append({'case': case['case'], 'status': status, 'artifact': str(result.path) if result else None, 'app_ready': False})
        save(output / 'summary.json', reports)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
