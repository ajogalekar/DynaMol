"""Bounded research seed/refinement comparison, preserving every rejected case."""
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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def bounded(command, folder, label, environment):
    """Supervise one process group, recording its actual outcome and sampled RSS."""
    started = time.monotonic()
    peak, stop = 0, None
    with (folder / (label + '.log')).open('w') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   env=environment, start_new_session=True)
        monitor = psutil.Process(process.pid)
        try:
            while process.poll() is None:
                try:
                    rss = sum(p.memory_info().rss for p in [monitor, *monitor.children(recursive=True)]
                              if p.is_running())
                    peak = max(peak, rss)
                except psutil.NoSuchProcess:
                    pass
                if time.monotonic() - started > 120 or peak > 4 * 1024**3:
                    stop = 'wall_time_limit' if time.monotonic() - started > 120 else 'memory_limit'
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                    break
                time.sleep(.2)
        except BaseException:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            raise
        code = process.wait()
    record = {'command': command, 'returncode': code, 'stop_reason': stop,
              'elapsed_seconds': time.monotonic() - started, 'sampled_peak_rss_bytes': peak,
              'limits': {'seconds': 120, 'rss_bytes': 4 * 1024**3, 'threads': 2}}
    write(folder / (label + '-supervisor.json'), record)
    if code or stop:
        raise RuntimeError(label + ' failed; see preserved log and supervisor')
    return record


def main(base, output, repo, native_python, rama_runtime):
    scripts = Path(__file__).resolve().parent
    worker = base / '6jxt-context-v1/worker.py'
    first = output / 'native/context/output.json'
    # The first enumeration defines the entire trial plan before refinement.
    # It has already been saved by compare_fragment_seeds.py in a fresh folder.
    baseline = json.loads(first.read_text())
    searches = baseline['construction']['attempts'][0]['gap_searches']
    plan = [{'name': 'native', 'choices': {'default': {'ranking': 'native', 'rank': 0}}},
            {'name': 'anchor', 'choices': {'default': {'ranking': 'anchor', 'rank': 0}}}]
    for gap_index, gap in enumerate(searches):
        for rank in range(1, min(3, gap['candidate_count'])):
            plan.append({'name': f'gap-{gap_index + 1}-native-{rank}',
                         'choices': {'default': {'ranking': 'native', 'rank': 0},
                                     gap['gap']: {'ranking': 'native', 'rank': rank}}})
    sources = [Path(__file__), scripts / 'compare_fragment_seeds.py', scripts / 'check_context_seed.py',
               scripts / 'rama_reference.py', worker, base / '6jxt-context-v1/input.json',
               base / '6jxt-context-v1/input.pdb', base / '6jxt-protein-v4/before-refinement.pdb',
               base / '6jxt-protein-v4/retained-environment.pdb',
               repo / 'backend/loop_refinement.py', repo / 'backend/loop_geometry.py',
               repo / 'backend/preparation_worker.py', first]
    source_hashes = {str(p): digest(p) for p in sources}
    if (output / 'plan.json').exists():
        raise ValueError('The comparison folder already has a frozen trial plan')
    write(output / 'plan.json', {'trials': plan, 'source_sha256': source_hashes,
        'purpose': 'Compare a bounded set of native fragment alternatives under the unchanged previous conformation-protection refinement. No parameter fitting or acceptance cutoff changes.',
        'app_ready': False, 'physical_model_validated': False})
    snapshot = output / 'implementation-snapshot'
    snapshot.mkdir()
    for p in sources:
        if p.suffix == '.py':
            (snapshot / p.name).write_bytes(p.read_bytes())
    env = os.environ.copy()
    env.update(PYTHONPATH=str(rama_runtime), OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2',
               DYNAMOL_CPU_THREADS='2', PM3_OPENMM_CPU_THREADS='2')
    native_env = {k: v for k, v in env.items() if k != 'PYTHONPATH'}
    records = []
    started = time.monotonic()
    for trial in plan:
        if time.monotonic() - started > 1200:
            raise RuntimeError('Overall 20-minute comparison budget exhausted')
        folder = output / trial['name']
        folder.mkdir(exist_ok=trial['name'] == 'native')
        choice_file = folder / 'trial-choices.json'
        write(choice_file, trial['choices'])
        write(output / 'progress.json', {'active_trial': trial['name'], 'completed': records})
        try:
            if trial['name'] != 'native':
                bounded([str(native_python), str(scripts / 'compare_fragment_seeds.py'),
                    '--worker', str(worker), '--source', str(base / '6jxt-context-v1/input.json'),
                    '--output', str(folder / 'context'), '--choices', str(choice_file)], folder, 'worker', native_env)
            bounded([sys.executable, str(scripts / 'check_context_seed.py'), '--repo', str(repo),
                '--original', str(base / '6jxt-protein-v4'), '--context', str(folder / 'context'),
                '--output', str(folder / 'refinement'), '--environment', str(base / '6jxt-protein-v4/retained-environment.pdb'),
                '--backbone-construction', '--preserve-backbone-conformation'], folder, 'refinement', env)
            result = json.loads((folder / 'refinement/result.json').read_text())
            for stage, filename in [('seed', 'coherent-context-seed.pdb'), ('refined', 'refined-candidate.pdb')]:
                bounded([sys.executable, str(scripts / 'rama_reference.py'),
                    '--pdb', str(folder / 'refinement' / filename),
                    '--selection-result', str(folder / 'refinement/result.json'),
                    '--output', str(folder / (stage + '-rama'))], folder, stage + '-rama', env)
            seed = json.loads((folder / 'seed-rama/result.json').read_text())
            refined = json.loads((folder / 'refined-rama/result.json').read_text())
            record = {'name': trial['name'], 'completed': True, 'result': str(folder / 'refinement/result.json'),
                      'source_sha256': {str(folder / 'refinement/result.json'): digest(folder / 'refinement/result.json')},
                      'geometry_accepted': result['geometry_accepted'], 'stereo_error': result['stereo_error'],
                      'maximum_observed_context_displacement_nm': result['maximum_observed_context_displacement_nm'],
                      'seed_outliers': [r['residue'] for r in seed['outliers']],
                      'refined_outliers': [r['residue'] for r in refined['outliers']],
                      'selected_residue_count': len(refined['rows']),
                      'passed_local_and_saved_reference_checks': bool(result['passed_local_checks'] and refined['all_selected_scored_without_outliers']),
                      'app_ready': False, 'physical_model_validated': False}
        except Exception as exc:
            record = {'name': trial['name'], 'completed': False, 'error': str(exc),
                      'passed_local_and_saved_reference_checks': False,
                      'app_ready': False, 'physical_model_validated': False}
        records.append(record)
        write(output / 'progress.json', {'active_trial': None, 'completed': records})
        print(json.dumps(record), flush=True)
    if source_hashes != {name: digest(Path(name)) for name in source_hashes}:
        raise ValueError('Comparison sources changed during execution')
    result = {'stage': 'bounded_native_fragment_comparison_complete', 'trials': records,
              'trial_count': len(records), 'all_trials_completed': all(r['completed'] for r in records),
              'passing_local_candidate_count': sum(r['passed_local_and_saved_reference_checks'] for r in records),
              'source_sha256': source_hashes, 'elapsed_seconds': time.monotonic() - started,
              'app_ready': False, 'physical_model_validated': False,
              'scope': 'Research local preparation checks only. Independent environment/PDB roundtrip, native assembly and dynamics remain necessary before admission.'}
    write(output / 'result.json', result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['base', 'output', 'repo', 'native-python', 'rama-runtime']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    main(args.base, args.output, args.repo, args.native_python, args.rama_runtime)
