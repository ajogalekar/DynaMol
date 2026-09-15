"""Run the frozen representative panel with two bounded preparations at a time."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

REPO = Path('/Users/ashujo/Documents/Science/DynaMol')
HERE = REPO / 'docs/audit/v1-loop-release/fresh-source-integration'
PYTHON = '/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python'
BASE = Path('/Users/ashujo/.cache/dynamol-research/v1-loop-release')
STAGED = BASE / 'frozen-panel-raw-inputs-v1'
ROOT = BASE / 'frozen-panel-fresh-validation-v1'
OVERLAY = BASE / 'fresh-preparation-protected-flanks-v1'
ORDER = ['1M17', '8K5R', '1HCK', '4HJO', '2HYY', '3HEG', '2CG9',
         '1FPU', '4HHY', '6A93', '2BEL', '2REN', '1RNE', '2ITY']
EXPECTED_REFUSALS = {'6A93', '2BEL', '2REN', '1RNE', '2ITY'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def stage(case, entry):
    root = ROOT / case
    root.mkdir()
    destination = root / 'workspace/datasets' / entry['dataset_id']
    original = Path(entry['dataset_path'])
    if sha(entry['staging_manifest']) != entry['staging_manifest_sha256']:
        raise ValueError('Staged case manifest changed: ' + case)
    records = json.loads(Path(entry['staging_manifest']).read_text())
    raw = []
    for row in records['files']:
        source = Path(row['copy'])
        target = destination / source.relative_to(original)
        if sha(source) != row['sha256']:
            raise ValueError('Staged raw source changed: ' + str(source))
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        if sha(target) != row['sha256']:
            raise ValueError('Raw copy mismatch')
        raw.append({'source': str(source), 'copy': str(target), 'sha256': row['sha256']})
    write(root / 'raw-source-provenance.json', {'staging_manifest': entry['staging_manifest'],
          'staging_sha256': entry['staging_manifest_sha256'], 'files': raw})
    overlay = []
    for row in json.loads((OVERLAY / 'overlay-provenance.json').read_text())['files']:
        source = Path(row['copy'])
        target = root / 'overlay' / source.relative_to(OVERLAY / 'overlay')
        if sha(source) != row['sha256']:
            raise ValueError('Qualified pure Python overlay changed')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        if sha(target) != row['sha256']:
            raise ValueError('Overlay copy mismatch')
        overlay.append({'source': str(source), 'copy': str(target), 'sha256': row['sha256']})
    write(root / 'overlay-provenance.json', {'purpose': 'Pinned pure Python dependencies only', 'files': overlay})
    config = {'case': case, 'root': str(root), 'dataset_id': entry['dataset_id'],
              'expected_raw_atoms': entry['n_atoms'], 'seed': records['settings']['seed']}
    write(root / 'input-configuration.json', config)
    return root, config


def one(case, entry):
    started = time.monotonic()
    print(json.dumps({'case': case, 'stage': 'fresh_preparation'}), flush=True)
    root, config = stage(case, entry)
    env = dict(os.environ, PYTHONPATH=str(REPO), OPENBLAS_NUM_THREADS='2',
               OPENMM_CPU_THREADS='2', DYNAMOL_CPU_THREADS='2', OMP_NUM_THREADS='2',
               PYTHONDONTWRITEBYTECODE='1')
    capture = [PYTHON, '-B', str(HERE / 'panel_capture.py'), '--config', str(root / 'input-configuration.json')]
    with (root / 'capture-launch.log').open('x') as stream:
        code = subprocess.run(capture, env=env, stdout=stream, stderr=subprocess.STDOUT).returncode
    result = {'case': case, 'expected': 'unsupported_input' if case in EXPECTED_REFUSALS else 'repair_candidate',
              'capture_returncode': code, 'full_preparation_published': False, 'app_default_enabled': False}
    if code == 0:
        captured = json.loads((root / 'capture-result.json').read_text())
        result['snapshot_sha256'] = captured['snapshot_sha256']
        print(json.dumps({'case': case, 'stage': 'complete_candidate_checks'}), flush=True)
        command = [PYTHON, '-B', str(HERE / 'run_snapshot_search.py'), str(root / 'snapshot'),
                   captured['snapshot_sha256'], str(root / 'candidate-search'), '--seed', str(config['seed'])]
        with (root / 'search-launch.log').open('x') as stream:
            code = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT).returncode
        result['search_returncode'] = code
        if (root / 'candidate-search/result.json').exists():
            result['search'] = json.loads((root / 'candidate-search/result.json').read_text())
            result['status'] = result['search']['status']
        else:
            result['status'] = 'search_harness_failure'
    else:
        result['status'] = 'preparation_failed'
        for name in ['capture-failure.json', 'early-failure.json']:
            if (root / name).is_file():
                result['failure'] = json.loads((root / name).read_text())
                break
    result['elapsed_seconds'] = time.monotonic() - started
    write(root / 'case-result.json', result)
    print(json.dumps({'case': case, 'status': result['status'], 'seconds': result['elapsed_seconds']}), flush=True)
    return result


def main():
    manifest_path = STAGED / 'staging-manifest.json'
    if sha(manifest_path) != '1d85de069830a7312311121405c63db53dc15c91bae93a4563c4d703056033fc':
        raise ValueError('Frozen input panel changed')
    manifest = json.loads(manifest_path.read_text())
    ROOT.mkdir(exist_ok=False)
    shutil.copyfile(__file__, ROOT / 'run_panel.py')
    write(ROOT / 'plan.json', {'cases': ['1UA2', *ORDER], 'parallel_cases': 2,
          'expected_unsupported_cases': sorted(EXPECTED_REFUSALS),
          'settings': 'Frozen per-case seed42, complete retained selected complex, default preparation settings',
          'first_1ua2_seed': 2026, 'first_1ua2_evidence': str(BASE / 'fresh-protected-loop-search-v1'),
          'capture_seconds_per_case': 600, 'search_seconds_per_case': 900,
          'scope': 'Fresh preparation and complete static loop validation; app publication remains gated',
          'MD': False, 'ordinary_ligand_AM1_BCC': True})
    started = time.monotonic()
    results = {'1UA2': {'case': '1UA2', 'expected': 'repair_candidate',
                'status': 'static_checks_passed', 'evidence': str(BASE / 'fresh-protected-loop-search-v1/result.json'),
                'snapshot': str(BASE / 'fresh-preparation-protected-flanks-v1/snapshot'),
                'seed': 2026, 'full_preparation_published': False}}
    write(ROOT / 'results.json', {'status': 'running', 'cases': results})
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(one, case, manifest['cases'][case]): case for case in ORDER}
        for future in as_completed(futures):
            case = futures[future]
            try:
                results[case] = future.result()
            except Exception as exc:
                results[case] = {'case': case, 'status': 'harness_failed',
                                 'error_type': type(exc).__name__, 'error': str(exc)}
                print(json.dumps(results[case]), flush=True)
            write(ROOT / 'results.json', {'status': 'running', 'cases': results,
                                         'elapsed_seconds': time.monotonic() - started})
    write(ROOT / 'results.json', {'status': 'complete', 'cases': results,
                                 'elapsed_seconds': time.monotonic() - started})
    print(json.dumps({'status': 'panel_complete', 'results': str(ROOT / 'results.json')}), flush=True)


if __name__ == '__main__':
    main()
