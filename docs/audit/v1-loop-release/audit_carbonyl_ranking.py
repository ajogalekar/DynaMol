"""Replay frozen loop inputs through the corrected native candidate ranking.

This is a construction-stage audit, not full preparation or release admission.
Historical artifacts are read only; each native worker gets private inputs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import numpy as np
import psutil
from openmm import app, unit

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'docs/audit/advanced-chemistry/covalent-v1-focused'))
from backend.loop_modeling import _transplant
from rama_reference import rama_report, residue_key


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def native(command, directory):
    env = {k: v for k, v in os.environ.items() if not k.startswith(('PYTHON', 'DYLD_', 'OPENMM', 'PM3_'))}
    env.update(PM3_OPENMM_CPU_THREADS='2', OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2')
    started = time.monotonic()
    peak = 0
    reason = None
    with (directory / 'worker.log').open('w') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   env=env, start_new_session=True)
        monitor = psutil.Process(process.pid)
        try:
            while process.poll() is None:
                try:
                    rss = sum(p.memory_info().rss for p in [monitor, *monitor.children(recursive=True)] if p.is_running())
                    peak = max(peak, rss)
                except psutil.NoSuchProcess:
                    pass
                if time.monotonic() - started > 180 or peak > 4 * 1024**3:
                    reason = 'wall_time' if time.monotonic() - started > 180 else 'memory'
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=3)
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
    result = {'returncode': code, 'stop_reason': reason,
              'elapsed_seconds': time.monotonic() - started, 'peak_sampled_rss_bytes': peak,
              'limits': {'seconds': 180, 'rss_bytes': 4 * 1024**3, 'threads': 2}, 'command': command}
    save(directory / 'supervisor.json', result)
    if code or reason:
        raise RuntimeError('Native construction failed; see retained supervisor and log')
    return result


def main(output, native_python):
    output.mkdir(parents=True, exist_ok=False)
    worker = ROOT / 'backend/promod_loop_worker.py'
    sources = [Path(__file__), worker, ROOT / 'backend/loop_modeling.py',
               ROOT / 'docs/audit/advanced-chemistry/covalent-v1-focused/rama_reference.py']
    snapshot = output / 'implementation'
    snapshot.mkdir()
    for p in sources:
        (snapshot / p.name).write_bytes(p.read_bytes())
    frozen = ROOT / 'docs/audit/loop-fallback'
    plan = []
    for case in ['8K5R', '1UA2']:
        report = frozen / 'runs/root-final' / case / 'result.json'
        record = json.loads(report.read_text())
        dataset = report.parent / 'workspace/datasets' / record['prepared_dataset_id']
        original_request = dataset / 'loop-model-input.json'
        request = json.loads(original_request.read_text())
        context = Path(request['input_pdb'])
        inputs = [report, original_request, context, dataset / 'loop-model-output.json',
                  dataset / 'loop-construction.pdb', dataset / 'loop-construction-before.npz',
                  dataset / 'loop-construction-after.npz']
        sources.extend(inputs)
        plan.append({'case': case, 'dataset': str(dataset), 'request': request,
                     'context': str(context)})
    hashes = {str(p): digest(p) for p in sources}
    save(output / 'plan.json', {'cases': plan, 'source_sha256': hashes,
        'scope': 'Native construction replay only; original observed coordinates are restored before static reference scoring. No minimization, full preparation, MD or release admission.'})
    rows = []
    for item in plan:
        case = item['case']
        directory = output / case
        directory.mkdir()
        try:
            context = directory / 'context.pdb'
            context.write_bytes(Path(item['context']).read_bytes())
            request = dict(item['request'], input_pdb=str(context), research_include_carbonyl_anchors=True)
            request.pop('_progress_path', None)
            save(directory / 'input.json', request)
            supervisor = native([str(native_python), '-I', '-B', str(snapshot / worker.name),
                                 str(directory / 'input.json'), str(directory / 'output.json')], directory)
            result = json.loads((directory / 'output.json').read_text())
            dataset = Path(item['dataset'])
            previous = json.loads((dataset / 'loop-model-output.json').read_text())
            old_keys = [tuple(a['identity']) for a in previous['atoms']]
            new_keys = [tuple(a['identity']) for a in result['atoms']]
            if len(set(old_keys)) != len(old_keys) or set(old_keys) != set(new_keys) or len(new_keys) != len(old_keys):
                raise ValueError('Candidate heavy-atom identities changed')
            scaffold = app.PDBFile(str(dataset / 'loop-construction.pdb'))
            with np.load(dataset / 'loop-construction-before.npz') as array:
                before = array['xyz_nm'].copy()
            with np.load(dataset / 'loop-construction-after.npz') as array:
                historical = array['xyz_nm'].copy()
            target = {(*residue_key(a.residue), a.name): a.index for a in scaffold.topology.atoms()
                      if (*residue_key(a.residue), a.name) in set(new_keys)}
            updated = _transplant(result, target, before)
            fixed = [a.index for a in scaffold.topology.atoms() if a.index not in target.values()]
            if not np.array_equal(updated[fixed], before[fixed]):
                raise ValueError('Construction transplant moved observed atoms')
            modeled = {key[:4] for key in new_keys}
            boundary = set()
            for a, b in scaffold.topology.bonds():
                if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
                    continue
                pair = {residue_key(a.residue), residue_key(b.residue)}
                if pair & modeled:
                    boundary |= pair - modeled
            selected = modeled | boundary
            references = {}
            for label, xyz in [('historical', historical), ('corrected', updated)]:
                report = rama_report(scaffold.topology, xyz * unit.nanometer, selected)
                save(directory / (label + '-reference.json'), report)
                references[label] = {'modeled_outliers': [r['residue'] for r in report['outliers'] if tuple(r['residue']) in modeled],
                                     'boundary_outliers': [r['residue'] for r in report['outliers'] if tuple(r['residue']) in boundary],
                                     'unavailable': report['unavailable']}
            with (directory / 'candidate-before-refinement.pdb').open('w') as handle:
                app.PDBFile.writeFile(scaffold.topology, updated * unit.nanometer, handle, keepIds=True)
            np.savez(directory / 'candidate-before-refinement.npz', xyz_nm=updated)
            searches = [search for attempt in result['construction']['attempts']
                        for search in attempt.get('gap_searches', [])]
            if any(s['candidate_count'] > 40 or s['anchor_atoms'] != ['N', 'CA', 'C', 'O'] for s in searches):
                raise ValueError('Corrected candidate ranking/budget was not used')
            rows.append({'case': case, 'completed': True, 'native_seconds': supervisor['elapsed_seconds'],
                         'candidate_heavy_atoms': len(new_keys), 'modeled_residues': sorted(modeled),
                         'boundary_residues': sorted(boundary), 'references': references,
                         'unchanged_observed_atoms': len(fixed), 'context_searches': searches,
                         'prepared_model_accepted': False})
        except Exception as exc:
            rows.append({'case': case, 'completed': False, 'error_type': type(exc).__name__, 'error': str(exc)})
        save(output / 'progress.json', {'completed': rows})
    if hashes != {name: digest(Path(name)) for name in hashes}:
        raise ValueError('A replay source changed during the audit')
    result = {'cases': rows, 'source_sha256': hashes, 'sources_unchanged': True,
              'app_ready': False, 'physical_model_validated': False,
              'scope': 'Native candidate construction and reference scoring only. Full-atom geometry, stereochemistry, environment and final refinement still required.'}
    save(output / 'result.json', result)
    print(json.dumps({'cases': [{k: v for k, v in row.items() if k != 'context_searches'} for row in rows]}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--native-python', type=Path, required=True)
    args = parser.parse_args()
    main(args.output.resolve(), args.native_python)
