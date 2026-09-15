"""Contract tests use explicit synthetic fixtures; no native calculation on import.

Optional native qualification:
  isolated-python test_geometric_adapter.py --native --root REPO --output NEWDIR
Two tiny H2O2 calculations share a 600-second total deadline. They exercise exact
Cartesian freezes plus a periodic dihedral; neither supports full-adduct accuracy.
"""
from copy import deepcopy
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location('qualification_geometric_adapter', Path(__file__).with_name('geometric_adapter.py'))
adapter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adapter)


def peroxide(angle=60.0, target=60.0):
    ids = ['H1', 'O1', 'O2', 'H2']
    return {'schema_version': 1, 'scope': adapter.SCOPE,
        'calculator': {'schema_version': 1, 'method': 'GFN2-xTB',
            'atom_ids': ids, 'elements': ['H', 'O', 'O', 'H'], 'charge': 0, 'spin': 0,
            'coords_bohr': [[0., 1.8, 0.], [0., 0., 0.], [2.7, 0., 0.],
                            [2.7, 1.8*math.cos(math.radians(angle)), 1.8*math.sin(math.radians(angle))]]},
        'identity': {'bonds': [{'atom_ids': ids[i:i+2], 'order': 1} for i in range(3)],
                     'stereocenters': []},
        'optimization': {'max_steps': 80, 'max_evaluations': 200, 'max_wall_seconds': 120,
            'freeze_atom_ids': ['O1', 'O2'], 'dihedral_constraints': [
                {'atom_ids': ids, 'target_degrees': target}], 'dihedral_tolerance_degrees': .1}}


def synthetic_evaluator(calc):
    return {'accepted': True, 'physical_acceptance': False,
        **{key: deepcopy(calc[key]) for key in ('atom_ids', 'elements', 'charge', 'spin', 'method', 'coords_bohr')},
        'energy_hartree': -1.25, 'gradient_hartree_per_bohr': [[.01, -.02, .03] for _ in calc['atom_ids']],
        'units': {'coordinates': 'Bohr', 'energy': 'Hartree', 'gradient': 'Hartree/Bohr; positive dE/dx'},
        'calculator_sha256': 'synthetic-test-only', 'runtime': {'synthetic': True},
        'native_version': ['synthetic'], 'settings': {'synthetic': True}}


def synthetic_runner(request, output, session):
    numerical, _ = session.evaluate(request['calculator']['coords_bohr'])
    # Explicitly confirm the adapter forwarded the raw E and positive gradient.
    assert numerical['energy'] == -1.25
    assert numerical['gradient'].reshape(-1, 3).tolist() == [[.01, -.02, .03]]*4
    return {'native_converged': True, 'iterations': 2,
            'coords_bohr': deepcopy(request['calculator']['coords_bohr'])}


def run(tmp_path, request=None, evaluator=synthetic_evaluator, runner=synthetic_runner):
    return adapter.run_synthetic_fixture(request or peroxide(), tmp_path/'case', evaluator, runner)


def test_forward_units_constraints_and_synthetic_label(tmp_path):
    result = run(tmp_path)
    assert result['accepted'] is True
    assert result['physical_acceptance'] is False
    assert result['simulation_ready'] is False
    assert result['evidence_kind'] == 'synthetic_engine_contract_fixture'
    assert result['evaluations'] == 2  # final coordinates explicitly re-evaluated
    assert result['added_restraint_energy_hartree'] == 0
    assert result['convergence_threshold_units'] == {'convergence_energy': 'Hartree',
        'convergence_grms': 'Hartree/Bohr', 'convergence_gmax': 'Hartree/Bohr',
        'convergence_drms': 'Angstrom', 'convergence_dmax': 'Angstrom'}
    request = adapter.validate(peroxide(target=-300))
    assert adapter.constraint_text(request) == '$freeze\nxyz 2,3\n$set\ndihedral 1 2 3 4 60\n'
    assert len(result['input_sha256']) == 64
    assert result['geometry_checks']['dihedrals'][0]['satisfied']


@pytest.mark.parametrize('measured,target', [(179., -181.), (-179., 181.), (60., -300.)])
def test_periodic_dihedral_comparison(measured, target):
    request = adapter.validate(peroxide(measured, target))
    checks = adapter.verify_geometry(request, request['calculator']['coords_bohr'])
    assert checks['passed']
    assert abs(checks['dihedrals'][0]['signed_periodic_error_degrees']) < 1e-12


@pytest.mark.parametrize('mutation,message', [
    (lambda r: r.update(extra=True), 'Unknown or missing'),
    (lambda r: r.update(scope='full_adduct_geometry'), 'Only explicit'),
    (lambda r: r['calculator'].update(method='GFN1-xTB'), 'Explicit GFN2'),
    (lambda r: r['calculator'].update(spin=2), 'closed-shell'),
    (lambda r: r['calculator'].update(charge=1), 'even electron'),
    (lambda r: r['calculator'].update(charge=True), 'integer charge'),
    (lambda r: r['calculator'].update(coords_bohr=[[0, 0, float('nan')]]*4), 'Finite explicit'),
    (lambda r: r['calculator'].update(coords_bohr=[[0., 0., 0.]]*4), 'Coincident'),
    (lambda r: r['calculator'].update(atom_ids=['H1', 'O1', 'O2', 'O1']), 'unique atom IDs'),
    (lambda r: r['calculator'].update(elements=['H', 'Zn', 'O', 'H']), 'H/C/N/O/F/S'),
    (lambda r: r['calculator'].update(max_iter=1000), 'iteration limit'),
    (lambda r: r['calculator'].update(accuracy=float('inf')), 'accuracy'),
    (lambda r: r['optimization'].update(max_wall_seconds=601), 'bounded'),
    (lambda r: r['optimization'].update(max_steps=True), 'bounded'),
    (lambda r: r['optimization'].update(convergence_grms=.5), 'Unknown optimization'),
    (lambda r: r['optimization'].update(freeze_atom_ids=['O1', 'O1']), 'uniquely'),
    (lambda r: r['optimization'].update(freeze_atom_ids=['H1', 'O1', 'O2', 'H2']), 'every atom'),
    (lambda r: r['optimization'].update(dihedral_tolerance_degrees=1.), 'at most 0.5'),
    (lambda r: r['optimization']['dihedral_constraints'].append(
        {'atom_ids': ['H2', 'O2', 'O1', 'H1'], 'target_degrees': 40}), 'Duplicate or reversed'),
    (lambda r: r['identity']['bonds'].append({'atom_ids': ['O1', 'H1'], 'order': 1}), 'Duplicate bond'),
    (lambda r: r['identity']['bonds'].pop(), 'connected'),
])
def test_malformed_input_rejected_before_evaluation(tmp_path, mutation, message):
    request = peroxide()
    mutation(request)
    def forbidden(_):
        pytest.fail('Malformed input reached evaluator')
    result = run(tmp_path, request, evaluator=forbidden)
    assert not result['accepted']
    assert message in result['error']
    assert result['evaluations'] == 0
    assert (tmp_path/'case/failure-traceback.txt').exists()


def test_full_parent_size_is_not_admitted():
    request = peroxide()
    request['calculator']['atom_ids'] = [str(i) for i in range(94)]
    with pytest.raises(ValueError, match='full-parent admission'):
        adapter.validate(request)


def test_singular_dihedral_rejected():
    request = peroxide()
    request['calculator']['coords_bohr'][0] = [-1.8, 0., 0.]
    with pytest.raises(ValueError, match='Singular dihedral'):
        adapter.validate(request)


@pytest.mark.parametrize('key,value', [
    ('accepted', False), ('physical_acceptance', True), ('atom_ids', ['O1', 'H1', 'O2', 'H2']),
    ('elements', ['H', 'S', 'O', 'H']), ('charge', 2), ('spin', 2), ('method', 'GFN1-xTB'),
    ('spin', False),
    ('coords_bohr', [[0, 0, 0]]*4), ('energy_hartree', float('nan')),
    ('gradient_hartree_per_bohr', [[0, 0, float('inf')]]*4),
    ('gradient_hartree_per_bohr', [1, 2, 3]),
    ('units', {'coordinates': 'Bohr', 'energy': 'Hartree', 'gradient': 'force'}),
    ('calculator_sha256', None),
])
def test_bad_calculator_response_preserved_and_rejected(tmp_path, key, value):
    def evaluator(calc):
        result = synthetic_evaluator(calc)
        result[key] = value
        return result
    result = run(tmp_path, evaluator=evaluator)
    assert not result['accepted']
    assert result['evaluations'] == 1
    assert list((tmp_path/'case/evaluations').glob('0001-*result*'))


def test_runtime_binding_cannot_change(tmp_path):
    calls = 0
    def evaluator(calc):
        nonlocal calls
        calls += 1
        result = synthetic_evaluator(calc)
        result['runtime']['invocation'] = calls
        return result
    result = run(tmp_path, evaluator=evaluator)
    assert not result['accepted']
    assert 'binding changed' in result['error']


@pytest.mark.parametrize('change,message', [
    (lambda result: result.update(native_converged=False), 'did not converge'),
    (lambda result: result.update(iterations=81), 'step budget'),
    (lambda result: result['coords_bohr'][1].__setitem__(0, .001), 'Independent final'),
    (lambda result: result.update(coords_bohr=peroxide(angle=70)['calculator']['coords_bohr']), 'Independent final'),
    (lambda result: result['coords_bohr'][0].__setitem__(1, 3.0), 'Independent final'),
])
def test_native_return_is_not_enough_for_acceptance(tmp_path, change, message):
    def runner(request, output, session):
        result = synthetic_runner(request, output, session)
        change(result)
        return result
    result = run(tmp_path, runner=runner)
    assert not result['accepted']
    assert message in result['error']
    assert (tmp_path/'case/optimizer-result.json').exists()


def test_evaluation_budget_includes_final_single_point(tmp_path):
    request = peroxide()
    request['optimization']['max_evaluations'] = 1
    result = run(tmp_path, request)
    assert not result['accepted']
    assert 'evaluation budget' in result['error']


def test_deadline_checked_on_each_evaluation(tmp_path):
    def runner(request, output, session):
        session.deadline = time.monotonic()-1
        return synthetic_runner(request, output, session)
    result = run(tmp_path, runner=runner)
    assert not result['accepted']
    assert result['error_type'] == 'TimeoutError'
    assert result['evaluations'] == 0


def test_fresh_output_required(tmp_path):
    run(tmp_path)
    original = (tmp_path/'case/result.json').read_bytes()
    with pytest.raises(FileExistsError):
        run(tmp_path)
    assert (tmp_path/'case/result.json').read_bytes() == original


def test_stereo_sign_verification():
    ids = ['C', 'H', 'F', 'N', 'O']
    xyz = [[0, 0, 0], [1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]]
    volume = np.linalg.det(np.array(xyz[1:4])-xyz[4])
    definitions = [{'center_atom_id': 'C', 'ordered_neighbor_ids': ids[1:], 'signed_volume_sign': int(np.sign(volume))}]
    assert adapter._stereo_rows(xyz, ids, definitions)[0]['satisfied']
    xyz[1], xyz[2] = xyz[2], xyz[1]
    with pytest.raises(ValueError, match='inverted or degenerate'):
        adapter._stereo_rows(xyz, ids, definitions)


def test_native_timeout_kills_process_group_and_preserves_files(tmp_path, monkeypatch):
    sources = tmp_path/'provider'; sources.mkdir()
    (sources/'calculator.py').write_text('synthetic supervisor fixture')
    (sources/'runtime-manifest.json').write_text('{}')
    monkeypatch.setattr(adapter, 'HERE', sources)
    killed = []
    class Process:
        pid = 123456789
        def wait(self, timeout=None):
            if timeout is not None and timeout <= .1:
                raise subprocess.TimeoutExpired('synthetic-worker', timeout)
            return -9
    monkeypatch.setattr(adapter.subprocess, 'Popen', lambda *args, **kwargs: Process())
    monkeypatch.setattr(adapter.os, 'killpg', lambda pid, sig: killed.append((pid, sig)))
    monkeypatch.setattr(adapter, '_group_is_alive', lambda pid: False)
    request = peroxide();request['optimization']['max_wall_seconds'] = 1
    result = adapter.run_native(request, tmp_path/'bounded', tmp_path)
    assert result['error_type'] == 'TimeoutError'
    assert killed == [(123456789, signal.SIGTERM), (123456789, signal.SIGKILL)]
    assert (tmp_path/'bounded/launch-input.json').exists()
    assert (tmp_path/'bounded/native.log').exists()


def test_supervisor_rejects_provider_change_after_worker_result(tmp_path, monkeypatch):
    sources = tmp_path/'provider'; sources.mkdir()
    (sources/'calculator.py').write_text('synthetic supervisor fixture')
    (sources/'runtime-manifest.json').write_text('{}')
    monkeypatch.setattr(adapter, 'HERE', sources)
    output = tmp_path/'binding'
    class Process:
        pid = 123456789
        def wait(self, timeout=None):
            adapter.write_json(output/'result.json', {'accepted': True, 'physical_acceptance': False})
            (sources/'calculator.py').write_text('changed provider')
            return 0
    monkeypatch.setattr(adapter.subprocess, 'Popen', lambda *args, **kwargs: Process())
    result = adapter.run_native(peroxide(), output, tmp_path)
    assert not result['accepted']
    assert 'binding changed' in result['error']
    assert (output/'worker-result-before-supervisor-rejection.json').exists()


def test_supervisor_preserves_native_rejection_reason(tmp_path, monkeypatch):
    sources = tmp_path/'provider'; sources.mkdir()
    (sources/'calculator.py').write_text('synthetic supervisor fixture')
    (sources/'runtime-manifest.json').write_text('{}')
    monkeypatch.setattr(adapter, 'HERE', sources)
    output = tmp_path/'native-refusal'
    class Process:
        pid = 123456789
        def wait(self, timeout=None):
            adapter.write_json(output/'result.json', {'accepted': False, 'physical_acceptance': False,
                'status': 'rejected', 'error': 'Synthetic native nonconvergence'})
            return 1
    monkeypatch.setattr(adapter.subprocess, 'Popen', lambda *args, **kwargs: Process())
    result = adapter.run_native(peroxide(), output, tmp_path)
    assert not result['accepted']
    assert result['error'] == 'Synthetic native nonconvergence'


def _wait_for_file(path, process, timeout=4):
    deadline = time.monotonic()+timeout
    while time.monotonic() < deadline:
        if path.exists():
            return
        if process.poll() is not None:
            raise AssertionError('Lifecycle supervisor exited before '+path.name)
        time.sleep(.01)
    raise TimeoutError('Lifecycle test did not produce '+path.name)


@pytest.mark.parametrize('mode,expected_error', [
    ('timeout', 'TimeoutError'), ('artifact-write-error', 'OSError'),
    ('repeated-cancel', 'Cancelled'), ('leader-exits', None),
    ('cancel-during-spawn', 'Cancelled'),
    ('repeat-at-cleanup-entry', 'Cancelled'),
    ('first-at-cleanup-entry-normal', 'Cancelled'),
    ('first-at-cleanup-entry-error', 'OSError'),
])
def test_real_process_lifetime_and_descendants(tmp_path, mode, expected_error):
    """Ordinary sleeping Python processes only: no numerical/native calculator."""
    # Each sleeping child and grandchild ignores SIGTERM, so the test also
    # exercises escalation to SIGKILL and a live group after its leader exits.
    grandchild = tmp_path/'grandchild.py'
    grandchild.write_text('import os,signal,time\nfrom pathlib import Path\n'
        'signal.signal(signal.SIGTERM,signal.SIG_IGN)\n'
        f'Path({str(tmp_path/"grandchild.pid")!r}).write_text(str(os.getpid()))\n'
        'time.sleep(20)\n')
    child = tmp_path/'child.py'
    child.write_text('import os,signal,subprocess,sys,time,json\nfrom pathlib import Path\n'
        'signal.signal(signal.SIGTERM,signal.SIG_IGN)\n'
        f'grand=subprocess.Popen([sys.executable,{str(grandchild)!r}])\n'
        f'Path({str(tmp_path/"children.json")!r}).write_text(json.dumps({{"leader":os.getpid(),"grandchild":grand.pid}}))\n'
        + ('sys.exit(0)\n' if mode in ('leader-exits', 'first-at-cleanup-entry-normal') else 'time.sleep(20)\n'))
    controller = tmp_path/'controller.py'
    controller.write_text('import importlib.util,json,os,signal,subprocess,sys,time\nfrom pathlib import Path\nfrom contextlib import contextmanager\n'
        f'spec=importlib.util.spec_from_file_location("adapter",{str(Path(adapter.__file__).resolve())!r})\n'
        'adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)\n'
        f'base=Path({str(tmp_path)!r});mode={mode!r}\n'
        'real_cleanup=adapter._terminate_owned_group\n'
        'def cleanup(child,grace_seconds):\n'
        '    (base/"cleanup-started").write_text("yes")\n'
        '    return real_cleanup(child,grace_seconds)\n'
        'adapter._terminate_owned_group=cleanup\n'
        'if mode in ("repeat-at-cleanup-entry","first-at-cleanup-entry-normal","first-at-cleanup-entry-error"):\n'
        '    real_guard=adapter._uninterrupted_cleanup\n'
        '    @contextmanager\n'
        '    def boundary_guard():\n'
        '        os.kill(os.getpid(),signal.SIGINT if mode=="repeat-at-cleanup-entry" else signal.SIGTERM)\n'
        '        (base/"boundary-signal-injected").write_text("before cleanup handler installation")\n'
        '        with real_guard():yield\n'
        '    adapter._uninterrupted_cleanup=boundary_guard\n'
        'if mode=="cancel-during-spawn":\n'
        '    real_popen=adapter.subprocess.Popen\n'
        '    def popen(*args,**kwargs):\n'
        '        child=real_popen(*args,**kwargs)\n'
        '        if "start_new_session" in kwargs:\n'
        '            deadline=time.monotonic()+3\n'
        '            while not (base/"grandchild.pid").exists() and time.monotonic()<deadline:time.sleep(.01)\n'
        '            os.kill(os.getpid(),signal.SIGTERM)\n'
        '        return child\n'
        '    adapter.subprocess.Popen=popen\n'
        'def on_start(pid):\n'
        '    deadline=time.monotonic()+3\n'
        '    while not (base/"grandchild.pid").exists():\n'
        '        if time.monotonic()>deadline:raise TimeoutError("Child startup")\n'
        '        time.sleep(.01)\n'
        '    (base/"ready").write_text(str(pid))\n'
        '    if mode=="repeat-at-cleanup-entry":os.kill(os.getpid(),signal.SIGTERM)\n'
        '    if mode in ("artifact-write-error","first-at-cleanup-entry-error"):raise OSError("Synthetic process.json write failure")\n'
        'try:\n'
        f'    code=adapter._run_owned_process([sys.executable,{str(child)!r}],base/"child.log",.1 if mode=="timeout" else 5,on_start,dict(os.environ))\n'
        '    result={"error_type":None,"returncode":code}\n'
        'except BaseException as error:result={"error_type":type(error).__name__,"error":str(error)}\n'
        'children=json.loads((base/"children.json").read_text())\n'
        'result["live_group"]=adapter._group_is_alive(children["leader"])\n'
        'try:os.waitpid(children["leader"],os.WNOHANG);result["leader_reaped"]=False\n'
        'except ChildProcessError:result["leader_reaped"]=True\n'
        '(base/"result.json").write_text(json.dumps(result))\n')
    supervisor = subprocess.Popen([sys.executable, str(controller)], start_new_session=True)
    try:
        if mode == 'repeated-cancel':
            _wait_for_file(tmp_path/'ready', supervisor)
            os.kill(supervisor.pid, signal.SIGTERM)
            _wait_for_file(tmp_path/'cleanup-started', supervisor)
            # The child ignores TERM, providing a real cleanup window.
            for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGINT):
                os.kill(supervisor.pid, sig)
                time.sleep(.03)
        assert supervisor.wait(timeout=8) == 0
        result = json.loads((tmp_path/'result.json').read_text())
        assert result['error_type'] == expected_error
        assert result['leader_reaped'] is True
        assert result['live_group'] is False
        if mode in ('artifact-write-error', 'first-at-cleanup-entry-error'):
            assert 'process.json' in result['error']
        if 'cleanup-entry' in mode:
            assert (tmp_path/'boundary-signal-injected').exists()
        if mode == 'repeat-at-cleanup-entry':
            assert result['error'] == 'Supervisor cancellation signal 15'
    finally:
        # Fail-safe for the tests themselves; never leave a failed fixture alive.
        if (tmp_path/'children.json').exists():
            pid = json.loads((tmp_path/'children.json').read_text())['leader']
            try:os.killpg(pid, signal.SIGKILL)
            except ProcessLookupError:pass
        if supervisor.poll() is None:
            os.killpg(supervisor.pid, signal.SIGKILL)
        supervisor.wait(timeout=3)


def native_qualification(root, output):
    """Two four-atom gas-phase optimizations, total wall <=600 seconds."""
    output.mkdir(parents=True, exist_ok=False)
    deadline = time.monotonic()+600
    results = []
    # The second target crosses the -180/+180 periodic representation boundary.
    for name, angle, target in [('freeze-and-dihedral', 80., 60.), ('periodic-target', 160., -181.)]:
        remaining = math.floor(deadline-time.monotonic())
        if remaining <= 0:
            raise TimeoutError('Native qualification total 600-second budget exhausted')
        request = peroxide(angle, target)
        request['optimization']['max_wall_seconds'] = min(240, remaining)
        result = adapter.run_native(request, output/name, root)
        results.append({'name': name, 'accepted': result['accepted'], 'error': result.get('error'),
                        'result': str(output/name/'result.json')})
        if not result['accepted']:
            break
    summary = {'accepted': len(results) == 2 and all(r['accepted'] for r in results),
               'physical_acceptance': False, 'scope': 'two four-atom H2O2 constrained optimizations only',
               'total_elapsed_seconds': 600-(deadline-time.monotonic()), 'runs': results,
               'adapter_sha256': adapter.file_sha(adapter.__file__), 'test_sha256': adapter.file_sha(__file__)}
    adapter.write_json(output/'result.json', summary)
    return summary


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--native', action='store_true', required=True)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(native_qualification(args.root.resolve(), args.output.resolve()), indent=2))
