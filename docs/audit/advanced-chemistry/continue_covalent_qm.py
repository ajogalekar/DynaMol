"""Durable, bounded DF-geometry -> exact-HF ESP -> constrained native RESP.

Numerical completion produces a research candidate, never an accepted physical
force field. A failed or already-started stage is retained, never overwritten.
"""
from __future__ import annotations
import argparse
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'docs/audit/advanced-chemistry'
RUNTIME = ROOT/'.tools/qm-parallel/bin/python'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_geometry(cap, request, result, arrays):
    if not result.get('accepted') or not result.get('optimization', {}).get('converged'):
        raise ValueError('Actual converged quantum optimization is required')
    if request['method'] != 'RHF' or request['basis'].lower() != '6-31g*' or request.get('ecp') or not request.get('density_fit'):
        raise ValueError('Expected explicitly labeled all-electron DF-RHF/6-31G* parent')
    source = json.loads((cap/'capped-adduct.json').read_text())
    molecule = Chem.SDMolSupplier(str(cap/'capped-adduct.sdf'), removeHs=False)[0]
    expected_ids = [f'{a["index"]:03d}:{a["name"]}' for a in source['atom_map']]
    expected_elements = [a.GetSymbol() for a in molecule.GetAtoms()]
    if request['atom_ids'] != expected_ids or result['atom_ids'] != expected_ids or request['elements'] != expected_elements or result['elements'] != expected_elements:
        raise ValueError('Quantum atom identities/order differ from the capped product')
    if request['charge'] != source['formal_charge'] or request['spin'] != 0:
        raise ValueError('Quantum electronic state differs from the capped product')
    xyz = np.asarray(arrays['coords_angstrom'])
    bohr = result['units']['bohr_to_angstrom']
    if xyz.shape != (len(expected_ids), 3) or not np.isfinite(xyz).all() or not np.allclose(xyz, arrays['coords_bohr']*bohr, atol=1e-10, rtol=0):
        raise ValueError('Quantum coordinates have invalid dimensions, units or values')
    cap_ids = [expected_ids[a['index']] for a in source['atom_map'] if a['role'] == 'temporary-cap']
    requested_frozen = request['optimization']['freeze_atom_ids']
    if set(cap_ids) != set(requested_frozen):
        raise ValueError('Exactly the declared temporary cap heavy atoms must be fixed')
    original_xyz = np.asarray(request.get('coords_angstrom', np.asarray(request.get('coords_bohr', []))*bohr))
    frozen = [expected_ids.index(x) for x in requested_frozen]
    frozen_motion = float(np.max(np.linalg.norm(xyz[frozen]-original_xyz[frozen], axis=1)))
    if frozen_motion > 1e-7:
        raise ValueError('Frozen cap heavy atoms moved beyond quantum constraint tolerance')
    for i, coordinate in enumerate(xyz):
        molecule.GetConformer().SetAtomPosition(i, coordinate)
    Chem.AssignStereochemistryFrom3D(molecule, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    names = {a['name']: a['index'] for a in source['atom_map']}
    stereo = []
    for name, expected in source['stereochemistry'].items():
        atom = molecule.GetAtomWithIdx(names[name])
        actual = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
        stereo.append({'name': name, 'expected': expected, 'observed': actual})
        if actual != expected:
            raise ValueError(f'Quantum optimization changed product stereochemistry at {name}')
    table = Chem.GetPeriodicTable()
    distances = []
    for bond in molecule.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        distance = float(np.linalg.norm(xyz[a.GetIdx()]-xyz[b.GetIdx()]))
        ratio = distance/(table.GetRcovalent(a.GetAtomicNum())+table.GetRcovalent(b.GetAtomicNum()))
        distances.append({'atoms': [expected_ids[a.GetIdx()], expected_ids[b.GetIdx()]], 'distance_angstrom': distance, 'covalent_radii_ratio': ratio})
        if not .7 < ratio < 1.35:
            raise ValueError('Implausible bonded distance in optimized product; ESP fit is not started')
    return molecule, {'scope': 'Identity, frozen-cap, product stereo and gross bonded-geometry checks, not chemical accuracy validation.',
        'cap_heavy_ids': cap_ids, 'maximum_cap_displacement_angstrom': frozen_motion,
        'stereochemistry': stereo, 'bonded_geometry': distances, 'passed': True}


def active_qm_processes():
    rows = subprocess.check_output(['ps', '-axo', 'pid=,command='], text=True).splitlines()
    return [{'pid': int(row.strip().split(None, 1)[0]), 'command': row.strip().split(None, 1)[1]}
            for row in rows if ('.tools/qm/bin/python' in row or '.tools/qm-parallel/bin/python' in row)
            and 'continue_covalent_qm.py' not in row and not row.strip().split(None, 1)[1].startswith(('ps ', '/bin/', 'bash '))]


def run(cap, optimizer, reference, workspace, optimizer_wait_seconds):
    for path in (cap, optimizer, reference):
        if not path.is_dir():
            raise ValueError('Expected existing source job directory: '+str(path))
    # A workspace lock alone does not prevent a second output directory from
    # repeating the same expensive parent calculation. Also guard by parent.
    lock_directory = ROOT/'build/advanced-chemistry/covalent/continuation-locks'
    lock_directory.mkdir(parents=True, exist_ok=True)
    parent_key = hashlib.sha256(str(optimizer.resolve()).encode()).hexdigest()
    parent_lock = (lock_directory/(parent_key+'.lock')).open('a')
    fcntl.flock(parent_lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    for row in subprocess.check_output(['ps', '-axo', 'pid=,command='], text=True).splitlines():
        fields = row.strip().split(None, 1)
        if len(fields) != 2 or int(fields[0]) == os.getpid():
            continue
        command = fields[1]
        executable = command.split(None, 1)[0]
        if ('python' in Path(executable).name and 'continue_covalent_qm.py' in command
                and str(optimizer.resolve()) in command):
            raise ValueError('A continuation for this exact parent is already running: PID '+fields[0])
    workspace.mkdir(parents=True, exist_ok=True)
    lock = (workspace/'controller.lock').open('a')
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (workspace/'result.json').exists():
        return json.loads((workspace/'result.json').read_text())
    prior = json.loads((workspace/'progress.json').read_text()) if (workspace/'progress.json').exists() else {}
    started = prior.get('started_unix', time.time())
    state = {'status': 'running', 'stage': 'waiting_for_optimization', 'controller_pid': os.getpid(),
        'started_unix': started, 'full_preparation_ready': False,
        'scope': 'Research candidate only; physical force-field acceptance requires full-adduct torsion and complete complex evidence.',
        'script_sha256': sha(__file__), 'cap_graph_sha256': sha(cap/'capped-adduct.json'),
        'cap_sdf_sha256': sha(cap/'capped-adduct.sdf'), 'optimizer_input_sha256': sha(optimizer/'input.json'),
        'optimizer': str(optimizer), 'exact_reference': str(reference), 'stages': {}}
    if prior:
        state['resumed_unix'] = time.time()
        state['previous_controller_pid'] = prior.get('controller_pid')
        state['previous_controller_script_sha256'] = prior.get('script_sha256')
    def update(**values):
        state.update(values)
        state['elapsed_seconds'] = time.time()-started
        write(workspace/'progress.json', state)
    try:
        deadline = prior.get('wait_deadline_unix', started+optimizer_wait_seconds)
        update(wait_deadline_unix=deadline)
        while not (optimizer/'result.json').exists():
            if time.time() > deadline:
                raise TimeoutError('Waiting for the existing optimizer exceeded its bounded continuation deadline')
            time.sleep(15)
        update(stage='checking_parent_geometry')
        if sha(cap/'capped-adduct.json') != state['cap_graph_sha256'] or sha(cap/'capped-adduct.sdf') != state['cap_sdf_sha256']:
            raise ValueError('Capped product source changed while waiting for optimization')
        request = json.loads((optimizer/'input.json').read_text())
        result = json.loads((optimizer/'result.json').read_text())
        if sha(optimizer/'input.json') != state['optimizer_input_sha256'] or result.get('input_sha256') != state['optimizer_input_sha256']:
            raise ValueError('Optimizer input changed or result does not bind to its original input')
        if not result.get('accepted'):
            raise ValueError('Existing optimization did not succeed: '+str(result.get('error')))
        if sha(optimizer/'arrays.npz') != result['arrays_sha256']:
            raise ValueError('Accepted optimizer coordinate archive hash differs')
        arrays = np.load(optimizer/'arrays.npz')
        optimized, geometry = check_geometry(cap, request, result, arrays)
        write(workspace/'geometry.json', geometry)
        (workspace/'optimized-adduct.sdf').write_text(Chem.MolToMolBlock(optimized)+'\n$$$$\n')
        if not (reference/'result.json').exists():
            raise ValueError('Exact initial-geometry reference is not complete; retained continuation needs explicit resume after review')
        comparison = json.loads((reference/'result.json').read_text())
        if not comparison.get('accepted'):
            raise ValueError('Exact initial-geometry reference did not converge; approximation needs review')
        state['stages']['parent'] = {'result_sha256': sha(optimizer/'result.json'), 'arrays_sha256': result['arrays_sha256'],
            'exact_reference_result_sha256': sha(reference/'result.json'), 'exact_reference_comparison': comparison,
            'interpretation': 'DF geometry approximation remains explicit; a converged reference and reported differences do not establish physical accuracy.'}
        surface_path = BASE/'vendor/resp-surface/vdw_surface.py'
        surface = load_module(surface_path, 'dynamol_resp_surface')
        shells = [surface.vdw_surface(arrays['coords_angstrom'], [x.upper() for x in request['elements']], scale, 1.0, {})[0]
                  for scale in (1.4, 1.6, 1.8, 2.0)]
        grid = np.concatenate(shells)/result['units']['bohr_to_angstrom']
        esp_request = {k: v for k,v in request.items() if k not in ('coords_angstrom', 'coords_bohr', 'optimization', 'initial_checkpoint', 'auxbasis')}
        esp_request.update(coords_bohr=arrays['coords_bohr'].tolist(), density_fit=False,
            operations=['esp'], esp_points_bohr=grid.tolist(), esp_batch_size=64,
            max_wall_seconds=14400, max_memory_mb=8000, threads=2,
            initial_checkpoint={'result_path': str(optimizer/'result.json'), 'result_sha256': sha(optimizer/'result.json'),
                'checkpoint_path': str(optimizer/'scf.chk'), 'checkpoint_sha256': sha(optimizer/'scf.chk')})
        write(workspace/'exact-esp-input.json', esp_request)
        espdir = workspace/'exact-esp'
        if espdir.exists() and not (espdir/'result.json').exists():
            raise ValueError('An earlier exact ESP directory exists without a terminal result; never duplicate or overwrite its process')
        if not espdir.exists():
            update(stage='waiting_for_qm_resource_slot')
            resource_deadline = time.time()+14400
            launch_lock = (ROOT/'build/advanced-chemistry/QM-LAUNCH.lock').open('a')
            while True:
                fcntl.flock(launch_lock.fileno(), fcntl.LOCK_EX)
                processes = active_qm_processes()
                free_disk_gib = shutil.disk_usage(workspace).free/(1024**3)
                if len(processes) < 3 and free_disk_gib >= 15:
                    break
                fcntl.flock(launch_lock.fileno(), fcntl.LOCK_UN)
                if time.time() > resource_deadline:
                    raise TimeoutError('No quantum process slot with at least 15 GiB free disk within four hours')
                update(active_qm_processes=processes, free_disk_gib=free_disk_gib,
                       minimum_free_disk_gib=15, resource_hold='Process count or disk reserve')
                time.sleep(15)
            state['stages']['exact_esp_launch'] = {'input_sha256': sha(workspace/'exact-esp-input.json'),
                'worker_sha256': sha(ROOT/'backend/qm_worker.py'), 'runtime': str(RUNTIME),
                'runtime_executable_sha256': sha(RUNTIME), 'active_qm_at_launch': processes,
                'free_disk_gib_at_launch': free_disk_gib, 'minimum_free_disk_gib': 15,
                'surface_source_sha256': sha(surface_path), 'native_wall_seconds': 14400, 'independent_timeout_seconds': 14520}
            with (workspace/'exact-esp-controller.log').open('w') as log:
                child = subprocess.Popen([str(RUNTIME), str(ROOT/'backend/qm_worker.py'), str(workspace/'exact-esp-input.json'), str(espdir)],
                    cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                update(stage='exact_hf_esp', child_pid=child.pid, stage_deadline_unix=time.time()+14520)
                fcntl.flock(launch_lock.fileno(), fcntl.LOCK_UN)
                try:
                    code = child.wait(timeout=14520)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGTERM)
                    try:
                        child.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid, signal.SIGKILL)
                        child.wait()
                    raise TimeoutError('Exact HF ESP exceeded independent four-hour deadline')
                if code:
                    raise RuntimeError('Exact HF ESP worker failed; its output is retained')
        esp = json.loads((espdir/'result.json').read_text())
        if not esp.get('accepted') or esp['input_sha256'] != sha(workspace/'exact-esp-input.json') or sha(espdir/'arrays.npz') != esp['arrays_sha256']:
            raise ValueError('Exact ESP did not pass numerical and provenance checks')
        esp_arrays = np.load(espdir/'arrays.npz')
        if esp['atom_ids'] != request['atom_ids'] or not np.array_equal(esp_arrays['coords_bohr'], arrays['coords_bohr']):
            raise ValueError('Exact ESP changed product atom identity or geometry')
        state['stages']['exact_esp_result'] = {'path': str(espdir), 'result_sha256': sha(espdir/'result.json'),
            'arrays_sha256': esp['arrays_sha256'], 'energy_hartree': esp['energy_hartree'], 'elapsed_seconds': esp['elapsed_seconds']}
        update(stage='constrained_native_resp', child_pid=None)
        nuclei, points, potentials = esp_arrays['coords_bohr'], esp_arrays['esp_points_bohr'], esp_arrays['esp_hartree_per_e']
        lines = [f'{len(nuclei):5d}{len(points):5d}{0:5d}']
        lines += [' '*17+''.join(f'{x:16.7E}' for x in row) for row in nuclei]
        lines += [' '+''.join(f'{x:16.7E}' for x in [potential,*point]) for potential,point in zip(potentials,points)]
        espfile = workspace/'hf-resp.esp'
        espfile.write_text('\n'.join(lines)+'\n')
        fitdir = workspace/'fitting'
        if fitdir.exists():
            raise ValueError('A prior RESP fitting directory exists; retained without duplicate fitting')
        fitdir.mkdir()
        shutil.copyfile(workspace/'optimized-adduct.sdf', fitdir/'capped-adduct.sdf')
        shutil.copyfile(cap/'capped-adduct.json', fitdir/'capped-adduct.json')
        fitter = load_module(BASE/'covalent_resp_spike.py', 'dynamol_covalent_resp')
        fit = fitter.run(fitdir, espfile)
        fitted = np.array(fit['charges'])
        # Bound memory by processing the grid in chunks for fit diagnostics.
        squared_errors, squared_reference, max_error = 0.0, 0.0, 0.0
        for first in range(0, len(points), 256):
            observed = potentials[first:first+256]
            predicted = (fitted/np.linalg.norm(points[first:first+256,None,:]-nuclei[None,:,:], axis=2)).sum(axis=1)
            error = predicted-observed
            squared_errors += float(error@error)
            squared_reference += float(observed@observed)
            max_error = max(max_error, float(np.max(np.abs(error))))
        state['stages']['resp'] = {'fitter_sha256': sha(BASE/'covalent_resp_spike.py'), 'fit': fit,
            'surface_shells': [1.4,1.6,1.8,2.0], 'surface_density_per_angstrom2': 1.0, 'point_count': len(points),
            'training_rmse_hartree_per_e': float(np.sqrt(squared_errors/len(points))),
            'training_relative_rmse': float(np.sqrt(squared_errors/squared_reference)), 'training_max_error_hartree_per_e': max_error,
            'held_out_accuracy': 'Not assessed; require parent-adduct profiles and separate conformer evidence.'}
        update(status='research_candidate', stage='awaiting_full_adduct_torsion_evidence', child_pid=None)
    except Exception as error:
        update(status='failed', error_type=type(error).__name__, error=str(error), child_pid=None)
    write(workspace/'result.json', state)
    return state


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cap', type=Path)
    parser.add_argument('optimizer', type=Path)
    parser.add_argument('reference', type=Path)
    parser.add_argument('workspace', type=Path)
    parser.add_argument('--optimizer-wait-seconds', type=int, default=21600)
    args = parser.parse_args()
    if not 1 <= args.optimizer_wait_seconds <= 21720:
        parser.error('Optimizer wait must be bounded to at most 21720 seconds')
    final = run(args.cap.resolve(), args.optimizer.resolve(), args.reference.resolve(), args.workspace.resolve(), args.optimizer_wait_seconds)
    print(json.dumps(final, indent=2))
    raise SystemExit(1 if final['status'] == 'failed' else 0)
