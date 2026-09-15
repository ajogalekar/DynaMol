"""Isolated research B3LYP-D3(BJ)/DZVP single-point provider.

Numerical completion is not force-field or chemical-state validation. This is
deliberately separate from the production neutral QM worker and its running jobs.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import time

import numpy as np
from pyscf import dft, gto, lib
from pyscf.gto.basis import parse_gaussian
from pyscf.data.nist import BOHR
from dftd3.interface import DispersionModel, RationalDampingParam

HERE = Path(__file__).resolve().parent
BASIS_FILE = HERE / 'reference-provider-sources/psi4-v1.9.1-dzvp.gbs'
BASIS_SHA = '629b5f1d8b7eb52ab43dd84e7c0c5e3e4d0f25cff8ed858edf407788d79a968c'
ELEMENTS = frozenset(('H', 'C', 'N', 'O', 'F', 'S'))
XC = 'B3LYPG'
D3_PARAMETERS = dict(s6=1., s8=1.9889, a1=.3981, a2=4.4211, s9=0.)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def write(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def basis_for(elements):
    if sha(BASIS_FILE) != BASIS_SHA:
        raise ValueError('Pinned Psi4 DZVP basis source changed')
    if not set(elements) <= ELEMENTS:
        raise ValueError('This research provider is qualified only for H/C/N/O/F/S')
    return {element: parse_gaussian.load(str(BASIS_FILE), element, optimize=False)
            for element in sorted(set(elements))}


def validate(request):
    if request.get('schema_version') != 1:
        raise ValueError('Expected schema_version 1')
    ids, elements = request.get('atom_ids'), request.get('elements')
    if (not isinstance(ids, list) or not ids or any(not isinstance(i, str) or not i for i in ids)
            or len(set(ids)) != len(ids) or not isinstance(elements, list) or len(elements) != len(ids)):
        raise ValueError('Unique nonempty atom IDs and one element per atom are required')
    if not set(elements) <= ELEMENTS:
        raise ValueError('Only explicitly qualified elements H/C/N/O/F/S are supported')
    if type(request.get('charge')) is not int or type(request.get('spin')) is not int or request.get('spin') != 0:
        raise ValueError('An explicit integer charge and closed-shell spin 0 are required')
    if request.get('method') != 'B3LYP-D3BJ/Psi4-DZVP':
        raise ValueError('Reference method must be explicitly selected')
    if ('coords_angstrom' in request) == ('coords_bohr' in request):
        raise ValueError('Supply exactly one coordinate unit')
    xyz = np.asarray(request.get('coords_bohr', request.get('coords_angstrom')), dtype=float)
    if 'coords_angstrom' in request:
        xyz = xyz / BOHR
    if xyz.shape != (len(ids), 3) or not np.isfinite(xyz).all():
        raise ValueError('Finite Nx3 coordinates are required')
    if len(ids) > 128:
        raise ValueError('Bounded research provider supports at most 128 atoms')
    if len(xyz) > 1:
        distances = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
        np.fill_diagonal(distances, np.inf)
        if distances.min() < .3:
            raise ValueError('Coincident or implausibly close atomic centers')
    for name, low, high, default in [('threads', 1, 2, 2), ('max_memory_mb', 256, 8000, 2048),
                                    ('max_wall_seconds', 1, 14400, 1200)]:
        value = request.get(name, default)
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f'Invalid bounded resource setting: {name}')
    if request.get('operations', ['energy', 'gradient']) not in (['energy'], ['energy', 'gradient']):
        raise ValueError('Only energy or energy plus gradient are qualified')
    allowed = {'schema_version', 'atom_ids', 'elements', 'charge', 'spin', 'method', 'coords_bohr',
               'coords_angstrom', 'threads', 'max_memory_mb', 'max_wall_seconds', 'operations', 'provenance'}
    if set(request) - allowed:
        raise ValueError('Unknown reference settings are not silently ignored: ' + str(sorted(set(request)-allowed)))
    return xyz


@contextmanager
def execution_guard(request):
    """Cooperate with active research jobs without changing any resident worker."""
    root = HERE.parents[2]
    with (root/'build/advanced-chemistry/QM-LAUNCH.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        active = []
        for row in subprocess.check_output(['ps', '-axo', 'pid=,command='], text=True).splitlines():
            fields = row.strip().split(None, 1)
            if len(fields) != 2 or int(fields[0]) == os.getpid():
                continue
            command = fields[1]
            if 'python' not in Path(command.split()[0]).name:
                continue
            if 'qm_worker.py' in command or 'covalent_reference_worker.py' in command:
                active.append({'pid': int(fields[0]), 'command': command})
        if len(active) >= 3:
            raise ValueError('Research concurrency limit would be exceeded')
        free = shutil.disk_usage(root).free
        reserve = (15 if len(request['atom_ids']) > 16 else 2)*1024**3
        if free < reserve:
            raise ValueError('Insufficient free disk reserve for the bounded reference')
        # Hold for this job: older resident launchers cannot see this new runtime
        # name in their process scanner, but all honor the shared launch lock.
        yield {'existing_qm_processes': active, 'free_disk_bytes_at_launch': free,
               'minimum_free_disk_bytes': reserve, 'launch_lock_held_through_completion': True}


def electronic_model(elements, xyz_bohr, charge, memory, log=None):
    molecule = gto.M(atom=list(zip(elements, xyz_bohr)), unit='Bohr', charge=charge, spin=0,
                     basis=basis_for(elements), cart=False, verbose=4 if log else 0,
                     output=str(log) if log else None, max_memory=memory)
    if molecule.nelectron % 2:
        raise ValueError('Closed-shell reference requires an even electron count')
    mf = dft.RKS(molecule)
    mf.xc = XC
    mf.grids.level = 4
    mf.conv_tol = 1e-11
    mf.conv_tol_grad = 1e-7
    mf.max_cycle = 150
    return molecule, mf


def correction(numbers, xyz_bohr, gradient=True):
    return DispersionModel(np.asarray(numbers), np.asarray(xyz_bohr)).get_dispersion(
        RationalDampingParam(**D3_PARAMETERS), grad=gradient)


def run(request, output):
    xyz = validate(request)
    started = time.time()
    lib.num_threads(request.get('threads', 2))
    write(output/'progress.json', {'status': 'running', 'stage': 'electronic_scf', 'pid': os.getpid(), 'started_unix': started})
    molecule, mf = electronic_model(request['elements'], xyz, request['charge'],
                                     request.get('max_memory_mb', 2048), output/'native.log')
    electronic_energy = float(mf.kernel())
    if not mf.converged or not np.isfinite(electronic_energy):
        raise ValueError('Electronic SCF did not converge to a finite energy')
    want_gradient = 'gradient' in request.get('operations', ['energy', 'gradient'])
    write(output/'progress.json', {'status': 'running', 'stage': 'gradient_and_dispersion', 'pid': os.getpid(), 'started_unix': started})
    disp = correction(molecule.atom_charges(), xyz, want_gradient)
    dispersion_energy = float(disp['energy'])
    arrays = {'coords_bohr': xyz, 'coords_angstrom': xyz*BOHR,
              'atomic_numbers': molecule.atom_charges(),
              'energy_hartree': np.asarray(electronic_energy + dispersion_energy),
              'electronic_energy_hartree': np.asarray(electronic_energy),
              'dispersion_energy_hartree': np.asarray(dispersion_energy)}
    if want_gradient:
        grad = mf.nuc_grad_method()
        grad.grid_response = True
        arrays['electronic_gradient_hartree_per_bohr'] = grad.kernel()
        arrays['dispersion_gradient_hartree_per_bohr'] = np.asarray(disp['gradient'])
        arrays['gradient_hartree_per_bohr'] = arrays['electronic_gradient_hartree_per_bohr'] + arrays['dispersion_gradient_hartree_per_bohr']
    if not all(np.isfinite(v).all() for v in arrays.values()):
        raise ValueError('Nonfinite reference output')
    np.savez(output/'arrays.npz', **arrays)
    result = {'schema_version': 1, 'accepted': True, 'status': 'numerically_complete_research_reference',
        'acceptance_scope': 'Converged finite single-point energy/gradient; not physical force-field acceptance.',
        'atom_ids': request['atom_ids'], 'elements': request['elements'], 'charge': request['charge'], 'spin': 0,
        'energy_hartree': electronic_energy + dispersion_energy,
        'electronic_energy_hartree': electronic_energy, 'dispersion_energy_hartree': dispersion_energy,
        'method': {'name': request['method'], 'electronic': 'conventional RKS', 'density_fitting': False,
            'xc': XC, 'libxc_id': 402, 'libxc_version': dft.libxc.__version__,
            'basis': 'Psi4 v1.9.1 DZVP (DFT orbital), Godbout et al. 1992',
            'basis_source_sha256': BASIS_SHA, 'expanded_basis_sha256': canonical_hash(molecule._basis),
            'expanded_basis': molecule._basis, 'spherical': True, 'ecp': None,
            'dispersion': 'D3(BJ)', 'dispersion_parameters': D3_PARAMETERS, 'atm': False,
            'grid_level': 4, 'grid_response_in_gradient': True, 'grid_pruning': mf.grids.prune.__name__,
            'conv_tol': mf.conv_tol, 'conv_tol_grad': mf.conv_tol_grad, 'scf_cycles': mf.cycles},
        'units': {'coordinates': 'Bohr and Angstrom explicitly named in arrays', 'energy': 'Hartree',
                  'gradient': 'Hartree/Bohr', 'bohr_to_angstrom': BOHR},
        'runtime': {'python': sys.version, 'executable': sys.executable,
            'pyscf': importlib.metadata.version('pyscf'), 'dftd3': importlib.metadata.version('dftd3'),
            'numpy': np.__version__, 'threads': lib.num_threads(), 'max_memory_mb': molecule.max_memory},
        'source': {'worker_sha256': sha(__file__), 'input_sha256': sha(output/'input.json'),
                   'provenance': request.get('provenance', {}),
                   'primary_source_manifest_sha256': sha(HERE/'reference-provider-sources/manifest.json')},
        'arrays_file': 'arrays.npz', 'arrays_sha256': sha(output/'arrays.npz'),
        'started_unix': started, 'elapsed_seconds': time.time()-started}
    write(output/'result.json', result)
    write(output/'progress.json', {'status': result['status'], 'accepted': True, 'elapsed_seconds': result['elapsed_seconds']})
    molecule.stdout.close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    request = json.loads(args.input.read_text())
    write(args.output/'input.json', request)
    def expire(signum, frame):
        raise TimeoutError('Bounded reference wall time expired')
    signal.signal(signal.SIGALRM, expire)
    try:
        validate(request)
        signal.alarm(request.get('max_wall_seconds', 1200))
        with execution_guard(request) as resources:
            write(args.output/'resource-check.json', resources)
            run(request, args.output)
    except Exception as error:
        write(args.output/'result.json', {'schema_version': 1, 'accepted': False, 'status': 'failed',
            'error_type': type(error).__name__, 'error': str(error), 'worker_sha256': sha(__file__),
            'input_sha256': sha(args.output/'input.json')})
        raise
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
