"""Bounded small-molecule implementation checks, never biological model data."""
from __future__ import annotations

import ast
import fcntl
import json
import os
from pathlib import Path
import signal
import sys
import time
import tomllib

import numpy as np
from pyscf import dft, gto, lib
from pyscf.data.nist import BOHR
from dftd3 import pyscf as d3pyscf
from dftd3.interface import DispersionModel, RationalDampingParam

import covalent_reference_worker as worker

ROOT = Path(__file__).resolve().parents[3]
SOURCES = worker.HERE/'reference-provider-sources'
OUT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT/'build/advanced-chemistry/covalent/reference-provider/qualification-v2'


def independent_basis(element):
    """Small G94 reader, independently implemented from the PySCF parser."""
    text = worker.BASIS_FILE.read_text()
    blocks = text.split('****')
    for block in blocks:
        lines = [x.strip() for x in block.splitlines() if x.strip() and not x.lstrip().startswith('!')]
        if not lines or lines[0].split() != [element, '0']:
            continue
        shells = []
        i = 1
        while i < len(lines):
            angular, count, scale = lines[i].split()
            assert scale == '1.00', 'Reader intentionally rejects nonunit scale'
            assert angular in 'SPDFGH' and len(angular) == 1
            n = int(count)
            primitives = [[float(x.replace('D', 'E')) for x in row.split()] for row in lines[i+1:i+n+1]]
            assert len(primitives) == n and all(len(p) == 2 for p in primitives)
            shells.append(['SPDFGH'.index(angular), *primitives])
            i += n + 1
        return shells
    raise ValueError('Element missing from exact source')


def d3_fixture():
    tree = ast.parse((SOURCES/'simple-dftd3-v1.6.0-test_pyscf.py').read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'test_energy_b3lyp_d3')
    call = next(n for n in ast.walk(function) if isinstance(n, ast.Call) and any(k.arg == 'atom' for k in n.keywords))
    coordinates = ast.literal_eval(next(k.value for k in call.keywords if k.arg == 'atom'))
    mol = gto.M(atom=coordinates, basis='sto-3g', verbose=0)
    return mol


def correction_checks():
    params = tomllib.loads((SOURCES/'simple-dftd3-v1.6.0-parameters.toml').read_text())
    source = params['parameter']['b3lyp']['d3']['bj']
    for name in ('a1', 'a2', 's8'):
        assert source[name] == worker.D3_PARAMETERS[name]
    mol = d3_fixture()
    xyz = mol.atom_coords()
    model = DispersionModel(mol.atom_charges(), xyz)
    explicit = worker.correction(mol.atom_charges(), xyz)
    named = model.get_dispersion(RationalDampingParam(method='b3lyp', atm=False), grad=True)
    adapter_energy, adapter_gradient, _ = d3pyscf.DFTD3Dispersion(mol, xc='B3LYPG', version='d3bj', atm=False).kernel()
    expected = -.034889411100056264
    assert abs(float(explicit['energy']) - expected) < 1e-10
    assert abs(float(named['energy'] - explicit['energy'])) < 1e-14
    assert np.max(np.abs(named['gradient']-explicit['gradient'])) < 1e-14
    assert abs(float(adapter_energy-explicit['energy'])) < 1e-14
    assert np.max(np.abs(adapter_gradient-explicit['gradient'])) < 1e-14
    h = 1e-4
    numeric = np.zeros_like(xyz)
    for atom in range(len(xyz)):
        for axis in range(3):
            plus, minus = xyz.copy(), xyz.copy()
            plus[atom, axis] += h
            minus[atom, axis] -= h
            numeric[atom, axis] = (float(worker.correction(mol.atom_charges(), plus, False)['energy'])-
                                   float(worker.correction(mol.atom_charges(), minus, False)['energy']))/(2*h)
    error = float(np.max(np.abs(numeric-explicit['gradient'])))
    assert error < 1e-8
    return {'source_fixture': 'upstream simple-dftd3 v1.6.0 test_energy_b3lyp_d3',
        'energy_hartree': float(explicit['energy']), 'published_expected_energy_hartree': expected,
        'all_54_gradient_components_finite_difference_max_error_hartree_per_bohr': error,
        'finite_difference_step_bohr': h, 'named_vs_explicit_vs_adapter_match': True,
        'independence_limit': 'Named/direct/adapter paths share s-dftd3; published fixture and numerical differentiation add independent checks, not a second dispersion implementation.'}


FIXTURES = {
    'water': [('O', [0., 0., 0.]), ('H', [.96, .05, 0.]), ('H', [-.24, .94, .07])],
    'ammonia': [('N', [0., 0., .12]), ('H', [.94, 0., -.2]), ('H', [-.47, .814, -.2]), ('H', [-.47, -.814, -.18])],
    'hydrogen-sulfide': [('S', [0., 0., 0.]), ('H', [1.33, .08, 0.]), ('H', [-.12, 1.34, .1])],
    'fluoromethane': [('C', [0., 0., 0.]), ('F', [0., 0., 1.38]), ('H', [1.03, 0., -.37]), ('H', [-.52, .90, -.37]), ('H', [-.52, -.90, -.37])],
}


def independent_electronic(elements, xyz_bohr, xc='B3LYPG'):
    basis = {e: independent_basis(e) for e in sorted(set(elements))}
    mol = gto.M(atom=list(zip(elements, xyz_bohr)), unit='Bohr', basis=basis, cart=False, spin=0, charge=0,
                max_memory=2048, verbose=0)
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = 4
    mf.conv_tol = 1e-11
    mf.conv_tol_grad = 1e-7
    mf.max_cycle = 150
    energy = mf.kernel()
    assert mf.converged
    grad = mf.nuc_grad_method()
    grad.grid_response = True
    return mol, mf, float(energy), grad.kernel()


def small_molecule_checks(name, fixture):
    request = {'schema_version': 1, 'atom_ids': [f'{i}:{e}' for i, (e, _) in enumerate(fixture)],
        'elements': [e for e, _ in fixture], 'coords_angstrom': [xyz for _, xyz in fixture],
        'charge': 0, 'spin': 0, 'method': 'B3LYP-D3BJ/Psi4-DZVP', 'threads': 2,
        'max_memory_mb': 2048, 'max_wall_seconds': 600, 'operations': ['energy', 'gradient'],
        'provenance': {'qualification_fixture': name, 'biological_model_data': False}}
    directory = OUT/name
    directory.mkdir()
    worker.write(directory/'input.json', request)
    result = worker.run(request, directory)
    arrays = np.load(directory/'arrays.npz')
    xyz = arrays['coords_bohr']
    mol, mf, electronic, gradient = independent_electronic(request['elements'], xyz)
    corr_energy, corr_gradient, _ = d3pyscf.DFTD3Dispersion(mol, xc='B3LYPG', version='d3bj', atm=False).kernel()
    energy_error = abs(result['energy_hartree'] - electronic - float(corr_energy))
    gradient_error = float(np.max(np.abs(arrays['gradient_hartree_per_bohr']-gradient-corr_gradient)))
    assert energy_error < 1e-9
    assert gradient_error < 1e-8
    # All water coordinates and one displacement in every remaining element fixture.
    probes = [(i, a) for i in range(len(xyz)) for a in range(3)] if name == 'water' else [(0, 2), (1, 0)]
    fd_errors = []
    h = 2e-4
    for atom, axis in probes:
        values = []
        for sign in (1, -1):
            coordinates = xyz.copy()
            coordinates[atom, axis] += sign*h
            displaced = mol.copy().set_geom_(coordinates, unit='Bohr', inplace=False)
            model = dft.RKS(displaced)
            model.xc = worker.XC
            model.grids.level = 4
            model.conv_tol = 1e-12
            model.conv_tol_grad = 1e-8
            value = model.kernel(dm0=mf.make_rdm1())
            assert model.converged
            values.append(float(value)+float(worker.correction(displaced.atom_charges(), coordinates, False)['energy']))
        fd_errors.append(abs((values[0]-values[1])/(2*h)-arrays['gradient_hartree_per_bohr'][atom, axis]))
    assert max(fd_errors) < 2e-6, (name, fd_errors)
    return {'fixture': name, 'atoms': len(fixture), 'energy_hartree': result['energy_hartree'],
        'independent_construction_energy_error_hartree': energy_error,
        'independent_construction_gradient_max_error_hartree_per_bohr': gradient_error,
        'finite_difference_components': len(probes), 'finite_difference_step_bohr': h,
        'total_gradient_finite_difference_max_error_hartree_per_bohr': max(fd_errors),
        'result_path': str(directory/'result.json'), 'result_sha256': worker.sha(directory/'result.json'),
        'independence_limit': 'Separately built basis/RKS plus upstream PySCF D3 adapter; same PySCF electronic engine and underlying D3 library.'}


def xc_checks():
    source = ast.parse((SOURCES/'psi4-v1.9.1-libxc_functionals.py').read_text())
    definition = next(n for n in ast.walk(source) if isinstance(n, ast.Dict)
        and any(isinstance(v, ast.Constant) and v.value == 'B3LYP' for v in n.values))
    value = ast.literal_eval(definition)
    assert value['xc_functionals'] == {'HYB_GGA_XC_B3LYP': {}}
    assert dft.libxc.parse_xc('B3LYPG')[1] == ((402, 1),)
    xyz = np.asarray([v for _, v in FIXTURES['water']])/BOHR
    elements = [e for e, _ in FIXTURES['water']]
    _, _, e1, g1 = independent_electronic(elements, xyz)
    _, _, e2, g2 = independent_electronic(elements, xyz, '.2*HF + .08*SLATER + .72*B88, .81*LYP + .19*VWNRPA')
    _, _, e5, _ = independent_electronic(elements, xyz, 'B3LYP5')
    assert abs(e1-e2) < 1e-9
    assert np.max(np.abs(g1-g2)) < 1e-8
    assert abs(e1-e5) > 1e-5
    return {'psi4_source_definition': value, 'libxc_id': 402,
        'equivalent_explicit_mixture': '.2 HF + .08 Slater + .72 B88 exchange; .81 LYP + .19 VWN-RPA correlation',
        'mixture_energy_error_hartree': abs(e1-e2), 'mixture_gradient_max_error_hartree_per_bohr': float(np.max(np.abs(g1-g2))),
        'different_VWN5_energy_delta_hartree': e5-e1}


def input_rejection_checks():
    base = {'schema_version': 1, 'atom_ids': ['O', 'H1', 'H2'], 'elements': ['O', 'H', 'H'],
        'coords_angstrom': [x for _, x in FIXTURES['water']], 'charge': 0, 'spin': 0,
        'method': 'B3LYP-D3BJ/Psi4-DZVP'}
    variants = [dict(base, atom_ids=['O', 'H', 'H']), dict(base, spin=2), dict(base, spin=False), dict(base, charge=True),
        dict(base, elements=['Zn', 'H', 'H']), dict(base, threads=3), dict(base, max_memory_mb=16000),
        dict(base, coords_bohr=base['coords_angstrom']), dict(base, density_fit=True),
        dict(base, method='B3LYP'), dict(base, coords_angstrom=[[0, 0, 0]]*3)]
    for request in variants:
        try:
            worker.validate(request)
        except ValueError:
            continue
        raise AssertionError('Invalid/ambiguous input accepted')
    return {'rejected_cases': len(variants), 'passed': True}


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    started = time.time()
    report = {'status': 'running', 'started_unix': started, 'pid': os.getpid(), 'tests': {},
        'scope': 'Implementation and numerical derivative qualification only; not independent physical validation of covalent parameters.',
        'worker_sha256': worker.sha(worker.__file__), 'qualification_script_sha256': worker.sha(__file__)}
    def update(stage):
        report['stage'] = stage
        worker.write(OUT/'progress.json', report)
    def expire(signum, frame):
        raise TimeoutError('Small-reference qualification exceeded bounded 1800-second wall time')
    signal.signal(signal.SIGALRM, expire)
    signal.alarm(1800)
    try:
        # Hold the cooperative launch lock for the entire small test process: the
        # older resident controller does not recognize this new runtime name.
        with (ROOT/'build/advanced-chemistry/QM-LAUNCH.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            lib.num_threads(2)
            update('source_and_input_checks')
            report['tests']['input_rejections'] = input_rejection_checks()
            basis = worker.basis_for(worker.ELEMENTS)
            for element, expanded in basis.items():
                assert expanded == independent_basis(element)
            report['tests']['basis'] = {'independent_parser_exact_match': True,
                'raw_basis_sha256': worker.BASIS_SHA, 'expanded_basis_sha256': worker.canonical_hash(basis),
                'elements': sorted(basis), 'shell_counts': {e: len(b) for e, b in basis.items()}, 'spherical': True}
            update('dispersion_checks')
            report['tests']['d3'] = correction_checks()
            update('xc_convention_checks')
            report['tests']['xc'] = xc_checks()
            report['tests']['small_molecules'] = []
            for name, fixture in FIXTURES.items():
                update('small_molecule_'+name)
                report['tests']['small_molecules'].append(small_molecule_checks(name, fixture))
        report['status'] = 'passed_numerical_qualification'
    except Exception as error:
        report['status'] = 'failed'
        report['error_type'] = type(error).__name__
        report['error'] = str(error)
        raise
    finally:
        report['elapsed_seconds'] = time.time()-started
        worker.write(OUT/'result.json', report)
        signal.alarm(0)


if __name__ == '__main__':
    main()
