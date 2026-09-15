"""Read-only failed-parent review and a new, explicitly unconverged restart input.

No worker, SCF, optimizer, checkpoint import, or source mutation is performed.
Coordinates are recovered from the numbered 12-decimal native log records.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
CAP = ROOT/'build/advanced-chemistry/covalent'
RUN = CAP/'df-optimization-parallel-v1'
JOB = RUN/'optimization'


def sha(path):
    raw = Path(path).read_bytes()
    if not raw:
        raise ValueError('Empty source read: '+str(path))
    return hashlib.sha256(raw).hexdigest()


def write(name, value):
    path = OUT/name
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def geometry_checks(xyz, request, source, template):
    molecule = Chem.Mol(template)
    ids = request['atom_ids']
    expected = [f'{atom["index"]:03d}:{atom["name"]}' for atom in source['atom_map']]
    assert ids == expected
    assert request['elements'] == [atom.GetSymbol() for atom in molecule.GetAtoms()]
    assert request['charge'] == source['formal_charge'] == 0 and request['spin'] == 0
    assert xyz.shape == (94, 3) and np.isfinite(xyz).all()
    frozen = [ids.index(name) for name in request['optimization']['freeze_atom_ids']]
    cap_indices = [atom['index'] for atom in source['atom_map'] if atom['role'] == 'temporary-cap']
    assert set(frozen) == set(cap_indices) and len(frozen) == 5
    original = np.asarray(request['coords_angstrom'])
    max_cap = float(np.linalg.norm(xyz[frozen]-original[frozen], axis=1).max())
    for index, position in enumerate(xyz):
        molecule.GetConformer().SetAtomPosition(index, position)
    Chem.AssignStereochemistryFrom3D(molecule, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    names = {row['name']: row['index'] for row in source['atom_map']}
    stereo = []
    for name, expected in source['stereochemistry'].items():
        atom = molecule.GetAtomWithIdx(names[name])
        observed = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
        stereo.append({'name': name, 'expected': expected, 'observed': observed,
                       'passed': observed == expected})
    table = Chem.GetPeriodicTable()
    distances = []
    for bond in molecule.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        distance = float(np.linalg.norm(xyz[a.GetIdx()]-xyz[b.GetIdx()]))
        ratio = distance/(table.GetRcovalent(a.GetAtomicNum())+table.GetRcovalent(b.GetAtomicNum()))
        distances.append({'atom_ids': [ids[a.GetIdx()], ids[b.GetIdx()]],
            'order': float(bond.GetBondTypeAsDouble()), 'distance_angstrom': distance,
            'covalent_radii_ratio': ratio, 'passed': .7 < ratio < 1.35})
    return {'geometry_identity_checks_passed': max_cap <= 1e-7 and
                all(row['passed'] for row in stereo+distances),
            'maximum_frozen_displacement_angstrom': max_cap,
            'frozen_tolerance_angstrom': 1e-7, 'stereochemistry': stereo,
            'bond_geometry': distances, 'physical_acceptance': False,
            'meaning': 'Declared reacted graph/order, frozen anchors, two expected CIP labels, and gross bond geometry only; not optimization convergence or chemical-model accuracy.'}


def main():
    source_files = [JOB/'input.json', JOB/'resolved-input.json', JOB/'result.json',
        JOB/'native.log', JOB/'error-traceback.txt', JOB/'constraints.txt',
        JOB/'initial-checkpoint/result.json', RUN/'controller.log', RUN/'progress.json',
        CAP/'capped-adduct.json', CAP/'capped-adduct.sdf', ROOT/'backend/qm_worker.py',
        ROOT/'.tools/qm-parallel/lib/python3.12/site-packages/pyscf/geomopt/geometric_solver.py']
    pins = {str(path): sha(path) for path in source_files}
    request = json.loads((JOB/'input.json').read_text())
    outcome = json.loads((JOB/'result.json').read_text())
    source = json.loads((CAP/'capped-adduct.json').read_text())
    initial = json.loads((JOB/'initial-checkpoint/result.json').read_text())
    template = Chem.SDMolSupplier(str(CAP/'capped-adduct.sdf'), removeHs=False)[0]
    assert template is not None and template.GetNumAtoms() == 94
    assert outcome['accepted'] is False and outcome['optimization']['converged'] is None
    assert outcome['input_sha256'] == pins[str(JOB/'input.json')]
    assert outcome['worker_sha256'] == pins[str(ROOT/'backend/qm_worker.py')]
    assert request['initial_checkpoint']['result_sha256'] == pins[str(JOB/'initial-checkpoint/result.json')]
    bohr = initial['units']['bohr_to_angstrom']
    text = (JOB/'native.log').read_text()
    pieces = re.split(r'Geometry optimization cycle (\d+)\n', text)
    number = r'([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)'
    pattern = re.compile(r'^\s*(\d+)\s+([A-Za-z]+)\s+'+r'\s+'.join([number]*3)+
        r'\s+AA\s+'+r'\s+'.join([number]*3)+r'\s+Bohr\s*$', re.M)
    frames = {}
    frame_checks = []
    for index in range(1, len(pieces), 2):
        cycle, block = int(pieces[index]), pieces[index+1]
        rows = pattern.findall(block)
        assert len(rows) == 94
        assert [int(row[0]) for row in rows] == list(range(1, 95))
        assert [row[1] for row in rows] == request['elements']
        xyz = np.array([[float(x) for x in row[2:5]] for row in rows])
        xyz_bohr = np.array([[float(x) for x in row[5:8]] for row in rows])
        discrepancy = float(np.max(np.abs(xyz-xyz_bohr*bohr)))
        assert discrepancy < 1e-12
        checks = geometry_checks(xyz, request, source, template)
        frames[cycle] = {'cycle': cycle, 'atom_ids': request['atom_ids'],
            'elements': request['elements'], 'coords_angstrom': xyz.tolist(),
            'coords_bohr': xyz_bohr.tolist(), 'source_log_sha256': pins[str(JOB/'native.log')],
            'native_log_start_line': text[:text.index('Geometry optimization cycle '+str(cycle)+'\n')].count('\n')+1,
            'coordinate_precision': '12 digits after decimal in each separately printed unit; not unrounded checkpoint coordinates',
            'unit_roundtrip_max_component_difference_angstrom': discrepancy,
            'geometry_checks': checks, 'scf_and_gradient_completed': cycle <= 57,
            'optimization_converged': False, 'accepted': False, 'physical_acceptance': False}
        frame_checks.append({'cycle': cycle, 'all_94_numbered_atoms_and_element_order_match': True,
            'unit_roundtrip_max_component_difference_angstrom': discrepancy,
            'maximum_frozen_displacement_angstrom': checks['maximum_frozen_displacement_angstrom'],
            'stereochemistry': checks['stereochemistry'],
            'geometry_identity_checks_passed': checks['geometry_identity_checks_passed']})
    assert sorted(frames) == list(range(1, 59))
    first_difference = float(np.max(np.abs(np.asarray(frames[1]['coords_angstrom'])-request['coords_angstrom'])))
    # Initial adapter unit roundtrip adds a small, directly measured numerical
    # difference; ordered source/input identity is also verified independently.
    assert first_difference < 1e-8
    assert all(row['geometry_identity_checks_passed'] for row in frame_checks)
    steps = outcome['optimization']['steps']
    assert len(steps) == 57 and [step['cycle'] for step in steps] == list(range(1, 58))
    assert all(step['scf_converged'] for step in steps)
    controller = re.sub(r'\x1b\[[0-9;]*m', '', (RUN/'controller.log').read_text())
    last_step_line = [line for line in controller.splitlines() if re.match(r'Step\s+56\s*:', line)][-1]
    assert '-2794.0693928687' in last_step_line and 'Quality = 1.135' in last_step_line
    assert 'Optimization converged' not in controller
    frames[57]['completed_native_observation'] = steps[-1]
    frames[57]['geometric_step'] = 56
    frames[57]['geometric_step_line'] = last_step_line
    frames[58]['reason_not_selected'] = 'Proposed geometry printed before failing HDF5 save_mol check; no completed SCF, gradient, or geomeTRIC step evaluation.'
    write('recovered-cycle-057.json', frames[57])
    write('recovered-cycle-058-unscored.json', frames[58])
    write('all-58-frame-checks.json', frame_checks)

    proposed = deepcopy(request)
    proposed.pop('initial_checkpoint')
    proposed['max_memory_mb'] = 5000
    proposed['max_wall_seconds'] = 21600
    proposed['optimization']['max_steps'] = 100
    proposed['coords_angstrom'] = deepcopy(frames[57]['coords_angstrom'])
    corrections = []
    frozen_ids = request['optimization']['freeze_atom_ids']
    for atom_id in frozen_ids:
        index = request['atom_ids'].index(atom_id)
        recovered = proposed['coords_angstrom'][index]
        original = request['coords_angstrom'][index]
        delta = (np.asarray(original)-recovered).tolist()
        corrections.append({'atom_id': atom_id, 'recovered_angstrom': recovered,
            'restored_original_angstrom': original, 'correction_angstrom': delta,
            'norm_angstrom': float(np.linalg.norm(delta)),
            'interpretation': 'Numerical native constraint residual, not merely 12-decimal text rounding. Explicit restoration to original exact frozen anchor.'})
        proposed['coords_angstrom'][index] = deepcopy(original)
    proposal_checks = geometry_checks(np.asarray(proposed['coords_angstrom']), request, source, template)
    assert proposal_checks['geometry_identity_checks_passed']
    assert proposal_checks['maximum_frozen_displacement_angstrom'] == 0
    free = [i for i, atom_id in enumerate(request['atom_ids']) if atom_id not in frozen_ids]
    assert np.array_equal(np.asarray(proposed['coords_angstrom'])[free], np.asarray(frames[57]['coords_angstrom'])[free])
    changed = [key for key in request if request[key] != proposed.get(key)]
    assert set(changed) == {'coords_angstrom', 'initial_checkpoint', 'max_memory_mb'}
    write('proposed-cold-restart-input.json', proposed)
    write('proposed-geometry-checks.json', proposal_checks)
    checkpoint_stat = (JOB/'scf.chk').stat()
    report = {'status': 'unconverged_restart_geometry_review_complete',
        'accepted_quantum_result': False, 'physical_acceptance': False, 'qm_launched': False,
        'source_files_sha256': pins, 'source_files_modified': False,
        'rdkit_version_for_declared_stereo_checks': rdBase.rdkitVersion,
        'all_native_frames_checked': 58, 'first_frame_vs_source_max_component_difference_angstrom': first_difference,
        'recommended_seed': 'PySCF cycle 57 / geomeTRIC step 56, with five explicitly restored original frozen anchors',
        'last_completed_energy_hartree': steps[-1]['energy_hartree'],
        'last_completed_gradient_norm_hartree_per_bohr': steps[-1]['gradient_norm_hartree_per_bohr'],
        'last_geometric_step_line': last_step_line,
        'remaining_convergence_failures_from_native_log': {
            'energy_change_abs_hartree': {'observed_rounded': 6.040e-5, 'threshold': 1e-6},
            'max_tangent_gradient_hartree_per_bohr': {'observed_rounded': 6.889e-4, 'threshold': 4.5e-4},
            'rms_displacement_angstrom': {'observed_rounded': .03525, 'threshold': .0012},
            'max_displacement_angstrom': {'observed_rounded': .1043, 'threshold': .0018}},
        'temporary_native_trajectory_directory': '/var/folders/cv/_3qd6yhx7lqdggx4cs3hlg5c0000gn/T/tmp72tvhhle',
        'temporary_native_trajectory_directory_present': False,
        'trajectory_recovery_provenance': 'PySCF numbered full-precision New geometry blocks, not a surviving native XYZ file. Pinned adapter source shows mol.set_geom_ preserving ordered atoms before each evaluation.',
        'failed_checkpoint': {'path': str(JOB/'scf.chk'), 'size_bytes': checkpoint_stat.st_size,
            'mtime_ns': checkpoint_stat.st_mtime_ns, 'usable_checkpoint_not_established': True,
            'read_probe': 'Read-only 25-second original-runtime probe timed out without data; import/file-read stage not resolved. Further probes stopped per root instruction.',
            'import_into_new_run': False,
            'reason': 'Failed worker result is unaccepted. No fabricated accepted checkpoint wrapper; fresh cold SCF required.'},
        'initial_checkpoint_provenance': outcome['initial_checkpoint'],
        'frozen_anchor_corrections': corrections,
        'proposal_path': str(OUT/'proposed-cold-restart-input.json'),
        'proposal_sha256': sha(OUT/'proposed-cold-restart-input.json'),
        'input_difference_summary': {'only_changed_top_level_keys': changed,
            'coords_angstrom': '89 free atoms at cycle-57 logged positions; five anchors exactly restored to original input positions',
            'initial_checkpoint': 'removed', 'max_memory_mb': {'old': request['max_memory_mb'], 'new': proposed['max_memory_mb']},
            'all_other_settings': 'exactly unchanged, including atom IDs/elements, charge/spin, method, basis, auxiliary basis, frozen IDs and all convergence thresholds'},
        'original_threshold_units': {'energy': 'Hartree', 'gradients': 'Hartree/Bohr', 'displacements': 'Angstrom'},
        'new_optimization_requires': 'Independent cold SCF and native optimization convergence at the original thresholds, then final geometry/identity checks and fresh bound result/arrays; downstream RESP and DFT must bind that new accepted parent.',
        'scope': 'Recoverable input geometry and provenance only; does not validate a force field or chemically accurate dynamics.'}
    assert all(sha(path) == expected for path, expected in pins.items())
    write('result.json', report)
    print(json.dumps({key: report[key] for key in ['status','recommended_seed','proposal_path','proposal_sha256']}, indent=2))


if __name__ == '__main__':
    main()
