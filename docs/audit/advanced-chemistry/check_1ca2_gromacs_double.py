"""Two native double-precision zero-step evaluations of an unchanged TPR."""
from pathlib import Path
import hashlib
import json
import subprocess
import time
import numpy as np
from read_trr_precision import read_first_trr

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
AUDIT = ROOT / 'build/advanced-chemistry/gromacs-double-v1'
P = BASE / 'prototype-1ca2-periodic-diagnostic'
provenance = json.loads((AUDIT / 'build-provenance.json').read_text())
if provenance['status'] != 'built':
    raise ValueError('The isolated double-precision build must finish first')
exe = Path(provenance['executable'])
tpr = P / 'matched-mesh/diagnostic.tpr'
digest = lambda x: hashlib.sha256(Path(x).read_bytes()).hexdigest()
if digest(tpr) != '5aa9285424526503ae95aebc3c786baf5c318ea86fbbdac3a12fe276b2ac3ac2':
    raise ValueError('Original fixed diagnostic TPR changed')
version = subprocess.run([str(exe), '--version'], capture_output=True, text=True, timeout=60)
(AUDIT / 'native-version.txt').write_text(version.stdout + version.stderr)
if version.returncode or 'Precision:           double' not in version.stdout + version.stderr:
    raise ValueError('Native executable did not identify double precision')
reference = np.load(P / 'openmm-reference.npz')
reference_energy = json.loads((P / 'openmm-reference.json').read_text())['energy_kj_mol']
original = read_first_trr(P / 'matched-mesh/diagnostic.trr')
fields = ['Bond', 'Angle', 'Proper-Dih.', 'Per.-Imp.-Dih.', 'LJ-14', 'Coulomb-14',
          'LJ-(SR)', 'Disper.-corr.', 'Coulomb-(SR)', 'Coul.-recip.', 'Potential']
records = []
arrays = []
started = time.monotonic()
for threads in [2, 1]:
    folder = AUDIT / f'static-thread-{threads}'
    folder.mkdir(exist_ok=False)
    command = [str(exe), 'mdrun', '-s', str(tpr), '-deffnm', 'double',
               '-ntmpi', '1', '-ntomp', str(threads), '-nb', 'cpu', '-notunepme']
    with (folder / 'mdrun.log').open('w') as log:
        run = subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=180)
    if run.returncode:
        raise RuntimeError('Native static evaluation failed; exact log retained')
    energy_command = [str(exe), 'energy', '-f', 'double.edr', '-o', 'energy.xvg', '-xvg', 'none', '-dp']
    with (folder / 'energy.log').open('w') as log:
        run = subprocess.run(energy_command, cwd=folder, input='\n'.join(fields) + '\n0\n',
                             text=True, stdout=log, stderr=subprocess.STDOUT, timeout=60)
    if run.returncode:
        raise RuntimeError('Native energy extraction failed')
    energy = np.loadtxt(folder / 'energy.xvg')
    if energy.shape != (len(fields) + 1,) or not np.isfinite(energy).all():
        raise ValueError('Expected one finite zero-step row with every selected energy component')
    frame = read_first_trr(folder / 'double.trr')
    if frame['stored_precision_bytes'] != 8 or frame['step'] != 0 or frame['time_ps'] != 0:
        raise ValueError('Expected a genuine double-precision zero-step frame')
    if frame['coords_nm'].shape != reference['coords_nm'].shape:
        raise ValueError('Native atom inventory changed')
    delta = frame['forces_kj_mol_nm'] - reference['forces_kj_mol_nm']
    coord_delta = frame['coords_nm'] - original['coords_nm']
    mic_delta = coord_delta - np.rint(coord_delta @ np.linalg.inv(frame['box_nm'])) @ frame['box_nm']
    box_delta = frame['box_nm'] - original['box_nm']
    if np.max(np.abs(mic_delta)) > 1e-6 or np.max(np.abs(box_delta)) > 1e-6:
        raise ValueError('Double runtime changed the physical input coordinates/box beyond stored rounding')
    np.savez_compressed(folder / 'double-frame.npz', coords_nm=frame['coords_nm'], box_nm=frame['box_nm'],
                        forces_kj_mol_nm=frame['forces_kj_mol_nm'])
    result = {'threads': threads, 'command': command, 'energy_command': energy_command,
              'stored_trr_precision_bytes': frame['stored_precision_bytes'],
              'components_kj_mol': dict(zip(fields, map(float, energy[1:]))),
              'energy_difference_from_openmm_reference_kj_mol': float(energy[-1] - reference_energy),
              'maximum_coordinate_difference_from_original_nm': float(np.max(np.abs(coord_delta))),
              'maximum_minimum_image_coordinate_difference_nm': float(np.max(np.abs(mic_delta))),
              'maximum_box_difference_nm': float(np.max(np.abs(box_delta))),
              'maximum_force_vector_difference_kj_mol_nm': float(np.linalg.norm(delta, axis=1).max()),
              'rms_force_vector_difference_kj_mol_nm': float(np.sqrt(np.mean(np.sum(delta * delta, axis=1)))),
              'finite_full_system_forces': bool(np.isfinite(frame['forces_kj_mol_nm']).all()),
              'input_tpr_sha256': digest(tpr),
              'output_hashes': {n: digest(folder / n) for n in ['double.trr', 'double.edr', 'energy.xvg', 'double-frame.npz']}}
    (folder / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    records.append(result)
    arrays.append(frame)
    print(json.dumps({k: result[k] for k in ['threads', 'energy_difference_from_openmm_reference_kj_mol',
                                            'maximum_force_vector_difference_kj_mol_nm', 'rms_force_vector_difference_kj_mol_nm']}, indent=2), flush=True)
force_delta = arrays[1]['forces_kj_mol_nm'] - arrays[0]['forces_kj_mol_nm']
report = {'status': 'completed', 'scope': 'Native double-precision zero-step discriminator only; no MD, changed TPR, parameter edits or arbitrary energy offsets.',
          'elapsed_seconds': time.monotonic() - started, 'source_archive_sha256': provenance['source_archive_sha256'],
          'native_executable': str(exe), 'native_executable_sha256': digest(exe), 'native_compiler_difference': provenance['compiler'],
          'openmm_reference_energy_kj_mol': reference_energy, 'cases': records,
          'thread_energy_difference_1_minus_2_kj_mol': records[1]['components_kj_mol']['Potential'] - records[0]['components_kj_mol']['Potential'],
          'thread_force_rms_vector_difference_kj_mol_nm': float(np.sqrt(np.mean(np.sum(force_delta * force_delta, axis=1)))),
          'thread_coordinate_arrays_exact': bool(np.array_equal(arrays[0]['coords_nm'], arrays[1]['coords_nm'])),
          'original_installed_binary_unchanged': digest(ROOT / '.gromacs/bin.ARM_NEON_ASIMD/gmx') == provenance['original_installed_binary_sha256'],
          'raw_periodic_differences_retained': True,
          'interpretation': 'Evaluate the measured precision/thread effects and remaining explicit energy conventions separately. This runtime comparison is not chemical-model validation.'}
if not report['original_installed_binary_unchanged']:
    raise ValueError('Installed native binary changed during isolated diagnostic')
(AUDIT / 'double-static-comparison.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: report[k] for k in ['status', 'thread_energy_difference_1_minus_2_kj_mol',
                                        'thread_force_rms_vector_difference_kj_mol_nm', 'thread_coordinate_arrays_exact']}, indent=2), flush=True)
