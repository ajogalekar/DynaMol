"""Read-only aggregation of the bounded, zero-step periodic energy diagnostic."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

import mdtraj
import numpy as np

HERE = Path(__file__).resolve().parent
REVIEW = HERE / 'independent-review'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trr(path):
    with mdtraj.formats.TRRTrajectoryFile(str(path)) as handle:
        a = handle._read(1, None, get_forces=True)
    return a[0][0].astype(float), a[3][0].astype(float), a[-1][0].astype(float)


reference = np.load(HERE / 'openmm-reference.npz')
eref = json.loads((HERE / 'openmm-reference.json').read_text())['energy_kj_mol']
original = trr(HERE / 'matched-mesh/diagnostic.trr')
cases = [
    ('SIMD, 2 threads (original)', HERE / 'matched-mesh/diagnostic.trr',
     HERE / 'matched-mesh/native-initial.xvg', [-3, -2, -1]),
    ('SIMD, 1 thread', REVIEW / 'thread-1/single.trr', REVIEW / 'thread-1/energy.xvg', [1, 2, 3]),
    ('SIMD, 2 threads (repeat)', REVIEW / 'thread-2-repeat/repeat.trr',
     REVIEW / 'thread-2-repeat/energy.xvg', [1, 2, 3]),
    ('Plain C, 2 threads', REVIEW / 'scalar-thread-2/scalar.trr',
     REVIEW / 'scalar-thread-2/selected-energy.xvg', [1, 2, 3]),
]
records = []
for name, traj, energies, columns in cases:
    xyz, box, forces = trr(traj)
    e = np.loadtxt(energies)[columns]
    delta = forces - reference['forces_kj_mol_nm']
    delta_original = forces - original[2]
    assert np.isfinite(xyz).all() and np.isfinite(forces).all() and np.isfinite(e).all()
    records.append({
        'case': name, 'coulomb_sr_kj_mol': float(e[0]),
        'coulomb_reciprocal_kj_mol': float(e[1]), 'potential_kj_mol': float(e[2]),
        'potential_minus_openmm_reference_kj_mol': float(e[2] - eref),
        'max_abs_coordinate_difference_from_original_nm': float(np.max(np.abs(xyz - original[0]))),
        'max_abs_box_difference_from_original_nm': float(np.max(np.abs(box - original[1]))),
        'max_abs_coordinate_difference_from_openmm_nm': float(np.max(np.abs(xyz - reference['coords_nm']))),
        'force_rms_vector_difference_from_openmm_kj_mol_nm': float(np.sqrt(np.mean(np.sum(delta * delta, axis=1)))),
        'force_max_vector_difference_from_openmm_kj_mol_nm': float(np.max(np.linalg.norm(delta, axis=1))),
        'force_rms_vector_difference_from_original_kj_mol_nm': float(np.sqrt(np.mean(np.sum(delta_original**2, axis=1)))),
        'trajectory': str(traj.relative_to(HERE)), 'trajectory_sha256': digest(traj),
        'energy_file': str(energies.relative_to(HERE)), 'energy_file_sha256': digest(energies),
    })
assert all(r['max_abs_coordinate_difference_from_original_nm'] == 0 for r in records)
assert all(r['max_abs_box_difference_from_original_nm'] == 0 for r in records)
sources = json.loads((REVIEW / 'source-retrieval.json').read_text())
for source in sources:
    if 'file' in source:
        assert digest(REVIEW / source['file']) == source['sha256']
report = {
    'created_utc': datetime.now(timezone.utc).isoformat(),
    'scope': 'Three fresh nsteps=0 static evaluations of the exact same pre-existing TPR; no MD/QM, model edits, arbitrary offsets, or backend changes.',
    'absolute_periodic_energy_parity_validated': False,
    'finding': 'Thread-sensitive direct-space Coulomb energy dominates this numerical discrepancy. A reciprocal-mesh-only explanation is contradicted by the component and thread comparisons. The exact residual cause has not been isolated to a particular accumulation operation.',
    'tpr_sha256': digest(HERE / 'matched-mesh/diagnostic.tpr'),
    'openmm_reference_energy_kj_mol': eref,
    'cases': records,
    'reciprocal_bookkeeping_estimate': json.loads((REVIEW / 'reciprocal-self-check.json').read_text()),
    'sources': sources,
    'next_discriminator': 'If further absolute-energy parity is required, repeat these static inputs using a pinned double-precision GROMACS build with the same settings and compare direct-space components and forces. This has not been run.',
}
(HERE / 'review-periodic-localization.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(records, indent=2))
