"""Compare saved peptide geometry in rejected loops and completed controls.

Read-only diagnostic, with no admission-limit changes or parameter fitting.
Whole coordinates were saved every ps by these research runners.
"""
import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
from openmm import app, unit

from prepare_adduct import digest


def audit(inputs, run, output):
    output.mkdir(exist_ok=False)
    paths = [inputs/'solvated.prmtop', inputs/'solvated.inpcrd', inputs/'source-mapping.json', run/'trajectory.dcd']
    hashes = {str(p): digest(p) for p in paths}
    native = app.AmberPrmtopFile(str(inputs/'solvated.prmtop'))
    xyz = np.asarray(app.AmberInpcrdFile(str(inputs/'solvated.inpcrd')).positions.value_in_unit(unit.nanometer))
    trajectory = md.load_dcd(str(run/'trajectory.dcd'), top=md.Topology.from_openmm(native.topology))
    atom_sources = {row['native_index']: row['source'] for row in json.loads((inputs/'source-mapping.json').read_text())}
    definitions = []
    for a, b in native.topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        left, right = (a, b) if a.name == 'C' else (b, a)
        before = {atom.name: atom for atom in left.residue.atoms()}
        after = {atom.name: atom for atom in right.residue.atoms()}
        if 'CA' not in before or 'CA' not in after:
            continue
        members = [before['CA'], left, right, after['CA']]
        definitions.append({'indices': [atom.index for atom in members],
            'left': atom_sources.get(left.index), 'right': atom_sources.get(right.index),
            'residue_names': [left.residue.name, right.residue.name]})
    indices = [row['indices'] for row in definitions]
    angles = np.degrees(md.compute_dihedrals(trajectory, indices, periodic=False)).astype(float)
    original = md.Trajectory(xyz[None, :, :], trajectory.topology)
    initial = np.degrees(md.compute_dihedrals(original, indices, periodic=False))[0].astype(float)
    if not np.isfinite(angles).all() or not np.isfinite(initial).all():
        raise ValueError('Undefined peptide geometry')
    deviations = np.minimum(abs(angles), 180-abs(angles))
    initial_deviations = np.minimum(abs(initial), 180-abs(initial))
    per_bond = []
    for i, row in enumerate(definitions):
        per_bond.append({**row, 'initial_omega_degrees': float(initial[i]),
            'initial_deviation_degrees': float(initial_deviations[i]),
            'sampled_mean_deviation_degrees': float(deviations[:, i].mean()),
            'sampled_max_deviation_degrees': float(deviations[:, i].max()),
            'frames_over_35_degrees': int((deviations[:, i] > 35).sum()),
            'first_7_ps_omega_degrees': angles[:7, i].tolist(),
            'first_7_ps_deviation_degrees': deviations[:7, i].tolist()})
    np.savez_compressed(output/'angles.npz', omega_degrees=angles, deviation_degrees=deviations,
                        initial_omega_degrees=initial, indices=np.asarray(indices))
    (output/'bonds.json').write_text(json.dumps(per_bond, indent=2)+'\n')
    report = {'stage': 'saved_peptide_planarity_diagnostic', 'physical_model_validated': False,
        'app_ready': False, 'admission_limits_changed': False, 'frames': len(trajectory),
        'peptide_bonds': len(definitions), 'bond_frame_samples': int(deviations.size),
        'samples_over_35_degrees': int((deviations > 35).sum()),
        'bonds_ever_over_35_degrees': int(np.any(deviations > 35, axis=0).sum()),
        'maximum_deviation_degrees': float(deviations.max()),
        'median_deviation_degrees': float(np.median(deviations)),
        'percentile_99_deviation_degrees': float(np.percentile(deviations, 99)),
        'initial_bonds_over_35_degrees': int((initial_deviations > 35).sum()),
        'source_sha256': hashes,
        'output_sha256': {n: digest(output/n) for n in ['angles.npz', 'bonds.json']},
        'scope': 'All saved backbone peptide omega angles, including observed residues and any modeled residues. No acceptance is inferred from control excursions; static preparation bounds and trajectory statistics must be interpreted separately.'}
    if hashes != {str(p): digest(p) for p in paths}:
        raise ValueError('Source changed during audit')
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['input', 'run', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.input, args.run, args.output)
