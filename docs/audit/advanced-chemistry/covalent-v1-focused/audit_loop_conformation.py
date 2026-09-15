"""Saved local backbone angles and nearby atoms; no conformation admission."""
import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import parmed

from prepare_adduct import digest


def audit(inputs, run, output, chain='A', left='750', right='751'):
    paths = [inputs/'solvated.prmtop', inputs/'solvated.inpcrd', inputs/'source-mapping.json',
             run/'result.json', run/'trajectory.dcd']
    hashes = {str(p): digest(p) for p in paths}
    result = json.loads((run/'result.json').read_text())
    for name, expected in result['input_sha256'].items():
        if digest(inputs/name) != expected:
            raise ValueError('Run model input changed')
    native = parmed.load_file(str(inputs/'solvated.prmtop'), xyz=str(inputs/'solvated.inpcrd'))
    traj = md.load_dcd(str(run/'trajectory.dcd'), top=str(inputs/'solvated.prmtop'))
    if len(traj) != 100 or not np.isfinite(traj.xyz).all():
        raise ValueError('Expected completed, finite 100-frame trajectory')
    labels = {}
    for row in json.loads((inputs/'source-mapping.json').read_text()):
        source = row['source'].get('residue')
        if source:
            labels.setdefault(native.atoms[row['native_index']].residue.idx, set()).add(tuple(source))

    def residue(number):
        hits = [native.residues[i] for i, sources in labels.items()
                if any(s[0] == chain and s[1] == number for s in sources)]
        if len(hits) != 1 or len(labels[hits[0].idx]) != 1:
            raise ValueError('Selected residue source identity is ambiguous')
        return hits[0]

    def atom(res, name):
        hits = [a for a in res.atoms if a.name == name]
        if len(hits) != 1:
            raise ValueError('Expected one named atom: '+name)
        return hits[0]

    def identity(a):
        sources = labels.get(a.residue.idx, set())
        return {'source_residues': sorted(sources), 'native_residue_index': a.residue.idx,
                'native_residue_name': a.residue.name, 'name': a.name, 'native_index': a.idx}

    def previous(res):
        return next(a for a in atom(res, 'N').bond_partners if a.name == 'C' and a.residue is not res)

    def following(res):
        return next(a for a in atom(res, 'C').bond_partners if a.name == 'N' and a.residue is not res)

    first, second = residue(left), residue(right)
    if following(first).residue is not second:
        raise ValueError('Requested peptide is not connected in the native topology')
    region = [previous(first).residue, first, second, following(second).residue]
    definitions = []
    for res in region:
        n, ca, c = (atom(res, name) for name in ['N', 'CA', 'C'])
        for name, atoms in [('phi', [previous(res), n, ca, c]),
                            ('psi', [n, ca, c, following(res)]),
                            ('omega', [ca, c, following(res), atom(following(res).residue, 'CA')])]:
            definitions.append({'kind': name, 'residue': sorted(labels[res.idx]),
                                'indices': [a.idx for a in atoms]})
    indices = [d['indices'] for d in definitions]
    angles = md.compute_dihedrals(traj, indices, periodic=False)
    initial = md.Trajectory(np.asarray(native.coordinates)[None, :, :]*.1, traj.topology)
    original = md.compute_dihedrals(initial, indices, periodic=False)[0]
    rows = []
    for i, definition in enumerate(definitions):
        windows = []
        for lo, hi in [(1, 10), (11, 30), (31, 100)]:
            values = angles[lo-1:hi, i]
            mean = np.mean(np.exp(1j*values))
            windows.append({'ps': [lo, hi], 'circular_mean_degrees': float(np.degrees(np.angle(mean))),
                            'circular_resultant_length': float(abs(mean))})
        rows.append(definition | {'initial_degrees': float(np.degrees(original[i])), 'windows': windows})
    contact_rows = []
    targets = [atom(first, 'O'), atom(second, 'N')]
    targets += [a for a in atom(second, 'N').bond_partners if a.atomic_number == 1]
    targets += [a for a in second.atoms if a.atomic_number == 8 and a.name != 'O']
    for target in targets:
        excluded = {target}
        edge = {target}
        for _ in range(3):
            edge = {n for a in edge for n in a.bond_partners}-excluded
            excluded |= edge
        others = [a for a in native.atoms if a.atomic_number > 1 and a not in excluded]
        pairs = [[target.idx, a.idx] for a in others]
        distances = md.compute_distances(traj, pairs, periodic=True)*10
        threshold = 2.5 if target.atomic_number == 1 else 3.0
        ranked = sorted(range(len(others)), key=lambda i: (-int((distances[:, i] < threshold).sum()),
                                                         float(distances[:, i].min())))[:12]
        nearest = []
        for i in ranked:
            values = distances[:, i]
            nearest.append({'atom': identity(others[i]), 'minimum_angstrom': float(values.min()),
                            'mean_angstrom': float(values.mean()),
                            'saved_frames_within_cutoff': int((values < threshold).sum()),
                            'final_70_ps_frames_within_cutoff': int((values[30:] < threshold).sum())})
        contact_rows.append({'target': identity(target), 'cutoff_angstrom': threshold, 'nearest': nearest})
    if hashes != {str(p): digest(p) for p in paths}:
        raise ValueError('Source changed during conformation audit')
    report = {'stage': 'saved_local_conformation_and_proximity', 'case': result['case'],
              'source_diagnostic_only': result.get('diagnostic_only', False),
              'physical_model_validated': False, 'app_ready': False,
              'frames': len(traj), 'angles': rows, 'contacts': contact_rows, 'source_sha256': hashes,
              'scope': 'Measured dihedrals and periodic proximity only. Neighbors within three covalent bonds are excluded from contact lists. Cutoffs define descriptive proximity, not accepted hydrogen bonds or a Ramachandran quality criterion. No native conformation or force-field accuracy is inferred.'}
    output.mkdir(exist_ok=False)
    np.savez_compressed(output/'backbone-angles.npz', degrees=np.degrees(angles), initial_degrees=np.degrees(original))
    report['angle_array_sha256'] = digest(output/'backbone-angles.npz')
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['angles', 'contacts', 'source_sha256']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['input', 'run', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--chain', default='A')
    parser.add_argument('--left', default='750')
    parser.add_argument('--right', default='751')
    args = parser.parse_args()
    audit(args.input, args.run, args.output, args.chain, args.left, args.right)
