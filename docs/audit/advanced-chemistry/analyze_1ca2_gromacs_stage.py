"""Independent saved-frame site geometry checks, without changing native outputs."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import mdtraj
import parmed

BASE = Path(__file__).resolve().parent
p = argparse.ArgumentParser()
p.add_argument('folder', type=Path)
a = p.parse_args()
folder = a.folder.resolve()
native = parmed.load_file(str(BASE / 'prototype-1ca2-zaff-motion/solvated.prmtop'))
exported = parmed.load_file(str(folder / 'system.top'))
identity = lambda model: [(x.name, x.residue.name, x.atomic_number) for x in model.atoms]
if identity(native) != identity(exported):
    raise ValueError('Native/export atom identity or order mismatch')
zn = next(x for x in native.atoms if x.residue.name == 'ZN6')
bonds = [b for b in native.bonds if zn in (b.atom1, b.atom2)]
angles = [x for x in native.angles if x.atom2.idx == zn.idx]
if len(bonds) != 4 or len(angles) != 6:
    raise ValueError('Expected native four-coordinate site terms absent')
adj = [[] for _ in native.atoms]
for bond in native.bonds:
    i, j = bond.atom1.idx, bond.atom2.idx
    adj[i].append(j)
    adj[j].append(i)
with mdtraj.formats.TRRTrajectoryFile(str(folder / 'stage.trr')) as stream:
    xyz, times, steps, boxes, *extra = stream._read(10000, None, get_forces=True)
forces = extra[-1]
if not np.isfinite(xyz).all() or not np.isfinite(forces).all():
    raise ValueError('Nonfinite saved frame')
status = json.loads((folder / 'stage-status.json').read_text())
if len(xyz) != status['stage_duration_ps'] + 1 or abs(float(times[-1]) - status['stage_duration_ps']) > 1e-5:
    raise ValueError('Saved frames do not span the complete requested stage')
ca = [x.idx for x in native.atoms if x.name == 'CA' and x.residue.idx < 256]
metrics = []
reference = None
for coords, box, t, force in zip(xyz.astype(float), boxes.astype(float), times, forces, strict=True):
    inv = np.linalg.inv(box)
    def image(delta):
        return delta - np.rint(delta @ inv) @ box
    whole = np.full_like(coords, np.nan)
    for anchor in range(len(coords)):
        if np.isfinite(whole[anchor]).all():
            continue
        whole[anchor] = coords[anchor]
        queue = [anchor]
        for i in queue:
            for j in adj[i]:
                expected = whole[i] + image(coords[j] - coords[i])
                if np.isnan(whole[j]).any():
                    whole[j] = expected
                    queue.append(j)
                elif not np.allclose(whole[j], expected, atol=1e-6, rtol=0):
                    raise ValueError('Inconsistent bonded cycle under periodic reconstruction')
    center = whole[ca] - whole[ca].mean(0)
    if reference is None:
        reference = center
    u, _, v = np.linalg.svd(center.T @ reference)
    sign = np.eye(3)
    sign[-1, -1] = np.linalg.det(u @ v)
    rmsd = float(np.sqrt(np.mean(np.sum((center @ u @ sign @ v - reference)**2, axis=1))) * 10)
    distances = [float(np.linalg.norm(whole[b.atom1.idx] - whole[b.atom2.idx]) * 10) for b in bonds]
    values = []
    for angle in angles:
        x = whole[angle.atom1.idx] - whole[angle.atom2.idx]
        y = whole[angle.atom3.idx] - whole[angle.atom2.idx]
        values.append(float(np.degrees(np.arccos(np.clip(np.dot(x, y) / np.linalg.norm(x) / np.linalg.norm(y), -1, 1)))))
    metrics.append({'time_ps': float(t), 'Zn_distances_angstrom': distances,
                    'Zn_angles_degrees': values, 'protein_CA_aligned_rmsd_angstrom': rmsd,
                    'maximum_force_kj_mol_nm': float(np.linalg.norm(force, axis=1).max())})
d = np.array([x['Zn_distances_angstrom'] for x in metrics])
angle_values = np.array([x['Zn_angles_degrees'] for x in metrics])
donor = lambda b: b.atom2 if b.atom1.idx == zn.idx else b.atom1
bounds = [(1.6, 2.8) if donor(b).name == 'O' else (1.6, 2.6) for b in bonds]
report = {'status': 'completed', 'frame_count': len(xyz), 'finite_saved_full_forces': True,
          'native_atom_order_and_identity_exact': True,
          'periodic_coordinates': 'Whole molecules reconstructed from native bonded graph; every bonded cycle checked',
          'meaning': 'Broad saved-frame geometry alarms and motion diagnostics, not chemical-model accuracy or thermodynamic convergence',
          'native_bonds': [{'atoms': [b.atom1.idx, b.atom2.idx], 'donor': donor(b).name,
                            'equilibrium_angstrom': b.type.req} for b in bonds],
          'native_angles': [{'atoms': [x.atom1.idx, x.atom2.idx, x.atom3.idx],
                             'equilibrium_degrees': x.type.theteq} for x in angles],
          'screening_bounds_angstrom': bounds,
          'broad_geometry_alarms_clear': all(np.all((d[:, i] >= lo) & (d[:, i] <= hi)) for i, (lo, hi) in enumerate(bounds)),
          'Zn_min_angstrom': d.min(0).tolist(), 'Zn_max_angstrom': d.max(0).tolist(),
          'Zn_mean_angstrom': d.mean(0).tolist(), 'Zn_angle_min_degrees': angle_values.min(0).tolist(),
          'Zn_angle_max_degrees': angle_values.max(0).tolist(), 'Zn_angle_mean_degrees': angle_values.mean(0).tolist(),
          'maximum_CA_rmsd_angstrom': max(x['protein_CA_aligned_rmsd_angstrom'] for x in metrics),
          'metrics': metrics, 'hashes': {n: hashlib.sha256((folder / n).read_bytes()).hexdigest()
                                       for n in ['system.top', 'stage.trr', 'stage-status.json']}}
(folder / 'independent-site-geometry.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: report[k] for k in ['frame_count', 'broad_geometry_alarms_clear', 'Zn_min_angstrom', 'Zn_max_angstrom', 'maximum_CA_rmsd_angstrom']}, indent=2))
