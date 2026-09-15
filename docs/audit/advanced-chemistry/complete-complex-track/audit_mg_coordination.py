"""Read-only complete-complex coordination diagnostic for saved research runs.

Reports donor identity, shell membership and all donor-Mg-donor angles. These
geometric summaries are not force-field validation or preparation admission.
No coordinates, topology parameters or source artifacts are modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import platform
import shutil
import time

import gemmi
import mdtraj as md
import numpy as np
import parmed


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def minimum_image(delta, box):
    """Double-precision image vectors; fail outside this audit's orthogonal scope."""
    delta = np.asarray(delta, dtype=float)
    box = np.asarray(box, dtype=float)
    if not np.isfinite(delta).all() or not np.isfinite(box).all():
        raise ValueError('Non-finite geometry')
    if box.shape != (3, 3) or np.max(abs(box - np.diag(np.diag(box)))) > 1e-7:
        raise ValueError('This audit requires orthogonal saved boxes')
    if np.min(np.diag(box)) <= 0:
        raise ValueError('Non-positive box')
    lengths = np.diag(box)
    return delta - np.rint(delta / lengths) * lengths


def angle_values(vectors):
    pairs = list(itertools.combinations(range(len(vectors)), 2))
    v = np.asarray(vectors, dtype=float)
    lengths = np.linalg.norm(v, axis=1)
    if np.min(lengths) < 1e-8:
        raise ValueError('Zero-length coordination vector')
    u = v / lengths[:, None]
    return np.array([np.degrees(np.arccos(np.clip(np.dot(u[i], u[j]), -1, 1)))
                     for i, j in pairs])


def statistics(values):
    x = np.asarray(values, dtype=float)
    if not x.size or not np.isfinite(x).all():
        raise ValueError('Empty/non-finite summary')
    return {k: float(v) for k, v in {
        'mean': x.mean(), 'std_population': x.std(), 'min': x.min(),
        'q05': np.quantile(x, .05), 'median': np.median(x),
        'q95': np.quantile(x, .95), 'max': x.max()}.items()}


def source_record(item):
    s = item['source']
    if 'residue' in s:
        return (*s['residue'], s['atom'])
    if 'source' in s:
        identity = s['source']['identity']
        return (identity['author_chain'], identity['author_residue_id'],
                identity.get('insertion_code') or '', s['resname'], s['source']['atom'])
    return (s['chain'], s['resid'], '', s['resname'], s['atom'])


def cif_atoms(path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    columns = ['label_asym_id', 'auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code',
               'label_comp_id', 'label_atom_id', 'type_symbol', 'Cartn_x',
               'Cartn_y', 'Cartn_z', 'label_alt_id', 'occupancy']
    result = []
    for row in block.find('_atom_site.', columns):
        result.append(dict(label_chain=row[0], author_chain=row[1], resid=row[2],
                           insertion='' if row[3] in ('.', '?') else row[3],
                           resname=row[4], atom=row[5], element=row[6].upper(),
                           xyz=np.array([float(row[i]) for i in (7, 8, 9)]),
                           altloc=row[10], occupancy=float(row[11])))
    return result


def resolve_source(item, atoms):
    chain, resid, ins, resname, atom = source_record(item)
    hits = [x for x in atoms if chain in (x['label_chain'], x['author_chain'])
            and (resid, ins, resname, atom) ==
            (x['resid'], x['insertion'], x['resname'], x['atom'])]
    if len(hits) != 1:
        raise ValueError(f'Ambiguous or absent deposited identity: {source_record(item)}')
    return hits[0]


def software_checks():
    xyz = np.array([[1., 0, 0], [-1., 0, 0], [0, 1., 0],
                    [0, -1., 0], [0, 0, 1.], [0, 0, -1.]])
    a = angle_values(xyz)
    assert np.count_nonzero(np.isclose(a, 90)) == 12
    assert np.count_nonzero(np.isclose(a, 180)) == 3
    assert np.allclose(angle_values(xyz * np.arange(1, 7)[:, None]), a)
    boundary = minimum_image([[2.9, -.2, 3.1]], np.diag([3., 3., 3.]))
    assert np.allclose(boundary, [[-.1, -.2, .1]])
    errors = 0
    for box in (np.array([[3, .2, 0], [0, 3, 0], [0, 0, 3]]), np.zeros((3, 3))):
        try:
            minimum_image([[0, 0, 0]], box)
        except ValueError:
            errors += 1
    assert errors == 2
    try:
        angle_values([[0, 0, 0], [1, 0, 0]])
    except ValueError:
        errors += 1
    assert errors == 3
    return {'ideal_octahedron': True, 'angles_invariant_to_radial_scaling': True,
            'periodic_boundary': True, 'invalid_geometry_refusals': 3}


def audit(root, output):
    started = time.monotonic()
    inputs = root / '6oim-complex-v3'
    run = root / '6oim-stability-v2'
    gmx = root / '6oim-gromacs-continuation-v2'
    cif = root / 'inputs/6OIM.cif'
    used = [inputs / n for n in ('solvated.prmtop', 'solvated.inpcrd',
                                 'source-mapping.json', 'result.json')]
    used += [run / n for n in ('result.json', 'trajectory.dcd', 'observations.csv',
                               'metal-sites.json', 'saved-run-audit.json')]
    used += [gmx / n for n in ('result.json', 'stage.trr', 'handoff.json', 'input.gro')]
    used += [cif]
    hashes = {str(p): digest(p) for p in used}
    result = json.loads((run / 'result.json').read_text())
    assembly = json.loads((inputs / 'result.json').read_text())
    gmx_result = json.loads((gmx / 'result.json').read_text())
    handoff = json.loads((gmx / 'handoff.json').read_text())
    if not handoff['source_atom_order_and_bonds_preserved']:
        raise ValueError('Native handoff does not confirm source atom order/bonds')
    for name, expected in result['input_sha256'].items():
        if digest(inputs / name) != expected:
            raise ValueError('OpenMM source artifact changed: ' + name)
    for path, expected in gmx_result['source_hashes'].items():
        if digest(path) != expected:
            raise ValueError('GROMACS source artifact changed: ' + path)
    if digest(gmx / 'stage.trr') != gmx_result['output_sha256']['stage.trr']:
        raise ValueError('GROMACS trajectory changed')
    saved = json.loads((run / 'saved-run-audit.json').read_text())
    if digest(run / 'trajectory.dcd') != saved['output_sha256']['trajectory.dcd']:
        raise ValueError('OpenMM trajectory changed')
    if not assembly['all_deposited_environment_retained'] or assembly['excluded_environment']:
        raise ValueError('Expected complete retained deposited environment')
    native = parmed.load_file(str(inputs / 'solvated.prmtop'), str(inputs / 'solvated.inpcrd'))
    mapping = json.loads((inputs / 'source-mapping.json').read_text())
    mapped = {r['native_index']: r for r in mapping}
    if len(mapped) != len(mapping):
        raise ValueError('Duplicate mapped native atom')
    if any(native.atoms[i].name != item['name'] or
           native.atoms[i].element_name.upper() != item['element'].upper()
           for i, item in mapped.items()):
        raise ValueError('Mapped native atom name/element mismatch')
    xyz0 = np.asarray(native.coordinates) / 10
    deposited = cif_atoms(cif)
    mg = [a.idx for a in native.atoms if a.atomic_number == 12]
    if len(mg) != 1:
        raise ValueError('Expected one retained Mg')
    mg = mg[0]
    mg_cif = resolve_source(mapped[mg], deposited)
    candidates = np.array([a.idx for a in native.atoms if a.atomic_number in (7, 8)])
    dist0 = np.linalg.norm(xyz0[candidates] - xyz0[mg], axis=1) * 10
    donors = candidates[dist0 < 2.8]
    if len(donors) != 6:
        raise ValueError('Expected six initial donors')
    recorded = json.loads((run / 'metal-sites.json').read_text())
    if set(donors) != {d['index'] for d in recorded[0]['donors']}:
        raise ValueError('Initial donor identity disagrees with saved native monitor')
    pairs = list(itertools.combinations(range(6), 2))
    source_vectors = []
    labels = []
    for i in donors:
        observed = resolve_source(mapped[int(i)], deposited)
        if observed['element'] != native.atoms[int(i)].element_name.upper():
            raise ValueError('Source/native donor element mismatch')
        source_vectors.append((observed['xyz'] - mg_cif['xyz']) / 10)
        labels.append({'native_index': int(i), 'source': list(source_record(mapped[int(i)])),
                       'label_chain': observed['label_chain'], 'author_chain': observed['author_chain'],
                       'occupancy': observed['occupancy'], 'altloc': observed['altloc']})
    source_vectors = np.asarray(source_vectors)
    if np.max(abs(source_vectors - (xyz0[donors] - xyz0[mg]))) > .00011:
        raise ValueError('Assembly changed the observed Mg coordination vectors')
    crystal_distances = np.linalg.norm(source_vectors, axis=1) * 10
    crystal_angles = angle_values(source_vectors)
    trans = crystal_angles > 135
    if trans.sum() != 3:
        raise ValueError('Observed geometry does not define three opposite donor pairs')
    ideals = np.where(trans, 180., 90.)

    def identity(i):
        if int(i) in mapped:
            return list(source_record(mapped[int(i)]))
        a = native.atoms[int(i)]
        return ['added-solvent-or-counterion', a.residue.idx, a.residue.name, a.name]

    summaries, frame_rows, crosschecks = {}, [], {}
    trajectories = [('OpenMM', md.load_dcd(str(run / 'trajectory.dcd'), top=str(inputs / 'solvated.prmtop'))),
                    ('GROMACS', md.load_trr(str(gmx / 'stage.trr'), top=str(inputs / 'solvated.prmtop')))]
    omm_endpoint = trajectories[0][1]
    gmx_start = trajectories[1][1]
    gro = md.load(str(gmx / 'input.gro'))
    # GROMACS wraps with its actual GRO box, whose text precision differs from
    # the source DCD. Compare through that saved input rather than incorrectly
    # removing periodic translations with the pre-export box.
    handoff_coordinate_error = float(np.max(abs(minimum_image(
        gro.xyz[0] - omm_endpoint.xyz[-1], gro.unitcell_vectors[0]))))
    native_start_error = float(np.max(abs(minimum_image(
        gmx_start.xyz[0] - gro.xyz[0], gro.unitcell_vectors[0]))))
    native_box_error = float(np.max(abs(gmx_start.unitcell_vectors[0] - gro.unitcell_vectors[0])))
    if max(handoff_coordinate_error, native_start_error, native_box_error) > 2e-6:
        raise ValueError('Saved OpenMM/GRO/GROMACS coordinate chain does not match')
    box_fields = (gmx / 'input.gro').read_text().splitlines()[-1].split()
    if len(box_fields) != 3:
        raise ValueError('Expected three orthogonal GRO box lengths')
    source_lengths = np.diag(np.asarray(handoff['source_box_nm']))
    box_text_error = np.abs(np.array([float(x) for x in box_fields]) - source_lengths)
    half_quantum = np.array([.5 * 10 ** -len(x.split('.')[1]) for x in box_fields])
    if np.any(box_text_error > half_quantum + 1e-12):
        raise ValueError('GRO box changed beyond its written decimal precision')
    observations = list(csv.DictReader((run / 'observations.csv').open()))
    for engine, trajectory in trajectories:
        expected = 100 if engine == 'OpenMM' else 21
        if len(trajectory) != expected or trajectory.n_atoms != len(native.atoms):
            raise ValueError('Unexpected saved trajectory dimensions')
        if not np.isfinite(trajectory.xyz).all():
            raise ValueError('Non-finite saved coordinates')
        times = np.array([float(r['time_ps']) for r in observations]) if engine == 'OpenMM' else trajectory.time
        phases = [r['stage'] for r in observations] if engine == 'OpenMM' else ['unrestrained continuation'] * expected
        if len(times) != expected or not np.allclose(times, np.arange(1, 101) if engine == 'OpenMM' else np.arange(21)):
            raise ValueError('Saved frame time mapping failed')
        vectors = np.array([minimum_image(x[donors] - x[mg], box)
                            for x, box in zip(trajectory.xyz, trajectory.unitcell_vectors)])
        distances = np.linalg.norm(vectors, axis=2) * 10
        angles = np.array([angle_values(v) for v in vectors])
        all_distances = md.compute_distances(trajectory, [[mg, int(i)] for i in candidates], periodic=True) * 10
        md_distances = md.compute_distances(trajectory, [[mg, int(i)] for i in donors], periodic=True) * 10
        md_angles = np.degrees(md.compute_angles(trajectory,
            [[int(donors[i]), mg, int(donors[j])] for i, j in pairs], periodic=True))
        d_error = float(np.max(abs(distances - md_distances)))
        a_error = float(np.max(abs(angles - md_angles)))
        if d_error > 2e-5 or a_error > .003:
            raise ValueError('Independent distance/angle comparison failed')
        octa_rms = np.sqrt(np.mean((angles - ideals) ** 2, axis=1))
        source_rms = np.sqrt(np.mean((angles - crystal_angles) ** 2, axis=1))
        nearest = candidates[np.argsort(all_distances, axis=1)[:, :6]]
        donor_sets_retained = np.array([set(row) == set(donors) for row in nearest])
        cutoffs = (2.4, 2.8, 3.2)
        shell_counts = {str(c): (all_distances < c).sum(axis=1) for c in cutoffs}
        outer_min = np.min(all_distances[:, [i not in donors for i in candidates]], axis=1)
        encountered = set(candidates[np.any(all_distances < max(cutoffs), axis=0)])
        for n in range(expected):
            row = {'engine': engine, 'time_ps': float(times[n]), 'phase': phases[n],
                   'octahedral_angle_rms_degrees': float(octa_rms[n]),
                   'deposited_angle_rms_difference_degrees': float(source_rms[n]),
                   'original_six_are_nearest_six': bool(donor_sets_retained[n]),
                   'nearest_nonoriginal_donor_angstrom': float(outer_min[n])}
            row.update({f'coordination_number_{c}_angstrom': int(v[n]) for c, v in shell_counts.items()})
            row.update({f'donor_{int(i)}_angstrom': float(distances[n, j]) for j, i in enumerate(donors)})
            row.update({f'angle_{int(donors[i])}_{int(donors[j])}_degrees': float(angles[n, k])
                        for k, (i, j) in enumerate(pairs)})
            frame_rows.append(row)
        group_reports = []
        groups = [('all saved frames', np.ones(expected, dtype=bool))]
        groups += [(p, np.array([x == p for x in phases])) for p in dict.fromkeys(phases)]
        for phase, mask in groups:
            group_reports.append({'phase': phase, 'frames': int(mask.sum()),
                'time_ps': [float(times[mask].min()), float(times[mask].max())],
                'donors': [dict(label, distance_angstrom=statistics(distances[mask, j]),
                               mean_change_from_deposited_angstrom=float(distances[mask, j].mean() - crystal_distances[j]))
                           for j, label in enumerate(labels)],
                'angles': [{'donors': [int(donors[i]), int(donors[j])],
                            'observed_angle_degrees': float(crystal_angles[k]),
                            'assigned_reference': 'trans' if trans[k] else 'cis',
                            'angle_degrees': statistics(angles[mask, k])} for k, (i, j) in enumerate(pairs)],
                'octahedral_angle_rms_degrees': statistics(octa_rms[mask]),
                'deposited_angle_rms_difference_degrees': statistics(source_rms[mask]),
                'nearest_nonoriginal_donor_angstrom': statistics(outer_min[mask]),
                'original_six_nearest_six_frames': int(donor_sets_retained[mask].sum()),
                'coordination_number_cutoff_sensitivity': {
                    c: {str(int(v)): int(np.count_nonzero(counts[mask] == v)) for v in np.unique(counts[mask])}
                    for c, counts in shell_counts.items()}})
        summaries[engine] = {'phases': group_reports,
                             'atoms_ever_within_3_2_angstrom': [dict(native_index=int(i), identity=identity(i)) for i in sorted(encountered)]}
        crosschecks[engine] = {'maximum_distance_difference_angstrom': d_error,
                               'maximum_angle_difference_degrees': a_error,
                               'frames': expected}
    gdp = next(r for r in native.residues if r.name == 'GDP')
    cov = next(r for r in native.residues if r.name == 'COV')
    if len(gdp.atoms) != 40 or any(a.atomic_number <= 0 for a in gdp.atoms):
        raise ValueError('Corrected GDP element inventory lost')
    if any(digest(path) != sha for path, sha in hashes.items()):
        raise ValueError('A source artifact changed while audited')
    report = {'stage': 'read_only_complete_complex_coordination_audit',
              'case': '6OIM', 'physical_model_validated': False, 'app_ready': False,
              'new_dynamics_or_parameters': False, 'input_sha256': hashes,
              'complete_complex': {'atoms': len(native.atoms), 'dry_atoms': assembly['dry_atoms'],
                 'retained_deposited_environment': assembly['retained_environment'],
                 'excluded_deposited_environment': assembly['excluded_environment'],
                 'all_deposited_environment_retained': True,
                 'covalent_adduct_atoms': len(cov.atoms), 'covalent_adduct_charge_e': sum(a.charge for a in cov.atoms),
                 'gdp_atoms': len(gdp.atoms), 'gdp_charge_e': sum(a.charge for a in gdp.atoms),
                 'mg_charge_e': native.atoms[mg].charge, 'mg_source': list(source_record(mapped[mg]))},
              'source_geometry': {'donors': [dict(label, distance_angstrom=float(crystal_distances[j]))
                                           for j, label in enumerate(labels)],
                                  'octahedral_angle_rms_degrees': float(np.sqrt(np.mean((crystal_angles - ideals) ** 2)))},
              'candidate_donor_atoms_searched': len(candidates), 'engines': summaries,
              'saved_engine_boundary': {'source_atom_order_and_bonds_preserved': True,
                                        'openmm_to_gro_coordinate_difference_nm': handoff_coordinate_error,
                                        'gro_to_native_first_frame_coordinate_difference_nm': native_start_error,
                                        'gro_to_native_first_frame_box_difference_nm': native_box_error,
                                        'source_to_gro_box_rounding_nm': box_text_error.tolist(),
                                        'gro_written_box_half_quantum_nm': half_quantum.tolist()},
              'independent_geometry_crosschecks': crosschecks, 'software_checks': software_checks(),
              'runtime': {'python': platform.python_version(), 'numpy': np.__version__, 'mdtraj': md.__version__,
                          'parmed': parmed.__version__, 'gemmi': gemmi.__version__},
              'limits': ['Saved-frame, short-run descriptive analysis; no exchange kinetics or equilibrium ensemble claim.',
                         'GROMACS continues the same model; it is not an independent physical replicate.',
                         'All O/N atoms are searched, including bulk water; cutoff counts are sensitivity diagnostics, not bond assignments.',
                         'Octahedral angular RMS is a descriptive ideal-geometry comparison, not an acceptance test.',
                         'Deposited coordinates are a reference pose, not an exact solution-phase distance target.',
                         'No force-field selection or adjustment is made by this analysis.'],
              'elapsed_seconds': time.monotonic() - started}
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, output / 'implementation.py')
    with (output / 'frames.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(frame_rows[0]))
        writer.writeheader()
        writer.writerows(frame_rows)
    (output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
    (output / 'manifest.json').write_text(json.dumps({p.name: digest(p) for p in output.iterdir() if p.is_file()}, indent=2) + '\n')
    print(json.dumps({'output': str(output), 'elapsed_seconds': report['elapsed_seconds'],
                      'independent_geometry_crosschecks': crosschecks}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.root, args.output)
