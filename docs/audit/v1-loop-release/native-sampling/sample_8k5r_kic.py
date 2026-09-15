"""One bounded native closure comparison; no full-atom model is admitted."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

from ost import io, geom
from promod3 import loop, modelling
import promod3


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(v):
    return [v.x, v.y, v.z]


def wrapped(v):
    return (v + math.pi) % (2 * math.pi) - math.pi


def residues(path):
    entity = io.LoadPDB(str(path))
    rows = [r for c in entity.chains if c.name == 'A' for r in c.residues]
    result = {r.number.num: r for r in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate chain A numbering prevents exact probe mapping')
    return entity, result


def coordinates(backbone, numbers, names):
    return [{'identity': ['A', str(n), '', names[n], atom],
             'xyz_angstrom': point(getattr(backbone, 'Get' + atom)(i))}
            for i, n in enumerate(numbers) for atom in ('N', 'CA', 'C', 'O')]


def basin(omega):
    return 'cis' if -math.pi / 2 < omega < math.pi / 2 else 'trans'


def compare(source, seeded, start, end, modeled, draws, pivot_sets, sampler, deadline):
    numbers = list(range(start, end + 1))
    if any(n not in seeded for n in numbers):
        raise ValueError('Missing seed residue in exact requested range')
    if any(n not in source for n in numbers if n not in modeled):
        raise ValueError('Missing observed source residue in requested range')
    if any(source[n].name != seeded[n].name for n in numbers if n not in modeled):
        raise ValueError('Observed source/seed sequence mismatch')
    if any(any(not seeded[n].FindAtom(a).IsValid() for a in ('N', 'CA', 'C', 'O')) for n in numbers):
        raise ValueError('Seed backbone inventory is incomplete')
    names = {n: seeded[n].name for n in numbers}
    sequence = ''.join(seeded[n].one_letter_code for n in numbers)
    indices = sampler.GetHistogramIndices(sequence)
    if len(indices) != len(numbers) - 2:
        raise ValueError('Native coil histogram inventory differs from internal residues')
    base = loop.BackboneList([seeded[n] for n in numbers])
    initial = coordinates(base, numbers, names)
    target = geom.DihedralAngle(*(source[start].FindAtom(a).pos for a in ('N', 'CA', 'C', 'O')))
    old = geom.DihedralAngle(base.GetN(0), base.GetCA(0), base.GetC(0), base.GetO(0))
    rotation = wrapped(target - old)
    base.SetPsiTorsion(0, base.GetPsiTorsion(0) + rotation)
    error = abs(wrapped(geom.DihedralAngle(base.GetN(0), base.GetCA(0), base.GetC(0), base.GetO(0)) - target))
    if error > 1e-5:
        raise ValueError('Observed n-stem carbonyl alignment failed')
    observed = {(n, a): point(source[n].FindAtom(a).pos)
                for n in numbers if n not in modeled for a in ('N', 'CA', 'C', 'O')}
    observed_omega = []
    for i, (n, following) in enumerate(zip(numbers, numbers[1:])):
        if n in modeled or following in modeled:
            continue
        omega = geom.DihedralAngle(source[n].FindAtom('CA').pos, source[n].FindAtom('C').pos,
                                   source[following].FindAtom('N').pos, source[following].FindAtom('CA').pos)
        observed_omega.append({'first': n, 'second': following, 'index': i,
                               'omega_radians': omega, 'basin': basin(omega)})
    closer = modelling.KIC(source[start], source[end])
    trials = []
    for pivots in pivot_sets:
        if len(pivots) != 3 or len(set(pivots)) != 3 or any(i < 1 or i > len(numbers)-2 for i in pivots):
            raise ValueError('Only three distinct internal KIC pivots are permitted')
        nonpivots = sorted(set(range(1, len(numbers)-1)) - set(pivots))
        # Reset to the aligned original seed for each independent draw. Draw
        # non-pivot torsions; changing solved pivots cannot expand the solutions.
        for draw in range(draws + 1):
            if time.monotonic() > deadline:
                raise TimeoutError('Native construction comparison exceeded its internal deadline')
            trial = {'pivots': pivots, 'draw_index': draw, 'sampling': [], 'solutions': []}
            try:
                candidate = base.Copy()
                if draw:
                    for index in nonpivots:
                        phi, psi = sampler.Draw(indices[index-1])
                        candidate.SetPhiPsiTorsion(index, phi, psi)
                        trial['sampling'].append({'index': index, 'residue': numbers[index],
                                                  'histogram_index': indices[index-1],
                                                  'phi_radians': phi, 'psi_radians': psi})
                trial['input_backbone_atoms'] = coordinates(candidate, numbers, names)
                # Native torsion rotations may introduce numerical roundoff;
                # they are never permitted to change an observed omega basin.
                trial['observed_omega_input_checks'] = [dict(row, candidate_basin=basin(candidate.GetOmegaTorsion(row['index'])))
                                                       for row in observed_omega]
                if any(r['basin'] != r['candidate_basin'] for r in trial['observed_omega_input_checks']):
                    raise ValueError('Torsion sampling changed an observed omega cis/trans basin')
                solutions = closer.Close(candidate, *pivots)
                if len(solutions) > 16:
                    raise ValueError('Native KIC exceeded its 16-solution bound')
                for s in solutions:
                    atoms = coordinates(s, numbers, names)
                    if not all(math.isfinite(v) for a in atoms for v in a['xyz_angstrom']):
                        raise ValueError('Native closure produced nonfinite coordinates')
                    displacements = [dict(identity=a['identity'], displacement_angstrom=math.dist(a['xyz_angstrom'], observed[int(a['identity'][1]), a['identity'][4]]))
                                     for a in atoms if int(a['identity'][1]) not in modeled]
                    checks = [dict(row, candidate_omega_radians=s.GetOmegaTorsion(row['index']),
                                   candidate_basin=basin(s.GetOmegaTorsion(row['index']))) for row in observed_omega]
                    maximum = max(r['displacement_angstrom'] for r in displacements)
                    trial['solutions'].append({'backbone_atoms': atoms,
                        'observed_backbone_displacements': displacements,
                        'maximum_observed_N_CA_C_O_displacement_angstrom': maximum,
                        'within_original_1A_observed_backbone_cap': maximum <= 1.0,
                        'observed_omega_checks': checks,
                        'observed_omega_basins_preserved': all(r['basin'] == r['candidate_basin'] for r in checks)})
            except Exception as exc:
                trial['error'] = str(exc)
            trials.append(trial)
    return {'start': start, 'end': end, 'sequence': sequence, 'modeled': sorted(modeled),
            'initial_backbone_atoms': initial,
            'carbonyl_alignment': {'source_dihedral_radians': target, 'initial_dihedral_radians': old,
                                   'native_psi_rotation_radians': rotation, 'error_radians': error},
            'trials': trials, 'raw_solutions': sum(len(t['solutions']) for t in trials),
            'solutions_within_observed_cap_and_omega_basins': sum(s['within_original_1A_observed_backbone_cap'] and s['observed_omega_basins_preserved'] for t in trials for s in t['solutions'])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--seed', type=Path, required=True)
    parser.add_argument('--control', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    files = [args.source, args.seed, args.control, Path(__file__)]
    hashes = {str(p): sha(p) for p in files}
    (args.output/'implementation.py').write_bytes(Path(__file__).read_bytes())
    source_entity, source = residues(args.source)
    seed_entity, seeded = residues(args.seed)
    control_entity, control = residues(args.control)
    sampler = loop.LoadTorsionSamplerCoil(seed=2026)
    started = time.monotonic()
    deadline = started + 120
    cases = []
    for name, original, seed, start, end, modeled, draws, pivots in [
        ('intact-control', control, control, 1, 6, set(), 0, [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]),
        ('original-stems', source, seeded, 176, 180, {177, 178, 179}, 0, [(1, 2, 3)]),
        ('two-observed-context-residues', source, seeded, 175, 181, {177, 178, 179}, 3,
         [(1, 3, 5), (2, 3, 4), (1, 2, 5), (1, 4, 5)]),
    ]:
        result = compare(original, seed, start, end, modeled, draws, pivots, sampler, deadline)
        result['name'] = name
        cases.append(result)
        (args.output / (name + '.json')).write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        print(json.dumps({'case': name, 'trials': len(result['trials']), 'solutions': result['raw_solutions'],
                          'within_cap_and_observed_omega': result['solutions_within_observed_cap_and_omega_basins']}), flush=True)
    if hashes != {str(p): sha(p) for p in files}:
        raise ValueError('Source or implementation changed during bounded native comparison')
    result = {'cases': cases, 'source_sha256': hashes, 'promod3_version': promod3.__version__,
              'coil_seed': 2026, 'sampled_conformer_draw_count': 12,
              'sampled_nonpivot_phi_psi_pair_count': 24,
              'elapsed_seconds': time.monotonic()-started,
              'bond_lengths_angles_or_force_constants_modified_by_workflow': False,
              'physical_model_validated': False, 'app_ready': False,
              'scope': 'Closure feasibility only. All native coordinates and failed checks retained; full-atom sidechains and environment have not been screened.'}
    (args.output/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
