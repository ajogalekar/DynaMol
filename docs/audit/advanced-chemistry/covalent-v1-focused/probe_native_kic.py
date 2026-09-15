"""Bounded native KIC feasibility probe; no protein candidate is accepted here."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

from ost import io, geom
from promod3 import modelling, loop
import promod3


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(vector):
    return [vector.x, vector.y, vector.z]


def run(source_path, seed_path, output, start, end, modeled, align_n_stem_carbonyl=False):
    output.mkdir(exist_ok=False)
    sources = [source_path, seed_path, Path(__file__)]
    hashes = {str(p): digest(p) for p in sources}
    (output / 'implementation.py').write_bytes(Path(__file__).read_bytes())
    source = io.LoadPDB(str(source_path))
    seed = io.LoadPDB(str(seed_path))
    original = {r.number.num: r for chain in source.chains if chain.name == 'A' for r in chain.residues}
    seeded = {r.number.num: r for chain in seed.chains if chain.name == 'A' for r in chain.residues}
    numbers = list(range(start, end + 1))
    if len(numbers) < 3 or len(numbers) > 8 or any(n not in original or n not in seeded for n in numbers):
        raise ValueError('Probe requires 3–8 explicitly mapped consecutive source residues')
    if any(original[n].name != seeded[n].name for n in numbers):
        raise ValueError('Seed and source sequences differ')
    backbone = loop.BackboneList([seeded[n] for n in numbers])
    carbonyl_alignment = None
    if align_n_stem_carbonyl:
        target = geom.DihedralAngle(*(original[start].FindAtom(atom).pos for atom in ('N', 'CA', 'C', 'O')))
        before = geom.DihedralAngle(backbone.GetN(0), backbone.GetCA(0), backbone.GetC(0), backbone.GetO(0))
        change = (target - before + math.pi) % (2 * math.pi) - math.pi
        backbone.SetPsiTorsion(0, backbone.GetPsiTorsion(0) + change)
        after = geom.DihedralAngle(backbone.GetN(0), backbone.GetCA(0), backbone.GetC(0), backbone.GetO(0))
        error = abs((after - target + math.pi) % (2 * math.pi) - math.pi)
        carbonyl_alignment = {'original_dihedral_radians': target, 'seed_dihedral_radians': before,
                              'psi_rotation_radians': change, 'final_angular_error_radians': error,
                              'method': 'Native torsion rotation aligns the observed N/CA/C/O dihedral before KIC; no Cartesian atom displacement or energy parameter is fitted.'}
        if error > 1e-5:
            raise ValueError('Native psi rotation did not preserve the observed carbonyl orientation')
    closer = modelling.KIC(original[start], original[end])
    reference = {(i, atom): point(original[n].FindAtom(atom).pos)
                 for i, n in enumerate(numbers) for atom in ('N', 'CA', 'C', 'O')}
    trials = []
    started = time.monotonic()
    # Three pivots are required. Include endpoint pivots; record a native error
    # explicitly if this runtime disallows a combination. No runtime fallback.
    for pivots in itertools.combinations(range(len(numbers)), 3):
        entry = {'pivots': pivots, 'solutions': []}
        try:
            solutions = closer.Close(backbone, *pivots)
            if len(solutions) > 16:
                raise ValueError('Native KIC exceeded its documented solution bound')
            for solution in solutions:
                atoms = [{'residue': numbers[i], 'name': atom,
                          'xyz_angstrom': point(getattr(solution, 'Get' + atom)(i))}
                         for i in range(len(numbers)) for atom in ('N', 'CA', 'C', 'O')]
                if not all(math.isfinite(v) for row in atoms for v in row['xyz_angstrom']):
                    raise ValueError('Nonfinite KIC coordinates')
                stem = [math.dist(point(getattr(solution, 'Get' + atom)(i)), reference[i, atom])
                        for i in (0, len(numbers) - 1) for atom in ('N', 'CA', 'C')]
                observed = [math.dist(row['xyz_angstrom'], reference[numbers.index(row['residue']), row['name']])
                            for row in atoms if row['residue'] not in modeled]
                entry['solutions'].append({'backbone_atoms': atoms,
                    'maximum_stem_N_CA_C_displacement_angstrom': max(stem),
                    'maximum_observed_backbone_displacement_angstrom': max(observed),
                    'within_original_observed_backbone_cap': max(observed) <= 1.0})
        except Exception as exc:
            entry['error'] = str(exc)
        trials.append(entry)
    if hashes != {str(p): digest(p) for p in sources}:
        raise ValueError('KIC probe source changed during execution')
    result = {'stage': 'native_kic_feasibility_probe_complete', 'trials': trials,
              'source_sha256': hashes, 'promod3_version': promod3.__version__,
              'residue_numbers': numbers, 'modeled_residue_numbers': sorted(modeled),
              'n_stem_carbonyl_alignment': carbonyl_alignment,
              'solution_count': sum(len(t['solutions']) for t in trials),
              'solutions_within_observed_backbone_cap': sum(s['within_original_observed_backbone_cap'] for t in trials for s in t['solutions']),
              'elapsed_seconds': time.monotonic() - started,
              'app_ready': False, 'physical_model_validated': False,
              'scope': 'Native closure feasibility only. Backbone displacements do not check sidechains, source stereochemistry, full retained environment or final peptide/backbone quality. No seed or prepared output is admitted.'}
    (output / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['trials', 'source_sha256']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['source', 'seed', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--start', type=int, required=True)
    parser.add_argument('--end', type=int, required=True)
    parser.add_argument('--modeled', type=int, nargs='*', default=[])
    parser.add_argument('--align-n-stem-carbonyl', action='store_true')
    args = parser.parse_args()
    run(args.source, args.seed, args.output, args.start, args.end, set(args.modeled), args.align_n_stem_carbonyl)
