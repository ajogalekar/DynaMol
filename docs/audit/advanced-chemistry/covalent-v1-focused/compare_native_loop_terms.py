"""Read-only comparison of native terms touching selected protein residues."""
import argparse
import json
from pathlib import Path

import parmed

from prepare_adduct import digest


def inventory(folder, residue_ids):
    topology = parmed.load_file(str(folder/'solvated.prmtop'))
    mapping = json.loads((folder/'source-mapping.json').read_text())
    labels = {}
    for row in mapping:
        source = row['source'].get('residue')
        if source is None:
            continue
        residue = topology.atoms[row['native_index']].residue.idx
        labels.setdefault(residue, set()).add(tuple(source))
    selected = set()
    found = set()
    for i, names in labels.items():
        for source in names:
            if source[0] == 'A' and source[1] in residue_ids:
                if len(names) != 1:
                    raise ValueError('Ambiguous source residue mapping at selected site')
                selected.add(i)
                found.add(source[1])
    if found != set(residue_ids):
        raise ValueError('Every requested chain A residue must be present')

    def key(atom):
        names = labels.get(atom.residue.idx, set())
        if len(names) != 1:
            raise ValueError('Ambiguous source identity at a local bonded term')
        return list(next(iter(names)))+[atom.name]

    atoms = [a for a in topology.atoms if a.residue.idx in selected]
    records = {'atoms': [{'identity': key(a), 'type': a.type, 'atomic_number': a.atomic_number,
                         'mass': a.mass, 'charge': a.charge, 'rmin': a.rmin, 'epsilon': a.epsilon}
                        for a in atoms], 'bonds': [], 'angles': [], 'proper_torsions': [], 'improper_torsions': []}
    for collection, count in [('bonds', 2), ('angles', 3), ('dihedrals', 4)]:
        for term in getattr(topology, collection):
            members = [getattr(term, 'atom'+str(i)) for i in range(1, count+1)]
            if not any(a.residue.idx in selected for a in members):
                continue
            identities = [key(a) for a in members]
            kind = collection
            if count == 2:
                params = {'k': term.type.k, 'req': term.type.req}
            elif count == 3:
                params = {'k': term.type.k, 'theteq': term.type.theteq}
            else:
                kind = 'improper_torsions' if term.improper else 'proper_torsions'
                types = term.type if isinstance(term.type, list) else [term.type]
                params = [{'phi_k': t.phi_k, 'per': t.per, 'phase': t.phase,
                           'scee': t.scee, 'scnb': t.scnb} for t in types]
            if kind != 'improper_torsions':
                identities = min(identities, identities[::-1])
            records[kind].append({'identities': identities, 'parameters': params})
    for rows in records.values():
        rows.sort(key=lambda row: json.dumps(row, sort_keys=True))
    return records


def compare(baseline, comparison, residues, output):
    paths = [folder/name for folder in [baseline, comparison]
             for name in ['solvated.prmtop', 'source-mapping.json']]
    hashes = {str(p): digest(p) for p in paths}
    a = inventory(baseline, residues)
    b = inventory(comparison, residues)
    equality = {kind: a[kind] == b[kind] for kind in a}
    if hashes != {str(p): digest(p) for p in paths}:
        raise ValueError('Input topology or atom mapping changed during comparison')
    report = {'stage': 'native_local_parameter_comparison', 'residue_ids_chain_A': residues,
              'exactly_equal': equality, 'all_compared_terms_equal': all(equality.values()),
              'baseline_counts': {k: len(v) for k, v in a.items()},
              'comparison_counts': {k: len(v) for k, v in b.items()},
              'source_sha256': hashes, 'physical_model_validated': False, 'app_ready': False,
              'scope': 'All native atom properties and bonded terms touching the selected residues, including neighboring endpoints. Identical parameters do not establish native coordinates, surrounding nonbonded environment, or a physically accurate loop ensemble.'}
    output.mkdir(exist_ok=False)
    for name, data in [('baseline.json', a), ('comparison.json', b), ('result.json', report)]:
        (output/name).write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['baseline', 'comparison', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--residues', nargs='+', required=True)
    args = parser.parse_args()
    compare(args.baseline, args.comparison, args.residues, args.output)
