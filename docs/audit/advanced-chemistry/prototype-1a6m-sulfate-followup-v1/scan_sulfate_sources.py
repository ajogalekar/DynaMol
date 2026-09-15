"""Read-only chemical-inventory search of the pinned official CHARMM archives.

This search screen is not a CHARMM parser or a parameter-acceptance engine.
MASS element/mass evidence is used so residue aliases cannot hide sulfate.
Drude subarchives are inspected separately and never called additive models.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile


def documents(archive):
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        return {m.name: tar.extractfile(m).read() for m in tar.getmembers()
                if m.isfile() and m.size < 20_000_000 and
                m.name.lower().endswith(('.rtf', '.str', '.inp', '.prm'))}


def scan(archive, label):
    files = documents(archive)
    masses = defaultdict(set)
    local_masses = {}
    for name, raw in files.items():
        local = defaultdict(set)
        for line in raw.decode(errors='replace').splitlines():
            f = line.split('!')[0].split()
            if not f or f[0].upper() != 'MASS' or len(f) < 4:
                continue
            try:
                mass = float(f[3])
            except ValueError:
                continue
            explicit = f[4].upper() if len(f) > 4 else ''
            element = explicit if explicit in ('H', 'C', 'N', 'O', 'S', 'P', 'F', 'CL', 'NA', 'K', 'CA', 'FE', 'LP') else (
                'S' if 31.5 < mass < 32.5 else 'O' if 15.5 < mass < 16.5 else
                'H' if .9 < mass < 1.2 else 'LP' if mass < .9 else 'OTHER')
            local[f[2].upper()].add(element)
            masses[f[2].upper()].add(element)
        local_masses[name] = local
    candidates, total = [], 0
    for name, raw in files.items():
        lines = raw.decode(errors='replace').splitlines()
        starts = [i for i, line in enumerate(lines) if re.match(r'^\s*(RESI|PRES)\s', line, re.I)]
        for index, first in enumerate(starts):
            total += 1
            last = starts[index+1] if index+1 < len(starts) else len(lines)
            head = lines[first].split('!')[0].split()
            atoms, bonds, flags = [], [], set()
            for line in lines[first+1:last]:
                f = line.split('!')[0].split()
                if not f:
                    continue
                if f[0].upper() in ('END', 'RETURN', 'READ'):
                    break
                if f[0].upper() == 'ATOM' and len(f) >= 4:
                    try:
                        charge = float(f[3])
                    except ValueError:
                        continue
                    elements = local_masses[name].get(f[2].upper(), masses.get(f[2].upper(), set()))
                    atoms.append({'name': f[1], 'type': f[2], 'charge': charge, 'elements': sorted(elements)})
                    if any(word.upper() in ('ALPHA', 'THOLE') for word in f):
                        flags.add('explicit_polarizability')
                if f[0].upper() in ('BOND', 'DOUBLE', 'TRIPLE'):
                    bonds += [f[i:i+2] for i in range(1, len(f)-1, 2)]
                if f[0].upper().startswith('LONE'):
                    flags.add('lonepairs')
            counts = Counter(a['elements'][0] if len(a['elements']) == 1 else 'AMBIGUOUS' for a in atoms)
            if counts['S'] < 1 or counts['O'] < 4:
                continue
            isolated = counts['S'] == 1 and counts['O'] == 4 and not any(
                count for element, count in counts.items() if element not in ('S', 'O', 'LP'))
            row = {'file': name, 'line': first+1, 'kind': head[0], 'residue': head[1],
                   'declared_charge': head[2], 'atom_count': len(atoms), 'element_counts': dict(counts),
                   'flags': sorted(flags), 'sulfate_only_atoms_including_possible_extra_points': isolated}
            if isolated:
                row.update(atoms=atoms, bonds=bonds, body='\n'.join(lines[first:last]),
                           source_sha256=hashlib.sha256(raw).hexdigest())
            candidates.append(row)
    return {'archive': label, 'sha256': hashlib.sha256(archive).hexdigest(), 'text_files_scanned': len(files),
            'residue_patch_blocks_scanned': total, 'sulfur_oxygen_candidates': candidates,
            'scope': 'Candidate search only; alias, extra-point, composition and explicit-polarizability screen. No model acceptance.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    parser.add_argument('archives', type=Path, nargs='+')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Keep prior source audit results; choose a new output')
    results = []
    for path in args.archives:
        raw = path.read_bytes()
        if not raw or len(raw) != path.stat().st_size:
            raise ValueError('Incomplete archive read: '+str(path))
        results.append(scan(raw, str(path)))
        with tarfile.open(fileobj=io.BytesIO(raw)) as tar:
            for member in tar.getmembers():
                if member.isfile() and 'drude/' in member.name.lower() and member.name.endswith('.tgz'):
                    nested = tar.extractfile(member).read()
                    results.append(scan(nested, str(path)+'::'+member.name))
    args.output.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps([{'archive': r['archive'], 'files': r['text_files_scanned'],
                      'blocks': r['residue_patch_blocks_scanned'], 'candidate_residues': [
                          {'file': c['file'], 'residue': c['residue'], 'flags': c['flags']}
                          for c in r['sulfur_oxygen_candidates'] if c['sulfate_only_atoms_including_possible_extra_points']]}
                     for r in results], indent=2))


if __name__ == '__main__':
    main()
