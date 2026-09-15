"""Preserve native numeric precision in selected ParmEd text-export fields.

Only round-trippable native values are written. Ambiguous atom mappings reject;
the result still requires independent complete energy/force comparison.
"""
import json

import parmed


def write_precise_topology(native, output):
    model = parmed.gromacs.GromacsTopologyFile.from_structure(native)
    raw = output.with_name(output.stem+'-parmed.top')
    model.write(str(raw))
    atomtypes = {}
    charges = {}
    for atom in model.atoms:
        value = (atom.sigma*.1, atom.epsilon*4.184)
        if atom.type in atomtypes and atomtypes[atom.type] != value:
            raise ValueError('Ambiguous native Lennard-Jones type')
        atomtypes[atom.type] = value
        key = (atom.type, atom.residue.name, atom.name, format(atom.charge, '.8f'))
        value = (atom.charge, atom.mass)
        if key in charges and max(abs(x-y) for x, y in zip(charges[key], value)) > 1e-12:
            raise ValueError('Ambiguous native charge/mass in a repeated molecule template')
        charges[key] = value
    section = None
    lines = []
    counts = {'defaults': 0, 'atomtypes': 0, 'atoms': 0}
    for line in raw.read_text().splitlines():
        text = line.split(';', 1)[0].strip()
        if text.startswith('['):
            section = text.strip('[] ')
        parts = text.split()
        changed = False
        if parts and not text.startswith(('[', '#')):
            if section == 'defaults':
                if len(parts) != 5:
                    raise ValueError('Unexpected GROMACS defaults layout')
                parts[3:] = [format(model.defaults.fudgeLJ, '.17g'), format(model.defaults.fudgeQQ, '.17g')]
                changed = True
            elif section == 'atomtypes':
                if parts[0] not in atomtypes:
                    raise ValueError('Exported atom type is not present in the native model')
                parts[-2:] = [format(x, '.17g') for x in atomtypes[parts[0]]]
                changed = True
            elif section == 'atoms':
                if len(parts) != 8:
                    raise ValueError('Alchemical or unexpected atom layout requires a separate exporter')
                key = (parts[1], parts[3], parts[4], format(float(parts[6]), '.8f'))
                if key not in charges:
                    raise ValueError('Exported atom has no unambiguous native charge/mass match')
                parts[6:] = [format(x, '.17g') for x in charges[key]]
                changed = True
        if changed:
            counts[section] += 1
            line = ' '.join(parts)
        lines.append(line)
    if counts['defaults'] != 1 or counts['atomtypes'] != len(atomtypes) or not counts['atoms']:
        raise ValueError('Incomplete precision-preserving topology export')
    output.write_text('\n'.join(lines)+'\n')
    report = {'fields_restored_from_native': counts, 'parameter_fitting_or_rescaling': False,
        'original_parmed_export': str(raw),
        'scope': 'Preserve native LJ, charge, mass and global 1-4 scaling values at 17 significant digits. Bonded fields remain ParmEd output and require complete energy/force validation.'}
    output.with_name(output.stem+'-precision.json').write_text(json.dumps(report, indent=2)+'\n')
    return report
