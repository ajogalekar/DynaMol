"""Research helpers for exact observed torsion support during local repair."""
import numpy as np


def residue_key(residue):
    return (residue.chain.id, residue.id, (residue.insertionCode or '').strip(), residue.name)


def peptide_neighbors(topology):
    previous, following = {}, {}
    for a, b in topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        left, right = (a, b) if a.name == 'C' else (b, a)
        if left.residue in following or right.residue in previous:
            raise ValueError('Ambiguous peptide connectivity')
        following[left.residue] = right.residue
        previous[right.residue] = left.residue
    return previous, following


def outer_torsion_anchor_indices(topology, modeled, flanks):
    """Freeze observed atoms supporting torsions outside the mobile region.

    At the left outer join, N/CA of the mobile flank define the preceding
    observed psi/omega. At the right outer join, CA/C define the following
    observed omega/phi. No residue number or target-specific exception is used.
    """
    previous, following = peptide_neighbors(topology)
    selected = set(modeled) | set(flanks)
    indices = set()
    for residue in topology.residues():
        if residue_key(residue) not in flanks:
            continue
        names = {a.name: a.index for a in residue.atoms()}
        required = set()
        if residue in previous and residue_key(previous[residue]) not in selected:
            required |= {'N', 'CA'}
        if residue in following and residue_key(following[residue]) not in selected:
            required |= {'CA', 'C'}
        if not required <= names.keys():
            raise ValueError('Incomplete observed boundary support atoms')
        indices.update(names[name] for name in required)
    return indices


def defining_indices(topology, reference_row):
    """Include preceding omega support, which can select the proline class."""
    atoms = list(topology.atoms())
    indices = set(reference_row['phi_atoms'] + reference_row['psi_atoms'])
    preceding = atoms[reference_row['phi_atoms'][0]].residue
    ca = [a.index for a in preceding.atoms() if a.name == 'CA']
    if len(ca) != 1:
        raise ValueError('Ambiguous preceding alpha-carbon identity')
    return indices | set(ca)


def affected_reference_keys(topology, initial, final, modeled, preliminary):
    """Retain new or actually changed definitions; never exempt changed outliers.

    A boundary depending on a modeled atom is new even if it happens to match
    the archived seed. Unavailable reference rows remain conservative failures.
    """
    atoms = list(topology.atoms())
    result = set(modeled)
    for row in preliminary['rows']:
        indices = sorted(defining_indices(topology, row))
        new_definition = any(residue_key(atoms[i].residue) in modeled for i in indices)
        changed = not np.array_equal(initial[indices], final[indices])
        if new_definition or changed:
            result.add(tuple(row['residue']))
    result.update(tuple(row['residue']) for row in preliminary['unavailable'])
    return result


def outside_observed_torsions(topology, initial, final, modeled, flanks):
    """Verify every fully observed external phi/psi/omega defining atom."""
    previous, following = peptide_neighbors(topology)
    selected = set(modeled) | set(flanks)
    names = {r: {a.name: a for a in r.atoms()} for r in topology.residues()}
    checked, changed = [], []
    for residue in topology.residues():
        if residue_key(residue) in selected:
            continue
        definitions = {}
        if residue in previous:
            before = previous[residue]
            definitions['phi'] = [(before, 'C'), (residue, 'N'), (residue, 'CA'), (residue, 'C')]
            definitions['omega'] = [(before, 'CA'), (before, 'C'), (residue, 'N'), (residue, 'CA')]
        if residue in following:
            definitions['psi'] = [(residue, 'N'), (residue, 'CA'), (residue, 'C'), (following[residue], 'N')]
        for kind, definition in definitions.items():
            if any(residue_key(r) in modeled for r, _ in definition):
                continue
            if any(name not in names[r] for r, name in definition):
                raise ValueError('Incomplete defining atoms for observed external backbone torsion')
            indices = [names[r][name].index for r, name in definition]
            row = {'residue': list(residue_key(residue)), 'torsion': kind, 'indices': indices}
            checked.append(row)
            if not np.array_equal(initial[indices], final[indices]):
                changed.append(row)
    return {'definitions_checked': len(checked), 'changed_definitions': changed,
            'all_defining_coordinates_exactly_preserved': not changed}


class AnchoredConstructionForceField:
    """Only the temporary minimization System receives additional fixed masses."""
    def __init__(self, forcefield, indices):
        self.forcefield = forcefield
        self.indices = frozenset(indices)

    def createSystem(self, *args, **kwargs):
        system = self.forcefield.createSystem(*args, **kwargs)
        for index in self.indices:
            system.setParticleMass(index, 0)
        return system
