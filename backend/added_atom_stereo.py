"""Correct inverted newly constructed branches without moving observed atoms."""
import numpy as np
from openmm import unit


def correct_added_branches(topology, positions, templates, observed_keys):
    """Reflect only wholly new pendant branches about their attachment plane.

    Standard residue identities define the intended stereoisomer. Reflection
    preserves all bonds within the branch and its attachment length. An
    observed atom, ring, crosslink or planar center cannot be corrected here.
    The caller must still enforce the ordinary geometry and chirality checks.
    """
    from .preparation_worker import atom_key, stereochemistry_report
    from .residue_identity import STANDARD_PROTEINS
    atoms = list(topology.atoms())
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float).copy()
    original = xyz.copy()
    observed = {a.index for a in atoms if atom_key(a) in observed_keys}
    adjacency = [set() for _ in atoms]
    for a, b in topology.bonds():
        adjacency[a.index].add(b.index)
        adjacency[b.index].add(a.index)
    report = {'method': 'Correct only newly built standard sidechain branches by attachment-plane reflection; no observed atom movement or atom renaming.',
              'corrections': [], 'observed_coordinates_preserved_exactly': True}
    for residue in topology.residues():
        if residue.name not in STANDARD_PROTEINS or residue.name not in templates or residue.name in {'GLY', 'PRO'}:
            continue
        members = {a.index for a in residue.atoms()}
        if not members & observed:
            continue
        named = {a.name: a.index for a in residue.atoms()}
        definitions = [('CA', ('N', 'C', 'CB'), ('CB',))]
        if residue.name == 'THR':
            definitions.append(('CB', ('CA', 'OG1', 'CG2'), ('CG2', 'OG1')))
        elif residue.name == 'ILE':
            definitions.append(('CB', ('CA', 'CG1', 'CG2'), ('CG2', 'CG1')))
        template = templates[residue.name]
        ref_named = {a.name: a.index for a in template.topology.atoms()}
        ref = np.asarray(template.positions.value_in_unit(unit.nanometer))
        for center, neighbors, mobile_choices in definitions:
            if not all(n in named and n in ref_named for n in (center, *neighbors)):
                continue
            ci = named[center]
            volume = float(np.linalg.det([xyz[named[n]] - xyz[ci] for n in neighbors]))
            expected = float(np.linalg.det([ref[ref_named[n]] - ref[ref_named[center]] for n in neighbors]))
            if abs(volume) < 1e-4 or volume * expected > 0:
                continue
            for branch_name in mobile_choices:
                branch = named[branch_name]
                moving, pending = set(), [branch]
                while pending:
                    index = pending.pop()
                    if index in moving:
                        continue
                    moving.add(index)
                    pending.extend(j for j in adjacency[index]
                                   if {index, j} != {ci, branch} and j not in moving)
                if ci in moving or not moving.issubset(members) or moving & observed:
                    continue
                plane = [named[n] for n in neighbors if n != branch_name]
                normal = np.cross(xyz[plane[0]] - xyz[ci], xyz[plane[1]] - xyz[ci])
                length = np.linalg.norm(normal)
                if length < 1e-10:
                    continue
                normal /= length
                indices = np.array(sorted(moving), dtype=int)
                displaced = xyz[indices] - xyz[ci]
                candidate = xyz[indices] - 2 * np.outer(displaced @ normal, normal)
                if not np.isfinite(candidate).all():
                    continue
                old = xyz[indices].copy()
                xyz[indices] = candidate
                after = float(np.linalg.det([xyz[named[n]] - xyz[ci] for n in neighbors]))
                if abs(after) < 1e-4 or after * expected <= 0:
                    xyz[indices] = old
                    continue
                report['corrections'].append({'identity': list(atom_key(atoms[ci])[:4]),
                    'center': center, 'moved_atom_names': [atoms[i].name for i in indices],
                    'moved_atom_indices': indices.tolist(), 'before_signed_volume_nm3': volume,
                    'after_signed_volume_nm3': after})
                break
    report['observed_coordinates_preserved_exactly'] = bool(np.array_equal(xyz[sorted(observed)], original[sorted(observed)]))
    if not report['observed_coordinates_preserved_exactly']:
        raise ValueError('Missing-atom chirality correction changed an observed coordinate.')
    report['stereochemistry_after'] = stereochemistry_report(topology, xyz * unit.nanometer, templates)
    return xyz * unit.nanometer, report


def corrected_branch_contacts(topology, positions, corrections):
    """Check corrected heavy atoms against the final retained environment."""
    from scipy.spatial import cKDTree
    from .preparation_worker import atom_key
    atoms = list(topology.atoms())
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
    if xyz.shape != (len(atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError('Corrected-branch geometry needs finite coordinates matching its topology.')
    indices = {}
    for a in atoms:
        indices.setdefault(atom_key(a), []).append(a.index)
    selected = set()
    for record in corrections:
        for name in record['moved_atom_names']:
            key = (*record['identity'], name)
            matches = indices.get(key, [])
            if len(matches) != 1:
                raise ValueError('A corrected branch atom was lost or has ambiguous identity.')
            selected.add(matches[0])
    heavy = [a.index for a in atoms if a.element and a.element.symbol not in {'H', 'D'}]
    selected &= set(heavy)
    adjacency = [set() for _ in atoms]
    for a, b in topology.bonds():
        adjacency[a.index].add(b.index)
        adjacency[b.index].add(a.index)
    tree = cKDTree(xyz[heavy])
    collisions, seen = [], set()
    for i in sorted(selected):
        excluded = {i} | adjacency[i]
        for j in adjacency[i]:
            excluded |= adjacency[j]
        for local in tree.query_ball_point(xyz[i], .13):
            j = heavy[local]
            pair = tuple(sorted((i, j)))
            if j in excluded or pair in seen:
                continue
            seen.add(pair)
            distance = float(np.linalg.norm(xyz[i] - xyz[j]) * 10)
            if distance < 1.3:
                collisions.append({'corrected_atom': list(atom_key(atoms[i])),
                                   'other_atom': list(atom_key(atoms[j])), 'distance_angstrom': distance})
    return {'accepted': not collisions, 'checked_corrected_heavy_atoms': len(selected),
            'minimum_nonbonded_distance_angstrom': 1.3, 'gross_collisions': collisions,
            'scope': 'Corrected new heavy atoms versus the complete retained complex; self, bonds and 1–3 pairs excluded. Unrelated observed contacts are not tested here.'}
