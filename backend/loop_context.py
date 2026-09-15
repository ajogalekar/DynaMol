"""Bind fresh loop proposals to an explicit complete preparation snapshot.

Admission here produces a seed, never a validated structure. Sequence/length
eligibility, hydrogen rebuilding, full-complex refinement, and all final quality
checks remain separate. No application path enables this module by importing it.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import json
import re
from types import MappingProxyType

import numpy as np
from openmm import unit

from .complex_topology import metal_environment
from .residue_identity import STANDARD_PROTEINS, backbone, residue_key


BACKBONE = frozenset({'N', 'CA', 'C', 'O'})
FINAL_OBSERVED_CAP_NM = 0.1


def _atom_key(atom):
    return (*residue_key(atom.residue), atom.name)


def _hydrogen(atom):
    return atom.element is not None and atom.element.symbol in {'H', 'D'}


def _peptide_bond(a, b):
    # A ligand/cross-chain C--N connection is not a polymer neighbor merely
    # because its atom names happen to match a peptide bond.
    return (a.residue is not b.residue and a.residue.chain is b.residue.chain
            and {a.name, b.name} == {'C', 'N'}
            and all(x.element is not None and x.element.symbol == x.name for x in (a, b))
            and backbone(a.residue) is not None and backbone(b.residue) is not None)


def peptide_neighbors(topology):
    previous, following = {}, {}
    for a, b in topology.bonds():
        if not _peptide_bond(a, b):
            continue
        left, right = (a, b) if a.name == 'C' else (b, a)
        if left.residue in following or right.residue in previous:
            raise ValueError('Ambiguous peptide connectivity')
        following[left.residue] = right.residue
        previous[right.residue] = left.residue
    return previous, following


def outer_torsion_anchor_indices(topology, modeled, flanks):
    """Fix flank atoms supporting fully observed external torsions."""
    previous, following = peptide_neighbors(topology)
    modeled, flanks = set(map(tuple, modeled)), set(map(tuple, flanks))
    selected, indices = modeled | flanks, set()
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
    """Include preceding omega support, which controls the proline class."""
    atoms = list(topology.atoms())
    phi, psi = reference_row['phi_atoms'], reference_row['psi_atoms']
    if len(phi) != 4 or len(psi) != 4 or any(type(i) is not int or not 0 <= i < len(atoms) for i in [*phi, *psi]):
        raise ValueError('Invalid backbone reference atom indices')
    preceding = atoms[phi[0]].residue
    ca = [a.index for a in preceding.atoms() if a.name == 'CA']
    if len(ca) != 1:
        raise ValueError('Ambiguous preceding alpha-carbon identity')
    return set(phi) | set(psi) | set(ca)


def affected_reference_keys(topology, initial, final, modeled, preliminary):
    """Keep modeled or actually changed definitions; retain unavailable rows."""
    atoms = list(topology.atoms())
    modeled = set(map(tuple, modeled))
    result = set(modeled)
    for row in preliminary['rows']:
        indices = sorted(defining_indices(topology, row))
        if (any(residue_key(atoms[i].residue) in modeled for i in indices)
                or not np.array_equal(initial[indices], final[indices])):
            result.add(tuple(row['residue']))
    result.update(tuple(row['residue']) for row in preliminary['unavailable'])
    return result


def outside_observed_torsions(topology, initial, final, modeled, flanks):
    """Verify every fully observed external phi/psi/omega defining atom."""
    previous, following = peptide_neighbors(topology)
    modeled, flanks = set(map(tuple, modeled)), set(map(tuple, flanks))
    selected = modeled | flanks
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
    """Only a fresh temporary construction System receives fixed masses."""
    def __init__(self, forcefield, indices):
        self.forcefield = forcefield
        self.indices = frozenset(indices)

    def createSystem(self, *args, **kwargs):
        system = self.forcefield.createSystem(*args, **kwargs)
        for index in self.indices:
            system.setParticleMass(index, 0)
        return system


@dataclass(frozen=True)
class LoopProposalSeed:
    seeded_xyz_nm: np.ndarray
    modeled_keys: frozenset
    proposed_flank_keys: frozenset
    mobile_flank_keys: frozenset
    extra_fixed_indices: frozenset
    provenance: Mapping


def _identity(value, length):
    if (not isinstance(value, (list, tuple)) or len(value) != length
            or any(not isinstance(v, str) for v in value)
            or any(not value[i] for i in (1, 3))
            or (length == 5 and not value[4])):
        raise ValueError('Invalid exact residue/atom identity')
    if any(v != v.strip() for v in value):
        raise ValueError('Identity fields must already be normalized')
    return tuple(value)


def _keys(values, label):
    rows = [_identity(v, 4) for v in values]
    if len(set(rows)) != len(rows):
        raise ValueError(f'Duplicate {label} residue identity')
    return frozenset(rows)


def _coordinates(value, shape):
    array = np.asarray(value)
    if array.shape != shape or array.dtype.kind not in {'i', 'u', 'f'}:
        raise ValueError('Coordinates must be numeric nanometer arrays of the exact shape')
    array = np.asarray(array, dtype=np.float64)
    if not np.isfinite(array).all():
        raise ValueError('Nonfinite proposal or source coordinates')
    return array.copy()


def _frame(named, xyz):
    origin = xyz[named['CA']]
    x = xyz[named['C']] - origin
    if np.linalg.norm(x) < 1e-12:
        raise ValueError('Degenerate observed residue frame')
    x = x / np.linalg.norm(x)
    y = xyz[named['N']] - origin
    y = y - np.dot(y, x) * x
    if np.linalg.norm(y) < 1e-12:
        raise ValueError('Degenerate observed residue frame')
    y = y / np.linalg.norm(y)
    return origin, np.column_stack((x, y, np.cross(x, y)))


def _scope(topology, modeled, flanks, residues):
    previous, following = peptide_neighbors(topology)
    remaining, component_for, components = set(modeled), {}, []
    while remaining:
        first = min(remaining)
        remaining.remove(first)
        segment, pending = set(), [first]
        while pending:
            key = pending.pop()
            segment.add(key)
            r = residues[key]
            for neighbor in (previous.get(r), following.get(r)):
                if neighbor is not None and residue_key(neighbor) in remaining:
                    nk = residue_key(neighbor)
                    remaining.remove(nk)
                    pending.append(nk)
        number = len(components)
        for key in segment:
            component_for[key] = number
        # The caller's sequence/length policy is separate; each admitted gap
        # still needs two observed stems in the actual connected topology.
        left = [previous.get(residues[k]) for k in segment if previous.get(residues[k]) is None or residue_key(previous[residues[k]]) not in segment]
        right = [following.get(residues[k]) for k in segment if following.get(residues[k]) is None or residue_key(following[residues[k]]) not in segment]
        if len(left) != 1 or len(right) != 1 or left[0] is None or right[0] is None:
            raise ValueError('Loop proposal requires an internal unambiguous peptide gap')
        components.append(segment)
    owners = {}
    for key in flanks:
        r = residues[key]
        owners[key] = {component_for[residue_key(n)] for n in (previous.get(r), following.get(r))
                       if n is not None and residue_key(n) in modeled}
        if not owners[key]:
            raise ValueError('Observed context is not an immediate peptide flank')
        if len(owners[key]) != 1:
            raise ValueError('Conflicting multi-gap observed context')
    for key in flanks:
        n = following.get(residues[key])
        if n is not None and residue_key(n) in flanks and owners[key] != owners[residue_key(n)]:
            raise ValueError('Conflicting adjacent multi-gap observed context')
    return components


def admit_loop_proposal(proposal, topology, source_xyz_nm, modeled_keys, allowed_flanks, *,
                        source_sha256, observed_context_policy='preserve_outer_torsions',
                        protected_residues=()):
    """Validate identity/scope and seed a fresh complete snapshot, without H work.

    ``source_sha256`` is the caller-verified exact native context PDB hash, not a
    hash of the complete preparation snapshot. ``source_xyz_nm`` is the separate
    explicit complete-topology preparation reference (including any documented
    prior sidechain preparation). The caller must bind it to its own snapshot
    hashes and retain the original observed source for final source checks.
    This function never substitutes proposal coordinates as that reference.

    The 1 Å observed-heavy cap is a FINAL refinement gate; raw seeds may exceed
    it. Context sidechains omitted from native rows receive a recorded rigid
    residue-frame initialization. This is not rotamer prediction or validation.
    """
    if (not isinstance(proposal, Mapping) or type(proposal.get('schema_version')) is not int
            or proposal.get('schema_version') != 1 or proposal.get('status') != 'candidate'):
        raise ValueError('Expected a schema_version=1 fresh candidate proposal')
    if not isinstance(source_sha256, str) or not re.fullmatch('[0-9a-f]{64}', source_sha256):
        raise ValueError('Expected an explicit native source SHA256')
    if proposal.get('source_sha256') != source_sha256:
        raise ValueError('Proposal source SHA256 differs from the bound native source')
    if observed_context_policy not in {'fixed_observed', 'mobile_flanks', 'preserve_outer_torsions'}:
        raise ValueError('Unknown observed-context refinement policy')
    modeled = _keys(modeled_keys, 'requested modeled')
    declared = _keys(proposal.get('modeled_residue_keys', []), 'proposal modeled')
    if not modeled or declared != modeled:
        raise ValueError('Proposal must cover the exact complete modeled residue set')
    allowed = _keys(allowed_flanks, 'allowed flank')
    if modeled & allowed:
        raise ValueError('Modeled residues cannot be observed flanks')
    atoms = list(topology.atoms())
    residues_list = list(topology.residues())
    residues = {residue_key(r): r for r in residues_list}
    index = {_atom_key(a): a.index for a in atoms}
    if (len(residues) != len(residues_list) or len(index) != len(atoms)
            or [a.index for a in atoms] != list(range(len(atoms)))):
        raise ValueError('Complete topology contains ambiguous residue/atom identities')
    if not (modeled | allowed) <= residues.keys():
        raise ValueError('Unknown requested modeled or allowed flank residue')
    if any(k[3] not in STANDARD_PROTEINS or backbone(residues[k]) is None for k in modeled):
        raise ValueError('Only complete canonical modeled residues are supported')
    source = _coordinates(source_xyz_nm, (len(atoms), 3))
    updated = source.copy()
    expected = {_atom_key(a) for a in atoms if residue_key(a.residue) in modeled and not _hydrogen(a)}
    rows_by_group, seen = {}, set()
    for group in ('modeled_heavy_atoms', 'observed_context_atoms'):
        rows = proposal.get(group)
        if not isinstance(rows, list):
            raise ValueError(f'Proposal requires explicit {group} rows')
        identities = set()
        for row in rows:
            if not isinstance(row, Mapping) or set(row) != {'identity', 'element', 'xyz_nm'}:
                raise ValueError('Atom rows require exact identity, element, and xyz_nm fields')
            identity = _identity(row['identity'], 5)
            if identity in seen:
                raise ValueError('Duplicate proposal atom identity')
            seen.add(identity)
            if identity not in index:
                raise ValueError('Unknown proposal atom identity')
            atom = atoms[index[identity]]
            if _hydrogen(atom) or row['element'] in {'H', 'D'}:
                raise ValueError('Unexpected hydrogen in heavy-atom proposal')
            if atom.element is None or row['element'] != atom.element.symbol:
                raise ValueError('Proposal atom element differs from complete topology')
            if (group == 'modeled_heavy_atoms') != (identity[:4] in modeled):
                raise ValueError('Modeled and observed atom groups disagree with requested identities')
            updated[atom.index] = _coordinates(row['xyz_nm'], (3,))
            identities.add(identity)
        rows_by_group[group] = identities
    if rows_by_group['modeled_heavy_atoms'] != expected:
        raise ValueError('Proposal must cover every modeled heavy atom exactly')
    flanks = frozenset(k[:4] for k in rows_by_group['observed_context_atoms'])
    if not flanks <= allowed:
        raise ValueError('Proposal changes an unapproved observed flank')
    for key in flanks:
        if key[3] not in STANDARD_PROTEINS or backbone(residues[key]) is None:
            raise ValueError('Unknown or modified observed context is not supported')
        required = BACKBONE | ({'OXT'} if any(a.name == 'OXT' for a in residues[key].atoms()) else set())
        if not all((*key, name) in seen for name in required):
            raise ValueError('Observed context must include complete N/CA/C/O backbone')
    components = _scope(topology, modeled, flanks, residues)
    protected = set(_keys(protected_residues, 'protected'))
    protected.update(metal_environment(topology, source * unit.nanometer)['protected_keys'])
    if flanks & protected:
        raise ValueError('Observed context touches a protected metal-binding residue')
    for a, b in topology.bonds():
        if a.residue is not b.residue and not _peptide_bond(a, b):
            if {residue_key(a.residue), residue_key(b.residue)} & flanks:
                raise ValueError('Observed context touches a chemical crosslink')
    context_heavy = [a.index for a in atoms if residue_key(a.residue) in flanks and not _hydrogen(a)]
    proposed_max = max((float(np.linalg.norm(updated[i]-source[i])) for i in context_heavy), default=0.)
    mobile, fixed = set(flanks), set()
    if observed_context_policy == 'fixed_observed':
        for a in atoms:
            if residue_key(a.residue) not in modeled:
                updated[a.index] = source[a.index]
        mobile.clear()
    elif observed_context_policy == 'preserve_outer_torsions':
        fixed = outer_torsion_anchor_indices(topology, modeled, flanks)
        for i in fixed:
            updated[i] = source[i]
    transferred = []
    for key in sorted(mobile):
        r = residues[key]
        names = {a.name: a.index for a in r.atoms()}
        missing = [a for a in r.atoms() if not _hydrogen(a) and _atom_key(a) not in seen]
        if missing:
            old_origin, old_frame = _frame(names, source)
            new_origin, new_frame = _frame(names, updated)
            for a in missing:
                updated[a.index] = new_origin + new_frame @ old_frame.T @ (source[a.index]-old_origin)
                transferred.append(_atom_key(a))
    if observed_context_policy == 'preserve_outer_torsions':
        if not outside_observed_torsions(topology, source, updated, modeled, mobile)['all_defining_coordinates_exactly_preserved']:
            raise ValueError('External observed torsion support was changed during seeding')
    seeded_max = max((float(np.linalg.norm(updated[i]-source[i])) for i in context_heavy), default=0.)
    provenance = MappingProxyType({
        'schema_version': 1, 'source_sha256': source_sha256,
        'proposal_sha256': hashlib.sha256(json.dumps(proposal, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
        'complete_reference_coordinates_sha256': hashlib.sha256(source.astype('<f8').tobytes()).hexdigest(),
        'reference_kind': 'Explicit caller-supplied complete preparation reference; native source hash is separate',
        'observed_context_policy': observed_context_policy,
        'modeled_components': tuple(tuple(sorted(c)) for c in components),
        'transferred_context_sidechain_atoms': tuple(transferred),
        'additional_fixed_observed_atoms': tuple(_atom_key(atoms[i]) for i in sorted(fixed)),
        'maximum_proposed_context_displacement_nm': proposed_max,
        'maximum_seeded_context_displacement_nm': seeded_max,
        'required_final_observed_heavy_cap_nm': FINAL_OBSERVED_CAP_NM,
        'hydrogen_coordinates_unchanged': True, 'complete_candidate_accepted': False,
    })
    immutable_xyz = np.frombuffer(updated.tobytes(), dtype=np.float64).reshape(updated.shape)
    return LoopProposalSeed(immutable_xyz, modeled, flanks, frozenset(mobile), frozenset(fixed), provenance)
