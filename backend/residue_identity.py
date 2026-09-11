"""Polymer identity is independent of whether a force-field template is supported.

Modified amino acids remain part of the protein.  This module assigns identity,
never a charge state, replacement residue, or ligand force field.
"""
from __future__ import annotations

import numpy as np

STANDARD_PROTEINS = frozenset({"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"})


def members(owner, name):
    value = getattr(owner, name)
    return list(value() if callable(value) else value)


def residue_key(residue):
    chain = getattr(residue.chain, 'id', None)
    if chain is None:
        chain = residue.chain.chain_id or str(residue.chain.index + 1)
    resid = getattr(residue, 'id', None)
    if resid is None:
        resid = str(residue.resSeq)
    return (chain, str(resid), (getattr(residue, 'insertionCode', '') or '').strip(), residue.name)


def backbone(residue):
    atoms = {atom.name: atom for atom in members(residue, 'atoms')}
    wanted = {'N': 'N', 'CA': 'C', 'C': 'C', 'O': 'O'}
    return atoms if all(name in atoms and atoms[name].element is not None and atoms[name].element.symbol == element for name, element in wanted.items()) else None


def _source_polymer_keys(dataset_id, topology):
    if not dataset_id:
        return set()
    from .ligands import _source_cif
    path = _source_cif(dataset_id)
    if not path:
        return set()
    import gemmi
    block = gemmi.cif.read_file(str(path)).sole_block()
    site = block.get_mmcif_category('_atom_site.')
    components = block.get_mmcif_category('_chem_comp.')
    entities = block.get_mmcif_category('_entity_poly.')
    peptide_components = {name for name, kind in zip(components.get('id', []), components.get('type', [])) if 'PEPTIDE' in kind.upper()}
    peptide_entities = {name for name, kind in zip(entities.get('entity_id', []), entities.get('type', [])) if 'polypeptide' in kind.lower()}
    rows = len(site.get('id', []))
    groups = {}
    for i in range(rows):
        if site.get('pdbx_PDB_model_num', ['1'] * rows)[i] != '1':
            continue
        seq = site.get('label_seq_id', [None] * rows)[i]
        component = site['label_comp_id'][i]
        entity = site.get('label_entity_id', [None] * rows)[i]
        if seq in {None, '', '.', '?', False} or not (component in peptide_components or entity in peptide_entities):
            continue
        if site['type_symbol'][i].upper() in {'H', 'D'}:
            continue
        insertion = site.get('pdbx_PDB_ins_code', [None] * rows)[i] or ''
        insertion = '' if insertion in {'.', '?'} else insertion
        key = (site['label_asym_id'][i], site['auth_asym_id'][i], str(site['auth_seq_id'][i]), insertion, component)
        groups.setdefault(key, set()).add(site['auth_atom_id'][i])
    result = set()
    for residue in members(topology, 'residues'):
        key = residue_key(residue)
        names = {atom.name for atom in members(residue, 'atoms') if atom.element and atom.element.symbol not in {'H', 'D'}}
        candidates = [(source, source_names) for source, source_names in groups.items() if source[2] == key[1] and source[3] == key[2] and source_names.issubset(names)]
        direct = [source for source, _ in candidates if key[0] in source[:2]]
        # Chain renaming is allowed only when the residue/atom identity is unique.
        selected = direct if direct else [source for source, _ in candidates]
        if selected and len({source[4] for source in selected}) == 1:
            result.add(key)
    return result


def protein_residue_keys(dataset_id, topology, positions=None):
    """Identify protein residues using source records and connected peptide backbones."""
    residues = members(topology, 'residues')
    known = {residue_key(residue) for residue in residues if residue.name in STANDARD_PROTEINS or bool(getattr(residue, 'is_protein', False))}
    if dataset_id:
        from . import storage
        try:
            metadata = storage.get_dataset(dataset_id)
        except FileNotFoundError:
            metadata = {}
        categorized = {(atom['chain'], str(atom['resid']), atom['residue']) for atom in metadata.get('atoms', []) if atom['category'] == 'protein'}
        known.update(residue_key(residue) for residue in residues if (residue_key(residue)[0], residue_key(residue)[1], residue.name) in categorized)
        known.update(_source_polymer_keys(dataset_id, topology))
    # A recognized amino-acid name remains polymer even when a few backbone atoms
    # are absent. MDTraj's amino-acid inventory provides identity, not MD support.
    from mdtraj.core.residue_names import _PROTEIN_RESIDUES
    known.update(residue_key(residue) for residue in residues if residue.name in _PROTEIN_RESIDUES)
    candidate = {residue: backbone(residue) for residue in residues}
    candidate = {residue: atoms for residue, atoms in candidate.items() if atoms is not None}
    links = []
    for first, second in members(topology, 'bonds'):
        if first.residue == second.residue or first.residue not in candidate or second.residue not in candidate:
            continue
        if {first.name, second.name} == {'C', 'N'}:
            links.append((residue_key(first.residue), residue_key(second.residue)))
    if positions is not None:
        if hasattr(positions, 'value_in_unit'):
            from openmm import unit
            positions = positions.value_in_unit(unit.nanometer)
        xyz = np.asarray(positions)
        for chain in members(topology, 'chains'):
            possible = [residue for residue in members(chain, 'residues') if residue in candidate]
            for first, second in zip(possible, possible[1:]):
                c, n = candidate[first]['C'], candidate[second]['N']
                separation = np.linalg.norm(xyz[c.index] - xyz[n.index])
                if 0.1 <= separation <= 0.2:
                    links.append((residue_key(first), residue_key(second)))
    # Only extend an established polymer. An isolated amino-acid-shaped ligand
    # with an unknown residue name must not be swept into protein preparation.
    changed = True
    while changed:
        previous = len(known)
        for first, second in links:
            if first in known or second in known:
                known.update((first, second))
        changed = len(known) != previous
    return known
