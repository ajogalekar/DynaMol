"""Explicit source/display identity maps for native parameterized Systems.

Native templates can rename residues or combine a protein residue and ligand.
Display topology is therefore separate and coordinates are permuted explicitly.
No atom, bond, parameter, chemical state or coordinate value is invented here.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json

import numpy as np


FIELDS = ("chain", "resid", "insertion", "resname", "atomname")


def _source_identity(record):
    identity = record.get("source_identity")
    legacy = record.get("source_or_added_id")
    legacy_identity = None
    if legacy is not None:
        if not isinstance(legacy, str) or len(legacy.split(":")) != 5:
            raise ValueError("Every native atom needs an explicit source or added-atom identity")
        legacy_identity = dict(zip(FIELDS, legacy.split(":")))
    if identity is None:
        identity = legacy_identity
    elif legacy_identity is not None and identity != legacy_identity:
        raise ValueError("Explicit and legacy source atom identities conflict")
    if (not isinstance(identity, dict) or set(identity) != set(FIELDS)
            or any(not isinstance(identity[key], str) for key in FIELDS)
            or any(not identity[key] for key in ("resid", "resname", "atomname"))):
        raise ValueError("Atom identities require chain, resid, insertion, resname and atomname strings")
    return tuple(identity[key] for key in FIELDS)


def _indices(values, count):
    if not isinstance(values, (list, tuple, np.ndarray)):
        raise ValueError("Atom selections must be a list of integer indices")
    if isinstance(values,np.ndarray) and values.ndim!=1:
        raise ValueError("Atom selection indices must be a one-dimensional list")
    if any(isinstance(i, (bool, np.bool_)) or not isinstance(i, (int, np.integer)) or not 0 <= i < count for i in values):
        raise ValueError("Atom selection index is outside the mapped particle inventory")
    if not len(values) or len(set(values))!=len(values):
        raise ValueError("Measurement atom indices must be nonempty and unique")
    return np.asarray(values, dtype=np.int64)


@dataclass(frozen=True)
class NativeAtomMap:
    display_topology: object
    display_to_native: tuple[int, ...]
    native_to_display: tuple[int, ...]
    source_identities_native_order: tuple[tuple[str, ...], ...]
    mapping_sha256: str

    def _coordinates(self, coordinates_nm, order):
        coords = np.asarray(coordinates_nm)
        count = len(order)
        if (coords.ndim not in (2, 3) or coords.shape[-2:] != (count, 3)
                or not np.issubdtype(coords.dtype, np.number) or not np.isrealobj(coords) or not np.isfinite(coords).all()):
            raise ValueError("Coordinates must be finite N×3 or F×N×3 arrays in nm for the complete mapped inventory")
        return np.take(coords, order, axis=-2)

    def to_display(self, coordinates_nm):
        return self._coordinates(coordinates_nm, self.display_to_native)

    def to_native(self, coordinates_nm):
        return self._coordinates(coordinates_nm, self.native_to_display)

    def measurement_to_native(self, display_indices, *, expected_mapping_sha256):
        """Map a selection pinned when picked/saved, never stale bare indices.

        The caller must persist the expected hash with the selection. Supplying
        the current hash for a previously saved unbound selection is not a
        valid recovery strategy: it defeats the identity check.
        """
        if expected_mapping_sha256!=self.mapping_sha256:
            raise ValueError("Measurement mapping checksum differs; selection belongs to another atom mapping")
        selected = _indices(display_indices, len(self.display_to_native))
        return np.asarray(self.display_to_native, dtype=np.int64)[selected].tolist()

    def manifest(self):
        return {"schema_version": 1, "mapping_sha256": self.mapping_sha256,
                "display_to_native": list(self.display_to_native),
                "native_to_display": list(self.native_to_display),
                "source_identities_native_order": [dict(zip(FIELDS, row)) for row in self.source_identities_native_order],
                "scope": "Coordinate/index permutation and display labels only; native System parameters remain separate."}


def build_native_atom_map(native_topology, records):
    """Build a bijection from an explicit complete per-native-particle map.

    The caller binds this map to its exact native topology/parameter snapshot.
    Every generated atom needs an explicit identity as well as observed atoms.
    Grouping source chains/residues may change display atom order; the returned
    permutation must accompany every trajectory frame and selected measurement.
    """
    from openmm import app
    atoms = list(native_topology.atoms())
    if not isinstance(records, list) or len(records) != len(atoms):
        raise ValueError("Source identity map must contain every native atom exactly once")
    if ([row.get("native_index") for row in records if isinstance(row, dict)] != list(range(len(atoms)))
            or any(type(row.get("native_index")) is not int for row in records)):
        raise ValueError("Source identity records must explicitly follow complete native particle order")
    identities = tuple(_source_identity(row) for row in records)
    if len(set(identities)) != len(identities):
        raise ValueError("Source atom identities are duplicated; refusing an ambiguous mapping")
    native_identities=tuple((str(a.residue.chain.id),str(a.residue.id),str(a.residue.insertionCode or ''),
                             str(a.residue.name),str(a.name)) for a in atoms)
    for index,record in enumerate(records):
        if 'native_identity' in record:
            declared=record['native_identity']
            if (not isinstance(declared,dict) or set(declared)!=set(FIELDS) or
                tuple(declared[k] for k in FIELDS)!=native_identities[index]):
                raise ValueError("Explicit native identity does not match the stated native particle index")
    # Preserve first-seen chain and residue order, keeping all atoms of a source
    # residue together even when a native adduct template interleaves them.
    chains = {}
    for index, identity in enumerate(identities):
        residues = chains.setdefault(identity[0], {})
        residues.setdefault(identity[1:4], []).append(index)
    display = app.Topology()
    display_to_native, imported = [], {}
    for chain_name, residues in chains.items():
        chain = display.addChain(chain_name)
        for (resid, insertion, resname), indices in residues.items():
            residue = display.addResidue(resname, chain, id=resid, insertionCode=insertion)
            for native_index in indices:
                atom = atoms[native_index]
                imported[native_index] = display.addAtom(identities[native_index][4], atom.element, residue,
                                                        id=str(len(display_to_native) + 1))
                display_to_native.append(native_index)
    bonds = []
    for bond in native_topology.bonds():
        a, b = bond[0].index, bond[1].index
        display.addBond(imported[a], imported[b], type=bond.type, order=bond.order)
        bonds.append(sorted((a, b)))
    display.setPeriodicBoxVectors(native_topology.getPeriodicBoxVectors())
    inverse = np.argsort(display_to_native).tolist()
    payload = {"source_identities_native_order": identities,
               "native_identities_native_order": native_identities,
               "display_to_native": display_to_native,
               "native_elements": [a.element.symbol if a.element else None for a in atoms],
               "native_bonds": sorted(bonds)}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return NativeAtomMap(display, tuple(display_to_native), tuple(inverse), identities, digest)
