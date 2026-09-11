"""Preserve observed metal environments while assembling a nonbonded complex."""
from __future__ import annotations

import numpy as np
from openmm import app, unit

from . import storage


def residue_key(residue):
    return (residue.chain.id, residue.id, (residue.insertionCode or "").strip(), residue.name)


def metal_environment(topology, positions):
    """Inventory close observed donors; distances are evidence, not new bonds."""
    atoms = list(topology.atoms())
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    metals = [a for a in atoms if a.element and a.element.symbol in {"Mg", "Ca", "Zn", "Fe", "Mn", "Cu", "Co", "Ni"}
              and len(list(a.residue.atoms())) == 1]
    waters, protected, contacts = set(), set(), []
    for metal in metals:
        cutoff = .32 if metal.element.symbol == "Ca" else .28
        for atom in atoms:
            if not atom.element or atom.element.symbol not in {"O", "N", "S"}:
                continue
            distance = float(np.linalg.norm(xyz[atom.index] - xyz[metal.index]))
            if .10 < distance <= cutoff:
                key = residue_key(atom.residue)
                if atom.residue.name.upper() in storage.WATERS:
                    waters.add(key)
                else:
                    protected.add(key)
                contacts.append({"metal": list(residue_key(metal.residue)), "donor": [*key, atom.name],
                                 "distance_angstrom": distance * 10})
    return {"water_keys": waters, "protected_keys": protected,
            "report": {"method": "Observed metal–N/O/S contacts within 2.8 Å (3.2 Å for Ca); these distances do not establish bond order or validate a metal model.",
                       "retained_coordinating_waters": [list(key) for key in sorted(waters)],
                       "protected_sidechain_residues": [list(key) for key in sorted(protected)], "contacts": contacts}}


def subset(topology, positions, residue_keys, *, remove_hydrogens=False):
    """Copy selected residues, retaining only their covalent internal connections.

    Monatomic ions use a nonbonded model: CONECT metal contacts are not bonds.
    """
    output = app.Topology()
    mapping, indices = {}, []
    for chain in topology.chains():
        selected = [r for r in chain.residues() if residue_key(r) in residue_keys]
        if not selected:
            continue
        new_chain = output.addChain(chain.id)
        for residue in selected:
            new_residue = output.addResidue(residue.name, new_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                if remove_hydrogens and atom.element == app.element.hydrogen:
                    continue
                mapping[atom] = output.addAtom(atom.name, atom.element, new_residue, atom.id)
                indices.append(atom.index)
    for bond in topology.bonds():
        a, b = bond
        if a not in mapping or b not in mapping:
            continue
        if any(atom.residue.name.upper() in storage.IONS and len(list(atom.residue.atoms())) == 1 for atom in (a, b)):
            continue
        output.addBond(mapping[a], mapping[b], bond.type, bond.order)
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    return app.Modeller(output, xyz[indices] * unit.nanometer)
