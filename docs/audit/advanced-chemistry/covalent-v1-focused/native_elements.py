"""Bind scoped native residue atom elements to an explicit source topology."""
import re

from openmm import app


def inventory(prmtop):
    topology = app.AmberPrmtopFile(str(prmtop)).topology
    result = []
    for atom in topology.atoms():
        if atom.element is None:
            raise ValueError(f'Unspecified element/extra point in all-atom research model: {atom.residue.name}:{atom.name}')
        result.append((atom.name, atom.element.atomic_number))
    return result


def write_unit_elements(source, unit_name, output):
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', unit_name):
        raise ValueError('Unsafe native unit identifier')
    topology = app.AmberPrmtopFile(str(source)).topology
    if topology.getNumResidues() != 1:
        raise ValueError('Element declarations require a single-residue source unit')
    expected = inventory(source)
    # Native integer atom selectors avoid interpreting * in Amber sugar names.
    lines = [f'set {unit_name}.1.{i+1} element "{atom.element.symbol}"'
             for i, atom in enumerate(topology.atoms())]
    output.write_text('\n'.join(lines)+'\n')
    return expected
