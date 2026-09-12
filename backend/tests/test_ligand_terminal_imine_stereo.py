"""CCD terminal-imine limitations must not weaken observable E/Z checks."""
import json
from pathlib import Path

import numpy as np
import pytest
from openmm import app, unit
from rdkit import Chem

from backend import ligands


class Reference:
    def __init__(self, molecule, descriptor):
        names = [atom.GetSymbol() + str(atom.GetIdx()) for atom in molecule.GetAtoms()]
        self.atoms = {'atom_id': names, 'type_symbol': [a.GetSymbol() for a in molecule.GetAtoms()],
                      'charge': [str(a.GetFormalCharge()) for a in molecule.GetAtoms()], 'pdbx_stereo_config': ['N'] * len(names)}
        self.bonds = {'atom_id_1': [], 'atom_id_2': [], 'value_order': [], 'pdbx_stereo_config': []}
        for bond in molecule.GetBonds():
            self.bonds['atom_id_1'].append(names[bond.GetBeginAtomIdx()])
            self.bonds['atom_id_2'].append(names[bond.GetEndAtomIdx()])
            self.bonds['value_order'].append('DOUB' if bond.GetBondType() == Chem.BondType.DOUBLE else 'SING')
            self.bonds['pdbx_stereo_config'].append(descriptor if bond.GetBondType() == Chem.BondType.DOUBLE else 'N')

    def get_mmcif_category(self, name):
        return self.atoms if name == '_chem_comp_atom.' else self.bonds


def fixture(smiles, descriptor, coordinates=None):
    molecule = Chem.MolFromSmiles(smiles)
    block = Reference(molecule, descriptor)
    topology = app.Topology()
    residue = topology.addResidue('LIG', topology.addChain('A'), '1')
    for name, element in zip(block.atoms['atom_id'], block.atoms['type_symbol']):
        topology.addAtom(name, app.element.Element.getBySymbol(element), residue)
    xyz = np.asarray(coordinates if coordinates is not None else [[1.3 * i, i % 2, .1 * (i % 3)] for i in range(molecule.GetNumAtoms())], dtype=float)
    return residue, xyz, block


@pytest.mark.parametrize('smiles,descriptor', [('CC=N', 'E'), ('CCC=N', 'Z')])
def test_terminal_imine_descriptor_is_unverified_without_moving_atoms(smiles, descriptor):
    residue, xyz, block = fixture(smiles, descriptor)
    before = xyz.copy()
    molecule, _, _, checked = ligands._graph(residue, xyz, block=block)
    record = json.loads(molecule.GetProp('_DynaMolUnverifiedCCDBondStereo'))[0]
    assert record['status'] == 'unverified' and record['ccd_descriptor'] == descriptor
    assert checked == {}
    assert all(b.GetStereo() == Chem.BondStereo.STEREONONE for b in molecule.GetBonds())
    np.testing.assert_array_equal(molecule.GetConformer().GetPositions(), before)
    np.testing.assert_array_equal(xyz, before)


@pytest.mark.parametrize('smiles', ['CC=CCl', 'CC=NC'])
@pytest.mark.parametrize('descriptor,side,expected_pass', [('E', -1, True), ('Z', -1, False), ('Z', 1, True), ('E', 1, False)])
def test_observable_alkene_and_substituted_imine_keep_strict_stereo(smiles, descriptor, side, expected_pass):
    residue, xyz, block = fixture(smiles, descriptor, [[-.7, 1, 0], [0, 0, 0], [1.3, 0, 0], [2, side, 0]])
    if expected_pass:
        molecule, *_ = ligands._graph(residue, xyz, block=block)
        assert not molecule.HasProp('_DynaMolUnverifiedCCDBondStereo')
    else:
        with pytest.raises(ValueError, match='E/Z'):
            ligands._graph(residue, xyz, block=block)


def test_unresolved_observable_double_bond_and_explicit_opposite_override_still_fail():
    residue, xyz, block = fixture('CC=CCl', 'E', [[i, 0, 0] for i in range(4)])
    with pytest.raises(ValueError, match='E/Z'):
        ligands._graph(residue, xyz, block=block)
    residue, xyz, block = fixture('CC=CCl', 'E', [[-.7, 1, 0], [0, 0, 0], [1.3, 0, 0], [2, -1, 0]])
    molecule, *_ = ligands._graph(residue, xyz, block=block)
    with pytest.raises(ValueError, match='E/Z'):
        ligands._mapped_variant(molecule, 'C/C=C\\Cl')


def test_cached_benzamidine_preserves_pose_and_selects_protonated_ph7_state():
    root = Path(__file__).resolve().parents[2] / 'data/integration-checks/five-system-openmm'
    ccd = root / 'chemistry/ccd/BEN.cif'
    pdb = root / 'datasets/59ff9ca69b454df1/topology.pdb'
    if not ccd.exists() or not pdb.exists():
        pytest.skip('Optional retained 3PTB/BEN audit fixture is unavailable; generic controls remain active')
    import gemmi
    structure = app.PDBFile(str(pdb))
    residue = next(r for r in structure.topology.residues() if r.name == 'BEN')
    xyz = np.asarray(structure.positions.value_in_unit(unit.angstrom))
    molecule, atoms, names, checked = ligands._graph(residue, xyz, block=gemmi.cif.read_file(str(ccd)).sole_block())
    selected, state = ligands._select_state(molecule, names, 'BEN', 7, None)
    assert state['original_formal_charge'] == 0 and state['formal_charge'] == 1
    assert checked == {} and molecule.HasProp('_DynaMolUnverifiedCCDBondStereo')
    np.testing.assert_array_equal(selected.GetConformer().GetPositions(), xyz[[atom.index for atom in atoms]])
