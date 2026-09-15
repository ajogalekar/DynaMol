"""Provisional loop joins may refine; observed breaks and bad final geometry may not."""
import numpy as np
from openmm import app, unit

from backend import config
from backend.loop_geometry import loop_geometry_report
from backend.preparation_worker import construction_backbone_gaps
from backend.residue_identity import residue_key


def test_only_requested_loop_joins_defer_and_final_geometry_still_rejects():
    pdb = app.PDBFile(str(config.ROOT / 'docs/audit/preparation-fixtures/six_residues_intact.pdb'))
    top = pdb.topology
    xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer)).copy()
    residues = list(top.residues())
    modeled = {residue_key(residues[2])}
    for atom in residues[2].atoms():
        xyz[atom.index] += [1, 0, 0]
    for atom in residues[5].atoms():
        xyz[atom.index] += [0, 1, 0]
    before = xyz.copy()
    fixed, pending = construction_backbone_gaps(top, xyz * unit.nanometer, modeled)
    assert {(g['after_index'], g['before_index']) for g in pending} == {(1, 2), (2, 3)}
    assert {(g['after_index'], g['before_index']) for g in fixed} == {(4, 5)}
    assert not loop_geometry_report(top, xyz * unit.nanometer, modeled)['accepted']
    np.testing.assert_array_equal(xyz, before)
    fixed, pending = construction_backbone_gaps(top, xyz * unit.nanometer, set())
    assert len(fixed) == 3 and not pending
