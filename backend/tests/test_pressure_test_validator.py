"""Independent edge controls for the offline pressure-test auditor; no MD."""
import importlib.util
from pathlib import Path

import mdtraj as md
import numpy as np

spec = importlib.util.spec_from_file_location("pressure_validator", Path(__file__).resolve().parents[2] / "scripts/validate_pressure_test.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_wrapped_torsions_match_but_large_changes_fail():
    assert audit.series_comparison([180, -179, None], [-180, 181, None], circular=True)["passed"]
    assert not audit.series_comparison([180], [-170], circular=True)["passed"]


def test_null_masks_and_nonfinite_geometry_cannot_pass():
    assert not audit.series_comparison([None], [0])["passed"]
    assert not audit.series_comparison([float("nan")], [float("nan")])["passed"]
    assert not audit.series_comparison([0], [0, 1])["passed"]


def test_protein_bond_screen_uses_periodic_physical_distances():
    topology = md.Topology()
    residue = topology.add_residue("GLY", topology.add_chain())
    a = topology.add_atom("CA", md.element.carbon, residue)
    b = topology.add_atom("C", md.element.carbon, residue)
    topology.add_bond(a, b)
    trajectory = md.Trajectory(np.asarray([[[.95, 0, 0], [.10, 0, 0]]]), topology)
    trajectory.unitcell_vectors = np.eye(3)[None]
    assert audit.protein_bond_screen(trajectory, [0, 1])["passed"]
    trajectory.unitcell_vectors = None
    assert not audit.protein_bond_screen(trajectory, [0, 1])["passed"]


def test_bond_screen_detects_a_single_late_gross_break():
    topology = md.Topology()
    residue = topology.add_residue("GLY", topology.add_chain())
    a = topology.add_atom("CA", md.element.carbon, residue)
    b = topology.add_atom("C", md.element.carbon, residue)
    topology.add_bond(a, b)
    trajectory = md.Trajectory(np.asarray([[[0, 0, 0], [.15, 0, 0]], [[0, 0, 0], [.4, 0, 0]]]), topology)
    result = audit.protein_bond_screen(trajectory, [0, 1])
    assert not result["passed"] and result["maximum_frame"] == 1
