"""Actual CPU OpenMM minimization on a tiny, analytically known potential."""
import json

import numpy as np
import openmm as mm
import pytest
from openmm import app, unit

from backend import config
from backend.loop_refinement import minimize_modeled_loops
from backend.residue_identity import residue_key


class AnalyticalForceField:
    """Fresh temporary Systems; no protein template/native model is involved."""

    def __init__(self):
        self.calls = []
        self.systems = []

    def createSystem(self, topology, **options):
        self.calls.append(options)
        assert options == {"nonbondedMethod": app.NoCutoff, "constraints": None, "rigidWater": False}
        system = mm.System()
        for atom in topology.atoms():
            system.addParticle(atom.element.mass)
        bonds = mm.HarmonicBondForce()
        # Independent mobile heavy atom and H attached to the fixed reference.
        bonds.addBond(0, 1, .10, 20000)
        bonds.addBond(0, 2, .15, 10000)
        system.addForce(bonds)
        # A substantial force on the observed atom demonstrates that zero mass
        # freezes coordinates during minimization, not just that force is zero.
        reference = mm.CustomExternalForce("1000*(x*x+y*y+z*z)")
        reference.addParticle(0, [])
        system.addForce(reference)
        self.systems.append(system)
        return system


class BackboneForceField:
    """Bond-only analytical model for checking actual temporary torsion scope.

    Residue labels exercise target-selection policy, not force-field chemistry.
    """

    def __init__(self, xyz):
        self.xyz = xyz.copy()
        self.systems = []

    def createSystem(self, topology, **options):
        assert options["constraints"] is None
        system = mm.System()
        for atom in topology.atoms():
            system.addParticle(atom.element.mass)
        force = mm.HarmonicBondForce()
        for a, b in topology.bonds():
            force.addBond(a.index, b.index, float(np.linalg.norm(self.xyz[a.index] - self.xyz[b.index])), 10000)
        system.addForce(force)
        self.systems.append(system)
        return system


@pytest.fixture
def analytical_system(monkeypatch):
    monkeypatch.setattr(config, "CPU_THREADS", 2)
    topology = app.Topology()
    chain = topology.addChain("A")
    observed = topology.addResidue("ALA", chain, "1")
    fixed = topology.addAtom("C", app.element.carbon, observed)
    hydrogen = topology.addAtom("H", app.element.hydrogen, observed)
    modeled = topology.addResidue("ALA", chain, "2", "A")
    mobile = topology.addAtom("CA", app.element.carbon, modeled)
    topology.addBond(fixed, mobile)
    topology.addBond(fixed, hydrogen)
    center = np.array([1.234567891, .234567891, -.345678912])
    xyz = np.array([center, center + [0, .30, 0], center + [.45, 0, 0]])
    return topology, xyz, {residue_key(modeled)}, AnalyticalForceField()


def test_actual_cpu_minimization_freezes_observed_atom_with_nonzero_force(analytical_system):
    topology, xyz, keys, forcefield = analytical_system
    original = xyz.copy()
    result, report = minimize_modeled_loops(topology, xyz * unit.nanometer, forcefield, keys)
    refined = result.value_in_unit(unit.nanometer)
    np.testing.assert_array_equal(refined[0], original[0])
    np.testing.assert_array_equal(xyz, original)
    assert np.linalg.norm(refined[1] - refined[0]) == pytest.approx(.10, abs=.001)
    assert np.linalg.norm(refined[2] - refined[0]) == pytest.approx(.15, abs=.002)
    assert not np.array_equal(refined[1:], original[1:])
    assert report["fixed_heavy_coordinates_preserved_exactly"]
    assert report["fixed_heavy_atoms"] == 1
    assert report["mobile_atoms"] == 2
    assert report["energy_after_kj_mol"] < report["energy_before_kj_mol"] - 100
    assert report["mobile_force_rms_kj_mol_nm"] <= 10
    system = forcefield.systems[0]
    assert system.getNumConstraints() == 0
    assert system.getParticleMass(0).value_in_unit(unit.dalton) == 0
    assert all(system.getParticleMass(i).value_in_unit(unit.dalton) > 0 for i in (1, 2))
    assert len(forcefield.calls) == 1
    assert report["max_iterations"] == 1000
    json.dumps(report, allow_nan=False)


def test_iteration_limited_result_reports_residual_force_without_convergence_claim(analytical_system):
    topology, xyz, keys, forcefield = analytical_system
    result, report = minimize_modeled_loops(topology, xyz * unit.nanometer, forcefield, keys, max_iterations=1)
    np.testing.assert_array_equal(result.value_in_unit(unit.nanometer)[0], xyz[0])
    assert report["max_iterations"] == 1
    assert np.isfinite(report["mobile_force_rms_kj_mol_nm"])
    assert report["energy_after_kj_mol"] <= report["energy_before_kj_mol"]
    assert "converged" not in report
    assert "no convergence" in report["scope"]


@pytest.mark.parametrize("iterations", [0, -1, 1001])
def test_outside_iteration_budget_rejected_before_native_context(analytical_system, iterations):
    topology, xyz, keys, forcefield = analytical_system
    with pytest.raises(ValueError, match="1–1000"):
        minimize_modeled_loops(topology, xyz * unit.nanometer, forcefield, keys, max_iterations=iterations)
    assert forcefield.calls == []


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_invalid_input_coordinates_rejected_before_parameter_creation(analytical_system, value):
    topology, xyz, keys, forcefield = analytical_system
    xyz[1, 0] = value
    with pytest.raises(ValueError, match="finite coordinates"):
        minimize_modeled_loops(topology, xyz * unit.nanometer, forcefield, keys)
    assert forcefield.calls == []


def test_exact_insertion_code_is_required_before_native_context(analytical_system):
    topology, xyz, _, forcefield = analytical_system
    with pytest.raises(ValueError, match="exact modeled residue identities"):
        minimize_modeled_loops(topology, xyz * unit.nanometer, forcefield, {("A", "2", "", "ALA")})
    assert forcefield.calls == []


def test_duplicate_exact_modeled_identity_is_rejected(analytical_system):
    topology, xyz, keys, forcefield = analytical_system
    duplicate = topology.addResidue("ALA", topology.addChain("A"), "2", "A")
    topology.addAtom("CA", app.element.carbon, duplicate)
    xyz = np.concatenate((xyz, [[10, 10, 10]]))
    with pytest.raises(ValueError, match="unique|exact modeled residue identities"):
        minimize_modeled_loops(topology, xyz * unit.nanometer, forcefield, keys)
    assert forcefield.calls == []


def test_empty_modeled_residue_is_rejected_before_native_context():
    topology = app.Topology()
    residue = topology.addResidue("GLY", topology.addChain("A"), "2")
    forcefield = AnalyticalForceField()
    with pytest.raises(ValueError, match="atom|exact modeled residue identities"):
        minimize_modeled_loops(topology, np.empty((0, 3)) * unit.nanometer, forcefield, {residue_key(residue)})
    assert forcefield.calls == []


def peptide_fixture():
    pdb = app.PDBFile(str(config.ROOT / 'docs/audit/preparation-fixtures/six_residues_intact.pdb'))
    return pdb.topology, np.asarray(pdb.positions.value_in_unit(unit.nanometer)).copy()


def test_construction_torsions_only_touch_modeled_residues_and_stay_out_of_fresh_system(monkeypatch):
    monkeypatch.setattr(config, "CPU_THREADS", 2)
    topology, xyz = peptide_fixture()
    residues = list(topology.residues())
    keys = {residue_key(residues[2])}
    forcefield = BackboneForceField(xyz)
    refined, report = minimize_modeled_loops(topology, xyz * unit.nanometer, forcefield, keys)
    temporary = forcefield.systems[0]
    torsions = [force for force in temporary.getForces() if isinstance(force, mm.CustomTorsionForce)]
    assert len(torsions) == 1
    assert torsions[0].getNumTorsions() == 2
    atoms = list(topology.atoms())
    actual_pairs = set()
    for index in range(torsions[0].getNumTorsions()):
        a, b, c, d, parameters = torsions[0].getTorsionParameters(index)
        assert [atoms[i].name for i in (a, b, c, d)] == ["CA", "C", "N", "CA"]
        actual_pairs.add((atoms[b].residue.index, atoms[c].residue.index))
        assert parameters[0] == pytest.approx(200)
        assert parameters[1] == pytest.approx(np.pi)
    # Both anchors are included; the three observed–observed links are absent.
    assert actual_pairs == {(1, 2), (2, 3)}
    assert report["energy_includes_construction_restraints"]
    assert report["peptide_construction_restraints"]["exported_to_simulation"] is False
    targets = report["peptide_construction_restraints"]["targets"]
    assert len(targets) == 2 and all(np.isfinite(row["initial_omega_degrees"]) for row in targets)
    fixed = [a.index for a in atoms if a.residue != residues[2]]
    np.testing.assert_array_equal(refined.value_in_unit(unit.nanometer)[fixed], xyz[fixed])
    # The factory is reusable and is not a container for the added construction force.
    fresh = forcefield.createSystem(topology, constraints=None)
    assert fresh.getNumForces() == 1
    assert not any(isinstance(force, mm.CustomTorsionForce) for force in fresh.getForces())


@pytest.mark.parametrize("name,cis,target", [("PRO", True, 0), ("PRO", False, 180),
                                          ("HYP", True, 0), ("PHE", True, 180)])
def test_target_basin_policy_is_recorded_for_new_link_to_observed_neighbor(monkeypatch, name, cis, target):
    monkeypatch.setattr(config, "CPU_THREADS", 2)
    topology, xyz = peptide_fixture()
    residues = list(topology.residues())
    residues[3].name = name
    if cis:
        before = {a.name: a.index for a in residues[2].atoms()}
        after = {a.name: a.index for a in residues[3].atoms()}
        pivot = xyz[after['N']].copy()
        axis = pivot - xyz[before['C']]
        axis /= np.linalg.norm(axis)
        downstream = [a.index for a in topology.atoms() if a.residue.index >= 3]
        vectors = xyz[downstream] - pivot
        xyz[downstream] = pivot - vectors + 2 * np.outer(vectors @ axis, axis)
    _, report = minimize_modeled_loops(topology, xyz * unit.nanometer, BackboneForceField(xyz),
                                      {residue_key(residues[2])})
    row = next(row for row in report['peptide_construction_restraints']['targets']
               if row['modeled_bond_before_residue'][1] == '4')
    assert row['target_degrees'] == target
    assert (abs(row['initial_omega_degrees']) < 90) is cis
