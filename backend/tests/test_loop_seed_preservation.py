"""Exact target contracts and static analytical OpenMM force checks; no sampling."""
import copy
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit
import pytest

from backend import config
from backend.loop_geometry import _dihedral
from backend.loop_refinement import _minimize_attempt, _validated_seed_torsion_targets
from backend.residue_identity import residue_key


class BondOnlyFactory:
    """Fresh analytical systems, deliberately without native protein parameters."""

    def __init__(self, xyz):
        self.xyz = xyz.copy()
        self.systems = []

    def createSystem(self, topology, **options):
        system = mm.System()
        for atom in topology.atoms():
            system.addParticle(atom.element.mass)
        bonds = mm.HarmonicBondForce()
        for a, b in topology.bonds():
            bonds.addBond(a.index, b.index, float(np.linalg.norm(self.xyz[a.index] - self.xyz[b.index])), 1000.0)
        system.addForce(bonds)
        self.systems.append(system)
        return system


def records(topology, xyz, selected):
    residues = list(topology.residues())
    names = [{a.name: a.index for a in r.atoms()} for r in residues]
    targets = []
    for i, residue in enumerate(residues):
        if residue_key(residue) not in selected:
            continue
        definitions = {'phi': [names[i-1]['C'], names[i]['N'], names[i]['CA'], names[i]['C']],
                       'psi': [names[i]['N'], names[i]['CA'], names[i]['C'], names[i+1]['N']]}
        for kind, indices in definitions.items():
            targets.append({'residue': list(residue_key(residue)), 'kind': kind, 'atoms': indices,
                            'target_radians': float(np.deg2rad(_dihedral(xyz[indices])))})
    return targets


@pytest.fixture
def sample(monkeypatch):
    monkeypatch.setattr(config, 'CPU_THREADS', 2)
    # Existing coordinates are only a convenient nondegenerate analytical test
    # geometry. There is no native model generation or conformational sampling.
    pdb = app.PDBFile(str(Path(config.ROOT) / 'docs/audit/preparation-fixtures/six_residues_intact.pdb'))
    top = pdb.topology
    xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer)).copy()
    residues = list(top.residues())
    residues[2].insertionCode = 'B'
    modeled = {residue_key(residues[2])}
    return top, xyz, residues, modeled, records(top, xyz, modeled)


def validate(sample, targets=None, xyz=None, flanks=()):
    top, original, residues, modeled, source = sample
    return _validated_seed_torsion_targets(top, original if xyz is None else xyz,
                                           modeled, set(flanks), source if targets is None else targets)


def test_exact_target_scope_and_input_records_are_preserved(sample):
    top, xyz, residues, modeled, _ = sample
    flanks = {residue_key(residues[1]), residue_key(residues[3])}
    targets = records(top, xyz, modeled | flanks)
    original = copy.deepcopy(targets)
    checked = validate(sample, targets, flanks=flanks)
    assert targets == original
    assert len(checked) == 6
    assert {tuple(row['residue']) for row in checked} == modeled | flanks
    assert {row['kind'] for row in checked} == {'phi', 'psi'}
    assert next(row for row in checked if row['residue'][1] == '3')['residue'][2] == 'B'
    assert all(np.deg2rad(row['initial_degrees']) == pytest.approx(row['target_radians']) for row in checked)


@pytest.mark.parametrize('value', [True, '0.0', float('nan'), float('inf'), -float('inf'), 4.0])
def test_target_requires_finite_canonical_numeric_radians(sample, value):
    targets = copy.deepcopy(sample[-1]); targets[0]['target_radians'] = value
    with pytest.raises(ValueError, match='finite canonical radians'):
        validate(sample, targets)


@pytest.mark.parametrize('value', [True, 1.0, -1, 100000])
def test_indices_require_exact_in_range_integers(sample, value):
    targets = copy.deepcopy(sample[-1]); targets[0]['atoms'][0] = value
    with pytest.raises(ValueError, match='integer atom indices'):
        validate(sample, targets)


def test_wrong_backbone_quadruple_and_insertion_code_are_rejected(sample):
    targets = copy.deepcopy(sample[-1]); targets[0]['atoms'][0] += 1
    with pytest.raises(ValueError, match='exact bonded peptide'):
        validate(sample, targets)
    targets = copy.deepcopy(sample[-1]); targets[0]['residue'][2] = ''
    with pytest.raises(ValueError, match='outside.*scope'):
        validate(sample, targets)


@pytest.mark.parametrize('kind', ['omega', 'chi1', None, []])
def test_only_phi_and_psi_are_supported(sample, kind):
    targets = copy.deepcopy(sample[-1]); targets[0]['kind'] = kind
    with pytest.raises(ValueError, match='outside.*scope'):
        validate(sample, targets)


def test_duplicate_missing_and_extra_target_records_are_rejected(sample):
    targets = sample[-1]
    for incomplete in [[], targets[:1], targets + targets]:
        with pytest.raises(ValueError, match='exactly phi and psi'):
            validate(sample, incomplete)
    with pytest.raises(ValueError, match='Duplicate'):
        validate(sample, [targets[0], copy.deepcopy(targets[0])])
    extra = copy.deepcopy(targets); extra[0]['strength_kj_mol'] = 1.0
    with pytest.raises(ValueError, match='target fields'):
        validate(sample, extra)


def test_targets_cannot_select_an_observed_residue_outside_allowed_flanks(sample):
    top, xyz, residues, _, _ = sample
    targets = records(top, xyz, {residue_key(residues[3])})
    with pytest.raises(ValueError, match='outside.*scope'):
        validate(sample, targets)


def test_every_allowed_flank_also_needs_both_qualified_targets(sample):
    flanks = {residue_key(sample[2][1])}
    with pytest.raises(ValueError, match='exactly phi and psi'):
        validate(sample, flanks=flanks)


def test_target_must_match_actual_seed_angle_before_system_creation(sample):
    top, xyz, _, modeled, targets = sample
    targets = copy.deepcopy(targets); targets[0]['target_radians'] += .01
    factory = BondOnlyFactory(xyz)
    with pytest.raises(ValueError, match='differs from.*initial'):
        _minimize_attempt(top, xyz * unit.nanometer, factory, modeled, seed_torsion_targets=targets)
    assert factory.systems == []


def test_equivalent_signed_pi_targets_match_circularly(sample):
    top, xyz, _, modeled, source = sample
    xyz = xyz.copy()
    xyz[source[0]['atoms']] = [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, -1, 0]]
    targets = records(top, xyz, modeled)
    assert abs(targets[0]['target_radians']) == pytest.approx(np.pi)
    targets[0]['target_radians'] *= -1
    assert validate(sample, targets, xyz=xyz)[0]['target_radians'] == targets[0]['target_radians']


def test_degenerate_seed_angle_and_nonfinite_coordinates_are_rejected(sample):
    xyz = sample[1].copy(); indices = sample[-1][0]['atoms']
    xyz[indices[2]] = xyz[indices[1]]
    with pytest.raises(ValueError, match='degenerate'):
        validate(sample, xyz=xyz)
    xyz[indices[2], 0] = np.nan
    with pytest.raises(ValueError, match='finite coordinates'):
        validate(sample, xyz=xyz)


def test_missing_peptide_connection_is_not_inferred_from_residue_number(sample):
    top, _, residues, _, _ = sample
    previous = next(a for a in residues[1].atoms() if a.name == 'C')
    current = next(a for a in residues[2].atoms() if a.name == 'N')
    top._bonds = [bond for bond in top.bonds() if set(bond) != {previous, current}]
    with pytest.raises(ValueError, match='complete peptide neighborhoods'):
        validate(sample)


def test_missing_intrarezidue_bond_invalidates_torsion_definition(sample):
    top, _, residues, _, _ = sample
    names = {a.name: a for a in residues[2].atoms()}
    top._bonds = [bond for bond in top.bonds() if set(bond) != {names['N'], names['CA']}]
    with pytest.raises(ValueError, match='exact bonded peptide'):
        validate(sample)


def test_ligand_crosslink_cannot_enter_seed_preservation(sample):
    top, xyz, residues, modeled, targets = sample
    ligand = top.addResidue('LIG', top.addChain('L'), '1')
    atom = top.addAtom('C', app.element.carbon, ligand)
    top.addBond(next(a for a in residues[2].atoms() if a.name == 'CB'), atom)
    xyz = np.concatenate((xyz, [[8, 8, 8]]))
    with pytest.raises(ValueError, match='crosslinked or nonprotein'):
        _validated_seed_torsion_targets(top, xyz, modeled, set(), targets)


def test_modified_target_and_ambiguous_peptide_branch_are_rejected(sample):
    top, xyz, residues, _, targets = sample
    original_name = residues[2].name
    residues[2].name = 'LIG'
    with pytest.raises(ValueError, match='canonical protein'):
        _validated_seed_torsion_targets(top, xyz, {residue_key(residues[2])}, set(), targets)
    residues[2].name = original_name
    top.addBond(next(a for a in residues[0].atoms() if a.name == 'C'),
                next(a for a in residues[2].atoms() if a.name == 'N'))
    with pytest.raises(ValueError, match='unambiguous peptide'):
        validate(sample)


def test_optional_prior_is_temporary_and_default_has_no_new_force(sample, monkeypatch):
    top, xyz, residues, modeled, targets = sample
    # Evaluate forces/energies only; do not run minimization or trajectory steps.
    monkeypatch.setattr(mm.LocalEnergyMinimizer, 'minimize', lambda *args: None)
    factory = BondOnlyFactory(xyz)
    _, default = _minimize_attempt(top, xyz * unit.nanometer, factory, modeled, max_iterations=1)
    assert 'initial_torsion_preservation' not in default
    assert not any(force.getName() == 'QualifiedSeedPhiPsiPreservation' for force in factory.systems[0].getForces())
    flanks = {residue_key(residues[1]), residue_key(residues[3])}
    targets = records(top, xyz, modeled | flanks)
    result, report = _minimize_attempt(top, xyz * unit.nanometer, factory, modeled,
        max_iterations=1, seed_torsion_targets=targets, mobile_flank_keys=flanks,
        reference_positions=xyz * unit.nanometer, preserve_context_peptides=True)
    temporary = factory.systems[-1]
    prior = next(force for force in temporary.getForces() if force.getName() == 'QualifiedSeedPhiPsiPreservation')
    assert prior.getNumTorsions() == 6
    assert all(prior.getTorsionParameters(i)[4][0] == 200 for i in range(6))
    provenance = report['initial_torsion_preservation']
    assert provenance['strength_kj_mol'] == 200
    assert provenance['seed_target_angles_verified']
    assert provenance['caller_reference_qualification_required']
    assert provenance['native_reference_check_performed_here'] is False
    assert provenance['exported_to_simulation'] is False
    assert np.isfinite(report['mobile_force_rms_kj_mol_nm'])
    np.testing.assert_array_equal(result.value_in_unit(unit.nanometer), xyz)
    assert report['max_iterations'] == 1
    assert report['energy_includes_construction_restraints']
    assert report['peptide_construction_restraints']['observed_context_peptides_preserved']
    assert report['flank_relaxation']['maximum_allowed_displacement_nm'] == .1
    assert report['flank_relaxation']['displacement_boundary_k_kj_mol_nm2'] == 1e6
    json.dumps(report, allow_nan=False)
    fresh = factory.createSystem(top, constraints=None)
    assert fresh.getNumForces() == 1
    assert not any(isinstance(force, mm.CustomTorsionForce) for force in fresh.getForces())
    assert all(fresh.getParticleMass(i).value_in_unit(unit.dalton) > 0 for i in range(fresh.getNumParticles()))

    # Isolate the actual added force, retaining its original atom indices.
    isolated = mm.System()
    for _ in top.atoms():
        isolated.addParticle(12.0)
    cloned = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(prior))
    isolated.addForce(cloned)
    integrator = mm.VerletIntegrator(.001)
    context = mm.Context(isolated, integrator, mm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(xyz * unit.nanometer)
        state = context.getState(getEnergy=True, getForces=True)
        assert state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) == pytest.approx(0, abs=1e-10)
        assert np.isfinite(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)).all()
        a, b, c, d, params = cloned.getTorsionParameters(0)
        cloned.setTorsionParameters(0, a, b, c, d, [params[0], params[1] + .2])
        cloned.updateParametersInContext(context)
        state = context.getState(getEnergy=True, getForces=True)
        assert state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) == pytest.approx(200 * (1 - np.cos(.2)), abs=1e-9)
        force = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))
        assert np.isfinite(force).all() and np.linalg.norm(force) > 0
        axis = int(np.argmax(abs(force[a])))
        step = 1e-6
        energies = []
        for direction in (-1, 1):
            moved = xyz.copy(); moved[a, axis] += direction * step
            context.setPositions(moved * unit.nanometer)
            energies.append(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        assert force[a, axis] == pytest.approx(-(energies[1] - energies[0]) / (2 * step), rel=1e-6, abs=1e-4)
    finally:
        del context, integrator
