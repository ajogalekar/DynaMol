from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pytest
from openmm import app, unit

from backend.stereo_monitor import StereoMonitor, require_valid
from backend.worker import Worker


SOURCE = Path(__file__).resolve().parents[2] / 'docs/audit/preparation-fixtures/six_residues_intact.pdb'


@pytest.fixture
def molecular_model():
    pdb = app.PDBFile(str(SOURCE))
    xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    return StereoMonitor(pdb.topology, xyz, SOURCE), xyz


def test_valid_residue_centers_pass_after_independent_periodic_wrapping(molecular_model):
    monitor, xyz = molecular_model
    box = np.diag([5., 5., 5.])
    wrapped = (xyz + [4.91, 3.94, 4.87]) % 5
    report = monitor.check(wrapped, box, 'Saved frame')
    assert report['passed'] and report['checked_centers'] > 5


@pytest.mark.parametrize('stage', ['After minimization', 'Initial relaxation step 100', 'Production step 500'])
def test_inversion_is_rejected_at_every_native_stage(molecular_model, stage):
    monitor, xyz = molecular_model
    changed = xyz.copy()
    ca, n, c, cb = monitor.indices[0]
    normal = np.cross(changed[n] - changed[ca], changed[c] - changed[ca])
    normal /= np.linalg.norm(normal)
    changed[cb] -= 2 * np.dot(changed[cb] - changed[ca], normal) * normal
    report = monitor.check(changed, stage=stage)
    assert not report['passed']
    assert report['violations'][0]['center'] == 'CA'
    with pytest.raises(ValueError, match='Dynamics stopped'):
        require_valid(report)
    assert monitor.check(xyz)['passed']


def test_flat_and_nonfinite_centers_are_not_accepted(molecular_model):
    monitor, xyz = molecular_model
    ca, n, c, cb = monitor.indices[0]
    for value in [xyz[ca], np.array([np.nan] * 3)]:
        changed = xyz.copy()
        changed[cb] = value
        assert not monitor.check(changed)['passed']


def make_worker(folder):
    worker = Worker.__new__(Worker)
    worker.folder = folder
    worker.provenance = {}
    return worker


def test_restarted_worker_preserves_native_phase_history(tmp_path):
    shutil.copy2(SOURCE, tmp_path / 'input.pdb')
    pdb = app.PDBFile(str(SOURCE))
    xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    make_worker(tmp_path).check_stereochemistry(pdb.topology, xyz, 'Before interrupted run')
    make_worker(tmp_path).check_stereochemistry(pdb.topology, xyz, 'After resume')
    history = json.loads((tmp_path / 'stereochemistry-checks.json').read_text())
    assert [entry['stage'] for entry in history['checks']] == ['Before interrupted run', 'After resume']


def test_failure_retains_matching_native_topology_and_each_attempt(tmp_path):
    shutil.copy2(SOURCE, tmp_path / 'input.pdb')
    pdb = app.PDBFile(str(SOURCE))
    model = app.Modeller(pdb.topology, pdb.positions)
    added = app.Topology()
    residue = added.addResidue('NA', added.addChain())
    added.addAtom('NA', app.element.sodium, residue)
    model.add(added, np.array([[9., 9., 9.]]) * unit.nanometer)
    xyz = np.asarray(model.positions.value_in_unit(unit.nanometer))
    monitor = StereoMonitor(model.topology, xyz, SOURCE)
    ca, n, c, cb = monitor.indices[0]
    normal = np.cross(xyz[n] - xyz[ca], xyz[c] - xyz[ca])
    normal /= np.linalg.norm(normal)
    xyz[cb] -= 2 * np.dot(xyz[cb] - xyz[ca], normal) * normal
    worker = make_worker(tmp_path)
    legacy_hash = None
    for attempt in range(2):
        with pytest.raises(ValueError, match='Dynamics stopped'):
            worker.check_stereochemistry(model.topology, xyz + attempt, f'Attempt {attempt}')
        report = json.loads((tmp_path / 'provenance.json').read_text())['stereochemistry_failure']
        for relative, expected in report['artifacts'].items():
            assert hashlib.sha256((tmp_path / relative).read_bytes()).hexdigest() == expected
        cif = next(tmp_path / key for key in report['artifacts'] if key.endswith('.cif'))
        saved = app.PDBxFile(str(cif))
        assert saved.topology.getNumAtoms() == model.topology.getNumAtoms() > pdb.topology.getNumAtoms()
        assert [(a.name, a.element) for a in saved.topology.atoms()] == [(a.name, a.element) for a in model.topology.atoms()]
        assert np.allclose(saved.positions.value_in_unit(unit.nanometer), xyz + attempt, atol=1e-5)
        current_hash = hashlib.sha256((tmp_path / 'stereochemistry-failure.npz').read_bytes()).hexdigest()
        if legacy_hash is not None:
            assert current_hash == legacy_hash
        legacy_hash = current_hash
    assert len(list((tmp_path / 'stereochemistry-failures').glob('*/coordinates.npz'))) == 2


def test_nonfinite_failure_coordinates_remain_in_npz_with_readable_topology(tmp_path, molecular_model):
    shutil.copy2(SOURCE, tmp_path / 'input.pdb')
    pdb = app.PDBFile(str(SOURCE))
    monitor, xyz = molecular_model
    xyz[monitor.indices[0, 3], 0] = np.nan
    worker = make_worker(tmp_path)
    with pytest.raises(ValueError, match='Dynamics stopped'):
        worker.check_stereochemistry(pdb.topology, xyz, 'Nonfinite coordinates')
    report = json.loads((tmp_path / 'provenance.json').read_text())['stereochemistry_failure']
    assert report['topology_coordinate_placeholders'] == 1
    npz = next(tmp_path / key for key in report['artifacts'] if key.endswith('.npz'))
    assert np.isnan(np.load(npz)['xyz_nm']).sum() == 1
    cif = next(tmp_path / key for key in report['artifacts'] if key.endswith('.cif'))
    assert np.isfinite(app.PDBxFile(str(cif)).positions.value_in_unit(unit.nanometer)).all()

def _fake_monitor(indices, reference_xyz):
    monitor = StereoMonitor.__new__(StereoMonitor)
    monitor.indices = np.asarray(indices, dtype=int).reshape(-1, 4)
    monitor.centers = [{'chain': 'A', 'resid': str(i + 1), 'residue': 'ALA', 'center': 'CA',
                        'atom_indices': list(map(int, idx)), 'insertion_code': '',
                        'template_signed_volume_nm3': _signed_volume(reference_xyz, idx)}
                       for i, idx in enumerate(monitor.indices)]
    monitor.expected = np.asarray([c['template_signed_volume_nm3'] for c in monitor.centers])
    return monitor


def _signed_volume(xyz, idx):
    p = np.asarray(xyz)[idx]
    return float(np.linalg.det([p[1] - p[0], p[2] - p[0], p[3] - p[0]]))


def test_chirality_guard_prevents_minimization_from_inverting_a_stereocenter():
    """A realistic distorting force that alone drives a center to its mirror
    image cannot invert it while the permanent guard is in the force field;
    without the guard the same force inverts it."""
    import openmm as mm
    from openmm import unit

    xyz = np.array([[0.0, 0.0, 0.05], [0.1, 0.0, 0.0],
                    [-0.05, 0.0866, 0.0], [-0.05, -0.0866, 0.0]])
    monitor = _fake_monitor([[0, 1, 2, 3]], xyz)
    v0 = _signed_volume(xyz, monitor.indices[0])
    assert abs(v0) > 1e-4

    def minimize(guard):
        system = mm.System()
        for mass in (12.0, 0.0, 0.0, 0.0):  # only the apex atom is free to move
            system.addParticle(mass)
        driver = mm.CustomExternalForce('0.5*K*(z-zt)^2')
        driver.addGlobalParameter('K', 4e4)   # ~4000 kJ/mol/nm at 0.1 nm: the real strained-loop scale
        driver.addGlobalParameter('zt', -0.05)  # pull the apex across the plane
        driver.addParticle(0, [])
        system.addForce(driver)
        if guard:
            force, guarded = monitor.chirality_guard_force(xyz, k=5000.0)
            assert len(guarded) == 1
            system.addForce(force)
        context = mm.Context(system, mm.VerletIntegrator(0.001 * unit.picosecond),
                             mm.Platform.getPlatformByName('Reference'))
        context.setPositions(xyz * unit.nanometer)
        mm.LocalEnergyMinimizer.minimize(context, 1e-6, 2000)
        return np.asarray(context.getState(getPositions=True)
                          .getPositions(asNumpy=True).value_in_unit(unit.nanometer))

    unguarded = minimize(False)
    guarded = minimize(True)
    assert np.sign(_signed_volume(unguarded, monitor.indices[0])) != np.sign(v0)
    assert np.sign(_signed_volume(guarded, monitor.indices[0])) == np.sign(v0)
    assert not monitor.check(unguarded)['passed']
    assert monitor.check(guarded)['passed']


def test_chirality_guard_is_flat_for_healthy_centers_and_rises_near_planarity():
    """The well is exactly zero while a center keeps its sign and >= 40% of its
    ideal volume, and only then rises -- so healthy centers feel no force."""
    import openmm as mm
    from openmm import unit

    xyz = np.array([[0.0, 0.0, 0.05], [0.1, 0.0, 0.0],
                    [-0.05, 0.0866, 0.0], [-0.05, -0.0866, 0.0]])
    monitor = _fake_monitor([[0, 1, 2, 3]], xyz)
    force, guarded = monitor.chirality_guard_force(xyz, k=5000.0, vmin_fraction=0.4)
    _, (s, vmin, vref) = force.getBondParameters(0)
    assert abs(vref - abs(_signed_volume(xyz, monitor.indices[0]))) < 1e-9
    assert abs(vmin - 0.4 * vref) < 1e-12

    system = mm.System()
    for _ in range(4):
        system.addParticle(12.0)
    system.addForce(force)
    context = mm.Context(system, mm.VerletIntegrator(0.001 * unit.picosecond),
                         mm.Platform.getPlatformByName('Reference'))
    # Healthy geometry: zero energy.
    context.setPositions(xyz * unit.nanometer)
    healthy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    assert abs(healthy) < 1e-6
    # Flatten the center toward the plane: energy must become positive.
    flattened = xyz.copy()
    flattened[0, 2] = 0.01  # apex nearly in the N/C/CB plane -> below the 40% floor
    context.setPositions(flattened * unit.nanometer)
    near_planar = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    assert near_planar > 1.0


def test_chirality_guard_skips_planar_and_inverted_centers(molecular_model):
    monitor, xyz = molecular_model
    # Every guarded bond keeps the center's own (correct) sign and targets 40% of ideal.
    force, guarded = monitor.chirality_guard_force(xyz, k=5000.0, vmin_fraction=0.4)
    names = [force.getGlobalParameterName(i) for i in range(force.getNumGlobalParameters())]
    assert names == ['k_chi'] and force.getGlobalParameterDefaultValue(0) == 5000.0
    assert 0 < len(guarded) == force.getNumBonds() <= len(monitor.centers)
    for b in range(force.getNumBonds()):
        atoms, (s, vmin, vref) = force.getBondParameters(b)
        assert vref > 0 and abs(vmin - 0.4 * vref) < 1e-12
        assert s == np.sign(_signed_volume(xyz, list(atoms)))
    # A planar center is skipped (nothing to hold), and an inverted one is left
    # for the before-minimization check to reject rather than pinned wrong-way.
    planar = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.05, 0.0866, 0.0], [-0.05, 0.0866, 0.0]])
    _, none_here = _fake_monitor([[0, 1, 2, 3]], planar).chirality_guard_force(planar)
    assert none_here == []
    changed = xyz.copy()
    ca, n, c, cb = monitor.indices[0]
    normal = np.cross(changed[n] - changed[ca], changed[c] - changed[ca]); normal /= np.linalg.norm(normal)
    changed[cb] -= 2 * np.dot(changed[cb] - changed[ca], normal) * normal  # mirror CB -> invert center 0
    _, after = monitor.chirality_guard_force(changed)
    assert not any(g['resid'] == monitor.centers[0]['resid'] and g['chain'] == monitor.centers[0]['chain']
                   and g['center'] == monitor.centers[0]['center'] for g in after)
