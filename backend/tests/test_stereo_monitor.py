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
