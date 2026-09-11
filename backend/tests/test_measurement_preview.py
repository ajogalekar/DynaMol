"""Live readouts share physical/periodic geometry with the saved plot contract."""
import time

import mdtraj as md
import numpy as np
import pytest
from fastapi.testclient import TestClient

from backend import analysis
from backend.main import app
from backend.models import MeasurementRequest


def trajectory(xyz, *, periodic=False):
    top = md.Topology()
    residue = top.add_residue('TST', top.add_chain('A'), resSeq=1)
    for name, element in [('N', md.element.nitrogen), ('H', md.element.hydrogen), ('O', md.element.oxygen), ('C', md.element.carbon)]:
        top.add_atom(name, element, residue)
    top.add_bond(top.atom(0), top.atom(1))
    traj = md.Trajectory(np.asarray(xyz, dtype=np.float32), top)
    if periodic:
        traj.unitcell_lengths = np.ones((traj.n_frames, 3), dtype=np.float32)
        traj.unitcell_angles = np.full((traj.n_frames, 3), 90, dtype=np.float32)
    return traj


@pytest.mark.parametrize('kind,atoms,expected', [
    ('distance', [0, 1], 1.0), ('angle', [0, 1, 2], 90.0), ('dihedral', [0, 1, 2, 3], 90.0),
])
def test_preview_equals_plot_geometry(kind, atoms, expected, monkeypatch):
    traj = trajectory([[[0, 0, 0], [.1, 0, 0], [.1, .1, 0], [.1, .1, .1]]])
    monkeypatch.setattr(analysis, 'load_physical', lambda _: traj)
    request = MeasurementRequest(kind=kind, atoms=atoms)
    value = analysis.preview('fixture', request)
    assert value['values'] == pytest.approx([expected], abs=1e-5)
    assert value['values'] == analysis.measure('fixture', request)['values']
    assert value['frame_errors'] == [None]
    assert 'id' not in value


def test_preview_uses_minimum_image_not_display_separation(monkeypatch):
    traj = trajectory([[[.05, 0, 0], [.95, 0, 0], [0, .2, 0], [0, 0, .3]]], periodic=True)
    monkeypatch.setattr(analysis, 'load_physical', lambda _: traj)
    result = analysis.preview('fixture', MeasurementRequest(kind='distance', atoms=[0, 1]))
    assert result['values'] == pytest.approx([1], abs=1e-5)
    assert result['warnings'] == []


def test_preview_preserves_valid_frames_when_other_frames_degenerate(monkeypatch):
    traj = trajectory([
        [[0, 0, 0], [.1, 0, 0], [.1, .1, 0], [.1, .1, .1]],
        [[0, 0, 0], [0, 0, 0], [.1, .1, 0], [.1, .1, .1]],
        [[0, 0, 0], [.1, 0, 0], [.2, 0, 0], [.3, 0, 0]],
    ])
    monkeypatch.setattr(analysis, 'load_physical', lambda _: traj)
    request = MeasurementRequest(kind='dihedral', atoms=[0, 1, 2, 3])
    result = analysis.preview('fixture', request)
    assert result['values'][0] == pytest.approx(90, abs=1e-5)
    assert result['values'][1:] == [None, None]
    assert 'coincident' in result['frame_errors'][1]
    assert 'collinear' in result['frame_errors'][2]
    with pytest.raises(ValueError, match='coincident'):
        analysis.measure('fixture', request)
    with TestClient(app) as client:
        response = client.post('/api/datasets/fixture/measurement-preview', json=request.model_dump())
    assert response.status_code == 200
    assert response.json()['values'][1:] == [None, None]


def test_hbond_preview_includes_current_geometry_without_claiming_chemistry(monkeypatch):
    traj = trajectory([
        [[0, 0, 0], [.1, 0, 0], [.3, 0, 0], [0, .4, 0]],
        [[0, 0, 0], [.1, 0, 0], [.4, 0, 0], [0, .4, 0]],
    ])
    monkeypatch.setattr(analysis, 'load_physical', lambda _: traj)
    result = analysis.preview('fixture', MeasurementRequest(kind='hbond', atoms=[0, 1, 2]))
    assert result['values'] == pytest.approx([3, 4], abs=1e-5)
    assert result['angle_values'] == pytest.approx([180, 180])
    assert result['geometry_passes'] == [True, False]
    assert 'occupancy' not in result
    assert any('do not establish chemical' in warning for warning in result['warnings'])
    with pytest.raises(ValueError, match='bonded to the donor'):
        analysis.preview('fixture', MeasurementRequest(kind='hbond', atoms=[2, 1, 0]))


@pytest.mark.parametrize('atoms', [[0, 0], [-1, 1], [0, 5], [0, 1, 2]])
def test_preview_rejects_invalid_selection(atoms, monkeypatch):
    monkeypatch.setattr(analysis, 'load_physical', lambda _: trajectory([[[0, 0, 0], [.1, 0, 0], [.1, .1, 0], [.1, .1, .1]]]))
    with TestClient(app) as client:
        response = client.post('/api/datasets/fixture/measurement-preview', json={'kind': 'distance', 'atoms': atoms})
    assert response.status_code == 422


def test_full_frame_limit_computes_only_selected_geometry_promptly(monkeypatch):
    # A reproducible throughput bound for 10,000 frames; excludes physical file I/O.
    traj = trajectory([[[0, 0, 0], [.1, 0, 0], [.1, .1, 0], [.1, .1, .1]]] * 10000)
    monkeypatch.setattr(analysis, 'load_physical', lambda _: traj)
    start = time.monotonic()
    result = analysis.preview('fixture', MeasurementRequest(kind='dihedral', atoms=[0, 1, 2, 3]))
    assert len(result['values']) == 10000
    assert time.monotonic() - start < 2.0
