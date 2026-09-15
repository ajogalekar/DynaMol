import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('continuation', BASE/'continue_covalent_qm.py')
continuation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(continuation)


@pytest.fixture
def parent():
    cap = continuation.ROOT/'build/advanced-chemistry/covalent'
    if not (cap/'capped-adduct.sdf').exists():
        pytest.skip('Isolated research cap not present')
    request = json.loads((cap/'df-optimization-parallel-v1-input.json').read_text())
    bohr = 0.529177210903
    xyz = np.array(request['coords_angstrom'])
    # Synthetic numerical metadata only to exercise geometry guards. No QM
    # result is written and this is not a claim of actual parent convergence.
    result = {'accepted': True, 'optimization': {'converged': True}, 'atom_ids': request['atom_ids'],
              'elements': request['elements'], 'units': {'bohr_to_angstrom': bohr}}
    return cap, request, result, {'coords_angstrom': xyz, 'coords_bohr': xyz/bohr}


def test_source_geometry_and_frozen_caps(parent):
    _, report = continuation.check_geometry(*parent)
    assert report['passed']
    assert len(report['cap_heavy_ids']) == 5
    assert report['maximum_cap_displacement_angstrom'] == 0


def test_failed_optimization_cannot_start_charge_fit(parent):
    parent[2]['accepted'] = False
    with pytest.raises(ValueError, match='converged'):
        continuation.check_geometry(*parent)


def test_atom_reordering_cannot_start_charge_fit(parent):
    parent[2]['atom_ids'] = list(reversed(parent[2]['atom_ids']))
    with pytest.raises(ValueError, match='identities/order'):
        continuation.check_geometry(*parent)


def test_frozen_cap_movement_is_rejected(parent):
    parent[3]['coords_angstrom'][0, 0] += .1
    parent[3]['coords_bohr'] = parent[3]['coords_angstrom']/parent[2]['units']['bohr_to_angstrom']
    with pytest.raises(ValueError, match='Frozen cap'):
        continuation.check_geometry(*parent)


def test_conflicting_declared_stereochemistry_is_rejected(parent):
    source = json.loads((parent[0]/'capped-adduct.json').read_text())
    source['stereochemistry']['CA'] = 'S'
    # Read substitute metadata via a temporary source copy, with original
    # molecule geometry untouched. This exercises product-state enforcement.
    import tempfile
    import shutil
    with tempfile.TemporaryDirectory() as temporary:
        folder = Path(temporary)
        shutil.copyfile(parent[0]/'capped-adduct.sdf', folder/'capped-adduct.sdf')
        (folder/'capped-adduct.json').write_text(json.dumps(source))
        with pytest.raises(ValueError, match='stereochemistry'):
            continuation.check_geometry(folder, *parent[1:])


def test_broken_link_geometry_is_rejected(parent):
    parent[3]['coords_angstrom'][18] += [4, 0, 0]
    parent[3]['coords_bohr'] = parent[3]['coords_angstrom']/parent[2]['units']['bohr_to_angstrom']
    with pytest.raises(ValueError, match='bonded distance'):
        continuation.check_geometry(*parent)
