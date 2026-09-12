"""Exact native lattice comparisons, boundary and bounded-allocation controls."""
import itertools
import math
from pathlib import Path

import numpy as np
import pytest
from openmm import app, unit
from openmm.app import modeller

from backend.solvent import _accepted_axis_tiles, tip3p_allocation_bound


@pytest.mark.parametrize('oxygen,period,side', [
    (0.0, 3.0, 3.0), (0.0, 3.0, math.nextafter(3.0, math.inf)),
    (0.0, 3.0, math.nextafter(3.0, -math.inf)),
    (-.001, 3.0, 3.0), (3.0, 3.0, 3.0), (.5, 3.0, 3.5),
    (3.1, 3.0, 3.0),
])
def test_axis_count_matches_literal_native_inclusive_boundary(oxygen, period, side):
    expected = sum(oxygen + index * period <= side for index in range(math.ceil(side / period)))
    assert _accepted_axis_tiles(oxygen, period, side) == expected


@pytest.mark.parametrize('side', [2., 3., math.nextafter(3., math.inf), 5.93, 9.6936])
def test_bound_contains_every_native_template_placement(side):
    source = Path(modeller.__file__).parent / 'data/tip3p.pdb'
    template = app.PDBFile(str(source))
    xyz = np.asarray(template.positions.value_in_unit(unit.nanometer))
    periods = template.topology.getUnitCellDimensions().value_in_unit(unit.nanometer)
    oxygens = [xyz[a.index] for a in template.topology.atoms() if a.element == app.element.oxygen]
    solute = np.array([[0., 0., 0.], [side - 1, 0., 0.]])
    before = solute.copy()
    bound = tip3p_allocation_bound(solute, 1)
    native_side = bound['box_side_nm']
    exact = sum(all(float(oxygen[i]) + tiles[i] * float(periods[i]) <= native_side for i in range(3))
                for oxygen in oxygens
                for tiles in itertools.product(*(range(math.ceil(native_side / periods[i])) for i in range(3))))
    assert bound['candidate_water_molecules'] >= exact
    assert bound['maximum_total_atoms'] == 2 + 3 * bound['candidate_water_molecules']
    assert bound['maximum_total_atoms'] >= 2 + 3 * exact
    assert bound['count_boundary_nm'] >= native_side
    assert np.array_equal(solute, before)
    assert len(bound['template_sha256']) == 64


def test_large_box_is_counted_without_allocating_a_solvent_grid():
    bound = tip3p_allocation_bound(np.array([[0., 0., 0.], [1e6, 0., 0.]]), 1)
    assert bound['maximum_total_atoms'] > 10**18
    assert math.isfinite(bound['box_side_nm'])


def test_padding_and_coordinates_validated_before_template_allocation():
    with pytest.raises(ValueError, match='coordinates'):
        tip3p_allocation_bound(np.array([[0., float('nan'), 0.]]), 1)
    with pytest.raises(ValueError, match='padding'):
        tip3p_allocation_bound(np.zeros((1, 3)), float('inf'))


def test_same_positions_quantity_and_array_have_identical_bound():
    coordinates = np.array([[1.3, -2.4, .8], [5.1, 3.4, 2.3], [-1.5, 2., .3]])
    assert tip3p_allocation_bound(coordinates, 1.2) == tip3p_allocation_bound(coordinates * unit.nanometer, 1.2)
