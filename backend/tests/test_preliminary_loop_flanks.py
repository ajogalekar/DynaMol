"""Preliminary chi sampling must preserve future loop-refinement references."""
import numpy as np
from openmm import app, unit

from backend.preparation_worker import adjust_sidechains, observed_loop_flanks
from backend.residue_identity import residue_key


def protein(chain_lengths):
    topology, xyz, residues, atoms = app.Topology(), [], {}, {}
    coordinates = {'N': [0, .145, 0], 'CA': [0, 0, 0], 'C': [0, -.14, .03],
                   'O': [0, -.25, .06], 'CB': [.15, 0, 0], 'OG': [.20, .10, .05]}
    for chain_id, length in chain_lengths.items():
        chain, previous = topology.addChain(chain_id), None
        for i in range(1, length+1):
            residue = topology.addResidue('SER', chain, str(100-i*7), 'B' if i == 3 else '')
            residues[chain_id, i] = residue
            named = {}
            origin = np.array([len(residues)*3., 0, 0])
            for name, point in coordinates.items():
                element = app.element.nitrogen if name == 'N' else app.element.oxygen if name in {'O', 'OG'} else app.element.carbon
                named[name] = topology.addAtom(name, element, residue)
                xyz.append(origin + point)
            for a, b in [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'OG')]:
                topology.addBond(named[a], named[b])
            if previous is not None:
                topology.addBond(previous, named['N'])
            previous = named['C']
            atoms[chain_id, i] = named
    return topology, np.array(xyz), residues, atoms


def test_multiple_gaps_and_chains_use_connectivity_with_exact_source_ids():
    top, xyz, residues, atoms = protein({'A': 8, 'B': 6})
    modeled = {residue_key(residues[k]) for k in [('A', 2), ('A', 3), ('A', 6), ('B', 4)]}
    expected = {residue_key(residues[k]) for k in [('A', 1), ('A', 4), ('A', 5), ('A', 7), ('B', 3), ('B', 5)]}
    assert observed_loop_flanks(top, modeled) == expected
    assert not observed_loop_flanks(top, modeled) & modeled
    assert not observed_loop_flanks(top, ())


def test_crosschain_connections_and_sidechain_crosslinks_are_not_peptide_neighbors():
    top, xyz, residues, atoms = protein({'A': 4, 'B': 4})
    modeled = {residue_key(residues['A', 2])}
    top.addBond(atoms['A', 2]['C'], atoms['B', 2]['N'])
    top.addBond(atoms['A', 2]['CB'], atoms['B', 3]['CB'])
    assert observed_loop_flanks(top, modeled) == {residue_key(residues['A', i]) for i in [1, 3]}


def test_missing_connection_is_not_inferred_from_order_or_coordinates():
    top, xyz, residues, atoms = protein({'A': 4})
    cut = {atoms['A', 1]['C'], atoms['A', 2]['N']}
    top._bonds = [bond for bond in top._bonds if set(bond) != cut]
    assert observed_loop_flanks(top, {residue_key(residues['A', 2])}) == {residue_key(residues['A', 3])}


def test_flanks_and_metals_stay_exact_while_other_chi_sampling_still_operates():
    top, xyz, residues, atoms = protein({'A': 4, 'M': 1})
    # Give every sidechain the same synthetic steric conflict. This tests
    # sampling behavior, not physical loop geometry or a native conformation.
    environment = top.addChain('E')
    points = list(xyz)
    for i, named in enumerate(atoms.values(), 1):
        obstacle = top.addResidue('LIG', environment, str(i))
        top.addAtom('O', app.element.oxygen, obstacle)
        points.append(xyz[named['OG'].index] + [0, 0, .01])
    xyz = np.array(points)
    source = xyz.copy()
    modeled = {residue_key(residues['A', 2])}
    metal = {residue_key(residues['M', 1])}
    old, old_report = adjust_sidechains(top, xyz * unit.nanometer, protected_residues=metal)
    new, report = adjust_sidechains(top, xyz * unit.nanometer, protected_residues=metal,
                                    modeled_residue_keys=modeled)
    before_fix = np.asarray(old.value_in_unit(unit.nanometer))
    after_fix = np.asarray(new.value_in_unit(unit.nanometer))
    flanks = {residue_key(residues['A', i]) for i in [1, 3]}
    for key in [('A', 1), ('A', 3)]:
        indices = [a.index for a in residues[key].atoms()]
        assert not np.array_equal(before_fix[indices], source[indices])
        np.testing.assert_array_equal(after_fix[indices], source[indices])
    protected_indices = [a.index for a in residues['M', 1].atoms()]
    np.testing.assert_array_equal(after_fix[protected_indices], source[protected_indices])
    changed_residues = {(r['chain'], r['resid']) for r in report['adjustments']}
    assert changed_residues == {('A', residues['A', i].id) for i in [2, 4]}
    assert {tuple(k) for k in report['preserved_loop_flanks']} == flanks
    assert report['residues_examined'] == 2
    assert old_report['preserved_loop_flanks'] == []
    assert any('observed metal coordination' in r for r in report['skipped'])
    assert sum('observed loop flank' in r for r in report['skipped']) == 2
    assert all('metal coordination' not in r for r in report['skipped'] if 'observed loop flank' in r)
    assert metal == {residue_key(residues['M', 1])}
    np.testing.assert_array_equal(xyz, source)


def test_overlapping_loop_and_metal_protection_reports_both_reasons():
    top, xyz, residues, atoms = protein({'A': 3})
    flank = residue_key(residues['A', 1])
    moved, report = adjust_sidechains(top, xyz * unit.nanometer,
                                     protected_residues={flank},
                                     modeled_residue_keys={residue_key(residues['A', 2])})
    record = next(row for row in report['skipped'] if f'A:{residues["A", 1].id}:SER' in row)
    assert 'observed metal coordination' in record and 'observed loop flank' in record
    ix = [a.index for a in residues['A', 1].atoms()]
    np.testing.assert_array_equal(np.asarray(moved.value_in_unit(unit.nanometer))[ix], xyz[ix])
