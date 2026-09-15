"""Protect original residue identities and native particle-index correspondence."""
import numpy as np
import pytest
from openmm import app

from backend.native_atom_mapping import build_native_atom_map


def specimen():
    native = app.Topology()
    residue = native.addResidue("COV", native.addChain("1"), id="1")
    # Native adduct heavy atoms followed by hydrogens; source residues interleave.
    atoms = [native.addAtom(name, element, residue) for name, element in
             [("S1", app.element.sulfur), ("C2", app.element.carbon),
              ("C3", app.element.carbon), ("H4", app.element.hydrogen),
              ("H5", app.element.hydrogen)]]
    for a, b in [(0, 1), (0, 2), (1, 3), (2, 4)]:
        native.addBond(atoms[a], atoms[b])
    ids = ["A:12:B:CYS:SG", "L:303::LIG:C25", "A:12:B:CYS:CB", "L:303::LIG:H1", "A:12:B:CYS:HB2"]
    records = [{"native_index": i, "source_or_added_id": identity} for i, identity in enumerate(ids)]
    return native, records


def test_interleaved_covalent_adduct_restores_labels_and_bonds_without_coordinate_change():
    native, records = specimen()
    mapping = build_native_atom_map(native, records)
    assert mapping.display_to_native == (0, 2, 4, 1, 3)
    assert [(r.chain.id, r.id, r.insertionCode, r.name) for r in mapping.display_topology.residues()] == [
        ("A", "12", "B", "CYS"), ("L", "303", "", "LIG")]
    assert mapping.measurement_to_native([0, 3], expected_mapping_sha256=mapping.mapping_sha256) == [0, 1]
    xyz = np.arange(5 * 3, dtype=float).reshape(5, 3)
    frames = np.stack([xyz, xyz + 0.25])
    shown = mapping.to_display(frames)
    assert np.array_equal(mapping.to_native(shown), frames)
    native_edges = {tuple(sorted((b[0].index, b[1].index))) for b in native.bonds()}
    mapped_edges = {tuple(sorted(mapping.measurement_to_native([b[0].index, b[1].index], expected_mapping_sha256=mapping.mapping_sha256)))
                    for b in mapping.display_topology.bonds()}
    assert native_edges == mapped_edges
    # The distance selected in the viewer uses the same actual engine atoms.
    assert np.array_equal(np.linalg.norm(shown[:, 0] - shown[:, 3], axis=1),
                          np.linalg.norm(frames[:, 0] - frames[:, 1], axis=1))


@pytest.mark.parametrize("selection", [[-1], [5], [True], [1.5]])
def test_measurements_reject_stale_or_noninteger_indices(selection):
    mapping = build_native_atom_map(*specimen())
    with pytest.raises(ValueError, match="index"):
        mapping.measurement_to_native(selection, expected_mapping_sha256=mapping.mapping_sha256)


def test_missing_and_duplicate_source_records_are_rejected():
    native, records = specimen()
    with pytest.raises(ValueError, match="every native atom"):
        build_native_atom_map(native, records[:-1])
    records[1]["source_or_added_id"] = records[0]["source_or_added_id"]
    with pytest.raises(ValueError, match="duplicated"):
        build_native_atom_map(native, records)


def test_atom_order_cannot_silently_change_without_a_new_mapping():
    native, records = specimen()
    with pytest.raises(ValueError, match="native particle order"):
        build_native_atom_map(native, records[::-1])
    original = build_native_atom_map(native, records)
    records[0]["source_or_added_id"] = "A:12:C:CYS:SG"
    changed = build_native_atom_map(native, records)
    assert original.mapping_sha256 != changed.mapping_sha256


@pytest.mark.parametrize("coordinates", [np.zeros((4, 3)), np.zeros((5, 2)), np.full((5, 3), np.nan)])
def test_incomplete_or_invalid_coordinate_inventory_is_rejected(coordinates):
    with pytest.raises(ValueError, match="complete mapped inventory"):
        build_native_atom_map(*specimen()).to_display(coordinates)


def test_in_range_saved_selection_is_rejected_after_mapping_identity_changes():
    native,records=specimen()
    old=build_native_atom_map(native,records)
    saved_indices=[0,3];saved_hash=old.mapping_sha256
    records[0]['source_or_added_id']='A:13:B:CYS:SG'
    new=build_native_atom_map(native,records)
    with pytest.raises(ValueError,match='another atom mapping'):
        new.measurement_to_native(saved_indices,expected_mapping_sha256=saved_hash)
    assert old.measurement_to_native(saved_indices,expected_mapping_sha256=saved_hash)==[0,1]


def test_native_labels_are_bound_even_when_elements_and_bonds_stay_equal():
    native,records=specimen();old=build_native_atom_map(native,records)
    list(native.atoms())[1].name='RENAMED'
    new=build_native_atom_map(native,records)
    assert new.mapping_sha256!=old.mapping_sha256
    with pytest.raises(ValueError,match='checksum'):
        new.measurement_to_native([0,3],expected_mapping_sha256=old.mapping_sha256)


def test_conflicting_source_identity_forms_and_wrong_native_identity_are_rejected():
    native,records=specimen()
    records[0]['source_identity']={'chain':'A','resid':'13','insertion':'B','resname':'CYS','atomname':'SG'}
    with pytest.raises(ValueError,match='identities conflict'):build_native_atom_map(native,records)
    records[0]['source_identity']['resid']='12'
    assert build_native_atom_map(native,records)
    records[0]['native_identity']={'chain':'1','resid':'1','insertion':'','resname':'COV','atomname':'C2'}
    with pytest.raises(ValueError,match='native identity'):build_native_atom_map(native,records)


@pytest.mark.parametrize('indices',[[],[0,0],np.array(1),np.array([[0,1]])])
def test_invalid_measurement_index_shapes_and_repeated_atoms_are_rejected(indices):
    mapping=build_native_atom_map(*specimen())
    with pytest.raises(ValueError,match='indices'):
        mapping.measurement_to_native(indices,expected_mapping_sha256=mapping.mapping_sha256)


def test_complex_coordinates_are_not_a_physical_trajectory():
    mapping=build_native_atom_map(*specimen())
    with pytest.raises(ValueError,match='finite'):
        mapping.to_display(np.ones((5,3),dtype=complex)*(1+1j))
