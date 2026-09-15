"""Serialization/identity fixtures only: no minimization, parameter fitting or MD."""
from dataclasses import FrozenInstanceError
import hashlib
import io
import json

import numpy as np
from openmm import app, unit
import pytest

from backend.loop_search import REQUIRED_CHECKS
from backend.loop_snapshot import (create_snapshot, load_snapshot,
                                   create_accepted_bundle, verify_accepted_bundle)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residue(topology, chain, name, resid, atoms):
    result = topology.addResidue(name, chain, resid)
    added = {name: topology.addAtom(name, element, result) for name, element in atoms}
    return result, added


@pytest.fixture
def case(tmp_path):
    # Resolve the macOS /var alias once at fixture creation. Production readers
    # correctly reject any explicitly supplied symlink path.
    root = tmp_path.resolve()
    source = app.Topology(); chain = source.addChain('A')
    backbone = [('N', app.element.nitrogen), ('CA', app.element.carbon),
                ('C', app.element.carbon), ('O', app.element.oxygen), ('H', app.element.hydrogen)]
    _, old = residue(source, chain, 'ALA', '1', backbone)
    for a, b in [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('N', 'H')]:
        source.addBond(old[a], old[b])
    _, water = residue(source, source.addChain('W'), 'HOH', '9', [('O', app.element.oxygen), ('H1', app.element.hydrogen)])
    source.addBond(water['O'], water['H1'])
    _, ligand = residue(source, source.addChain('B'), 'LIG', '3', [('OLD', app.element.carbon)])
    prepared = app.Topology(); chain = prepared.addChain('A')
    _, first = residue(prepared, chain, 'ALA', '1', backbone)
    _, modeled = residue(prepared, chain, 'GLY', '2', backbone)
    for atoms in [first, modeled]:
        for a, b in [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('N', 'H')]:
            prepared.addBond(atoms[a], atoms[b])
    prepared.addBond(first['C'], modeled['N'])
    _, renamed = residue(prepared, prepared.addChain('B'), 'LIG', '3', [('NEW', app.element.carbon)])
    parameters = root / 'inputs'; parameters.mkdir()
    (parameters / 'main.xml').write_text('<ForceField><Include file="nested/child.xml"/></ForceField>')
    (parameters / 'child.xml').write_text('<ForceField><Include file="../base.xml"/></ForceField>')
    (parameters / 'base.xml').write_text('<ForceField/>')
    (parameters / 'source.pdb').write_text('Explicit source-file fixture; topology supplied separately.\n')
    args = dict(source_topology=source, source_positions=np.arange(24).reshape(8, 3) * .01 * unit.nanometer,
                topology=prepared, positions=np.arange(33).reshape(11, 3) * .01 * unit.nanometer,
                modeled_keys=[('A', '2', '', 'GLY')], requested_modeled_keys=[('A', '2', '', 'GLY')],
                source_to_prepared=[0, 1, 2, 3, None, None, None, 10],
                preparation_provenance={'retention': {'remove_waters': True}, 'hydrogens': 'regenerated',
                                        'renames': [{'source': 'B:3::LIG:OLD', 'prepared': 'B:3::LIG:NEW'}]},
                parameter_files={'main.xml': parameters / 'main.xml', 'nested/child.xml': parameters / 'child.xml', 'base.xml': parameters / 'base.xml'},
                source_files={'input.pdb': parameters / 'source.pdb'})
    return root, args


def snapshot(case):
    root, args = case
    return create_snapshot(root / 'snapshot', **args)


def rewrite_snapshot(snap, change):
    data = json.loads(snap.manifest_path.read_text()); change(data)
    snap.manifest_path.write_text(json.dumps(data))
    return digest(snap.manifest_path)


def array_bytes(value, npz=False):
    stream = io.BytesIO()
    if npz:
        np.savez(stream, xyz_nm=value)
    else:
        np.save(stream, value, allow_pickle=False)
    return stream.getvalue()


def bundle_inputs(case, snap):
    root, _ = case
    (root / 'final.npz').write_bytes(array_bytes(snap.prepared_xyz_nm, npz=True))
    (root / 'prepared.cif').write_text('Explicit topology fixture; scientific check supplied by validator.\n')
    (root / 'validation.json').write_text(json.dumps({'checks': dict.fromkeys(REQUIRED_CHECKS, True)}))
    (root / 'atom-map.json').write_text(json.dumps({'source_to_prepared': snap.source_to_prepared}))
    return {'coordinates': {'final.npz': root / 'final.npz'},
            'topology': {'prepared.cif': root / 'prepared.cif'},
            'parameters': {f.name: f.path for f in snap.parameter_files},
            'provenance': {'validation.json': root / 'validation.json', 'atom-map.json': root / 'atom-map.json'}}


def test_roundtrip_preserves_different_source_prepared_and_explicit_mapping(case):
    snap = snapshot(case)
    restored = load_snapshot(snap.manifest_path.parent, expected_sha256=snap.sha256)
    assert len(restored.source.atoms) == 8 and len(restored.prepared.atoms) == 11
    assert restored.source_to_prepared == (0, 1, 2, 3, None, None, None, 10)
    assert restored.source.atoms[7].identity[-1] == 'OLD'
    assert restored.prepared.atoms[10].identity[-1] == 'NEW'
    assert restored.removed_source_identities == tuple(restored.source.atoms[i].identity for i in (4, 5, 6))
    assert restored.source.hydrogen_parents == ((4, (0,)), (6, (5,)))
    assert restored.prepared.hydrogen_parents == ((4, (0,)), (9, (5,)))
    assert restored.modeled_keys == (('A', '2', '', 'GLY'),)
    np.testing.assert_array_equal(restored.source_xyz_nm, np.arange(24).reshape(8, 3) * .01)
    np.testing.assert_array_equal(restored.prepared_xyz_nm, np.arange(33).reshape(11, 3) * .01)
    assert restored.preparation_provenance['renames'][0]['prepared'] == 'B:3::LIG:NEW'
    assert {f.name for f in restored.parameter_files} == {'main.xml', 'nested/child.xml', 'base.xml'}


def test_arrays_metadata_and_records_are_immutable(case):
    snap = snapshot(case)
    with pytest.raises(ValueError):
        snap.prepared_xyz_nm[0, 0] = 20
    with pytest.raises(ValueError):
        snap.prepared_xyz_nm.setflags(write=True)
    with pytest.raises(TypeError):
        snap.preparation_provenance['retention']['remove_waters'] = False
    with pytest.raises(FrozenInstanceError):
        snap.sha256 = '0' * 64


def test_copies_are_independent_of_later_original_file_changes(case):
    snap = snapshot(case)
    for path in case[1]['parameter_files'].values():
        path.write_text('original changed after capture')
    case[1]['source_files']['input.pdb'].unlink()
    assert load_snapshot(snap.manifest_path.parent, expected_sha256=snap.sha256).sha256 == snap.sha256


@pytest.mark.parametrize('field,value', [('positions', np.zeros((10, 3))),
    ('positions', np.full((11, 3), np.nan)), ('positions', np.full((11, 3), np.inf)),
    ('positions', np.zeros((11, 3), complex)), ('positions', np.ones((11, 3), bool)),
    ('source_positions', np.zeros((7, 3)))])
def test_reject_invalid_coordinates_before_freezing(case, field, value):
    root, args = case; args[field] = value
    with pytest.raises(ValueError):
        create_snapshot(root / 'snapshot', **args)


@pytest.mark.parametrize('field,value', [('modeled_keys', []),
    ('requested_modeled_keys', [('A', '2', '', 'GLY'), ('A', '3', '', 'ALA')]),
    ('modeled_keys', [('A', '2', '', 'GLY'), ('A', '2', '', 'GLY')]),
    ('source_to_prepared', [0]), ('source_to_prepared', [0, 0, 2, 3, None, None, None, 10]),
    ('source_to_prepared', [True, 1, 2, 3, None, None, None, 10]),
    ('source_to_prepared', [1, 0, 2, 3, None, None, None, 10]),
    ('preparation_provenance', {})])
def test_reject_partial_coverage_or_ambiguous_retention(case, field, value):
    root, args = case; args[field] = value
    with pytest.raises(ValueError):
        create_snapshot(root / 'snapshot', **args)


def test_duplicate_atom_identity_rejected(case):
    root, args = case
    atoms = list(args['topology'].atoms()); atoms[1].name = atoms[2].name
    with pytest.raises(ValueError, match='duplicate atom'):
        create_snapshot(root / 'snapshot', **args)


def test_missing_loop_backbone_rejected(case):
    root, args = case
    list(args['topology'].atoms())[8].name = 'WRONG'
    with pytest.raises(ValueError, match='missing prepared backbone'):
        create_snapshot(root / 'snapshot', **args)


def test_prepared_hydrogen_with_two_parents_rejected(case):
    root, args = case
    atoms = list(args['topology'].atoms()); args['topology'].addBond(atoms[9], atoms[6])
    with pytest.raises(ValueError, match='one exact heavy parent'):
        create_snapshot(root / 'snapshot', **args)


def test_raw_unbonded_hydrogen_is_retained_as_unknown_not_guessed(case):
    root, args = case
    args['source_topology']._bonds.pop(3)
    snap = create_snapshot(root / 'snapshot', **args)
    assert snap.source.hydrogen_parents[0] == (4, ())
    assert snap.prepared.hydrogen_parents[0] == (4, (0,))


def test_source_deuterium_element_is_not_silently_changed_to_hydrogen(case):
    root, args = case
    list(args['source_topology'].atoms())[4].element = app.element.deuterium
    snap = create_snapshot(root / 'snapshot', **args)
    assert snap.source.atoms[4].element == 'D'
    assert snap.prepared.atoms[4].element == 'H'
    assert snap.source_to_prepared[4] is None  # Explicitly removed, not equated.
    args['source_to_prepared'][4] = 4
    with pytest.raises(ValueError, match='changes an atom element'):
        create_snapshot(root / 'equated-isotopes', **args)


def test_missing_recursive_parameter_include_rejected(case):
    root, args = case; del args['parameter_files']['base.xml']
    with pytest.raises(ValueError, match='absent from the explicit file list'):
        create_snapshot(root / 'snapshot', **args)


@pytest.mark.parametrize('content', ['<ForceField><Include file="/tmp/external.xml"/></ForceField>',
    '<ForceField><Include file="../../external.xml"/></ForceField>',
    '<ForceField><Include file="main.xml"/></ForceField>'])
def test_escaping_or_cyclic_include_rejected(case, content):
    root, args = case; args['parameter_files']['main.xml'].write_text(content)
    with pytest.raises(ValueError):
        create_snapshot(root / 'snapshot', **args)


@pytest.mark.parametrize('name', ['../escape.xml', '/absolute.xml', 'a/../x.xml', 'a\\x.xml'])
def test_manifest_path_escape_rejected(case, name):
    root, args = case; args['source_files'] = {name: args['source_files']['input.pdb']}
    with pytest.raises(ValueError):
        create_snapshot(root / 'snapshot', **args)


def test_symlink_input_and_snapshot_root_rejected(case):
    root, args = case
    link = root / 'source-link'; link.symlink_to(args['source_files']['input.pdb'])
    original = args['source_files']; args['source_files'] = {'input.pdb': link}
    with pytest.raises(ValueError, match='Symlink'):
        create_snapshot(root / 'snapshot', **args)
    args['source_files'] = original
    snap = create_snapshot(root / 'snapshot', **args)
    link = root / 'snapshot-link'; link.symlink_to(snap.manifest_path.parent, target_is_directory=True)
    with pytest.raises(ValueError, match='Symlink'):
        load_snapshot(link, expected_sha256=snap.sha256)


@pytest.mark.parametrize('relative', ['parameters/main.xml', 'sources/input.pdb', 'coordinates/source.npy'])
def test_stale_copied_file_bytes_rejected(case, relative):
    snap = snapshot(case); (snap.manifest_path.parent / relative).write_bytes(b'changed')
    with pytest.raises(ValueError, match='SHA256'):
        load_snapshot(snap.manifest_path.parent, expected_sha256=snap.sha256)


def test_mutated_manifest_cannot_self_authorize(case):
    snap = snapshot(case)
    rewrite_snapshot(snap, lambda r: r.update(modeled_keys=[]))
    with pytest.raises(ValueError, match='SHA256'):
        load_snapshot(snap.manifest_path.parent, expected_sha256=snap.sha256)


@pytest.mark.parametrize('change', [lambda r: r['prepared']['atoms'][1].update(identity=r['prepared']['atoms'][2]['identity']),
    lambda r: r['prepared']['hydrogen_parents'].__setitem__(0, [4, [1]]),
    lambda r: r['prepared'].pop('hydrogen_parents'),
    lambda r: r['source_to_prepared'].__setitem__(0, 5),
    lambda r: r.update(modeled_keys=[]),
    lambda r: r.update(schema_version=True),
    lambda r: r.update(chemistry_validated=0),
    lambda r: r['parameter_files'][0].update(path='../escaped.xml')])
def test_reader_rejects_semantically_invalid_self_consistent_manifests(case, change):
    snap = snapshot(case); current_hash = rewrite_snapshot(snap, change)
    with pytest.raises(ValueError):
        load_snapshot(snap.manifest_path.parent, expected_sha256=current_hash)


@pytest.mark.parametrize('array', [np.full((11, 3), np.nan), np.zeros((10, 3)),
                                 np.zeros((11, 3), complex), np.zeros((11, 3), np.float32)])
def test_reader_rejects_wrong_arrays_even_with_matching_hashes(case, array):
    snap = snapshot(case); data = array_bytes(array)
    (snap.manifest_path.parent / 'coordinates/prepared.npy').write_bytes(data)
    def change(record):
        record['coordinates']['prepared'].update(sha256=hashlib.sha256(data).hexdigest(), size=len(data))
    current_hash = rewrite_snapshot(snap, change)
    with pytest.raises(ValueError):
        load_snapshot(snap.manifest_path.parent, expected_sha256=current_hash)


@pytest.mark.parametrize('symlink', [False, True])
def test_undeclared_or_symlink_payload_rejected(case, symlink):
    snap = snapshot(case); extra = snap.manifest_path.parent / 'extra'
    if symlink:
        extra.symlink_to(snap.manifest_path)
    else:
        extra.write_text('undeclared payload')
    with pytest.raises(ValueError):
        load_snapshot(snap.manifest_path.parent, expected_sha256=snap.sha256)


def test_accepted_bundle_binds_all_sidecars_and_remains_integrity_only(case):
    snap = snapshot(case); artifacts = bundle_inputs(case, snap)
    bundle = create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=artifacts,
                                    checks=dict.fromkeys(REQUIRED_CHECKS, True))
    restored = verify_accepted_bundle(bundle.manifest_path, expected_sha256=bundle.sha256)
    assert set(restored.artifacts) == {'coordinates', 'topology', 'parameters', 'provenance'}
    assert {f.name for f in restored.artifacts['provenance']} == {'validation.json', 'atom-map.json'}
    assert restored.snapshot.sha256 == snap.sha256
    data = json.loads(bundle.manifest_path.read_text())
    assert data['chemistry_validated'] is False and data['physical_model_validated'] is False
    with pytest.raises(TypeError):
        restored.validator_checks['geometry'] = False


@pytest.mark.parametrize('role', ['coordinates', 'topology', 'parameters', 'provenance'])
def test_handoff_refuses_changed_role_artifact(case, role):
    snap = snapshot(case)
    bundle = create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=bundle_inputs(case, snap),
                                    checks=dict.fromkeys(REQUIRED_CHECKS, True))
    bundle.artifacts[role][0].path.write_bytes(b'changed after acceptance')
    with pytest.raises(ValueError):
        verify_accepted_bundle(bundle.manifest_path, expected_sha256=bundle.sha256)


@pytest.mark.parametrize('value', [False, None, 1])
def test_missing_or_nontrue_scientific_check_cannot_create_bundle(case, value):
    snap = snapshot(case); checks = dict.fromkeys(REQUIRED_CHECKS, True); checks['geometry'] = value
    with pytest.raises(ValueError, match='validator'):
        create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=bundle_inputs(case, snap), checks=checks)


def test_changed_parameters_and_nonfinite_final_array_cannot_be_promoted(case):
    snap = snapshot(case); artifacts = bundle_inputs(case, snap)
    altered = case[0] / 'changed.xml'; altered.write_text('<ForceField/>different')
    artifacts['parameters']['main.xml'] = altered
    with pytest.raises(ValueError, match='parameter bytes differ'):
        create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=artifacts, checks=dict.fromkeys(REQUIRED_CHECKS, True))
    artifacts = bundle_inputs(case, snap)
    artifacts['coordinates']['final.npz'].write_bytes(array_bytes(np.full((11, 3), np.inf), npz=True))
    with pytest.raises(ValueError, match='Coordinates'):
        create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=artifacts, checks=dict.fromkeys(REQUIRED_CHECKS, True))


def test_existing_snapshot_and_bundle_are_never_overwritten(case):
    snap = snapshot(case)
    with pytest.raises(FileExistsError):
        create_snapshot(snap.manifest_path.parent, **case[1])
    artifacts = bundle_inputs(case, snap)
    bundle = create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=artifacts,
                                    checks=dict.fromkeys(REQUIRED_CHECKS, True))
    with pytest.raises(FileExistsError):
        create_accepted_bundle(bundle.manifest_path.parent, snapshot=snap, artifacts=artifacts,
                                checks=dict.fromkeys(REQUIRED_CHECKS, True))
    assert verify_accepted_bundle(bundle.manifest_path, expected_sha256=bundle.sha256)


def test_forged_huge_array_shape_rejected_before_allocation(case):
    snap = snapshot(case)
    stream = io.BytesIO()
    np.lib.format.write_array_header_1_0(stream, {'descr': '<f8', 'fortran_order': False, 'shape': (10**12, 3)})
    data = stream.getvalue()
    (snap.manifest_path.parent / 'coordinates/prepared.npy').write_bytes(data)
    def change(record):
        record['coordinates']['prepared'].update(sha256=hashlib.sha256(data).hexdigest(), size=len(data))
    current_hash = rewrite_snapshot(snap, change)
    with pytest.raises(ValueError, match='exact atom-count'):
        load_snapshot(snap.manifest_path.parent, expected_sha256=current_hash)


def test_final_npz_extra_array_and_missing_bundle_role_rejected(case):
    snap = snapshot(case); artifacts = bundle_inputs(case, snap)
    stream = io.BytesIO(); np.savez(stream, xyz_nm=snap.prepared_xyz_nm, unbound=np.ones(1))
    artifacts['coordinates']['final.npz'].write_bytes(stream.getvalue())
    with pytest.raises(ValueError, match='only xyz_nm'):
        create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=artifacts, checks=dict.fromkeys(REQUIRED_CHECKS, True))
    del artifacts['provenance']
    with pytest.raises(ValueError, match='roles are mandatory'):
        create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=artifacts, checks=dict.fromkeys(REQUIRED_CHECKS, True))


def test_same_bytes_symlink_replacement_is_not_a_valid_handoff(case):
    snap = snapshot(case); artifacts = bundle_inputs(case, snap)
    bundle = create_accepted_bundle(case[0] / 'accepted', snapshot=snap, artifacts=artifacts,
                                    checks=dict.fromkeys(REQUIRED_CHECKS, True))
    path = bundle.artifacts['topology'][0].path
    path.unlink(); path.symlink_to(artifacts['topology']['prepared.cif'])
    with pytest.raises(ValueError, match='Symlink'):
        verify_accepted_bundle(bundle.manifest_path, expected_sha256=bundle.sha256)
