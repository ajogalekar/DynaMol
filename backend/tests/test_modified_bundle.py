"""A prepared modification exports its actual portable parameter files and hashes."""
import hashlib
import json
from pathlib import Path

import pytest

from backend.modified_residues import DATA, modified_forcefield_provenance
from backend.prepared_system import copy_ligand_parameters, load_prepared_forcefield, snapshot_modified_parameters


def state():
    return {'modified_residues': [{'chain': 'A', 'resid': '170', 'residue': 'TPO', 'formal_charge': -2}], 'modified_residue_parameters': modified_forcefield_provenance()}


def test_modified_parameter_snapshot_and_copy_work_without_any_ligand_bundle(tmp_path):
    source, destination, continuation = [tmp_path / name for name in ['prepared-job', 'solvated-dataset', 'md-job']]
    for folder in [source, destination, continuation]:
        folder.mkdir()
    preparation = state()
    snapshot_modified_parameters(source, preparation)
    expected = {'manifest.json', *json.loads((DATA / 'manifest.json').read_text())['files']}
    assert {path.name for path in (source / 'residue-parameters').iterdir()} == expected
    for name in expected:
        assert (source / 'residue-parameters' / name).read_bytes() == (DATA / name).read_bytes()
    for source_folder, target in [(source, destination), (destination, continuation)]:
        copy_ligand_parameters(source_folder, target, preparation)
        for name in expected:
            assert hashlib.sha256((target / 'residue-parameters' / name).read_bytes()).digest() == hashlib.sha256((source / 'residue-parameters' / name).read_bytes()).digest()
    assert not (continuation / 'ligands').exists()
    field, names = load_prepared_forcefield(continuation, preparation)
    assert 'TPO' in field._templates
    assert any('phosaa14SB.xml' in name for name in names)


def test_modified_parameter_copy_rejects_tampered_xml(tmp_path):
    source, destination = tmp_path / 'source', tmp_path / 'destination'
    source.mkdir(); destination.mkdir()
    preparation = state()
    snapshot_modified_parameters(source, preparation)
    path = source / 'residue-parameters' / 'phosaa14SB.xml'
    path.write_bytes(path.read_bytes() + b'\n<!-- changed -->\n')
    with pytest.raises(ValueError, match='checksum|hash|match|changed|mismatch'):
        copy_ligand_parameters(source, destination, preparation)


def test_modified_parameter_copy_rejects_symlinked_source_file(tmp_path):
    source, destination = tmp_path / 'source', tmp_path / 'destination'
    source.mkdir(); destination.mkdir()
    preparation = state()
    snapshot_modified_parameters(source, preparation)
    path = source / 'residue-parameters' / 'phosaa14SB.xml'
    external = tmp_path / 'outside.xml'
    external.write_bytes(path.read_bytes())
    path.unlink()
    path.symlink_to(external)
    with pytest.raises(ValueError, match='[Ss]ymlink|regular|escape'):
        copy_ligand_parameters(source, destination, preparation)


def test_modified_parameters_reject_version_mismatch_and_implicit_solvent(tmp_path):
    preparation = state()
    snapshot_modified_parameters(tmp_path, preparation)
    with pytest.raises(ValueError, match='explicit'):
        load_prepared_forcefield(tmp_path, preparation, solvent='implicit')
    preparation['modified_residue_parameters']['registry_version'] = -1
    with pytest.raises(ValueError, match='match|version|state'):
        load_prepared_forcefield(tmp_path, preparation)
