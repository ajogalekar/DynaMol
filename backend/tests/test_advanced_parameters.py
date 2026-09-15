"""Native artifact continuity checks, not new chemistry/force-field validation."""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from backend.advanced_parameters import create_snapshot, copy_snapshot, read_snapshot, load_native_snapshot

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = Path(__file__).parent / "fixtures/ligand-6dbk"


@pytest.fixture
def snapshot(tmp_path):
    # This existing, independently generated native ligand is solely a transport
    # fixture. It does not claim advanced chemistry from its wrapper metadata.
    spec = importlib.util.spec_from_file_location("native_manifest_exporter", ROOT / "docs/audit/advanced-chemistry/export_native_parameter_manifest.py")
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    expected = tmp_path / "expected.json"
    exporter.export(FIXTURE / "ligand.prmtop", expected)
    state = tmp_path / "state.json"
    state.write_text(json.dumps({"purpose": "Native artifact transport test only", "model_accuracy_claimed": False}))
    folder = tmp_path / "snapshot"
    record = create_snapshot({"native_topology": FIXTURE / "ligand.prmtop",
                              "coordinates": FIXTURE / "ligand.inpcrd",
                              "expected_terms": expected, "chemical_state": state,
                              "source_structure": FIXTURE / "ligand.inpcrd"},
                             folder, model_id="transport-fixture", kind="covalent",
                             provenance={"fixture_only": True})
    return folder, record


def test_native_snapshot_survives_copy_with_exact_system_and_coordinates(snapshot, tmp_path):
    from openmm import XmlSerializer, unit
    source, record = snapshot
    first = load_native_snapshot(source, expected_manifest_sha256=record["manifest_sha256"])
    copy = tmp_path / "reopened"
    copy_snapshot(source, copy, expected_manifest_sha256=record["manifest_sha256"])
    second = load_native_snapshot(copy, expected_manifest_sha256=record["manifest_sha256"])
    assert first.mechanics_report["accepted"] and second.mechanics_report["accepted"]
    assert XmlSerializer.serialize(first.system) == XmlSerializer.serialize(second.system)
    assert np.array_equal(first.positions.value_in_unit(unit.nanometer), second.positions.value_in_unit(unit.nanometer))
    assert not first.simulation_ready and not second.simulation_ready
    assert first.chemical_state == second.chemical_state


@pytest.mark.parametrize("role", ["native_topology", "coordinates", "expected_terms", "chemical_state", "source_structure"])
def test_modified_artifact_fails_before_native_loading(snapshot, role):
    folder, record = snapshot
    contents = read_snapshot(folder)
    path = folder / contents["manifest"]["files"][role]["path"]
    path.write_bytes(path.read_bytes() + b"\nchanged")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_native_snapshot(folder, expected_manifest_sha256=record["manifest_sha256"])


def test_rehashed_false_native_parameters_fail_mechanics_comparison(snapshot):
    folder, _ = snapshot
    record = read_snapshot(folder)
    manifest = record["manifest"]
    entry = manifest["files"]["expected_terms"]
    expected = json.loads(record["contents"]["expected_terms"])
    expected["bonds"][0]["k_kj_mol_nm2"] *= 0.5
    data = json.dumps(expected).encode()
    (folder / entry["path"]).write_bytes(data)
    entry["sha256"] = hashlib.sha256(data).hexdigest()
    serialized = json.dumps(manifest).encode()
    (folder / "snapshot.json").write_bytes(serialized)
    with pytest.raises(ValueError, match="mechanics validation"):
        load_native_snapshot(folder, expected_manifest_sha256=hashlib.sha256(serialized).hexdigest())


def test_manifest_edit_cannot_promote_prototype_to_ready(snapshot):
    folder, _ = snapshot
    manifest = read_snapshot(folder)["manifest"]
    manifest["simulation_ready"] = True
    (folder / "snapshot.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="research prototypes only"):
        read_snapshot(folder)


def test_snapshot_cannot_overwrite_existing_result(snapshot, tmp_path):
    folder, record = snapshot
    with pytest.raises(ValueError, match="replace"):
        copy_snapshot(folder, folder, expected_manifest_sha256=record["manifest_sha256"])
    assert read_snapshot(folder)["manifest_sha256"] == record["manifest_sha256"]


@pytest.mark.parametrize("path", ["../elsewhere", "/tmp/elsewhere", "nested/file", "01/../file", "..\\elsewhere"])
def test_snapshot_rejects_escaping_paths(snapshot, path):
    folder, _ = snapshot
    manifest = read_snapshot(folder)["manifest"]
    manifest["files"]["coordinates"]["path"] = path
    (folder / "snapshot.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="relative filenames"):
        read_snapshot(folder)


def test_snapshot_rejects_symlinked_artifact(snapshot, tmp_path):
    folder, record = snapshot
    data = read_snapshot(folder)
    target = folder / data["manifest"]["files"]["native_topology"]["path"]
    outside = tmp_path / "outside.prmtop"
    outside.write_bytes(target.read_bytes())
    target.unlink()
    target.symlink_to(outside)
    with pytest.raises(ValueError, match="symlink"):
        load_native_snapshot(folder, expected_manifest_sha256=record["manifest_sha256"])


@pytest.mark.parametrize('digest',[None,'','0'*63,True])
def test_native_loading_and_copying_require_a_pinned_manifest_digest(snapshot,tmp_path,digest):
    folder,_=snapshot
    with pytest.raises(ValueError,match='explicit lowercase SHA-256'):
        load_native_snapshot(folder,expected_manifest_sha256=digest)
    with pytest.raises(ValueError,match='explicit lowercase SHA-256'):
        copy_snapshot(folder,tmp_path/'copy',expected_manifest_sha256=digest)
    assert not (tmp_path/'copy').exists()


@pytest.mark.parametrize('index',[False,0.0])
def test_rehashed_boolean_or_float_native_indices_are_not_integer_identity_maps(snapshot,index):
    folder,_=snapshot
    record=read_snapshot(folder);manifest=record['manifest'];entry=manifest['files']['expected_terms']
    terms=json.loads(record['contents']['expected_terms']);terms['atoms'][0]['native_index']=index
    data=json.dumps(terms).encode();(folder/entry['path']).write_bytes(data);entry['sha256']=hashlib.sha256(data).hexdigest()
    raw=json.dumps(manifest).encode();(folder/'snapshot.json').write_bytes(raw)
    with pytest.raises(ValueError,match='complete native atom index map'):
        load_native_snapshot(folder,expected_manifest_sha256=hashlib.sha256(raw).hexdigest())


def test_malformed_manifest_object_is_rejected_before_native_parsing(snapshot):
    folder,_=snapshot;(folder/'snapshot.json').write_text('[]')
    with pytest.raises(ValueError,match='research prototypes'):
        read_snapshot(folder)
