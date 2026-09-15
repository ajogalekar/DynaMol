"""Fresh, isolated release API audit; analytic fixtures are not simulation evidence."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import time
import traceback
import zipfile

ROOT = Path(__file__).resolve().parents[4]
DATA = ROOT / "data/integration-checks/final-release-2026-09-12/regression/api"
os.environ["DYNAMOL_DATA_DIR"] = str(DATA)
os.environ["OPENMM_CPU_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import mdtraj as md
import numpy as np
from fastapi.testclient import TestClient
from rdkit import Chem
from rdkit.Chem import AllChem

from backend import config, storage
from backend.main import app
from backend.tests.test_sources_solvent import ETHANOL_MOL2

HERE = Path(__file__).resolve().parent
FIXTURES = DATA / "fixtures"
FIXTURES.mkdir(parents=True, exist_ok=True)
parser = argparse.ArgumentParser()
parser.add_argument("--output", default="api-initial.json")
args = parser.parse_args()
report = {"started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "data_root": str(DATA),
          "method": "Actual FastAPI TestClient requests; generated analytic fixtures and live fixed-provider network fetches. No native MD or preparation jobs are submitted.",
          "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob("backend/**/*.py")},
          "checks": []}


def check(name, fn):
    started = time.monotonic()
    record = {"name": name}
    try:
        record["details"] = fn()
        record["passed"] = True
    except Exception as exc:
        record.update(passed=False, error=str(exc), traceback=traceback.format_exc())
    record["elapsed_seconds"] = round(time.monotonic() - started, 4)
    report["checks"].append(record)
    (HERE / args.output).write_text(json.dumps(report, indent=2))
    print(f"{'PASS' if record['passed'] else 'FAIL'} {name}: {record.get('error', '')}", flush=True)


def expect(response, status=200):
    assert response.status_code == status, (response.status_code, response.text[:1000])
    return response.json()


def topology():
    top = md.Topology()
    residue = top.add_residue("TST", top.add_chain("A"), resSeq=1)
    elements = [md.element.nitrogen, md.element.hydrogen, md.element.oxygen, md.element.carbon]
    atoms = [top.add_atom(name, e, residue) for name, e in zip(["N", "H", "O", "C"], elements)]
    top.add_bond(atoms[0], atoms[1])
    return top


def make(xyz, periodic=False):
    trajectory = md.Trajectory(np.asarray(xyz, dtype=np.float32), topology(), time=np.arange(len(xyz)) * .25)
    if periodic:
        trajectory.unitcell_lengths = np.ones((len(xyz), 3))
        trajectory.unitcell_angles = np.ones((len(xyz), 3)) * 90
    return trajectory


base = make([[[0, 0, 0], [.1, 0, 0], [.1, .1, 0], [.1, .1, .1]],
             [[0, 0, 0], [.2, 0, 0], [.2, .1, 0], [.2, .1, -.1]],
             [[0, 0, 0], [.3, 0, 0], [.3, .1, 0], [.3, .1, .1]]], periodic=True)
base[0].save_pdb(str(FIXTURES / "topology.pdb"))
base.save_cif(str(FIXTURES / "sample.cif"))
base[0].save_gro(str(FIXTURES / "sample.gro"))


with TestClient(app, base_url="http://127.0.0.1", raise_server_exceptions=False) as client:
    def validate_dataset(meta):
        assert meta["n_atoms"] > 0 and meta["n_frames"] > 0
        assert len(meta["atoms"]) == meta["n_atoms"]
        coordinates = client.get(meta["coordinates_url"])
        assert coordinates.status_code == 200
        assert len(coordinates.content) == meta["n_atoms"] * meta["n_frames"] * 12
        assert np.isfinite(np.frombuffer(coordinates.content, dtype="<f4")).all()
        assert client.get(meta["topology_url"]).status_code == 200
        physical = storage.load_physical(meta["id"])
        assert (physical.n_atoms, physical.n_frames) == (meta["n_atoms"], meta["n_frames"])
        return {"dataset_id": meta["id"], "atoms": meta["n_atoms"], "frames": meta["n_frames"], "time_unit": meta["time_unit"]}

    def source_import(path):
        original = path.read_bytes()
        meta = expect(client.post("/api/structures/upload", files={"file": (path.name, original)}))
        folder = storage.dataset_dir(meta["id"])
        provenance = json.loads((folder / "provenance.json").read_text())
        assert (folder / meta["source_file"]).read_bytes() == original
        assert provenance["source_sha256"] == hashlib.sha256(original).hexdigest()
        return validate_dataset(meta)

    structure_paths = [ROOT / "examples/ubiquitin-start/topology.pdb", FIXTURES / "sample.cif"]
    for suffix in ("mmcif", "pdbx"):
        target = FIXTURES / f"sample.{suffix}"
        shutil.copy2(FIXTURES / "sample.cif", target)
        structure_paths.append(target)
    (FIXTURES / "ethanol.mol2").write_text(ETHANOL_MOL2)
    structure_paths.append(FIXTURES / "ethanol.mol2")
    molecule = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)C(=O)[O-]"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=2026) == 0
    for suffix in ("mol", "sdf"):
        path = FIXTURES / f"lactate.{suffix}"
        path.write_text(Chem.MolToMolBlock(molecule) + ("\n$$$$\n" if suffix == "sdf" else ""))
        structure_paths.append(path)
    for suffix in ("smi", "smiles"):
        path = FIXTURES / f"lactate.{suffix}"
        path.write_text("C[C@H](O)C(=O)[O-] lactate\n")
        structure_paths.append(path)
    for path in structure_paths:
        check("structure-upload:" + path.suffix, lambda path=path: source_import(path))

    def smiles_roundtrip():
        meta = expect(client.post("/api/structures/smiles", json={"smiles": "C[C@H](O)C(=O)[O-]", "seed": 2026}))
        chemistry = json.loads((storage.dataset_dir(meta["id"]) / "chemistry.json").read_text())
        assert chemistry["total_formal_charge"] == -1 and "@" in chemistry["canonical_isomeric_smiles"]
        return validate_dataset(meta)
    check("smiles:charge-and-stereochemistry", smiles_roundtrip)

    methods = {"pdb": "save_pdb", "xtc": "save_xtc", "trr": "save_trr", "dcd": "save_dcd", "h5": "save_hdf5",
               "nc": "save_netcdf", "mdcrd": "save_mdcrd", "lammpstrj": "save_lammpstrj", "xyz": "save_xyz"}
    paths = {}
    for suffix, method in methods.items():
        path = FIXTURES / f"frames.{suffix}"
        getattr(base, method)(str(path))
        paths[suffix] = path
    for alias, native in (("hdf5", "h5"), ("netcdf", "nc"), ("crd", "mdcrd")):
        path = FIXTURES / f"frames.{alias}"
        shutil.copy2(paths[native], path)
        paths[alias] = path

    def trajectory_import(path):
        meta = expect(client.post("/api/datasets/upload", files={"topology": ("topology.pdb", (FIXTURES / "topology.pdb").read_bytes()), "trajectory": (path.name, path.read_bytes())}, data={"stride": "2", "frame_interval_ps": ".25"}))
        actual = storage.load_physical(meta["id"])
        assert actual.n_frames == 2
        np.testing.assert_allclose(actual.xyz, base.xyz[::2], atol=2e-4)
        np.testing.assert_allclose(actual.time, [0, .5], atol=1e-7)
        folder = storage.dataset_dir(meta["id"]) / "originals"
        assert (folder / ("trajectory" + path.suffix)).read_bytes() == path.read_bytes()
        return validate_dataset(meta)
    for suffix, path in paths.items():
        check("trajectory-upload:" + suffix + ":stride-and-time", lambda path=path: trajectory_import(path))

    def topology_import(path):
        meta = expect(client.post("/api/datasets/upload", files={"topology": (path.name, path.read_bytes())}))
        return validate_dataset(meta)
    for path in [FIXTURES / "topology.pdb", FIXTURES / "sample.cif", FIXTURES / "sample.mmcif", FIXTURES / "sample.pdbx", FIXTURES / "sample.gro", paths["h5"], paths["hdf5"]]:
        check("topology-only:" + path.suffix, lambda path=path: topology_import(path))

    def invalid_upload(endpoint, field, filename, data):
        before = {p.name for p in config.DATASETS_DIR.iterdir()}
        response = client.post(endpoint, files={field: (filename, data)})
        result = expect(response, 422)
        assert isinstance(result["detail"], str) and result["detail"]
        assert {p.name for p in config.DATASETS_DIR.iterdir()} == before, "invalid import left a partial dataset"
        return {"status": response.status_code, "detail": result["detail"]}
    for endpoint, field in [("/api/structures/upload", "file"), ("/api/datasets/upload", "topology")]:
        for filename, data in [("empty.pdb", b""), ("broken.pdb", b"not a pdb\n"), ("broken.cif", b"not a CIF\n"), ("unknown.exe", b"arbitrary bytes")]:
            check("invalid:" + endpoint + ":" + filename, lambda ep=endpoint, field=field, filename=filename, data=data: invalid_upload(ep, field, filename, data))
    for filename, data in [("broken.sdf", b"bad\n"), ("broken.mol2", b"bad\n"), ("many.smi", b"CCO\nCCC\n"), ("nonutf8.smiles", b"\xff\xff")]:
        check("invalid:structure:" + filename, lambda filename=filename, data=data: invalid_upload("/api/structures/upload", "file", filename, data))

    def invalid_smiles(value):
        return expect(client.post("/api/structures/smiles", json={"smiles": value}), 422)
    for value in ("bad smiles", "C.C", "CCO\nCCC"):
        check("invalid:smiles:" + repr(value), lambda value=value: invalid_smiles(value))

    def analytic_measurements():
        meta = storage.save_dataset(base, "Analytic signed geometry", "Generated audit fixture", "Not simulation data")
        for kind, atoms, expected in [("distance", [0, 1], [1, 2, 3]), ("angle", [0, 1, 2], [90, 90, 90]), ("dihedral", [0, 1, 2, 3], [90, -90, 90])]:
            body = {"kind": kind, "atoms": atoms}
            result = expect(client.post(f"/api/datasets/{meta['id']}/measurements", json=body))
            np.testing.assert_allclose(result["values"], expected, atol=1e-4)
            preview = expect(client.post(f"/api/datasets/{meta['id']}/measurement-preview", json=body))
            np.testing.assert_allclose(preview["values"], expected, atol=1e-4)
            assert preview["frame_errors"] == [None] * 3
        bad = expect(client.post(f"/api/datasets/{meta['id']}/measurements", json={"kind": "distance", "atoms": [0, 0]}), 422)
        return {"dataset_id": meta["id"], "distance_angstrom": [1, 2, 3], "angle_degrees": [90] * 3, "signed_dihedral_degrees": [90, -90, 90], "invalid": bad}
    check("numeric:distance-angle-signed-dihedral-preview", analytic_measurements)

    def hbond_and_periodic():
        hbond = make([[[0, 0, 0], [.1, 0, 0], [.3, 0, 0], [0, .4, 0]], [[0, 0, 0], [.1, 0, 0], [.4, 0, 0], [0, .4, 0]]])
        meta = storage.save_dataset(hbond, "Analytic H-bond", "Audit fixture", "Geometry only")
        result = expect(client.post(f"/api/datasets/{meta['id']}/measurements", json={"kind": "hbond", "atoms": [0, 1, 2]}))
        assert result["occupancy"] == .5
        np.testing.assert_allclose(result["values"], [3, 4], atol=1e-5)
        np.testing.assert_allclose(result["angle_values"], [180, 180], atol=1e-5)
        pbc = make([[[.05, 0, 0], [.95, 0, 0], [0, .2, 0], [0, 0, .3]]], True)
        other = storage.save_dataset(pbc, "Periodic distance", "Audit fixture", "Analytic 1nm box")
        minimum = expect(client.post(f"/api/datasets/{other['id']}/measurements", json={"kind": "distance", "atoms": [0, 1]}))
        np.testing.assert_allclose(minimum["values"], [1], atol=1e-5)
        return {"hbond": result, "periodic_distance": minimum}
    check("numeric:hbond-occupancy-and-minimum-image", hbond_and_periodic)

    def rms_controls():
        shifted = base.xyz[0] + np.array([.2, 0, 0])
        data = make([base.xyz[0], shifted])
        meta = storage.save_dataset(data, "Rigid translation", "Audit fixture", "Known 2 angstrom shift")
        body = {"kind": "rmsd", "selection": "custom", "atoms": [0, 1, 2, 3], "alignment": "none", "periodic": "cartesian"}
        result = expect(client.post(f"/api/datasets/{meta['id']}/structural-analysis", json=body))
        np.testing.assert_allclose(result["values"], [0, 2], atol=1e-5)
        fitted = expect(client.post(f"/api/datasets/{meta['id']}/structural-analysis", json={**body, "alignment": "selection"}))
        np.testing.assert_allclose(fitted["values"], [0, 0], atol=1e-5)
        rmsf = expect(client.post(f"/api/datasets/{meta['id']}/structural-analysis", json={**body, "kind": "rmsf", "residue_average": False}))
        np.testing.assert_allclose(rmsf["values"], [1] * 4, atol=1e-5)
        return {"unfitted_rmsd_angstrom": result["values"], "fitted_rmsd_angstrom": fitted["values"], "unfitted_rmsf_angstrom": rmsf["values"]}
    check("numeric:RMSD-RMSF-independent-analytic-controls", rms_controls)

    def fetch(provider, identifier):
        meta = expect(client.post("/api/structures/fetch", json={"provider": provider, "identifier": identifier}))
        details = validate_dataset(meta)
        provenance = json.loads((storage.dataset_dir(meta["id"]) / "provenance.json").read_text())
        assert provenance["url"].startswith("https://") and len(provenance["source_sha256"]) == 64
        return {**details, "provenance": provenance}
    for provider, identifier in [("pdb", "1UBQ"), ("pdb", "pdb_00001ubq"), ("pubchem", "2244"), ("pubchem", "aspirin")]:
        check("live-fetch:" + provider + ":" + identifier, lambda provider=provider, identifier=identifier: fetch(provider, identifier))

    def workspace_roundtrip():
        meta = storage.save_dataset(base, "Portable geometry", "Audit fixture", "Known generated frames")
        values = expect(client.post(f"/api/datasets/{meta['id']}/measurements", json={"kind": "distance", "atoms": [0, 1]}))
        state = {"version": 1, "dataset_id": meta["id"], "frame": 2, "camera": np.eye(4).reshape(-1).tolist(), "selected_atoms": [0, 1],
                 "measurements": [{**values, "id": "test-distance", "label": "N-H", "color": "#73d8c0", "visible": True}], "active_measurement": "test-distance",
                 "named_selections": [{"id": "two-atoms", "name": "Selection", "atoms": [0, 1]}]}
        saved = expect(client.post("/api/workspace", json={"state": state}))
        assert len(saved["trajectory_signature"]) == 64
        restored = expect(client.get("/api/workspace"))["state"]
        assert restored["frame"] == 2 and restored["camera"] == state["camera"]
        assert restored["measurements"][0]["values"] == values["values"]
        assert restored["color_scheme"] == "element"
        project = expect(client.post("/api/projects", json={"name": "API portability audit", "state": state}))
        archive = client.get(f"/api/projects/{project['id']}/export")
        assert archive.status_code == 200 and archive.headers["content-type"] == "application/zip"
        (FIXTURES / "portable-project.zip").write_bytes(archive.content)
        with zipfile.ZipFile(io.BytesIO(archive.content)) as z:
            assert "dynamol-project.json" in z.namelist()
        old_root = (config.DATA_ROOT, config.DATASETS_DIR, config.JOBS_DIR)
        destination = DATA / "portable-destination"
        destination.mkdir(exist_ok=True)
        config.DATA_ROOT, config.DATASETS_DIR, config.JOBS_DIR = destination, destination / "datasets", destination / "jobs"
        config.DATASETS_DIR.mkdir(exist_ok=True); config.JOBS_DIR.mkdir(exist_ok=True)
        try:
            imported = expect(client.post("/api/projects/import", files={"file": ("portable.zip", archive.content, "application/zip")}))
            result = expect(client.get(f"/api/projects/{imported['project']['id']}"))
            assert result["state"]["color_scheme"] == "element"
            assert result["state"]["measurements"][0]["values"] == values["values"]
            np.testing.assert_array_equal(storage.load_physical(meta["id"]).xyz, base.xyz)
            assert not list(config.JOBS_DIR.iterdir())
        finally:
            config.DATA_ROOT, config.DATASETS_DIR, config.JOBS_DIR = old_root
        other = storage.save_dataset(base, "Unrelated geometry", "Audit fixture", "Fresh scene")
        expect(client.post("/api/workspace", json={"state": {"version": 1, "dataset_id": other["id"]}}))
        fresh = expect(client.get("/api/workspace"))["state"]
        assert fresh["color_scheme"] == "element" and fresh["measurements"] == [] and fresh["selected_atoms"] == []
        expect(client.post("/api/workspace", json={"state": {"version": 1, "dataset_id": other["id"], "color_scheme": "chain"}}))
        custom = expect(client.get("/api/workspace"))["state"]
        assert custom["color_scheme"] == "chain"
        return {"source_dataset": meta["id"], "project_id": project["id"], "archive_bytes": len(archive.content), "fresh_default": fresh["color_scheme"], "explicit_saved_color_preserved": custom["color_scheme"]}
    check("workspace:save-export-fresh-import-color-defaults-and-selection-identity", workspace_roundtrip)

    def stability_reads():
        datasets = expect(client.get("/api/datasets"))
        for i in range(50):
            assert client.get("/api/datasets").status_code == 200
            assert client.get("/api/workspace").status_code == 200
            assert client.get("/api/library").status_code == 200
            assert client.get("/api/jobs").json() == []
        return {"iterations": 50, "requests": 200, "datasets": len(datasets), "submitted_jobs": 0}
    check("stability:200-library-workspace-and-job-reads", stability_reads)

report["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
report["passed"] = sum(row["passed"] for row in report["checks"])
report["failed"] = sum(not row["passed"] for row in report["checks"])
(HERE / args.output).write_text(json.dumps(report, indent=2))
print(json.dumps({"passed": report["passed"], "failed": report["failed"]}), flush=True)
raise SystemExit(bool(report["failed"]))
