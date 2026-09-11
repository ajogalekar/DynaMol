"""Independent analytic fixtures for the public DynaMol storage/analysis boundary.

Run from repository root: .venv/bin/python docs/audit/check_numerics.py
Temporary structures are synthetic numerical fixtures, not simulated molecules.
"""
import datetime
import importlib.metadata
import json
import pathlib
import sys
import tempfile

import mdtraj as md
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from backend import config
from backend.analysis import measure
from backend.models import MeasurementRequest
from backend.storage import load_physical, load_uploaded, save_dataset


def topology():
    top = md.Topology()
    chain = top.add_chain("A")
    res = top.add_residue("LIG", chain, resSeq=1)
    atoms = [top.add_atom(n, e, res) for n, e in zip(
        ["N1", "H1", "O1", "C1"],
        [md.element.nitrogen, md.element.hydrogen, md.element.oxygen, md.element.carbon],
    )]
    top.add_bond(atoms[0], atoms[1])
    return top


def close(actual, expected, tolerance=1e-4):
    if not np.allclose(actual, expected, atol=tolerance, rtol=0):
        raise AssertionError(f"Expected {expected!r}; received {actual!r}")


def must_reject(call):
    try:
        call()
    except ValueError as exc:
        return str(exc)
    raise AssertionError("Invalid input was accepted")


results = []


def check(name, function):
    try:
        evidence = function()
        results.append({"name": name, "status": "pass", "evidence": evidence})
    except Exception as exc:
        results.append({"name": name, "status": "fail", "error": repr(exc)})


with tempfile.TemporaryDirectory(prefix="dynamol-review-") as temporary:
    workspace = pathlib.Path(temporary)
    config.DATASETS_DIR = workspace / "datasets"
    config.DATASETS_DIR.mkdir()

    def dataset(xyz, ident, periodic=False):
        coordinates = np.asarray(xyz, dtype=np.float32)
        if coordinates.ndim == 2:
            coordinates = coordinates[None, :, :]
        trajectory = md.Trajectory(coordinates, topology(), time=np.arange(len(coordinates)) * .25)
        if periodic:
            trajectory.unitcell_lengths = np.ones((len(coordinates), 3))
            trajectory.unitcell_angles = np.full((len(coordinates), 3), 90.)
        return save_dataset(trajectory, ident, "analytic fixture", "Not MD", dataset_id=ident)

    def measurement(ident, kind, indices):
        return measure(ident, MeasurementRequest(kind=kind, atoms=indices))

    def coordinate_units():
        meta = dataset([[0, .1, 0], [0, 0, 0], [.1, 0, 0], [.1, 0, .1]], "units")
        physical = load_physical("units")
        close(physical.xyz[0, 2] - physical.xyz[0, 1], [.1, 0, 0])
        display = np.fromfile(config.DATASETS_DIR / "units" / "coordinates.bin", dtype="<f4").reshape(1, 4, 3)
        close(display[0, 2] - display[0, 1], [1, 0, 0])
        assert [a["name"] for a in meta["atoms"]] == ["N1", "H1", "O1", "C1"]
        assert [a.name for a in physical.topology.atoms] == ["N1", "H1", "O1", "C1"]
        return {"stored_distance_nm": .1, "display_distance_angstrom": 1, "atom_order": "preserved"}

    check("nm to Angstrom boundary and canonical atom order", coordinate_units)

    def nonperiodic_geometry():
        distance = measurement("units", "distance", [0, 1])
        angle = measurement("units", "angle", [0, 1, 2])
        torsion = measurement("units", "dihedral", [0, 1, 2, 3])
        close(distance["values"], [1.])
        close(angle["values"], [90.])
        close(torsion["values"], [90.])
        return {"distance_angstrom": distance["values"], "angle_degrees": angle["values"], "dihedral_degrees": torsion["values"]}

    check("analytic Cartesian distance, angle and signed dihedral", nonperiodic_geometry)

    def periodic_geometry():
        dataset([[.95, 0, 0], [.05, 0, 0], [.05, .1, 0], [.05, .1, .1]], "periodic", True)
        d = measurement("periodic", "distance", [0, 1])
        a = measurement("periodic", "angle", [0, 1, 2])
        close(d["values"], [1.])
        close(a["values"], [90.])
        skew = md.Trajectory(np.asarray([[[0, 0, 0], [.45, .9, 0], [.1, .1, 0], [.1, .1, .1]]], dtype=np.float32), topology(), time=[0.])
        skew.unitcell_vectors = np.asarray([[[1, 0, 0], [.5, 1, 0], [0, 0, 1]]], dtype=np.float32)
        save_dataset(skew, "triclinic", "analytic fixture", "Not MD", dataset_id="triclinic")
        triclinic = measurement("triclinic", "distance", [0, 1])
        close(triclinic["values"], [np.sqrt(1.25)])
        return {"raw_separation_angstrom": 9, "minimum_image_angstrom": d["values"], "angle_degrees": a["values"], "triclinic_minimum_image_angstrom": triclinic["values"]}

    check("periodic minimum-image distance and angle", periodic_geometry)

    def hydrogen_bond():
        dataset([
            [[0, 0, 0], [.1, 0, 0], [.3, 0, 0], [.2, .2, .2]],
            [[0, 0, 0], [.1, 0, 0], [.4, 0, 0], [.2, .2, .2]],
        ], "hydrogen")
        result = measurement("hydrogen", "hbond", [0, 1, 2])
        close(result["values"], [3., 4.])
        close(result["angle_values"], [180., 180.])
        close(result["occupancy"], .5)
        rejection = must_reject(lambda: measurement("hydrogen", "hbond", [2, 1, 0]))
        return {"distances_angstrom": result["values"], "occupancy": result["occupancy"], "unbonded_donor_rejected": rejection}

    check("explicit donor-H connectivity and geometric occupancy", hydrogen_bond)

    def degeneracy():
        dataset([[0, 0, 0], [.1, 0, 0], [.2, 0, 0], [.3, 0, 0]], "collinear")
        return {"dihedral": must_reject(lambda: measurement("collinear", "dihedral", [0, 1, 2, 3])),
                "duplicate_indices": must_reject(lambda: measurement("collinear", "angle", [0, 1, 1])),
                "outside_topology": must_reject(lambda: measurement("collinear", "distance", [0, 4]))}

    check("undefined geometry and invalid indices rejected", degeneracy)

    def import_formats():
        trajectory = load_physical("units")
        trajectory = md.join([trajectory, trajectory])
        trajectory.time = np.array([0., .25])
        pdb = workspace / "topology.pdb"
        trajectory[0].save_pdb(str(pdb))
        evidence = {}
        for suffix in ["pdb", "gro", "xtc", "trr", "nc", "h5", "xyz", "dcd"]:
            path = workspace / f"trajectory.{suffix}"
            trajectory.save(str(path))
            independent = suffix in {"pdb", "gro", "h5"}
            loaded, warnings = load_uploaded(path if independent else pdb, None if independent else path, 1, None)
            close(loaded.xyz, trajectory.xyz, tolerance=1e-3)
            if suffix in {"xtc", "trr", "nc", "h5"}:
                close(loaded.time, [0., .25])
            if suffix == "dcd":
                assert any("time" in item.lower() for item in warnings), "DCD synthesized frame times need an explicit warning"
                explicit, _ = load_uploaded(pdb, path, 1, .25)
                close(explicit.time, [0., .25])
            evidence[suffix] = {"frames": loaded.n_frames, "coordinate_roundtrip": "pass", "times": loaded.time.tolist()}
        return evidence

    check("supported format round trips and DCD timestamp override", import_formats)

    def self_describing_identity():
        trajectory = load_physical("units")
        pdb = workspace / "identity-reference.pdb"
        trajectory[0].save_pdb(str(pdb))
        mismatch = workspace / "identity-mismatch.pdb"
        # Edit fixed PDB fields directly. Mutating MDTraj Atom names changes
        # their hashes and invalidates its internal bond lookup dictionary.
        lines = []
        for line in pdb.read_text().splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                name = line[12:16].strip()
                if name in {"N1", "O1"}:
                    replacement = "O1" if name == "N1" else "N1"
                    line = line[:12] + f"{replacement:>4}" + line[16:]
            lines.append(line)
        mismatch.write_text("\n".join(lines) + "\n")
        return {"same_count_reordered_names": must_reject(lambda: load_uploaded(pdb, mismatch, 1, .25))}

    check("self-describing trajectory identity mismatch rejected", self_describing_identity)

    def mmcif_import():
        from openmm.app import PDBxFile
        trajectory = load_physical("units")
        path = workspace / "structure.cif"
        with path.open("w") as handle:
            PDBxFile.writeFile(trajectory.topology.to_openmm(), trajectory.openmm_positions(0), handle)
        loaded, warnings = load_uploaded(path, None, 1, None)
        close(loaded.xyz, trajectory.xyz, tolerance=1e-4)
        return {"atoms": loaded.n_atoms, "coordinate_roundtrip": "pass", "warnings": warnings}

    check("OpenMM-written mmCIF import", mmcif_import)

report = {
    "recorded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "scope": "Analytic software fixtures only; these do not establish scientific validity or simulation convergence.",
    "packages": {name: importlib.metadata.version(name) for name in ["numpy", "mdtraj", "openmm"]},
    "checks": results,
    "passed": sum(item["status"] == "pass" for item in results),
    "failed": sum(item["status"] == "fail" for item in results),
}
output = pathlib.Path(__file__).with_name("numerical-checks.json")
output.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
raise SystemExit(1 if report["failed"] else 0)
