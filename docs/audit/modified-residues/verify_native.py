"""Compare bundled modified-residue support with native AmberTools24.8 fixtures."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
from openmm import app, Context, Platform, VerletIntegrator, unit
from backend.modified_residues import SUPPORTED_MODIFIED, modified_forcefield_files, modified_forcefield_provenance, register_modified_forcefield, register_topology_definitions

NATIVE = Path(__file__).parent / "native"


def state(system, positions):
    integrator = VerletIntegrator(.001)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    result = context.getState(getEnergy=True, getForces=True)
    return result.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), result.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)


def run():
    register_topology_definitions()
    rows = []
    for name in sorted(SUPPORTED_MODIFIED):
        pdb = app.PDBFile(str(NATIVE / f"{name}.pdb"))
        native = app.AmberPrmtopFile(str(NATIVE / f"{name}.prmtop"))
        xyz = np.asarray(app.AmberInpcrdFile(str(NATIVE / f"{name}.inpcrd")).positions.value_in_unit(unit.nanometer))
        forcefield = register_modified_forcefield(app.ForceField("amber14/protein.ff14SB.xml", *modified_forcefield_files()))
        options = {"nonbondedMethod": app.NoCutoff, "constraints": None, "rigidWater": False, "removeCMMotion": False}
        native_system, generated_system = native.createSystem(**options), forcefield.createSystem(pdb.topology, **options)
        rng = np.random.default_rng(2026)
        poses = []
        for index in range(3):
            coordinates = xyz if index == 0 else xyz + rng.normal(0, .002, xyz.shape)
            native_energy, native_forces = state(native_system, coordinates)
            energy, forces = state(generated_system, coordinates)
            poses.append({"pose": index, "native_energy_kj_mol": native_energy, "generated_energy_kj_mol": energy, "energy_difference_kj_mol": abs(native_energy - energy), "max_force_component_difference_kj_mol_nm": float(np.max(np.abs(native_forces - forces)))})
        passed = all(row["energy_difference_kj_mol"] < .002 and row["max_force_component_difference_kj_mol_nm"] < .02 for row in poses)
        rows.append({"residue": name, "atoms": native_system.getNumParticles(), "poses": poses, "passed": passed})
    return {"passed": all(row["passed"] for row in rows), "reference": "Native AmberTools24.8 TLEAP ACE-residue-NME systems, generated from leaprc.protein.ff14SB and leaprc.phosaa14SB", "platform": "OpenMM Reference", "registry": modified_forcefield_provenance(), "thresholds": {"energy_kj_mol": .002, "max_force_component_kj_mol_nm": .02}, "threshold_explanation": "Native TLEAP prmtop serializes angular constants with finite-precision pi; published XML uses mathematical pi. Differences are serialization-level numerical behavior, not uncertainty estimates or scientific accuracy validation.", "fixtures": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(NATIVE.iterdir()) if path.is_file()}, "rows": rows}


if __name__ == "__main__":
    report = run()
    (Path(__file__).parent / "native-parity.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"passed": report["passed"], "residues": len(report["rows"]), "poses": sum(len(row["poses"]) for row in report["rows"])}))
    raise SystemExit(0 if report["passed"] else 1)
