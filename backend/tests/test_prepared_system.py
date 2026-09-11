"""Parameter-file integrity and state continuity, not a force-field benchmark."""
import hashlib
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import mdtraj as md
import numpy as np
import pytest

from backend import config, jobs, storage
from backend.models import SimulationConfig
from backend.prepared_system import copy_ligand_parameters, ligand_parameter_files, load_prepared_forcefield


# Synthetic methane parameters test actual template loading without native
# charge calculation. These are test data, not parameters offered to users.
METHANE_XML = """<ForceField>
<AtomTypes><Type name="test-C" class="test-C" element="C" mass="12.011"/><Type name="test-H" class="test-H" element="H" mass="1.008"/></AtomTypes>
<Residues><Residue name="LIG"><Atom name="C1" type="test-C" charge="-0.4"/>
<Atom name="H1" type="test-H" charge="0.1"/><Atom name="H2" type="test-H" charge="0.1"/><Atom name="H3" type="test-H" charge="0.1"/><Atom name="H4" type="test-H" charge="0.1"/>
<Bond atomName1="C1" atomName2="H1"/><Bond atomName1="C1" atomName2="H2"/><Bond atomName1="C1" atomName2="H3"/><Bond atomName1="C1" atomName2="H4"/>
</Residue></Residues>
<HarmonicBondForce><Bond class1="test-C" class2="test-H" length="0.109" k="284512"/></HarmonicBondForce>
<HarmonicAngleForce><Angle class1="test-H" class2="test-C" class3="test-H" angle="1.910633" k="276"/></HarmonicAngleForce>
<NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5"><UseAttributeFromResidue name="charge"/>
<Atom type="test-C" sigma="0.34" epsilon="0.45"/><Atom type="test-H" sigma="0.25" epsilon="0.12"/>
</NonbondedForce></ForceField>"""


def bundle(folder, xml=METHANE_XML):
    path = Path(folder) / "ligands" / "LIG" / "parameters.xml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(xml)
    (path.parent / "native.log").write_text("Synthetic continuity fixture; no native calculation claimed.\n")
    return {"ph": 7, "simulation_ready": True, "ligand_parameters": {
        "forcefield": "GAFF2", "charge_method": "AM1-BCC", "requires_explicit_solvent": True,
        "ligands": [{"name": "LIG", "test_fixture": True}],
        "files": [{"path": "ligands/LIG/parameters.xml", "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}],
    }}


@pytest.fixture
def isolated_data(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(config, "DATASETS_DIR", tmp_path / "datasets")
    monkeypatch.setattr(config, "JOBS_DIR", tmp_path / "jobs")
    config.DATASETS_DIR.mkdir()
    config.JOBS_DIR.mkdir()
    return tmp_path


def test_bundle_copies_native_artifacts_and_rejects_modified_parameters(tmp_path):
    state = bundle(tmp_path / "source")
    copy_ligand_parameters(tmp_path / "source", tmp_path / "destination", state)
    assert (tmp_path / "destination/ligands/LIG/native.log").read_bytes() == (tmp_path / "source/ligands/LIG/native.log").read_bytes()
    assert len(ligand_parameter_files(tmp_path / "destination", state)) == 1
    (tmp_path / "destination/ligands/LIG/parameters.xml").write_text(METHANE_XML.replace('charge="-0.4"', 'charge="0.4"'))
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_prepared_forcefield(tmp_path / "destination", state)


@pytest.mark.parametrize("path", ["../outside.xml", "/tmp/outside.xml", "ligands/../../outside.xml", "ligands\\outside.xml", "ligands//outside.xml", "ligands/./outside.xml", "ligands/outside.json"])
def test_parameter_paths_cannot_escape_bundle(tmp_path, path):
    state = bundle(tmp_path)
    state["ligand_parameters"]["files"][0]["path"] = path
    with pytest.raises(ValueError, match="path"):
        ligand_parameter_files(tmp_path, state)


@pytest.mark.parametrize("content", ["<Include file='outside.xml'/>", "<Script>raise RuntimeError('should not run')</Script>", "<InitializationScript>print('should not run')</InitializationScript>"])
def test_parameter_xml_cannot_execute_or_include_unverified_content(tmp_path, content):
    state = bundle(tmp_path, "<ForceField>" + content + "</ForceField>")
    with pytest.raises(ValueError, match="self-contained"):
        load_prepared_forcefield(tmp_path, state)


def test_parameter_archive_symlinks_are_rejected(tmp_path):
    source = tmp_path / "source"
    state = bundle(source)
    (source / "ligands/linked.log").symlink_to(tmp_path / "outside.log")
    with pytest.raises(ValueError, match="symlink"):
        copy_ligand_parameters(source, tmp_path / "destination", state)
    assert not (tmp_path / "destination/ligands").exists()


def test_missing_bundle_and_implicit_ligand_model_fail_clearly(tmp_path):
    assert ligand_parameter_files(tmp_path, None) == []
    with pytest.raises(ValueError, match="no verified ligand"):
        ligand_parameter_files(tmp_path, None, required=True)
    state = bundle(tmp_path)
    with pytest.raises(ValueError, match="explicit TIP3P"):
        load_prepared_forcefield(tmp_path, state, solvent="implicit")
    (tmp_path / "ligands/LIG/parameters.xml").unlink()
    with pytest.raises(ValueError, match="missing"):
        load_prepared_forcefield(tmp_path, state)


def prepared_complex():
    from openmm import app, unit, Vec3

    pdb = app.PDBFile(str(config.ROOT / "examples/demo/topology.pdb"))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    top = app.Topology()
    residue = top.addResidue("LIG", top.addChain("Z"), id="1")
    carbon = top.addAtom("C1", app.element.carbon, residue)
    for index in range(1, 5):
        hydrogen = top.addAtom(f"H{index}", app.element.hydrogen, residue)
        top.addBond(carbon, hydrogen)
    center = np.asarray(pdb.positions.value_in_unit(unit.nanometer)).max(axis=0) + 0.4
    offsets = np.array([[0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) * 0.109 / np.sqrt(3)
    modeller.add(top, [Vec3(*row) for row in center + offsets] * unit.nanometer)
    traj = md.Trajectory(np.asarray(modeller.positions.value_in_unit(unit.nanometer))[None], md.Topology.from_openmm(modeller.topology), time=[0])
    metadata = storage.save_dataset(traj, "Synthetic parameter continuity fixture", "test", "No force-field accuracy claim.")
    folder = storage.dataset_dir(metadata["id"])
    metadata["preparation"] = bundle(folder)
    with (folder / "prepared.pdb").open("w") as output:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, output, keepIds=True)
    storage.atomic_json(folder / "metadata.json", metadata)
    return metadata


def assert_ligand_charge(folder, preparation, pdb_name="prepared.pdb"):
    import openmm as mm
    from openmm import app, unit

    pdb = app.PDBFile(str(folder / pdb_name))
    ff, files = load_prepared_forcefield(folder, preparation)
    assert files[-1] == "ligands/LIG/parameters.xml"
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)
    force = next(force for force in system.getForces() if isinstance(force, mm.NonbondedForce))
    ligand = [atom for atom in pdb.topology.atoms() if atom.residue.name == "LIG"]
    assert len(ligand) == 5
    charges = [force.getParticleParameters(atom.index)[0].value_in_unit(unit.elementary_charge) for atom in ligand]
    assert charges == pytest.approx([-0.4, 0.1, 0.1, 0.1, 0.1])


def test_real_solvation_job_and_output_preserve_parameter_bundle(isolated_data, monkeypatch):
    from backend.solvent import solvate_dataset
    from backend.worker import Worker

    parent = prepared_complex()
    parent_xyz = storage.load_physical(parent["id"]).xyz[0].copy()
    with pytest.raises(ValueError, match="explicit TIP3P"):
        jobs.submit_job(SimulationConfig(dataset_id=parent["id"], solvent="implicit"))
    preview = solvate_dataset(parent["id"], padding_nm=1, seed=29)
    assert preview["preparation"] == parent["preparation"]
    assert preview["solvation"]["forcefield"][-1] == "ligands/LIG/parameters.xml"
    assert np.allclose(storage.load_physical(preview["id"]).xyz[0, :len(parent_xyz)], parent_xyz, atol=0.000051)
    assert_ligand_charge(storage.dataset_dir(preview["id"]), preview["preparation"])
    with pytest.raises(ValueError, match="GROMACS adapter"):
        jobs.submit_job(SimulationConfig(dataset_id=preview["id"], engine="gromacs", solvent="explicit"))
    monkeypatch.setattr(jobs, "health", lambda: {"engines": [{"id": "openmm", "available": True}]})
    monkeypatch.setattr(jobs.subprocess, "Popen", lambda *args, **kwargs: SimpleNamespace(pid=424242))
    job = jobs.submit_job(SimulationConfig(dataset_id=preview["id"], solvent="explicit", duration_ps=0.002, report_interval=1, minimize=False, equilibration_steps=0))
    job_folder = config.JOBS_DIR / job["id"]
    input_state = json.loads((job_folder / "input-state.json").read_text())
    assert input_state["preparation"] == parent["preparation"]
    assert_ligand_charge(job_folder, input_state["preparation"], "input.pdb")
    # Test persistence without launching another expensive dynamics test.
    worker = Worker(job["id"])
    trajectory = storage.load_physical(preview["id"])
    shutil.copyfile(job_folder / "input.pdb", job_folder / "prepared.pdb")
    monkeypatch.setattr(worker, "run_openmm", lambda: trajectory)
    worker.run()
    assert worker.job["status"] == "completed", worker.job
    output = storage.get_dataset(worker.job["dataset_id"])
    assert output["preparation"] == parent["preparation"]
    assert output["solvation"] == preview["solvation"]
    assert_ligand_charge(storage.dataset_dir(output["id"]), output["preparation"])
    assert "ligands/LIG/parameters.xml" in worker.provenance["outputs"]
    assert "ligands/LIG/native.log" in worker.provenance["outputs"]
    assert (storage.dataset_dir(output["id"]) / "ligands/LIG/native.log").is_file()
