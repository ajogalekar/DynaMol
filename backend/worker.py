"""One independent, bounded CPU simulation process per job.

Run only through jobs.submit_job. All preparation, logs and actual output remain
in a persistent job directory for inspection and download.
"""
import copy
import hashlib
import importlib.metadata
import json
import math
import os
import re
import random
import signal
import shutil
import uuid
import subprocess
import sys
import time
import tomllib
import traceback
from pathlib import Path

from . import config
from .jobs import gromacs_executable
from .prepared_system import copy_ligand_parameters, load_prepared_forcefield
from .storage import atomic_json, save_dataset, safe_id, dataset_dir


_termination_requested = False


class Cancelled(Exception):
    pass


class Worker:
    def __init__(self, job_id: str):
        self.folder = config.JOBS_DIR / safe_id(job_id)
        for _ in range(50):
            self.job = json.loads((self.folder / "status.json").read_text())
            if self.job.get("worker_pid") == os.getpid():
                break
            time.sleep(0.02)
        self.settings = self.job["config"]
        state_path = self.folder / "input-state.json"
        self.input_state = json.loads(state_path.read_text()) if state_path.exists() else {}
        self.started = time.monotonic()
        self.previous_elapsed = self.job.get("elapsed_seconds", 0) if self.job.get("resume_requested") else 0
        application_version = tomllib.loads((config.ROOT / "pyproject.toml").read_text())["project"]["version"]
        self.provenance = {"application": f"DynaMol {application_version}", "purpose": "Short exploratory molecular dynamics/software demonstration, not converged scientific validation.", "config": self.settings, "input_sha256": hashlib.sha256((self.folder / "input.pdb").read_bytes()).hexdigest(), "cpu_threads": config.CPU_THREADS, "versions": {name: importlib.metadata.version(name) for name in ("openmm", "mdtraj", "numpy")}, "commands": [], "preparation": [], "worker_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "source_dataset_id": self.settings["dataset_id"]}

        if self.job.get("resume_requested") and (self.folder / "provenance.json").exists():
            self.provenance = json.loads((self.folder / "provenance.json").read_text())
            self.provenance.pop("failure", None)
        self.provenance["restarts"] = self.job.get("restarts", [])
        self.provenance["resource_limits"] = {"maximum_atoms": config.MAX_ATOMS, "maximum_coordinate_bytes": config.MAX_COORD_BYTES,
                                              "maximum_frames": config.MAX_FRAMES}

    def check_cancel(self):
        if _termination_requested or (self.folder / "cancel.request").exists():
            raise Cancelled("Cancelled by user.")

    def check_stereochemistry(self, topology, xyz_nm, stage, box_nm=None):
        from .stereo_monitor import StereoMonitor, require_valid
        if not hasattr(self, '_stereo_monitor'):
            self._stereo_monitor = StereoMonitor(topology, xyz_nm, self.folder / 'input.pdb')
        report = self._stereo_monitor.check(xyz_nm, box_nm, stage)
        history_path = self.folder / 'stereochemistry-checks.json'
        if not hasattr(self, '_stereo_checks'):
            # A new worker may be continuing an interrupted run. Keep the
            # preparation and earlier production evidence across that restart.
            self._stereo_checks = json.loads(history_path.read_text())['checks'] if history_path.exists() else []
        history = self._stereo_checks
        history.append({k: v for k, v in report.items() if k != 'method'})
        atomic_json(history_path, {'method': report['method'], 'checks': history})
        if not report['passed']:
            import numpy as np
            from openmm import app, unit
            evidence = self.folder / 'stereochemistry-failures' / uuid.uuid4().hex[:16]
            evidence.mkdir(parents=True)
            coordinates = evidence / 'coordinates.npz'
            topology_path = evidence / 'topology.cif'
            np.savez_compressed(coordinates, xyz_nm=xyz_nm,
                                **({'box_nm': box_nm} if box_nm is not None else {}))
            # Native preparation can add atoms before prepared.pdb exists. Save
            # its matching topology too. NaNs remain in the NPZ; finite zero
            # placeholders only make the topology file readable in that case.
            xyz = np.asarray(xyz_nm)
            nonfinite = int((~np.isfinite(xyz)).sum())
            with topology_path.open('w') as stream:
                app.PDBxFile.writeFile(topology, np.where(np.isfinite(xyz), xyz, 0) * unit.nanometer, stream, keepIds=True)
            report['artifacts'] = {str(path.relative_to(self.folder)): hashlib.sha256(path.read_bytes()).hexdigest()
                                   for path in (coordinates, topology_path)}
            report['topology_coordinate_placeholders'] = nonfinite
            atomic_json(evidence / 'report.json', report)
            # Preserve the first failure's legacy path for existing download
            # consumers, while every failure has its own immutable evidence.
            legacy = self.folder / 'stereochemistry-failure.npz'
            if not legacy.exists():
                shutil.copy2(coordinates, legacy)
            self.provenance['stereochemistry_failure'] = report
            atomic_json(self.folder / 'provenance.json', self.provenance)
        require_valid(report)

    def update(self, stage=None, completed=None, message=None, **fields):
        self.check_cancel()
        if stage is not None:
            self.job["stage"] = stage
        if completed is not None:
            self.job["completed_steps"] = int(completed)
            self.job["progress"] = min(100, 100 * completed / self.job["total_steps"])
        if message:
            self.job["logs"] = (self.job["logs"] + [message])[-150:]
            print(message, flush=True)
        self.job.update(fields)
        self.job["elapsed_seconds"] = round(self.previous_elapsed + time.monotonic() - self.started, 2)
        atomic_json(self.folder / "status.json", self.job)

    def record_preparation(self, message):
        self.provenance["preparation"].append(message)
        atomic_json(self.folder / "provenance.json", self.provenance)
        self.update(message=message)

    def verify_measurement_identity(self, topology, *, coordinate_tolerance=2e-4):
        from .live_measurements import verify_native_identity
        if not self.settings.get("measurements"):
            return
        proof = {"input_sha256": hashlib.sha256((self.folder / "input.pdb").read_bytes()).hexdigest(), "verified": False}
        try:
            mapping = verify_native_identity(self.folder, self.input_state, topology, coordinate_tolerance=coordinate_tolerance)
            proof.update(verified=True, input_to_native=mapping, coordinate_tolerance_nm=coordinate_tolerance)
        except Exception as exc:
            proof["error"] = str(exc)
        atomic_json(self.folder / "measurement-identity.json", proof)

    def start_measurements(self, topology):
        from .live_measurements import LiveMeasurements, empty_snapshot
        self.live_measurements = None
        try:
            state = dict(self.input_state)
            if self.settings.get("measurements"):
                proof = json.loads((self.folder / "measurement-identity.json").read_text())
                if not proof.get("verified"):
                    raise ValueError(proof.get("error", "Tracked native atom identities could not be verified."))
                if proof["input_sha256"] != hashlib.sha256((self.folder / "input.pdb").read_bytes()).hexdigest():
                    raise ValueError("The tracked atom identity proof belongs to a different input structure.")
                state["measurement_verified_output_atoms"] = proof["input_to_native"]
            self.live_measurements = LiveMeasurements(self.folder, self.job, state, topology)
        except Exception as exc:
            snapshot = empty_snapshot(self.job)
            snapshot["errors"] = [f"Live measurements unavailable: {exc}"]
            atomic_json(self.folder / "measurements.json", snapshot)
            self.update(message=f"Live measurements unavailable: {exc}. Dynamics can continue.")

    def report_measurements(self, method, *args, **kwargs):
        monitor = getattr(self, "live_measurements", None)
        if monitor is not None and monitor.active:
            try:
                getattr(monitor, method)(*args, **kwargs)
            except Exception as exc:
                monitor.stop(exc)
                self.update(message=f"Live measurements stopped: {exc}. Dynamics can continue.")

    def run_command(self, arguments, stage, stdin=None, production=False, log_name=None):
        self.update(stage=stage, message="$ " + " ".join(arguments))
        self.provenance["commands"].append(arguments)
        atomic_json(self.folder / "provenance.json", self.provenance)
        output_path = self.folder / (log_name or f"command-{len(self.provenance['commands']):02d}.log")
        with output_path.open("w") as output:
            child = None
            try:
                child = subprocess.Popen(arguments, cwd=self.folder, env=os.environ.copy(), stdin=subprocess.PIPE if stdin is not None else subprocess.DEVNULL, stdout=output, stderr=subprocess.STDOUT, text=True)
                self.job["native_pid"] = child.pid
                atomic_json(self.folder / "status.json", self.job)
                if stdin is not None:
                    child.stdin.write(stdin)
                    child.stdin.close()
                last_diagnostics = 0.0
                previous = -1
                while child.poll() is None:
                    self.check_cancel()
                    completed = None
                    if production:
                        mdlog = self.folder / "production.log"
                        if mdlog.exists():
                            content = mdlog.read_text(errors="replace")[-20000:]
                            matches = re.findall(r"Step\s+Time\s*\n\s*(\d+)\s+[\d.eE+-]+", content)
                            if matches:
                                completed = min(int(matches[-1]), self.job["total_steps"])
                    if production and time.monotonic() - last_diagnostics > 3:
                        try:
                            self.export_gromacs_energies(arguments[0])
                        except (OSError, subprocess.TimeoutExpired):
                            pass
                        # This is the current command's log, not the appended
                        # production log. Read only after native startup, when
                        # checkpoint append has restored its retained prefix.
                        if "starting mdrun" in output_path.read_text(errors="replace")[-20000:]:
                            self.report_measurements("read_xtc")
                        last_diagnostics = time.monotonic()
                    if completed is not None and completed != previous:
                        self.update(completed=completed, message=f"Production step {completed:,}/{self.job['total_steps']:,}")
                        previous = completed
                    else:
                        self.update()
                    time.sleep(0.5)
            finally:
                if child is not None and child.poll() is None:
                    child.send_signal(signal.SIGTERM)
                    try:
                        child.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        child.kill()
                        child.wait()
                self.job.pop("native_pid", None)
            text = output_path.read_text(errors="replace")
            if child.returncode:
                raise RuntimeError(f"{stage} failed (exit {child.returncode}). " + text[-5000:])
            self.update(message=f"{stage} finished.")
            return text

    def run_openmm(self):
        import mdtraj as md
        import numpy as np
        import openmm as mm
        from openmm import app, unit

        settings = self.settings
        from .recovery import create_manifest, validate_manifest
        from .modified_residues import register_topology_definitions
        implicit = settings["solvent"] == "implicit"
        total = self.job["total_steps"]
        register_topology_definitions()
        if self.job.get("resume_requested"):
            manifest = validate_manifest(self.folder, "openmm")
            pdb = app.PDBFile(str(self.folder / "prepared.pdb"))
            system = mm.XmlSerializer.deserialize((self.folder / "system.xml").read_text())
            integrator = mm.XmlSerializer.deserialize((self.folder / "integrator.xml").read_text())
            simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatformByName("CPU"), {"Threads": str(config.CPU_THREADS), "DeterministicForces": "true"})
            checkpoint = manifest["checkpoint"]
            simulation.loadCheckpoint(str(self.folder / "checkpoints" / checkpoint["directory"] / "checkpoint.chk"))
            if simulation.currentStep != checkpoint["step"]:
                raise ValueError("Checkpoint step does not match the saved frame manifest.")
            self.output_solvation = manifest.get("output_solvation")
            frame_records = checkpoint["frames"]
            self.update(status="running", stage="Resuming dynamics", completed=simulation.currentStep, message=f"Native checkpoint restored at production step {simulation.currentStep:,}. Coordinates, velocities and engine random state retained.")
        else:
            self.update(stage="Preparing protein", status="running")
            from .modified_residues import register_topology_definitions
            register_topology_definitions()
            pdb = app.PDBFile(str(self.folder / "input.pdb"))
            modeller = app.Modeller(pdb.topology, pdb.positions)
            prepared_state = self.input_state.get("preparation")
            solvent_state = self.input_state.get("solvation")
            ff, files = load_prepared_forcefield(self.folder, prepared_state, solvent=settings["solvent"])
            self.provenance.update(forcefield_files=files, integrator="LangevinMiddleIntegrator", ensemble="NVT", solvent="GBn2 implicit" if implicit else "TIP3P explicit", ph=prepared_state.get("ph", 7.0) if prepared_state else 7.0, platform="CPU", input_preparation=prepared_state, input_solvation=solvent_state)
            random.seed(settings["seed"])
            np.random.seed(settings["seed"])
            self.provenance["preparation_random_seed"] = settings["seed"]
            before = modeller.topology.getNumAtoms()
            if prepared_state:
                self.record_preparation(f"Preserving exact prepared topology and hydrogens at recorded pH {prepared_state['ph']:g}; no automatic hydrogen reassignment.")
                if prepared_state.get("ligand_parameters"):
                    if not solvent_state:
                        raise ValueError("Create and inspect the explicit-water preview for this prepared complex before starting simulation.")
                    self.record_preparation("Reusing verified GAFF2 ligand templates and AM1-BCC charges; ligand molecular states and atom identities are retained.")
                if prepared_state.get("modified_residues"):
                    if not solvent_state:
                        raise ValueError("Create and inspect explicit water for the prepared modified protein before starting simulation.")
                    self.record_preparation("Reusing recorded compatible modified-residue templates, with the covalent modifications and peptide bonds retained.")
            else:
                modeller.addHydrogens(ff, pH=7.0)
                self.record_preparation(f"Added {modeller.topology.getNumAtoms() - before} hydrogens using OpenMM templates at pH 7; retained all input atoms. Protonation uses heuristic residue defaults and must be reviewed for scientific studies.")
            self.verify_measurement_identity(md.Trajectory(
                np.asarray(modeller.positions.value_in_unit(unit.nanometer))[None], md.Topology.from_openmm(modeller.topology)
            ))
            if implicit:
                modeller.topology.setPeriodicBoxVectors(None)
                system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
            else:
                before = modeller.topology.getNumAtoms()
                if solvent_state:
                    if modeller.topology.getPeriodicBoxVectors() is None:
                        raise ValueError("Prepared solvent preview has no periodic box; regenerate the preview.")
                    self.record_preparation(f"Reusing the exact prepared TIP3P solvent preview ({before:,} atoms); no water or ions added and no new box generated.")
                else:
                    modeller.addSolvent(ff, model="tip3p", padding=settings["padding_nm"] * unit.nanometer, ionicStrength=0 * unit.molar, neutralize=True)
                    self.record_preparation(f"Added explicit TIP3P solvent and neutralizing counterions with {settings['padding_nm']:g} nm padding ({modeller.topology.getNumAtoms() - before:,} added atoms); fixed volume NVT, no pressure equilibration.")
                    self.output_solvation = {"parent_dataset_id": settings["dataset_id"], "padding_nm": settings["padding_nm"], "seed": settings["seed"], "water_model": "tip3p", "method": "Legacy OpenMM simulation setup; retained to prevent duplicate solvation on continuation", "added_atoms": modeller.topology.getNumAtoms() - before}
                system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1 * unit.nanometer, constraints=app.HBonds)
            frames = math.ceil(total / settings["report_interval"]) + 1
            if modeller.topology.getNumAtoms() > config.MAX_ATOMS or frames * modeller.topology.getNumAtoms() * 12 > config.MAX_COORD_BYTES:
                raise ValueError("Prepared solvated output exceeds viewer limits. Increase the report interval or use a smaller input.")
            integrator = mm.LangevinMiddleIntegrator(settings["temperature_k"] * unit.kelvin, settings["friction_ps"] / unit.picosecond, settings["timestep_fs"] * unit.femtosecond)
            integrator.setRandomNumberSeed(settings["seed"])
            platform = mm.Platform.getPlatformByName("CPU")
            simulation = app.Simulation(modeller.topology, system, integrator, platform, {"Threads": str(config.CPU_THREADS), "DeterministicForces": "true"})
            simulation.context.setPositions(modeller.positions)
            before_min_xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer))
            self.check_stereochemistry(modeller.topology, before_min_xyz, 'Before minimization')
            # Permanent chirality guard. The force field has no explicit term
            # keeping a Calpha/Cbeta stereocenter left-handed; chirality is held
            # only by the local bonded geometry. Where a strained prepared
            # geometry (a residue beside a low-confidence rebuilt loop) distorts
            # that geometry, the force field's own local minimum can be the
            # inverted (D) center -- free minimization drives it there and the
            # first MD steps racemize it. This gentle flat-bottom restraint is
            # part of the minimization and dynamics force field: it is exactly
            # zero while a center keeps its correct sign and at least 40% of its
            # ideal signed-volume magnitude and rises only as a center nears
            # planarity, so it never biases the many healthy centers and only
            # resists racemization of the few strained ones. Because it is
            # permanent there is no release into which a held center could invert.
            guard, guarded = self._stereo_monitor.chirality_guard_force(before_min_xyz)
            if guarded:
                system.addForce(guard)
                simulation.context.reinitialize(preserveState=True)
                self.record_preparation(f"Added a permanent flat-bottom chirality restraint over {len(guarded)} standard-residue stereocenters "
                                        "(zero force while a center keeps its sign and >= 40% of its ideal signed volume; it resists only the "
                                        "approach to planarity). It is part of the minimization and dynamics force field and prevents strained "
                                        "rebuilt-loop geometry from racemizing a center; healthy centers are unaffected.")
            if settings["minimize"]:
                self.update(stage="Energy minimization", message="Minimizing up to 1,000 iterations (tolerance 10 kJ/mol/nm).")
                simulation.minimizeEnergy(tolerance=10 * unit.kilojoule_per_mole / unit.nanometer, maxIterations=1000)
                minimized = simulation.context.getState(getPositions=True)
                self.check_stereochemistry(simulation.topology, minimized.getPositions(asNumpy=True).value_in_unit(unit.nanometer), 'After minimization')
            simulation.context.setVelocitiesToTemperature(settings["temperature_k"] * unit.kelvin, settings["seed"])
            if settings["equilibration_steps"]:
                self.update(stage="Initial relaxation", message=f"Running {settings['equilibration_steps']:,} initial relaxation steps; this does not establish equilibration.")
                remaining = settings["equilibration_steps"]
                while remaining:
                    count = min(remaining, 100)
                    simulation.step(count)
                    remaining -= count
                    relaxed = simulation.context.getState(getPositions=True)
                    self.check_stereochemistry(simulation.topology, relaxed.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
                                               f"Initial relaxation step {settings['equilibration_steps'] - remaining}")
                    self.update()
            simulation.currentStep = 0
            simulation.context.setTime(0 * unit.picosecond)
            state = simulation.context.getState(getPositions=True)
            with (self.folder / "prepared.pdb").open("w") as handle:
                app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle, keepIds=True)
            (self.folder / "system.xml").write_text(mm.XmlSerializer.serialize(system))
            (self.folder / "integrator.xml").write_text(mm.XmlSerializer.serialize(integrator))
            manifest = create_manifest(self.folder, "openmm")
            manifest["output_solvation"] = getattr(self, "output_solvation", None)
            frame_records = []
        self.update(stage="Running dynamics", message=f"Production started: {total:,} steps, {settings['timestep_fs']:g} fs, seed {settings['seed']}; OpenMM {mm.__version__} CPU.")
        self.start_measurements(md.load_topology(str(self.folder / "prepared.pdb")))
        # Degrees of freedom follow OpenMM StateDataReporter: massive particles,
        # constraints between massive particles, and removed centre-of-mass motion.
        dof = sum(3 for i in range(system.getNumParticles()) if system.getParticleMass(i) > 0 * unit.dalton)
        for i in range(system.getNumConstraints()):
            a, b, _ = system.getConstraintParameters(i)
            if system.getParticleMass(a) > 0 * unit.dalton and system.getParticleMass(b) > 0 * unit.dalton:
                dof -= 1
        if any(isinstance(force, mm.CMMotionRemover) for force in system.getForces()):
            dof -= 3
        self.provenance["diagnostics"] = {"temperature_method": "2 * kinetic energy / (gas constant * constrained degrees of freedom)", "degrees_of_freedom": dof, "checkpoint_note": "Binary CPU checkpoints retain native state; compatibility is checked before loading. Bitwise identity across runs or computers is not promised."}
        atomic_json(self.folder / "provenance.json", self.provenance)
        frames_folder = self.folder / "frames"
        frames_folder.mkdir(exist_ok=True)
        checkpoints = self.folder / "checkpoints"
        checkpoints.mkdir(exist_ok=True)

        def checkpoint_now():
            from .recovery import digest
            step = simulation.currentStep
            directory = f"step-{step:09d}-{uuid.uuid4().hex[:8]}"
            temporary = checkpoints / (directory + ".tmp")
            temporary.mkdir()
            simulation.saveCheckpoint(str(temporary / "checkpoint.chk"))
            # Commit native state and its exact saved-frame prefix together.
            temporary.rename(checkpoints / directory)
            manifest["checkpoint"] = {"directory": directory, "step": step, "time_ps": step * settings["timestep_fs"] / 1000, "sha256": digest(checkpoints / directory / "checkpoint.chk"), "frames": list(frame_records)}
            atomic_json(self.folder / "recovery.json", manifest)
            for previous in checkpoints.iterdir():
                if previous.name != directory and previous.is_dir():
                    shutil.rmtree(previous)

        def read_frames():
            return [dict(np.load(frames_folder / record["name"], allow_pickle=False)) for record in frame_records]

        with (self.folder / "energies.csv").open("w") as energies:
            energies.write("step,time_ps,potential_kj_mol,kinetic_kj_mol,temperature_k\n")
            def write_energy(frame):
                energies.write(f"{int(frame['step'])},{float(frame['time']):.8f},{float(frame['potential']):.8f},{float(frame['kinetic']):.8f},{float(frame['temperature']):.8f}\n")
                energies.flush()
            # Roll back any uncommitted energy rows to the native checkpoint.
            for frame in read_frames():
                write_energy(frame)
                self.report_measurements("add_frame", frame["xyz"], float(frame["time"]), int(frame["step"]), frame.get("box"))
            while True:
                step = simulation.currentStep
                is_frame = step == 0 or step % settings["report_interval"] == 0 or step == total
                if is_frame and (not frame_records or frame_records[-1]["step"] != step):
                    state = simulation.context.getState(getPositions=True, getEnergy=True)
                    xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                    self.check_stereochemistry(simulation.topology, xyz, f'Production step {step}',
                                               None if implicit else state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
                    potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                    kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                    temperature = 2 * kinetic / (0.00831446261815324 * dof) if dof > 0 else float("nan")
                    if not np.isfinite(xyz).all() or not all(np.isfinite(v) for v in (potential, kinetic, temperature)):
                        raise RuntimeError("The engine produced non-finite coordinates or energy; simulation stopped.")
                    name = f"frame-{step:09d}.npz"
                    frame = {"step": step, "time": step * settings["timestep_fs"] / 1000, "xyz": xyz.astype(np.float32), "potential": potential, "kinetic": kinetic, "temperature": temperature}
                    if not implicit:
                        frame["box"] = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
                    temporary = frames_folder / (name + ".tmp")
                    with temporary.open("wb") as output:
                        np.savez(output, **frame)
                    temporary.replace(frames_folder / name)
                    from .recovery import digest
                    frame_records.append({"name": name, "step": step, "sha256": digest(frames_folder / name)})
                    write_energy(frame)
                    self.report_measurements("add_frame", frame["xyz"], float(frame["time"]), int(frame["step"]), frame.get("box"))
                checkpoint_now()
                self.report_measurements("flush", force=step >= total)
                self.update(completed=step, message=f"Step {step:,}/{total:,} · {step * settings['timestep_fs'] / 1000:g} ps · checkpoint saved")
                if step >= total:
                    break
                next_frame = (step // settings["report_interval"] + 1) * settings["report_interval"]
                simulation.step(min(250, next_frame - step, total - step))
        state = simulation.context.getState(getPositions=True)
        with (self.folder / "final.pdb").open("w") as handle:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle, keepIds=True)
        saved = read_frames()
        times = np.asarray([float(frame["time"]) for frame in saved])
        traj = md.Trajectory(np.array([frame["xyz"] for frame in saved], dtype=np.float32), md.load_topology(str(self.folder / "prepared.pdb")), time=times)
        if not implicit:
            traj.unitcell_vectors = np.array([frame["box"] for frame in saved], dtype=np.float32)
        # Reconstruct once from the checkpoint's committed frame prefix. Never
        # append a repeated frame after resuming a partially written trajectory.
        traj.save_dcd(str(self.folder / "trajectory.dcd"))
        np.savetxt(self.folder / "frame_times_ps.csv", times, header="time_ps", comments="")
        traj.save_xtc(str(self.folder / "trajectory.xtc"))
        shutil.copy2(checkpoints / manifest["checkpoint"]["directory"] / "checkpoint.chk", self.folder / "checkpoint.chk")
        self.provenance["timestamp_note"] = "frame_times_ps.csv and XTC carry exact production times. DCD readers may discard timing; use the sidecar. Final frame can have a shorter interval."
        self.report_measurements("finish", traj)
        return traj

    def export_gromacs_energies(self, gmx):
        """Best-effort native EDR extraction, including during a running job."""
        import csv
        import numpy as np
        energy = self.folder / "production.edr"
        if not energy.is_file() or energy.stat().st_size == 0:
            return
        target = self.folder / ".energies-live.xvg"
        result = subprocess.run([gmx, "energy", "-f", str(energy), "-o", str(target), "-xvg", "none"], input="Potential\nKinetic-En.\nTemperature\n0\n", text=True, capture_output=True, env={**os.environ, "GMX_MAXBACKUP": "-1"}, timeout=10)
        if result.returncode or not target.is_file():
            return
        rows = []
        for line in target.read_text().splitlines():
            try:
                values = [float(value) for value in line.split()]
            except ValueError:
                continue
            if len(values) == 4 and np.isfinite(values).all():
                time_ps, potential, kinetic, temperature = values
                rows.append([round(time_ps * 1000 / self.settings["timestep_fs"]), time_ps, potential, kinetic, temperature])
        if rows:
            temporary = self.folder / "energies.csv.tmp"
            with temporary.open("w") as stream:
                writer = csv.writer(stream)
                writer.writerow(["step", "time_ps", "potential_kj_mol", "kinetic_kj_mol", "temperature_k"])
                writer.writerows(rows)
            temporary.replace(self.folder / "energies.csv")
        target.unlink(missing_ok=True)

    def gromacs_production(self, gmx, cpu, resume=False):
        import mdtraj as md
        from .recovery import create_manifest, validate_manifest
        if resume:
            validate_manifest(self.folder, "gromacs")
        else:
            create_manifest(self.folder, "gromacs")
        # Resume reconstructs from the native append output after GROMACS has
        # restored/truncated it, never from the previous worker's JSON series.
        self.start_measurements(md.load_topology(str(self.folder / "system.gro")))
        arguments = [gmx, "mdrun", "-deffnm", str(self.folder / "production"), "-cpt", "0.1", *cpu]
        if resume:
            arguments += ["-cpi", "production.cpt", "-append"]
        self.provenance["checkpoint_note"] = "GROMACS native checkpoints every 6 seconds and at a graceful stop. Resume reuses the exact TPR and native append validates output checksums. No bitwise reproducibility claim."
        self.provenance["diagnostics"] = {"temperature_method": "Native GROMACS Temperature term extracted from production.edr", "energy_extraction": "gmx energy; Potential, Kinetic-En., Temperature"}
        try:
            self.run_command(arguments, "Resuming dynamics" if resume else "Running dynamics", production=True)
        finally:
            try:
                self.export_gromacs_energies(gmx)
            except (OSError, subprocess.TimeoutExpired):
                pass
            self.report_measurements("read_xtc")
            self.report_measurements("flush", force=True)
        traj = md.load(str(self.folder / "production.xtc"), top=str(self.folder / "system.gro"))
        for frame in traj:
            self.check_stereochemistry(traj.topology.to_openmm(), frame.xyz[0], f'Production time {float(frame.time[0]):g} ps', frame.unitcell_vectors[0])
        self.report_measurements("read_xtc", final=True)
        self.report_measurements("finish", traj)
        traj[0].save_pdb(str(self.folder / "prepared.pdb"))
        return traj

    def run_gromacs(self):
        import mdtraj as md
        settings = self.settings
        gmx = gromacs_executable()
        if not gmx:
            raise ValueError("GROMACS executable is unavailable.")
        self.update(stage="Preparing protein", status="running")
        cpu = ["-ntmpi", "1", "-ntomp", str(config.CPU_THREADS), "-nb", "cpu", "-pme", "cpu", "-bonded", "cpu", "-update", "cpu", "-pin", "off"]
        if self.job.get("resume_requested"):
            self.update(message="Restoring native GROMACS checkpoint; reusing the original TPR, topology, velocities and random state. Native append will verify all prior output checksums.")
            return self.gromacs_production(gmx, cpu, resume=True)
        version = self.run_command([gmx, "--version"], "Engine probe")
        self.provenance.update(gromacs_version=version, forcefield="amber99sb-ildn", water_model="tip3p", integrator="GROMACS stochastic dynamics (sd)", ensemble="NVT", ph="pdb2gmx default protonation (review required)")
        self.record_preparation("pdb2gmx rebuilds hydrogen coordinates with Amber99SB-ILDN templates (-ignh); all input heavy atoms must be retained. Default termini and protonation require review. No ligand parameterization.")
        self.run_command([gmx, "pdb2gmx", "-f", "input.pdb", "-o", "processed.gro", "-p", "topol.top", "-ff", "amber99sb-ildn", "-water", "tip3p", "-ignh"], "Building topology")
        initial = md.load(str(self.folder / "input.pdb"))
        processed = md.load(str(self.folder / "processed.gro"))
        # GRO rounds Cartesian coordinates to 0.001 nm (PDB to 0.0001 nm).
        self.verify_measurement_identity(processed, coordinate_tolerance=6e-4)
        # GROMACS may reorder atoms and add terminal atoms. Detect removal of heavy atoms by per-residue element counts.
        from collections import Counter
        def heavy_counts(top):
            return Counter((a.residue.index, a.element.symbol) for a in top.atoms if a.element and a.element.symbol != "H")
        missing = heavy_counts(initial.topology) - heavy_counts(processed.topology)
        if missing:
            raise ValueError(f"GROMACS preparation removed input heavy atoms ({dict(missing)}); refusing to continue.")
        self.run_command([gmx, "editconf", "-f", "processed.gro", "-o", "boxed.gro", "-c", "-d", str(settings["padding_nm"]), "-bt", "cubic"], "Creating periodic box")
        self.run_command([gmx, "solvate", "-cp", "boxed.gro", "-cs", "spc216.gro", "-o", "solvated.gro", "-p", "topol.top"], "Adding TIP3P water")
        # Restrict ion replacement to newly added waters; preserve input crystallographic waters.
        boxed = md.load(str(self.folder / "boxed.gro"))
        solvated = md.load(str(self.folder / "solvated.gro"))
        import numpy as np
        if not np.allclose(solvated.xyz[0, :boxed.n_atoms], boxed.xyz[0], atol=1e-5):
            raise ValueError("Solvation changed the original coordinate prefix; cannot safely identify newly added solvent for neutralization.")
        added_water_atoms = list(range(boxed.n_atoms + 1, solvated.n_atoms + 1))
        if not added_water_atoms or len(added_water_atoms) % 3:
            raise ValueError("Could not identify a complete group of newly added three-site water molecules.")
        group_lines = ["[ SOL ]"] + [" ".join(map(str, added_water_atoms[i:i+15])) for i in range(0, len(added_water_atoms), 15)]
        (self.folder / "added-solvent.ndx").write_text("\n".join(group_lines) + "\n")
        # The coordinate template spc216 is used for coordinates only; topology selects TIP3P.
        minimization = "integrator = steep\nnsteps = 2000\nemtol = 1000\nemstep = 0.01\ncutoff-scheme = Verlet\nnstlist = 10\nrlist = 1.0\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n"
        (self.folder / "ions.mdp").write_text(minimization.replace("coulombtype = PME", "coulombtype = Cut-off"))
        (self.folder / "minimize.mdp").write_text(minimization)
        self.run_command([gmx, "grompp", "-f", "ions.mdp", "-c", "solvated.gro", "-p", "topol.top", "-o", "ions.tpr"], "Preparing neutralization")
        self.run_command([gmx, "genion", "-s", "ions.tpr", "-o", "system.gro", "-p", "topol.top", "-pname", "NA", "-nname", "CL", "-neutral", "-seed", str(settings["seed"]), "-n", "added-solvent.ndx"], "Neutralizing system", stdin="SOL\n")
        self.record_preparation(f"Added TIP3P solvent and neutralizing ions, replacing only newly added waters; retained input crystal waters. {settings['padding_nm']:g} nm cubic padding. NVT at fixed volume; no pressure equilibration or added salt.")
        prepared = md.load(str(self.folder / "system.gro"))
        if self.settings.get("measurements") and (
            prepared.n_atoms < boxed.n_atoms or not np.allclose(prepared.xyz[0, :boxed.n_atoms], boxed.xyz[0], atol=1e-5)
            or [(a.residue.index, a.name, a.element) for a in list(prepared.topology.atoms)[:boxed.n_atoms]] != [(a.residue.index, a.name, a.element) for a in boxed.topology.atoms]
        ):
            proof_path = self.folder / "measurement-identity.json"
            proof = json.loads(proof_path.read_text())
            proof.update(verified=False, error="Native solvation or ion addition changed the retained atom prefix; tracked atom identities could not be verified.")
            atomic_json(proof_path, proof)
        if prepared.n_atoms > config.MAX_ATOMS or (math.ceil(self.job["total_steps"] / settings["report_interval"]) + 1) * prepared.n_atoms * 12 > config.MAX_COORD_BYTES:
            raise ValueError("Prepared solvated output exceeds viewer limits. Increase report interval or use a smaller input.")
        cpu = ["-ntmpi", "1", "-ntomp", str(config.CPU_THREADS), "-nb", "cpu", "-pme", "cpu", "-bonded", "cpu", "-update", "cpu", "-pin", "off"]
        coordinates = "system.gro"
        self.check_stereochemistry(prepared.topology.to_openmm(), prepared.xyz[0], 'Before minimization', prepared.unitcell_vectors[0])
        if settings["minimize"]:
            self.run_command([gmx, "grompp", "-f", "minimize.mdp", "-c", coordinates, "-p", "topol.top", "-o", "minimize.tpr"], "Preparing minimization")
            self.run_command([gmx, "mdrun", "-deffnm", str(self.folder / "minimize"), *cpu], "Energy minimization")
            coordinates = "minimize.gro"
            minimized = md.load(str(self.folder / coordinates))
            self.check_stereochemistry(minimized.topology.to_openmm(), minimized.xyz[0], 'After minimization', minimized.unitcell_vectors[0])
        def mdp(steps, continuation):
            return f"""integrator = sd
nsteps = {steps}
dt = {settings['timestep_fs'] / 1000:.9f}
nstxout-compressed = {settings['report_interval']}
nstenergy = {settings['report_interval']}
nstlog = {min(settings['report_interval'], 100)}
cutoff-scheme = Verlet
nstlist = 10
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
constraints = h-bonds
constraint-algorithm = lincs
pbc = xyz
tc-grps = System
tau-t = {1 / settings['friction_ps']:.8f}
ref-t = {settings['temperature_k']}
ld-seed = {settings['seed']}
continuation = {'yes' if continuation else 'no'}
gen-vel = {'no' if continuation else 'yes'}
gen-temp = {settings['temperature_k']}
gen-seed = {settings['seed']}
pcoupl = no
"""
        if settings["equilibration_steps"]:
            (self.folder / "relax.mdp").write_text(mdp(settings["equilibration_steps"], False))
            self.run_command([gmx, "grompp", "-f", "relax.mdp", "-c", coordinates, "-p", "topol.top", "-o", "relax.tpr"], "Preparing initial relaxation")
            self.run_command([gmx, "mdrun", "-deffnm", str(self.folder / "relax"), *cpu], "Initial relaxation")
            coordinates = "relax.gro"
            relaxed = md.load(str(self.folder / coordinates))
            self.check_stereochemistry(relaxed.topology.to_openmm(), relaxed.xyz[0], 'After initial relaxation', relaxed.unitcell_vectors[0])
        (self.folder / "production.mdp").write_text(mdp(self.job["total_steps"], bool(settings["equilibration_steps"])))
        self.run_command([gmx, "grompp", "-f", "production.mdp", "-c", coordinates, "-p", "topol.top", "-o", "production.tpr"], "Preparing production")
        return self.gromacs_production(gmx, cpu)

    def run(self):
        try:
            trajectory = self.run_openmm() if self.settings["engine"] == "openmm" else self.run_gromacs()
            self.update(stage="Preparing trajectory for viewer", completed=self.job["total_steps"])
            self.provenance["outputs"] = {str(p.relative_to(self.folder)): hashlib.sha256(p.read_bytes()).hexdigest() for p in self.folder.rglob("*") if p.is_file() and p.name not in {"status.json", "worker.log", "provenance.json"} and p.suffix not in {".tmp", ".zip"}}
            atomic_json(self.folder / "provenance.json", self.provenance)
            warnings = ["Short exploratory simulation. This run does not establish equilibration, convergence, biological function, or binding stability.", "Protein force-field templates, terminal states and protonation require scientific review before a production study."]
            # Reserve the result identity before materializing it. A crash during
            # viewer conversion can then safely rebuild the same result instead
            # of leaving multiple apparent trajectories in the user's library.
            if not self.job.get("output_dataset_id"):
                self.update(output_dataset_id=uuid.uuid4().hex[:16])
            dataset = save_dataset(trajectory, self.settings["name"], f"{self.settings['engine']} simulation", f"{trajectory.time[-1]:g} ps of real {self.settings['engine']} dynamics at {self.settings['temperature_k']:g} K. Fixed-volume exploratory NVT trajectory; seed {self.settings['seed']}.", warnings=warnings, dataset_id=self.job["output_dataset_id"], provenance=self.provenance)
            # A simulation result is itself a valid starting dataset. Carry exact
            # preparation/solvent state forward so continuation cannot reset pH or
            # regenerate a second solvent box.
            prepared_state = copy.deepcopy(self.input_state.get("preparation"))
            solvent_state = copy.deepcopy(self.input_state.get("solvation") or getattr(self, "output_solvation", None))
            dataset["parent_dataset_id"] = self.settings["dataset_id"]
            if prepared_state:
                dataset["preparation"] = prepared_state
                copy_ligand_parameters(self.folder, dataset_dir(dataset["id"]), prepared_state)
            if solvent_state:
                dataset["solvation"] = solvent_state
            if prepared_state or solvent_state:
                import shutil
                exact_topology = self.folder / "prepared.pdb"
                if not exact_topology.is_file():
                    raise ValueError("Simulation output lacks the exact prepared first-frame topology required to preserve its state.")
                shutil.copy2(exact_topology, dataset_dir(dataset["id"]) / "prepared.pdb")
            atomic_json(dataset_dir(dataset["id"]) / "metadata.json", dataset)
            self.update(status="completed", stage="Complete", dataset_id=dataset["id"], message=f"Complete: {trajectory.n_frames} frames, {trajectory.n_atoms:,} atoms. Output files and reproducibility manifest saved.")
        except (Cancelled, KeyboardInterrupt):
            self.job.update(status="cancelled", stage="Cancelled", error="Cancelled by user. Partial outputs retained.", elapsed_seconds=round(self.previous_elapsed + time.monotonic() - self.started, 2))
            atomic_json(self.folder / "status.json", self.job)
            atomic_json(self.folder / "provenance.json", self.provenance)
        except Exception as exc:
            traceback.print_exc()
            self.job.update(status="failed", stage="Failed", error=str(exc)[-6000:], elapsed_seconds=round(self.previous_elapsed + time.monotonic() - self.started, 2))
            self.job["logs"] = (self.job["logs"] + ["ERROR: " + str(exc)[-1500:]])[-150:]
            atomic_json(self.folder / "status.json", self.job)
            self.provenance["failure"] = str(exc)
            atomic_json(self.folder / "provenance.json", self.provenance)


def main():
    def terminate(signum, frame):
        # Do not raise asynchronously between native Popen and child ownership.
        # The worker's regular checkpoints handle graceful cancellation.
        global _termination_requested
        _termination_requested = True
    signal.signal(signal.SIGTERM, terminate)
    Worker(sys.argv[1]).run()


if __name__ == "__main__":
    main()
