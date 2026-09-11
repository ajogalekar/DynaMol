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
import subprocess
import sys
import time
import traceback
from pathlib import Path

from . import config
from .jobs import gromacs_executable
from .prepared_system import copy_ligand_parameters, load_prepared_forcefield
from .storage import atomic_json, save_dataset, safe_id, dataset_dir


class Cancelled(Exception):
    pass


class Worker:
    def __init__(self, job_id: str):
        self.folder = config.JOBS_DIR / safe_id(job_id)
        for _ in range(50):
            self.job = json.loads((self.folder / "status.json").read_text())
            if self.job.get("worker_pid"):
                break
            time.sleep(0.02)
        self.settings = self.job["config"]
        state_path = self.folder / "input-state.json"
        self.input_state = json.loads(state_path.read_text()) if state_path.exists() else {}
        self.started = time.monotonic()
        self.provenance = {"application": "DynaMol 0.1.0", "purpose": "Short exploratory molecular dynamics/software demonstration, not converged scientific validation.", "config": self.settings, "input_sha256": hashlib.sha256((self.folder / "input.pdb").read_bytes()).hexdigest(), "cpu_threads": config.CPU_THREADS, "versions": {name: importlib.metadata.version(name) for name in ("openmm", "mdtraj", "numpy")}, "commands": [], "preparation": [], "worker_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "source_dataset_id": self.settings["dataset_id"]}

    def check_cancel(self):
        if (self.folder / "cancel.request").exists():
            raise Cancelled("Cancelled by user.")

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
        self.job["elapsed_seconds"] = round(time.monotonic() - self.started, 2)
        atomic_json(self.folder / "status.json", self.job)

    def record_preparation(self, message):
        self.provenance["preparation"].append(message)
        self.update(message=message)

    def run_command(self, arguments, stage, stdin=None, production=False, log_name=None):
        self.update(stage=stage, message="$ " + " ".join(arguments))
        self.provenance["commands"].append(arguments)
        output_path = self.folder / (log_name or f"command-{len(self.provenance['commands']):02d}.log")
        with output_path.open("w") as output:
            child = subprocess.Popen(arguments, cwd=self.folder, env=os.environ.copy(), stdin=subprocess.PIPE if stdin is not None else subprocess.DEVNULL, stdout=output, stderr=subprocess.STDOUT, text=True)
            if stdin is not None:
                child.stdin.write(stdin)
                child.stdin.close()
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
                if completed is not None and completed != previous:
                    self.update(completed=completed, message=f"Production step {completed:,}/{self.job['total_steps']:,}")
                    previous = completed
                else:
                    self.update()
                time.sleep(0.5)
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
        self.update(stage="Preparing protein", status="running")
        pdb = app.PDBFile(str(self.folder / "input.pdb"))
        modeller = app.Modeller(pdb.topology, pdb.positions)
        implicit = settings["solvent"] == "implicit"
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
        else:
            modeller.addHydrogens(ff, pH=7.0)
            self.record_preparation(f"Added {modeller.topology.getNumAtoms() - before} hydrogens using OpenMM templates at pH 7; retained all input atoms. Protonation uses heuristic residue defaults and must be reviewed for scientific studies.")
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
        total = self.job["total_steps"]
        frames = math.ceil(total / settings["report_interval"]) + 1
        if modeller.topology.getNumAtoms() > config.MAX_ATOMS or frames * modeller.topology.getNumAtoms() * 12 > config.MAX_COORD_BYTES:
            raise ValueError("Prepared solvated output exceeds viewer limits. Increase the report interval or use a smaller input.")
        integrator = mm.LangevinMiddleIntegrator(settings["temperature_k"] * unit.kelvin, settings["friction_ps"] / unit.picosecond, settings["timestep_fs"] * unit.femtosecond)
        integrator.setRandomNumberSeed(settings["seed"])
        platform = mm.Platform.getPlatformByName("CPU")
        simulation = app.Simulation(modeller.topology, system, integrator, platform, {"Threads": str(config.CPU_THREADS), "DeterministicForces": "true"})
        simulation.context.setPositions(modeller.positions)
        if settings["minimize"]:
            self.update(stage="Energy minimization", message="Minimizing up to 1,000 iterations (tolerance 10 kJ/mol/nm).")
            simulation.minimizeEnergy(tolerance=10 * unit.kilojoule_per_mole / unit.nanometer, maxIterations=1000)
        simulation.context.setVelocitiesToTemperature(settings["temperature_k"] * unit.kelvin, settings["seed"])
        if settings["equilibration_steps"]:
            self.update(stage="Initial relaxation", message=f"Running {settings['equilibration_steps']:,} initial relaxation steps; this does not establish equilibration.")
            remaining = settings["equilibration_steps"]
            while remaining:
                count = min(remaining, 100)
                simulation.step(count)
                remaining -= count
                self.update()
        simulation.currentStep = 0
        simulation.context.setTime(0 * unit.picosecond)
        state = simulation.context.getState(getPositions=True)
        with (self.folder / "prepared.pdb").open("w") as handle:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle, keepIds=True)
        (self.folder / "system.xml").write_text(mm.XmlSerializer.serialize(system))
        (self.folder / "integrator.xml").write_text(mm.XmlSerializer.serialize(integrator))
        self.update(stage="Running dynamics", message=f"Production started: {total:,} steps, {settings['timestep_fs']:g} fs, seed {settings['seed']}; OpenMM {mm.__version__} CPU.")
        coordinates = []
        times = []
        cell_vectors = []
        with (self.folder / "trajectory.dcd").open("wb") as output, (self.folder / "energies.csv").open("w") as energies:
            dcd = app.DCDFile(output, simulation.topology, settings["timestep_fs"] * unit.femtosecond, firstStep=0, interval=settings["report_interval"])
            energies.write("step,time_ps,potential_kj_mol,kinetic_kj_mol\n")
            while True:
                step = simulation.currentStep
                state = simulation.context.getState(getPositions=True, getEnergy=True)
                xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                if not np.isfinite(xyz).all() or not np.isfinite(potential) or not np.isfinite(kinetic):
                    raise RuntimeError("The engine produced non-finite coordinates or energy; simulation stopped.")
                box = state.getPeriodicBoxVectors()
                dcd.writeModel(state.getPositions(), periodicBoxVectors=None if implicit else box)
                coordinates.append(xyz)
                times.append(step * settings["timestep_fs"] / 1000)
                if not implicit:
                    cell_vectors.append(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
                energies.write(f"{step},{times[-1]:.8f},{potential:.8f},{kinetic:.8f}\n")
                energies.flush()
                self.update(completed=step, message=f"Step {step:,}/{total:,} · {times[-1]:g} ps · potential {potential:.2f} kJ/mol")
                if step >= total:
                    break
                simulation.step(min(settings["report_interval"], total - step))
        simulation.saveCheckpoint(str(self.folder / "checkpoint.chk"))
        with (self.folder / "final.pdb").open("w") as handle:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle, keepIds=True)
        traj = md.Trajectory(np.array(coordinates, dtype=np.float32), md.load_topology(str(self.folder / "prepared.pdb")), time=np.asarray(times))
        if cell_vectors:
            traj.unitcell_vectors = np.array(cell_vectors, dtype=np.float32)
        np.savetxt(self.folder / "frame_times_ps.csv", times, header="time_ps", comments="")
        # XTC stores explicit physical times, including a shorter final interval.
        traj.save_xtc(str(self.folder / "trajectory.xtc"))
        self.provenance["timestamp_note"] = "frame_times_ps.csv and XTC carry exact production times. DCD readers may discard timing; use the sidecar. Final frame can have a shorter interval."
        return traj

    def run_gromacs(self):
        import mdtraj as md
        settings = self.settings
        gmx = gromacs_executable()
        if not gmx:
            raise ValueError("GROMACS executable is unavailable.")
        self.update(stage="Preparing protein", status="running")
        version = self.run_command([gmx, "--version"], "Engine probe")
        self.provenance.update(gromacs_version=version, forcefield="amber99sb-ildn", water_model="tip3p", integrator="GROMACS stochastic dynamics (sd)", ensemble="NVT", ph="pdb2gmx default protonation (review required)")
        self.record_preparation("pdb2gmx rebuilds hydrogen coordinates with Amber99SB-ILDN templates (-ignh); all input heavy atoms must be retained. Default termini and protonation require review. No ligand parameterization.")
        self.run_command([gmx, "pdb2gmx", "-f", "input.pdb", "-o", "processed.gro", "-p", "topol.top", "-ff", "amber99sb-ildn", "-water", "tip3p", "-ignh"], "Building topology")
        initial = md.load(str(self.folder / "input.pdb"))
        processed = md.load(str(self.folder / "processed.gro"))
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
        if prepared.n_atoms > config.MAX_ATOMS or (math.ceil(self.job["total_steps"] / settings["report_interval"]) + 1) * prepared.n_atoms * 12 > config.MAX_COORD_BYTES:
            raise ValueError("Prepared solvated output exceeds viewer limits. Increase report interval or use a smaller input.")
        cpu = ["-ntmpi", "1", "-ntomp", str(config.CPU_THREADS), "-nb", "cpu", "-pme", "cpu", "-bonded", "cpu", "-update", "cpu", "-pin", "off"]
        coordinates = "system.gro"
        if settings["minimize"]:
            self.run_command([gmx, "grompp", "-f", "minimize.mdp", "-c", coordinates, "-p", "topol.top", "-o", "minimize.tpr"], "Preparing minimization")
            self.run_command([gmx, "mdrun", "-deffnm", "minimize", *cpu], "Energy minimization")
            coordinates = "minimize.gro"
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
            self.run_command([gmx, "mdrun", "-deffnm", "relax", *cpu], "Initial relaxation")
            coordinates = "relax.gro"
        (self.folder / "production.mdp").write_text(mdp(self.job["total_steps"], bool(settings["equilibration_steps"])))
        self.run_command([gmx, "grompp", "-f", "production.mdp", "-c", coordinates, "-p", "topol.top", "-o", "production.tpr"], "Preparing production")
        self.run_command([gmx, "mdrun", "-deffnm", "production", *cpu], "Running dynamics", production=True)
        traj = md.load(str(self.folder / "production.xtc"), top=str(self.folder / "system.gro"))
        traj[0].save_pdb(str(self.folder / "prepared.pdb"))
        return traj

    def run(self):
        try:
            trajectory = self.run_openmm() if self.settings["engine"] == "openmm" else self.run_gromacs()
            self.update(stage="Preparing trajectory for viewer", completed=self.job["total_steps"])
            self.provenance["outputs"] = {str(p.relative_to(self.folder)): hashlib.sha256(p.read_bytes()).hexdigest() for p in self.folder.rglob("*") if p.is_file() and p.name not in {"status.json", "worker.log", "provenance.json"} and p.suffix not in {".tmp", ".zip"}}
            atomic_json(self.folder / "provenance.json", self.provenance)
            warnings = ["Short exploratory simulation. This run does not establish equilibration, convergence, biological function, or binding stability.", "Protein force-field templates, terminal states and protonation require scientific review before a production study."]
            dataset = save_dataset(trajectory, self.settings["name"], f"{self.settings['engine']} simulation", f"{trajectory.time[-1]:g} ps of real {self.settings['engine']} dynamics at {self.settings['temperature_k']:g} K. Fixed-volume exploratory NVT trajectory; seed {self.settings['seed']}.", warnings=warnings, provenance=self.provenance)
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
            self.job.update(status="cancelled", stage="Cancelled", error="Cancelled by user. Partial outputs retained.", elapsed_seconds=round(time.monotonic() - self.started, 2))
            atomic_json(self.folder / "status.json", self.job)
        except Exception as exc:
            traceback.print_exc()
            self.job.update(status="failed", stage="Failed", error=str(exc)[-6000:], elapsed_seconds=round(time.monotonic() - self.started, 2))
            self.job["logs"] = (self.job["logs"] + ["ERROR: " + str(exc)[-1500:]])[-150:]
            atomic_json(self.folder / "status.json", self.job)
            self.provenance["failure"] = str(exc)
            atomic_json(self.folder / "provenance.json", self.provenance)


def main():
    def terminate(signum, frame):
        raise Cancelled("Cancelled by user.")
    signal.signal(signal.SIGTERM, terminate)
    Worker(sys.argv[1]).run()


if __name__ == "__main__":
    main()
