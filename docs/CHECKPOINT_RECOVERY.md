# Simulation checkpoints and recovery

DynaMol saves native production checkpoints automatically. Failed, interrupted and stopped simulation cards offer **Resume** when a complete checkpoint exists, and **Files** always remains available after a job stops. Resume continues the original job and configuration; it does not run preparation, hydrogen assignment, solvation, minimization or initial relaxation again.

## OpenMM

A binary CPU checkpoint is committed at the start of production, at every saved frame and at most every 250 production steps. Its manifest identifies the exact saved trajectory prefix. Frame coordinates and energies are persisted before the checkpoint becomes current, so a killed process cannot cause an incomplete frame to be treated as saved. Resume loads the native checkpoint, rolls diagnostics back to that committed prefix, then reconstructs DCD/XTC output with each frame represented once. XTC and `frame_times_ps.csv` contain physical production times, including a shorter final reporting interval. A checkpoint at the final step can also recover interrupted trajectory processing.

The checkpoint contains OpenMM's native simulation state, including velocities and random generator state. DynaMol checks CPU model, operating system, thread count, OpenMM version/native library fingerprints and worker code, plus hashes of the original input, configuration, serialized System/integrator and complete ligand/modified-residue parameter directories. Native checkpoint loading supplies an additional compatibility check. This is a conservative local recovery feature; it does not promise bitwise identity or portability to another computer or upgraded runtime.

## GROMACS

`mdrun` saves its native checkpoint every six seconds and on graceful termination. Resume uses the exact `production.tpr` and `production.cpt` with native `-append`; GROMACS verifies checksums of its preceding outputs. DynaMol additionally checks input/configuration/topology/parameter hashes, runtime and native library fingerprints. Modified or missing append outputs cause the native engine to reject continuation rather than concatenate inconsistent data.

Stopping a job signals its supervisor, which lets the native child finish its current operation and save before exit. If the supervisor was killed, DynaMol identifies the recorded native child for that job, asks it to stop, and waits for it to exit before permitting Resume. GROMACS handles coordinates, velocities, random state and existing output positions through its checkpoint; bitwise reproducibility across process restarts is not claimed.

## Boundaries

Preparation and initial relaxation must complete before the first production checkpoint exists. A job stopped sooner retains inputs/logs but must be started again. Structural preparation/solvation jobs do not offer MD Resume. A binary checkpoint from an old runtime may remain available in Files even when compatibility checks prevent resuming it in an upgraded app.

Closing a browser does not stop simulation workers. The native **Quit** menu checks active jobs and offers **Quit viewer, keep running**, **Stop jobs and quit**, or **Keep DynaMol open**. Stopping waits for the jobs to become terminal before shutting down the viewer service. A computer reboot or forced engine termination can lose work since the last complete checkpoint.

## Diagnostics

OpenMM writes potential energy, kinetic energy and temperature at each stored frame. Temperature is calculated from kinetic energy using massive-particle degrees of freedom, subtracting constraints and removed centre-of-mass motion. GROMACS temperature and energies come directly from its native EDR terms, extracted during production when readable and once again after it stops. Early stages or an EDR buffer with no readable samples can legitimately have no chart points yet. These diagnostics are numerical health checks for exploratory NVT runs; they do not establish convergence or production readiness.

Native check script: `scripts/check_recovery.py`. Recorded evidence: `docs/audit/release-readiness/recovery-native.json`. The script deliberately interrupts and stops short isolated runs, rejects changed systems/runtimes, then checks retained frame prefixes, increasing unique times, native append, restart provenance and finite coordinates/energies/temperature.
