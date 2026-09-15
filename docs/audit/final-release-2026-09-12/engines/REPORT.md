# Engine, stability and recovery audit — 12 September 2026

**Pass for the bounded software audit.** Both current-source native engines completed fresh 10 ps ubiquitin simulations in explicit TIP3P water, including two interruption/restart cycles per engine. Their saved results open through the application API with intact live measurement definitions, and downloaded native files match the originals. This is evidence of local execution and artifact continuity, not of equilibrated or converged molecular dynamics.

## Fresh current-source runs

Both runs used seed 20260912, 300 K Langevin/stochastic NVT dynamics, a 2 fs timestep, 1 nm solvent padding, 500 initial relaxation steps, 5,000 production steps, and a frame every 100 steps. Only two CPU threads were allocated to the active engine. Initial preparation, minimization, restart settings, native runtime identities, input hashes and worker source hashes are retained in the isolated job folders.

| Engine | Water-box atoms | Saved frames | Production | Wall time including setup/restarts | Artifact checks |
|---|---:|---:|---:|---:|---:|
| OpenMM 8.6.0 | 17,242 | 51 | 10 ps | 161.90 s | 21 passed |
| GROMACS 2025.4-conda_forge | 26,842 | 51 | 10 ps | 62.12 s | 16 passed |

The two engines use their documented distinct force fields and independently built boxes. These runs do not compare force-field energies or establish equivalence between engines.

- Every stored coordinate, energy, temperature and time is finite. Frame times are unique and strictly increasing from 0 to 10 ps; periodic boxes have finite positive fixed volume.
- Native XTC and app physical trajectories agree within their documented format precision. OpenMM's DCD and time sidecar also agree. Display payload dimensions and recorded file hashes match.
- All 84 monitored protein stereocenters pass in each of 51 frames. All 608 protein heavy-atom bonds pass the gross length screen; maximum observed lengths are 1.886 Å (OpenMM) and 1.885 Å (GROMACS).
- The 300 K target is approached during these short runs: recorded temperatures span 243.4–302.3 K for OpenMM and 263.2–303.2 K for GROMACS. These are intentionally short relaxation/software tests; the broad numerical screen is not an equilibration test.
- The native logs contain no searched LINCS/SHAKE, fatal-error or nonfinite-coordinate/energy alarms. GROMACS's preparatory steepest-descent minimization reached its configured force threshold in 437 steps. The log scan is a gross alarm search, not exhaustive scientific review.

## Recovery and opening results

Each worker was deliberately killed after an established native checkpoint. GROMACS's orphaned native child was stopped before recovery. Both engines rejected deliberately changed system/TPR bytes and a modified runtime-thread identity. Restoring the exact originals enabled native checkpoint resume; a subsequent cooperative Stop and second Resume both worked.

OpenMM resumed after the first checkpoint at step 500 and the cooperative stop at step 1,100. GROMACS was first interrupted after reported step 800 and cooperatively stopped after reported step 1,300. The exact native checkpoint is authoritative for GROMACS. The retained saved prefix remains unchanged, including ten retained GROMACS native frames. Both runs finish with 51 frames and two restart provenance entries.

For each engine, **12 HTTP handoff/export checks pass**: completed job and dataset opening, topology and coordinate delivery, retained distance/angle/dihedral definitions, agreement between live and posthoc values, ZIP integrity, and byte-for-byte native/configuration/provenance/energy contents in downloaded output. This exercises the same application routes used by the UI through an isolated loopback test client. Actual browser mode switching is covered by the separate UI audit.

## Additional historical evidence

The current saved-artifact validator was run read-only over all five prior 100 ps systems: 1XQZ kinase (GROMACS), 3PTB protease with benzamidine/ion chemistry (OpenMM), 1SUG phosphatase (GROMACS), 1K4C ion channel (OpenMM), and 4N6H GPCR with its ligand (OpenMM). **All 95 checks pass.** No new dynamics were run for those systems in this audit. Their original preparation choices and limitations remain recorded in the previous pressure-test report. In particular, water-only membrane-protein cases are software exercises, not physical membrane simulations; the 145,913-atom GPCR also requires the documented larger resource profile.

## Fix and regression checks

The audit found stale simulation provenance: `backend/worker.py` exported the application as version 0.1.0. It now reads the project version from the `pyproject.toml` shipped in both source and application bundles. Fresh files correctly identify DynaMol 0.1.1; resumed runs retain their original application provenance.

The initial focused recovery/resource/validator/stereochemistry test set passed **56 tests**. After the small provenance correction, the new provenance test plus repeated recovery tests passed **34 tests**; these counts overlap and must not be added as unique coverage. The first API harness attempt used TestClient's default non-loopback hostname, which the app correctly rejected. The harness was corrected to `127.0.0.1`; no host-protection code changed.

## Runtime and practical limits

The source runtime probe and native runs used `PATH=/usr/bin:/bin:/usr/sbin:/sbin`. No `gmx` was available from PATH; DynaMol used its private `.gromacs/bin/gmx`, and OpenMM came from `.venv`. Thus these fresh runs did not rely on a separately installed system MD engine. They are not a new DMG installation or separate clean-machine test. Existing relocated packaged-runtime evidence remains separate and historical.

The fresh runs cover a small standard soluble protein and three live geometry types. Variety of ligand, cofactor, modified-residue and ion preparation is covered in the separate chemistry audit and historical trajectory revalidation. These results do not establish production sampling quality, native loop accuracy, binding-pose correctness, free energies, membrane stability, GPU operation, Windows/Linux support or cross-machine checkpoint portability.

## Evidence and reproduction

- `native-recovery.json`: exact fresh settings, identities, restart evidence and output hashes.
- `fresh-openmm-validation.json`, `fresh-gromacs-validation.json`: numerical and native-artifact checks.
- `fresh-handoff-validation.json`: actual application-route and ZIP checks.
- `native-log-screen.json`, `runtime-probe.json`, `focused-tests.json`: supporting checks.
- `historical-revalidation.json` and `historical-*.json`: clearly labeled prior-output revalidation.
- `run_native_audit.py`: isolated fresh-run driver; requires both private engines and starts new jobs.
- `validate_fresh.py`: read-only validation and API export checks for the saved fresh jobs.

Native data remain under `data/integration-checks/final-release-2026-09-12/engines/`, separate from user projects. Nothing was published or repackaged by this subtask.
