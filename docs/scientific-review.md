# DynaMol scientific and numerical review

This is a language-model review of implementation, executable checks, and cited documentation. It is not a review by an independent MD expert or evidence that a simulation has converged. The bundled short protein trajectory is an integration demonstration. It cannot establish equilibrium behavior, binding, stability, a mechanism, or agreement between engines.

## Review status

Preflight permits exploratory software development and a short real OpenMM demonstration, with unresolved version and seed records to be resolved from execution artifacts. The automated PatAgent checklist is inspired by public writings; it is not Pat Walters, is not endorsed by him, and does not speak for him. Its self-test detected its deliberately flawed fixtures. An automated checklist does not validate physical correctness.

Storage, analysis, renderer, and engine code have been inspected. Eight independent numerical fixture groups pass in `audit/numerical-checks.json`; six chart fixture groups pass in `audit/plot-checks.json`; 19 backend tests pass (the initial 14-test run is retained in `audit/backend-tests.txt`). The executable independent checks are `audit/check_numerics.py` and `audit/check_plots.cjs`. These counts describe software checks, not independent scientific validations.

The postflight checklist has four passes and three failures. It **blocks scientific-result promotion** because independent replicas/convergence and uncertainty estimates were not established. It also flags incomplete preparation RNG capture for the initial OpenMM demonstration. This does not negate the evidence that the app executes its tested exploratory workflows; it limits what can be claimed from those trajectories. The preflight and postflight manifests and machine-readable reports are `pipeline-manifest.preflight.json`, `pipeline-manifest.postflight.json`, `audit/pat-preflight.json`, and `audit/pat-postflight.json`.

## Numerical contracts to verify

| Surface | Required interpretation | Review criterion |
| --- | --- | --- |
| Coordinates | API and rendering use Å; MDTraj uses nm | Explicit conversion by 10 once at the boundary |
| Time | API and plots use ps; configured time step uses fs | Explicit conversion by 0.001; distinguish frame number from physical time |
| Picking | Zero-based canonical topology index | The PDB shown by the renderer and every coordinate frame must retain identical atom order |
| Distance | Å from saved frames | Minimum-image distances when valid unit cells are present |
| Angle | Degrees with the middle selected atom at the vertex | Convert MDTraj radians to degrees; reject coincident atom geometry |
| Dihedral | Signed degrees in a documented convention | Reject collinear planes; avoid ordinary arithmetic means across the ±180° branch cut |
| Display interpolation | Visual transition between saved conformations | Never use interpolated frames for measurements; handle periodic jumps explicitly |
| Hydrogen bond | User-picked D–H–A geometry | Require an explicit H and bonded D–H; report geometric occupancy and criteria |

MDTraj uses nm, degrees for cell angles, and ps internally, while its angle and dihedral functions return radians. Its distance, angle, and dihedral routines apply minimum-image conventions when requested and box information is available. These unit distinctions must be handled at the app boundary. Sources: [MDTraj trajectory units](https://mdtraj.org/1.9.3/api/generated/mdtraj.Trajectory.html), [distance API](https://www.mdtraj.org/1.9.7/api/generated/mdtraj.compute_distances.html), [angle API](https://mdtraj.readthedocs.io/en/stable/api/generated/mdtraj.compute_angles.html), [dihedral API](https://mdtraj.readthedocs.io/en/stable/api/generated/mdtraj.compute_dihedrals.html).

The implemented hydrogen-bond test is D–A ≤ 3.5 Å and D–H–A ≥ 150°. It is a chosen geometric definition, not evidence of hydrogen-bond energy, and not an exact reproduction of GROMACS hydrogen-bond analysis. GROMACS uses an angle at the donor, whereas MDAnalysis uses the D–H–A angle and defaults to a 3.0 Å donor–acceptor cutoff. N/O/S element screening alone does not determine protonation-specific acceptor chemistry. Bond-bearing topology is preferable to distance-only guesses for donor–H pairing. Sources: [GROMACS geometric definition](https://manual.gromacs.org/2021-current/reference-manual/analysis/hydrogen-bonds.html), [MDAnalysis hydrogen-bond API](https://docs.mdanalysis.org/2.8.0/documentation_pages/analysis/hydrogenbonds.html).

## Completed numerical checks

The analytic fixtures check the app's public storage and analysis boundaries, independently of a protein simulation. A 0.1 nm displacement becomes 1 Å in the display buffer while original coordinates remain in nm. Canonical atom names and indices survive the PDB storage round trip. Known test coordinates give 1 Å, 90° angle, and +90° dihedral. A periodic pair separated by 9 Å in a 10 Å cubic cell measures 1 Å under the minimum-image convention. A two-frame D–H–A fixture gives 3 Å then 4 Å and 50% geometric occupancy; a hydrogen assigned to the wrong donor is rejected. Collinear torsions, duplicate picks, and indices outside the topology are rejected.

PDB, GRO, XTC, TRR, NetCDF, MDTraj HDF5, XYZ, DCD, and OpenMM-written mmCIF coordinate round trips passed. These are small representative fixtures, not exhaustive parser validation. XTC/TRR/NetCDF/HDF5 physical times survived the fixture round trips. Review identified that the installed MDTraj DCD loader synthesizes frame indices instead of recovering physical time. DCD imports now carry an unavailable-time warning, label the horizontal axis as frame indices, and accept an explicit original-frame interval; that override was checked. The initial missing PyTables dependency for HDF5 and unsupported iterative mmCIF loading were corrected. Same-count atom-name mismatches in self-describing PDB trajectory files are rejected. Source for DCD reader behavior: [MDTraj DCD implementation](https://github.com/mdtraj/mdtraj/blob/master/mdtraj/formats/dcd/dcd.pyx).

A skewed triclinic cell fixture also returns its analytically expected minimum-image separation. Chart checks verify nonuniform timestamp spacing and nearest-time seeking. Torsions use a labeled circular mean: +179° and −179° average to 180°, not 0°; a zero circular resultant is reported as undefined. Curve segments break across ±180° wraps. Unordered timestamps explicitly fall back to a frame-index plot. These fixes prevent rendering conventions from being mistaken for physical changes.

Renderer inspection confirms that it validates atom count and atom/residue names against the canonical topology, writes interpolated coordinates to a separate buffer, and uses saved-frame backend measurements for numeric labels. Periodic trajectories currently play saved frames, avoiding unverified interpolation across cell boundaries. This is a deliberate limitation of the initial viewer. Interpolating Cartesian coordinates, even for nonperiodic data, is an animation operation and does not add physically sampled conformations.

## Real engine integration evidence

| Run | System and model | Recorded output |
| --- | --- | --- |
| OpenMM 8.6.0, job `1e5a626039224d3c` | Ubiquitin, ff14SB + GBn2 implicit solvent, Langevin-middle, 300 K, 2 fs, seed 2026, 0.2 ps initial relaxation | 101 frames, 1,231 atoms, 0–2 ps production |
| GROMACS 2025.4, job `b3898c9d78e94214` | Ubiquitin, Amber99SB-ILDN + TIP3P explicit solvent, stochastic dynamics, 300 K, 2 fs, seed 2026, 0.1 ps initial relaxation | 6 frames, 26,842 atoms, 0–0.1 ps production |
| OpenMM explicit-water check, job `a6d4e0ed2f124807` | Ubiquitin, ff14SB + TIP3P with periodic PME; isolated integration-test storage | 3 frames, 18,556 atoms, 0–0.02 ps production |
| GROMACS solvent-selection check, job `9963ef95444b447e` | Repeated short explicit-water workflow after restricting neutralization to newly added waters | 6 frames, 26,842 atoms, 0–0.1 ps production |

All four completed through the real local engines. Inspection verified finite coordinates, increasing timestamps, topology/coordinate counts, and all 111 recorded output hashes across the runs. The initial OpenMM serialized system has the expected particle count and integration settings, and both OpenMM runs have finite reported energies. The GROMACS commands completed without a `grompp -maxwarn` override. Exact configuration, logs, inputs, trajectories, checkpoints, and provenance remain in the corresponding local job directories. Evidence is in `audit/engine-smoke-checks.json`. These different force fields and solvent models do not constitute a controlled engine comparison or independent replicas for convergence assessment.

The neutral ubiquitin full runs did not exercise nonzero ion replacement. A separate GROMACS integration check forced one sodium and one chloride replacement within the new-solvent-only selection, preserving all 1,231 input atoms. That check caught and corrected GROMACS's requirement for the solvent group name `SOL`; the affected branch was then verified directly. Actual OpenMM worker cancellation was also exercised. Those supplementary result records are captured in `audit/engine-smoke-checks.json`, with underlying artifacts in `data/integration-checks/`.

## Preparation and interpretation

The application should preserve input files and record every preparation transformation. Unsupported residue chemistry should fail clearly rather than silently removing ligands or inventing parameters. Automatic hydrogen placement and protonation are assumptions: OpenMM's Modeller selects protonation using the requested pH and structural context, and does not replace a system-specific assessment. A trajectory topology after preparation can legitimately have new atoms; it must receive its own consistent metadata and coordinates. Source: [OpenMM model building](https://docs.openmm.org/latest/userguide/application/03_model_building_editing.html).

Fixed engine presets must disclose their actual protein force field, solvent and water parameters, electrostatics, constraints, integrator, thermostat, time step, seed, minimization, equilibration, and production duration. OpenMM supports Amber protein fields with specific compatible water/ion XML files or separately selected Generalized Born models; those choices describe different physical models. A short implicit-solvent preview and explicit-solvent GROMACS simulation are not a controlled engine comparison. Source: [OpenMM force-field and solvent documentation](https://docs.openmm.org/latest/userguide/application/02_running_sims.html).

## Reproducibility evidence

For each real run, preserve the original input and checksum, prepared topology, exact settings, seed, engine and package versions, engine platform, command/log records, trajectory, and output checksums. Seed recording supports provenance but does not promise bitwise repeatability across platforms. Downloadable artifacts should be enough to inspect what ran; a browser movie alone is insufficient evidence.

The bundled OpenMM demo was generated before Python/NumPy preparation RNG controls were added. Its integrator and velocity seed is known, but hydrogen-preparation randomness was uncontrolled. This limitation is recorded in its provenance. Prepared coordinates, serialized system/integrator, original source, and outputs are retained; exact replay of that initial preparation is not claimed. Subsequent runs seed preparation as well as engine RNGs. Original uploaded files and their names, hashes, stride, and user time override are retained for subsequent imports.

The complete indexed public Practical Cheminformatics corpus was refreshed during this review: 91 records, 74 archived Blogger posts and 17 current GitHub Pages posts, with guest authorship preserved. Searches used MD, trajectories, reproducibility, validation, and code-sharing terms. The directly relevant result is Pat Walters's [“Where's the code?” (2019)](https://practicalcheminformatics.blogspot.com/2019/05/wheres-code.html), which advocates accessible code and data so that others can reproduce and inspect computational methods. Applying that principle here means preserving executable code, exact inputs, configurations, and outputs. This application to DynaMol is review synthesis, not a statement by Walters about this project. Corpus metadata, attributed search results, and audit reports are in `docs/audit/`.

## Remaining scientific validation

Independent replicas, sensitivity analysis, equilibration and convergence assessment, and uncertainty estimates are outside the bundled smoke demonstration. They remain prerequisites for claims based on a scientific simulation campaign. No predictive model, training data, held-out set, docking ranking, compound design, or statistical comparison is being evaluated in this software task.

## Final packaged application check

The one-command launcher was exercised against the frozen dependency lock. Its final environment has OpenMM 8.6.0, MDTraj 1.11.1.post2, and NumPy 2.4.6. All 19 backend tests and all eight independent numerical groups passed again in that environment. Seven production-browser integration tests passed, including actual molecular rendering, live simulation completion, import errors, distance/hydrogen-bond measurement and CSV export, race handling, and mobile/tablet layouts. These checks validate tested software behavior; the scientific convergence/uncertainty limitations above remain. Local data created by browser tests were archived under `data/integration-checks/browser-tests/`.
