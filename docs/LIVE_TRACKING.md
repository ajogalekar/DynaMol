# Live simulation tracking

Simulation setup is organized as engine, structure/preparation, water/environment, and run settings. Tracking is optional. Expand **Track during simulation** to select existing measurements or add a distance (2 atoms), angle (3), dihedral (4), or geometric hydrogen bond (donor, bonded H, acceptor). Completed picks appear immediately in both the setup list and the bottom plot. Definitions remain with saved workspaces and projects.

A run accepts at most 12 observations. They do not apply forces or restraints. Selections are fixed for that run; checkboxes unlock when it stops. Live curves contain saved production frames, not integrator substeps or equilibration. The starting structure remains in view while curves grow. The same bottom panel stays present when switching Simulate and Explore. Live plot clicks do not seek the static starting structure; completed plots regain click-to-seek.

When a followed run completes, it opens in Explore in the same workspace. Saved plot IDs, labels, colors, visibility and tracked choices carry over, with indices verified against output topology. Final values are calculated once from the saved physical trajectory, in bounded sequence, so timestamps and frame counts agree with playback. Camera rotation and zoom persist; the scene translation compensates for water-box recentering. Display coordinates and physical analysis are never modified to make the handoff appear smoother.

Starting New simulation or choosing **Run in background** cancels automatic opening. If another molecule is displayed at completion, that scene, its plots, camera and visibility remain unchanged. The result notice and background job's **Open trajectory** open the completed run in Explore with its own tracked plots. A late network response cannot overwrite a new setup. The followed job is remembered in local browser storage; monitored series and definitions remain in the job files.

## Identity and geometry

Preparation and solvation remap retained selections using unique chain, residue number/name, atom name and element correspondence. Missing or ambiguous atoms produce a visible request to pick again. Hydrogen identities changed by protein preparation must be selected again. Engine setup pins correspondence before integration, including format-aware GROMACS coordinate precision, and hashes the identity proof into checkpoint recovery. GROMACS's current adapter rebuilds hydrogens, so it rejects hydrogen-containing monitor selections and supports heavy-atom distance/angle/dihedral tracking. OpenMM supports explicit donor–hydrogen–acceptor tracking with validated connectivity.

Analysis uses original physical coordinates and minimum-image geometry when valid periodic box data is available. Degenerate measurements retain their timestamp with a null value and per-frame reason; plots show gaps and CSV cells remain empty. Dihedral plots break at the angular branch cut and report circular means. Hydrogen-bond occupancy is a geometric cutoff summary over valid frames, not evidence of chemical donor/acceptor eligibility or binding.

Snapshots are atomic bounded JSON files. OpenMM updates from saved frame records; GROMACS reads only complete production XTC frames and tolerates a partial tail. Stop/resume rebuilds the retained prefix without duplicate boundary samples. Tracking failures are reported separately and do not terminate an otherwise healthy simulation.

## API

- `POST /api/jobs`: optional `measurements` array, each with `id`, `kind`, source-dataset `atoms`, `label`, and optional `color`.
- `GET /api/jobs/{id}/measurements`: definitions, source/output atom mappings, nullable values, physical times, frame errors, warnings and job status. Polling does not reread the trajectory.
- `POST /api/datasets/{target}/remap-measurements`: `{source_dataset_id, measurements}` returns mapped definitions plus errors for lost/ambiguous selections. Up to 100 existing workspace plots can be remapped; a run still accepts at most 12.

## Validation scope

The [native validation record](audit/live-simulation-tracking/native-validation.json) covers actual bounded OpenMM and GROMACS production runs, live observations, stop/resume and exact comparisons to offline saved-frame analysis. It includes an OpenMM hydrogen-bond case. Backend regression tests cover periodic geometry, undefined values, identity ambiguity, frame-prefix recovery and persistence. Browser tests cover setup controls and replay the delivery timing of those native results against the real API/viewer to check automatic/manual handoff, mode continuity and focus protection. These checks validate software consistency; these short runs do not establish equilibration, convergence, biological relevance or predictive accuracy.
