# DynaMol release-readiness review

This is a language-model review of source code and recorded tests, using an automated evidence checklist. It is not expert certification, independent experimental validation, or an endorsement by Pat Walters.

The reviewed scope is workspace restoration and portable projects, recoverable simulations, readiness checks, library management, live run diagnostics, RMSD/RMSF and named selections. Large-trajectory streaming is outside this change. Test results support the bounded software behavior below; they do not establish universal chemistry coverage, correct biological states or converged molecular dynamics.

## Evidence and interpretation

| Area | Recorded result | What it establishes |
| --- | --- | --- |
| Backend regression suite | [251 passed; 2 dependency deprecation warnings; 14.25 s](backend-tests.txt) | Included chemistry, identity, recovery, workspace, resource, geometry and diagnostic guards |
| UI coverage union | [82 distinct cases; latest results all passed, none failed/skipped](ui-coverage-union.md) | Isolated browser workflows across separate integration runs; affected cases were rerun after fixes |
| Structural-analysis UI | [4 passed; 9.779 s](analysis-ui.json) | Real RMSD/RMSF responses, CSV/JSON exports, seeking, atom selection, readiness transitions and desktop/narrow layout |
| Native chemistry | [43 positive cases](chemistry-matrix.json) and [6 boundary/pipeline checks](chemistry-boundaries.json) passed | Representative native parameter consistency, finite short integration and preserved/rejected chemistry |
| Native checkpoint recovery | [Both OpenMM and GROMACS passed](recovery-native.json) | Deliberate interruption, compatible continuation, unchanged retained frame prefixes, unique increasing times and finite observations |
| Recovery interface | [1 native OpenMM browser workflow passed; 31.7 s](recovery-validation.json) | Live temperature, Stop, Resume, partial Files and completed trajectory opening |
| Early cancellation | [4 immediate startup checks](cancellation-startup.json) and [2 original browser cancellation cases](workspace-cancellation-ui.json) passed | Persisted cancellation intent survives worker exit; stopped jobs retain usable records |

The [machine-readable union](ui-coverage-union.json) retains the earlier five failing legacy UI cases and their later passing resolutions. It is not a claim that all 82 cases ran in one command or on one final source revision. Historical product hashes that were not recorded are explicitly null; the final source hashes are not assigned retroactively. Every referenced run-report hash was checked against its actual bytes for this review. The [aggregate manifest](aggregate-manifest.postflight.json) records source/evidence hashes and the exact scope.

The individual test groups overlap; their counts must not be summed into an independent benchmark size. Chemistry checks used seed 2026 with recorded versions and native inputs/outputs. The matrix includes all 20 standard amino acids across 28 residue/state/terminus/disulfide fixtures, five ion models, and ten organic ligand/cofactor fixtures. The initial failing matrix is retained alongside the corrected matrix, so this is an iterative regression audit, not a blind benchmark.

The five-step chemistry integrations used 0.2-fs steps and nonperiodic Reference-platform systems. Finite numbers and native parameter consistency at these geometries demonstrate numerical continuity only. They do not validate aqueous thermodynamics, equilibrium populations, coordination chemistry, redox behavior, binding affinity, or long-run stability. Checkpoint tests likewise establish state/file continuity on short local runs, not identical future trajectories across hardware or runtimes.

## Numerical and interface review

RMSD uses original saved coordinates and an explicit reference, atom selection, alignment selection and frame window. The fit is a proper Kabsch rotation with a determinant correction, so mirror images cannot be aligned through reflection. Float64 fitting covers the tested planar three-atom limiting case; fewer than three non-collinear fit atoms produce a clear error. Values use equal atom weighting and convert nanometers to angstroms.

RMSF is the square root of population mean squared displacement about each atom's window mean after the stated fit. Residue aggregation takes the square root of mean atomic variance. The implementation uses stable per-atom Welford accumulators and retains original coordinates. Periodic processing makes bonded molecules whole and chooses nearby images; it is explicitly not continuous diffusion unwrapping. Turning alignment off retains rigid translation and rotation in the result. The independent analytic fixtures and comparison with MDTraj constrain implementation errors; neither provides biological validation.

Recorded diagnostic temperatures are actual kinetic/native observations. Missing temperature does not become zero or the thermostat target. The reader drops partial final rows, reports malformed observations and downsamples while retaining extrema. Estimated remaining time is labeled approximate and includes setup. Final diagnostics and the original files remain available after a run stops.

Readiness uses the same structure/configuration validation and resource limits as submission. The resource estimate is a conservative planning bound, not a guaranteed memory reservation or runtime benchmark. A green readiness state means the current implementation accepts the selected input/model; it is not a scientific quality certificate.

The final integration resolved concrete findings rather than suppressing their symptoms:

- A planar quaternion-fit failure was replaced with the proper SVD/Kabsch fit and explicit degeneracy checks.
- Infinite/overflowed, negative and fractional diagnostic steps, and negative/non-finite times, now become reported invalid observations instead of causing an uncaught server error.
- Non-increasing timestamps now use an explicit saved-frame axis with a warning instead of collapsing the plot domain.
- A companion **Analysis record** JSON now retains the exact calculation request, measured/fit atoms, frame mapping, method and warnings; the separate CSV remains a simple table.
- Browser checks verified the analysis pane's corrected flex layout and full visible controls, and workspace restoration now preserves expanded analysis state because it changes canvas dimensions.
- An early Stop race now respects the persisted cancellation marker after the worker exits; original native cancellation assertions and added backend guards pass.

See the [analysis contract](../../ANALYSIS.md), [workspace contract](../../WORKSPACES.md) and [recovery contract](../../CHECKPOINT_RECOVERY.md).

Workspace identity checks include atom ordering and physical-coordinate hashes. Imported project checksums detect changed bytes, but do not authenticate authorship or certify molecular meaning. Backups preserve source ancestry and preparation parameters; historical runs imported from a backup are not silently queued as live work. Coordinate and archive size limits remain explicit.

## Chemistry applicability and exclusions

The supported path includes ff14SB standard protein states/caps, curated SEP/TPO/PTR/HYP states, compatible noncovalent organic GAFF2/AM1-BCC ligands, and the tested TIP3P-compatible Na+, Cl−, K+, Mg2+ and Ca2+ models. The ATP/ADP/NAD/FAD fixtures passed their fixed-state native checks. Cofactor names alone do not establish a transferable redox/protonation model, and monatomic Mg/Ca parameters do not establish correct coordination energetics.

Unsupported modified amino acids remain identifiable and retained. MSE/ALY and other unregistered covalent residues require dedicated compatible parameterization; arbitrary parent-residue substitution or disconnected small-molecule parameters are not a valid substitute. Heme, organometallics, covalent ligands/cofactors, radicals, incomplete/ambiguous ligand graphs and unsupported crosslinks receive explicit blockers. The 1UA2 TPO fix does not resolve its four 12-residue internal sequence gaps: these exceed the local six-residue loop-builder limit. Prepared complex/PTM/ion state transfer to GROMACS remains unsupported; those supported preparation paths use OpenMM and explicit water. See [chemistry support](../../CHEMISTRY_SUPPORT.md) for the full matrix and state assumptions.

## Automated postflight and claim boundary

The [aggregate automated postflight](aggregate-postflight.json) reports **7 passes, 1 failure and 0 unresolved gates**, with decision **`block`**. The failing critical gate is `uncertainty`: scientifically meaningful uncertainty intervals are explicitly absent. The tool's own audit/parser smoke test passed, including deliberately flawed fixtures, but this only checks the checklist implementation.

The absence of uncertainty intervals is deliberate and disclosed. Deterministic implementation checks do not supply meaningful confidence intervals for biological accuracy. The failed gate therefore remains a block on accuracy/convergence or predictive-model promotion; it is not relabeled as a pass to obtain a release badge. No statistical efficacy comparison, prospective experiment, production replicate/convergence analysis, or exhaustive applicability-domain study was performed.

The cached full public Practical Cheminformatics corpus contains 91 records at the recorded sync: 74 archived Blogger posts and 17 GitHub Pages posts, with guest attribution retained. [The Trouble With Tautomers](https://patwalters.github.io/The-Trouble-With-Tautomers/) supports distinguishing a chosen representation from experimentally populated states. [We Need Better Benchmarks for Machine Learning in Drug Discovery](https://practicalcheminformatics.blogspot.com/2023/08/we-need-better-benchmarks-for-machine.html) supports preserving structure identity, provenance, relevant endpoints and explicit benchmark boundaries. Applying these general principles to this MD software audit is this review's inference. [Retrieval evidence](release-review-blog-search.json) records the searched terms and canonical sources.

Public distribution still requires separate packaging acceptance, including platform signing/notarization, clean-machine launch validation and complete redistribution materials as applicable. Bundle-specific checks are recorded separately and do not expand the chemistry or convergence evidence. Native Quit-menu interactions were not exercised because the desktop was locked; Swift validation and local engine/browser tests do not replace that UI check. This review does not certify a public release artifact.
