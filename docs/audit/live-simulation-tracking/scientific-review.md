# Live simulation tracking: evidence-linked review

This is a language-model code and evidence review, separate from the implementation task, plus an automated PatAgent checklist. It is not a domain-expert certification, independent experimental validation, or Pat Walters endorsement. PatAgent retrieves public recommendations; it does not represent Pat Walters.

## Preflight decision

`review-preflight.json` reports **proceed: six passes, no failures or unresolved gates**. `review-manifest.preflight.json` records parsed finite input fixtures, SHA-256 hashes, installed versions, seed 2026, two CPU threads, and bounded native tests. This permits the specified software checks. The frozen backend and final UI have passed the scoped checks below. A preflight pass does not establish scientific convergence.

The endpoint is a live display of distance (Å), angle/torsion (degrees), and an explicitly defined geometric hydrogen-bond screen on saved physical frames. Display interpolation is excluded. No affinity, biological function, equilibrium occupancy, relative engine accuracy, or converged dynamics claim is intended. This task introduces no new chemistry, force-field parameters, protonation choices, or production MD settings.

## Required correctness checks

| Risk | Required evidence before software promotion |
| --- | --- |
| Wrong atoms after preparation | Source-to-native mapping retains atom identity despite added/reordered atoms; missing, duplicate, renamed or rebuilt identities are rejected or disclosed. GROMACS rebuilt hydrogens must not inherit a selected input hydrogen by index. |
| Periodic artifacts | Compare analytic boundary-crossing controls and live/offline calculations using the same saved unit-cell vectors. No inferred box or interpolated coordinates. |
| Undefined geometry | Coincident bonds and collinear torsion planes remain null at their original time positions, with a reason. JSON must not contain NaN/Infinity or a plausible fabricated zero. |
| Hydrogen-bond interpretation | Validate D,H,A order, explicit donor–hydrogen connectivity, and the geometric thresholds. Distances/angles are geometric observations; element/connectivity checks alone do not prove donor/acceptor chemistry. |
| Interrupted/resumed runs | Read only complete samples, retain step/time alignment, ignore incomplete XTC/JSON tails, and remove samples outside the committed checkpoint prefix. Resume does not duplicate its boundary sample. |
| Monitor failure changes MD | A missing mapping, snapshot read failure or undefined measurement cannot alter forces, engine configuration, or completion status. Native samples and reproducibility records remain accessible. |
| Completion replaces another view | Completion is scoped to the originating dataset/session and load generation. Another dataset currently being inspected retains camera, selections, plots and water visibility. Explicit Explore remains available for the completed job. |
| Promoted plot/camera mismatch | Remap selections to final topology, retain monitor IDs and appearance, carry the actual camera state, and suppress default-distance injection during restored-workspace loading. No positional-array copying across unrelated topology. |

The review found and verified fixes for chain identity, ordinary GRO coordinate rounding, identity-proof persistence through Resume, and completion-state restoration. The final UI passes a complete remapped workspace, suppresses default-plot injection, scopes late visibility and warning updates to the current load/session, and restores tracked plots for automatic and explicitly opened results. Plot restoration is sequential to bound physical-trajectory memory use. Camera translation compensates for current engine box recentering while preserving rotation and zoom. `review-findings.json` records the original observations and their resolution evidence.

## Sources and scope of inference

The cached Practical Cheminformatics corpus contains 74 archived Blogger and 17 GitHub Pages posts; the manifest checksum and three task-specific searches are retained in `review-corpus-search.json`. This records coverage of that cache, not a claim that no newer post exists.

Pat Walters argues that benchmark structure/representation quality and well-defined endpoints require explicit review in [We Need Better Benchmarks for Machine Learning in Drug Discovery](https://practicalcheminformatics.blogspot.com/2023/08/we-need-better-benchmarks-for-machine.html) (2023, Pat-authored; local corpus used because the live Blogger page did not fetch). [The Trouble With Tautomers](https://patwalters.github.io/The-Trouble-With-Tautomers/) (2025, Pat-authored) distinguishes algorithmic molecular representations from physically relevant states. Applying these principles to persistent atom identity and avoiding silent hydrogen reassignment is this reviewer's inference, not a quoted MD prescription.

For the actual geometry convention, primary [MDTraj distance documentation](https://mdtraj.readthedocs.io/en/stable/api/generated/mdtraj.compute_distances.html), [angle documentation](https://mdtraj.readthedocs.io/en/stable/api/generated/mdtraj.compute_angles.html), and [dihedral documentation](https://mdtraj.readthedocs.io/en/stable/api/generated/mdtraj.compute_dihedrals.html) state that periodic computation uses available unit-cell information and minimum images. Angles/torsions are returned in radians and require explicit conversion for degree displays. These documented semantics must be corroborated by the installed-version tests, especially for degeneracies and truncated trajectory reads.

## Verified backend evidence

The native report records 115 passing backend regressions. A 400-step OpenMM job produced 17 saved frames, and a 3,000-step GROMACS job produced 31; both were stopped and resumed, showed values while running, and matched offline distance, angle and dihedral calculations exactly on the same saved frames. A separate 50-step OpenMM run retained an explicit donor–hydrogen bond over six saved frames and matched offline donor–acceptor distances, donor–hydrogen–acceptor angles, times and geometric occupancy. These are consistency checks of stored observations, not independent scientific accuracy estimates. The job worker hashes match the reviewed source. See `native-validation.json` and the retained isolated job directories.

Independent review checks reproduced and then verified fixes for same-sequence chain permutation, ordinary GRO rounding, and preservation of a failed identity verdict on a new/resumed worker. Four focused identity/truncation tests and three failure/reconstruction/workspace tests passed. Separately, `review-analytic-geometry.json` specifies four known-coordinate frames and independent expected values for all four measurement kinds, including periodic crossings, undefined torsions, coincident bonds, and hydrogen-bond pass/fail/null cases. Thus the check is not limited to comparing a function with itself. `review-identity-regressions.json` retains the original identity reproductions.

Identity verification is persisted, checked against the input hash and native atom indices, and included in checkpoint dependencies. Runtime tracking failure reports an unavailable monitor while MD continues. OpenMM reconstructs from committed checkpoint frames; GROMACS verifies its last consumed physical-frame signature and rebuilds when an append operation truncates or regenerates the tail.

## Verified browser and presentation evidence

`browser-regression.json` records **30 passing tests, no failures and no skips** against the isolated real API. This includes all four OpenMM completion-handoff cases, tracking persistence across mode changes and reload, undefined plot values and torsion branch cuts, existing pinch/zoom/fullscreen controls, and a real six-residue preparation plus explicit water build. The tracked heavy-atom selection survived preparation and solvation; 3,648 water atoms were created, and hide/show changed the canvas. Water setup is visible before run settings and describes a ready saved box without asserting that waters are currently visible.

Four GROMACS completion-handoff cases also passed. These browser handoff controls replay only delivery timing of actual native result data through the real API; they are not additional live MD runs. Actual live observations and Stop/Resume were verified separately in `native-validation.json`. The browser controls verify all three mapped native plots, physical timestamps, final frame, retained camera scale/rotation, reload persistence, manual result opening, and that another scene's camera, plots and visibility remain unchanged. A completion already fetching coordinates cannot overwrite a new setup.

`review-camera-control.json` independently checks a translated two-frame structure against the installed NGL 2.5.0 orientation convention: molecular center, rotation and zoom are preserved, and input camera state is not mutated. It records explicit input coordinates, output camera and the final helper hash. This supports the current engines' recentering behavior; arbitrary rigid reorientation or preserving biological focus after unvalidated correspondence changes is outside this check. Final reviewed source hashes match `browser-validation.json`.

## Evidence limit and review decision

No unresolved software finding remains within this task's scope. The final postflight checklist records **seven passes and one failed uncertainty gate** (`review-postflight.json`). No independent-replica or scientific confidence interval is supplied by these bounded software checks, so the generic scientific-promotion gate remains unsatisfied. The PatAgent result is therefore **block for scientific promotion**, not an all-gates pass. The separate scoped software checks above passed.

A stable trace or complete software test run must not be presented as equilibrium, convergence, accurate population/affinity prediction, universal chemical coverage, or proof that one engine is scientifically superior. Hydrogen-bond occupancy is the geometric pass fraction of defined saved frames. The current GROMACS adapter rebuilds input hydrogens and explicitly refuses live monitors selecting those hydrogens. These limits remain visible in the implementation and evidence.
