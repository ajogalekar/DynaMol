# DynaMol UI functional audit

Audit date: September 10–11, 2026. Application: local DynaMol at `http://127.0.0.1:8765/`.

The workspace/recovery/analysis pass has **82 distinct browser cases with passing latest results and zero skips**, assembled from documented isolated runs, plus **251 passing backend tests**. See [the coverage union](audit/release-readiness/ui-coverage-union.md) and [release additions](RELEASE_READINESS.md). Earlier counts below are historical snapshots; their results have not been relabeled as current tests. Current browser regressions require a disposable data root and refuse user-facing ports; use [frontend/TESTING.md](../frontend/TESTING.md).

## Scope and method

This audit inventories distinct interactive controls from `App`, `ImportDialog`, `StructureWorkbench`, `SimulationPanel`, `StructureJobMonitor`, `MolecularViewer`, and `PlotPanel`. It combines the existing browser regressions with additional real Chrome/Playwright interaction tests. Molecular rendering uses the actual NGL/WebGL viewer; structure imports, measurements, preparation, solvation, and dynamics use the local API and installed engines. Camera, representation, hydrogen-context, and group-visibility checks inspect actual canvas output rather than trusting selected-button state alone.

New audit tests are in `frontend/tests/ui-audit-*.spec.ts`. Small format fixtures, a regeneration script, and provenance are in `docs/audit/ui-audit-fixtures/`. These contain five original saved frames of the bundled ubiquitin demo, restricted to six residues. `mixed-groups.pdb` includes arbitrarily placed ethanol, water, and sodium for visibility testing only; it is not a simulation-ready scientific model. MOL/SDF/SMILES fixtures contain ethanol. The LAMMPS fixture uses an artificial cubic cell for decoding tests. The schematic nucleic-acid fixture is a four-residue backbone without bases or validated stereochemistry, used only to test polymer visibility.

The audit uses short CPU jobs and never runs multiple application jobs concurrently. Downloads are opened and inspected locally. Invalid-input, delayed-response, unavailable-fullscreen, and context-loss cases are explicitly identified as fault injection; they do not substitute simulated molecular results for real backend calculations.

## Defects found and corrected

| Finding | Reproduction | Correction and regression |
| --- | --- | --- |
| mmCIF import returned HTTP 500 | Open files → choose a valid mmCIF → Open in DynaMol. MDTraj retained string atom IDs and its PDB exporter attempted numeric modulo on them. | Canonical export now gives a copied topology integer serials. The backend regression checks exact physical coordinates, atom identities and ordering, preserved original CIF bytes, and original-file SHA-256. Browser format coverage includes mmCIF. |
| Escape did not exit picking after pressing the pick button | Click Pick atoms in the scene; press Escape while that button retains focus. | Escape is handled before the form-control shortcut guard. A browser regression exercises that exact focus state. |
| Keyboard picking left the trajectory playing | Start playback, focus the scene, press M. | Entering picking via M now pauses playback, matching the visible pick button. A browser regression verifies the movie stops. |
| Hidden atoms retained gold selection spheres and labels | Select an atom, then hide its molecule group or hydrogen class. | Selection representations now use the same visibility predicate as molecular geometry. The inspector retains the selection, while a six-case real-canvas regression verifies its hidden marker disappears. |
| Ligand-state overrides were submitted with the wrong media type | Enter a valid explicit-state SMILES in a ligand card. | The inspection POST now declares JSON. The source audit verifies actual updated ligand chemistry and requested pH. |
| Bundled-demo inspection failed during complex integration | Open Studio on the bundled trajectory, or Retry inspection after a network error. | Parent-source traversal now accepts the demo's metadata shape. The inspection-retry browser regression passes against the real service. |
| LAMMPS trajectories were absent from the chooser filter | Open files → Add trajectory did not list `.lammpstrj` as an accepted extension. | Added the extension and a visible LAMMPS format caption; a real file import verifies its coordinates, stride, and timing. |

Native GROMACS also produced a clearly reported failure when minimization was deliberately disabled for a newly solvated system. Repeating the smoke run with normal minimization completed. This was a scientifically unsuitable test configuration, not a claim that every user-selected MD configuration must converge or remain stable.

## Control inventory and evidence map

Each row groups parallel controls; named actions were individually exercised unless the notes explicitly limit the claim.

| Area | Controls / behavior checked | Evidence |
| --- | --- | --- |
| Workspace navigation | Home link; Explore; Simulate; New simulation; Studio close; library open, dataset selection, Escape close; structure details expand/collapse | `ui-audit-controls`, existing `e2e` dataset-switch race regressions |
| Help | Header quick guide, scene help, footer keyboard shortcuts; close button and Escape; focus trap and restored focus; four documentation links | `ui-audit-controls`. Each link is clicked and its new-tab target checked; external documentation availability is outside the local-app pass criterion. |
| Representation | Ribbons, Ball & stick, Sticks, Surface; sequence, chain, element color modes | `ui-audit-controls`, actual scene screenshots for every option |
| Molecular visibility | Protein/nucleic-acid toggle; ligand toggle; water toggle; ion toggle; Hidden, Polar only, All hydrogens; hidden selection highlights | `ui-audit-imports` mixed and schematic nucleic fixtures, `hydrogen-context`, existing preparation/water pixel comparisons. |
| Camera | Zoom in/out; fit; auto rotate on/off; drag rotation; right-drag pan; wheel zoom; residue number/name/chain queries and no-match feedback; real atom hover and double-click focus; full screen enter/exit | `ui-audit-controls`, existing `e2e`. Progressive edge antialiasing is allowed to settle before judging stopped rotation. |
| Snapshots | Save snapshot; PNG download; disabled while viewer unavailable; saved notification and dismissal | Existing `e2e`, `ui-audit-recovery` |
| Playback | Play/pause; first/previous/next/last frame; boundary clamping; timeline drag/seek; each speed (0.25, 0.5, 1, 2, 4×); looping and stop-at-end; single-frame disabled playback | `ui-audit-controls`, existing `e2e` and source-loading tests |
| Keyboard | Space, arrow keys, F, M, Escape, ?; form controls retain their own arrow-key behavior; Escape with button focus; M while movie playing | `ui-audit-controls` |
| Inspector | Toggle/close inspector; Pick atoms button; Stop picking; real canvas atom click and repeat-click deselection; atom search, no results, clear search; selected-atom remove; clear all; plot disabled for incomplete selection | `ui-audit-controls` |
| Measurements | Distance, angle, dihedral and donor–H–acceptor hydrogen-bond analysis; invalid hydrogen-bond error; current numeric values; occupancy; saved-frame count | Existing `e2e` (distance/H bond), `ui-audit-controls` (angle/dihedral/errors). Geometry correctness also has existing backend tests. |
| Plots | Active tab switch; per-measurement show/hide; global show/hide without deleting plot; remove; CSV export; hover readout; click-to-seek; add/Measure actions; disabled export with no measurements | `ui-audit-controls`, `preparation.e2e`, existing `e2e` |
| General import | Open files; choose topology and trajectory file dialogs; drag/drop; no-file disabled submit; close/cancel; invalid file and atom-count mismatch; stride bounds; explicit timing override | `ui-audit-controls`, `ui-audit-imports`, existing `e2e` |
| Molecular formats | PDB, mmCIF, GRO, MDTraj HDF5 structures; XTC, DCD, TRR, NetCDF, XYZ, MDCRD, LAMMPS trajectories; multi-frame stride and timing preservation | `ui-audit-imports`; real files, decoded atom/frame counts, and saved-frame distance checks against original source coordinates. Filename aliases such as `.mmcif`, `.pdbx`, `.hdf5`, `.netcdf`, and `.crd` share loaders and are not separately certified. |
| Studio sources | Upload/Fetch/SMILES tabs; file chooser and source loading in view; PDB, MOL2; invalid and valid SMILES; molecule name; RCSB and PubChem database selection and real fetches; invalid fetch and empty-button states | `preparation.e2e`, `ui-audit-controls`; supplementary file variants covered by the audit source tests |
| Preparation settings | Requested pH; add missing heavy atoms; side-chain refinement; water-removal selection; explicit heterogen-removal selection; known short-gap build opt-in and disabled unknown-gap state | Existing preparation suites and `ui-audit-controls`. This validates UI/configuration behavior; molecular preparation limitations remain documented separately. |
| Preparation progress | Immediate submitting feedback; active spinner/stage/progress/elapsed time; persistent monitor across closing/reopening Studio; logs; completed structure auto-load; success dismissal; actual rejection; failed display followed by same-result retry; native Cancel and prepared-PDB download | `preparation-progress`, `preparation-jobs`, `ui-audit-jobs` |
| Solvent | Implicit/explicit selection; real visible TIP3P box; water toggle changes actual scene; implicit returns to dry prepared parent; saved state tied to correct source; Update box at changed padding; setup Cancel/log/dismiss controls | `preparation-jobs`, `ui-audit-jobs` |
| Complex-specific controls | Ligand detail expand/collapse; current charge and chosen SMILES; explicit-state override; pH reinspection; removal opt-in disables override; inspection failure Retry; prepared parameters download | `ui-audit-sources`; native complex release evidence is recorded separately. |
| Simulation configuration | OpenMM/GROMACS engine buttons and installed state; name; duration; temperature; solvent restrictions; advanced panel; integration step; save interval; seed; friction; equilibration; box padding; minimization checkbox; step/frame summary | `ui-audit-controls`, existing `e2e`, native job audit |
| Background jobs | Real OpenMM and GROMACS Start; progress to 100%; view/hide log; background footer reopens monitor; Stop cancellation; disabled duplicate start; completed output loading; Files ZIP download inspected for native trajectory and configuration/status files | `ui-audit-jobs`, existing `e2e` |
| Errors and recovery | Invalid input retains scene; dismiss error and notification; failed initial example offers fallback/import; viewer context-loss feedback and reload recovery; unavailable fullscreen notification | `ui-audit-controls`, `ui-audit-recovery`, existing viewer-retry regression |
| Responsive use | Desktop 1440×1000, tablet 820×1180, phone 390×844; reachable representation/inspector/source/advanced controls; vertical Studio split; no horizontal page overflow | `ui-audit-recovery`, existing responsive `e2e` |

## Validation status

The combined release run passed **59 of 59 browser tests in 145.6 seconds**, with no failures, skips, retries, or flaky results. This includes 17 existing regressions and 42 new audit tests across 10 files. The 22 inventory rows above cover the implemented interface, including all 11 listed trajectory/structure format families, preparation and water cancellation, a real GROMACS completion, a real OpenMM completion and cancellation, and the prepared whole-complex downloads. The complex-download case used the actual successfully prepared dataset `4cf5bcedf38847c1`; it was not skipped.

The command, per-test results, 88 attachment records, download manifests, and hashes of the tested working source are preserved in [the final results](audit/ui-audit-final-results.json). The full local Playwright report is at `frontend/test-results/e2e-results.json`, with a retained copy under `data/integration-checks/ui-functional-audit/playwright-full-results.json`. [The collector](audit/ui-audit-collect.py) extracts this compact record from actual test output. The earlier [first-pass record](audit/ui-audit-first-pass.json) retains exploratory native-job evidence, including the deliberately unsuitable GROMACS configuration described above.

```sh
cd frontend
DYNAMOL_BASE_URL=http://127.0.0.1:8765 \
DYNAMOL_AUDIT_COMPLEX_DATASET=4cf5bcedf38847c1 \
npx playwright test
```

Afterward, 103 identified UI test datasets and 21 test jobs were moved into the local audit archive, with no deletion. [The archive manifest](audit/ui-audit-archive-manifest.json) records every original path, archived path, provenance rule, and metadata checksum. The active library was checked afterward: all original 9AX6 inputs, the prepared whole complex, the bounded native-validation system and trajectory, and other non-test datasets remain present; there are no active jobs. The prepared-complex parameter archive remains available in the application.

Representative captured views: [desktop ribbons](audit/ui-audit-desktop-ribbons.png), [surface](audit/ui-audit-surface.png), [polar hydrogens with protein context](audit/ui-audit-polar-hydrogens.png), [real explicit water](audit/ui-audit-explicit-water.png), [preparation completed](audit/ui-audit-preparation-complete.png), and [phone Studio](audit/ui-audit-phone-studio.png).

This is functional coverage of the implemented local interface on macOS Chrome. It is not exhaustive browser/device certification, formal WCAG conformance, a GPU benchmark, or validation of every chemical species, force-field choice, protein conformation, and simulation duration. Native full-screen behavior inside an embedded app browser can differ from Chrome; the app's rejection fallback is checked separately. Long trajectories, production convergence, alternate operating systems, and hardware-specific GPU paths need their own validation.


## Follow-up: pinch, fullscreen, live values and modified residues

The earlier 59-test audit missed the user-reported trackpad pinch, clipped fullscreen toolbar and delayed measurement-value behaviors. Those defects are now corrected and covered by focused regressions:

| Behavior | Fix and evidence |
| --- | --- |
| Trackpad/touch pinch | Captured Ctrl-wheel and WebKit gesture scaling zoom the molecular camera; continuous two-finger touch ratios also preserve center panning. Four canvas-geometry tests cover zoom in/out, slow touch motion, duplicate-event suppression and ordinary wheel behavior. See [pinch results](audit/pinch-zoom-results.json). WebKit gesture dispatch is synthetic; physical Safari/trackpad certification remains outstanding. |
| Fullscreen toolbar | Fullscreen contains the entire app shell, toolbar and dialogs, with reachable actions at desktop and narrow/short viewports. Five tests include entering/exiting fullscreen, Open files, New simulation and missing fullscreen API recovery. |
| Immediate measurements | Complete distance, angle, dihedral and donor–H–acceptor selections show the current original-frame value in the inspector and scene before plotting. Three browser tests exercise values, playback, selection changes, invalid geometry and plot creation; backend tests verify geometry and periodic behavior. |
| Modified-residue inspection | A real 1UA2 browser test checks four TPO records, separate ATP ligand cards, fixed-state disclosures, genuine gaps and modification retention when free heterogen removal is selected. Unsupported-residue blockers have separate backend coverage. Native preparation and parameter validation are recorded separately in the [modified-residue review](modified-residue-review.md). |

The combined first three files passed **12/12 browser tests in 22.9 seconds**, followed by **1/1 modified-residue UI test**. [The combined record](audit/ui-interaction-fixes.json) retains actual test results. The backend suite passed **110/110** after the molecular changes; [its output](audit/ui-modified-backend-tests.txt) is retained. These are follow-up checks, not a claim that the original 59-test run was repeated against every subsequent change.

After the package-origin correction and missing-modified-sequence guard, the complete backend suite passed **131/131** in 12.28 seconds; [final output](audit/final-backend-tests.txt) is retained.

Package launch testing subsequently exposed an origin allowlist tied to fixed development ports. The launcher now supplies one exact, validated localhost origin for its allocated port; middleware and CORS use the same list. Foreign origins and other unconfigured ports remain rejected. `backend/tests/test_packaged_origin.py` covers dynamic-port requests, CORS preflight and invalid configuration. Packaged browser and job-lifecycle validation is recorded under `packaging/validation/`.

Final preparation review also found that PDBFixer substituted ordinary parent amino acids for missing modified sequence residues. The corrected alignment retains exact residue identities; ten regression cases cover SEP/TPO/PTR/HYP/MSE/ALY/unknown identities, submission and worker guards, and unchanged ordinary short-loop support. The real TPO preparation → solvation → MD continuity check was repeated successfully on the final worker.

## Workspaces, recovery and informative analysis

The final additions extend the inventory with automatic scene/camera restoration; named project save/update/open/export/import; library search/rename/archive/trash/restore; named atom selections; preparation/simulation readiness with shared admission limits; live native energy/temperature plots; checkpoint Resume; RMSD/RMSF settings, numeric curves, click-to-seek/select, CSV and JSON analysis records. Fullscreen coverage now includes Projects and its modal, and workspace restoration includes expanded analysis layout.

The [82-case union](audit/release-readiness/ui-coverage-union.json) retains each test's command, source report hash, available build evidence and run history. The first remaining legacy run had five failures: two assertions expected invalid preparation to remain clickable, one isolated fixture lacked its native preparation-job folder, and two exposed a real startup-cancellation classification bug. Corrected readiness assertions and the complete fixture passed on retry. The supervisor now distinguishes requested cancellation from unexpected interruption after the worker exits; both original browser cancellation regressions and four additional native startup checks pass. Earlier failure records remain available.

All listed molecular format families were exercised again. The native whole-9AX6 prepared PDB and parameter archive download was checked without skipping, as was real 1UA2 inspection with four TPO and four ATP records plus unresolved sequence gaps. These checks do not claim that the full 1UA2 structure was prepared or that the whole 9AX6 complex fits the viewer's explicit-solvent size cap. The broader [chemistry support matrix](CHEMISTRY_SUPPORT.md) records molecular models and expected unsupported cases separately.

New views: [structural analysis](audit/release-readiness/structural-analysis-desktop.png), [narrow analysis](audit/release-readiness/structural-analysis-narrow.png), [ready-to-run estimates](audit/release-readiness/simulation-readiness.png), and [restored workspace](audit/release-readiness/workspace-restored.png).


## Follow-up: automatic bottom-panel measurements and zoom buttons

On September 11, the user reported two defects that earlier final-state checks missed. The bottom **Add a measurement / Measure** actions previously enabled picking without actually adding the completed selection to the plot panel. They now open a visible draft with atom-count progress and automatically create and activate its plot when enough atoms are selected. Pending calculations expose progress, cancellation and retry. Selection/type/dataset changes invalidate stale results, and stopping a draft returns to ordinary inspection. Inspector-only selections retain immediate values and explicit **Plot over time**.

The +/− camera animation passed a positive distance to an NGL animation that starts from signed negative camera Z. Its path crossed zero, causing a transient oversized molecule and a blank frame. Button zoom now interpolates positive distances through the same camera control as pinch. Repeated clicks accumulate a target and cancel competing button animations; pinching immediately takes control. Initial camera spacing and resize use the direct distance control as well.

The production frontend passed **30/30 focused Chrome browser tests**, zero retries/skips, in **58.3 seconds**, using a separate API and data directory. Four new measurement cases cover automatic distance/angle/dihedral plots, exact saved values and reload, visible errors/retry, cancel/type-change races, and stopping/Escape. Three new zoom cases inspect the real WebGL image on each animation frame, including rapid clicks and pinch interruption. Existing live measurements, gestures, fullscreen and general control regressions passed in the same command. TypeScript and the Vite production build passed. No chemistry or simulation implementation changed, so the native engine/chemistry suite was not repeated.

[Per-test results and source hashes](audit/measurement-zoom-fixes/browser-results.json), [selection progress](audit/measurement-zoom-fixes/measurement-draft.png), and [three automatically added plots](audit/measurement-zoom-fixes/three-measurements.png) retain the evidence. The [original zoom frames](audit/measurement-zoom-fixes/original-in-frames.json) reproduce a 109→372→130 pixel width and a blank frame; the [fixed frames](audit/measurement-zoom-fixes/in-frames.json) stay within the start/end molecular size with the test's rasterization tolerance. Gesture tests use synthetic browser events; this is not physical-device or Safari certification. The 30 cases are a focused follow-up, not a repeat of every prior release check.
