# Chrome major-function audit

All 31 grouped checks below passed within the stated scope. The water-padding defect was corrected and rechecked. Screenshots and machine-readable receipts are retained beside this report.

| Area | Fresh observation |
|---|---|
| Fresh scene | Ribbons render; Element is selected by default. |
| Representations | Ribbons, Sticks, Ball & stick and Surface render actual geometry; Element, Sequence and Chain colors were exercised. |
| Hydrogen context | Hidden, Polar only and All change hydrogen visibility; Polar only retains the protein heavy-atom context. |
| Visibility | Protein, ligand, retained Zn ion and explicit-water show/hide controls exercised with nonempty relevant structures; measurement overlay show/hide retains plots. |
| Camera | Zoom in/out settle; Fit molecule, Auto rotate, direct drag rotation and scroll zoom respond. |
| Viewport | Keyboard resize increased 280→360 px; native Chrome fullscreen retains Open files and New simulation and exits correctly. |
| Playback | Play/pause, first/last/previous/next frame, direct midpoint timeline click, 2× speed and no-loop stop at last frame exercised on real trajectories. A measurement-chart click also seeks to the matching midpoint frame. |
| Distance | Atom search MET1:N–CA gives 1.46 Å immediately; Add Measurement automatically creates plot. |
| Direct picking | Two actual canvas clicks picked LEU15:CA and GLN2:CA and automatically created a 5.60 Å distance plot. |
| Angle | MET1:N–CA–C gives 109.82° and auto-commits via Add Measurement. |
| Dihedral | MET1:N–CA–C–GLN2:N gives 158.67° immediately; explicit Plot over time creates its trace. |
| Hydrogen bond | Connected donor/H plus acceptor selection yields distance, angle criterion and 0.0% geometric occupancy for the deliberately distant pair. |
| Residue focus and selections | GLY35 focus and named 7-atom selection save successfully. |
| Structural analysis | RMSD and RMSF compute on the demo; RMSD reference frame equals 0.000 Å; RMSF export has 76 residue rows. |
| Client exports | Distance/angle pane CSV, PNG molecular snapshot, RMSF CSV and project ZIP downloads verified on disk; checksums and file structure recorded. |
| Native file chooser | PDB topology plus genuine XTC selected in macOS chooser; stride 2 and 0.25 ps original interval import 113 atoms and 3 frames at 0, 0.5, 1 ps. |
| SMILES | Caffeine built from SMILES and rendered as 24-atom ligand; ligand hide/show works. |
| Web fetch | 1UBQ and 1UA2 fetched in Simulate; PubChem CID 2244 fetched and rendered. Fresh datasets select Element. |
| Engine first | OpenMM and GROMACS engine choices precede prep; GROMACS displays its native standard-protein limitations. |
| Repair options | Build supported missing loops checkbox is selectable before prep; 1UA2 reports 12-residue internal gaps and current 12-per-gap/96-total cap. |
| Monomer | 1UA2 four chains→chain A with ATP/TPO retained; Restore full structure returns all four chains. |
| Preparation monitor | Real 1UBQ prep progresses 5/7→7/7 stages; view log exposes executed stages; completion loads the prepared structure. |
| Explicit water | Selecting Explicit produces a real visible 17,395-atom box with 16,164 water atoms; completion and water visibility controls work. |
| Padding regression | After fix, 0.7 nm shows inline 1–3 nm error and disables Update box/Start; 1 nm removes error and enables valid actions. |
| Pre-run tracking | Add measurement from Simulate tracks MET1:CA–GLN2:CA before starting; unfinished selection blocks Start. |
| Run UI and handoff | Real 2 ps OpenMM run submitted through Chrome; progress/log accessible; completion automatically opens Explore with 21 saved frames and retained measured trace. |
| Run artifact integrity | 19 strict native-output checks pass on the Chrome-created run, including live/posthoc measurement identity and values. |
| Persistence and reset | Reload preserves completed dataset, frame 21, Element color and plot; New simulation resets name, duration and transient prep/water monitors. |
| Named projects | Saved scene opens; project export is a valid 88-entry ZIP including dataset/parameter/provenance files. |
| Help and console | Quick guide opens/closes. Inspected console logs contain only NGL/Three.js legacy-light deprecation warnings, no application errors. |

## Scope and limits

- No fresh physical two-finger trackpad gesture or right-button drag was synthesized. Scroll zoom and direct rotation were exercised; pinch/right-pan historical regression evidence is separate.
- Native chooser was exercised with PDB+XTC, not once for every extension. All advertised formats have fresh API contract tests.
- Fresh native numerical stability scope is two 10 ps engine runs plus this 2 ps UI smoke run; no convergence or biological-state claim.
- The full browser walkthrough began on the initial audit bundle, then reloaded the final bundle for padding regression and workspace/results checks. This is not represented as a single uniform-build automated UI suite.

## Harness observations

- DOM locator fullscreen activation lacked a trusted user gesture; native Chrome click succeeded and header visibility was visually verified.
- Locator.fill on a range input did not trigger the app seek in this automation. Direct pointer midpoint seek was verified at frame 11/21 (1 ps).
- Choosing Explicit automatically built water; the following Build water locator no longer existed because it had become Update box.
- Inspector preview without Add Measurement requires Plot over time; Add Measurement flow auto-commits at a complete selection. Both paths exercised.

## Final-build water error test

Copied the freshly prepared 6A93 audit fixture into the isolated UI library; 1 nm water request refuses above 100,000 atoms. Final UI correctly says Explicit-water setup failed and retains the actionable size reason. No caps changed.

The 6A93 final-UI resource rejection uses a copied prepped audit dataset, not a second full preparation. The new failed-job header was exercised; fast request-failure operation retention was reviewed in source.
