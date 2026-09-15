# Final DynaMol release audit · 12 September 2026

**Source verdict: suitable for a v1 release candidate within the documented support limits.** No unresolved functional blocker was found in the exercised major workflows. The fresh local macOS DMG candidate also passed isolated runtime, application and Chrome smoke acceptance from the copy made from that exact DMG. Publication is separate from this audit.

This is an audit of application behavior and representative numerical preparation/execution. It is not scientific validation of predicted loop conformations, protonation populations, metal coordination, thermodynamic properties, or convergence. The reviewer is an AI assistant using recorded tests and an automated evidence checklist, not an independent professional certification.

## Fresh evidence

| Area | Result | Receipt |
|---|---|---|
| Backend regression | **601 passed**, zero failures/skips in the final full run | [Regression](regression/README.md) |
| File/API contracts | **53 passed**: all advertised import/trajectory extensions, live PDB/PubChem fetches, malformed inputs, numerical measurement controls, workspaces and repeated reads | [Verification](regression/verification.json) |
| Chemistry matrix | **44 passed**: 28 protein templates/states/termini/disulfides, six ions, ten organic ligands/cofactors | [Chemistry](chemistry/REPORT.md) |
| Chemistry boundaries | **Six checks passed**, including explicit rejection/review paths and retained-ion solvation | [Chemistry summary](chemistry/summary.json) |
| Real complex preparation | **6DBK, 6A93 observed chain A, and 1UA2 observed chain A passed** with retained identities and native ligand parameter parity | [Complex details](chemistry/summary.json) |
| Fresh native runs | OpenMM and GROMACS each completed **10 ps**, 51 frames and two restarts; **61 numerical/artifact/API checks passed** | [Engines](engines/REPORT.md) |
| Chrome UI | **31 grouped major-function checks**, including direct picking, representations, water, file chooser, monomers, preparation monitor, exports, project persistence and completion handoff | [UI walkthrough](ui/REPORT.md) |
| Source-material coverage | **22 checks passed** binding existing source archives, notices, native dependencies and actual packaged source to this exact build | [Binding receipt](package/source-materials-binding.json) |
| Final macOS DMG | Build **589ff2dddea916f5**: container/copy integrity, bundled-runtime self-tests, real workers in both engines, startup/restart persistence and six packaged Chrome smoke checks passed | [Package receipts](package/acceptance-status.json) |
| Chrome-created simulation | Additional **2 ps** explicit-water OpenMM run, 21 frames, **19 artifact checks passed**; automatically opened in Explore with its plot and survived reload | [Native validation](ui/native-ui-run-validation.json) |

Counts above describe different levels of evidence and must not be summed into a single number of independent scientific cases. The UI walkthrough spans the initial and final frontend builds; the final-build fixes and affected persistence paths were rechecked. Source hashes and exact changes are in [source-final.json](source-final.json) and [source-changes.patch](source-changes.patch).

Five earlier **100 ps** trajectories also pass the current validator (**95 checks**). Those kinase/GPCR/channel/protease/phosphatase results are explicitly **historical, read-only revalidation**, not five new runs. The two fresh 10 ps tests and 2 ps UI run demonstrate execution and continuity; they do not establish equilibration or stability on production timescales.

## Concrete chemistry outcomes

- **6DBK:** PTR and G5D retained; preparation and a 50,033-atom explicit-water preview completed.
- **6A93 observed chain A:** charged 8NU, cholesterol, three 1PE molecules and Zn retained. Explicit repair built nine missing heavy atoms in two 1PE molecules. Preparation passed; its 1 nm water box exceeded the default 100,000-atom cap and was correctly refused. The final Chrome UI reproduced that refusal with the correct water-specific error heading. No partner was silently removed and no cap was relaxed.
- **1UA2 observed chain A:** TPO and ATP retained; its 12-residue internal loop was built with the documented ProMod3 workflow. The 58,392-atom water preview completed. Modeled loop coordinates remain provisional, with no native-conformation accuracy claim.
- **1,013 protein/PTM stereocenters** passed the applied screens. Original inputs remained unchanged; the chemistry audit rechecked **1,171 artifact hashes**.

## Fixes and corrections

1. Water-box padding outside **1–3 nm** now receives an inline message and cannot be submitted. Water request errors retain the correct operation heading instead of appearing as protein-prep errors.
2. New simulation provenance reads the shipped application version instead of hardcoding `DynaMol 0.1.0`. Resuming a historical run preserves its original provenance.
3. Added regression coverage for default Element coloring and preservation of saved user choices. The requested Element defaults are present in fresh loads and unrelated result handoffs; Ribbons remains the initial geometry.
4. Corrected current documentation that still described the older six-residue loop limit, universal rejection of incomplete ligands, and a narrower named-ion list. Historical reports were preserved.

The first audit-harness failures were retained: one API test initially expected saved state rather than a save receipt, and one TestClient used a forbidden non-loopback hostname. Correcting those test assumptions required no application change. Chrome fullscreen needed a native user gesture; a native click verified the header is not clipped. The final frontend build and `git diff --check` passed.

## Release boundaries

- The default profile caps structures/solvation at **100,000 atoms**. The complete retained 6A93 case cannot finish the default water setup at its minimum padding; users must choose a smaller appropriate system or an explicitly configured larger profile outside this test.
- Automatic prepared protein–ligand/PTM workflows are **OpenMM** workflows. GROMACS retains its separate standard-protein topology path; cross-engine complex transfer is not implied.
- Unsupported covalent ligands/cofactors, radicals, ambiguous graphs, unregistered amino acids such as MSE/ALY, and specialized metal chemistry still require additional modeling. Named Zn support is a fixed nonbonded ion model, not an inferred coordination/oxidation-state model.
- Internal-loop modeling is bounded at 12 residues per gap and 96 total; terminal sequence extensions are not built. Fixed protonation/PTM choices are recorded approximations.
- Fresh UI testing was on **Chrome/macOS**. No fresh physical two-finger trackpad gesture or right-button drag was synthesized. Scroll zoom, direct rotation, buttons and native fullscreen were tested. This is not a Windows/Linux or every-browser certification.
- The existing GitHub DMG predates the current Element-default and audit fixes. This source audit must not be presented as acceptance of that old download. A new candidate has now passed isolated bundled-runtime acceptance; use that candidate rather than the old download when replacing release assets. The bundle remains without Developer ID notarization, as previously chosen; this audit does not establish clean-machine Gatekeeper acceptance.

## Evidence-linked review

[Preflight](review.preflight.json) passed six gates. [Postflight](review.postflight.json) retains seven passes and an **uncertainty block**: no independent replica confidence intervals or predictive scientific validation were established. This remains a block on promoting molecular-model accuracy/convergence claims; it is not relabeled as a full scientific pass. The software-release verdict above is based on the specific deterministic, native and UI evidence, with those claims excluded.

The cached Practical Cheminformatics corpus was searched for structure representation and reproducibility guidance; [retrieval receipt](corpus-search.json). Relevant public discussions include [The Trouble With Tautomers](https://patwalters.github.io/The-Trouble-With-Tautomers/) and [We Need Better Benchmarks for Machine Learning in Drug Discovery](https://practicalcheminformatics.blogspot.com/2023/08/we-need-better-benchmarks-for-machine.html). Their role here is methodological context for preserving representation assumptions and test scope, not certification or endorsement of DynaMol.

Independent read-only review checked the aggregate against the sibling receipts and corrected the chain-subset labels for 6A93/1UA2. No other count or scope mismatch was identified.

## Accepted download candidate

- Build ID: `589ff2dddea916f5`
- Artifact: `build/releases/final-audit-2026-09-12/DynaMol-0.1.1-macos-arm64-candidate.dmg`
- Size: 1,692,901,364 bytes
- SHA-256: `29adffe128d249f7040abcf05c65a6206a70da2d59c9375393b6c0d63b196e7a`
- Tested by mounting read-only, copying into an isolated Applications folder, and running that copied app with bundled Python and a fresh private home. Both native engines completed application-worker runs. OpenMM, GROMACS, AmberTools and ProMod3 dependency checks passed with restricted PATH and loaded-library provenance.
- The copied app opened in Chrome with the real bundled trajectory, Element color default and both engines marked Installed. [Packaged UI receipt](package/chrome-smoke.json).
- No GitHub asset was uploaded or replaced during this audit; the user’s Applications directory was unchanged.

The final source-material review passed 22 checks for this exact candidate, including its three engine archives and 1,052 manifest-matched packaged source files. It records the working-tree snapshot separately from the baseline Git commit. Existing source-material archives remain applicable. No concrete source/notice coverage blocker was found.
