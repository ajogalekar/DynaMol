# DynaMol missing-loop and preparation audit

**Final result: 10 of 15 structures prepared successfully; 5 remain blocked by explicit chemistry or loop-size limits.** The 10 completed preparations rebuilt 16 internal gaps comprising 90 residues. All 20 independent checks passed for each completed result. No unexpected final workflow failure remains in this fixed panel.

The original 8K5R failure is resolved in the final run: all four sequence-supported gaps were built, the complex was retained, and the saved model passed the geometry and force-field checks. This supports the revised preparation workflow within its stated limits; it does not establish reliable preparation for every PDB entry or the native conformations of reconstructed loops.

## Final 15-entry panel

The entries and chosen chains were frozen before this panel's native preparation runs; the previously reported 8K5R failure was deliberately included. Sequence strings below are the actual missing internal amino acids, with their lengths in parentheses. Original source files remain unchanged.

| PDB / author chain | Target | Missing internal sequence | Chemistry retained or requiring attention | Final outcome | Walltime (s) |
|---|---|---|---|---|---:|
| [8K5R A](https://www.rcsb.org/structure/8K5R) | Human CDK9 | GQGTFG (6); YNRC (4); AKN (3); LYEK (4) | KB-0742 (VQE); phosphothreonine (TPO) | [Prepared](runs/root-final/8K5R/result.json) | 28.07 |
| [1UA2 A](https://www.rcsb.org/structure/1UA2) | Human CDK7 | KLGHRSEAKDGI (12) | ATP; phosphothreonine (TPO) | [Prepared](runs/root-final/1UA2/result.json) | 41.73 |
| [1HCK A](https://www.rcsb.org/structure/1HCK) | Human CDK2 | LDTE (4) | ATP; Mg²⁺ | [Prepared](runs/review-final/1HCK/result.json) | 31.45 |
| [4HJO A](https://www.rcsb.org/structure/4HJO) | Human EGFR | EGEK (4) | Erlotinib (AQ4) | [Prepared](runs/review-final/4HJO/result.json) | 30.80 |
| [1M17 A](https://www.rcsb.org/structure/1M17) | Human EGFR | LPSPTDSNFYRA (12) | Erlotinib (AQ4) | [Prepared](runs/root-final/1M17/result.json) | 59.87 |
| [2HYY C](https://www.rcsb.org/structure/2HYY) | Human ABL1 | D (1); TGD (3) | Imatinib (STI) | [Prepared](runs/root-final/2HYY/result.json) | 77.18 |
| [3HEG A](https://www.rcsb.org/structure/3HEG) | Human p38α | GSGA (4); LARHTDDEMTG (11) | Sorafenib (BAX) | [Prepared](runs/root-final/3HEG/result.json) | 48.62 |
| [6A93 B](https://www.rcsb.org/structure/6A93) | 5-HT2A receptor fusion | IHHSRFN (7); SGS (3) | Risperidone; cholesterol; covalent PLM | [Blocked: covalent model](runs/identities-final/6A93/result.json) | 2.76 |
| [2CG9 X → C](https://www.rcsb.org/structure/2CG9) | Yeast Sba1 cochaperone | PHVGDEN (7); QH (2) | Protein only | [Prepared](runs/identities-final/2CG9/result.json) | 54.40 |
| [2BEL C](https://www.rcsb.org/structure/2BEL) | Human 11β-HSD1 | IVH (3) | NADP; carbenoxolone (CBO); Cl⁻ | [Blocked: reference conflict](runs/review-final/2BEL/result.json) | 3.11 |
| [1FPU B](https://www.rcsb.org/structure/1FPU) | Mouse Abl | LMTGD (5) | Abl inhibitor (PRC) | [Prepared](runs/review-final/1FPU/result.json) | 23.39 |
| [4HHY B](https://www.rcsb.org/structure/4HHY) | Human PARP1 | GGSDDSSK (8) | PARP inhibitor (15R); sulfate | [Prepared](runs/review-final/4HHY/result.json) | 122.95 |
| [2REN A](https://www.rcsb.org/structure/2REN) | Human renin | RLY (3); SENSQ (5); QESYSSKKL (9) | Covalent NAG glycan | [Blocked: covalent model](runs/identities-final/2REN/result.json) | 1.08 |
| [1RNE A](https://www.rcsb.org/structure/1RNE) | Human renin | SENSQS (6) | Renin inhibitor (C60); covalent NAG | [Blocked: covalent model](runs/identities-final/1RNE/result.json) | 1.26 |
| [2ITY A](https://www.rcsb.org/structure/2ITY) | Human EGFR | EKEYHAEGGK (10); PSPTDSNFYRALM (13) | Gefitinib (IRE) | [Blocked: 13-residue gap](runs/identities-final/2ITY/result.json) | 1.27 |

Walltimes are observed full workflow durations with up to three concurrent isolated workers, two CPU threads each. They include import, inspection, parameter assignment, construction and output checks under shared load; blocked cases stopped at eligibility. No MD or water-box generation was performed. All ten final successful cases used the fixed-heavy-atom refinement path. Fixed-heavy here refers to the loop-refinement stage; separately recorded optional sidechain preparation can adjust protein rotamers beforehand.

## v1.0.0 recheck on the released source

Re-run on 14 September 2026 (local time) against the v1.0.0 source, with the
same frozen inputs, driver and settings (seed 42, three concurrent isolated
workers, two CPU threads each). The nine metal-free repair candidates ran
through the real preparation worker; 1HCK stays deferred as a metal case and
the five chemistry/size refusals were not rerun. Results and job-artifact
hashes are in [runs/v1-recheck-20260914](runs/v1-recheck-20260914/); per-run
`workspace/` copies are kept locally and identified by the recorded hashes.

| PDB | Gaps (residues) | Outcome | Refinement path | Walltime (s) | Modeled-residue Ramachandran outliers |
|---|---|---|---|---:|---|
| 1UA2 | 12 | Prepared, 20/20 checks | fixed-heavy | 44.1 | 3 of 12 (LYS44, GLY46, ARG48) |
| 1M17 | 12 | Prepared, 20/20 checks | fixed-heavy | 63.4 | 1 of 12 (PRO966) |
| 8K5R | 6, 4, 3, 4 | Prepared, 20/20 checks | restrained flanks | 33.8 | 5 of 17 |
| 4HHY | 8 | Prepared, 20/20 checks | fixed-heavy | 115.5 | none |
| 1FPU | 5 | Prepared, 20/20 checks | fixed-heavy | 24.6 | none |
| 4HJO | 4 | Prepared, 20/20 checks | fixed-heavy | 33.0 | none |
| 2HYY | 1, 3 | Prepared, 20/20 checks | fixed-heavy | 75.1 | 1 of 4 (ASP276) |
| 2CG9 | 7, 2 | Rejected: modeled heavy-atom bond outside 1.06–1.90 Å after refinement | — | 106.0 | — |
| 3HEG | 4, 11 | Rejected: non-finite coordinates during loop refinement | — | 49.9 | — |

Both rejections are refusals, not accepted models. 3HEG originally surfaced
OpenMM's raw "Particle coordinate is NaN" text as the job error; the worker
now reports it as an explicit loop rejection with unchanged behaviour
([rerun](runs/v1-recheck-20260914-nan-message/3HEG/result.json)). 2CG9 and
3HEG had prepared on the 13 September source, but the
[saved-output audit](../v1-loop-release/SAVED-BASELINE.md) had already found
backbone outliers in both of those models, so the earlier passes are not
evidence of a regression in model quality.

The Ramachandran column uses
[score_ramachandran_recheck.py](score_ramachandran_recheck.py) with the native
CCTBX 2025.11 tables from the separate research runtime, scoring the requested
modeled residues of each accepted `prepared.pdb`
([results](runs/v1-recheck-20260914/ramachandran-recheck.json)). Boundary
residues bonded to a modeled residue are also recorded there; the script does
not separate newly defined boundary torsions from pre-existing deposited
outliers, so the boundary counts are informational only. Ramachandran
classification is a static reference check, not evidence of the native loop.

This recheck confirms the panel's software behaviour on the released source
and quantifies the known conformational limitation stated in
[LOOP_REPAIR.md](../../LOOP_REPAIR.md). It does not change the earlier
verdicts or artifacts.

## Remaining blockers

- **2BEL:** CBO has missing heavy atoms and a stereochemistry conflict in the CCD reference coordinates. The supplied ligand was neither inverted nor removed to force a pass.
- **6A93:** the selected chain retains a covalent PLM attachment requiring a specialized force-field model.
- **2REN and 1RNE:** covalent NAG glycosylation requires a specialized force-field model.
- **2ITY:** one internal gap contains 13 residues, exceeding the 12-residue per-gap limit. The larger-gap refusal is intentional.

The glycan and 13-residue boundaries were anticipated at panel selection. The palmitate and CBO blockers became clear during inspection. All five remain in the 15-entry denominator; none counts as a successful preparation.

## General fixes exercised

- Bounded structural-fragment and sampling fallbacks, followed by validation against original anchors and staged local refinement. The final acceptance gates remain active.
- Source label/author chain mapping, exact ligand connection identities, insertion-code preservation, and correct handling of connections to absent crystallographic images.
- Correction of inverted newly built standard sidechain branches without moving observed atoms, followed by a complete-environment gross-contact check for corrected atoms.
- A declared sulfate reference state, correct handling of ligand preparation with zero newly added hydrogens, and native-validated 1–4 handling for molecules without applicable torsion pairs.

## Validation and provenance

All 15 final runs used the same frozen backend and benchmark driver, with matching start/end hashes and unchanged source inputs. The 10 completed preparations each passed 20 independent checks covering gap preservation, identity and bond preservation, fixed-atom/refinement constraints, expanded and serialized geometry, ligand/ion retention, native ligand conversion, complete force-field reconstruction, and finite static energy and forces. The full backend regression suite passed 647 tests; the frontend production build also passed.

See the [independent final verification](root-final-verification.json), [source checkpoint](final-source-checkpoint.json), [machine-readable summary](summary.json), and [backend test log](backend-final-regression.log). Each table outcome links to its complete result and retained workspace artifacts.

## Initial failures and additional seed checks

The original inventory and all common-runner development attempts remain listed in summary.json and in their original run folders; intermediate successes with changing source files are not substituted for final evidence. Preserved failures include 8K5R fragment/anchor failures, 2CG9 peptide-omega rejection, THR sidechain inversions in 2HYY/1FPU, source-identity/selection failures, and the successive sulfate preparation failures in 4HHY.

An earlier 8K5R seed 2026 attempt failed the 1 Å flank-motion cap at 1.46 Å; that failure remains recorded. Additional final-candidate checks subsequently passed at seed 42 (29.37 s, explicitly recorded local flank movement up to 0.781 Å) and seed 2026 (26.34 s, fixed-heavy path). The final 15-case 8K5R run used the fixed-heavy path. These additional checks exercise both refinement paths, not universal or cross-platform reproducibility. The restrained-flank fallback is a provisional local remodeling step with protected residues and retained molecules fixed, an unchanged 1 Å acceptance cap, and geometry/chirality checks.

See [final 8K5R repeat artifacts](8k5r-full-prep-final/summary.json), [independent flank verification](review-8k5r-flank-verification.json), and the additional seed history in summary.json.

## Interpretation limits

This is a convenience panel of host/model-organism structures, enriched for kinases and chosen for known missing sequence and ligand variety. It is not a random sample of the PDB, and 10/15 is not an estimate of general preparation success. Only one observed chain was selected per entry; that choice is not a biological-assembly assignment. Missing terminal extensions remain outside this loop workflow.

Rebuilt coordinates are provisional starting models. Passing geometry and force-field checks does not validate native loop conformations, conformational accuracy, thermodynamic convergence or MD stability. Those claims require separate evidence.
