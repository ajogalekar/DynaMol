# Monomer setup: completed software review

The scoped software review is complete, with no unresolved finding. This is a language-model code and evidence review, not biological assembly validation or expert certification. The feature selects one observed protein coordinate chain; it does not establish a functional biological monomer. No new MD, parameterization, predictions or network retrieval was performed for this task.

The PatAgent preflight recorded six passes and no failures or unresolved gates. Representation and provenance principles reuse the evidence-linked review in `../live-simulation-tracking/scientific-review.md`; this deterministic selection task makes no new scientific modeling claim. `review-results.json` records the final reviewed source and evidence hashes.

## Verified behavior

The final combined backend run passed **31 tests**; the reviewer independently reran all **nine extraction regressions**, which passed in 0.63 seconds. These cover unchanged source data, exact retained physical coordinates and atom correspondence, whole-ligand retention, fresh preparation/solvation state, opt-out behavior, authoritative cross-chain covalent bonds, cross-chain SSBOND records with missing SG atoms, selected versus excluded dangling CONECT endpoints, retained sequence gaps, excluded-chain ligand-link scope, and periodic safety refusal. Final build validation also passed.

All **three browser workflow tests passed with no skips or failures**. Missing-atom and supported-loop controls are visible before Prep, with advanced options closed. Suggested and alternate chain selection produce new datasets, preserve the original, transfer retained measurements, explain excluded selections, and restore the full source on request. The UI clears the old preparation banner, locks Start/engine changes during selection, and prevents a delayed result from replacing a new setup. The chain/biological-assembly distinction is visible before selection.

The final isolated **6N7A** API check preserved original bytes and exact mapped first-frame coordinates. It selected chain B, retained two associated residues, removed the old periodic/preparation/solvation state, and retained both original five-residue internal sequence gaps. This was extraction and sequence inspection, not loop building or ligand parameterization. See `current-structure-extraction.json`, `browser-checks.json`, and `validation.json`.

## Review findings resolved

- Archived LINK/struct_conn records on excluded chains could incorrectly block an unlinked retained ligand with the same component name. The derived selection now records retained connection scope, and later ligand screening respects that scope while keeping original records archived unchanged.
- Covalent source identity is checked at residue level, so missing SG atoms do not erase a cross-chain SSBOND. Unresolved links affect extraction only when their known endpoint is retained. The reviewer reproduced a dangling CONECT that was initially dropped; the final fix records that unresolved endpoint, with separate passing selected/excluded regressions.
- Cartesian-only proximity could omit a molecule across a periodic boundary, and removing a box could split retained molecules. Saved-box checks now detect these cases and refuse extraction requiring unwrapping. An independent abstract geometry control put a ligand 0.6 Å away by minimum image but 9.4 Å in direct coordinates; extraction refused, source hashes stayed unchanged, and no derived dataset was created.
- The retained-bond pass now constructs its atom set once, avoiding avoidable quadratic work.

Original risk reproductions remain explicitly labeled in `review-source-link-risk.json` and `review-dangling-conect.json`; they are not final failures. `review-findings.json` and `review-focused-tests.txt` record their closure.

## Scope limits

Nearby molecule retention uses a documented 5 Å heuristic and covalent connected components, with shared-interface notes. It does not identify proven binding partners or validate the isolated biological state. Retained/excluded inventory entries are residue records; connected multi-residue molecules remain whole. Periodic input requiring unwrapping must be supplied in a suitable unwrapped form; the feature does not silently move atoms. Sequence records preserve missing-region identity, while existing length, terminal-region and modified-residue restrictions remain in force. Making repair choices visible does not validate a built loop or expand force-field coverage.
