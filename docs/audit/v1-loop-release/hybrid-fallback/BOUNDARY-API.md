# Native boundary conditioning and stem preservation

**There is a concrete, bounded intervention available: supply the actual external boundary torsions and neighbor identities to `PhiPsiSampler`, and compare `CCDCloser` with the existing `DirtyCCDCloser`. This conditions construction on the source structure; it does not hold every observed stem atom fixed or establish a repaired loop as acceptable.** No protein search was performed for this audit.

## What the current path does

[native_mc_adapter.py](/Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/hybrid-fallback/native_mc_adapter.py:125) uses the sampler's default boundary values: n-stem phi approximately −60°, c-stem psi approximately −45°, and alanine neighbors. It then uses `DirtyCCDCloser` and reduced/cb_packing/clash scoring. The installed native `FillLoopsByMonteCarlo` does the same at `_closegaps.py:1079–1090`. The app's `construct_loops` calls that native function, so this behavior is inherited by the app fallback; it is not unique to the research adapter. The app calls MC only if gaps remain after its earlier construction stages.

The exact release source distinguishes initialization from later proposals. `PhiPsiSampler.Initialize` draws unconditional phi/psi pairs using sequence-dependent histograms. During `ProposeStep`, the first residue's psi is drawn conditional on `n_stem_phi`, and the last residue's phi conditional on `c_stem_psi`. The supplied neighbors also determine histogram identities, including during initialization. Merely supplying the constructor fields does not make initial closure conditioned on those angles; the constrained closer supplies that additional step. [Official sampler documentation](https://openstructure.org/promod3/3.5/modelling/monte_carlo/), [release 3.6 sampler source](https://git.scicore.unibas.ch/schwede/ProMod3/-/blob/3.6.0/modelling/src/monte_carlo_sampler.cc).

## Exact available calls and identities

These calls were successfully constructed with the installed runtime, without drawing or proposing any conformations:

```python
sampler = modelling.PhiPsiSampler(
    gap.full_seq, torsions,
    n_stem_phi=source_phi_radians,
    c_stem_psi=source_psi_radians,
    prev_aa=previous.one_letter_code,
    next_aa=following.one_letter_code,
    seed=seed,
)
closer = modelling.CCDCloser(
    gap.before, gap.after, gap.full_seq, torsions, seed
)
```

For the sampler, arguments after the first two must be keyword arguments. `CCDCloser` above uses the single `TorsionSampler` object overload. Internally it uses 1,000 CCD iterations and a 0.1 Å N/CA/C endpoint RMSD criterion, as does `DirtyCCDCloser`; these are native closure settings, not the complete-model acceptance thresholds. [Release Python bindings](https://git.scicore.unibas.ch/schwede/ProMod3/-/blob/3.6.0/modelling/pymod/export_monte_carlo.cc), [release closer source](https://git.scicore.unibas.ch/schwede/ProMod3/-/blob/3.6.0/modelling/src/monte_carlo_closer.cc).

| Case | Original missing residues | Original stems → native model numbering | Full native span sequence | External neighbors | Source n-stem phi / c-stem psi (degrees, native calculation) |
|---|---|---|---|---|---:|
| 1UA2 | A:44–55 | ILE43/ASN56 → 31/44 | IKLGHRSEAKDGIN | LYS42 / ARG57, K/R | −123.942867 / 62.367154 |
| 8K5R | A:177–179 | LEU176/SER180 → 172/176 | LAKNS | SER175 / GLN181, S/Q | −67.999123 / −52.704670 |

Both external neighbor links passed native `ost.mol.InSequence`. Exact source identities came from the frozen alignment maps; native residue numbers must not be mistaken for original PDB residue IDs. The small difference from the double-precision 1UA2 values in [SOURCE-BOUNDARIES.md](/Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/hybrid-fallback/SOURCE-BOUNDARIES.md) is native coordinate/geometry precision, not a different residue selection. For production-quality binding, calculate the fields from exact source atom identities and preserve their provenance.

`CCDCloser` obtains these external angles and residue classes itself from its stem handles. If the neighbor is absent, not in sequence, or lacks its defining atom, it falls back to its default angle; noncanonical neighbor classes can fall back to alanine. Such fallback must be explicit in a research attempt record. Neither target examined here required it.

## What closure holds, and what it can change

| Native component | Available boundary behavior | Limit relevant to this task |
|---|---|---|
| `DirtyCCDCloser` | Aligns the first N/CA/C frame and minimizes last N/CA/C RMSD with phi/psi rotations. | First psi and last phi are free closure variables; no torsion-probability checks. |
| `CCDCloser` | Uses the same coordinate target with a torsion-probability Metropolis check, including external source phi/psi at the boundaries. | Still rotates first psi and last phi. It does not impose exact source O coordinates or a hard CCTBX pass. |
| `KICCloser` | Uses three strictly internal pivots and independent N/CA/C stem frame fits; chooses its first returned KIC solution. | No explicit oxygen target or reference check; at least three internal residues are needed. Direct `KIC.Close` exposes all solutions for a prescribed pivot set. |
| `SoftSampler` | Small-angle proposals with the same external conditional fields. | Smaller proposals do not freeze observed geometry. |
| `FragmentSampler` | Can restrict insertion positions through `sampling_start_index`, selected fraggers, and `init_bb_list`. | Insertions transform the downstream tail; avoiding a stem as an insertion site does not make its Cartesian position fixed. The subsequent closer can still alter boundaries. |

The CCD target contains **only N, CA and C**. Before CCD, `FitCStem` explicitly adjusts the final residue's N–CA/CA–C lengths and N–CA–C angle to those of the target. Both CCD modes then reconstruct c-stem O using the following observed residue if it is in sequence. Thus the algorithm neither preserves all incoming stem geometry unchanged nor matches original O explicitly. N-stem N/CA/C are rigidly fitted, so differences in the incoming versus observed triangle can also leave small residuals. [Official closure documentation](https://openstructure.org/promod3/3.5/modelling/loop_closing/), [release CCD source](https://git.scicore.unibas.ch/schwede/ProMod3/-/blob/3.6.0/modelling/src/ccd.cc).

Source hooks for a future exact-geometry implementation are `SetTargetPositions_`, `FitCStem`, the first-residue branch of `GetSimpleCloseRotations_`, the endpoint iterations in `GetConstrainedCloseRotations_`, and final `ReconstructCStemOxygen`. There is no exposed constructor flag to freeze O or exclude first psi from CCD. Editing these hooks would be a distinct algorithmic intervention requiring controls, not a sampler configuration change.

The earlier [KIC endpoint audit](/Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/native-sampling/KIC-ENDPOINT-AUDIT.md) already verified actual fitted endpoints against returned coordinates. It did not reveal a mapping or transform bug. For KIC, pre-aligning the n-stem carbonyl with a native psi rotation can persist through internal-pivot closure, but incoming stem geometry must already be coherent and matching the complete observed frame must still be checked. This is not a guarantee of loop or complete-complex acceptance.

An intact five-residue fixed-angle control confirms an important implementation detail: `SetPsiTorsion(0, old + 0.1)` defaults to `sequential=False`, moving first-residue N by 0.140863 Å while its O stays put. With `sequential=True`, N/CA/C stay put and O moves 0.107742 Å. CCD applies accumulated rotations with `sequential=True`. A claim that default sampler proposals directly rotate O would therefore be imprecise; the combined proposal/frame-fit/closure process changes its orientation. Last-residue `SetPhiTorsion` in the control moved C and O by 0.149057 and 0.233146 Å. None preserves all four observed stem atoms.

Holding a complete observed carbonyl also restricts the newly constructed peptide direction: the new peptide N must have compatible local geometry with the observed O. A freely sampled first psi cannot independently guarantee this. The control report includes an explicitly schematic planar direction derived from N–CA–C–O; it is not a measured missing-atom conformation or an accepted structural reconstruction.

## Scoring and custom-hook limits

`LinearScorer` accepts named native backbone scores and an explicit environment. Native `TorsionScorer` / `LoadTorsionScorer` exist, but the current MC score subset omits that term. A torsion score is not an exact oxygen target, the independent CCTBX reference, or a full-complex force field. Changing weights would confound the proposed default-versus-conditioned comparison and is not recommended for this next test.

Installed abstract native bases are named `MonteCarloSampler`, `MonteCarloCloser` and `MonteCarloScorer`; the documentation's `SamplerBase`/`CloserBase`/`ScorerBase` names are absent. The native bases are not directly instantiable, and their Boost.Python exports do not provide Python virtual-override wrappers. Do not assume a Python subclass can be passed to native `LoopCandidates.FillFromMonteCarloSampler`. A separate installed Python `SampleMonteCarlo` helper uses duck typing, but its source has local-variable reassignment at acceptance/return instead of reliable caller output mutation; it was not executed or validated here and is not a recommended substitute. An explicit research driver calling native proposal/closure/scoring methods would need its own acceptance and state-update controls.

## Bounded next comparison and limits

The clean next comparison is the parent-owned paired default versus source-conditioned sampler/closer run, with frozen original input, unchanged seeds, step counts, scoring weights, context cap, and complete-model gates. Record exact boundary fields and native neighbor availability before generation. Preserve all coherent context moves, raw closure results and complete N/CA/C/O residuals. Reject or report failures normally; do not restore observed coordinates afterward to manufacture a passing stem.

This comparison changes boundary treatment without changing chemistry or force-field parameters. It remains a hypothesis until candidate and full-complex checks pass. No evidence here validates a 1–6 versus 7–12 length threshold, guarantees either target can be repaired, or establishes robustness over the release panel.

## Evidence and provenance

All new evidence is in `/Users/ashujo/.cache/dynamol-research/v1-loop-release/boundary-api-v1`:

- `introspection_control.py` and `introspection-control.json`: native constructors, exact source mapping, neighbor checks, and three deterministic intact-geometry operations. Native activity took 0.056 seconds with 171,966,464 bytes peak RSS; two-thread environment. No initialization, random draw, proposal, MC, closure, scoring, refinement, MD or QM ran.
- `fetch.json`: official release 3.6.0 source URLs and SHA256 hashes. Key files: `ccd.cc` `e6b2714c5d3fa8afe9fab1db541964ca3f39e60f0c2a3fadb3e6858da50b5ecc`; `monte_carlo_sampler.cc` `81e36f295f6f1147148e00318a9d634767e50a4280a9701200f8d641524340a3`; `monte_carlo_closer.cc` `e84d5790e71337a3e5b1f3594947a1932214ec6a52bf681314a1dd0bdb430e79`.
- `local-provenance.json`: snapshots/hashes of installed headers and native Python integration. Installed `_closegaps.py` SHA256 `0333c158e09e6bd2a11abcb0ceba8b0995642f1134686ccca6662864ff6ae5e6`.
- Introspection output SHA256 `9cbc23bc8a53d5480be8291a82c254bd34f5611d08648991b1bc7283adcea3de`; script SHA256 `2d7a2a617f3dea09b0df4aea28c76c90542598e7978a5fc93b2e59ae1619dfce`.

Runtime: `/Users/ashujo/Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2`, ProMod3 3.6.0 / OST 2.11.1. Existing audit coordinates, input artifacts, backend, shared status and thresholds were not modified.
