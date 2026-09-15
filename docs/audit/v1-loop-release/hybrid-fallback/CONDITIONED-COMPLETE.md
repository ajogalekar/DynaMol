# Source-conditioned construction: complete-complex replay

**Supplying the real boundary context is implemented and tested, but none of
the 12 complete-complex refinements passes every gate.** One source-conditioned
torsion/CCD candidate now passes both the modeled-residue reference and the
newly defined loop boundaries after refinement. Its remaining reference
failure is the pre-existing ARG57 outlier whose defining coordinates moved.
No candidate is admitted, and no acceptance policy changes.

## Intervention and preserved controls

The installed native default pipeline, inherited by DynaMol's fallback,
constructs the torsion sampler with default flanking phi/psi and alanine
neighbors. The research adapter now has explicit source-conditioned options.
For 1UA2, these use observed phi43 = −123.943°, psi56 = 62.367° and K/R
neighbors. The source-conditioned sampler and native torsion-aware CCD closer
are separate interventions. FragmentSampler has no equivalent boundary
fields; that branch changes the closer only. Exact API and source inspection
are in [BOUNDARY-API.md](BOUNDARY-API.md).

The [independent native comparison](CONDITIONED-NATIVE.md) completed 15 fixed
attempts. Nine candidates were returned; all three masked controls passed
the local backbone screen. All six 8K5R attempts failed initialization. All
six 1UA2 candidates still exceeded the observed-backbone cap because of
**ILE43 O**, despite small N/CA/C fitting residuals. Native CCD's boundary
probabilities do not preserve that carbonyl orientation.

Every returned 1UA2 proposal was evaluated in the complete archived complex,
without selecting on its local screen. Each was checked twice:

- Existing mobile-flank refinement, retaining coherent native context as the
  starting proposal and the original 1 Å final observed-heavy cap.
- All observed heavy atoms fixed at their archived coordinates. This restores
  source anchors before the same bounded Cartesian refinement and can strain
  the seed's joins; initial failures remain recorded. Restoration itself is
  never treated as a successful repair.

The same original 1UA2 archive, ATP/TPO chemistry, native force field, hydrogen
identities, 1,000-iteration refinement and acceptance thresholds apply. The
plan was frozen before full-complex scores. There were no stronger
construction forces, extra random draws or outcome-driven retries.

## Results

Each cell lists the number of modeled reference outliers followed by affected
observed outliers. An empty observed list does not imply a passing model.

| Native proposal / seed | Mobile flanks: modeled; observed | Fixed observed: modeled; observed |
|---|---|---|
| Fragment + CCD / 2026 | 2; ARG57 | 3; none |
| Fragment + CCD / 2027 | 4; ARG57 | 4; none |
| Torsion + CCD / 2026 | 2; ARG57 | 3; none; stereo also fails |
| Torsion + CCD / 2027 | **0; ARG57** | 1; none |
| Torsion + DirtyCCD / 2026 | 3; ARG57 | 4; ASN56 |
| Torsion + DirtyCCD / 2027 | 2; ARG57 | 2; none |

All 12 final outputs pass the broad geometry, original displacement,
identity, retained-environment and current parameter-integrity checks.
Eleven pass stereo and the saved-file checks. None passes the complete
backbone reference. All six fixed-observed seeds initially fail broad
geometry; Cartesian refinement repairs those gross defects but does not
remove their modeled reference failures. One fixed-observed output also
fails stereochemistry and remains rejected.

The best mobile-flank candidate has zero modeled outliers and valid newly
defined ILE43/ASN56 boundary pairs. Its remaining ARG57 failure is not a newly
invented source defect: ARG57 was already an outlier, but movement of ASN56 C
changes its phi. Fixing every observed atom preserves that existing outside
conformation, yet the matched repair then leaves LYS44 as a modeled outlier.
The [source-boundary distinction](SOURCE-BOUNDARIES.md) remains in force.
Neither result can substitute for the other or be combined into a claimed pass.

All outputs retain the complete 4,825-atom archive and exact atom/state/bond
inventory. Retained heavy atoms outside each explicitly modeled/mobile region
remain fixed. The fresh native systems are identical across all 12 attempts;
construction forces do not enter them. Historical base-XML equivalence remains
unverified, as in the earlier replay; current base files and the archived
ligand/modified-residue bundles are hash-bound. Source, parameter and output
hashes were rechecked after completion.

The 12 validations took 45.74 seconds total, with maximum sampled RSS236.06
MiB, within 180 seconds/4 GiB per attempt. These are replay timings, not full
app preparation latency. All native jobs are complete. No MD/QM, new chemical
parameters, app default, package, video or publication changed.

## Concrete next construction test

The remaining observed motion is localized well enough to test a specific
closure change. Align the observed n-stem and rotate its carbonyl together
with the downstream chain; then move the numerical closure anchor to the
first modeled residue. Native suffix closure can vary modeled torsions
without rotating the preceding observed carbonyl. An intact control verifies
the necessary transforms and first-join preservation; see
[CARBONYL-LOCK-DESIGN.md](CARBONYL-LOCK-DESIGN.md).

That design is not yet a protein result. A fixed virtual anchor per attempt
freezes the first modeled phi; per-proposal anchor updates need a tested
driver. The far-end carbonyl and full environment still require independent
checks. Evaluate a fixed small set of coherent saved proposals before any
larger search. Preserve full original coordinates and all failed outputs.

Refinement must also avoid unnecessary changes to fully observed backbone
torsions outside the repair region. If partial anchor mobility is studied,
freeze the complete set of atoms defining those outside torsions and verify
actual unchanged coordinates; do not simply ignore an affected pre-existing
outlier. Do not loosen the existing geometry/reference or movement thresholds.

Artifacts: [machine-readable complete summary](conditioned-complete-summary.json),
[frozen comparison driver](validate_conditioned_candidates.py), and
[validator](validate_complete_candidate.py).
Cache: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/boundary-conditioned-complete-v1`.
