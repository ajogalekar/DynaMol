# Hybrid fallback: first complete-complex comparison

The acceptance-driven search prototype works, but **no candidate from this
bounded comparison clears all preparation checks**. It remains disconnected
from the app's default preparation path. The useful result is narrower: one
native fragment proposal for the 12-residue 1UA2 gap passes the modeled-residue
backbone reference after complete-complex refinement; its affected joins still
fail. DiSGro is not selected as a release dependency.

## What was implemented

`backend/loop_search.py` orchestrates a predeclared bounded attempt plan. A
closed candidate rejected after refinement reaches the next generator.
Unavailable candidates also advance; cancellation, timeouts, and integrity or
programming exceptions stop the search. Every attempt keeps its own files,
checks, rejection reasons and progress stages. Missing checks cannot become a
pass. Successful decisions must identify the final validated artifact and its
SHA256, verified within that attempt's directory; the search never returns an
unrefined proposal as the accepted model. Validation mappings are immutable.

Twenty focused lifecycle tests pass, and 69 tests pass when combined with the
existing refinement and native-worker suites. The real comparison also
demonstrates the intended transition: each archived closed fragment fails the
backbone reference, and the next scheduled proposal is actually validated.
Both case searches end exhausted with no accepted artifact.

The UI remains one proposed automatic workflow. Length guides the intended
search budget; this work does not establish a reliable six-residue algorithm
boundary or justify exposing algorithm choices to users.

## Frozen proposals and results

The separate [native comparison](NATIVE-MC-COMPARISON.md) tested 20 attempts,
returning 16 coherent heavy-atom candidates, with one local masked-control
pass. The [private DiSGro evaluation](../disgro-evaluation/REPORT.md) retained
40 proposals across two cases; none initially met the observed-anchor cap.

Before full-complex refinement, the root plan fixed these eight attempts:
the archived seed for each case, all four original-stem native 1UA2 proposals,
and one DiSGro representative per case selected from the previous local
screens. This is a feasibility probe, not a population success-rate estimate.
There was no resampling after these final outcomes.

| Case / proposal | Modeled backbone outliers | Other affected outliers | Maximum observed heavy motion | Broad geometry | All checks |
|---|---:|---|---:|---|---|
| 1UA2 archived fragment | 3 | ASN56 | 0 Å | Pass | Fail |
| 1UA2 torsion MC, seed2026 | 2 | ASN56, ARG57 | 0.751 Å | Pass | Fail |
| 1UA2 torsion MC, seed2027 | 2 | ARG57 | 0.739 Å | Pass | Fail |
| 1UA2 fragment MC, seed2026 | 3 | ASN56, ARG57 | 0.709 Å | Pass | Fail |
| 1UA2 fragment MC, seed2027 | 0 | ASN56, ARG57 | 0.757 Å | Pass | Fail |
| 1UA2 DiSGro representative003 | 2 | ARG57 | 0.787 Å | Pass | Fail |
| 8K5R archived fragment | 5 | PRO91 | 0 Å | Pass | Fail |
| 8K5R DiSGro representative005 | 5 | PRO91, GLN181 | 0.759 Å | Fail | Fail |

The 8K5R proposal replaces only the three-residue AKN loop; the complete
archived structure contains four gaps and 17 modeled residues. All those
modeled residues remain checked and refined. Its five outliers therefore do
not mean five outliers in the three-residue target. The geometry failure is
an out-of-range C–N–CA angle. Seven of eight saved models also pass the saved
geometry/stereo check; the remaining one preserves the same rejected geometry.

ARG57 in 1UA2 is a **pre-existing source outlier**, but its phi changes when
ASN56 C moves. ASN56 is different: its phi becomes defined only after the
missing ILE55 is built. The [source-boundary audit](SOURCE-BOUNDARIES.md)
records this distinction and confirms that these selections do not
unnecessarily include a wholly unchanged source-supported outlier. No
threshold was relaxed and no source defect was silently reclassified.

## Complete-complex and parameter checks

The validator binds each proposal to its hashed source and verifies the exact
observed protein-backbone identities and coordinate frame against the archive.
It inserts only the explicit modeled/context coordinates into the complete
archived topology. ATP/TPO for 1UA2 and VQE/TPO for 8K5R remain present.
Unprovided canonical sidechain coordinates are initialized by a proper
residue-frame transform, not asserted to be predicted rotamers. Only affected
hydrogens are rebuilt with the exact archived hydrogen names and parent bonds;
subsequent hydrogen relaxation preserves those identities and protonation
states. Explicit context cannot move crosslinked or metal-coordinating
residues. This does not qualify a new metal model.

All eight checks of exact identity, retained heavy environment, stereo and
the original 1 Å observed-heavy displacement cap pass. The CCTBX reference
covers the modeled residues and affected boundaries. Observed-only peptide
cis/trans basins are checked independently after refinement. PDB round trips
compare atoms/elements, full connectivity, coordinate rounding, geometry,
stereo, reference classifications and observed peptide basins.

The archived ligand and modified-residue parameter bundles pass their existing
checksums. Current protein/water XML files and includes are recorded and frozen
for this replay. Fresh native systems are identical across all attempts within
each case, and temporary construction forces do not enter them. The historical
archive did not hash the base protein/water distribution: the parameter check
establishes current replay integrity, not unverified historical file equivalence.
No new parameters, charges, MD or QM were generated.

All eight validations completed within their 180-second/4-GiB limits:
79.81 seconds total native validation, with maximum sampled RSS234.27 MiB.
Those timings exclude the earlier generator evaluations and do not estimate
full application preparation latency. All recorded input, implementation,
parameter and output hashes are unchanged. There are no remaining native jobs.

## Next step and release boundary

Retain the native fragment route as the leading fallback prototype. Work on
the geometry of the observed joins during candidate construction, including
their carbonyl direction and the torsions defined across the gap. Preserve
the distinction between newly modeled geometry and existing source outliers.
Do not just resample the same malformed seed, increase construction forces,
or relax the movement/reference gates. These failures do not establish that
the missing conformation is unrecoverable.

DiSGro's parser omits heterogens, modified residues and some source identities;
its sampler did not see that omitted environment. Reinsertion and later
validation preserve the full complex but do not retroactively make its search
environment-aware. Its inherited code/data redistribution provenance also
remains unresolved. Keep the current evaluation isolated and private.

A successful candidate still needs a full preparation rerun, the frozen loop
panel, UI/engine checks and exact bundle/install verification. No output here
is admitted to the app, used for dynamics or considered a validated physical
model. The inaugural-video update remains queued after the loop release gates.

Evidence: [machine-readable summary](complete-comparison-summary.json),
[frozen-plan runner](run_complete_comparison.py),
[complete validator](validate_complete_candidate.py), and
[summary verifier](summarize_complete_comparison.py).
Full artifacts: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/hybrid-complete-validation-v1`.
