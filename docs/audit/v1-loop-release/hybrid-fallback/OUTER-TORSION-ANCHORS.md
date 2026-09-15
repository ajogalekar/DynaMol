# Outer-torsion anchors: first complete static passes

Two of eleven frozen **1UA2 12-residue-loop** candidates pass all eight private
archived-complex checks. Both retain the 4,825-atom topology, ATP and TPO. This
is a static construction and saved-artifact milestone, not fresh-source app
preparation, a panel success rate, native-loop accuracy or release qualification.
An [independent saved-artifact review](independent-complete-audit/REPORT.md)
confirms both endpoints and separately confirms the handoff-test endpoint.

## General correction

The previous source-conditioned candidate had no modeled or newly defined
boundary outliers, but moved ASN56 C, which changed the defining geometry of
the pre-existing ARG57 outlier. Exempting that changed residue would have been
incorrect. The new policy instead holds the observed outer torsion-support
atoms exactly fixed during temporary refinement. Peptide connectivity selects
left-flank N/CA and right-flank CA/C; no PDB ID, residue number or outlier label
selects an exception.

All 852 fully observed external phi/psi/omega definitions remain exactly
unchanged in every comparison. Modeled definitions and all actually changed
or newly defined boundary torsions still receive the same CCTBX reference
check. Preceding CA is included because omega can affect the proline class.
Unavailable reference rows remain failures. A replay of the new selection
rule on all 20 previous complete-complex comparisons reproduces every old
selection, so none of their failures was reclassified.

Only the temporary construction System receives additional fixed masses.
The original force constants, 1 Å observed-heavy movement limit, iteration
budget and other acceptance checks are unchanged. Every attempt uses the same
fresh physical System hash. The historical archive did not pin its original
protein/water XML distribution; this proves current replay consistency and
archived ligand/modified-parameter integrity, not historical base-XML identity.

## Frozen results

The plan includes all six returned source-conditioned 1UA2 proposals and all
five mechanically converged shifted-anchor proposals. Failed native closures
remain rejected and are not full-refinement candidates. See the separate
[shifted-anchor comparison](shifted-anchor/REPORT.md), including its failed
masked-control closures. Native proposal and final refinement outcomes must
not be conflated.

| Attempt | Proposal | Maximum observed-heavy movement | Failed final checks |
|---|---|---:|---|
| 01 | Fragment MC, source CCD, seed 2026 | 0.754 Å | Backbone reference |
| 02 | Fragment MC, source CCD, seed 2027 | 0.758 Å | Backbone reference |
| 03 | Torsion MC, source CCD, seed 2026 | 0.754 Å | Backbone reference |
| 04 | Torsion MC, source CCD, seed 2027 | 0.754 Å | None |
| 05 | Torsion MC, source-dirty closure, seed 2026 | 0.740 Å | Backbone reference |
| 06 | Torsion MC, source-dirty closure, seed 2027 | 0.753 Å | Backbone reference |
| 07 | Shifted fragment/source CCD, seed 2026 | 0.709 Å | Geometry, reference, saved geometry |
| 08 | Shifted fragment/source CCD, seed 2027 | 0.602 Å | None |
| 09 | Shifted torsion/source CCD, seed 2027 | 0.785 Å | Geometry, reference, saved geometry |
| 10 | Shifted torsion/source-dirty closure, seed 2026 | 0.744 Å | Geometry, reference, saved geometry |
| 11 | Shifted torsion/source-dirty closure, seed 2027 | 0.592 Å | Backbone reference |

Identity, retained environment, stereochemistry, observed displacement and
parameter-integrity checks pass in all eleven. Geometry and saved-output
checks pass in eight. Two pass the affected backbone reference and therefore
all eight checks. All eleven children exit normally. The complete-refinement
stage takes 45.0 seconds total, 3.26–4.79 seconds per candidate; peak sampled
RSS is 242,499,584 bytes. These timings exclude candidate generation and
initial chemistry preparation and do not predict total UI latency.

The seven new context-integrity tests and twenty existing search-runner tests
pass. Tests cover graph-derived anchors, changed external definitions,
preceding-omega support, missing reference data, ambiguous peptide branches
and isolation of temporary masses. Full evidence is in
[outer-anchor-summary.json](outer-anchor-summary.json), with source and output
hashes. All 52 distinct recorded source/output files were verified unchanged.

## Actual accepted-coordinate handoff

A separate known-outcome integration regression runs the generic search
orchestrator with an archived fast-fragment seed, then two conditioned fragment
MC and two conditioned torsion MC proposals. The first four fail complete
quality checks. The fifth passes; the runner returns the verified final
`refined.npz` coordinate artifact, not the proposal or a result JSON. Fifteen
progress events record generation, validation and completion. The original
inputs stay unchanged. This regression deliberately uses known candidates;
it does not estimate independent prediction performance.

The repeat starts with bitwise-identical heavy-atom coordinates but slightly
different hydrogen coordinates (maximum 0.00108 Å). Its refined coordinates
therefore are not a bitwise replay: the maximum difference is 0.579 Å at a
modeled Asp sidechain oxygen. Its own saved endpoint independently passes all
checks. This demonstrates accepted-artifact handoff, not deterministic native
minimization or equivalence of different final coordinate files.

The independent audit verifies all atom identities, 4,886 indexed bonds,
ATP/TPO heavy coordinates, 14 affected reference residues and outside ARG57
definitions across source, archive, refined arrays and reloaded PDB. The fresh
physical System has 4,825 positive particle masses and no construction forces.
Reloading the PDB changes harmonic-bond term order in serialized XML, but all
bond parameters match as exact atom-indexed multisets and every other System
field/force term is canonically identical. No physical parameter discrepancy
was found.

Driver: [run_accepted_handoff.py](run_accepted_handoff.py). Full output:
`/Users/ashujo/.cache/dynamol-research/v1-loop-release/accepted-coordinate-handoff-v1/`.

## Next integration boundary

The [independent integration review](integration-review/REVIEW.md) identifies
the seam after creation of the immutable, sequence-restored scaffold. Retry
must wrap complete candidate assembly, state-preserving H reconstruction,
refinement, reference and saved-file checks. The initial fragment must pass
the same gates as every fallback. Chemical parameter bundles should be fixed
once and reused; only a fully bound accepted output bundle may be published.

Fresh preparation needs stronger complete-source and multi-gap binding than
this intentionally archived, partial-loop replay. Coherent native context,
both coordinate baselines, cancellation during native validation and all
publication sidecars must be explicit. First run this as an internal opt-in
with fresh-source 1UA2, then the frozen multi-length panel including 8K5R and
the second 12-residue case. Keep one automatic repair control in the UI;
length guides the search budget, and quality failure can trigger fallback at
any supported length. App enablement, bundle checks and video update follow
those gates.
