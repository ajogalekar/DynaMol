# Carbonyl-aware candidate ranking: bounded comparison

The new backbone baseline identifies a refinement problem that ranking alone
does not solve. This experiment addresses one narrower construction defect:
the old observed-stem fit omitted carbonyl oxygen. A candidate could fit N/CA/C
while placing O more than 2 Å from its observed position.

The worker now provides an opt-in research ranking over N/CA/C/O, reports the
maximum anchor displacement and comparison count, copies native candidates to
avoid dependence on the original collection's lifetime, and enforces the
40-candidate retained budget inside each collection. The old loop checked the
budget only between collections and could retain more than 40 candidates.
The normal application request does **not** enable the new ranking.

The worker and adapter suites pass **42 tests**. Controls distinguish an exact
N/CA/C fit with a displaced O from a coherent backbone, reject nonfinite
carbonyl coordinates under the opt-in, enforce the budget, and verify the
experimental ranking remains disabled in normal application input. These are
software checks, not scientific qualification.

## Native replay

Two fresh native construction runs used the original saved 8K5R and 1UA2
requests. Inputs and source code were hashed before and after; the workers
received private copied inputs and progress files. Every requested missing
heavy-atom identity was preserved, and transplant left all observed scaffold
coordinates bitwise unchanged. Reference scoring included the newly defined
torsions of observed loop boundaries.

| Case | Native construction time | Requested heavy atoms | Result before refinement |
| --- | ---: | ---: | --- |
| 8K5R | 1.30 s | 136 | No modeled outliers; new boundary outlier at LEU265 |
| 1UA2 | 0.66 s | 91 | GLU50 modeled outlier and ASN56 boundary outlier persist |

8K5R's context searches retained 28, 12 and 40 candidates. Even after including
oxygen, the selected temporary contexts have maximum observed-backbone
displacements of 1.64, 2.54 and 1.81 Å. Those contexts are discarded on
transplant; these are **not accepted observed-coordinate movements**. Their
mismatch explains why complete join geometry needs attention before a coherent
candidate can be refined.

1UA2 closes in the initial database search and never reaches context ranking,
so this correction cannot address its current defect. Its earlier refinement
also creates new outliers from initially favored residues. The two paths need
boundary-inclusive conformation checks and candidate retry through refinement.

No full preparation, full-environment acceptance or new dynamics was performed.
Neither candidate qualifies the v1 release. Native internal-coordinate closure
with complete observed endpoint geometry was then tested separately on the
smallest failing 8K5R gap. Starting from the transplanted scaffold preserves
already malformed joins. Starting from one coherent native fragment removes
that defect, but none of its closed solutions remains within the unchanged
observed-backbone displacement cap. Both comparisons remain research-only;
multiple coherent fragments and full-atom validation are still required.

## Artifacts

- [Replay implementation](audit_carbonyl_ranking.py)
- [Final explicit opt-in comparison](/Users/ashujo/.cache/dynamol-research/v1-loop-release/carbonyl-ranking-v2/result.json)
- [Plan and source hashes](/Users/ashujo/.cache/dynamol-research/v1-loop-release/carbonyl-ranking-v2/plan.json)
- The earlier `carbonyl-ranking-v1` experiment and implementation snapshot
  remain unchanged. Its native results were repeated after making the ranking
  an explicit research opt-in; no failed candidate was admitted by that change.

Workers had a 180-second/4-GB supervisor limit and two numerical-library
threads. No bounds, physical force-field parameters or geometry/reference
thresholds were relaxed.
