# Representative loop panel: fresh preparation and validation

Updated 14 September 2026. The frozen panel has 15 structures, deliberately
covering representative short, medium and long gaps plus multiple-gap cases.
Every integer loop length is not required. Both 12-residue examples passed.

The first two-at-a-time sweep took **5 minutes 17 seconds**. That excludes
subsequent debugging and candidate comparisons. All ten repair candidates
reached a fresh common-preparation snapshot with their retained complex;
five unsupported inputs failed before that boundary without deleting chemistry.

Per the user's release direction, **metal-containing complexes are deferred**.
1HCK is kept as evidence but is not a v1 gate. Among nine metal-free repair
candidates, five pass complete static checks and four remain unresolved.
The five separate unsupported-input checks are not successful preparations.

| Structure | Gaps (residues) | Current outcome |
|---|---|---|
| 1UA2 | 12 | Passed; ATP and TPO retained |
| 1M17 | 12 | Passed; erlotinib retained |
| 4HJO | 4 | Passed; erlotinib retained |
| 1FPU | 5 | Passed; ligand retained |
| 4HHY | 8 | Passed; ligand and sulfate retained |
| 8K5R | 6, 4, 3, 4 | Unresolved backbone reference checks |
| 2CG9 | 7, 2 | Unresolved complete checks |
| 3HEG | 4, 11 | Unresolved reference/movement checks |
| 2HYY | 1, 3 | Unresolved complete checks |
| 1HCK | 4 | ATP–magnesium complex; deferred for v1 |
| 6A93 | 7, 3 | Expected refusal: covalent PLM |
| 2BEL | 3 | Expected refusal: incomplete ligand and inconsistent stereochemical reference |
| 2REN | 3, 5, 9 | Expected refusal: covalent NAG |
| 1RNE | 6 | Expected refusal: covalent NAG |
| 2ITY | 10, 13 | Expected refusal: gap exceeds 12-residue limit |

## Corrections and bounded follow-ups

- Native parser precision is compared at the source PDB's three-decimal
  coordinate precision, then exact parsed coordinates must survive raw-model
  construction. A real 1HCK parser difference of 0.0000105 Å had triggered
  the previous fixed tolerance; actual subsequent model movement still fails.
- The exact documented native Monte Carlo initialization failure is now an
  unavailable candidate, so another bounded candidate can run. Other native
  errors and identity/parameter failures remain fatal. The installed native API was also inspected directly; its explicit
  initial-backbone overload was exercised in the separate closed-fragment prototype.
- A private bounded fragment comparison uses one or two extra native context
  residues while exporting only original missing residues and immediate flanks.
  It produces closed models for difficult cases but has not produced a new
  accepted case. It does not weaken the final 1 Å caps, geometry, stereochemistry,
  backbone references, atom identity, retained environment or parameter checks.
- A further private 8K5R comparison seeds torsion sampling with complete closed
  fragments. All four candidates were rejected; this is not an app feature.

The focused worker/search/process regressions passed 62 tests; the subsequent
worker/context changes passed 75 tests. These are software tests, not proof of
native loop conformations. Saved bundles for all five accepted cases were
re-verified: all bound coordinate/topology/parameter/provenance bytes match,
all eight complete checks pass, and original-source flank caps pass.

## Remaining release work

The fresh full-worker 1UA2 handoff passed in 61.6 seconds,
including its final retained-complex, stereochemistry and geometry checks,
private dataset publication, and a copied/verified complete validation and
parameter bundle. The dataset has 4,825 atoms and all owned children exited.
The fresh 1M17 handoff also passed in 69.5 seconds, including verified dataset and parameter-bundle publication. Both private workflows completed and all owned children exited. The app's automatic fallback
remains disabled. UI/engine regressions, exact bundled installation
verification and video capture still follow that integration step.

The target is public v1 on 15 September morning, but the unresolved loop cases
and unfinished generic app integration prevent an unconditional readiness claim. Modeled
loops remain provisional; successful static checks do not establish native
conformation, ensemble accuracy or MD stability.

Detailed decisions and evidence hashes are in [panel-summary.json](panel-summary.json).
All failed/superseded results remain preserved. No build, release, MD or advanced
covalent/metal parameterization was launched by this panel sweep.

A later [targeted Abl–imatinib comparison](ABL-IMATINIB.md) now passes complete
static checks for both 2HYY gaps, with imatinib retained and 0.767 Å maximum
observed-heavy movement. This does not retroactively change the automated sweep
above: general selection and fresh full-worker qualification are still pending.

## What distinguishes the unresolved cases

All four have multiple missing segments, including very short gaps. The latest
complete checks primarily reject backbone conformations in reconstructed residues
or their immediate joins. 3HEG also exceeds the source movement limit and sometimes
fails gross geometry. Atom identity, retained chemistry and parameter integrity
are not the failing checks. 2CG9 uses the selected Sba1 protein chain without a
ligand, so ligand chemistry is not a shared cause.

RCSB reports 8K5R at 3.75 Å; that is relatively low resolution and might contribute
to difficult endpoints, but this has not been established as the cause. 2HYY is a
2.4 Å Abl–imatinib structure with one- and three-residue gaps. These are relevant
inputs for the intended app, not grounds to label all four abnormal. Source pages:
[8K5R](https://www.rcsb.org/structure/8K5R),
[2HYY](https://www.rcsb.org/structure/2HYY),
[2CG9](https://www.rcsb.org/annotations/2CG9),
[3HEG](https://www.rcsb.org/experimental/3HEG).
