# Backbone conformation exposed a refinement defect

The 6JXT research refinement can turn a statistically favored fragment
conformation into a Ramachandran outlier while passing its existing gross
bond, angle, peptide-plane and stereo checks. The current research export path
now checks the saved candidate with native CCTBX before creating a protein
intermediate. It rejects the previously screened candidate; no old artifact is
modified and no app integration is claimed.

## Evidence from saved coordinates

`audit_loop_conformation.py` measured local backbone angles and nearby atoms in
the two completed 6JXT diagnostics and the independent 6JX0 control. THR751
starts at approximately phi/psi +54.9/+171.7 degrees in the refined 6JXT model,
versus -116.3/+26.7 degrees in the control. The 6JXT threonine sidechain oxygen
stays within 3 Å of ALA750 carbonyl oxygen in 93/100 and 96/100 saved frames;
the respective mean distances are 2.76 and 2.72 Å. The control instead has
persistent contacts with GLU872. These are descriptive distances, not assigned
hydrogen bonds or proof that the control is the missing loop's native state.

Comparing refinement stages locates the change more precisely. In the coherent
fragment seed, THR751 has phi/psi -102.6/+119.1 degrees and is favored by the
native CCTBX reference. After refinement it becomes an outlier. The complete
selected region contains 15 modeled/context residues:

| Saved structure | Native-reference outliers |
| --- | --- |
| Coherent fragment seed | PRO753 |
| Previously screened refined candidate | THR751, SER752, ALA866, GLY874 |
| Conformation-protection experiment | SER752, GLY874; also fails peptide omega |

Thus the seed itself is not wholly validated. The specific good-to-outlier
THR751 transition and the additional context outliers identify a missing
quality check; they do not establish that one constraint or interaction is the
sole cause of the dynamics behavior. Source E865A/E866A/K867A identities remain
unchanged.

Evidence folders are `6jx0-conformation-v1`, `6jxt-conformation-v3`,
`6jxt-conformation-v4`, `6jxt-seed-rama-v2`, `6jxt-refined-rama-v2` and
`6jxt-protected-rama-roundtrip-v2` in the focused cache.

## Native reference and independent integration check

The separate cache runtime contains CCTBX 2025.11 and six 1.17.0; the main
research runtime and application bundle were not changed. `rama_reference.py`
uses the compiled `mmtbx_validation_ramachandran_ext.rama_eval`, including its
six residue classes and existing cutoffs. For general residues, the native
outlier boundary is 0.0005; the refined THR751 score is approximately 0.000238.
[CCTBX scoring implementation](https://github.com/cctbx/cctbx_project/blob/dca90aab7a53cd1a747d425542a5724eea1b17e9/mmtbx/validation/ramachandran/rama_eval.h).

The installed evaluator/table headers are byte-identical to pinned source
commit `dca90aab7a53cd1a747d425542a5724eea1b17e9`. Native binary SHA-256 is
`b7a531377948dc6c1534ffe2f4041cfe28c38eac6817e12e0808e828da28935c`.
Source headers, retrieval hashes and the source license are retained under
`ramachandran-reference-v1`. No redistribution decision is implied.

The adapter initially used MDTraj float32 dihedrals. All classifications agreed
with full CCTBX, but a steep part of the table amplified a 0.00023-degree angle
difference to a 0.0000286 score difference. That exceeded the predeclared numeric
comparison tolerance. The original reports remain preserved; the adapter now
uses double-precision angles. The complete CCTBX PDB parser and `ramalyze`
independently agree on all 45 selected residue checks across three structures,
with maximum angular difference below 5e-13 degrees and score difference below
3.5e-14. Evidence: `cctbx-reference-parity-v2`.

Full-reference verification needs `PYTHONPATH` set to
`/Users/ashujo/.cache/dynamol-runtimes/cctbx-rama-reference-v1` and `LIBTBX_BUILD`
set to its `libtbx/core/share/cctbx` directory. The light native evaluator itself
does not need the full parser environment. These are research-runtime settings,
not application setup instructions.

## First corrective experiment remains rejected

`6jxt-conformation-protected-v1` adds temporary phi/psi tethers for allowed seed
conformations during the existing anchor/environment refinement. Existing seed
outliers are left free to improve. Tether strength is 1,000 kJ/mol, matching the
existing final peptide construction strength; this is a geometric construction
choice, not a physical torsion parameter or a fitted statistical energy.
The native MD force field is unchanged and these forces are never exported.

The bounded experiment finished in 4.02 seconds, with sampled peak memory
181 MB. THR751 remains favored at phi/psi approximately -78.7/+107.0 degrees,
but its peptide omega departs 42.24 degrees from planarity, exceeding the
unchanged 35-degree screen. SER752 and GLY874 remain Ramachandran outliers.
Stereo and the 1 Å observed-context movement cap pass (maximum 0.813 Å).
The result is rejected and supplies no coordinates to another complex or run.

`export_loop_intermediate.py` now independently scores the saved PDB and
requires all selected modeled/context residues to be scorable without outliers,
in addition to the original geometry, stereo and environment checks. The known
four-outlier candidate is refused before any export is created; the control's
ALA750/THR751 reference check passes. Evidence: `rama-export-admission-v1`.
This is a static preparation check, not a new per-frame MD cutoff.

## Bounded native alternatives: completed, no candidate admitted

`compare_fragment_seeds.py` retains the coordinates and identities of up to
40 native closed, ring-filtered fragments per gap, considering up to two context
residues. The preserved worker continues to enforce source sequence and atom
mapping. The first enumeration retained 40, two and 19 alternatives for the
three gaps; the first count reaches the explicit cap and is not an exhaustive
inventory. Native backbone scores guide candidate selection; they are not an
accuracy or acceptance score.

`run_fragment_comparison.py` freezes seven trial choices before refinement:
the native-score and anchor-displacement leaders, then up to two additional
native-ranked choices at each gap with the other choices held at their native
leaders. Every trial uses the existing conformation-protection, anchor,
environment and peptide construction strengths and all original final limits.
Every seed and refined PDB receives a fresh native Ramachandran check.

All seven trials finished in 42.15 seconds. Four seeds have no selected
backbone-reference outliers. One refined candidate also has no such outliers,
but still fails peptide planarity; all seven remain rejected. Stereo and the
1 Å observed-context displacement cap pass in all seven. This is not a new
accepted protein intermediate, native complex or MD run.

| Trial | Seed outliers | Refined outliers | Geometry accepted |
| --- | --- | --- | --- |
| Native-score leaders | 0 | 1 | No |
| Anchor-displacement leaders | 1 | 2 | No |
| Gap 1, next native choice | 0 | 1 | No |
| Gap 1, third native choice | 0 | 0 | No |
| Gap 2, next native choice | 1 | 1 | No |
| Gap 3, next native choice | 1 | 1 | No |
| Gap 3, third native choice | 0 | 1 | No |

Evidence is preserved under `6jxt-fragment-alternatives-v1`, including the trial
plan, all native alternatives, source hashes, implementation snapshots, process
logs, bounded process supervision and both sets of saved-PDB reference scores.
Source hashes and the original frozen 15+15 panel hashes remain unchanged.

The resulting coverage audit identifies a specific next correction to test.
`backend.loop_refinement._minimize_attempt` applies temporary peptide-planarity
restraints only to peptide links touching rebuilt residues. Observed context
residues are also mobile, but peptide links solely between them are omitted.
SER752–PRO753, which exceeds the planarity screen in five of these trials, is
one such omitted link. All native physical peptide force-field terms are still
present; this is a construction-restraint coverage gap, not missing simulation
parameters. Coverage alone does not prove the cause of every failure.

## Context peptide coverage: paired experiment remains rejected

The backend helper now has a research opt-in `preserve_context_peptides` option,
off by default. It adds temporary construction targets for peptide links
touching explicitly mobile observed context. Observed cis/trans basins come
from the original reference coordinates, including observed non-Pro cis links;
modeled-link policy and the existing 200/1,000 kJ/mol strength schedule remain
unchanged. The UI and normal preparation path do not enable the option.

All 27 `test_loop_refinement.py` checks pass in the separate research runtime.
They exercise actual CPU minimization, default-path behavior, full mobile-link
coverage, original cis/trans preservation even when the seed has the opposite
state, PRO/HYP/non-Pro identity policies, fixed atoms, and exclusion of temporary
forces from fresh simulation systems. This is software validation on analytical
potentials, not validation of protein chemistry. The repository virtualenv
pytest launch stalled before collection and was terminated; pytest 8.4.2 was
installed in a separate cache target for the successful run.

`replay_fragment_refinement.py` applies both the original and corrected coverage
to each of the seven saved seeds, using identical initial coordinates and
hydrogens within every pair. All 14 refinements completed in 61.16 seconds;
none passes all final checks. Corrected coverage produces three gross-geometry
passes, but all three fail SER752 stereochemistry and reference backbone checks.
Another candidate eliminates peptide-plane failures yet stretches the SER752
N–CA bond beyond its unchanged bound. Other candidates still fail peptide or
backbone-reference checks. All 14 remain within the original 1 Å context cap.
No geometry-only success is promoted.

Evidence: `6jxt-context-peptide-replay-v1`, including `plan.json`, `result.json`,
`summary.json`, `implementation-checks.json`, frozen implementation snapshots,
all original/refined coordinates and per-process supervision. Original seeds
and earlier failed outputs remain intact. No native assembly, MD or QM followed.

## Native closure probe and a more specific boundary condition

The next direction is native torsion-based closure rather than further Cartesian
restraint tuning. ProMod3 exposes both CCD and KIC and supports manipulating
backbones through internal coordinates. Its published workflow also supports
full-atom candidate scoring and sidechain reconstruction; these are established
library capabilities, not evidence that they solve this case.
[ProMod3 methods paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008667).

The installed ProMod3 3.6.0 KIC interface was checked on the intact six-residue
fixture before an EGFR probe. The control returns native closures that reproduce
all observed backbone coordinates within 0.000007 Å. The initial adapter's
ResidueView/ResidueHandle error is preserved in `6jxt-kic-feasibility-v1`; the
corrected control and probe are in `6jxt-kic-feasibility-v2`. Native runtime
errors establish that pivot residues must be internal; future probes should
exclude the two stem indices rather than repeat those invalid combinations.

For one saved EGFR seed covering residues 749–753, KIC returns eight raw
solutions, none within the existing observed-backbone cap. The best maximum
displacement is 2.0705 Å at observed GLU749 carbonyl O; SER752 N moves 1.8806 Å.
The N/CA/C stem atoms alone match within 0.0195 Å. Thus matching those three
stem atoms does not ensure agreement with the observed carbonyl orientation.
The current fragment ranking counts N/CA/C displacement but excludes O.

In `6jxt-kic-carbonyl-v1`, a native psi rotation first aligns the seed's
N/CA/C/O dihedral to the observed carbonyl orientation. The needed change is
about 165.36 degrees, with native angular error below 0.000004 radians. KIC then
returns no solutions for this one seed and its fixed internal geometry. This
does not establish that the gap is impossible to model: it identifies an
additional observed-stem condition that a candidate search must respect.
No KIC output is an accepted full-atom protein or preparation intermediate.

Next test bounded native torsion sampling/closure with the complete observed
stem geometry, including carbonyl orientation, represented before scoring and
refinement. Compare multiple native fragments and full-atom context, preserving
source sequence, all retained molecules and every current admission limit.
Do not keep increasing construction forces to compress an incompatible seed.
Broader loop/panel and physical covalent/Mg–GDP validation remain required.
