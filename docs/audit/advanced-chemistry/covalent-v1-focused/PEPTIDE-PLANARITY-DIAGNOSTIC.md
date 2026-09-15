# Peptide planarity: control comparison and unresolved loop bias

The original 6JXT full-complex tests remain rejected under their predeclared
checks. The first stopped at 1 ps; the run freeing the remodeled/context region
from warm-up position restraints stopped at 7 ps. Both failures concern the
rebuilt ALA750–THR751 peptide. The second sampled a 38.383° deviation from the
nearest cis/trans plane; the unchanged static screen permits at most 35°.
No finite-force, standard-protein stereo, clash or bond-angle failure accompanied
that event. Neither failed run is reclassified as successful.

## Read-only comparison with completed controls

`audit_peptide_planarity.py` independently reads the saved trajectories and
complete native topologies. It computes every protein backbone peptide omega,
with whole coordinates and no wrapping. The original inputs and saved files are
hash-bound and unchanged. All initial models had zero peptides beyond 35°.

| Case | Saved frames | Peptide bonds | Bond-frame samples beyond 35° | Bonds ever beyond 35° | Maximum deviation |
| --- | ---: | ---: | ---: | ---: | ---: |
| Rejected 6JXT v2 | 7 | 287 | 1 / 2,009 | 1 | 38.38° |
| EGFR control 6JX0 | 100 | 288 | 7 / 28,800 | 7 | 41.16° |
| BTK 5P9J | 100 | 265 | 25 / 26,500 | 12 | 49.58° |
| Corrected KRAS 6OIM | 100 | 169 | 8 / 16,900 | 4 | 42.11° |

These observations show that applying the static preparation cutoff as a
per-frame MD acceptance rule can also flag completed controls. They do not prove
that those excursions are physically accurate, nor supply a replacement cutoff.

The specific EGFR loop remains concerning: in the first seven ps its mean
deviation is 27.96°, versus 2.65° for the corresponding sequence-identical peptide
in the 6JX0 control during the same heating interval. The latter averages 7.34°
over its whole 100 ps and never exceeds 21.17° at that peptide. The control is a
different deposited structure and genotype, not an experimental trajectory or
proof of the native conformation of the missing 6JXT loop.

Evidence folders in the focused cache are `6jxt-peptide-planarity-diagnostic-v1`,
`6jx0-peptide-planarity-diagnostic-v1`, `5p9j-peptide-planarity-diagnostic-v1` and
`6oim-peptide-planarity-diagnostic-v1`. Each retains per-bond statistics, the
complete angle arrays and source hashes.

## Sequence and template check

The original deposited CIFs confirm ALA750/THR751 in both structures. 6JXT
declares E865A/E866A/K867A; those alanines are preserved, not replaced with the
glutamates/lysine in the control. [6JXT entry](https://www.rcsb.org/structure/6JXT)
and [6JX0 entry](https://www.rcsb.org/structure/6JX0) describe WT gatekeeper and
T790M complexes, respectively. Titles do not replace exact source sequence checks.

A proper rigid fit using backbone atoms of observed residues 748, 749, 752 and
753 gives a local reference-to-target RMSD of 1.166 Å and maximum displacement
2.288 Å. Extending the fit to four observed residues on each side gives 1.773 Å
RMSD and maximum displacement 3.369 Å. Thus copying the control coordinates
unchanged would not satisfy the current 1 Å observed-context bound. No reference
coordinates or sequence substitutions were applied.

Published work supports interpreting peptide nonplanarity in context, but does
not establish a DynaMol runtime threshold: [Berkholz et al., structural study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3258596/)
reports conserved nonplanarity, while [Rick and Cachau, MD study abstract](https://scholarworks.uno.edu/chem_facpubs/7/)
describes an environment-dependent peptide model different from our Amber
candidate. These sources are context, not replacement parameters or validation
of this repaired loop.

## Bounded diagnostic sampling, without admission

`6jxt-stability-diagnostic-v3` tests whether the specific loop's deformation
relaxes or persists through heating and unrestrained equilibration. The starting
model, native force field, seed, schedule and free-remodeled-region choice match
v2. No extra planarity force or attachment constraint is added.

The optional diagnostic mode records each finite peptide excursion as a failed
geometry frame, retaining the original 35° threshold. Initial and minimized
frames must still pass. Stereo failures, nonfinite geometry, other geometry
errors and cis/trans basin changes remain stopping conditions. This mode can
finish only as `diagnostic_trajectory_complete`; `diagnostic_only=true` and
`qualified_for_handoff=false` prevent it becoming a qualified GROMACS endpoint.
All excursions remain in both frame reports and a collected event list.

Six checks using the actual rejected frame and intentionally corrupted copies
confirm strict-mode refusal, failed-frame recording without promotion, refusal
to waive preparation, and rejection of distorted bonds, inverted stereochemistry
and nonfinite coordinates. Two handoff checks reject diagnostic/unqualified
endpoints before invoking any engine. Evidence:
`loop-monitor-diagnostic-admission-v1`. Future hard failures also save the engine
state and binary checkpoint for diagnosis. No application acceptance rule has
changed, and no loop or force-field accuracy is claimed.

## Completed 100 ps diagnostic and matched comparison

`6jxt-stability-diagnostic-v3` completed 100 ps in 2,782 seconds with 42,900
atoms. Its mean temperature during the final 70 ps is 300.63 K, final density
1.01768 g/mL, unconstrained attachment range 1.741–1.945 Å and maximum aligned
C-alpha RMSD 1.177 Å. All 100 saved frames and the final box/coordinates agree
with the engine state within DCD precision. The portable restart preserves
Reference energy and forces exactly and explicitly remains unqualified.

Twenty loop/context geometry frames remain flagged: 6, 7, 13, 20, 21, 22, 23,
24, 28, 29, 30, 43, 45, 52, 63, 71, 85, 92, 93 and 98 ps. These are finite
peptide-omega excursions with retained stereo, other geometry and cis/trans
basins. Independent replay of all 18 monitored peptides at each saved ps agrees
with the live angle reports within 0.00055 degrees and confirms every flagged
frame. Across all 287 protein peptides, 85 of 28,700 bond-frame samples exceed
35 degrees, involving 19 peptides; the maximum is 59.13 degrees. This is more
frequent than in the existing controls and is not interpreted as validation.

The originally problematic ALA750–THR751 peptide shows a time-dependent change:

| Saved interval | Restraints on remaining protein | Mean planarity departure | Maximum | Samples over 35 degrees |
| --- | --- | ---: | ---: | ---: |
| 1–10 ps heating | 1,000 kJ/mol/nm² | 30.02° | 50.27° | 2 / 10 |
| 11–30 ps NPT | 100 kJ/mol/nm² | 35.03° | 55.36° | 9 / 20 |
| 31–100 ps NPT | None | 11.20° | 35.20° | 1 / 70 |

The rebuilt and adjusted-context region is free of those position restraints
throughout. Reduced distortion coincides with release of the remaining
restraints, but elapsed relaxation time is a confounder. This observation does
not establish the original missing loop's native conformation.

Evidence: `6jxt-peptide-planarity-diagnostic-v2` retains every protein peptide
angle; `6jxt-completed-loop-diagnostic-v1` cross-checks monitored angles. The
additional `6jxt-completed-loop-diagnostic-v2` verifies phase labels against
each saved positional-restraint value and preserves identical statistics.
The earlier first-70-ps analysis remains in `6jxt-live-loop-trend-70ps-v1` as a
partial-run observation. `saved-run-audit.json` and `portable-restart-v1` under
the completed run retain diagnostic and unqualified flags; neither clears the
20 recorded failures.

The matched `6jxt-stability-diagnostic-v4` comparison releases position
restraints after the same initial minimization, before heating. It uses the
same original complete model, mapped intermediate, native parameters, seed and
temperature/pressure schedule. Two CPU threads and sampled two-hour/4 GB bounds
are unchanged. It is diagnostic-only, with the same geometry/stereo/basin
stopping rules, and cannot qualify a handoff even if it finishes. The runner
snapshot matches its launch hash. CPU runs are not assumed bitwise identical.
For comparison, use the same time windows and final-70-ps temperature
mean, since the all-unrestrained mean in v4 includes its heating interval.

### Completed v4 comparison: mixed effects, no warm-up fix accepted

The v4 run completed all 100 ps in 2,791 seconds. Its mean temperature over
31–100 ps is 300.78 K, final density 1.01570 g/mL and unconstrained drug attachment
range 1.749–1.957 Å. The saved trajectory, final coordinates/box and zero-restraint
portable restart are independently verified; Reference reload energy and force
differences are zero. Diagnostic and unqualified flags remain explicit.

The 13 flagged loop/context frames are 4, 5, 14, 22, 26, 28, 33, 34, 43, 62, 87,
92 and 97 ps. All 18 monitored peptide angles agree with the independent DCD
replay within 0.00043 degrees. Across all protein peptides, 66 of 28,700 samples
exceed 35 degrees, involving 24 peptides with a maximum of 46.37 degrees. Fewer
total flags but more affected peptides is not a uniform improvement.

| ALA750–THR751 interval | v3 mean departure | v4 mean departure | v3 / v4 samples over 35 degrees |
| --- | ---: | ---: | ---: |
| 1–10 ps | 30.02° | 27.48° | 2 / 2 |
| 11–30 ps | 35.03° | 17.00° | 9 / 0 |
| 31–100 ps | 11.20° | 23.14° | 1 / 4 |

The native input coordinates match exactly. Initial minimization uses the same
protocol but is not bitwise identical: target omega differs by 0.416 degrees
between minimized poses. A single pair of short trajectories does not establish
causality or native loop accuracy. In particular, the final phase is worse at
the target peptide when restraints are released earlier. No new default is
accepted and no further identical random rerun is planned.

Evidence: `6jxt-peptide-planarity-diagnostic-v3`,
`6jxt-completed-loop-diagnostic-v3`, and `6jxt-warmup-comparison-v1` in the focused
cache. `compare_loop_diagnostics.py` verifies the original model hashes, mapped
intermediate, launch-bound runner snapshots, recorded restraint schedule and
all peptide identities/time windows before comparing. Both failures and all
diagnostic outputs remain preserved.

### Native local terms match the independently prepared EGFR control

`compare_native_loop_terms.py` compares every atom property and bonded term
touching chain A residues 750 and 751 in native `6jxt-complex-v1` against
`6jx0-complex-v1`, with source residue/atom identities including neighboring
term endpoints. Exact equality holds for all 24 atom records (types, elements,
masses, charges and Lennard-Jones values), 25 bonds, 47 angles, 137 proper torsion
terms and six improper terms. Evidence is
`6jxt-target-native-terms-comparison-v1`, with both full inventories and input
hashes retained.

This rules out a mismatch in those assigned local terms between the two
assembled models. It does not establish that the loop coordinates or surrounding
nonbonded environment are correct, nor that the shared force field is physically
accurate. Local backbone conformation and its environment are the next diagnostic
targets; app acceptance and the frozen preparation limits remain unchanged.

The subsequent [native Ramachandran audit](RAMACHANDRAN-REFINEMENT.md) identifies
the favored-to-outlier THR751 transition during construction refinement, along
with three other selected backbone outliers. This is a missing preparation
quality check rather than evidence supporting a relaxed per-frame omega limit.
The research exporter now refuses that candidate before producing an
intermediate; the original completed diagnostic trajectories remain unchanged.

## Separate unrestrained minimization did not improve the starting pose

`6jxt-native-relaxation-v1` is a separate experiment with the original native
complex, no position/construction restraints, and a 1,500-iteration limit.
It finished in 242 seconds and passed its static loop/stereo checks. Its target
peptide nevertheless ends 33.41 degrees from planarity, worse than the 23.32
degrees after the previous restrained minimization. The maximum remaining
force is 2,644 kJ/mol/nm; completion at the iteration cap is not evidence of
convergence.

Relative to the assembled MD input, observed context moves up to 4.391 Å and
other observed protein up to 4.016 Å. Those values do not satisfy the previous
1 Å local preparation cap and are not cumulative deviations from the original
deposition. The experiment preserves original source files and its complete
displacement breakdown, but its coordinates are not promoted or used for v4.
