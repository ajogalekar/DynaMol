# Complete-complex development track

This track continues covalent and metal chemistry independently of the next
DynaMol v1 release. The user split these tracks on 14 September 2026: the v1
priority is loop repair through 12 residues per gap; advanced chemistry remains
active research and is not a release blocker. Nothing here enables covalent
preparation in the app or promotes a research candidate.

Support is assessed for the complete selected complex, including required
cofactors, ions, covalent adducts and retained waters. The original frozen
15-covalent and 15-metal panels are unchanged.

## First completed diagnostic: KRAS–sotorasib–GDP–Mg

`audit_mg_coordination.py` independently replays the corrected existing
`6oim-stability-v2` OpenMM trajectory and `6oim-gromacs-continuation-v2` native
GROMACS trajectory. It ran in 0.76 seconds; it launches no MD or QM.

The final evidence directory is
`/Users/ashujo/.cache/dynamol-research/complete-complex-track/6oim-coordination-v4`.
It contains `result.json`, 121 frame rows in `frames.csv`, the exact implementation
snapshot and an output manifest. Source trajectories, inputs and reports are
hash-checked and unchanged after analysis. The original 100 ps OpenMM test has
10 ps heating, 20 ps restrained NPT and 70 ps unrestrained NPT; the 20 ps GROMACS
continuation contributes 21 saves including its starting frame. These are not
121 independent samples or two independent physical replicates.

The complete native model has 27,477 atoms and retains sotorasib, GDP, Mg and
all 207 deposited waters, with no deposited environmental molecule excluded.
The joined cysteine–drug residue has 82 atoms and a charge of approximately
−0.000007 e; GDP has 40 explicitly identified atoms and charge −3 e; Mg is +2 e.
These inventory/charge checks do not establish an accurate charge distribution.
[RCSB identifies 6OIM](https://www.rcsb.org/structure/6OIM) as the 1.65 Å human
KRAS G12C–AMG 510 structure with GDP and magnesium.

All six donor identities are resolved back to the retained mmCIF. Label and
author chain identifiers are both recorded: label chains B/C/E for Mg/GDP/water
correspond to author chain A. The donor atoms have full deposited occupancy and
unambiguous atom identities. Their vectors relative to Mg are unchanged by the
initial assembly, within recorded coordinate precision.

| Mg donor | Deposited distance, Å | OpenMM final 70 ps mean, Å | GROMACS continuation mean, Å |
| --- | ---: | ---: | ---: |
| Ser17 OG | 2.164 | 2.098 | 2.112 |
| GDP302 O2B | 2.157 | 1.891 | 1.899 |
| Water405 O | 2.161 | 2.026 | 2.048 |
| Water415 O | 2.166 | 2.054 | 2.041 |
| Water440 O | 2.166 | 2.066 | 2.038 |
| Water447 O | 2.164 | 2.122 | 2.127 |

The audit searches **all 8,746 oxygen/nitrogen atoms**, including added bulk
water. Every frame has the same six nearest donors and a count of six at each
of the descriptive 2.4, 2.8 and 3.2 Å cutoffs. No other atom enters 3.2 Å; the
closest nonoriginal donor reaches 3.221 Å in OpenMM and 3.251 Å in GROMACS.
This is a saved-frame result, not evidence about exchange between saves or
long-timescale water-exchange kinetics.

All 15 donor–Mg–donor angles are reported. The mean Ser17–Mg–GDP O2B angle is
98.87° during unrestrained OpenMM and 97.39° in GROMACS, versus 90.51° deposited.
The angular RMS difference from an ideal octahedron is 3.96° deposited,
6.59° in unrestrained OpenMM and 6.17° in GROMACS. The original three opposite
donor pairs define the ideal-angle comparison; no fit or acceptance threshold
is introduced.

**Interpretation:** the persistent GDP contact contraction is accompanied by
changes in the retained coordination geometry. An unnoticed donor replacement
does not explain it at the saved frames. This strengthens the reason to assess
the joint ion/nucleotide/protein/water interaction model. It does not prove the
force field is incorrect from one crystal pose, establish the correct
solution-phase distance, or validate the covalent force field.

## Verification and preserved failures

Double-precision periodic vector calculations agree with independent MDTraj
distance/angle calculations across all 121 saves to within 0.00000029 Å and
0.000435°. Analytic checks cover ideal octahedral geometry, invariance of angles
to radial changes, periodic boundaries and invalid-geometry refusals. Only
orthogonal boxes are supported; other boxes fail explicitly.

The native handoff's atom-order/bond check is required, all mapped native atom
names/elements are verified, and the saved OpenMM endpoint is compared through
the actual exported GRO input to the first native GROMACS frame. Coordinate
differences at both steps are below 0.000001 nm. The native starting box agrees
with the GRO box within 0.00000023 nm; the GRO box differs from the source by
the explicitly measured five-decimal text rounding.

The failed first adapter attempt (`v1`, unsupported Gemmi row slicing) remains
preserved. The completed first geometry audit (`v2`) remains unchanged. A
stricter boundary comparison (`v3`) incorrectly used the source box to undo
translations made using the rounded GRO box; its failure and original script
remain preserved. The corrected `v4` verifies both saved boundary steps with
the same 0.000002 nm coordinate tolerance and explicitly checks decimal box
rounding. No chemistry tolerance was changed and no source failure was
reclassified.

Frozen panel SHA256 values remain:

- Covalent: `5fdb1dba8cfa43e90c3396c87bde1cdc2a1a8a7091ec813cb120255ad70a613a`
- Metal: `28db88c10ce5a80dde8513d9c02542c2aa7ffda4910d8bbaef64ce9f312af694`

## Completed applicability milestone

The pinned microMg compatibility inventory has now been completed; see
[its report](MICROMG-APPLICABILITY.md) and `micromg-applicability.json`.
The complete current GDP/protein/drug/water model is not eligible for direct
substitution. The report identifies exact source aliases, omitted named pairs,
charge differences, a small TIP3P LJ variant and a licensed-paper implementation
route distinct from copying unlicensed repository files.

The original milestone below records the assessment's intended scope.

Construct a component-by-component compatibility inventory for a published
alternative Mg model against the **actual retained GDP and Ser17 donor atom
types/charges**, including every pair override and its water model. The
[microMg author implementation](https://github.com/bio-phys/Magnesium-FFs)
explicitly targets TIP3P and named RNA/ion combinations; its RNA success does
not validate direct transfer to this GDP/protein model. The existing
[Mg–GDP assessment](../covalent-v1-focused/MG-GDP-ASSESSMENT.md) records those
applicability and redistribution gaps.

Only after that inventory is coherent should a bounded alternative be tested
with explicit native pair-energy/force transport checks in both engines and
complete-complex comparisons. Do not tune a radius to this single crystal
distance, remove GDP/Mg or constrain the covalent drug bond to obtain a pass.
Broader covalent charge/torsion validation and the frozen panels remain open.

Reproduce the diagnostic with the isolated existing research Python runtime:

```sh
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
  /Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python \
  docs/audit/advanced-chemistry/complete-complex-track/audit_mg_coordination.py \
  --root /Users/ashujo/.cache/dynamol-research/covalent-v1-focused \
  --output /Users/ashujo/.cache/dynamol-research/complete-complex-track/NEW-UNUSED-DIRECTORY
```

Existing output directories are never overwritten.
