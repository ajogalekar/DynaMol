# Focused covalent v1 validation

**Direction update, 14 September:** this research now continues separately
from the next v1 release, which is gated on the loop fix and release audit.
See [the two-track plan](../../v1-loop-release/PLAN.md). The historical directory
and artifact names are preserved. Support must cover each complete complex,
including its required ions, ligands and cofactors.

Scope approved by the user: known human KRAS–sotorasib (6OIM),
BTK–ibrutinib (5P9J), and EGFR–osimertinib (6JXT) adducts.
Unknown-chemistry QM controls remain a v2 item. These research outputs have
not been integrated into the app or approved for release.

**Latest preparation check:** the opt-in context-peptide coverage correction
passes 27 software tests, but none of 14 paired refinements passes every final
geometry, stereo and backbone-reference check. A native KIC probe identifies
an additional source-stem condition: observed carbonyl orientation is not fixed
by matching N/CA/C alone. Next test native torsion sampling/closure with complete
observed stem geometry and full-atom screening, preserving all current limits.
See `RAMACHANDRAN-REFINEMENT.md`, `6jxt-context-peptide-replay-v1` and the native
KIC probe artifacts. The normal app workflow does not enable the option; no
candidate was admitted or used for a new MD/QM run.

**Latest correction:** the original KRAS run omitted 12 GDP hydrogen constraints
because scoped atom types lost native element metadata. It is preserved but
superseded as a protocol qualification. Explicit source-element declarations
are fixed in `gdp-scoped-v6` / `6oim-complex-v3`; the corrected 100 ps run completed
in `6oim-stability-v2`, including its saved-trajectory and restart audits.
KRAS, BTK and the EGFR control also completed 20 ps GROMACS continuations.
The original 6JXT loop-water clash is resolved by including retained molecules
as fixed steric obstacles during construction. The new v3 candidate passes
independent geometry/stereo/environment checks, and its complete native complex
contains 42,900 atoms with all original drug, chloride and waters retained.
Its first two OpenMM runs stopped at 1 and 7 ps on peptide nonplanarity.
Control analysis also finds occasional static-cutoff excursions, but this
specific loop is more persistently distorted during warm-up. The diagnostic-only
v3 run completed 100 ps, with all saved frames and its restart independently
verified; 20 loop-geometry flags remain failures. The target peptide's mean
deviation falls from 35.03 degrees during restrained NPT to 11.20 degrees after
restraint release. The completed v4 diagnostic releases restraints before heating:
it has 13 flagged frames and improves early equilibration, but the target peptide
is more distorted in the final 70 ps (mean 23.14 versus 11.20 degrees). This is
not an accepted warm-up fix. Its saved trajectory and restart are verified.
Native atom/bond/angle/torsion terms at ALA750/THR751 match the EGFR control
exactly; local conformation and environment still need assessment. Neither
diagnostic qualifies a handoff or app admission. Wholly unrestrained minimization did not supply an improved starting
pose and was not promoted. See [planarity diagnostic evidence](PEPTIDE-PLANARITY-DIAGNOSTIC.md)
and [handoff and loop evidence](HANDOFF-AND-LOOPS.md).

**New preparation finding:** refinement turns the favored THR751 fragment
conformation into a native CCTBX Ramachandran outlier. The previously screened
candidate has four selected backbone outliers that gross geometry/stereo checks
missed. The research exporter now checks for these before creating an
intermediate. A first corrective experiment remains rejected; a validated
replacement and app integration are still pending. See
[backbone refinement evidence](RAMACHANDRAN-REFINEMENT.md).

## Actual results

All three capped drug–cysteine candidates have native joint AM1-BCC charges
and native PREPGEN residue templates. Each assembled between two peptide
neighbors with its complete retained atom/bond inventory, two peptide links,
eight boundary angles and 44 boundary torsion terms. Final residue charges
agree with the declared neutral/neutral/+1 states within 7e-6 e. These are
AM1-BCC/PREPGEN candidates, not constrained RESP results. Canonical ff14SB
backbone types are retained; the backbone charges differ from standard
ff14SB by approximately 0.15 e at the most affected atom. Physical assessment
of the resulting attachment flexibility is still required.

KRAS and BTK passed the existing ProMod3 loop-repair workflow for their
three-residue internal gaps. The complete native BTK complex contains 36,093
atoms with water/counterions. The complete KRAS complex contains 27,474
atoms in the superseded build; the corrected build has 27,477 atoms, including
GDP, Mg(II), and all original waters. Both cases have verified
source-to-native atom mappings after solvation.

Three 100 ps OpenMM CPU tests completed, including the additional 6JX0 EGFR
control, with results in their respective
`result.json`, `saved-run-audit.json` and `observations.csv`. The protocol is 10 ps of NVT heating,
20 ps of restrained NPT equilibration, and 70 ps of unrestrained NPT. There
is no constraint or special bond restraint holding the drug attachment.
Completed output requires `result.json`; a live progress file is not success.
The checks track finite energies/forces, attachment geometry, adduct bonds,
mapped stereocenters, temperature, density and aligned C-alpha RMSD. KRAS
also records initial magnesium-donor distances. These short tests cannot
establish thermodynamic convergence or validate the force field physically.

| Check | BTK–ibrutinib | KRAS–sotorasib | EGFR–osimertinib (6JX0) |
|---|---:|---:|---:|
| Total / unrestrained time | 100 / 70 ps | 100 / 70 ps | 100 / 70 ps |
| Mean unrestrained temperature | 300.22 K | 300.38 K | 300.51 K |
| Final density | 1.029 g/mL | 1.024 g/mL | 1.016 g/mL |
| Attachment bond sampled range | 1.720–1.953 Å | 1.762–1.956 Å | 1.763–1.935 Å |
| Maximum aligned C-alpha RMSD | 0.976 Å | 0.738 Å | 1.124 Å |
| Attachment torsion sampled span | 37.6° | 42.4° | 37.4° |

The KRAS column uses corrected `6oim-stability-v2`; the original v1 output is
preserved as superseded evidence. Its verified restart is `portable-restart-v1`
within the corrected run; BTK and the EGFR control use their `portable-restart-v2`.

All 100 saved frames per run were independently checked. Final coordinates and
boxes match the engine state within DCD precision; sampled stereocenters retain
their signs. The exported Systems have no constraint on the drug attachment.
KRAS retains all six initial magnesium donors, but its Mg–GDP oxygen contact
shortened from 2.157 Å to a sampled mean of 1.887 Å. That is an unresolved
physical-model issue to assess, despite the numerical checks passing.

The original 6JXT EGFR adduct template and corrected complete protein are now
assembled for research dynamics. Four earlier preserved attempts, including an experimental
two-neighbor relaxation and a restart from the pre-relaxation candidate, failed.
Fragment construction uses expanded observed
context but exports only the missing residues, leaving incompatible fixed
anchors in some cases. The geometry failures remain rejected; do not loosen
the angle/planarity checks or repeatedly rerun random candidates to claim success.
The report-only native context audit in `6jxt-context-v1` found backbone
displacements up to 2.101 Å and a serine sidechain displacement of 6.084 Å
before the current adapter discards those context coordinates. This explains
the incompatible starting joins; it does not authorize those changes or prove
the native candidate is accurate. The bounded coherent-context correction and
retained-environment screen are documented in [the loop follow-up](HANDOFF-AND-LOOPS.md);
the remaining native warm-up distortion is being diagnosed separately. The
unmodified original 6JXT and all failed models remain available.

## General defects identified and addressed in the research path

- SQM prints individual Mulliken charges to three decimal places in this
  bundle. BCC preserved those rounded sums, which initially failed an overly
  strict intermediate sum check. Admission now verifies completion, declared
  molecular total, printed precision, atom inventory, conservation through
  BCC, and unchanged native output. It does not modify charges. The final
  native PREPGEN total remains strictly checked. Six regression tests pass,
  including wrong-state, altered-charge and incomplete-calculation rejection.
- EGFR's global ligand dictionary describes the pre-reaction acrylamide.
  The working product graph explicitly applies the bond-order change only
  after confirming the deposited cysteine linkage and the supported graph.
  BTK's dictionary already describes the reacted bond and is not changed again.
  EGFR's distal aliphatic amine is an explicit +1 state hypothesis; a bound-state
  pKa has not been calculated.
- Native solvation centers the coordinates and reorders water/counterions.
  Retained atoms are matched by coordinates, identity and common translation,
  rather than assumed index stability.
- Context parameter changes do not update serialized System defaults. The
  original research runner saved a nonzero warm-up-restraint default and omitted
  parameters from State XML. All three completed runs now have independently
  reloaded `portable-restart-v2` exports with zero restraint and explicit Context
  parameters; Reference energies and forces are unchanged. Use these exports
  for a portable physical restart. The original binary checkpoints and raw
  exports are preserved; exact RNG continuation is a different operation.
- Native counterion placement replaced a deposited KRAS water in the first
  attempt. The corrected path chooses only newly added bulk waters, at least
  5 angstrom from retained heavy atoms and other selected ions under periodic
  distances. The second KRAS assembly preserves every original retained atom.
- The published GDP model is isolated into its own atom-type namespace so
  its Amber99-based source parameters cannot override protein ff14SB terms.
  General improper-type sorting changed one improper's ordering. Preserving
  native ordering and native phase conventions restored energy/force agreement
  on three coordinate probes to roughly 1e-12 in kJ/mol and kJ/mol/nm.
  GDP extracted from the complete KRAS topology also retains those mechanics.
  Failed export/format attempts are preserved. GDP parameter transport is
  validated here; GDP/Mg site physics and redistribution permission remain
  separate issues.

## Files and continuation

Native GROMACS peptide comparisons now pass for BTK, KRAS and the 6JX0 EGFR
control on three poses each. The original 6JXT pose exposes default export
rounding and a separate mixed-precision/high-force discrepancy. The preserved
precision diagnostic fixes the Reference-reader mismatch without changing
source parameters. See [cross-engine evidence](CROSS-ENGINE.md).
These capped-peptide checks are distinct from the completed full-complex
Reference-reader comparisons and native 20 ps continuations for corrected KRAS,
BTK and 6JX0. The 6JXT full-complex handoff remains unqualified.

The short KRAS Mg–GDP contact is not caused by a missing magnesium parameter.
The actual loaded CM/TIP3P values were verified and a published microMg source
was pinned for assessment. Its RNA/ion pair overrides require a compatibility
review before use with the current GDP/protein model. See
[Mg–GDP assessment](MG-GDP-ASSESSMENT.md). No replacement was accepted.

Actual outputs are under
`/Users/ashujo/.cache/dynamol-research/covalent-v1-focused/`:

- `6oim-prepgen-v2`, `5p9j-prepgen-v2`, `6jxt-prepgen-v2` — native residue/peptide candidates.
- `6oim-protein-v1`, `5p9j-protein-v1` — successful protein intermediates.
- `6jxt-protein-v1`, `v2`, `v3`, `v4` — rejected loop candidates/diagnostics.
- `gdp-scoped-v6` — GDP model with corrected native element metadata; v5 remains superseded evidence.
- `5p9j-complex-v2`, `6oim-complex-v3`, `6jx0-complex-v1` — complete solvated models used in the current completed numerical tests.
- `5p9j-stability-v1`, `6oim-stability-v2`, `6jx0-stability-v1` — completed 100 ps numerical tests, independently audited saved trajectories and portable restarts.
- `6jxt-protein-v5`, `6jxt-complex-v1` — screened research intermediate and complete model with all retained source molecules.
- `6jxt-stability-v1`, `6jxt-stability-v2` — rejected qualification attempts, unchanged.
- `6jxt-stability-diagnostic-v3` — completed 100 ps diagnostic; recorded loop failures and verified restart remain unqualified.
- `6jxt-stability-diagnostic-v4` — completed, audited 100 ps comparison with restraints released before heating; 13 flagged frames remain unqualified.
- `6jxt-warmup-comparison-v1` — source-bound comparison of all 18 monitored peptides and identical time windows; no default change accepted.
- `6jxt-target-native-terms-comparison-v1` — exact native local parameter comparison against 6JX0, including charges, proper and improper torsions.
- `6jxt-native-relaxation-v1` — separate completed native minimization experiment, not promoted or used for v4 dynamics.

Use `/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python`.
Do not restart the previous 94-atom HF optimization or held metal-site
optimization. Further attachment-parameter and metal-site assessment,
6JXT loop/warm-up validation, its full-complex cross-engine checks, broader frozen
panels and app lifecycle integration remain.

BTK's process loaded the runner before the final-PDB box export fix
and magnesium-monitor additions. Its saved state XML has the evolving box;
`final-with-current-box.pdb` was regenerated from that state and verified after
completion. The original export is preserved. KRAS used the corrected exporter
and has `launch.json` provenance.
The completed EGFR run loaded the runner before the Context-parameter export
fix. Its saved trajectory and corrected `portable-restart-v2` were independently
verified after completion; the Reference reload preserves energies and forces
exactly with zero warm-up restraint. Original exports remain preserved.

## Additional EGFR control

[6JX0](https://www.rcsb.org/structure/6JX0), from the same primary study as 6JXT,
is a 2.53 Å human EGFR T790M–osimertinib complex. Its deposited CIF explicitly
links Cys797 SG to YY3 C9 and has no internal sequence gaps. Its protein
intermediate passed. Native joint charge calculation and peptide assembly
completed under `6jx0-baseline-v1` and `6jx0-prepgen-v1`. The complete selected
model in `6jx0-complex-v1` has 53,791 atoms; its 100-ps test completed in
`6jx0-stability-v1`, including all 100 saved frames and a verified portable
restart. The control isolates the known
adduct workflow while 6JXT remains a rejected loop-repair case, and does not
replace or remove a frozen-panel case. The original structure contains three
drug copies. A focused single-adduct control explicitly selects the linked
copy and records exclusion of the two noncovalently associated copies while
retaining ions/waters and the full original files.
