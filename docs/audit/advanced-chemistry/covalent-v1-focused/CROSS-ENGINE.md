# Focused adduct transport checks, 14 September 2026

Actual GROMACS 2025.4 mixed-precision static energies and forces were compared
with OpenMM Reference reading the original Amber peptide and the exported
GROMACS topology. Each case uses the source pose plus two deterministic 0.003 Å
perturbations, with identical recorded coordinates in each comparison. These
are capped modified residues between peptide neighbors, not the full solvated
complexes. No dynamics or chemistry-accuracy claim follows from these checks.

| Case | Actual native GROMACS comparison | Double-precision topology transport |
| --- | --- | --- |
| 5P9J, BTK–ibrutinib | 3/3 pass | 3/3 pass |
| 6OIM, KRAS–sotorasib | 3/3 pass | 3/3 pass |
| 6JX0, EGFR–osimertinib control | 3/3 pass | 3/3 pass |
| 6JXT, original EGFR case | 0/3 pass; retained for diagnosis | Default export fails; precision diagnostic passes |

Native tolerances were fixed before the runs at 0.001 kJ/mol energy and
0.01 kJ/mol/nm maximum force-component difference. The separate Reference
reader comparison uses 0.0001 and 0.001 in those units. No thresholds were
changed to make the original EGFR case pass.

The 6JXT default ParmEd export rounds Lennard-Jones parameters to eight
significant digits and atomic charges to eight decimal places. Restoring those
fields directly from the unchanged native Amber values reduces its maximum
Reference-reader force discrepancy from 0.0130 to 0.0000123 kJ/mol/nm, and
energy discrepancy from 0.000166 to 0.000000500 kJ/mol across the three probes.
This is a demonstrated export-precision issue. It does not resolve the separate
actual mixed-precision GROMACS discrepancy (up to 0.204 kJ/mol/nm; relative RMS
errors around 1.5e-6 at this high-force source pose), which remains recorded and
needs a native precision/geometry assessment. The protein loop-repair rejection
also remains independent and unresolved.

The new research harness is `check_peptide_gromacs.py`. Native outputs, TRR
forces, source hashes and comparisons are under the focused cache:

- `5p9j-peptide-gromacs-v2`
- `6oim-peptide-gromacs-v1`
- `6jxt-peptide-gromacs-v1`
- `6jx0-peptide-gromacs-v1`
- `6jxt-export-precision-diagnostic-v1` (exact before/after parameter lines)

The first BTK invocation rejected an unused FLEXIBLE macro in this water-free
peptide. Removing the unused macro corrected the harness; no GROMACS warning
override was used. The rejected invocation remains in `5p9j-peptide-gromacs-v1`.
No source chemistry, original topology or application converter was changed.
