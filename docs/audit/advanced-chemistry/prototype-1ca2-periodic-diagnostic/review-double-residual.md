# Independent review of the pinned double-precision discriminator

The large mixed-build/thread energy discrepancy is absent in the new double-precision build. With the same TPR, exactly identical atom coordinates and box, the one-versus-two-thread energy difference falls from **41.71875 to 0.000000256 kJ/mol**, and its force RMS vector difference is 3.47 × 10⁻¹³ kJ/mol/nm. This strongly supports a numerical execution/build origin for the previous large offset. The new compiler also differs (AppleClang16 versus conda Clang19), so precision alone is not isolated as the sole build variable. Actual TRRs contain eight-byte reals; the independent parser/synthetic precision check is retained separately.

The two-thread double result differs from OpenMM Reference by **−0.408748287 kJ/mol**. Its force RMS/max vector differences are **0.001243445 / 0.007062159 kJ/mol/nm**. No native output or force-field parameter has been altered.

| Recorded component accounting | kJ/mol |
|---|---:|
| Raw GROMACS double − OpenMM Reference potential | −0.408748287 |
| GROMACS double dispersion correction | −2395.302924765 |
| OpenMM recorded correction, same-system on/off difference | −2394.845913265 |
| Difference between recorded correction components | −0.457011501 |
| Remaining budget after identifying that component difference | +0.048263213 |

This is **component accounting, not an applied offset or a normalized pass criterion**. The residual remains unassigned. The existing component decomposition was recorded on OpenMM CPU, while the comparison total is Reference; the same analytical volume-dependent correction can be identified, but other CPU components must not be substituted for exact Reference components. No reciprocal-only explanation follows from these numbers.

GROMACS explicitly removes excluded pairs from its averaged dispersion coefficients ([pinned source](https://github.com/gromacs/gromacs/blob/v2025.4/src/gromacs/mdlib/dispersioncorrection.cpp#L160), [manual](https://manual.gromacs.org/documentation/2025.4/reference-manual/functions/long-range-vdw.html)). Exact tail averaging, whether repulsive C12 tails are included, switching functions, cutoff shifts, box volume and exception conventions must therefore be declared when comparing engines. A volume-dependent correction also affects pressure; treating the observed value as an arbitrary constant across boxes would be wrong. This review has not derived every fraction of the tail difference or the remaining +0.048263213 kJ/mol from first principles.

The reusable validation contract should keep separate results for mapped physical parameter integrity, raw same-coordinate energy and forces, explicit analytical conventions, and numerical precision/repeatability. Close forces or disappearance of the thread effect does not establish native site accuracy, realistic metal dynamics, or general acceptance for the full metal panel. Existing full-case motion results and chemistry boundaries remain separate.

`review-double-residual.json` pins the completed native comparison and source code and reproduces this arithmetic. The native artifacts are in `build/advanced-chemistry/gromacs-double-v1/`; no further native calculation was launched for this review.
