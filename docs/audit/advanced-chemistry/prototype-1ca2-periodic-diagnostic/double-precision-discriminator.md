# Native double-precision discriminator

The large thread-dependent energy discrepancy collapses in an isolated GROMACS 2025.4 double-precision build. This is evidence of a native numerical-accumulation issue in the earlier comparison, **not acceptance of exact periodic energy parity or validation of the metal chemistry**.

## Fixed inputs and scope

Only two native `nsteps=0` evaluations were run, with one and two CPU threads. Both read the original `matched-mesh/diagnostic.tpr`, SHA256 `5aa9285424526503ae95aebc3c786baf5c318ea86fbbdac3a12fe276b2ac3ac2`. Its 34,414 atoms, coordinates, box, force-field parameters, PME coefficient/order/mesh and cutoffs are unchanged. Both native logs explicitly report zero steps and a single-precision input TPR. The new output TRRs store eight-byte reals; the independently checked audit reader preserves that precision.

All coordinates and box entries match the original mixed-precision frame **bitwise**. The one- and two-thread coordinate arrays also match exactly. There was no new MD, replacement of an installed engine, topology rewrite, parameter edit or applied energy offset.

## Measured results

OpenMM Reference potential energy at this geometry is −429191.2916377305 kJ/mol.

| Native runtime | Threads | Potential energy, kJ/mol | Difference from OpenMM Reference, kJ/mol |
|---|---:|---:|---:|
| Existing mixed precision | 2 | −429166.656250 | +24.635388 |
| Existing mixed precision | 1 | −429208.375000 | −17.083362 |
| Isolated double precision | 2 | −429191.700386018 | −0.408748287 |
| Isolated double precision | 1 | −429191.700386274 | −0.408748543 |

The native one-minus-two-thread energy change shrinks from **−41.71875** to **−2.5594e−7 kJ/mol**. The double-runtime full-force RMS/max vector differences from OpenMM Reference are **0.00124345 / 0.00706216 kJ/mol/nm**. The one-versus-two-thread force RMS vector difference is **3.4682e−13 kJ/mol/nm**. All full-system force values are finite. The two static evaluations and their extraction took 3.65 seconds after compilation.

The double-runtime short-range Coulomb and reciprocal energies are −565353.718491989 and +3962.082392056 kJ/mol. The previous independent source/component audit placed the large mixed-precision offset mainly in direct-space Coulomb/exclusion/self-energy accumulation, rather than the reciprocal mesh. These new observations are consistent with that localization; they do not identify a single exact floating-point operation as the sole cause.

## Smaller residual kept explicit

The double-native dispersion correction is −2395.302924765479 kJ/mol. The preserved OpenMM calculation with and without its dispersion correction gives −2394.845913264900 kJ/mol. Their **−0.457011501 kJ/mol** difference is a separately measured energy-component convention difference. After identifying that component, **+0.048263213 kJ/mol** remains in the raw total comparison. This is explanatory bookkeeping, not an offset applied to either engine, and the remaining difference has not been fully apportioned among serialized parameter precision and other implementation conventions.

The unchanged input TPR was written in single precision. Running it with double arithmetic cannot recover parameter precision already lost during serialization. The raw −0.408748287 kJ/mol total difference therefore remains visible, and the earlier strict absolute gate is not declared passed by changing its threshold.

## Build provenance and limits

The build uses the preserved official 2025.4 archive with SHA256 `ca17720b4a260eb73649211e9f6a940ee7543452129844213c3accb0a927a5c3`, ARM64 CPU SIMD, OpenMP, and double FFTW 3.3.11. Native source files were not edited. Apple Clang 16 differs from the bundled conda Clang 19, and unused optional muparser, Colvars and lmfit components were disabled through supported CMake options. **Precision was therefore not literally the sole build variable.** The large reduction strongly supports a numerical runtime cause; it is not a controlled compiler-independent proof of one arithmetic instruction's behavior.

Earlier optional-component build failures are preserved. In particular, the failed lmfit consumer dependency file shows an unrelated `/usr/local/include/lmmin.h` shadowing the bundled header. The final successful build took 936 seconds at two compiler jobs; all source/build files remain in the isolated cache. The installed GROMACS executable's SHA256 is unchanged.

The machine-readable result is `double-precision-discriminator.json`. Exact native logs, eight-byte TRRs, energy files, force arrays, executable/core-library hashes, compiler flags, linked-library provenance, failed-build evidence and commands are retained under `build/advanced-chemistry/gromacs-double-v1/`. The helper is `docs/audit/advanced-chemistry/check_1ca2_gromacs_double.py`.

No broad metal-model, equilibrium, production-stability or application-release claim follows from this single-endpoint numerical diagnostic.
