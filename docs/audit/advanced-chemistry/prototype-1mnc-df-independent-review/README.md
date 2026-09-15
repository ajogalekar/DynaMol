# 1MNC fixed-geometry density-fitting review

The saved conventional and density-fitted calculations use the same 97 ordered atoms, coordinates, charge +1, singlet state, resolved B3LYPG functional, spherical 6-31G* orbital basis, and level-4 grid. Both SCF calculations converged. Independent analysis reproduced the controller's differences exactly; no native calculation was rerun.

| Quantity, DF minus conventional | Saved result |
|---|---:|
| Energy, Hartree | −0.000487164787046 |
| Maximum gradient-component difference, Hartree/Bohr | 0.0000459547007168 |
| RMS gradient-component difference, Hartree/Bohr | 0.00000978536199290 |
| Maximum gradient-vector difference, Hartree/Bohr | 0.0000507321490533 |
| RMS gradient-vector difference, Hartree/Bohr | 0.0000169487441422 |

The DF auxiliary basis is `def2-universal-jkfit`. It started from the converged conventional density and then performed nine SCF cycles. Its 734.93 seconds versus the original 6,505.75 seconds is **not a controlled speedup measurement**: initial guesses, worker, and native runtime builds differ. Their declared library versions and resolved scientific method agree. The copied runtime declaration is hash checked; this review does not re-attest every loaded native binary.

This is one fixed, unoptimized geometry. The DF maximum absolute gradient component is 0.120622 Hartree/Bohr. The observed approximation differences do not establish relative conformer-energy accuracy, optimized-geometry agreement, Hessian accuracy, protonation-state correctness, or a validated site force field. Neither full 1MNC preparation nor dynamics is ready.

The source, array, checkpoint, worker, method, and request hashes are preserved in [review.json](review.json). Actual conventional and DF coordinates are bitwise identical in both Bohr and Angstrom. The original decimal input matches the native atom map exactly. The saved Angstrom values differ from those decimals by at most 7.1×10⁻¹⁵ Å, entirely reproduced by the exact PySCF unit-conversion operations. The first independent analysis stopped on an overly strict decimal-input identity assertion; the corrected analysis verifies that exact roundtrip, with no coordinate-relaxation tolerance. It did not alter either calculation.

The 12 native caps are explicitly mapped. In each of His218, His222 and His228, the original CA becomes a CH3 carbon at the same coordinates; N/HA/C become methyl hydrogens. The three carbon boundary identities are `cap:small:15:CH3`, `cap:small:54:CH3`, and `cap:small:92:CH3`. This matches native `write_sc` and `build_small_model`; the reviewed source is copied and hashed here.

The original native small-model optimization input requests unconstrained optimization. It supplies **no cap-freeze policy**. [optimization-readiness.json](optimization-readiness.json) records the exact boundary candidates and the unresolved modeling decisions, without selecting or launching a new request. Freezing the three carbon positions while allowing their hydrogens to relax would be a distinct constrained-cluster approximation. Freezing all twelve cap atoms imposes a different constraint and preserves unrelaxed cap-hydrogen guesses. Neither follows automatically from native input.

The explicit model hypothesis remains Zn(II), three NE2-coordinating HID residues, and a full PLH hydroxamate anion with O1 deprotonated and N1 protonated. Structural Zn, Ca, nonlocal protein, and solvent are omitted. Geometry optimization, checked Hessian extraction, larger-cluster ESP/RESP, other-site models, full-complex mapping, native force validation, and motion assessment remain separate unfinished stages.

Reproduction command: `.venv/bin/python -B docs/audit/advanced-chemistry/review_1mnc_density_fit.py`. The script is read-only with respect to source calculations and refuses to overwrite its review directory; inspect it before choosing a fresh review output directory. It invokes no native chemistry program.
