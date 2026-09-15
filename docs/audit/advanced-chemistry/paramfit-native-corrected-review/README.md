# Corrected native Paramfit qualification review

Both previously executed native cases returned zero and produced valid **four-column**, finite energy tables with the requested **500 / 1,000 rows**. Their current upstream fixture assertions are satisfied within their stated scope. The original runner marked them failed because it expected three columns; that failed analysis and every native input/output remain unchanged. **No native calculation was rerun.**

This qualifies these small synthetic program fixtures. It does not validate a ligand, covalent adduct, metal site, QM target dataset or generally robust parameter fitting.

| Saved native case | Rows / exact IDs | Fitted−target RMSE, kcal/mol | Maximum absolute error, kcal/mol | Current upstream comparison |
|---|---|---:|---:|---|
| Dihedral least squares, NMA | 500 / 0–499 | 0.000242033 | 0.000420 | Target-column reference exact; generated `frcmod` byte-identical to current saved reference |
| Simplex, NMA | 1,000 / 0–999 | 0.000058798 | 0.000180 | Fitted-energy max difference 0.000100; native `ndiff` max relative error 3.95667e−6, below upstream 5e−6 |

All four columns are finite, row IDs are unique and sequential, and the target column equals the sorted first requested entries of the corresponding 20,000-entry input target file exactly. Sorting is requested by the original controls. Residuals above use **Amber+K minus Quantum**, without centering or fitting an offset. Native energies are printed to five decimal places, so these are residuals of the saved decimal output.

## Which upstream references apply

The least-squares `run_test.x` selects tab field 3: the **target** energy column, not the fitted energy. Its 500-row `energy.out.saved` matches exactly. This assertion alone checks target reading/sorting; the byte-identical `frcmod.saved` supplies the fitted-parameter regression evidence. The nearby three-column `energy.dat.saved` is an older **100-row** artifact and is not silently truncated or promoted to the current reference. Its historical log also describes a different fit.

The simplex `run_test.x` selects the fitted energy column and specifies a 5e−6 relative limit. The calculation here follows the actual installed `ndiff.awk` rule, `abs(a−b)/min(abs(a),abs(b))` for these nonzero values. The current control reads `prms.in`, selecting **six torsion amplitudes**. It requests no `frcmod`; the final parameters are preserved in `native.log`, so absence of that file is not a native failure. The historical saved log instead describes a **19-dimensional** selection including other parameters. Its printed parameters are compared diagnostically, but cannot establish exact current-configuration parameter parity.

The review records 30 finite printed parameter terms per case: seven bonds, ten angles and thirteen torsion terms. The largest initial-to-final printed torsion-amplitude changes are 0.0004 kcal/mol for least squares and 0.0019 kcal/mol for simplex. Full identities, initial/final numbers and historical differences are in the separate parameter JSON files. Wrapped phase differences are reported alongside raw differences; no parameter file is changed.

## Fourth-column source issue

The native header reads `Num, Amber+K, Quantum, Initial Amber+K`. However, `write_input.c::write_energy` writes `init_energy+K`. In the sum-of-squares implementation, `eval_amber_std.c` initializes the accumulator with `K−target`, adds Amber terms, and stores that **residual** in `init_energy`. With this fixture's K=0, the fourth column therefore is not an initial Amber energy curve. Its values range from −0.00017 to +0.00018 kcal/mol. The least-squares path records zero fitness calls and its fourth column is all zero.

The corrected review preserves the original column and validates its finiteness, but does not use it as fitted energy or as an initial-energy baseline. The source files and hashes identifying this behavior are recorded. A future adapter should compute an actual initial energy explicitly if that curve is required.

## Warnings and limits

The least-squares log has **16 warning lines**: six torsion-coverage warnings, four angle-coverage warnings, two bond-coverage warnings, two aggregate insufficient-data warnings, and two notifications that all torsion amplitudes are fitted while phases are fixed. The angle and bond warnings remain despite successful completion; `CHECK_BOUNDS=WARN` was explicitly requested by the unmodified upstream control. The simplex log has one warning: seven missing sampling bins for an HC–CT–C–O torsion. Its optimizer reports 150 cycles, 269 evaluations at convergence and 270 total fitness calls, with convergence ratio 2.0978e−6 against the requested 0.1 criterion.

The fixture descriptions state that targets were generated from **Sander/parm99 energies**, with increased printed precision and a corrected trajectory/energy indexing offset. These are synthetic self-consistency targets, not quantum calculations. They are training samples, not an independent held-out dataset. The source also deliberately starts near the desired parameters: the simplex initial sum of squares was already 3.4784e−6 kcal²/mol². Its small final residual does not demonstrate recovery from a poor starting model or exploration of every relevant torsion. Neither fixture tests fitted forces, chemical-state assignment, broad ligand coverage or protein motion.

The original runner omitted the simplex selection file `prms.in` from its pre-run hash inventory. The current copy matches the pinned source and the native six-dimensional selection, but this does not fabricate the missing historical pre-run hash. Future qualification manifests should include every parameter-selection file.

## Evidence preservation

`review.json` contains all measurements, warnings, source references, input/output hashes and original errors. Original native data are under `~/.cache/dynamol-native-audits/paramfit-v1/`; the untouched failed report remains `build/advanced-chemistry/paramfit-native-audit-v1/result.json`. The review verifies those original hashes and the native executable hash. The separate script is `docs/audit/advanced-chemistry/review_native_paramfit.py`; it invokes no external executable. An initial review-parser attempt that mistook an explanatory parenthetical line for a parameter row is retained in the sibling `paramfit-native-corrected-review-attempt1` directory.
