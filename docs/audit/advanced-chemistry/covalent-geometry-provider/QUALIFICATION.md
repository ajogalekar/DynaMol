# GFN2-xTB calculator: small-molecule numerical qualification

The isolated energy/gradient interface passed the tests below. This qualifies its numerical plumbing on small neutral closed-shell examples; it does **not** validate a covalent adduct, attachment torsion, metal state, optimized geometry, or force field. No large adduct, optimization, scan, charge replacement, or biological target was calculated in this qualification.

The reusable module is `calculator.py`. The installed runtime is `/Users/ashujo/.cache/dynamol-runtimes/covalent-geometry-v1`. Current accepted evidence is `build/advanced-chemistry/covalent/geometry-provider-qualification-v2/result.json`; the first failed attempt remains in `geometry-provider-qualification-v1`.

## Method and units

The exact method is the built-in **GFN2-xTB** parameterization in tblite 0.7.0. It includes GFN2's own self-consistent D4 dispersion; no additional D3 correction, solvent, restraint, or penalty potential is added. This is a semiempirical electronic method, not the independently qualified B3LYP-D3(BJ)/DZVP reference. The [primary GFN2 paper](https://doi.org/10.1021/acs.jctc.8b01176) and [tblite Python API](https://tblite.readthedocs.io/en/latest/api/python.html) describe the native method and interface.

Coordinates are explicitly **Bohr**, native energy is **Hartree**, and the returned gradient is **Hartree/Bohr**, with positive derivative dE/dx. A force requires the negative of this gradient. Units follow the pinned native interface and are checked by finite differences. Absolute GFN2 and ab initio energies have different model conventions; their raw totals are not interchangeable reference origins.

The default SCC accuracy is **0.001**, maximum iterations **250**, and electronic occupation temperature **0.00095 Hartree**. Temperature is not specified in kelvin. Fresh calculations are used for each evaluation; no density or charge restart is silently reused. Native SCC nonconvergence raises an error and is not accepted as a last-iterate result.

## Numerical results

All tolerances were fixed before calculation. The original upstream regression limits are 10⁻⁶ Hartree and 10⁻⁶ Hartree/Bohr. The derivative limit is 2×10⁻⁶ Hartree/Bohr, using central steps 2×10⁻⁴ and 1×10⁻⁴ Bohr. Permutation/rigid-transformation limits are 10⁻⁸ Hartree and 10⁻⁷ Hartree/Bohr.

| Check | Observed result |
|---|---|
| Exact upstream 22-atom alanine-dipeptide restart sequence, two conformers | Maximum energy difference 1.43×10⁻¹⁴ Hartree |
| All 66 stored upstream gradient components | Maximum difference 4.88×10⁻¹¹ Hartree/Bohr |
| Water, ammonia, H₂S, fluoromethane, methanethiol | All 63 gradient components checked at two finite-difference steps |
| Worst derivative difference, across both steps/all five molecules | 1.40×10⁻⁸ Hartree/Bohr |
| Atom permutation | Maximum energy difference 6.22×10⁻¹⁵ Hartree; gradient difference 9.81×10⁻¹⁶ Hartree/Bohr |
| Rotation plus translation | Maximum energy difference 7.11×10⁻¹⁵ Hartree; gradient difference 2.13×10⁻⁹ Hartree/Bohr |
| Input/state/method/unit/resource errors | 15 invalid requests rejected |
| Deliberately insufficient one-cycle SCC | Native `TBLiteRuntimeError` rejected |
| Missing enclosing execution guard | Rejected |
| Shared lock held by another process | Actual second process rejected before a native calculation |

The accepted run attempted 272 tiny native evaluations and took **1.73 seconds inside the runner** (1.91 seconds including its outer process). These timings do not establish performance for a 94-atom optimization. Base inputs/results, analytic and finite-difference arrays, deterministic perturbation-generation code, summary residuals, source/runtime hashes and failure outputs are retained. Individual displaced energies and transformed/permuted input/result arrays were not separately saved.

The upstream fixtures come from the pinned [tblite v0.7.0 test source](https://raw.githubusercontent.com/tblite/tblite/v0.7.0/python/tblite/test_interface.py). The finite-difference checks independently differentiate the returned energy, but they still use the same electronic implementation. This is not agreement between independent quantum engines or experimental validation.

## Preserved first failure and SCC interpretation

The initial test evaluated the second upstream conformer from a fresh default-accuracy calculation. The upstream gradient actually comes from updating the first conformer and reusing its result. The fresh default calculation agreed in energy but differed by **5.17×10⁻⁶ Hartree/Bohr** in a gradient component, failing the original 10⁻⁶ threshold. That failed report, traceback and executed script are preserved.

The corrected regression follows the upstream restart sequence exactly and passes without widening the threshold. The adapter's tighter fresh SCC calculation is checked separately. Its gradient differs by up to **4.60×10⁻⁶ Hartree/Bohr** from the loose restarted upstream result; this numerical convergence difference is reported, not treated as a matching reference. The five complete small-molecule finite-difference tests validate the adapter's tighter default directly.

## Runtime provenance and resource contract

The downloaded PyPI wheel is `tblite-0.7.0-cp312-cp312-macosx_15_0_arm64.whl`, SHA256 **29f4838f63bd22da973462c8912f2dfa6d260d39e1db2eb4bda455ec562ccfe4**. All **18 tblite package members** match the installed files byte for byte. The wheel, PyPI metadata, pinned upstream interface/test files, and their license notices are retained. tblite's upstream source carries LGPL-3.0-or-later terms; runtime/dependency redistribution needs its own packaging review.

`runtime-manifest.json` pins **1,590 files / 129,841,031 bytes**, covering the Python executable and relevant tblite/geomeTRIC/numpy/scipy Python and native-library files. Package versions are tblite 0.7.0, geomeTRIC 1.1.1, numpy 2.5.3, scipy 1.18.1 and pytest 9.0.2. The installed runtime was not modified during qualification.

`execution_guard(root)` verifies these hashes, holds `build/advanced-chemistry/QM-LAUNCH.lock` for the enclosing operation, refuses a fourth recognized quantum worker, and requires 2 GiB free for this small-system qualification. Native tblite OpenMP and OpenBLAS thread counts were explicitly set and observed as **one**. Existing running quantum processes and their runtimes were untouched.

The guard deliberately does not provide a native wall-clock interrupt. The caller must supply an outer process-group timeout. Qualification used bounded outer controllers and finished in under three seconds across both attempts; a future optimizer must enforce its own reviewed limit. The whole optimization should remain inside one guard so that each callback cannot race another launcher.

## Callable contract and scope

```python
with calculator.execution_guard(project_root):
    result = calculator.evaluate({
        "schema_version": 1,
        "method": "GFN2-xTB",
        "atom_ids": atom_ids,
        "elements": elements,
        "charge": 0,
        "spin": 0,
        "coords_bohr": coordinates,
    })
```

Inputs require unique stable atom IDs, explicit H/C/N/O/F/S elements, an integer charge and closed-shell spin zero, finite coordinates, and a positive even electron count. Unknown settings are refused. The current limit is **32 atoms**. The positive numerical fixtures are all neutral; support for explicit nonzero integer charge in the interface is not evidence that charged chemistry has been qualified. Outputs preserve atom order and state and expose energy/gradient, exact numerical settings, versions and hashes. No atomic charges are exported or substituted into Amber.

The module is only a calculator. It does not build structures, check a reaction graph, enforce cap freezes/dihedral constraints, determine protonation, or validate stereochemistry. A separate admitted parent workflow and independently checked optimizer constraints are required before any larger covalent use. The planned cheaper-geometry/DFT-energy strategy still needs a matched constrained DFT comparison and held-out conformational evidence.

## Execution-guard follow-up

A small optimizer test exposed a conservative process-counting error: an idle `test_geometric_adapter.py` supervisor was counted as native work because of a substring match. The guard now classifies the executed Python script basename, and counts `geometric_adapter.py` only in `--worker` mode. Twelve no-compute classification tests cover actual workers and idle supervisors. The `validate`, `evaluate`, and `verify_runtime` implementations are unchanged; the original numerical source snapshot is retained, and `guard-classification-review.json` links its hash to the corrected guard. No scientific fixture was rerun or reclassified to conceal the first failed launch.
