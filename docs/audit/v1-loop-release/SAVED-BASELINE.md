# Saved loop baseline: native Ramachandran audit

**The old geometry-only acceptance is insufficient for the v1 loop release.** Five of ten originally accepted panel preparations contain modeled-residue Ramachandran outliers. Including fixed boundary residues whose torsions become newly defined by the reconstructed loop, only **2 of 10** pass this additional reference check. Both extra accepted 8K5R repeats fail it. No original result was changed or relabeled, and no app admission was changed.

This read-only audit evaluates the frozen panel’s 90 modeled residues, plus the recorded observed context that actually moved in the repeat. All 12 coordinate audits completed without mapping or artifact-hash errors. Native CCTBX 2025.11 scores and classifications are used directly; all selected residues were scorable, and selected classifications matched between full-precision NPZ coordinates and the saved PDB files.

## Results

| Saved case | Modeled residues | Moved observed residues | Entry selected outliers | Final selected outliers | Additional boundary outliers | Unchanged observed outliers |
|---|---:|---:|---:|---:|---:|---:|
| 8K5R | 17 | 0 | 0 | 5 | 1 | 5 |
| 1UA2 | 12 | 0 | 1 | 3 | 1 | 18 |
| 1HCK | 4 | 0 | 0 | 0 | 1 | 3 |
| 4HJO | 4 | 0 | 0 | 0 | 1 | 0 |
| 1M17 | 12 | 0 | 2 | 1 | 0 | 7 |
| 2HYY | 4 | 0 | 0 | 1 | 2 | 1 |
| 3HEG | 15 | 0 | 1 | 2 | 0 | 10 |
| 2CG9 | 9 | 0 | 0 | 0 | 2 | 8 |
| 1FPU | 5 | 0 | 0 | 0 | 0 | 6 |
| 4HHY | 8 | 0 | 0 | 0 | 0 | 5 |
| 8K5R-repeat-42 | 17 | 4 | 0 | 8 | 2 | 4 |
| 8K5R-repeat-2026 | 17 | 0 | 0 | 3 | 1 | 5 |

“Selected” means requested modeled residues plus observed residues whose heavy atoms actually moved during loop refinement. “Boundary” identifies a fixed observed residue whose phi/psi includes a newly modeled or changed atom. A boundary failure is not counted as a modeled residue failure, but it is relevant to the finished local model. The unaffected observed outliers are reported separately and do not count against either local criterion. For all ten main cases, the phi/psi-defining atoms of those unaffected observed outliers are present and unchanged in the original imported selected structure.

## Findings that guide the fix

- **8K5R:** all 17 modeled residues were non-outliers at refinement entry; five become outliers after fixed-heavy refinement: GLY28, TYR92, ARG94, CYS95 and LYS178 in chain A. Fixed boundary PRO91 also fails. The smallest affected gap is A:177–179 (AKN), between observed LEU176 and SER180.
- **1UA2:** LYS44, GLY46 and ARG48 become outliers although all three were favored at refinement entry. Its 12-residue loop therefore needs the stronger conformational check even though the earlier bond, angle, chirality and gross-contact checks passed.
- **1M17:** the other 12-residue example also has one modeled outlier. A length limit alone cannot establish an acceptable conformation.
- **1HCK, 4HJO and 2CG9:** modeled residues pass, but one or more newly defined fixed boundary torsions fail. A modeled-residues-only gate would miss these.
- **8K5R repeat seed 42:** the declared movable flanks match measured motion exactly; SER180 and GLU260 are native-reference outliers, alongside six modeled outliers. The motion is within the existing 1 Å cap, showing that a displacement cap alone is insufficient.

## Scope and provenance

The five originally blocked panel entries (6A93, 2BEL, 2REN, 1RNE and 2ITY) remain in the frozen 15-entry panel and were not reparsed as successful models. They have no accepted coordinate output eligible for this baseline audit. 1HCK is included for its saved loop geometry; this does not validate its Mg–ATP interaction model or change the separate complete-complex chemistry requirements.

Unchanged refers to the saved refinement-entry coordinates for the two repeat jobs because their copied job artifacts do not provide a verified selected-source dataset link. This limitation is recorded explicitly; original source coordinates were not inferred. The 10 main cases do have that verified source comparison.

Ramachandran classifications describe a static structural reference distribution. They do not recover the experimental missing coordinates, establish native loop accuracy, validate force-field dynamics, or provide an MD per-frame acceptance rule. No new minimization, force-field calculation, MD, QM, backend edit or publication was performed.

- [Audit implementation](audit_saved_loop_baseline.py)
- [Machine-readable compact results](saved-baseline-summary.json)
- [Full results](/Users/ashujo/.cache/dynamol-research/v1-loop-release/saved-baseline-v1/result.json)
- [Frozen plan and hashes](/Users/ashujo/.cache/dynamol-research/v1-loop-release/saved-baseline-v1/plan.json)
- Full result SHA-256: `658186a3260f162a620573fb10a3c78c55e07f506be2dff6a1f5848d23c549f0`.
- Frozen 15-structure manifest SHA-256: `a0ebcf95c9f1ea27299ea8fe4c717e55d6fdab596c566c2869f5e775c9b2c587`.
- Native reference extension SHA-256: `b7a531377948dc6c1534ffe2f4041cfe28c38eac6817e12e0808e828da28935c`.

Execution used the isolated research Python with the separate CCTBX target and at most two numerical-library threads. The source program and reference implementation are copied into the new cache run; original artifact hashes were checked before and after scoring. The audit script also passed a Python AST parse and a targeted whitespace check.

```bash
PYTHONPATH=/Users/ashujo/.cache/dynamol-runtimes/cctbx-rama-reference-v1 \
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python -u \
  docs/audit/v1-loop-release/audit_saved_loop_baseline.py \
  --output /Users/ashujo/.cache/dynamol-research/v1-loop-release/NEW-UNUSED-RUN
```
