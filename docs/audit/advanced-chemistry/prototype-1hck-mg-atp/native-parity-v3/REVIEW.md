# 1HCK native static parity: residual explained

The complete 5,216-atom Mg–ATP/protein/water model agrees between independent Sander and OpenMM Reference evaluations after matching two documented native numerical conventions **in diagnostic copies only**. The original raw differences and failed comparisons remain preserved. This does not establish chemical transferability, equilibration, or production application support.

## Native lifecycle failure

Each of three poses now runs in a fresh native Python/Sander subprocess, with a 120-second per-pose timeout and separate immutable request/result/log files. All three complete successfully. The previous repeated-setup failure remains consistent with native CMAP allocation state persisting across setup/cleanup; the pinned source and failure are retained. No installed engine was edited.

The general helper rejects reused output paths and wrong-pose native results. Three focused regression tests pass; the actual three native runs supply additional integration evidence. The complete original structure and two deterministic 0.003 Å perturbations are identical between evaluators. No dynamics or quantum calculations were run here.

## Numerical localization

Original tolerances remain **0.001 kJ/mol total energy** and **0.01 kJ/mol/nm maximum force-component difference**, over every atom. Raw comparisons fail. The first diagnostic copy matches Amber's source-backed legacy charge scale 18.2223 to OpenMM's independently measured electrostatic constant; both particle charge and 1–4 products are handled consistently. This is a known evaluator convention, not a charge fit or an energy offset chosen from this molecule.

A second diagnostic copy additionally implements Sander `dihpar`'s documented near-pi rule: phases within 0.001 radians of pi evaluate at exact mathematical pi. There are 6,225 affected terms, each stored 1.3464102069×10⁻⁶ radians above pi. Their coordinate-independent energy-error bound is 0.07175103 kJ/mol. Eight distinct fitted ATP phases outside that rule are retained exactly. The first diagnostic attempt stopped on these general phases; its script and error log are preserved. The corrected audit explicitly verifies that they trigger neither the near-pi rule nor native small-trigonometric-value clamping.

| Pose | Raw ΔE (kJ/mol) | Coulomb-matched ΔE | Coulomb + phase ΔE | Final max component ΔF (kJ/mol/nm) |
|---|---:|---:|---:|---:|
| Bound | −1.435244729 | +0.002559612 | +9.211×10⁻⁹ | 2.532×10⁻⁹ |
| Perturbation 2027 | −1.434939168 | +0.002550723 | +7.683×10⁻⁹ | 2.094×10⁻⁹ |
| Perturbation 2028 | −1.435120836 | +0.002556246 | +5.472×10⁻⁹ | 6.446×10⁻⁹ |

Differences are OpenMM minus Sander. Raw maximum component force discrepancies reach 0.12869042 kJ/mol/nm; Coulomb matching reduces them below 0.003804 before phase matching. With both conventions matched, all three cases pass the original tolerances by a wide margin. No constant subtraction is used.

Per-component comparisons separately cover all bonds, angles, periodic torsions, CMAP, pair-specific Lennard-Jones interactions, electrostatics, and 1–4 terms. At the bound pose the CMAP energy difference is 2.36×10⁻¹¹ kJ/mol and custom LJ difference −1.98×10⁻⁹ kJ/mol. Component energies and force vectors sum back to the full OpenMM totals. Every raw and diagnostic force component, and every native total force vector, is retained in `full-force-components.npz`; native energies/components/forces remain in the individual JSON results. The earlier 576-value CMAP and 361-entry pair-table transport audit remains independent evidence.

## What this does and does not resolve

This localizes the earlier −1.435 kJ/mol discrepancy to general native electrostatic and torsion-phase evaluation conventions. It supplies evidence that the full CMAP and custom pair model is transported numerically, without discarding special terms. It does not make the raw evaluators identical or silently change the app's force field.

The original published global Mg–oxygen corrections also affect protein oxygen classes. Their transferability to this protein-bound site still needs scientific justification; ATP/Mg aqueous benchmarks alone do not establish it. The SI's CC BY-NC 4.0 license remains a distribution limitation. The crystal/template coordinates remain unrelaxed and have high forces, so this is not a stability result. Application support remains disabled.

## Evidence

`component-convention-diagnostic.json` records detailed components, original limits, source/hash evidence, force statistics, and unchanged model hashes. `sander-parity.json` preserves the failed raw/Coulomb-only comparison. `provenance.json` pins the native runtime extension, helper, diagnostic, tests, and generated artifacts. Prior status and review snapshots, the old helper, the failed first diagnostic, and the original repeated-native-setup failure are all retained.
