# Amber ↔ GROMACS periodic validation contract

**Decision:** validate parameter conversion and native execution as separate, mandatory evidence lanes. Preserve raw energies. A known convention difference may make an absolute-energy comparison inapplicable; it must never turn a failed comparison into a pass by fitting or subtracting a constant. A force match is useful evidence, but cannot substitute for a complete parameter inventory or for volume/pressure checks.

This is the general development contract for declared, supported classical models. It does not establish a metal oxidation state, protonation state, chemical model accuracy or catalytic mechanism. It introduces **no production code changes, replacement engines or new numerical tolerances**.

## 1. Native parameter and identity fidelity

Before comparing numbers, require a full-system manifest derived independently from the original native parameter artifact. Bind it to the exact source topology, exported topology, imported System, native executable/core libraries, importer/validator version and atom-order map. Preserve chemical identifiers and record native template aliases explicitly.

Require every atom, element, mass, charge, pair interaction, bond, angle, ordered torsion term, exclusion and 1–4 interaction. Do not reconstruct the expected values from the System being tested. Bond graphs and physical bond terms are separate requirements. Record virtual sites and constraints as explicit protocol inventories; a restraint or undeclared constraint cannot stand in for a native bond. Unknown CustomForce, CMAP, NBFIX-like pair overrides or other unsupported terms prevent full-system acceptance until independently covered. A passing site-only report remains partial.

Narrow mathematical no-ops already supported by the validator remain narrow: an exactly zero torsion amplitude may be omitted only with valid positive periodicity, while its independent 1–4 graph must still match; sigma is irrelevant only when both compared epsilon values are exactly zero. These rules do not permit removal of small nonzero terms. Every validator change requires new evidence bound to its new hash; old reports do not silently become current.

## 2. Source-to-export functional check

Evaluate both the original Amber parameters and exported GROMACS parameters using the **same double-precision evaluator**, coordinates, units and explicit mathematical conventions. Use the existing unmodified engineering gates. This separates a conversion error from differences between native arithmetic implementations.

Use several preselected configurations, including small distortions of the actual ligand, covalent attachment and/or metal site, with finite energies and full forces. An equilibrium pose alone can hide an incorrect bond force constant or a missing term. Exclude no atom because its force is inconvenient. A solvated case must retain solvent and ions in its identity inventory; a deliberately isolated parameter test must identify its narrower scope. Retain native Amber/Sander comparisons as an independent check where available, including documented legacy phase and electrostatic conventions rather than silent parameter normalization.

## 3. Native periodic comparison at one physical state

Compare native **zero-step** evaluations at the same saved coordinates and box, not independently minimized structures or later chaotic trajectory frames. If native serialization changes coordinates, recompute the reference at those exact saved values. Do not treat a geometry tolerance as permission to compare forces at different geometries. An intentional whole-molecule image transformation requires its explicit graph/map and equivalence check; arbitrary atom-by-atom wrapping is not automatically valid for bonded terms.

Record the effective settings actually used by each engine, including:

- Cutoffs, LJ combination/pair rules, Coulomb and LJ shift/switch functions, switching distance, dielectric and 1–4/exclusion treatment.
- PME/Ewald method, actual alpha, actual mesh/order, boundary and net-charge/surface conventions, precision, backend, thread/SIMD settings and any automatic tuning.
- Dispersion-correction form and averaging, excluded-pair handling, C6 versus C6/C12 tails, switching contributions and box volume.
- Bonded/exception periodicity, virtual-site evaluation, constraints and any extra forces used by the evaluation protocol.

Matching the text of an error-tolerance setting is insufficient: the engines can choose different effective PME parameters. OpenMM also describes its error tolerance as empirical, not a rigorous force-error bound; making it extremely small in lower precision can make roundoff worse. [OpenMM standard-force definitions](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html).

## 4. Numerical precision and force assessment

Report raw total and semantic component energies, all atomic force vectors, absolute RMS/max vector differences, reference force scale, and local-site metrics. A relative RMS is descriptive; near-zero reference forces need an absolute measure. A solvent-dominated global RMS alone is insufficient. Preserve full precision in native outputs and readers; a double executable followed by a float32 reader is not a double-precision check.

A native runtime profile must be frozen **before applying it to new candidates**, with exact versions/builds and an independent calibration set spanning system sizes, solvent/charge content and the supported interactions. It needs absolute maximum and RMS force limits, negative controls with missing/altered site terms and exclusions, and explicit scalar-energy criteria where energies are comparable. Repeats across threads/precision help localize numerical behavior; their observed spread is not a rigorous upper bound or a license to enlarge a threshold. Neither relative total energy nor an atom-count divisor may hide localized defects.

There is no new universal native force/energy cutoff in this contract. The existing small-system conversion gates remain unchanged. **1CA2 cannot serve simultaneously as the case that defines a relaxed threshold and the independent case proving that threshold is adequate.** An uncalibrated profile yields measurements plus `not_assessed`, not a pass. Different conventions yield `not_comparable`, not success.

## 5. Scalar energy and volume conventions

First decide whether the two raw energies describe the same Hamiltonian and conventions. Match supported settings in a separately documented diagnostic when possible, while retaining the original run settings and failed outputs. Never silently replace production settings to obtain agreement.

Energy components must be mapped by meaning: one engine can include Ewald self contributions in a different reported component. Equal component names do not establish equal definitions. Small residuals remain unassigned until supported by a reproducible component calculation or bounded error analysis; an explanation must reproduce more than a fitted number at one pose.

GROMACS documents distinct LJ combination rules and cutoff shifts; these settings must be explicit. [GROMACS nonbonded definitions](https://manual.gromacs.org/2025.4/reference-manual/functions/nonbonded-interactions.html). Its dispersion treatment also depends on excluded-pair averaging and produces energy **and pressure** corrections. [GROMACS long-range dispersion](https://manual.gromacs.org/2025.4/reference-manual/functions/long-range-vdw.html).

A volume-dependent tail contribution can leave fixed-box atomic forces unchanged while affecting pressure and an NPT volume distribution. Therefore atomic-force agreement cannot establish NPT equivalence. Where such equivalence is claimed, require consistent native virial/pressure conventions or a separately verified volume derivative for the intended Hamiltonian. Do not disable tails, self terms or other constants merely to pass. Per-configuration analytical accounting may explain differences, but raw totals and the remaining budget stay visible and are never replaced with a normalized pass value.

## 6. Report and release decision

Keep separate outcomes for **full parameter fidelity**, **same-evaluator conversion**, **native force/profile checks**, **raw absolute-energy comparability/agreement**, and **volume/pressure treatment**. Use `pass`, `fail`, `not_assessed`, `not_comparable` and `unsupported` precisely. Missing/unvalidated terms or mismatched identity are hard failures. A convention difference is not evidence of missing ligand chemistry, but also is not evidence of equivalence.

Engine-specific operation may be valid under its documented Hamiltonian even when cross-engine absolute energies are not directly comparable. That must be stated as engine-specific support with an explicit limitation, never “identical physics” or “conversion validated” on the strength of close forces alone. No overall approval is inferred from the isolated analyzer; deployment/integration must consume all required lanes and their exact evidence bindings.

Short motion, temperature/density monitoring and bond survival remain separate engineering checks. They cannot certify chemical-state choice, equilibrium, force-field accuracy, or every member of the metal/covalent panels.

## Present 1CA2 evidence and remaining boundary

The original full-system parameter and dry source-to-export checks support the declared neutral-water ZAFF6 template. The unchanged periodic input then revealed a native mixed-build/thread energy difference of 41.72 kJ/mol. The isolated double build reduced it to 2.56e−7 kJ/mol with identical coordinates and box; its compiler also differs, so precision was not the sole build variable.

The raw double-native difference from OpenMM Reference is **−0.408748287 kJ/mol**, with force RMS/max **0.00124345/0.00706216 kJ/mol/nm**. The separately recorded dispersion-component difference is **−0.457011501 kJ/mol**, leaving an algebraic **+0.048263213 kJ/mol** budget unassigned. No exact periodic-parity pass or general runtime profile follows from this single endpoint. See the original [discriminator](prototype-1ca2-periodic-diagnostic/double-precision-discriminator.md) and [independent residual review](prototype-1ca2-periodic-diagnostic/review-double-residual.md).

## Isolated executable analysis

`periodic_comparison_analysis.py` implements only the read-only numerical lane. It requires explicit units, unique atom IDs in saved order, exactly matching coordinates/box, finite arrays, declared conventions and evidence metadata. It computes unmodified energy/force differences, optional named-site metrics and source-array/metadata hashes. It has no default policy and rejects energy-offset policy fields. Its convention metadata are upstream assertions requiring source verification; it does not infer their scientific completeness.

`test_periodic_comparison_analysis.py` uses synthetic closed-form and deliberately corrupted witnesses: a single bad site hidden by global RMS, constant energy error, altered conventions, atom-order/coordinate/box/unit corruption, nonfinite forces and malformed policies. The synthetic test limits are **not** native-engine acceptance limits. No production validator or installed engine is changed by this helper.
