# Explicit bonded-site validation

`backend/bonded_site_validation.py` checks a declared native parameter manifest against an actual OpenMM `Topology` and `System`. It does not assign chemistry, infer bonds from proximity, validate a metal oxidation state, or enable application preparation for unsupported sites.

## Parameter contract

Call `validate_bonded_site(topology, system, manifest)`. The version-1 manifest contains:

| Field | Content |
| --- | --- |
| `schema_version` | `1` |
| `model_id` | Stable model identifier |
| `scope` | `full_system` or `site` |
| `native_reference` | Independent native artifact metadata, including `sha256` |
| `atoms` | `id`, `identity`, `element`, `mass_da`, `charge_e`, `sigma_nm`, `epsilon_kj_mol` |
| `bonds` | `atoms:[2 IDs]`, `length_nm`, `k_kj_mol_nm2` |
| `angles` | `atoms:[3 IDs]`, `angle_radian`, `k_kj_mol_rad2` |
| `torsions` | Ordered `atoms:[4 IDs]`, `periodicity`, `phase_radian`, `k_kj_mol` |
| `exceptions` | `atoms:[2 IDs]`, `chargeprod_e2`, `sigma_nm`, `epsilon_kj_mol` |
| `constraints` | Independent protocol inventory: `atoms:[2 IDs]`, `length_nm`; omission means no constraints |

Each atom identity explicitly contains the strings `chain`, `resid`, `insertion`, `resname`, and `atomname`. No author/label alias, insertion code, or residue name is guessed. Missing or duplicate identities fail. `chemical_state` may carry provenance but remains explicitly **unvalidated** by this mechanics comparison.

All four term lists are required, including when empty. Values must originate from an independent native artifact rather than the System under test. For Amber harmonic parameters, OpenMM’s half-k convention requires bond k × `2×4.184×100`, angle k × `2×4.184`, Å × `0.1`, and degrees × `π/180`. Amber torsion amplitude converts by `4.184`; preserve the native divisor, ordering, phase and 1–4 scaling conventions through the native parser.

The validator matches a multiset of actual force terms, catching missing, duplicate and undeclared internal terms. Every declared chemical bond requires both a topological bond and a positive-k `HarmonicBondForce` term. Constraints and custom restraints cannot substitute. Torsion periodicities must be positive: even a zero-amplitude periodicity-zero term prevents OpenMM from creating a Context. Sigma is physically irrelevant for a particle or exception when both expected and observed epsilon are exactly zero; all other charge/LJ values are checked. Nonbonded parameter offsets are not supported and fail rather than silently changing the tested charges.

Native Amber water H–H geometry terms are an explicit exception to the chemical-edge requirement, because AmberPrmtopFile omits their topological edges. The native exporter must label the term `topology_role: native_water_hh_geometry`. The validator requires two H atoms in the same residue, both with `native_residue_name: WAT`, sharing exactly one topological oxygen within that residue. The physical harmonic term must still match. This exception cannot cover a Zn–donor or other chemical bond.

An absent **exactly zero-amplitude** periodic torsion is recorded under `zero_amplitude_torsions_omitted_by_reader`, because OpenMM’s GROMACS reader intentionally omits such null potentials. Positive periodicity is still required in the native declaration. Missing nonzero terms, even 1e-9 kJ/mol, fail this omission rule; every native Coulomb/LJ exception is still compared independently. The actual 1CA2 GROMACS-reader comparison passes with 4,265 explicitly recorded zero-term omissions in [the reviewed report](prototype-1ca2-gromacs-static/gromacs-reader-native-mechanics-null-terms-reviewed.json). Its earlier failed report is retained. This does not claim GROMACS dynamic or force-field accuracy from parameter agreement alone.

`full_system` requires complete atom coverage and rejects any additional potential-energy force without a supported parameter comparison, including `CustomBondForce` and `CMAPTorsionForce`. A correct physical bond plus an unvalidated restraint cannot pass. `CMMotionRemover` is the sole non-potential force excluded from this inventory. This standard-term validator is not a validator for arbitrary custom potentials, CMAP surfaces, or a complete MD protocol.

Every constraint within the declared atom inventory must match an independently declared atom pair and positive length. An omitted `constraints` field means zero constraints; it never authorizes deriving the expected constraints from the tested System. Missing, additional, changed, or duplicate constraints fail. A declared constraint still cannot replace a physical harmonic bond in this native-parameter check. Validate a model before any engine protocol intentionally removes bonded energy terms. Constraints remain listed as `constraints_not_used_as_bond_evidence` even when their inventory matches.

Virtual sites are not yet supported by this manifest. Any virtual site inside the validated scope fails, even if its mass and all force parameters match. Reports retain the virtual particle's identity, class, parent indices, and scope under `virtual_sites_not_validated`. A virtual site outside a partial `site` inventory remains visible as unvalidated context; that partial report cannot approve the whole System. This prevents undeclared coordinate and force redistribution from passing a term-only comparison.

`site` checks all standard terms and constraints entirely within the declared atom set and reports boundary terms without pretending to validate them. An accepted partial-scope comparison with additional forces is explicitly incomplete: inspect `scope`, `all_potential_energy_forces_covered`, and `additional_forces_not_used_as_parameter_evidence`. OpenMM does not retain a general proper/improper semantic label; the physical ordered quartet and periodic terms are what this function can check. Reports include `validator_sha256`, `system_xml_sha256`, `manifest_canonical_sha256`, and `topology_identity_sha256`. JSON hashes use sorted keys, compact separators and no NaN. The topology hash covers ordered five-field atom identities and sorted undirected bond index pairs. Artifact provenance must also be verified by the surrounding immutable-bundle workflow; trajectory reports reject a changed expected manifest.

## Trajectory contract

`audit_site_trajectory(frame_atom_ids, frames_nm, manifest, *, parameter_report, bond_bounds, motion_atom_ids, ...)` requires the exact saved frame atom-ID order and a matching accepted parameter report. Coordinates are F×N×3 in nm. Every monitored bound must refer to an explicitly parameterized bond and supply `minimum_nm` and `maximum_nm`; no bounds are inferred from instantaneous geometry.

Optional `boxes_nm` contains one 3×3 row-vector box per frame in the reduced OpenMM/GROMACS triclinic convention. Nearby lattice images are searched rather than relying only on fractional rounding. Unsupported box conventions fail explicitly. `coordination_sites` can list a `metal_atom_id` and its declared `donor_atom_ids`; the diagnostic counts only those parameterized donors within the supplied bounds and never assigns a new chemical bond from distance.

The caller selects 2–64 `motion_atom_ids`. Meaningful internal motion is measured by changes in pair distances, so a repeated frame or rigid translation alone does not pass. The default required change is 1e-4 nm and can be explicitly configured. Pair-distance monitoring is a small-site diagnostic, not a full-protein unwrapping/alignment method. Forces must be supplied in matching F×N×3 order (`forces_kj_mol_nm`); missing force evidence is reported as incomplete. `check_finite_forces(system, positions_nm)` provides a separate one-configuration OpenMM energy/force evaluation.

A passing trajectory report means the supplied frames passed the declared checks. It is not a claim of force-field accuracy, equilibrium sampling, production stability, or an unbiased simulation. Constraints and extra forces remain visible and never become accuracy evidence.

## Evidence

The [focused tests](bonded-site-tests.log) cover identity and insertion-code failures, half-k errors, missing/duplicate/reordered torsions, altered 1–4 electrostatic and LJ scaling, constraints/restraints that substitute for physical bonds, scope boundaries, periodic wrapping, skewed reduced cells, repeated/rigidly translated frames, nonfinite/missing forces, and undeclared coordination monitors.

The preserved [initial 94-atom native-cap check](native-covalent-mechanics-check.json) matched serialized terms, but is **superseded as usable-System evidence**: a later actual Context check found zero-periodicity native torsions, which OpenMM cannot instantiate. The current validator rejects those terms even when their amplitude is zero. The native covalent builder is correcting their source; the failed prototype and its earlier comparison remain recorded. This distinction prevents serialized agreement from being mistaken for an executable or accurate chemical model.
