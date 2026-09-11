# Ligand conversion review

Status: the atom-mapping fix passes the scoped software checks. The raw postflight retains one unmet generic uncertainty gate, discussed below. This is an independent automated evidence review, not an endorsement by Pat Walters or expert certification.

## Question and evidence boundary

DynaMol must preserve the energy function of the native Amber topology when it writes an OpenMM force-field template. The retained 6DBK failure has 54 ligand atoms and a bound-pose discrepancy of −0.0034635454183 kJ/mol in energy and 2.8057523302 kJ/mol/nm in the largest force component. Both exceed the existing limits (0.0001 kJ/mol and 0.001 kJ/mol/nm). These differences establish a translation discrepancy; they do not establish that the selected chemical state or Amber parameters are scientifically correct.

The original job and force-field artifacts remain unchanged under `data/jobs/96b4ab9d156a4961`. Their hashes, package versions, input policy, seeds, and planned checks are recorded in `review-manifest.preflight.json`. The automated preflight reports six passes, zero failures, and zero unresolved gates in `review-preflight.json`.

## Verified results

The defect was graph-template matching: eight symmetric ligand atoms were exchanged, changing four improper torsions. The native Amber-to-XML converter itself preserved the intended terms. `review-original-mapping.json` independently records the eight assignments. The fix keeps the verified native atom identities, including hydrogens, through preparation and saved-bundle loading; it leaves force-field constants and validation tolerances unchanged.

| Check | Result |
| --- | --- |
| Original 6DBK ligand, bound pose plus two fixed perturbations | Maximum energy error 1.3713 × 10⁻⁸ kJ/mol; maximum force component error 1.4099 × 10⁻⁷ kJ/mol/nm |
| 12 additional retained ligand bundles, 36 poses | All pass; maximum energy error 1.3223 × 10⁻⁷ kJ/mol; force error 2.6395 × 10⁻⁶ kJ/mol/nm |
| Backend regressions | 94 pass, including 27 new identity tests, with zero failures or skips |
| Full 6DBK preparation | Completes in 64.86 s; 4,713 atoms, G5D ligand and PTR retained |
| Explicit TIP3P preview | Completes in 3.52 s; 50,015 atoms |
| Fresh PDB and factory reload of both systems | All 729 native ligand parameter terms match at six-decimal canonicalization; separate energy/force checks use the original production tolerances |

Evidence is in `matrix-results.json`, `full-prep-results.json`, `backend-tests.txt`, and `backend-results.xml`. The additional matrix covers ADP, ATP, FAD, NAD, acetate, a chiral amide, ethanol, a halogenated aromatic, methylammonium, phosphate, and the 9AX6 GNP/A1AHB ligands, spanning formal charges −4 through +1. The fixed tests also cover reordered atoms/bonds, two ligand copies, hydrogen addition, saved provenance and PDB reload, and deliberate corruption of numeric parameters or identity metadata. Original input artifacts remain unchanged.

The full worker began before the final constructor-consistency guard was added. Final regressions, the matrix, and fresh prepared/solvated-system reloads use the current code; source hashes attached to the reload report do not prove the bytes loaded by the earlier worker. The last guard did not alter numeric mapping behavior.

Diagnosis narrowed the initial suspected converter problem to atom identity. No torsion conversion algorithm was rewritten; identity regressions and exact native-term comparison replaced speculative synthetic converter expansion. Existing limitations on unequal duplicate improper terms remain.

## Review criteria

- Compare the retained native prmtop against the new XML on OpenMM's Reference platform, with the same coordinates, no constraints, and NoCutoff. Keep existing absolute energy and force limits. Preserve all force terms, charges, masses, atom mapping, exclusions, and 1–4 interactions; no chemical substitutions or tolerance relaxation.
- Test the original bound pose plus deterministic Cartesian perturbations. Planarity can conceal a wrong improper atom order, so forces and nonplanar poses are required alongside total energy.
- Exercise reversed Amber quartet encoding, multiple torsion terms, aromatic and amide centers, equivalent outer atoms, charged molecules, distinct ligand copies, and template matching under atom reordering where supported. Include a deliberately corrupted parameter file that must fail validation.
- Confirm the full 6DBK preparation path with the bound ligand and modified residue retained. Report any failure separately from the isolated conversion check. Parameter equivalence does not demonstrate equilibration, accurate free energies, a preferred protonation state, or universal ligand/cofactor support.

These checks are software equivalence tests, not a predictive modeling benchmark. There is no training/test split, parameter fitting, significance test, or experimental uncertainty estimate. Report all fixed poses' errors and the maximum, rather than implying that a few numerical comparisons estimate chemical or statistical uncertainty.

## Supporting evidence

OpenMM's force-field guide states that residue templates match elements and connectivity, independently of atom names and order. It also distinguishes proper and improper ordering and allows multiple periodic terms. This means template matching and quartet ordering must both be inspected; a visually unchanged molecule alone is insufficient. See [OpenMM 8.6 force-field creation](https://docs.openmm.org/latest/userguide/application/06_creating_ffs.html).

ParmEd documents Amber impropers with the third atom as the center. Its exporter moves the center to the first XML position; the eventual physical quartet is determined by OpenMM's improper matcher. The installed versions are ParmEd 4.3.1 and OpenMM 8.6.0; the implementation under test, rather than current online source alone, determines actual behavior. See [ParmEd Dihedral](https://parmed.github.io/ParmEd/html/topobj/parmed.topologyobjects.Dihedral.html) and [ParmEd OpenMM exporter source](https://parmed.github.io/ParmEd/html/_modules/parmed/openmm/parameters.html). The [OpenMM PeriodicTorsionForce API](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.PeriodicTorsionForce.html) provides direct access to particle indices, phase, periodicity, and force constant for term-level comparison.

The complete local public Practical Cheminformatics corpus was refreshed September 10, 2026 and includes 91 posts: 74 archived Blogger posts and 17 GitHub Pages posts. Full-corpus searches for representation/stereochemistry/validation, reproducibility/code sharing, and force-field/improper/torsion are saved in `review-corpus-search.json`. No blog result is treated as direct authority for an Amber/OpenMM implementation detail.

The review applies three general recommendations as an inference to this conversion problem: validate representations and benchmark inputs ([Pat Walters, “We Need Better Benchmarks for Machine Learning in Drug Discovery,” 2023](https://practicalcheminformatics.blogspot.com/2023/08/we-need-better-benchmarks-for-machine.html)); avoid silently changing a molecular state while standardizing a representation ([Pat Walters, “The Trouble With Tautomers,” 2025](https://patwalters.github.io/The-Trouble-With-Tautomers/)); and retain executable code and procedures for reproducibility ([Pat Walters, “Where's the code?”, 2019](https://practicalcheminformatics.blogspot.com/2019/05/wheres-code.html)). These are all Pat-authored posts; no guest writing is attributed to him.

## Postflight disposition

`review-postflight.json` reports seven passes and one failed gate, `uncertainty`, with raw decision `block`. That result is retained without alteration: no statistical or biophysical uncertainty intervals were calculated. A few deterministic error values are not confidence intervals. This gate blocks claims of predictive or physical accuracy; it is not evidence of a failed atom-mapping equivalence test. The scoped software correction is supported by the numerical and integration results above. There is no new claim of force-field accuracy, experimental protonation correctness, simulation convergence, or universal ligand/cofactor coverage.

`review-findings.json` records no unresolved blocking implementation defect. The small adapter uses two private OpenMM hooks, so upgrading OpenMM requires these regressions. Exact identity checks reject missing or inconsistent mappings rather than falling back to symmetric graph matching for generated DynaMol ligand templates.
