# Protein–ligand preparation: evidence and review

This is an implementation and evidence review by a language model, supported by the local PatAgent checklist. PatAgent is not Pat Walters, does not speak for him, and is not expert certification. A successful parameterization or short MD run establishes software operability, not an experimentally correct molecular state, equilibrium sampling, affinity, or accuracy of metal coordination.

## Preflight scope and decision

The preflight manifest and audit are `docs/audit/complex-manifest.preflight.json` and `docs/audit/complex-preflight.json`. The automated checklist reports **proceed: six passes, no blocking gates**, for implementation and bounded parameterization/CPU smoke checks. This was the design-stage decision; completed implementation and postflight findings are recorded below.

The intended behavior is to retain protein, ligands, cofactors, ions, and coordinated waters by default. Ligand preparation must restore a chemically defined graph, add hydrogens in a documented molecular state, produce actual force-field parameters, and keep those parameters and the bound pose through solvation and simulation. A label saying “prepared” is insufficient evidence.

## The reported 9AX6 case

Both local datasets, `9fb26d98962c453b` and `0bc7b7de87ec4e1d`, have identical source mmCIF and physical coordinates. The source hash is recorded in the preflight manifest. The imported topology has 5,936 atoms, including two 32-heavy-atom GNP residues, two Mg ions, and two 58-heavy-atom inhibitors.

**The topology PDB truncates the inhibitor identifier A1AHB to A1A. These are different CCD components.** A1A has 43 heavy atoms and unrelated chemistry. Resolve component identity from the preserved mmCIF, with an explicit atom mapping; never look up the truncated residue name blindly. The correct inhibitor is A1AHB/RMC-6236, with six CCD stereo flags, including a pyramidal nitrogen flag. The CCD descriptor providers differ in whether they retain that nitrogen stereo tag. [RCSB 9AX6](https://www.rcsb.org/structure/9AX6), [A1AHB CCD](https://www.rcsb.org/ligand/A1AHB), [A1A CCD](https://www.rcsb.org/ligand/A1A)

The source contains 12 `struct_conn` records, all marked `metalc`; it reports no covalent protein–ligand link. Each Mg coordinates a serine oxygen, a threonine oxygen, two GNP phosphate oxygens, and two waters. The recorded distances span 1.961–2.181 Å. Those waters are part of the initial coordination environment and should not be silently stripped by a generic “remove waters” default. `docs/audit/complex-source-inventory.json` preserves these records.

An independent, non-optimizing RDKit graph/coordinate check matches all five A1AHB carbon stereo flags and all four GNP carbon flags in both copies. The neutral GNP graph also matches its two phosphorus flags; phosphate deprotonation can remove distinctions between equivalent oxygens. RDKit does not retain the A1AHB N3 tetrahedral tag, so this check is not validation of every CCD stereochemical feature or of macrocycle atropisomerism. Results are in `docs/audit/complex-bound-stereochemistry.json`.

The CCD formal charges are A1AHB 0, GNP 0, and Mg +2. **CCD charge is a component reference, not a pH prediction.** GNP's neutral dictionary form contains protonated phosphate oxygens. Spectroscopy and QM/MM work on Ras supports deprotonated GppNHp phosphates, while discussing earlier contradictory neutron-crystal evidence. A GNP −4 state with the imido NH retained can be a documented starting assumption for this Ras-like complex; it must not be presented as a universal ligand pKa calculation. [GNP CCD](https://www.rcsb.org/ligand/GNP), [Mann et al., 2018](https://doi.org/10.1074/jbc.RA117.001110)

## Chemical-state and parameter requirements

- Resolve the ligand from authoritative atom-named CCD data or explicit user chemistry, validate elements, connectivity, bond orders, formal charges and defined stereochemistry, and retain all bound heavy-atom coordinates. Reject mismatches instead of guessing a different component or re-embedding a macrocycle.
- Expose and record the ligand molecular state. A protonation enumerator returns possible states; output order is not a population ranking. Avoid automatic canonical-tautomer substitution. Protein pH preparation does not set a ligand's state by itself.
- Use a compatible protein/water/ligand force-field combination and actual ligand partial charges. Record package versions, exact force-field identifier, charge method, atom-name mapping, parameter files, seeds, input hashes and logs. Sum the partial charges and verify agreement with the selected formal charge within numerical tolerance.
- Carry the exact ligand templates and selected molecular state into solvation and simulation. The prepared component graph, atom count, atom identities and charges must survive every serialization step. A protein-only `ForceField` reconstruction must not quietly replace the prepared system.
- Treat covalent inhibitors, unsupported elements, radicals, incomplete ligand heavy atoms, ambiguous identity and unsupported chemistry as explicit, component-specific failures. Keep the original complex intact when a step fails. A browser viewer can still display a structure that the selected MD route cannot parameterize.

`openmmforcefields` documents GAFF small-molecule templates and AM1-BCC charge generation through AmberTools. It requires explicit protonation/tautomer state and stereochemistry because an OpenMM topology lacks full chemical identity. A charge-assignment conformer must not overwrite the bound simulation pose. OpenFF documents that unmapped SMILES do not guarantee atom indexing; mapped inputs or verified graph mappings are needed. [OpenMM force fields](https://github.com/openmm/openmmforcefields), [OpenFF molecule construction](https://docs.openforcefield.org/projects/toolkit/en/latest/users/molecule_cookbook.html)

Dimorphite-DL is an open-source rule-based protonation enumerator with a pH range and bounded variants. An upstream report currently identifies erroneous amide/anilide deprotonation with version 2.0.2; that report is a reason to add a local regression check, not independent proof that every version fails. Do not take the first returned state as dominant, or hide enumerated ambiguity. [Original methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6689865/), [Upstream issue 12](https://github.com/durrantlab/dimorphite_dl/issues/12)

## Metal scope

A standard Mg nonbonded model is an explicit approximation. Merely loading a Mg parameter or observing finite energy does not validate bound phosphate coordination. Published comparisons show sensitivity to parameter choice, water model, and coordination energetics. Retain the crystal coordination environment, report the selected ion/water model, and make any restraints explicit. Do not add fabricated covalent bonds to force a metal site through an organic-ligand parameterizer. [Optimized Mg parameters](https://pmc.ncbi.nlm.nih.gov/articles/PMC8047801/), [ATP–Mg force-field comparison](https://pubmed.ncbi.nlm.nih.gov/33616388/)

## Public-writing retrieval

The existing complete public Practical Cheminformatics corpus contains 91 records: 74 archived Blogger and 17 GitHub Pages posts, meeting the skill's minimum baseline. This task used targeted retrieval rather than claiming a fresh full-corpus sync; results are retained in `docs/audit/complex-blog-search.json`.

Pat Walters's [Assigning Bond Orders to PDB Ligands](https://practicalcheminformatics.blogspot.com/2018/09/assigning-bond-orders-to-pdb-ligands.html) motivates using known component chemistry rather than inferring bond orders solely from imperfect PDB geometry. [The Trouble With Tautomers](https://patwalters.github.io/The-Trouble-With-Tautomers/) explains that algorithmically generated tautomers need not be physically populated. Our implementation checks extend those general recommendations to this app; the resulting gates are reviewer inferences, not statements attributed to the author.

## Required postflight evidence

1. A real protein–ligand fixture completes preparation without losing the ligand, and preserved heavy-atom coordinate changes are quantified.
2. Atom-named CCD mapping distinguishes A1AHB from A1A; charges, hydrogen counts and stereocenters are validated for the actual selected states.
3. Parameterization and a finite-energy/force check succeed with recorded native outputs. Any bounded MD smoke is labeled a software check.
4. Prepared-system parameters and molecular state survive solvent preview, job submission, worker startup and trajectory output; no duplicate solvent or silent reparameterization.
5. Failure fixtures exercise wrong CCD identity, unsupported/covalent chemistry, missing parameters and native-tool errors without deleting source components.
6. User controls show progress, completion, components retained, molecular state and actionable errors.

## Implementation review and corrected findings

The reviewer implemented `backend/prepared_system.py` and continuity changes in the solvent, job and simulation modules. Review of those contributions is self-review, supplemented by separate artifact and integration checks by other agents. This is not independent expert certification.

The loader verifies relative XML paths and SHA-256 hashes, excludes symlinks, executable XML and unverified includes, and feeds verified bytes to OpenMM. ff14SB + TIP3P + GAFF2/AM1-BCC parameters and complete native ligand archives persist through solvent previews, simulation inputs and trajectory outputs. Implicit ligand models and prepared-state GROMACS conversion remain unavailable.

The initial A1AHB Amber-to-XML conversion exposed an improper-torsion ordering mismatch. The corrected conversion preserves native Amber improper quartets with atom-specific classes and Amber ordering. A mandatory energy/force comparison at the bound pose and two deterministic perturbations now guards every new ligand template. The independently recomputed maximum difference across A1AHB and GNP was 8.46×10⁻⁸ kJ/mol in energy and 1.26×10⁻⁶ kJ/mol/nm in a force component, within the recorded 10⁻⁴ and 10⁻³ tolerances. These quantify conversion agreement, not force-field physical accuracy. [Native artifact checks](audit/complex-preparation-checks.json)

The requested source-chemistry guards are implemented: original PDB LINK/mmCIF covalent-link checks, rejection of ambiguous unmapped overrides, defined carbon and E/Z stereochemistry comparisons, fixed-host CCD fetches without redirects and with size bounds, and native-tool/data/version/command provenance. The ligand suite includes regression cases for those mapping/stereo/link paths, unexplained charge residuals, native template equivalence, hydrogen preservation and identical-copy provenance. The full backend run reports **77 passed** in [the test record](audit/complex-backend-tests.txt). Intentional failure fixtures do not imply a chemistry accuracy benchmark.

The earlier focused archive tests used a synthetic methane model with real OpenMM template checks and water construction, mocking only the final dynamics step to isolate persistence. Native results below provide the separate real-chemistry and MD evidence.

## Completed molecular and native software checks

The full uploaded 9AX6 dataset was prepared as dataset `4cf5bcedf38847c1` by job `57083a2415754858`: 10,654 atoms, all four ligands, both magnesium ions and four coordinating waters. Independent comparison after the PDB round trip found **exactly zero displacement** of all 180 ligand heavy atoms, both Mg atoms, all four water oxygens and all 26 heavy atoms in the protected Ser17/Thr35 residues. All twelve metal-donor distances were unchanged. Per-atom ligand charges in the assembled OpenMM system exactly matched the native Amber values, with GNP −4 and A1AHB +1 to numerical precision. [Full-structure checks](audit/complex-preparation-checks.json)

A source-version limitation remains explicit: the full prep job ran while development edits were occurring. Its recorded completion-time implementation hash is not a frozen snapshot of loaded Python modules. Final source adds guards and a repeated-copy provenance fix after that job started. The original shared-folder annotation incorrectly identified the last identical ligand copy as the parameterization source; the annotation was corrected to identify the first copy, with original annotations and old/new hashes preserved. All 104 checked native/parameter files remained byte-identical. The independent audit rechecks the actual molecular artifacts, and current-code regressions cover the corrected guards and provenance behavior. This does not retroactively turn the earlier run into an immutable build.

The complete structure remains intact. Its full explicit-water cube would exceed the unchanged 100,000-atom alpha cap. For a bounded native test, the documented fixture selects **whole KRAS chain A, whole PPIA chain C, GNP E201, Mg F202, inhibitor I201 (full CCD A1AHB), and waters K307/K348** from the prepared parent. In source identifiers, PPIA/inhibitor author chain D maps to canonical C/I. Original CIF sequence/CCD evidence and parent-to-fixture atom maps are preserved. This selection removes contacts with the second crystallographic copy and is not a separately validated biological-assembly determination. [Fixture builder and rationale](audit/complex-native-fixture.py), [selection manifest](audit/complex-native-fixture.json)

The fixture has 5,327 prepared atoms and an estimated 79,526-atom upper bound at 1 nm cube padding. The live API produced a **64,494-atom TIP3P system** (`865380bba4a54a53`) with zero displacement of the prepared solute. Job `f66ad52d18834ec1` performed minimization, ten initial dynamics steps and ten production steps (0.02 ps production), saving six frames to `7549b750c03f4817`. XML hashes and every ligand atom's charge survived through the saved native system and output dataset; energy records and coordinates were finite. These are a real CPU software continuity check, not equilibration or a useful sampling trajectory. [Native run report](audit/complex-native-run.json)

The UI audit separately covers state overrides, progress, actual canvas visibility, preparation/water completion and downloads. Its own report carries the current test totals and device/browser limits. [UI functional audit](UI_FUNCTIONAL_AUDIT.md)

## Postflight decision and unresolved scientific scope

The postflight manifest records current source/build hashes, exact artifact hashes, native versions, seeds, exclusions and results. The automated PatAgent checklist reports **seven passes and one failed gate: uncertainty**, with decision **block**. No replicate distributions or confidence intervals were generated, and the manifest explicitly sets `uncertainty_intervals_reported` to false. [Manifest](audit/complex-manifest.postflight.json), [checklist output](audit/complex-postflight.json)

The deterministic software contracts and bounded native run passed. The failed general modeling gate blocks promotion to scientific accuracy, uncertainty or convergence claims; it has not been relabeled a pass or bypassed with invented intervals. The deliverable is tested exploratory application support for the documented noncovalent chemistry and explicit-water OpenMM path. Scientific use still requires appropriate state selection, metal modeling, equilibration, replicates and convergence analysis for the actual system. Nitrogen inversion, macrocycle atropisomerism, general bound-state pKa accuracy, untested chemistry and packaging/installer portability remain unvalidated.
