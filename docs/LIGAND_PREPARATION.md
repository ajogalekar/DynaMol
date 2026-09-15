# Ligand preparation in DynaMol

DynaMol prepares supported **noncovalent organic ligands** alongside the protein. It retains each bound ligand heavy atom, assigns a fixed protonation model, adds hydrogens without moving the heavy atoms, and generates real GAFF2/AM1-BCC parameters with AmberTools. The resulting OpenMM XML files travel with the prepared dataset, solvated system, and simulation result.

## Install the native runtime

Run `uv sync --frozen` for the Python dependencies, then:

```bash
./scripts/install_ligand_tools.sh
```

The script uses an available Miniforge/Conda installation and installs AmberTools 24.8 in DynaMol's private `.tools/ambertools` directory. It saves an exact platform-specific Conda package list in `.tools/ambertools-explicit.txt`. It does not modify `.gromacs` or global MD environments. An existing complete installation can be selected with `DYNAMOL_AMBERTOOLS=/absolute/path/to/environment`.

The [local standalone macOS prototype](../packaging/README.md) bundles this private runtime, including native libraries, `sqm`, `antechamber`, `parmchk2`, `tleap`, and GAFF2/BCC data. The commands above apply to the source checkout. Package signing, clean-machine validation and complete public redistribution materials remain release work.

## Chemical identity and protonation

- PDB/mmCIF chemistry is recovered from the wwPDB Chemical Component Dictionary, matching every bound heavy-atom name and element. Original mmCIF records restore full component identifiers that the PDB-based display format truncates; for example, the 58-heavy-atom **A1AHB** must not be confused with **A1A**.
- Imported SDF/SMILES and chemically supported MOL2 files use their retained authoritative chemical graph. MOL2 partial charges are source information, not a substitute for AM1-BCC charges.
- Explicit user SMILES can define a state or provide a missing chemical reference. A unique graph mapping is required; otherwise an atom-mapped SMILES must identify the bound heavy-atom order, starting at one. Heavy-atom mismatch, conflicting carbon/E–Z stereochemistry, missing atoms, and covalent ligand/polymer links are rejected explicitly.
- The default ligand pH model uses Dimorphite-DL's nominal site rules, with zero pKa uncertainty width. Changes at amide nitrogen and N–N sites are filtered and the source states retained, with a visible warning. Original heavy-atom connectivity and stereotags are preserved rather than using Dimorphite's output geometry or atom ordering.
- If the protonation tool omits an atom map, recovery requires exactly one complete element/connectivity match constrained by every surviving map. Conflicting or ambiguous maps remain blocked. This addresses the observed NAD phosphate-oxygen case and records the recovery warning.
- GNP/GppNHp uses a documented −4 reference state at pH 6–8, deprotonating its four terminal phosphate hydroxyls while retaining the imido N–H. This is a model assumption informed by the Ras nucleotide study [doi:10.1074/jbc.RA117.001110](https://doi.org/10.1074/jbc.RA117.001110), not a site-specific pKa calculation.
- The named [CCD SO4 sulfate dianion](https://www.rcsb.org/ligand/SO4) is retained at pH 6–8 only when its existing five-atom graph, bonds and charges already match that reference. This explicit fixed-state policy avoids a protonation-tool normalization that neutralized sulfate and lost two oxygen maps. It preserves all source atom identities and coordinates; it does not recover a symmetric map or predict local sulfate/bisulfate populations. Other source graphs or pH values require an explicit atom-mapped state, and a supplied state always takes precedence. This narrowly supports SO4, not arbitrary inorganic ligands.

Every ligand reports its reference and selected SMILES, formal charge, protonation method, warnings, atom map, and parameterization provenance. These are fixed-state exploratory models. They do not establish the dominant bound protonation state or rank tautomers. Carbon tetrahedral and specified E–Z configurations are checked; nitrogen inversion, atropisomerism, and phosphate stereodescriptors after ionization are not independently validated. Bound heavy-atom coordinates remain unchanged during ligand preparation.

## Parameter validation and retained evidence

AmberTools assigns GAFF2 atom types and AM1-BCC charges. Missing parameters, failed native steps, charge inconsistencies, and non-finite coordinates fail preparation. There is no Gasteiger, MMFF, uniform-charge, or zero-charge fallback for MD parameters. MMFF/UFF, when available, only optimize hydrogen placement with every bound heavy atom fixed.

When the selected state requires no added hydrogens, hydrogen relaxation is skipped because there are no movable atoms. This prevents an all-fixed RDKit minimizer failure; it does not skip native MD parameter assignment or its validation.

SQM prints individual Mulliken charges to three decimal places. A small residual can therefore remain after bond-charge corrections. DynaMol only normalizes this residual when all evidence agrees: SQM reports geometry convergence and the requested integer charge, its printed atomic-charge sum matches the native AM1-BCC sum within 0.00005 e, and the total correction is no more than 0.005 e. The tiny correction per atom and both pre/post sums are recorded; the unmodified native MOL2 is retained.

ParmEd conversion uses unique atom classes and the exact native Amber improper-torsion quartets. This avoids the force differences caused by treating all outer-atom improper permutations as interchangeable. **Every newly generated XML is compared with the native Amber topology** at the original bound pose and two deterministic perturbed poses. Energy differences must be below 0.0001 kJ/mol and the maximum force-component difference below 0.001 kJ/mol/nm.

Native Amber can reverse a quartet to avoid a sign-coded negative-zero index. The converter restores the central-third convention by reversing all four atoms, which preserves the torsion exactly; the halogenated aromatic regression checks this against native energies and forces.

Small ligands such as sulfate can have no 1–4 interactions. ParmEd's generic parameter container then emits unused global scaling defaults that conflict with Amber protein templates. DynaMol uses Amber's compatible defaults only after confirming both that there are no proper torsion parameter types and that no atom pair has shortest bond-graph distance three. No interaction in that ligand uses these changed defaults. The proof and before/after defaults are saved in `nonbonded-scale-compatibility.json` and included in conversion provenance. The same native energy/force comparison remains mandatory; existing torsion scales and any graph with a 1–4 pair are left unchanged.

Each ligand artifact folder retains input SDF, native/normalized charged MOL2, SQM output, GAFF2 supplemental parameters, Amber topology/coordinates, XML, command arguments, logs, charge-normalization evidence, and conversion-validation results. Provenance records AmberTools build, GAFF2/BCC data hashes, Python dependency versions, CCD source/hash/retrieval time where available, and the DynaMol implementation hash.

Supported organic elements are H, C, N, O, F, P, S, Cl, Br, and I, with a 200-heavy-atom limit per ligand. Each native command has a bounded ten-minute timeout and cooperates with cancellation. Metal coordination, covalent inhibitors, radicals, and other specialized chemistries require appropriate specialized models. Preserving Mg and nearby waters with a standard nonbonded ion model does not validate a coordination-shell model.

Primary implementation references: [RCSB CCD downloads](https://www.rcsb.org/downloads/ligands), [OpenMM small-molecule parameterization](https://github.com/openmm/openmmforcefields), and [Dimorphite-DL](https://github.com/durrantlab/dimorphite_dl).

See the [representative chemistry support matrix](CHEMISTRY_SUPPORT.md) for ATP, ADP, NAD, FAD, ligand functional classes, ions and expected unsupported cases. These bounded checks do not establish arbitrary cofactor or metal-model coverage.
