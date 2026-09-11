# DynaMol structure sources and preparation contract

All endpoints below are relative to `/api`. Preparation, solvation and simulation share the persistent background-job queue, cancellation controls and progress monitor. A completed operation creates a new dataset and records its parent; the uploaded source remains available.

## Structure sources

- `POST /structures/upload`, multipart `file`: supported PDB/mmCIF or MOL2/SDF/SMI/SMILES input.
- `POST /structures/smiles`, JSON `{smiles, name?, seed:2026}`: generate a small-molecule conformer.
- `POST /structures/fetch`, JSON `{provider:'pdb'|'pubchem', identifier, name?, seed:2026}`: retrieve from fixed RCSB/PubChem providers. Arbitrary URLs are not accepted.

Imports retain original files, sequence evidence and relevant chemical metadata. A SMILES conformer is a generated structure; importing a separate molecule does not dock it or merge it into a protein. Loading alone does not assign MD parameters.

## Inspection and preparation

`GET /datasets/{id}/inspection?ph=7` inspects without geometry optimization or native charge assignment. `POST` to the same path additionally accepts `ph` and `ligand_overrides`. Inspection reports missing atoms, sequence-supported missing residues, numbering gaps, ligand identity/state choices, warnings and blockers. Numbering gaps do not establish missing sequences.

`POST /preparations` returns a job and accepts:

```json
{
  "dataset_id": "dataset-id",
  "ph": 7,
  "add_missing_atoms": true,
  "build_missing_residues": false,
  "optimize_sidechains": true,
  "remove_waters": true,
  "remove_heterogens": false,
  "ligand_overrides": {},
  "seed": 2026
}
```

Optional `name` labels the result. The interface displays an active spinner, current stage, elapsed time, logs and completion or failure. Missing-loop building requires the explicit option and available sequence evidence. Rebuilt geometry remains a model, not an experimentally resolved loop.

PDBFixer repairs supported heavy atoms; OpenMM assigns protein hydrogen states at the requested pH. A bounded steric chi-angle search records its accepted changes and stereochemistry checks. These are not residue pKa calculations or exhaustive rotamer packing.

Noncovalent complexes retain supported ligands and ions by default. Observed metal-coordinating waters survive removal of other water, and coordinating protein donors are protected during chi sampling. Explicit heterogen removal prepares the protein alone. Unsupported chemistry produces a specific failure without silently deleting molecules.

## Ligand identity and parameter bundle

The automatic path supports connected, closed-shell organic ligands of at most 200 heavy atoms with elements H/C/N/O/F/P/S/Cl/Br/I and a validated chemical graph. CCD matching requires the observed heavy-atom names/elements and reference bonds. Original mmCIF ancestry preserves full component identifiers: 9AX6's `A1AHB` must not be mistaken for the unrelated CCD component `A1A` when a display PDB truncates its name.

Supplied graphs and explicit state overrides must match the observed heavy graph and defined stereochemistry. Override keys identify a residue (`chain:resid:insertion_code:residue`) or its full component ID. Atom maps refer to the bound ligand's one-based heavy-atom order. Ambiguous or incompatible mappings are rejected. Bound heavy coordinates remain fixed while hydrogens are placed around them.

Default ligand states are recorded heuristics at nominal pH. Dimorphite-DL is not a binding-site pKa calculation or tautomer ranking. Unsupported amide/N–N charge proposals are filtered and reported. GNP uses a documented fixed −4 phosphate model at pH 6–8 with imido N–H retained; outside that range the general heuristic applies and may require an override if ambiguous. Review warnings and provide an explicit SMILES override when appropriate.

Native AmberTools 24.8 runs Antechamber/SQM AM1-BCC charges, GAFF2 typing, Parmchk2 and LEaP. ParmEd converts the native Amber model to OpenMM XML. There is no guessed-charge fallback when native parameterization fails. Every new template must match native Amber energies and forces at bound and perturbed geometries. Bounded charge-rounding correction requires evidence from converged SQM output and preserves the original charges and its correction record.

The `ligands/` archive retains mapped identities, structures, native inputs/outputs/logs, charges, parameters, conversion validation and provenance. Preparation metadata has this shape (the `ligands` array contains objects, abbreviated below):

```json
{
  "ligand_parameters": {
    "forcefield": "GAFF2",
    "charge_method": "AM1-BCC",
    "requires_explicit_solvent": true,
    "ligands": [],
    "files": [{"path": "ligands/ligand-id/ligand.xml", "sha256": "64 hexadecimal characters"}]
  }
}
```

The actual bundle requires nonempty ligand and file lists. `backend.prepared_system.ligand_parameter_files` checks safe relative paths, hashes and limits, rejecting symlinks and executable/include directives. `load_prepared_forcefield` loads ff14SB, TIP3P and verified ligand templates. `copy_ligand_parameters` carries the full archive into solvent datasets, simulation inputs and trajectory outputs. Prepared hydrogen and ligand states are reused.

Install optional native tools with `scripts/install_ligand_tools.sh`; its private environment defaults to `.tools/ambertools`, overridable with `DYNAMOL_AMBERTOOLS`. Python dependencies include RDKit, ParmEd, Gemmi and Dimorphite-DL. CCD references come from the fixed official RCSB endpoint and are cached under `DYNAMOL_DATA_DIR/chemistry/ccd`. An uncached component requires network access or a sufficient supplied graph; prepared bundles do not need to re-fetch or reassign charges for simulation.

## Solvation and engine continuity

`POST /datasets/{id}/solvate`, JSON `{padding_nm:1, seed:2026, ph:7}`, returns a job building a real TIP3P cube and neutralizing ions. The preview contains actual coordinates and is reused by OpenMM. Existing coordinating waters remain present. A pH change requires re-preparation. The 100,000-atom local limit includes a conservative preallocation check; a complete prepared structure may still be too large to solvate in this alpha.

The implemented complex path is **ff14SB protein + GAFF2/AM1-BCC ligand + TIP3P explicit-water OpenMM**. Complexes require explicit solvent: bundled implicit GBn2 has no validated ligand parameters. Complex prep skips implicit relaxation; minimize the assembled explicit system before dynamics. GROMACS currently rejects prepared-state inputs because its conversion path does not preserve this parameter/state bundle.

Supported simple ions include Na, Cl, K, Mg and Ca using the bundled Amber/TIP3P nonbonded model. Preserving metal-donor coordinates does not validate coordination energetics or create metal bonds. Covalent ligands, metal-organic bonding, unsupported metals/cofactors, radicals and incomplete or ambiguous graphs require a specialized externally parameterized workflow. Nucleic-acid complexes are outside this preparation path.

Dataset `preparation` and `solvation` metadata record summaries, warnings, state choices, parent IDs, atom maps and parameter provenance. `prepared.pdb` is the exact downstream input. Successful preparation establishes parameterization support; water construction and a short MD smoke establish software continuity. Neither demonstrates equilibration, convergence, native protonation or binding affinity. See [the complex review](complex-preparation-review.md) for validation evidence.

## Viewer behavior

Simulate remains beside the live canvas. Imports and completed prep/water results load into view; solvation makes water and ions visible. Global and per-measurement visibility controls hide both lines and labels while retaining plots. Polar-hydrogen mode preserves heavy atoms and hides nonpolar hydrogens.
