# DynaMol structure sources and preparation contract

All API under `/api`, root owns main.py route wiring and Pydantic models. Functions own validation too when sensible. Existing Job and Dataset API unchanged except additive metadata.

## Structure sources

- POST `/structures/upload` multipart `file`: protein PDB/mmCIF or small molecule MOL2/SDF/SMI/SMILES. Sources helper `import_structure(path, name=None, provenance=None)` → Dataset.
- POST `/structures/smiles` JSON `{smiles, name?, seed:2026}` calls `create_smiles(smiles,name=None,seed=2026)` → Dataset.
- POST `/structures/fetch` JSON `{provider:'pdb'|'pubchem', identifier, name?, seed:2026}` calls `fetch_structure(provider,identifier,name=None,seed=2026)` → Dataset.
- All imports/fetches retain original files/sequence evidence and chemically relevant metadata. Sources agent owns modules; root wires run_in_threadpool endpoints. RCSB and PubChem only; no user URLs.
- Small molecules are loaded and visualized; automatic protein preparation and current protein MD presets must reject arbitrary ligands clearly. SMILES is a generated 3D conformer, not an experimentally resolved structure.

## Inspection and preparation

- GET `/datasets/{id}/inspection` calls `preparation.inspect_preparation(dataset_id)` → Inspection schema below. Must inspect quickly without optimization. Missing sequence identities must not be inferred from numbering gaps.
- POST `/preparations` body `{dataset_id, name?:string, ph:7, add_missing_atoms:true, build_missing_residues:false, optimize_sidechains:true, remove_waters:true, remove_heterogens:false, seed:2026}` → existing Job. Function `preparation.submit_preparation(settings:dict)`; persistent background worker shares jobs lock and current status/cancel API. Name defaults `Protein preparation`.
- POST `/datasets/{id}/solvate` body `{padding_nm:1,seed:2026,ph:7}` → existing Job. Function `preparation.submit_solvation(dataset_id,settings:dict)` invokes sources agent's `solvent.solvate_dataset(dataset_id,padding_nm=1,seed=2026,ph=7)` in background. Name `Explicit water preview`. Uses engine='preparation'/'solvation' for type. `config.operation` distinguishes. One worker at a time. Cancellation must match new worker module and stop native work.
- Completion creates immutable new Dataset id and `Job.dataset_id` as usual. Frontend polls jobs, loads resulting dataset automatically while keeping Simulate pane open. During preparation progress from real stages (no fake percentage); solvation turns water/ions visibility on.

Inspection fields (backend supply all):
```
{
 dataset_id:string, protein_atoms:number, hydrogen_atoms:number, water_atoms:number,
 heterogen_residues:string[], can_prepare:boolean, has_sequence: boolean,
 missing_atoms:[{chain:string,resid:string,residue:string,atoms:string[]}],
 missing_residues:[{chain:string,position:number,residues:string[],count:number,terminal:boolean,buildable:boolean}],
 gaps:[{chain:string,after:string,before:string,message:string}],
 warnings:string[], blockers:string[]
}
```

Metadata additive optional `preparation` containing `{ph,method,summary:string[],warnings:string[],...}` and `solvation` containing `{parent_dataset_id,padding_nm,water_model,added_water_atoms,...}`. Root stores these in Types as flexible fields with summary access. Preserve parent_dataset_id in metadata/provenance for ancestry; prepared geometry and protonation must actually persist into simulation. If an engine cannot honor a state, reject clearly before launching. Exact helper contracts/inspection may evolve through messages; coordinate root early.

## UI behavior

Simulate is a dock next to the live canvas, not a blurred modal. It has Upload / Fetch PDB / SMILES inputs, inspection warnings, pH controls, preparation button and explicit opt-in missing-loop building. Source/preview load doesn't close studio. Choosing explicit solvent initiates real solvation preview when input is usable; changing the box may request a refresh. No invented solvent spheres. A protein-only input can need prep first; actionable error and button.

Measurement visibility is independent of plot existence. Global scene measurements switch plus per-measurement eye controls; hide both lines and labels when disabled, retain plotted values.
