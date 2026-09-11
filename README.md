# DynaMol

**Molecules in motion.** A local, open source molecular dynamics workbench: create a trajectory, explore it in 3D, and follow the interactions that matter.

![DynaMol molecular workspace](docs/dynamol-desktop.png)

DynaMol is a working **v0.1 alpha**, built with React, NGL, FastAPI, MDTraj, OpenMM, GROMACS, PDBFixer, and RDKit. It runs on your own computer. Uploads and simulations stay local; optional Fetch contacts RCSB or PubChem for the requested identifier. The application needs no cloud account or API key. Fonts and the starter trajectory are bundled for offline use after dependencies are installed.

## Start exploring

Install [Node.js 22+](https://nodejs.org/) and [uv](https://docs.astral.sh/uv/getting-started/installation/), then run:

```bash
./start.sh
```

Open **http://127.0.0.1:8765**. The launcher installs locked Python dependencies, builds the frontend, and restores the bundled example into a fresh workspace. On macOS, double-click **DynaMol.command** to launch and open the browser. Keep the terminal running while using the app. Simulation workers are separate processes and can continue when you close the browser or stop the web server; reopen the app to see their progress or cancel them.

The initial scene contains a **real 2 ps OpenMM trajectory of ubiquitin**: 101 saved frames, 1,231 atoms, Amber ff14SB/GBn2, 300 K. Press Play or Space. A real Cα distance trace is ready to explore below the timeline. Its 3D annotation starts hidden: click the eye on its plot tab or the viewer toolbar to show it. The eye preserves the plot; × removes the measurement. This short run demonstrates the software; it does not establish equilibration, convergence, or a biological conclusion.

### Enable GROMACS

OpenMM installs with the Python dependencies. GROMACS is a separate optional native dependency. Install a compatible `gmx` executable and put it on PATH, or use a project-local Conda environment:

```bash
conda create -p .gromacs -c conda-forge gromacs=2025.4 -y
```

DynaMol also accepts `DYNAMOL_GMX=/absolute/path/to/gmx`. The current workspace already has **OpenMM 8.6 and GROMACS 2025.4** installed and tested. Availability and engine versions appear in Simulation Studio.

## Inside the workbench

**Explore.** Ribbons by default, with ball-and-stick, sticks, and molecular surface views. Rotate by dragging, zoom with the wheel or explicit buttons, pan with right-drag, and focus on a residue by name or number. Show or hide proteins/nucleic acids, ligands, waters, ions, all hydrogens, or only polar hydrogens. Save the scene as a PNG. Short or interrupted protein fragments are shown as atoms when a ribbon cannot be drawn.

In ribbon view, **Polar only** retains a light trace of the full heavy-atom structure with attached polar hydrogens; non-polar hydrogens are hidden. **All** adds the remaining hydrogens, and **Hidden** returns to the ribbon view.

**Play.** Frame stepping, scrubbing, looped playback, variable speed, and a linked timeline. Nonperiodic trajectories use interpolated visual transitions without changing saved coordinates. Periodic trajectories currently play saved frames, avoiding interpolation across box discontinuities. The camera is preserved throughout playback.

**Measure.** Pick two atoms for distance, three for an angle, or four for a signed dihedral. Pick donor–hydrogen–acceptor for a geometric hydrogen-bond trace and occupancy. Atom search supplements direct 3D picking. Plots use actual timestamps, support click-to-seek, and export CSV; torsion summaries use circular statistics. Measurements come from original saved coordinates with minimum-image periodic geometry when valid box data exists. Interpolation and display alignment do not affect the calculated values.

**Prepare.** Open Simulate to upload a PDB, mmCIF, MOL2, SDF or SMILES file, fetch a PDB accession or PubChem molecule, or paste SMILES to generate a local 3D conformer. Each result loads immediately beside the controls. Click **Prep protein** to rebuild hydrogens for the selected pH, repair missing heavy atoms and refine side-chain rotamers/clashes with a restrained backbone. Inspection warns about missing sequence and backbone gaps; known short internal loops can be built with an explicit checkbox. The resulting structure, settings and warnings are saved as a separate dataset.

The Prep button spins immediately. A monitor above the scrolling controls shows the current stage, completed-stage percentage, elapsed time, cancellation and logs. It reconnects when the studio is reopened and distinguishes a completed preparation from a structure still loading in the viewer. Errors stay visible, with an Open result action to retry a failed display load.

**Solvate.** After preparation, select **Explicit water**. A background job builds an actual TIP3P box, displays water and ions, and saves the system. OpenMM reuses those coordinates and protonation states. Update box applies the chosen padding; switching back to implicit returns to the unsolvated parent. A solvent preview has not been equilibrated.

**Simulate.** Choose OpenMM or GROMACS, configure temperature, duration, solvent and optional advanced settings, then start a real local job. Follow stages, step counts, progress and logs while exploring another trajectory. Cancel a running job, open the completed result, or download its files with configuration, seed, versions, provenance, and engine output.

### File support

| Input | Supported formats |
| --- | --- |
| Simulation structure | PDB, mmCIF, MOL2, SDF/MOL, SMILES; RCSB/PubChem fetch |
| Trajectory topology | PDB, mmCIF, GRO, MDTraj HDF5 |
| Trajectory | XTC, DCD, TRR, NetCDF, multi-model PDB, MDTraj HDF5, MDCRD, LAMMPS text, XYZ |
| Exports | PNG scene, CSV measurements, simulation output ZIP |

Topology and trajectory must have identical atom ordering. Self-describing trajectories undergo atom-identity checks; binary formats only allow atom-count validation. Original uploads and hashes are retained. If a reader cannot supply trustworthy timing, DynaMol labels its axis as frame indices; enter an interval in picoseconds during import to provide physical time. The interval is per original frame, before stride.

Imports are limited to **250 MiB per file**, **100,000 atoms**, **10,000 retained frames**, and **256 MiB of retained coordinates**. The browser receives the retained coordinate buffer in full. Use import stride or trim large trajectories externally. This alpha does not yet stream multi-gigabyte trajectories.

The Simulate structure loader accepts one molecule per SDF/MOL2/SMILES file, with a 25 MiB source limit and a 500-atom limit for RDKit small molecules. MOL2 requires explicit supported Tripos atom/bond types and contiguous residue atom blocks. Chemical connectivity, charges and original files are retained; importing or generating a ligand does not assign MD parameters.

### Keyboard shortcuts

| Key | Action |
| --- | --- |
| Space | Play / pause |
| ← / → | Previous / next frame |
| F | Fit molecule |
| M | Toggle atom picking |
| Esc | Clear picking / close dialog |
| ? | Quick guide |

## Simulation scope

Automatic preparation currently supports **standard amino-acid proteins**. Unsupported ligands, cofactors, nonstandard residues, and nucleic acids are rejected with an explanation. DynaMol does not silently strip user-uploaded chemistry or parameterize arbitrary ligands.

- **OpenMM:** ff14SB with GBn2 implicit solvent, or TIP3P explicit solvent with PME. Langevin-middle, hydrogen-bond constraints, up to 2 fs integration steps. Implicit mode requires a protein-only input.
- **GROMACS:** Amber99SB-ILDN/TIP3P, PME, stochastic dynamics, hydrogen-bond constraints, explicit solvent. Hydrogen reconstruction is recorded; existing heavy atoms are retained. Neutralizing ions replace only newly added solvent.
- Both presets use **fixed-volume NVT**, optional minimization, and a short initial relaxation. They do not include pressure equilibration, site-specific pKₐ prediction, production convergence assessment, ligand parameterization, or GPU selection. Different engine defaults are not equivalent physical protocols.
- Prepared structures and solvent previews currently run with **OpenMM**. GROMACS transfer is explicitly blocked because its preset rebuilds hydrogen states and uses a different force field; unprepared standard-protein inputs retain the original GROMACS workflow.
- Missing-loop building is limited to sequence-supported internal gaps of at most 6 residues each and 12 residues total. Terminal extensions and larger gaps require external modeling. PDBFixer loop coordinates and bounded side-chain sampling are starting models, not validated native conformations. Ionization uses pH/template heuristics; existing hydrogens are removed before rebuilding at a new pH. Rebuilt/relaxed structures undergo template-based stereochemistry checks; invalid models are rejected.
- One simulation, preparation or solvation job runs at a time, using two CPU threads by default. Set `DYNAMOL_CPU_THREADS` to 1–4. Progress reflects completed production steps; preparation stages do not have a fabricated percentage.

The hydrogen-bond tool reports a **custom geometric criterion**: donor–acceptor ≤ 3.5 Å and donor–H–acceptor ≥ 150°, with explicit hydrogen connectivity required. This is not proof of chemical donor/acceptor eligibility or equivalence to other packages' default definitions.

See [preparation review](docs/preparation-review.md) and [scientific review](docs/scientific-review.md) for numerical validation and limitations, and [backend notes](backend/README.md) for exact preparation behavior. The first bundled demo has an explicitly recorded preparation-randomness limitation; new jobs seed preparation as well as dynamics. No convergence or uncertainty claims are made for the demonstration runs.

## Develop and contribute

```bash
./start.sh --dev                         # API 8765 + Vite 5173
.venv/bin/python -m pytest -q            # Backend tests
npm --prefix frontend run build         # TypeScript + production build
cd frontend && npm run test:e2e         # Live Chrome + real local API/engines
```

Browser tests require Google Chrome installed and the local app running. Set `DYNAMOL_BASE_URL=http://127.0.0.1:8765` to test the production build. They use actual molecular data and a tiny real OpenMM job, and create test datasets/jobs in the local workspace.

```text
frontend/src/        React interface, NGL viewer, plots, import & simulation panes
backend/             FastAPI API, persistent jobs, native engine workers, analysis
examples/            Bundled real demo and its provenance
scripts/             Bootstrap and optional demo regeneration
docs/                API contract, numerical audits, review and screenshots
data/                Local uploads, trajectories and jobs (gitignored)
```

The production UI is served by FastAPI at the same localhost origin. The API reference is available at `/docs`. Use one backend process. `DYNAMOL_DATA_DIR` changes the data location; environment settings are illustrated in `.env.example`.

To regenerate the example with the current preparation pipeline:

```bash
.venv/bin/python scripts/generate_demo.py --force
```

This downloads public PDB 1UBQ if necessary and starts a short real CPU simulation. It explicitly selects protein for the implicit-solvent demo and records the excluded crystallographic waters. Regeneration replaces the local `demo` dataset; ordinary uploads are unaffected.

Good next contributions: bounded trajectory streaming, verified continuous periodic display paths, validated prepared-state transfer to GROMACS, ligand parameterization, GPU controls, reusable projects and selections, RMSD/RMSF/contact maps, and accessible atom tables for very large structures. Please include representative molecular files or reproducible fixtures with format/analysis changes, and preserve the separation between display coordinates and physical analysis.

## Credits and license

DynaMol application code is [MIT licensed](LICENSE). Dependencies retain their own licenses; see [third-party notices](THIRD_PARTY.md). The starter structure is [RCSB PDB 1UBQ](https://www.rcsb.org/structure/1UBQ), from Vijay-Kumar, Bugg & Cook, *Journal of Molecular Biology* (1987), [doi:10.1016/0022-2836(87)90679-6](https://doi.org/10.1016/0022-2836(87)90679-6).
