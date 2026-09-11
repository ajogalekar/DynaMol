# DynaMol local API

Start from the project root with `.venv/bin/python -m uvicorn backend.main:app --host 127.0.0.1 --port 8765`. The interactive API reference is at `http://127.0.0.1:8765/docs`. `DYNAMOL_DATA_DIR` changes the default `data/` directory; `DYNAMOL_CPU_THREADS` selects 1–4 threads (default 2). Use one API process. The runner starts one independent simulation process at a time, persists progress in JSON, and retains intermediate files after failure or cancellation. Engine availability probes are shared and cached for 30 seconds across browser polling.

## Molecular files and measurements

Topologies: PDB, mmCIF/PDBx, GRO, MDTraj HDF5. Trajectories: DCD, XTC, TRR, NetCDF, PDB, MDTraj HDF5, MDCRD, LAMMPS text and XYZ. Matching topology and trajectory atom order is mandatory. Atom signatures are compared for self-describing trajectories; binary trajectories can only be checked for matching atom count. Uploaded original files, source names, SHA256 hashes, stride and time overrides are retained with the canonical dataset.

Each file is limited to 250 MiB. Retained coordinate arrays are limited to 100,000 atoms, 10,000 frames and 256 MiB. Use stride or an externally trimmed trajectory for larger data. Small-trajectory coordinates are served as a complete float32 little-endian `[frame,atom,xyz]` array in Ångströms. MDTraj HDF5 requires PyTables, included in the Python dependency set.

Original physical coordinates in nanometers and times are preserved in `physical.npz`. The viewer receives separate protein-aligned coordinates, with periodic molecules made whole when topology permits. Measurements always use the original saved frames; display interpolation and alignment do not alter analysis. Distances, angles and dihedrals use MDTraj minimum-image geometry when periodic box vectors are present. Cartesian geometry is explicitly identified when they are absent. Torsions use the signed −180° to 180° convention. Coincident atoms and collinear torsion planes are rejected.

DCD reader timestamps, structure-only files and trajectory formats without trusted physical times are represented as original frame indices (`time_unit="frame"`). A user-supplied `frame_interval_ps` is per original frame, before stride. Internally generated OpenMM trajectories receive exact step-derived times; XTC and a CSV sidecar retain them. A DCD reader may discard timing metadata. The final OpenMM frame may use a shorter report interval to include the final simulated step.

A hydrogen-bond selection is **donor, explicit hydrogen, acceptor**. Hydrogen connectivity is validated. The custom geometric occupancy is the fraction of saved frames with D–A ≤ 3.5 Å and D–H–A ≥ 150° (angle at H). N/O/S element screening does not establish donor/acceptor chemical eligibility, protonation or electronic state. This criterion is not a claim of equivalence with another analysis tool's defaults.

## Real simulation workflows

OpenMM: Amber ff14SB protein parameters; GBn2 implicit solvent or TIP3P explicit water with PME and neutralizing ions. Langevin-middle dynamics, hydrogen-bond constraints, up to 2 fs timestep, fixed-volume NVT. Hydrogens are added with pH 7 template heuristics. Optional minimization uses tolerance 10 kJ/mol/nm and at most 1,000 iterations. Python/NumPy preparation RNG, integrator and velocities are seeded in current jobs.

GROMACS: Amber99SB-ILDN protein parameters, TIP3P water, PME and stochastic dynamics (`sd`), hydrogen-bond constraints, fixed-volume NVT. `pdb2gmx -ignh` explicitly rebuilds hydrogen coordinates; input heavy-atom retention is checked. Protein protonation/termini use GROMACS defaults. Solvent coordinates use the standard `spc216.gro` packing template, with TIP3P parameters selected in topology. An explicit index restricts neutralizing-ion replacement to newly added solvent; existing input waters are retained. Optional steepest-descent minimization uses at most 2,000 steps and an emtol of 1,000 kJ/mol/nm. No grompp warnings are suppressed. Set `DYNAMOL_GMX` or put `gmx` on PATH; project-local `.gromacs/bin/gmx` and `.tools/gromacs/bin/gmx` are also discovered.

RNA and DNA residues are identified as nucleic-acid polymers for viewing and analysis. Pure nucleic-acid systems and protein–nucleic-acid complexes are explicitly rejected by the standard-protein simulation preset.

These are deliberately constrained protein workflows. Unknown ligands, cofactors and nonstandard residue templates are rejected; no ligand parameterization is implied. Explicit solvent is required for GROMACS. Implicit OpenMM rejects existing water/ions instead of silently removing them. No pressure equilibration is performed, and a short relaxation stage is not evidence of equilibration. The engines use different forcefields and integrators and are not presented as identical simulation protocols.

Each completed job includes its config, input, preparation record, engine versions, seed, source hash, stable output hashes, topology, trajectory and engine logs. OpenMM also saves XML system/integrator files, an energy CSV and a checkpoint; GROMACS retains its MDP/topology/TPR/checkpoint files. Download the output ZIP after completion or cancellation. Native engine calls are terminated as a process group on cancellation; partial output may be incomplete.

The first bundled 1UBQ demonstration predates preparation RNG seeding. Its dynamics seed is recorded, but its preparation randomness was uncontrolled, as stated in its provenance. It has 2 ps of production after 0.2 ps of initial relaxation, and is a real software demonstration, not a converged scientific study. Its starting structure explicitly excludes 58 crystallographic waters for implicit solvent.

## Verification

Run `.venv/bin/python -m pytest -q` for API, unit/PBC geometry, explicit H-bond selection, identity-order checking, file timing, retained uploads and error handling. `docs/audit/check_numerics.py` supplies a separate analytical fixture audit. Actual OpenMM and GROMACS integration jobs and cancellation were also executed during this build; their evidence is recorded in the scientific audit and local `data/` artifacts.
