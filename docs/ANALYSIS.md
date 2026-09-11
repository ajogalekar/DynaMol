# Structural analysis, readiness and live diagnostics

Expand **Structural analysis** below the geometry plot. RMSD measures change from a chosen saved reference; RMSF measures fluctuations about each atom's mean position over a frame window. Defaults use protein Cα atoms, fit on protein Cα, and make periodic molecules whole. Choose another atom group or a named selection for a ligand, local region or other molecular component. Alignment needs three noncollinear atoms; use protein alignment or turn fitting off for smaller selections.

Reference, first and last frame controls are one-based and inclusive. The analyzed stride starts at the first chosen frame. The reference can be outside that window. RMSF requires at least two analyzed frames. One-frame RMSD remains defined. Changing controls clears a previous result. Analysis settings and the expanded panel are saved with the workspace/project; results recalculate when requested.

Click the RMSD curve to seek the corresponding saved trajectory frame. Click an RMSF point to select its measured atoms in the viewer. **Export CSV** writes plotted values and labels; **Analysis record** writes JSON containing those values, the original request, exact measured/fit atom indices, frame indices, methods and warnings. Keep the project backup as well: the analysis record does not contain the source coordinates. The zero-based array indices in the record correspond to the preserved dataset atom order.

## Numerical contract

- Calculations use saved physical coordinates, independent of the camera, display alignment, rendering and animation interpolation. Coordinates are stored in nm; reported RMSD/RMSF are in Å.
- Fitting uses an equal-weight proper Kabsch rotation in float64, disallowing reflection. The requested fit group determines the transformation applied to all measured atoms. Degenerate collinear fits fail explicitly.
- RMSD is the square root of the mean squared atom displacement from the specified reference after the chosen fit. Atomic masses are not used as weights.
- RMSF uses population mean squared displacement about each atom's window mean. Welford accumulation avoids full-trajectory float64 temporaries. A grouped residue value is the square root of the mean atomic variance, not the arithmetic mean of atomic RMSFs.
- With valid periodic cells, “Make molecules whole” uses bonded connectivity and nearest molecular images before fitting. The anchor is the connected component containing the most fit atoms (or measured atoms with fitting off). This is structural imaging, not continuous diffusion unwrapping; inter-molecular jumps or changing nearest images require interpretation. A missing usable bonded topology produces an error. “Saved Cartesian coordinates” explicitly retains box crossings and reports a warning.
- Time axes use saved timestamps. Unknown timing is labeled as frame indices. Non-finite or non-increasing timestamps use saved frame indices with an explicit warning; elapsed time is not invented.
- RMSF groups use topology residue identity, so repeated residue numbers do not merge atoms from distinct residues. Labels can repeat when an input's chain/residue identifiers repeat; the record retains distinct atom groups.

Analytic fixtures test rigid translation/rotation, a proper fit of planar atoms, reflection rejection, known RMSF, residue aggregation, periodic boundary crossings, fit degeneracy, selection/window validation, stride mapping and invalid timestamps. A nondegenerate fit is compared independently against MDTraj RMSD. These establish the stated numerical behavior, not physical model accuracy or sampling convergence.

## Before preparation or simulation

Readiness shows blocking issues and planning estimates alongside the controls. Preparation uses the same chemistry validator as submission. Simulation uses the same input/model and resource guards as submission; calling the API directly cannot bypass these resource checks. No readiness request launches a job or creates solvent. Inspection may cache a missing public CCD reference.

The estimates include system atoms, saved frames, coordinate storage, working memory and output disk use. Unbuilt solvent boxes use conservative geometric estimates, so an over-limit estimate can block a system whose exact solvated count would be smaller. Preparing a system and constructing its explicit-water preview provides an actual atom count. Available disk space must cover the estimate plus a 256 MiB reserve. The existing 100,000-atom, 10,000-frame and 256 MiB retained-coordinate limits remain in force; multi-gigabyte streaming is future work.

Model details identify the engine, fixed-volume NVT preset, solvent and recorded preparation assumptions. Component charges are model assignments. The [chemistry matrix](CHEMISTRY_SUPPORT.md) describes their evidence and limitations. Readiness means the implemented checks allow an attempt; it does not predict stable dynamics, an experimental protonation state or an accurate model.

## During a run

Simulation diagnostics polls native potential energy, kinetic energy and kinetic temperature as they become available. It never substitutes the requested thermostat temperature for an observation. OpenMM temperature accounts for constraints and center-of-mass removal; GROMACS uses its recorded native temperature. Older runs without temperature report it as unavailable. Incomplete CSV rows wait for the next poll, and invalid/nonmonotonic observations are counted and omitted.

The display is bounded to 600 samples by retaining channel minima/maxima in bins, plus endpoints. This preserves sampled extrema while reducing plot size; it is not a complete time-series export. The chart CSV exports these displayed samples. **Files** contains the full recorded `energies.csv`, native logs and run records. Elapsed time and estimated remaining time use observed overall throughput, including setup; the estimate can change substantially during startup or resume.

A flat-looking energy or temperature trace does not establish equilibration or convergence. Use these plots to notice gross changes and investigate logs, with an appropriate scientific protocol for any subsequent conclusions. [Checkpoint recovery](CHECKPOINT_RECOVERY.md) explains compatible resume and partial files.

Primary API references: [MDTraj RMSD](https://mdtraj.readthedocs.io/en/latest/api/generated/mdtraj.rmsd.html), [MDTraj trajectory and molecular imaging](https://mdtraj.readthedocs.io/en/latest/api/generated/mdtraj.Trajectory.html), and [OpenMM StateDataReporter](https://docs.openmm.org/latest/api-python/generated/openmm.app.statedatareporter.StateDataReporter.html). DynaMol's exact choices above, including residue RMSF aggregation and timestamp fallback, are its own explicit implementation contract.
