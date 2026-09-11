# Engine/preparation compatibility review

Scope: source and documentation review of DynaMol's engine-selection workflow, 2026-09-11. No molecular structures were newly prepared, no MD was run, and no scientific algorithm was changed by this review. This is language-model evidence synthesis and an automated PatAgent checklist, not independent expert certification or Pat Walters endorsement.

## Finding

Engine-first selection is appropriate for the current DynaMol adapters. A repaired molecular structure is not inherently tied to one MD engine, and supported complete simulation inputs can be transferred. For example, OpenMM documents reading GROMACS coordinate and topology files. That does not mean a coordinate PDB alone carries the force-field parameters or guarantees state-preserving transfer. [OpenMM: Using Gromacs Files](https://docs.openmm.org/latest/userguide/application/02_running_sims.html#using-gromacs-files).

DynaMol currently uses different preparation paths and protein force fields. The OpenMM path preserves its recorded protonation state and exact prepared topology, and uses ff14SB with validated supplementary ligand/modified-residue parameter files when applicable. The GROMACS worker calls `pdb2gmx` with Amber99SB-ILDN, TIP3P and `-ignh`, rebuilding hydrogens through its own templates and defaults before constructing its own solvent box. GROMACS documents both topology generation and ignoring input hydrogens with this option. [GROMACS: pdb2gmx](https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html).

## Evidence and UI implications

- `backend/jobs.py:148`: simulation validation rejects DynaMol preparation **or solvation** metadata for GROMACS (lines 158–159). It also rejects unverified ligand parameters, unsupported residue templates and unresolved backbone gaps. `backend/readiness.py:35` invokes that same validator; these guards need no algorithm change for an engine-first UI.
- `backend/prepared_system.py:104`: verified supplementary parameter files and ff14SB are loaded into an OpenMM ForceField. `backend/worker.py:171`: prepared hydrogen state is retained. Prepared complexes/modified residues/ions have explicit-solvent restrictions; visibility or preparation alone is not a universal simulation compatibility claim.
- `backend/worker.py:368`: native GROMACS setup runs at simulation start. It preserves checked input heavy-atom counts, adds its own water and neutralizing ions, and does not parameterize ligands or repair missing loops. A UI must not describe this as the OpenMM interactive preparation/solvent-preview workflow.
- Show the engine selector before structure preparation. For GROMACS, explain the native setup on start and require a compatible complete standard-protein source. OpenMM-prepared and solvent-preview datasets need a clear incompatibility notice or a route back to OpenMM/original source. Engine switching must not silently strip preparation metadata, molecular atoms or recorded state.
- General transfer language should describe topology, parameter and molecular-state preservation. Claiming that engines can *never* share prepared systems is incorrect; claiming smooth conversion in DynaMol is unsupported.

These are code-supported workflow conclusions. This review does not establish universal chemical coverage, native-engine equivalence, energy/force parity, pKa accuracy, convergence or state-conversion validity. Root-agent browser checks separately validate the changed controls.

## Evidence-linked review

The existing cached Practical Cheminformatics corpus was searched for molecular dynamics/preparation/protonation/reproducibility and tautomers. The cached manifest contains 74 Blogger and 17 GitHub Pages records; no claim of fresh exhaustive corpus coverage is made. The relevant Pat-authored post, [The Trouble With Tautomers](https://patwalters.github.io/The-Trouble-With-Tautomers/), discusses how tautomer representation changes chemical inputs. Applying the principle of recording and preserving chosen molecular representations to engine routing is this review's inference, not an MD-engine compatibility claim by the author. Official engine documentation and current project code provide the direct compatibility evidence.

The lightweight postflight checklist records source checksums and environment. Its unresolved uncertainty gate is retained because this UI/source review produced no scientific replicate distributions or uncertainty intervals; it is not a scientific validation pass. No expensive-compute preflight was required or performed.

## Final frontend review

Read-only review of the final `SimulationPanel.tsx`, `StructureWorkbench.tsx`, `ReadinessPanel.tsx` and `App.tsx` change confirms that engine cards precede structure loading/preparation; GROMACS hides OpenMM preparation/pH/water-preview controls; the native setup and lack of a separate GROMACS preview are described explicitly. Both prepared and solvation-only metadata activate the incompatibility notice and disable Start. Recovery explicitly switches to OpenMM without discarding dataset state. The selected engine is held by App across Simulation Studio closure/reopening (not claimed to persist across a full page refresh).

One control-state issue was identified and resolved before completion: clicking the already-selected engine invalidated `runReady` without changing the readiness request key. The final handler returns early for the current engine, preventing a Ready-to-run panel with a disabled Start button. Switching engines continues to clear readiness so the new engine is checked. No further actionable bypass or misleading engine/pH/solvent wording was found in this bounded code review. Browser behavioral tests are owned and reported separately by the implementation/test agents.
