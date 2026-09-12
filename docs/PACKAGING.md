# DynaMol 0.1 for macOS

DynaMol 0.1 is an early release for individual scientists to prepare supported
molecules, run local simulations and inspect trajectories. The intended GitHub
distribution is the source repository plus an unsigned, self-contained Apple
Silicon Mac archive. No release has been uploaded by this workflow. Use the matching acceptance
record beside the archive; prior package checks apply only to their recorded
builds.

## Open the app

The target is **Apple Silicon with macOS 14 or newer**. The archive is named
`DynaMol-0.1.0-macos-arm64.zip`. It contains Python, OpenMM, GROMACS, AmberTools,
the web interface and their needed runtime files. Recipients do not need to
install Python, Node, uv, Conda or either MD engine separately.

1. Unzip the archive and place `DynaMol.app` in Applications or another folder
   you can write to.
2. Open the app. First launch displays progress while unpacking the included
   engines; it does not download them or install global software.
3. Use the interface in your default browser. The DynaMol menu-bar item reopens
   it. The service chooses an available localhost port.
4. Load a structure, review preparation/readiness, then start a simulation.
   Progress and diagnostics appear on its card. A completed trajectory opens
   automatically in Explore only while you are following that run; otherwise
   use its Open trajectory button. Save named work in Projects and use portable
   backups when needed.

The candidate has an **ad-hoc signature, not an Apple Developer ID signature or
notarization**. A downloaded copy may be warned about or blocked by Gatekeeper.
Signing and notarization are optional distribution improvements that reduce
installation friction; they are not prerequisites for hosting an unsigned beta
on GitHub. See [Apple's explanation](https://support.apple.com/en-us/102445).

## Your data and background work

Structures, projects, jobs, private engine prefixes and logs live in
`~/Library/Application Support/DynaMol`, outside the app bundle. Replacing or
deleting the app does not delete this directory. Only the public Ubiquitin
examples are seeded; developer uploads and pressure-test jobs are not bundled.

Closing the browser leaves background work running. Quit offers a choice to
leave active jobs running or stop them first. Compatible dynamics can resume
from a saved production checkpoint; preparation and solvent construction must
be restarted if stopped. Geometry-rejected jobs require repaired input and a
new simulation. See [checkpoint recovery](CHECKPOINT_RECOVERY.md).

## Scope and limits

- CPU execution is the checked path. Intel Mac, Windows, Linux and additional
  GPU support require their own builds and validation.
- Fetching new structures needs internet. Ligand preparation needs cached CCD
  chemistry or a sufficient supplied graph; arbitrary uncached ligands are not
  promised to work offline.
- Supported prepared complexes and modified residues use OpenMM with explicit
  TIP3P. The GROMACS workflow currently handles compatible standard-protein
  inputs; packaging does not add complex parameter export or membrane setup.
- The default limit is **100,000 atoms**. The source service permits an explicit
  opt-in up to 200,000 using `DYNAMOL_MAX_ATOMS=200000 ./start.sh`. The packaged
  launcher currently keeps the 100,000 default and does not forward that
  override; there is no packaged UI toggle. The 145,913-atom GPCR benchmark used
  the larger source-service profile, so that case exceeds this package's limit.
- Chemistry and size checks remain meaningful blockers. Successful preparation
  or a short finite trajectory does not establish protonation, metal
  coordination, ligand pose accuracy or convergence.

## What has been checked

The [five-system software pressure test](audit/five-system-pressure-test/REPORT.md)
completed five 100 ps cases with three OpenMM and two GROMACS runs; all five real
trajectory browser audits passed. A failed rebuilt GPCR loop is preserved in the
report, and the successful replacement required no internal loop construction.
Membrane proteins used explicitly declared water-only boxes. These source-app
results are not a claim that every case ran in the rebuilt package.

Earlier [relocated-runtime](../packaging/validation/relocated-runtime.json),
[package application](../packaging/validation/workbench-app.json) and
[archive](../packaging/validation/workbench-archive.json) checks passed on the
**development Mac** with isolated data, paths containing spaces and a restricted
PATH. They exercised both native workers and workspace/project persistence
across a changed port. The [native launcher check](../packaging/validation/native-launch.json)
verified its process chain; menu clicks and automatic browser opening were not
checked while that desktop was locked. An independent clean-Mac check remains
unperformed. New candidate results must match its build ID and archive hash.

## Build and redistribution records

See [packaging instructions](../packaging/README.md) for build commands. The ZIP
has a SHA-256 companion file; its Resources manifest records source and engine
archive hashes. Included license texts, dependency metadata, recipes and source
URLs are catalogued in `Notices/DEPENDENCY_INVENTORY.json`.

**Public binary redistribution materials remain unresolved:** all corresponding
source archives needed for the bundled copyleft components have not yet been
fully assembled or reviewed. Existing notices and URLs alone do not establish
completion. AmberTools has component-specific terms; its inclusion is not a
license for proprietary Amber/pmemd. See [third-party notices](../THIRD_PARTY.md).
This is separate from whether the app is Apple signed and must remain visible
when deciding which binary assets to attach to a GitHub release.
