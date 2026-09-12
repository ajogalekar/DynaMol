# DynaMol 0.1 for macOS

DynaMol is a free early release for individual scientists to prepare supported
molecules, run local simulations and inspect trajectories. The Mac download is
an unsigned, self-contained DMG for **Apple Silicon with macOS 14 or newer**.
The GitHub source and third-party source-materials companion are separate from
the app download; neither is needed to use it.

## Install and open

1. Download `DynaMol-0.1.0-macos-arm64.dmg` and open it.
2. Drag **DynaMol** onto the **Applications** shortcut in the disk-image window.
3. Eject the disk image and open DynaMol from Applications.
4. First launch shows progress while unpacking the included engines into your
   private application data directory. The interface then opens in your default
   browser, using an available localhost port.

Python, OpenMM, GROMACS, AmberTools and the interface are included. You do not
need a separate Python, Node, uv, Conda or MD-engine installation. First-run
engine setup uses the files already in the download; it does not download them
or install global software.

**This release is not Developer ID signed or notarized.** If macOS blocks the
first opening because the developer cannot be verified, follow
[Apple's per-app opening guidance](https://support.apple.com/en-us/102445): after
trying to open it, use System Settings → Privacy & Security → Open Anyway if
you trust the download. Managed Macs may restrict this option. No Apple
Developer enrollment is part of the installation or release path.

Load a structure, review preparation and readiness, then start a simulation.
Progress and diagnostics appear on its card. A completed trajectory opens in
Explore automatically while you are following that run; otherwise use its
Open trajectory button. Save named work in Projects and use portable backups
when needed.

## Your data and background work

Structures, projects, jobs, private engine prefixes and logs live in
`~/Library/Application Support/DynaMol`, outside the app bundle. Replacing or
deleting the app does not delete this directory. Only the public Ubiquitin
examples are seeded; developer uploads and pressure-test jobs are not bundled.

Closing the browser leaves background work running. The DynaMol menu-bar item
reopens the interface. Quit offers a choice to leave active jobs running or stop
them first. Compatible dynamics can resume from a saved production checkpoint;
preparation and solvent construction must be restarted if stopped.
Geometry-rejected jobs require repaired input and a new simulation.
See [checkpoint recovery](CHECKPOINT_RECOVERY.md).

## Scope and limits

- CPU execution is the checked path. Intel Mac, Windows, Linux and additional
  GPU support require their own builds and validation.
- Fetching new structures needs internet. Ligand preparation needs cached CCD
  chemistry or a sufficient supplied graph; arbitrary uncached ligands are not
  promised to work offline.
- Supported prepared complexes and modified residues use OpenMM with explicit
  TIP3P. The GROMACS workflow currently handles compatible standard-protein
  inputs; packaging does not add complex parameter export or membrane setup.
- The packaged limit is **100,000 atoms**, with no UI toggle. The source service
  permits an explicit `DYNAMOL_MAX_ATOMS=200000 ./start.sh` profile; the packaged
  launcher does not forward that override. The 145,913-atom GPCR benchmark used
  the larger source profile and exceeds the package's limit.
- Successful preparation or a short finite trajectory does not establish
  protonation, metal coordination, ligand pose accuracy or convergence.

## Installation evidence

The release-side `first-release-acceptance.json` identifies the exact ZIP, build
ID and checksums used for full bundled-app acceptance. Both native MD workers
completed six-frame installation checks; native AmberTools parameterization,
local Host/Origin checks, restart with saved datasets/projects/workspace, and
101-frame structural RMSD passed. Seven real Chrome pinch/zoom controls passed.
Native app launch and automatic default-browser opening also passed. Native
Quit-menu clicking was not verified; authenticated shutdown and restart were.

`isolated-installation-validation.json` adds a stronger dependency check. The
same extracted payload ran with reads denied from the developer's entire home,
including the checkout, venv and native engine installations, plus common
external package prefixes. Both a bundled-Python parent and inherited child
received permission errors on existing developer files. Both MD workers,
AmberTools and dynamic-library checks passed with these restrictions intact.

That isolation test required one **external test-harness adjustment**: macOS
`sandbox-exec` refuses to execute the system's setuid `ps`, which the app uses
when polling active jobs. The harness waited on saved worker status, then
verified completed jobs through the normal API. The original failure and exact
harness diff are retained; app code and filesystem-denial rules were unchanged.
Normal API progress had already passed outside this test sandbox.

The final DMG's checksum and `.dmg.validation.json` identify its mounted/copy
checks. Release acceptance must link its native-payload byte parity to the
validated ZIP and record the copied app's launch check. A documentation/notices
or container-only update does not mean the full native simulations were rerun.
These are tests on the development Mac, not an independent clean-Mac,
minimum-OS, physical Safari-trackpad or Gatekeeper certification.

The separate [five-system software pressure test](audit/five-system-pressure-test/REPORT.md)
completed three OpenMM and two GROMACS 100 ps source-app runs. It preserves a
failed reconstructed GPCR loop and states the water-only membrane-protein and
larger-resource-profile limits. Those results are distinct from package
installation checks.

## Build and third-party materials

See [packaging instructions](../packaging/README.md). Release assets have
SHA-256 companion files; the app's Resources manifest records source and engine
archive hashes. Keep acceptance records alongside the assets they identify.

`Notices/DEPENDENCY_INVENTORY.json` includes dependency notices, metadata and
available recipes. The separate [source-materials companion](../packaging/source-materials/README.md)
adds verified upstream archives, recipe/patch mappings and explicitly recorded
provenance limitations. It is available for inspection and rebuilding and is
not needed to run DynaMol. See [third-party notices](../THIRD_PARTY.md) for its
scope, compiler-runtime exceptions and AmberTools component-specific terms.
