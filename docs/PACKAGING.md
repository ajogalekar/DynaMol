# Self-contained DynaMol distribution

This is a packaging proposal. The current source launcher is not a portable installer: it requires uv and Node, builds the frontend, and creates a local Python environment. OpenMM is installed with the Python dependencies; GROMACS is discovered separately. Native ligand preparation additionally uses the private AmberTools environment installed by `scripts/install_ligand_tools.sh`. No standalone release or DMG has been built or validated yet.

## Recommended first release

Build an offline Apple Silicon macOS application, distributed in a signed and notarized DMG. Opening DynaMol would start its local service and open the existing browser interface. A native window can be added later without changing the engine distribution.

Include inside the application:

- The prebuilt frontend, fonts and example data. Node is a build dependency only.
- A private Python interpreter and pinned production dependencies, including OpenMM, PDBFixer, RDKit, MDTraj, ParmEd, Gemmi and Dimorphite-DL.
- A private, pinned AmberTools 24.8 runtime for Antechamber/SQM AM1-BCC, Parmchk2 and LEaP, including GAFF2, BCC and LEaP data and required libraries. The development installer writes a native package export; packaging still needs relocation and clean-machine validation.
- A private CPU GROMACS installation, including dependent libraries and force-field data, with paths that remain valid after relocation.
- A curated CCD cache for bundled examples and an application-data cache for later components, preserving reference bytes, source URLs, hashes and retrieval metadata.
- Dependency license notices and required redistribution materials, including the actual AmberTools environment's component and transitive licenses.

Recipients would not need to install Python, Node, uv, Conda, AmberTools, OpenMM or GROMACS. Local viewing, protein preparation and simulations with available models would work offline. New CCD-based ligand preparation requires a cached component or sufficient supplied graph; an unknown, uncached component cannot be promised offline support. Database fetching still requires a connection. Prepared complex archives contain their parameters and do not need another charge calculation or CCD fetch to simulate.

Use the bundled engines by default. Keep external engine paths as an advanced override. Write uploads, jobs and logs to the user's application-data directory, such as `~/Library/Application Support/DynaMol`, outside the application bundle.

The source override for ligand tools is `DYNAMOL_AMBERTOOLS`. A packaged launcher must resolve its private prefix and data paths from the installed application location. Users should not need to activate Conda. The implemented complex path uses GAFF2/AM1-BCC with explicit-water OpenMM; bundling GROMACS does not add GROMACS complex parameter export.

## Implementation constraints

Do not copy the current `.venv` and `.gromacs` directories into a ZIP: their Python symlink and GROMACS wrapper contain paths specific to the development machine. GROMACS supports relocatable installation trees, but libraries and data paths must be packaged together. [GROMACS relocation documentation](https://manual.gromacs.org/current/dev-manual/relocatable-binaries.html)

Apply the same rule to `.tools/ambertools`: use a supported relocatable package strategy and test Antechamber, SQM, Parmchk2 and LEaP from a path containing spaces. Include the GAFF2 and BCC data actually used. An OpenMM import or a native tool's help output does not establish that parameterization can find its data.

Collect notices and applicable source/redistribution materials for every binary and dependency shipped. AmberTools has several component licenses: its default is GPLv3, with exceptions and separately licensed third-party components, while force-field files under `dat/leap` are public-domain data. Do not include proprietary Amber/pmemd under an assumption that all Amber software is freely redistributable. This proposal is not a completed license audit. [AmberTools distribution](https://ambermd.org/AmberTools.php), [Amber project and force-field information](https://ambermd.org/index.php)

Preserve a real private Python executable initially: background workers currently launch with `sys.executable -m backend.worker` or `backend.preparation_worker`. A frozen executable requires a worker dispatcher and explicit inclusion of native libraries, package metadata and scientific resources. PyInstaller can bundle Python, but builds are specific to their operating system. [PyInstaller operating mode](https://www.pyinstaller.org/en/stable/operating-mode.html)

Ship CPU support first. Optional GPU support requires compatible platform libraries and vendor drivers; drivers cannot be assumed on a recipient's machine. [OpenMM installation documentation](https://docs.openmm.org/latest/userguide/application/01_getting_started.html)

Build and test other OS/architecture packages separately. Windows also needs replacement of the current Unix process-control code and a validated GROMACS distribution strategy.

## Release validation

Test the relocated application on a clean machine with no development tools or MD engines. Exercise structure loading, protein preparation, native complex preparation with a cached CCD ligand and a supplied graph, supported paths in both engines, cancellation and reopening a running job. Verify native-to-XML energy/force checks, parameter archive reuse through explicit-water OpenMM, resource discovery, writable storage, cache behavior with networking disabled, license materials and the signed installer. The 100,000-atom alpha limit remains independent of packaging.

A smaller bootstrap installer could download verified engine bundles on first use, with installation progress and integrity checks. That is an alternative distribution mode, not a fully offline package.
