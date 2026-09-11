# Self-contained DynaMol distribution

This is a packaging proposal. The current source launcher is not a portable installer: it requires uv and Node, builds the frontend, and creates a local Python environment. OpenMM is installed with the Python dependencies; GROMACS is discovered separately. No standalone release has been built or validated yet.

## Recommended first release

Build an offline Apple Silicon macOS application, distributed in a signed and notarized DMG. Opening DynaMol would start its local service and open the existing browser interface. A native window can be added later without changing the engine distribution.

Include inside the application:

- The prebuilt frontend, fonts and example data. Node is a build dependency only.
- A private Python interpreter and pinned production dependencies, including OpenMM, PDBFixer, RDKit and MDTraj.
- A private CPU GROMACS installation, including dependent libraries and force-field data, with paths that remain valid after relocation.
- Dependency license notices and required redistribution materials.

Recipients would not need to install Python, Node, uv, Conda, OpenMM or GROMACS. Local preparation, viewing and simulation would work offline. Database fetching would still require a connection.

Use the bundled engines by default. Keep external engine paths as an advanced override. Write uploads, jobs and logs to the user's application-data directory, such as `~/Library/Application Support/DynaMol`, outside the application bundle.

## Implementation constraints

Do not copy the current `.venv` and `.gromacs` directories into a ZIP: their Python symlink and GROMACS wrapper contain paths specific to the development machine. GROMACS supports relocatable installation trees, but libraries and data paths must be packaged together. [GROMACS relocation documentation](https://manual.gromacs.org/current/dev-manual/relocatable-binaries.html)

Preserve a real private Python executable initially: background workers currently launch with `sys.executable -m backend.worker` or `backend.preparation_worker`. A frozen executable requires a worker dispatcher and explicit inclusion of native libraries, package metadata and scientific resources. PyInstaller can bundle Python, but builds are specific to their operating system. [PyInstaller operating mode](https://www.pyinstaller.org/en/stable/operating-mode.html)

Ship CPU support first. Optional GPU support requires compatible platform libraries and vendor drivers; drivers cannot be assumed on a recipient's machine. [OpenMM installation documentation](https://docs.openmm.org/latest/userguide/application/01_getting_started.html)

Build and test other OS/architecture packages separately. Windows also needs replacement of the current Unix process-control code and a validated GROMACS distribution strategy.

## Release validation

Test the relocated application on a clean machine with no development tools or MD engines. Exercise structure loading, protein preparation, both simulation engines, cancellation and reopening a running job. Verify resource discovery, writable data storage, offline operation, license materials and the signed installer.

A smaller bootstrap installer could download verified engine bundles on first use, with installation progress and integrity checks. That is an alternative distribution mode, not a fully offline package.
