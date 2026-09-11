# DynaMol macOS prototype package

The local Apple Silicon application bundles its own Python 3.12 interpreter, Python dependencies (including OpenMM), prebuilt web interface, GROMACS 2025.4, and AmberTools 24.8. The recipient does not need Python, Node, uv, Conda, OpenMM, GROMACS, or AmberTools installed. Opening the app launches DynaMol in the default browser and adds a DynaMol menu in the macOS menu bar. This is a local prototype; it is not Apple Developer signed or notarized.

Unzip `DynaMol-0.1.0-macos-arm64.zip`, put `DynaMol.app` in a writable folder such as Applications, and open it. First launch shows setup progress while unpacking the included engines. It does not install global software or download engines. Data, jobs, engine prefixes and logs live in `~/Library/Application Support/DynaMol`; deleting or replacing the app does not delete that directory. The server binds an available localhost port, so it can coexist with the source-development server. Reopening the app reuses its running service. Closing the browser leaves simulations running; use the menu-bar item to reopen the workspace.

The runtime requires Apple Silicon and macOS 14 or newer (the bundled NumPy and SciPy wheels require macOS 14); platform-specific wheels and native engines cannot be reused for Windows or Linux. The app includes only the public Ubiquitin examples, not the developer's uploaded structures or jobs. Fetching new database structures still needs internet. Native ligand preparation needs either cached CCD chemistry or an explicit supplied graph; arbitrary uncached ligands are not promised to work offline. Prepared complex simulation remains OpenMM with explicit TIP3P.

The updated workbench includes durable scene restoration, Projects and portable backups, named selections, readiness estimates, live native diagnostics and RMSD/RMSF. Quit offers a choice to leave active jobs running or stop them first; compatible dynamics can resume from a saved production checkpoint. Preparation and solvent construction must be restarted if stopped. See the [release additions](../docs/RELEASE_READINESS.md) for current functional and chemistry coverage. Guides, cited audit evidence and frontend/backend test sources are included with this prototype.

## Build

On the prepared source checkout, install the build-only packing tool:

```sh
uv venv --python .venv/bin/python .tools/packaging
uv pip install --python .tools/packaging/bin/python conda-pack==0.8.1 setuptools==80.9.0
.venv/bin/python scripts/package_macos.py --prepare-runtimes
cd frontend && npm run build && cd ..
.venv/bin/python scripts/package_macos.py
```

The builder copies the real standalone Python distribution underlying the venv and overlays installed site packages. It does not copy a venv's external Python symlink or `pyvenv.cfg`. [python-build-standalone](https://github.com/astral-sh/python-build-standalone) supplies the interpreter distribution. Native engines are archived using [conda-pack](https://conda.github.io/conda-pack/) and relocated once into a versioned user prefix on first launch. No package manager is called at runtime. Prefixes are keyed by archive hashes; incomplete extraction is retried, and changed engine versions receive separate prefixes. Runtime archives are verified with SHA-256 before extraction.

Build outputs are under `build/releases`, with the ZIP SHA-256 next to the archive. `Contents/Resources/manifest.json` records included application source and engine archive hashes. The bundled app source, Python distribution metadata/licenses, native package license files, and available native build recipes/source URLs are included under Resources. This preserves evidence of what was shipped.

## Release limits

The launcher receives an ad-hoc signature; there is no Apple Developer certificate or notarization. A downloaded, quarantined copy may be blocked by Gatekeeper. This local prototype does not establish clean-machine compatibility or permission to bypass that protection. A public release needs signing/notarization and testing on a separate clean Mac.

`Notices/DEPENDENCY_INVENTORY.json` records the actual dependency versions, licenses and collected materials. All corresponding source archives required for redistribution of every copyleft binary have not yet been assembled or legally reviewed. Package recipes/source URLs and license texts are included, but this is not a claim that public binary redistribution requirements are complete. No package is published or uploaded by the builder.

## Validation performed

The [relocated runtime check](validation/relocated-runtime.json) passed using the copied Python distribution and separately unpacked engine archives from paths containing spaces. With a sanitized system-only PATH, it loaded OpenMM's force-field data and ran CPU minimization/integration, generated ethanol GAFF2/AM1-BCC parameters with the native AmberTools pipeline, and ran GROMACS Ubiquitin topology generation plus five minimization steps. Dynamic-loader traces for Python, SQM, and GROMACS showed only the relocated bundled libraries and macOS system libraries. Database/network functions were not used; proxy endpoints were disabled for native subprocesses. This is stronger than import/help probes, but remains a test on the development Mac rather than a fresh computer.

The standalone validation entry point is `packaging/macos/validate_runtime.py`. Run it with the bundled interpreter and Resources directory after building; it creates a separate validation workspace. The [application check](validation/packaged-app.json) also passed: startup progress, both engines, service reuse, authenticated shutdown/restart and data retention. Actual background workers completed OpenMM Ubiquitin (1,231 atoms) and explicit-water GROMACS (26,842 atoms), with six trajectory frames each. The [four browser zoom tests](validation/packaged-ui-tests.txt) passed against the packaged service on its assigned port. No convergence claim follows from these short functional runs.

A [native LaunchServices check](validation/native-launch.json) also started the relocated `.app` through macOS, confirmed the Swift launcher → bundled Python → local API process chain, and verified both engines plus same-origin frontend asset access. The generated app/test data were placed outside protected Documents folders, with data under a dedicated Application Support directory. Native menu clicks and automatic browser opening were not verified: the Mac was locked, so computer-use inspection could not proceed; browser opening was disabled in that isolated check. No lock or OS permission protection was bypassed. Those native interactions remain release-check items.
