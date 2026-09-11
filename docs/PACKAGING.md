# Self-contained DynaMol distribution

A local Apple Silicon macOS prototype has been built as `DynaMol.app` and
`DynaMol-0.1.0-macos-arm64.zip` under `build/releases/`. It bundles the engines
and application runtime: recipients do not need Python, Node, uv, Conda,
OpenMM, GROMACS or AmberTools installed. This is a local prototype, not a public
release. It has no Apple Developer signature or notarization, and clean-Mac
compatibility and the full redistribution audit remain unverified.

## What the package contains

The app includes a standalone Python 3.12 interpreter, the installed Python
scientific dependencies including OpenMM, the built frontend and fonts,
GROMACS 2025.4, and AmberTools 24.8 with the native programs, libraries and data
needed by the implemented preparation pipeline. Only the public Ubiquitin
examples are seeded; developer uploads and jobs are excluded.

The supported target is **Apple Silicon with macOS 14 or newer**. The minimum
OS version follows the bundled NumPy and SciPy wheels. Intel macOS, Windows and
Linux require separately built and tested distributions. CPU execution has been
checked; bundling the app does not establish additional GPU compatibility.

Unzip the archive, place `DynaMol.app` in a writable location such as Applications,
and open it. The launcher opens the local interface in the default browser and
provides a DynaMol menu in the macOS menu bar. First launch shows progress while
unpacking the included engine archives. It does not download engines, invoke a
package manager, or install global software. Engine prefixes are versioned by
archive hashes, and SHA-256 verification precedes extraction.

Uploads, jobs, engine prefixes and logs live in
`~/Library/Application Support/DynaMol`, outside the app bundle. Replacing the
app does not remove that directory. The service selects an available localhost
port and can coexist with the source-development server. Reopening the app
reuses its service; closing the browser leaves background work running.

Database fetching still needs internet. Ligand preparation needs cached CCD
chemistry or a sufficient supplied graph; arbitrary uncached ligands are not
promised to work offline. Prepared complexes and supported modified residues
currently require OpenMM with explicit TIP3P water. Packaging does not add
GROMACS complex parameter export or change the 100,000-atom alpha limit.

## Build and runtime design

See [the package build instructions](../packaging/README.md) and
[`scripts/package_macos.py`](../scripts/package_macos.py) for the current
commands. Builds require the prepared source checkout and its scientific and
native runtimes, frontend build tools, the macOS build toolchain, and the
build-only `conda-pack` environment. These requirements belong to the builder,
not the recipient.

The builder copies the real interpreter distribution underlying the development
venv and overlays its installed site packages. It does not ship an external
venv Python symlink or `pyvenv.cfg`. The interpreter comes from
[python-build-standalone](https://github.com/astral-sh/python-build-standalone).
Native engines are packed with [conda-pack](https://conda.github.io/conda-pack/)
and relocated into a user prefix once at launch. This preserves a real private
`sys.executable` for the existing background worker processes. GROMACS libraries
and data must travel with its executable; see the
[GROMACS relocation documentation](https://manual.gromacs.org/current/dev-manual/relocatable-binaries.html).

The ZIP has a neighboring SHA-256 file. `Contents/Resources/manifest.json`
records application source and engine archive hashes. Application source,
Python metadata/licenses, native package notices and available package recipes
and source URLs accompany the prototype. `Notices/DEPENDENCY_INVENTORY.json`
records the collected versions, licenses and materials.

The source checkout remains a separate way to run DynaMol. Its `start.sh`
requires uv and Node/npm (Node 22 or newer is documented), creates the Python
environment, builds the frontend and starts the local service. OpenMM comes
from its Python dependencies; the source workflow discovers GROMACS separately
and uses `scripts/install_ligand_tools.sh` for the private AmberTools runtime.
Those source prerequisites do not apply to the built app.

## Validation and remaining release work

The [relocated runtime check](../packaging/validation/relocated-runtime.json)
passed using the copied interpreter and unpacked engine archives in paths
containing spaces, with a system-only PATH. It loaded OpenMM force-field data,
ran CPU minimization and integration, generated native GAFF2/AM1-BCC ethanol
parameters, and ran GROMACS Ubiquitin topology generation and five minimization
steps. Dynamic-loader traces for Python, SQM and GROMACS used the relocated
libraries and macOS system libraries. No network operations were used, and
proxy endpoints were disabled; this was not an OS-level network sandbox.

The [application-launch check](../packaging/validation/packaged-app.json) passed:
setup progress, both bundled engine availability probes, service reuse,
authenticated shutdown and restart with datasets and jobs retained. Real
background workers launched with the private interpreter completed a 1,231-atom
OpenMM Ubiquitin run and a 26,842-atom explicit-water GROMACS run, each producing
six frames. Four browser pinch/zoom regression tests also passed against the
packaged service on its dynamically assigned port; the exact local origin is
configured without accepting arbitrary foreign origins. These are short
functional checks, not convergence or accuracy validation. The entry points are
[`validate_runtime.py`](../packaging/macos/validate_runtime.py) and
[`validate_app.py`](../packaging/macos/validate_app.py). The checks used relocated
paths containing spaces on the development Mac, not a separate clean computer.

The [native LaunchServices check](../packaging/validation/native-launch.json)
started the relocated `.app` and verified the Swift launcher → private Python →
API process chain, both engines and frontend asset access. Native menu clicks
remain untested because the Mac was locked; automatic browser opening was
disabled in this isolated check. No lock or OS file-access protection was
bypassed. These native interactions remain release validation items.

The launcher has an ad-hoc signature only. A quarantined downloaded copy may be
blocked by Gatekeeper; Apple Developer signing, notarization and a separate
clean-Mac test are still required before a supported public release.

Corresponding source archives required for redistribution of all copyleft
binaries have not yet been fully assembled or legally reviewed. Included license
texts, recipes and source URLs preserve provenance but do not establish that
all public binary redistribution requirements are complete. AmberTools has
component-specific licenses; its installed license identifies public-domain
force-field data under `dat/leap`, separately from program licenses. This is not
a blanket redistribution license for proprietary Amber/pmemd. See
[AmberTools](https://ambermd.org/AmberTools.php) and
[DynaMol's third-party notices](../THIRD_PARTY.md).

No artifact has been publicly published or uploaded by this packaging workflow.
