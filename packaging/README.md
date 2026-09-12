# Build the DynaMol 0.1 Mac archive

The first-release target is a GitHub source repository and an unsigned,
self-contained **Apple Silicon / macOS 14+** archive for individual scientists.
The app bundles Python 3.12, scientific dependencies including OpenMM, the built
interface, GROMACS 2025.4 and AmberTools 24.8. Recipients do not install those tools.
See [opening the app, data locations and limits](../docs/PACKAGING.md).

Use the matching acceptance record beside the archive; prior checks apply only
to their recorded build IDs and source hashes. No artifact is uploaded by the
builder.

## Build

These are requirements for the prepared developer checkout, not for recipients:
scientific/native runtimes, frontend build tools, the macOS build toolchain and
the build-only packing environment.

```sh
uv venv --python .venv/bin/python .tools/packaging
uv pip install --python .tools/packaging/bin/python conda-pack==0.8.1 setuptools==80.9.0
.venv/bin/python scripts/package_macos.py --prepare-runtimes
cd frontend
npm run build
cd ..
.venv/bin/python scripts/package_macos.py
```

Outputs are `build/releases/DynaMol.app`,
`build/releases/DynaMol-0.1.0-macos-arm64.zip` and its SHA-256 file.
`Contents/Resources/manifest.json` records included source and engine hashes.
Only the public examples are seeded; user uploads, jobs and projects are excluded.

The builder copies the real [standalone Python distribution](https://github.com/astral-sh/python-build-standalone)
underlying the venv and overlays site packages, rather than shipping an external
Python symlink. [conda-pack](https://conda.github.io/conda-pack/) archives native
engines and relocates them once into hash-versioned private prefixes. SHA-256
verification precedes extraction. No package manager runs at recipient launch.

The service uses an allocated localhost port and stores data under
`~/Library/Application Support/DynaMol`. It reuses a running service. Its private
environment defaults to 100,000 atoms: the source backend's opt-in 200,000 limit
is currently not forwarded by the packaged launcher and has no packaged UI toggle.

## Validate the exact archive

Use `packaging/macos/validate_runtime.py` and `validate_app.py` with the bundled
interpreter, Resources directory and a fresh isolated data location. Match the
result's build ID and source hashes to the archive being offered. Keep the
acceptance record alongside the release; a record generated after packing is
not automatically inside that same archive.

Existing [runtime checks](validation/relocated-runtime.json) exercised native
OpenMM, GROMACS and ligand parameterization from relocated paths. The earlier
[workbench package](validation/workbench-app.json) ran both MD workers and
verified readiness, RMSD, and scene/named-project persistence on a new port;
[archive checks](validation/workbench-archive.json) verified its recorded files.
These passed on the developer's Mac with separate data and a restricted PATH,
not on an independent clean Mac. [Native launch](validation/native-launch.json)
was checked, but menu clicks and automatic browser opening were not verified on
the locked desktop. Those gaps are stated limits of the evidence.

The [five-system audit](../docs/audit/five-system-pressure-test/REPORT.md) records
five completed 100 ps source-app runs and passing real-browser checks. It also
preserves a failed modeled loop, water-only membrane-protein tests, and the
larger GPCR resource profile. It does not certify the rebuilt archive or physical
accuracy. Current feature/evidence scope is in [release readiness](../docs/RELEASE_READINESS.md).

## Distribution status

The launcher currently receives an ad-hoc signature, without Developer ID or
notarization. An unsigned GitHub beta is a possible distribution choice;
downloaded copies may encounter Gatekeeper warnings or blocking. Signing and
notarization can reduce that friction, and are not an absolute prerequisite for
a GitHub release. [Apple explains these checks](https://support.apple.com/en-us/102445).

`Notices/DEPENDENCY_INVENTORY.json` records collected versions, licenses and
materials. Full corresponding-source materials for bundled copyleft binaries
have not yet been assembled or reviewed. Recipes, URLs and license texts do not
by themselves establish that public binary redistribution requirements are
complete. Preserve this unresolved status separately from signing; see
[THIRD_PARTY.md](../THIRD_PARTY.md). No public-ready claim or upload is made here.
