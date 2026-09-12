# Build the DynaMol 0.1.1 Mac replacement

The 0.1.0 installer was withdrawn for a packaging signature error. These are
the replacement build instructions; downloaded-app opening validation is still
pending, and no replacement download is announced here.

The release format is an unsigned, self-contained **Apple Silicon / macOS 14+**
DMG: open it, drag DynaMol to Applications, then open the app. Python 3.12,
OpenMM, GROMACS 2025.4, AmberTools 24.8 and the interface are included.
Recipients do not install development tools or either MD engine separately.
See [installation, data locations and limits](../docs/PACKAGING.md).

## Build

These commands are for the prepared developer checkout. They require its
scientific/native runtimes, frontend tools, the macOS build toolchain and the
build-only packing environment:

```sh
uv venv --python .venv/bin/python .tools/packaging
uv pip install --python .tools/packaging/bin/python conda-pack==0.8.1 setuptools==80.9.0 ds_store==1.3.3 mac_alias==2.2.3
.venv/bin/python scripts/package_macos.py --prepare-runtimes
cd frontend
npm run build
cd ..
.venv/bin/python scripts/package_macos.py --output-dir build/releases
DYNAMOL_BUILT_APP=$(.venv/bin/python -c 'import json; print(json.load(open("build/releases/DynaMol.app.signature.json"))["bundle"])')
.tools/packaging/bin/python scripts/create_macos_dmg.py \
  --app "$DYNAMOL_BUILT_APP" --output build/releases/DynaMol-0.1.1-macos-arm64.dmg
```

The sealed `DynaMol.app` remains in a private system temporary directory recorded
in `DynaMol.app.signature.json`; do not assume it is under `build/releases`.
The requested output directory receives the ZIP, signature report and checksum.
The DMG command writes its finished DMG, checksum and validation report beside
the explicit `--output` path. Retain the staged app until acceptance is complete.
The DMG builder preserves its input app, uses built-in `hdiutil`/`ditto`, creates
an Applications shortcut, and checks a read-only mount and an isolated copy.
It refuses to overwrite an existing DMG; retain earlier evidence or choose a
new `--output` path. It does not install into the real `/Applications` directory
or launch simulations. No builder uploads release assets.

The app builder copies the real [standalone Python distribution](https://github.com/astral-sh/python-build-standalone)
underlying the venv and overlays site packages, avoiding an external interpreter
symlink. [conda-pack](https://conda.github.io/conda-pack/) archives native engines
and relocates them once into hash-versioned private prefixes. SHA-256 checks
precede extraction. No package manager runs at recipient launch.

Only public examples are seeded; user uploads, jobs and projects are excluded.
The app stores data under `~/Library/Application Support/DynaMol`, uses an
allocated localhost port and reuses a running service. The packaged atom limit
is 100,000; the source backend's optional larger profile is not forwarded by
the launcher and has no packaged UI toggle.

## Validate the payload and final container

The builder preserves valid inner signatures, repairs invalid inner signatures
when necessary, then seals the completed outer app last. A strict recursive
signature check must pass before packing. Verify the extracted ZIP, staged and
mounted DMG, and installed copy as well; the DMG builder enforces its copy checks.
Runtime bytecode writes into the sealed app are disabled. The signature report
records any repaired native files and their before/after hashes.

Use `packaging/macos/validate_app.py` and `validate_runtime.py` with the extracted
app's private interpreter, Resources directory and a fresh isolated home. For
example, after extracting the candidate ZIP into `/tmp/DynaMol Check`:

```sh
"/tmp/DynaMol Check/DynaMol.app/Contents/Resources/python/bin/python3.12" -I -B \
  "/tmp/DynaMol Check/DynaMol.app/Contents/Resources/app/packaging/macos/validate_app.py" \
  "/tmp/DynaMol Check/DynaMol.app/Contents/Resources" "/tmp/DynaMol Check/Private Home"
"/tmp/DynaMol Check/DynaMol.app/Contents/Resources/python/bin/python3.12" -I -B \
  "/tmp/DynaMol Check/DynaMol.app/Contents/Resources/app/packaging/macos/validate_runtime.py" \
  "/tmp/DynaMol Check/DynaMol.app/Contents/Resources" --home "/tmp/DynaMol Check/Private Home"
```

Release-side `first-release-acceptance.json` records the exact ZIP's hash/build
ID, both native worker completions, private runtime/AmberTools checks,
Host/Origin controls, workspace/project persistence across restart, native
launcher/default-browser startup and seven Chrome pinch/zoom tests. Native
Quit-menu clicking remains untested; authenticated shutdown passed.

`isolated-installation-validation.json` records the additional test denying
reads from the developer's home and common external package prefixes. Parent
and inherited child sentinel reads failed as required; both native workers and
all runtime/library checks passed. macOS's deprecated `sandbox-exec` cannot run
setuid `ps`, so an external copy of the test harness read saved worker status
while jobs ran and checked terminal results through the API. The original
failure and exact diff are saved. Neither app code nor file-denial rules changed.
This is test isolation, not a shipping App Sandbox or clean-Mac claim.

For the replacement DMG, keep the `.dmg.validation.json`, signature report,
code/data comparison and copied-app launch record with its checksum. If native
signatures were repaired, full binary hashes change: verify their code/data
identity and record signature differences separately. Link unchanged scientific
payload to the earlier exact-ZIP acceptance without claiming simulations were
repeated. Historical runtime checks did not validate the withdrawn outer seal;
the replacement also requires its own downloaded/quarantined opening check.
A rebuild that changes executable payload needs relevant new acceptance.
Acceptance files generated after packing live beside the release, not
necessarily inside the app being tested.

The [five-system audit](../docs/audit/five-system-pressure-test/REPORT.md) covers
separate 100 ps source-app runs and real trajectory browser checks. Its modeled
loop failure, water-only membrane-protein setup and larger GPCR profile remain
explicit. See [release scope](../docs/RELEASE_READINESS.md).

## Unsigned release and source materials

The workflow ad-hoc signs the complete app; it requires no paid Apple Developer
enrollment and performs no Developer ID signing or notarization. Signature
integrity and Gatekeeper trust are separate, as described in
[Apple TN2206](https://developer.apple.com/library/archive/technotes/tn2206/_index.html).
Downloaded copies may be blocked by Gatekeeper; link recipients
to [Apple's per-app opening instructions](https://support.apple.com/en-us/102445).
Do not ask them to disable system-wide protections.

`Notices/DEPENDENCY_INVENTORY.json` records bundled package versions and license
materials. The [source-materials companion](source-materials/README.md) adds
upstream archives, shipped recipes/patches and per-package provenance mappings.
The 0.1.1 candidate reuses the 0.1.0 collection because upstream component
versions are unchanged. Provide it beside the binary release with its index
and checksum. Acceptance must link both build IDs, unchanged engine archives
and native code/data; signature-only byte changes are recorded separately.
It is not required for app use.
See [THIRD_PARTY.md](../THIRD_PARTY.md) for the supplied scope and remaining
compiler-runtime provenance limitations; the collection makes no legal
certification or bit-identical-rebuild claim.
