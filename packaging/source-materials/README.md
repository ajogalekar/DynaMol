# DynaMol source-materials companion

This optional companion supports inspection and rebuilding of selected third-party components in the self-contained Apple Silicon app. Users of the app do not need it or a separate engine installation. The app DMG/ZIP is unchanged by source collection.

Run the collector from the DynaMol source checkout with its Python environment (the copied script under `materials/collector/` is an inspection reference):

```sh
.venv/bin/python scripts/collect_release_sources.py --bundle /path/to/DynaMol.app
.venv/bin/python scripts/collect_release_sources.py --verify-only
```

The collector uses PyYAML and Python's standard library. Output goes to `build/releases/source-materials`, which is gitignored. `--resume` retries interrupted downloads against the same app manifest; it does not silently switch to a newer package version. Each downloaded file must match the source checksum recorded by its shipped recipe, `uv.lock`, or official release checksum list. OpenMM is pinned to the Git revision embedded in the shipped `openmm/version.py`; its archive checksum is recorded on retrieval.

## Download parts and verification

The combined collection is distributed as two optional assets so each remains below the hosting service’s per-file limit:

1. `DynaMol-0.1.0-source-materials.tar.gz`, retained byte-for-byte as the base; SHA-256 `50c63be58fafcf63a7334a9d74a72b13a460d663c2faaf1efc38cff877365929`.
2. `DynaMol-0.1.1-loop-runtime-source-supplement.tar.gz`, containing new archives plus the current combined index, recipes, patches and notices.

Extract the base first, then the supplement into the same directory. Both use a `source-materials/` root; the supplement updates its index and guides. The old version in the base filename identifies reused source materials, not the version of the new application. Neither asset is required to run DynaMol. Keep both assets and their checksum records together.

```sh
mkdir dynamol-sources
tar -xzf DynaMol-0.1.0-source-materials.tar.gz -C dynamol-sources
tar -xzf DynaMol-0.1.1-loop-runtime-source-supplement.tar.gz -C dynamol-sources
python3 dynamol-sources/source-materials/materials/collector/collect_release_sources.py --verify-only --output dynamol-sources/source-materials
```

The copied verifier requires PyYAML. Verification checks the combined archive and material hashes; it does not reconstruct native binaries. The supplement is not standalone: the combined index also expects the reused archives from the base.

## What is recorded

- `INDEX.json` maps exact package name/version/build and native binary-package hashes to source archives, build recipes, patches, and license materials. Its source-collection ID and engine archive hashes identify the collected native payload; the final application binding is a separate acceptance record.
- `archives/` contains unmodified upstream archives and numbered AmberTools updates. Repeated source archives used by multiple binary packages are deduplicated by hash.
- `materials/native/` preserves shipped Conda recipes, resolved build configuration, patches, and notices. Use these with the corresponding source archive and the upstream build instructions. Original build tools, platform SDKs, and build dependencies are separate; no bit-identical rebuild was performed.
- `materials/python/` retains the shipped distribution metadata and hashes of compiled extensions. ParmEd's pinned source distribution and the recorded local wheel deployment-target override are included. A generic upstream source archive alone is not represented as a recreation of every wheel's original build environment.

The selected package set includes GROMACS, AmberTools, LGPL native libraries, the four direct copyleft Python distributions (ParmEd, MDTraj, Gemmi, certifi), and OpenMM's exact reported source commit because the wheel includes LGPL OpenCL platform plugins. OpenMM's CPU/core/application layers are MIT-licensed; the entire project is not labeled copyleft. The upstream [OpenMM license description](https://docs.openmm.org/latest/userguide/library/01_introduction.html#license) explains the split.

Packages with a recorded non-GPL alternative are listed separately rather than downloaded simply because their license expression contains `OR GPL`. GCC runtime sources used by the native environments are included conservatively, with their runtime exception distinguished. This is a materials collection, not a conclusion that each included package legally requires a source archive.

## AmberTools 24.8

The exact shipped recipe starts from `AmberTools24_rc5.tar.bz2`, applies the listed Conda source patches, and runs `update_amber --update-to=AmberTools.8` before compilation. All eight official numbered update files are included; their URLs and frozen hashes are in `materials/collector/ambertools-24.8-updates.json` and the index. Keep their upstream names `update.1` through `update.8` when using the upstream updater's local patch repository. The base updater identifies that repository as `.patches/AmberTools24_Unapplied_Patches` and the official source series as `https://ambermd.org/bugfixes/AmberTools/24.0/`.

Official `update.5` contains an `update.1` heading. Its original bytes are preserved; the download URL/filename and checksum identify it. No heading was rewritten. The full numbered patch application/rebuild has not been executed by this collector.

## Limits recorded explicitly

SciPy is permissively licensed but this wheel includes three compiler-runtime libraries. Their exact binary hashes, shipped license/exception text, and build configuration are preserved. GCC 13.4.0 source is included as a conservative supplement because it matches the reported compiler version; its official SHA-512 checksum is preserved and checked. The compiler version alone does not establish each runtime library's source revision. That uncertain provenance is recorded as a conditional component; the collector does not label it a universal release blocker or assume exception applicability.

This frozen collection is identified as `dynamol-three-runtime-sources-20260912`. Its input manifest records the exact staged Python and AmberTools, GROMACS, and ProMod3 runtime archives. The new loop runtime includes OpenStructure 2.11.1 and private OpenMM 8.5.1. Its recorded full OpenMM source commit is taken from the actual runtime archive and checked against the shipped recipe. A final application build has not yet been bound to this collection; engine archive and relevant native-file hashes must match before that binding is accepted.

Source/version and checksum correspondence is checked. Compiler signatures, complete reproducibility, every transitive license interpretation, and legal certification are outside this collection's claims. Availability of an archive is not by itself a legal determination. A release should keep this index, its checksums, and its source materials available together so that the mapping remains usable.
