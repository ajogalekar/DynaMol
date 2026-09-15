# Local quantum runtime recovery v1

The relocated runtime is ready for a bounded research recovery outside Documents. The original runtime was not edited.

- Python: `/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/bin/python`
- Resolved interpreter: `/Users/ashujo/.local/share/uv/python/cpython-3.12.13-macos-aarch64-none/bin/python3.12`
- Interpreter SHA-256: `f1c4ad59c4adf1e5dfb5e9f558f7b95e4d7a62387d4bb4572fcf84f42648fee3`
- Runtime manifest: `/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/dynamol-qm-runtime.json`
  SHA-256 `91a060defde8a3506861dffbb5277a4672bbfd041208c2cbbfeedb66dfe0c1b9`
- Complete non-generated file inventory: `/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/local-recovery-evidence/runtime-files.json`
  SHA-256 `683b958788605f35ce2fc9a5b7ac94edb03f5316bd0b166cfcb415ffaea1f416`
- Detailed report: `/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/local-recovery-evidence/recovery-report.json`

## Copy and relocation

All 12,444 original regular files and six symbolic links were copied. Source/copy SHA-256 values are retained in `copy-manifest.json`, with explicit revalidation of two scripts in `copy-revalidation.json`. The initial two script checks rejected changing source metadata during File Provider hydration; revalidation found stable, identical bytes and did not rewrite their content. The initial failure record is retained. The source inventory contained 3,042 initially dataless files; targeted Foundation downloads were requested before reads.

The original native OpenMP library was copied byte-for-byte from `.tools/resp-build/lib/libomp.dylib` into the local PySCF library directory (SHA-256 `8cd3b360f4fa84036a617cbfb7fc3ef22e659f3f7e6674930b467b4f8a6cd5e4`). Sixteen copied PySCF dylibs had two obsolete Documents RPATHs removed and received local ad-hoc signatures. **Every file-backed Mach-O section payload is identical before and after relocation.** Whole-file hashes change because loader metadata and signatures change. Original copied binaries, commands, signing checks, section hashes and pre/post binary hashes are retained in `native-relocation.json` and `pre-relocation-native/`.

Fifteen copied activation/CLI scripts had their runtime prefix changed to the local cache. Interpreter symlink targets and `pyvenv.cfg` remain unchanged; the latter retains a historical creation command, while the actual interpreter and home are outside Documents. No sitecustomize/usercustomize file was found, the sole `.pth` file has no Documents path, and no copied executable script still points to the original prefix.

The final inventory pins **9,051 regular files (332,079,540 bytes), six links, and all 177 native libraries**. Only audit evidence and generated Python bytecode caches are excluded. No final pinned file is dataless. `native-linkage-after.json` checks all 177 native libraries: no dependency or RPATH enters Documents. The post-calculation loaded module paths and native image paths likewise remain outside Documents.

## Memory accounting

The unmodified macOS PySCF `current_memory()` implementation falls back to `(0, 0)` when psutil cannot be imported. **psutil 7.2.2 was added only to the new runtime**, from the official PyPI arm64 wheel, SHA-256 `1a7b04c10f32cc88ab39cbf606e117fd74721c831c98a27dc04578deb0c16979`; package metadata, exact wheel and installation record are retained. After the native test, PySCF reports 126.451712 MB, matching both psutil and `ps` at 126,451,712 RSS bytes. This restores allocation accounting; it is not an operating-system memory limit. The parent supervisor still needs its external sampled RSS bound.

## Numerical relocation check

Exactly one new, three-atom water RHF/STO-3G single-point energy/gradient calculation was run. Its input coordinates and resolved basis/ECP/method are identical to the retained qualified original-runtime reference. It held the shared `QM-LAUNCH.lock`, checked fewer than three active workers and 25 GiB free disk, used two threads and a 2,000 MB allocation hint, with a 60-second native and 70-second process-group outer timeout. Inputs, scratch, worker copy and outputs all live in the local cache.

- Energy: `-74.96361400506848` Ha; reference `-74.96361400506842` Ha.
- Absolute energy difference: **5.68434188608e-14 Ha**.
- Maximum absolute gradient difference: **6.85229650799e-12 Ha/Bohr**.
- Both fixed comparison tolerances: `1e-10`; no tolerance was relaxed.
- SCF converged in eight cycles; native worker elapsed 0.268178 seconds.
- Actual native OpenMP and PySCF thread counts: two each.
- Frozen worker SHA-256: `ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1` (byte-identical to the original worker).

The first post-run method comparison rejected two newly recorded observation fields (`requested_threads` and `actual_pyscf_threads`) absent from the older reference. That failure is preserved in `initial-comparison-failure.json`. The recheck verifies all original method fields exactly and requires both additional thread fields to equal two. No numerical calculation was repeated.

## Scope

This qualifies runtime relocation and memory reporting, not a large adduct, geometry optimization, force field or release bundle. It preserves the existing external uv-managed interpreter and stdlib. It does not repair previously offloaded checkpoint bytes. No source runtime, central status/checkpoint or active worker was edited. The cache is outside Documents, and the parent must place all new checkpoint, output and scratch files there as well.
