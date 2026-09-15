# Local final-candidate package acceptance

**Passed.** Candidate build `589ff2dddea916f5` was built from the final application source, copied from its new DMG into a private Applications folder, and tested using the copied app's private Python and engine runtimes. No GitHub upload or change to the user's Applications folder occurred.

- Candidate: `build/releases/final-audit-2026-09-12/DynaMol-0.1.1-macos-arm64-candidate.dmg`
- Size: 1,692,901,364 bytes
- SHA-256: `29adffe128d249f7040abcf05c65a6206a70da2d59c9375393b6c0d63b196e7a`
- Summary: `acceptance-summary.json`

## Evidence

| Check | Result | Receipt |
| --- | --- | --- |
| Complete application seal | Verified | `app-signature.json` |
| DMG integrity, read-only mount and drag-style copy | All recorded checks passed | `dmg-validation.json` |
| Private runtime operations | OpenMM CPU integration, AmberTools/SQM GAFF2/AM1-BCC, native GROMACS operations and ProMod3 ABI/data lookup passed | `runtime-validation.json` |
| Application acceptance | Both native workers completed; startup progress, service reuse, Host/Origin checks, diagnostics, 101-frame RMSD, workspace/project restoration and restart on a new port passed | `app-validation.json` |
| Exact copied-app Chrome smoke | Six checks passed, including demo rendering, element coloring, installed engines and final setup controls | `chrome-smoke.json` and screenshots |
| Cleanup | Owned service stopped; API unreachable; private service record removed; no active jobs | `cleanup.json`, `owned-smoke-server.json` |

Both native worker records identify `DynaMol 0.1.1` and the worker code included in this package. Native worker checks are short installation tests; the separate source audit supplies longer numerical and chemistry evidence.

The initial stock app build was retained in its separate build receipt. The final candidate used `build_filtered_candidate.py` to omit this audit's project ZIP and package-output folder from the recursive documentation copy. This local build-harness filter does not change application source. All existing published release artifacts were preserved.

The cleanup endpoint correctly returned HTTP 204 and shut down. An initial harness assertion expected HTTP 200; `cleanup.json` records that assertion correction and the subsequent completed cleanup verification.

This is local development-Mac evidence, with private prefixes and native-library path checks. It does not establish Apple notarization, quarantined-download opening, independent clean-Mac or minimum-OS compatibility. The app remains ad-hoc signed. ProMod3's package self-test verifies imports, pinned ABI and bundled data lookup; it does not repeat full loop construction. Public release receipts and source-materials binding must identify this exact candidate before upload.

## Final source-material follow-up

After the application acceptance above, the exact candidate also passed all 22 source-material binding checks. See [source-materials-binding.json](source-materials-binding.json). The aggregate report records this as completed; public uploading/committing of the candidate and corresponding release records has not been performed.
