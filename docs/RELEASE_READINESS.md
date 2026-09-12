# DynaMol 0.1 early-release scope

**The 0.1.0 installer was withdrawn for a packaging signature error.** The 0.1.1
replacement remains pending its downloaded-app opening check. The scope and
historical software evidence below do not announce a replacement download.

DynaMol 0.1 serves individual scientists using a local Mac: load supported
molecules, review preparation, run a bounded simulation and explore its
trajectory. The unsigned Apple Silicon/macOS 14+ DMG contains the app and both
MD engines. Open the DMG, drag DynaMol to Applications and open it. The GitHub
source and source-materials companion are separate downloads for inspection or
rebuilding. See [installation, data and packaging details](PACKAGING.md).

Use the checksum and acceptance records supplied with each release asset.
Documentation and container changes are linked to prior native tests by payload
identity; necessary signature repairs require code/data comparison with the
changed signature bytes recorded separately. These are not fresh simulation runs.

| Feature | Where to use it | Evidence |
| --- | --- | --- |
| Persistent scenes, named Projects and backups | Projects toolbar | [Workspace contract](WORKSPACES.md), real scene/camera restoration and archive checks |
| Checkpoint recovery | Resume on a stopped/interrupted dynamics card | [Recovery contract](CHECKPOINT_RECOVERY.md), native OpenMM/GROMACS interruption tests; current resource and geometry rejection guards |
| Preparation/readiness and native diagnostics | Simulate | Shared submission guards, background progress, recorded energies and temperatures |
| Managed molecule library | Dataset name or Projects | Rename, trash, restore and source/dependency preservation tests |
| Measurements, RMSD/RMSF and selections | Explore | [Analysis contract](ANALYSIS.md), numerical controls and actual browser plot/export checks |
| Supported complex/modified-residue preparation | Preparation inspection | [Chemistry coverage](CHEMISTRY_SUPPORT.md), identity/parameter continuity and explicit unsupported-case blockers |

## Evidence and its limits

The [five-system pressure test](audit/five-system-pressure-test/REPORT.md)
completed five 100 ps trajectories: kinase, protease, phosphatase, channel and
GPCR; three OpenMM and two GROMACS. All five actual trajectory browser audits
passed, including movie controls, water/zoom, saved and posthoc plots and native
file downloads. All saved frames passed the recorded residue stereochemistry
checks; the GPCR's four defined ligand carbon centers were checked separately.

These are **software and gross numerical checks**, not evidence of physical
convergence, biological stability, force-field accuracy, protonation populations,
ligand pose accuracy or metal coordination quality. The initially rebuilt 4BVN
GPCR loop failed native stereochemistry checks; it was rejected and replaced by
a complete observed construct. Channel and GPCR runs were expressly water-only,
without membranes. Unsupported covalent chemistry and specialized ion/residue
models still produce explicit blockers.

The [final scoped backend run](audit/five-system-pressure-test/backend-regressions-final.json)
passed **265 tests with nine recorded warnings**. Earlier [251-test workbench
coverage](audit/release-readiness/backend-tests.txt) and [82 distinct UI cases](audit/release-readiness/ui-coverage-union.md)
remain separate historical evidence; counts are not added as though they were
one suite. The [final review](audit/five-system-pressure-test/review-final-postflight.json)
preserves generic checklist flags for missing replicate/convergence and
uncertainty evidence. They restrict scientific claims, not the stated smoke-test
results. Native runs retain their launch-time worker hashes; subsequent recovery
and diagnostic-history changes were tested separately.

## First-release boundaries

The default atom limit is 100,000. The source service supports an explicit
200,000 opt-in; the packaged launcher currently retains the default and offers
no toggle. The 145,913-atom GPCR pressure test therefore exceeds the packaged
profile. Large-trajectory streaming and membrane construction are outside 0.1.

The exact-ZIP acceptance record covers both bundled native workers, private
runtime and AmberTools checks, restart/workspace persistence, seven Chrome
pinch/zoom cases, and native launch with automatic default-browser opening.
The stronger installation-isolation record denies reads from developer
installations and common package prefixes; both engines and AmberTools still
passed. Its external harness waited on saved status while workers ran because
macOS `sandbox-exec` prohibits the system's setuid `ps`, then checked completed
jobs through the API. The original failure and harness diff are retained; the
application and filesystem-denial rules were unchanged.

The replacement requires strict signature-integrity checks before packaging,
after archive extraction and through the mounted/copied DMG, plus its own
downloaded/quarantined opening check. Earlier runtime acceptance did not detect
the withdrawn installer's invalid outer signature. Its replacement records must
bind the corrected package to the tested scientific code/data. These local
checks do not establish an independent clean-Mac, minimum-OS or Gatekeeper
result. Native Quit-menu clicking and physical Safari trackpad testing remain
unverified. A completed short installation check does not validate every
supported molecule or simulation setup.

The replacement seals the completed app last, preserving valid inner signatures
and repairing them only when necessary; bytecode writes into the bundle are
disabled. Its ad-hoc signature establishes integrity, not Apple developer trust
([TN2206](https://developer.apple.com/library/archive/technotes/tn2206/_index.html)).
It has no Developer ID signature, notarization or paid Apple enrollment step.
Downloaded copies may require
[Apple's per-app opening procedure](https://support.apple.com/en-us/102445).

The [source-materials companion](../packaging/source-materials/README.md) records
selected upstream archives, numbered AmberTools updates and exact shipped
recipe/patch mappings. It preserves conditional compiler-runtime provenance
and exception questions rather than presenting the collection as legal
certification. See [third-party notices](../THIRD_PARTY.md) and the release's
matching materials index.
The 0.1.0 source collection is reused for 0.1.1's unchanged upstream components;
the release binding must distinguish signature repairs from code/data changes.
