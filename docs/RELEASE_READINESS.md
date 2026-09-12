# DynaMol 0.1 early-release scope

The intended first release serves individual scientists using a local Mac:
load supported molecules, review preparation, run a bounded simulation and
explore its trajectory. GitHub source plus an unsigned self-contained Apple
Silicon archive is the proposed distribution. Use the matching acceptance record beside the archive; prior checks apply only
to their recorded builds. This page does not declare a particular archive
publicly ready. See [launch, data and packaging details](PACKAGING.md).

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

Previous relocated package tests passed on the development Mac with fresh
private data; they do not establish an independent clean-Mac launch or validate
newly rebuilt bits. The candidate presently has an ad-hoc signature only.
Developer ID signing/notarization is an optional way to reduce download/install
friction, not a condition for placing an unsigned beta on GitHub.

Public binary redistribution materials for all bundled components remain
incomplete. That unresolved work is separate from signing and must remain
visible before binary publication. See [packaging status](PACKAGING.md) and
[third-party notices](../THIRD_PARTY.md). No binary has been uploaded by this step.
