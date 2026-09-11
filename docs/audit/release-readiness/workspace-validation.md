# Workspace, library and UI regression evidence

Verified in the source application on 2026-09-11. No final application bundle was rebuilt. The browser used the real production frontend, molecular API, WebGL viewer and installed native engines against an isolated data root on port 8767. The user's port 8765 library was not reset or edited by these tests.

## Implemented

- Server-side autosave and restore of dataset identity, exact NGL camera orientation/zoom, frame, representations, visibility, atom selections, named selections, geometry plots and structural-analysis settings/layout. Playback resumes paused.
- Named project snapshots with explicit replacement/removal and portable ZIP export/import.
- Searchable molecular library with rename, archive/unarchive, guarded recoverable trash/restore, and storage figures.
- Bounded archive validation, checksums, atom/trajectory binding, ancestry and preparation parameter preservation, collision rejection and rollback.

See [WORKSPACES.md](../../WORKSPACES.md) for behavior and current limits.

## Backend

`backend/tests/test_workspaces.py`: **29 passed**, including cross-root portability, preservation of a real capped-TPO parameter bundle's serialized OpenMM System, corrupt-state recovery, immutable project snapshots, protected source/job references, traversal/symlink/case-collision rejection, malicious XML, inconsistent coordinates and oversized NPY headers. Repeated after integrating the temporary download response cleanup helper. Tests use pytest temporary directories.

## Real browser

- `workspace-ui.json`: 27 passed, no skips. Workspace, general viewer controls, fullscreen, immediate measurements and pinch zoom.
- `workspace-final-ui.json`: 10 passed, no skips. Expanded fullscreen/project-dialog checks and the added delayed-project-open race regression; overlaps the preceding run.
- `workspace-release-ui.json`: **13 passed, no skips**, on final frontend asset `index-YBadVVti.js`. Five workspace workflows, five fullscreen checks, two updated preparation readiness workflows, and a complete prepared 9AX6 protein/ligand PDB + parameter ZIP download.
- `workspace-restored.json` and `workspace-restored.png`: actual before/after camera, frame, plots, named-selection and analysis-setting evidence and screenshot from the final release check.
- `workspace-legacy-ui.json`: initial remaining 49-case regression run, **44 passed / 5 failed / 0 skipped**. Passing cases cover genuine PDB/CIF/GRO/H5/XTC/DCD/TRR/NC/XYZ/MDCRD/LAMMPSTRJ uploads, MOL/MOL2/SDF/SMILES imports, RCSB/PubChem fetches, protein/ligand/water/ion visibility, 1UA2 TPO-versus-ATP classification, background protein repair/solvation, real OpenMM simulation and real GROMACS completion/log/output ZIP.

Three initial legacy failures were resolved and passed in `workspace-release-ui.json`: two old assertions expected invalid preparation to remain clickable, whereas the new readiness UI correctly blocks missing heavy-atom repair and unresolved internal loops; the third fixture lacked its completed preparation-job directory in the isolated root, which has now been copied with its original parameter files. The two cancellation failures exposed a real early-worker-exit status bug. After the supervisor fix, both unchanged native cases passed in `workspace-cancellation-ui.json` (2 passed / 0 failed / 0 skipped, 12.7 seconds) on the final frontend. All five initial failures therefore have explicit passing retry evidence.

The final complex fixture was copied read-only from datasets `4cf5bcedf38847c1` and its ancestor `9fb26d98962c453b`, plus stopped preparation job `57083a2415754858`. The 1UA2 classification fixture was copied from `f491cb784d7f4fbf`. No native complex case was skipped. Simulation and preparation regressions created disposable outputs under `data/integration-checks/workspace-ui`.

Results from overlapping runs must not be added as distinct tests. There are 28 distinct focused viewer/workspace cases and 49 distinct remaining legacy cases. Root-owned structural analysis and recovery-specific tests are recorded separately.

The latest-per-case union, including root structural-analysis and recovery-control checks, is recorded in [ui-coverage-union.md](ui-coverage-union.md) and [ui-coverage-union.json](ui-coverage-union.json). It retains historical failures, exact commands and source/evidence hashes rather than presenting the separate runs as one monolithic test run.
