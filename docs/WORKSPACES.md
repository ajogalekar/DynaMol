# Workspaces, projects and portable backups

DynaMol saves the current scene to its local molecular service. This is server-side storage, so a packaged application's changing localhost port does not lose its workspace. Reopening DynaMol restores the dataset, camera orientation/zoom, frame, representation and coloring, visibility controls, selected atoms, named selections, geometry measurement plots, structural-analysis settings and expanded/collapsed panel, and playback speed/loop preference. Playback reopens paused. Structural RMSD/RMSF results recalculate when you click Calculate; their saved settings identify the requested analysis.

The toolbar displays the save state. Changes are written about every 700 ms while the viewer is ready; browser-close delivery is also attempted for small scenes. Wait for **Saved locally** before closing if a just-made change matters. A failed save is visible in the toolbar and Projects dialog. A saved scene is bound to the dataset's ordered atom identities and physical-trajectory checksum; it is rejected if that topology or trajectory has changed rather than applying selections or plots to the wrong coordinates.

**Projects** are named snapshots. Save a new project, open a previous snapshot, explicitly replace one with the current scene, or remove the named snapshot while retaining its molecular files. Ordinary autosaving does not overwrite named snapshots. Pick atoms or use **Go to a residue**, give the selection a name, and save it from the Scene pane. Clicking a named selection selects and focuses those exact atoms; these selections are also available to trajectory analysis.

## Molecular library

Click the current dataset name, or open Projects → Molecule library. Search by name/source/identifier, rename an entry, and switch between active, archived, and trashed structures. Storage figures show molecular datasets, retained run files, recoverable trash, and available disk space.

Archive hides an entry from the active list without changing its files or simulation provenance. Moving a structure to trash requires its exact name and preserves all files for restoration. Bundled examples, the current workspace, named-project inputs, derived structures' source datasets, and retained jobs' input/output datasets are protected from trash. Use Archive for those entries. Trash retains disk space; permanent erasure is deliberately not exposed in this alpha.

## Project backups

Use the download button beside a named project. Its ZIP contains:

- The project scene, measurements and named selections.
- The complete active dataset and all source ancestors referenced by metadata/provenance.
- Original uploaded files, physical and display trajectories, topology, and all retained preparation/ligand/modified-residue parameter files.
- Related stopped runs' output files, logs and provenance as historical records.
- A versioned manifest with file sizes, SHA-256 checksums and ordered-atom fingerprints.

Import the ZIP through **Import project backup**. DynaMol validates paths, file counts/sizes, archive checksums, topology identities/order, frame dimensions, finite coordinates, periodic-cell dimensions, and preparation parameter manifests before committing files. It rejects symlinks, case-colliding paths, malformed/nested oversized coordinate archives, unsafe XML scripts/includes, external dataset URLs and incomplete ancestry. No code from a backup is executed. SHA-256 detects changed bytes; it does not authenticate the person who created a backup or establish scientific validity.

Dataset identifiers and source bytes remain intact. Identical existing datasets are reused. An existing identifier with different files is a conflict: import stops without overwriting it. Import into a separate/empty workspace in that case. Historical run files are retained under the imported project and survive re-export; importing never launches or adds a historical job to the live queue. Active related jobs must stop before export so the backup represents a consistent state.

This is a **project backup**, including required sources, rather than an export of every unrelated library entry. Save/export separate projects for unrelated studies. Limits are 100 related datasets, 10,000 files, 500 MiB compressed and 2 GiB expanded; the existing viewer atom/frame/coordinate limits also apply on import. Oversized exports/imports fail visibly.

## Local storage

The workspace service stores scenes and project records under `DYNAMOL_DATA_DIR/workspaces`; molecular datasets and jobs remain in their existing sibling directories. In the macOS package, the data root is `~/Library/Application Support/DynaMol/Workspace`. Development uses the repository's `data` directory unless `DYNAMOL_DATA_DIR` is set. These records belong to the user and are not bundled into release archives.

Validation is covered by isolated backend tests in `backend/tests/test_workspaces.py` and real-browser workflows in `frontend/tests/workspace-library.spec.ts`.
