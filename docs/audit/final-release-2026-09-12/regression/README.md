# Final release regression and API audit

This audit uses the repository's Python environment and isolated data roots. It does not mutate the scientist's running workspace, publish a release, or test a rebuilt app bundle. Native simulation stability and Chrome interaction checks are recorded separately by the other audit workstreams.

## Results

| Check | Result | Evidence |
| --- | --- | --- |
| Initial ordinary backend regressions | 594 passed; 5 native cases deliberately deferred | `backend-initial.txt`, `backend-initial.xml` |
| Deferred native ligand / solvent regressions | 5 passed; 594 deselected | `native-deferred.txt`, `native-deferred.xml` |
| Final complete backend suite | **601 passed, 0 failed, 0 skipped** in 21.70 s | `backend-final.txt`, `backend-final.xml` |
| Fresh API contract audit | 53 passed, 0 failed | `api-final.json`, `api-final.txt` |
| Initial API audit receipt | 52 passed, 1 audit-harness assertion mistake | `api-initial.json`, `api-initial.txt` |

The initial API failure was in the audit script: it expected `POST /api/workspace` to return scene state. That endpoint intentionally returns a save receipt. The corrected script uses `GET /api/workspace` to verify the restored scene. No application change was needed. The original failure and traceback remain intact.

`runtime.json` records Python and dependency versions. `source-initial.json`, `source-final.json`, and the source hashes embedded in the API JSON bind these results to the tested source. The final complete suite includes the new element-color workspace regression and the worker application-version provenance regression. `verification.json` confirms unchanged backend source hashes through the completed final suite and records artifact hashes.

The nine pytest warnings comprise two dependency deprecation notices, six expected OpenMM duplicate-number warnings from retained insertion-code regression fixtures, and one intentional non-finite-coordinate negative-test warning. There were no test failures or skipped tests.

## API coverage

The script makes real requests through FastAPI's `TestClient`, including multipart file uploads, response downloads, independent numeric expectations, live network fetches from the app's supported providers, and portable-project import into a second isolated data root.

| Area | Verified |
| --- | --- |
| Simulation structure inputs | PDB, CIF, mmCIF, PDBx, MOL2, SDF, MOL, SMI, SMILES |
| Trajectory topology inputs | PDB, CIF, mmCIF, PDBx, GRO, MDTraj HDF5 (`.h5`, `.hdf5`) |
| Trajectory files | PDB, XTC, TRR, DCD, NetCDF (`.nc`, `.netcdf`), MDTraj HDF5 (`.h5`, `.hdf5`), MDCRD (`.mdcrd`, `.crd`), LAMMPS dump, XYZ |
| Import contracts | Finite coordinates, metadata/binary dimensions, retained original bytes and SHA-256, frame stride, explicit frame interval, stereo and formal charge from SMILES |
| Invalid inputs | Empty, malformed, unknown extensions, bad CIF/PDB/SDF/MOL2, multiple SMILES, disconnected or invalid SMILES, non-UTF8 input; readable 422 errors with no partial datasets |
| Web fetches | RCSB `1UBQ`, extended accession `pdb_00001ubq`, PubChem aspirin by CID `2244` and name |
| Measurements | Immediate preview and full distance, angle, signed dihedral; H-bond occupancy; periodic minimum-image distance |
| Structural analysis | Known 2 Å rigid translation, fitted RMSD of zero, independently known 1 Å RMSF |
| Workspace and export | Save/restore camera, frame, selections and plots; ZIP export/import into a fresh data root; unchanged physical coordinates; no imported executable jobs |
| Color defaults | Fresh or legacy scene with no specified scheme defaults to element coloring; explicit saved chain coloring remains intact; fresh scene clears old measurements and selection |
| Read stability | 200 consecutive library/workspace/dataset/job requests across 50 rounds, no failed responses or accidental jobs |

The existing backend suite additionally covers topology mismatch, undefined geometry, live-measurement atom remapping, job download snapshots, recovery/cancellation state, archive integrity and traversal, chemistry guards, ion parameters, modified residues, loop limits, and resource caps. These are contract and regression checks, not proof that arbitrary protein–ligand chemistry will parameterize successfully.

The numeric fixtures are deliberately generated analytic coordinates. Their values are known controls; they are not presented as molecular simulation results. This API pass does not establish physical model accuracy or browser rendering performance.

## Browser obligations

The separate Chrome pass must exercise actual representation rendering and atom-color changes, ligand/water/hydrogen visibility, real atom picking and immediate labels, trajectory movie controls, plot updates, client-side CSV/PNG exports, file chooser controls, web-fetch dialogs, resize/zoom/fullscreen, saved-scene restoration, and a fresh unrelated dataset's color reset. API success alone does not verify those interactions.

## Reproduction

Run from the repository root:

```bash
PYTHONPATH=. .venv/bin/python docs/audit/final-release-2026-09-12/regression/api_contract_audit.py --output api-repeat.json
```

The script exclusively uses `data/integration-checks/final-release-2026-09-12/regression/api`. It adds fresh audit datasets and saves synthetic import/export fixtures there. Each run records source hashes and live-provider provenance. The source dependencies are the already-installed `.venv` packages listed in `runtime.json`.

No new format was advertised or silently accepted. MOL2/SDF/SMILES are molecular structure inputs; the trajectory-upload endpoint intentionally requires one of its supported topology formats. HDF5 support refers to MDTraj's molecular trajectory schema, not arbitrary HDF5 files.
