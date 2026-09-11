# DynaMol API contract

All routes under `/api`; local FastAPI server on 8765. Vite frontend on 5173 proxies `/api`.

- GET `/health` → Health from frontend/src/types.ts. engines openmm/gromacs.
- GET `/datasets` → Dataset[] (metadata, may omit atoms/bonds for list; detail mandatory).
- GET `/datasets/demo` → complete Dataset; bundled real protein MD trajectory with clear provenance. If not yet generated, real PDB structure acceptable honestly labeled.
- GET `/datasets/{id}` → Dataset.
- GET `/datasets/{id}/topology` → canonical single-model PDB, same atom ordering as metadata/coordinates.
- GET `/datasets/{id}/coordinates` → raw little-endian float32 Å coordinates, layout [frame, atom, xyz]. Whole small trajectory for v0.1, enforce load caps, document.
- POST `/datasets/upload` multipart `topology` required, `trajectory` optional, `stride` optional (default 1), `frame_interval_ps` optional for missing timestamps → Dataset.
- POST `/datasets/{id}/measurements` JSON `{kind: distance|angle|dihedral|hbond, atoms: number[]}` → `{kind, atoms, unit, values, times_ps, occupancy?, warnings?, angle_values?}`; atom indices zero based. hbond selection D,H,A with D-A <=3.5 Å AND D-H-A >=150°, validate explicit H and connectivity, label geometric occupancy. PBC-aware physical-coordinate analysis (not interpolated display). Å and degrees. Reject degenerate geometry with clear error.
- GET `/jobs` → Job[].
- POST `/jobs` → Job; body SimulationConfig. Real subprocess-backed engines; reject unsupported chemistry/parameters, no silent atom removals. Defaults modest laptop smoke / exploration. GROMACS explicit only; OpenMM implicit or explicit. Fixed protein FF per engine documented.
- GET `/jobs/{id}` → Job.
- POST `/jobs/{id}/cancel` → Job.
- GET `/jobs/{id}/download` → output archive (configs, topology, trajectory, logs, engine versions/seed/provenance).

Persist jobs/datasets locally in `data/` (gitignored). Errors structured `{detail: string}`. No fake progress or invented simulation/analysis. Bind localhost only, no arbitrary shell/path endpoints.

Viewer component contract (frontend/src/components/MolecularViewer.tsx): forwardRef exposing `fit()`, `zoom(factor:number)`, `focus(indices:number[])`, `snapshot(): Promise<Blob | null>`. Props: `dataset: Dataset | null`, `coordinates: Float32Array | null`, `frame: number` (fractional), `visibility: Visibility`, `representation: Representation`, `colorScheme: 'chain'|'residue'|'element'`, `selectedAtoms: number[]`, `measurements: Measurement[]`, `picking: boolean`, `spin: boolean`, `onAtomPick(index:number)`, `onReady()`, `onError(message:string)`. Keep canvas initial/empty background #0b121b. CSS container fills parent. Agent owns component and molecular viewer helper files only.

### Passive run measurements

`SimulationConfig.measurements` optionally supplies up to 12 measurement definitions. `GET /api/jobs/{id}/measurements` returns the current atomic snapshot, including nullable geometry and verified output atom indices. `POST /api/datasets/{target}/remap-measurements` maps retained definitions across preparation/solvation, returning explicit errors for lost or ambiguous identities. See [the live tracking contract](LIVE_TRACKING.md).
