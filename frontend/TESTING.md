# UI regressions use an isolated workspace

Application startup now restores the user's previous scene. Legacy UI regressions explicitly seed the real bundled demo before each test through `tests/testWorkspace.ts`; they require `DYNAMOL_TEST_ISOLATED=1` and refuse the usual user-facing ports 8765, 5173 and 4173. Never point the regression suite at a user's data root.

Create a disposable data directory with `datasets/demo` and `datasets/ubiquitin-start` copied from the bundled examples, and an empty `jobs` directory. Start the real API against it on a dedicated port, for example:

```bash
DYNAMOL_DATA_DIR=/absolute/path/to/disposable-data DYNAMOL_ORIGIN=http://127.0.0.1:8917 .venv/bin/python -m uvicorn backend.main:app --host 127.0.0.1 --port 8917
```

Build the frontend, then run from `frontend/` in another terminal:

```bash
npm run build
DYNAMOL_TEST_ISOLATED=1 DYNAMOL_BASE_URL=http://127.0.0.1:8917 npx playwright test --output=/absolute/path/to/ui-test-artifacts --reporter=list
```

Use a different artifact directory for each concurrent test run. Playwright clears its output directory at startup, so sharing the default `test-results` directory between agents can erase another run's downloads or traces. The API serves the production build; no viewer or API mocks replace successful molecular workflows. Some failure tests deliberately intercept requests to exercise error paths.

The 1UA2 inspection regression fetches its public RCSB fixture if `DYNAMOL_MODIFIED_DATASET` is not supplied. Tests requiring a completed native complex explicitly skip unless their documented fixture environment variables are provided. Live engine tests create and stop real jobs in the isolated data root; verify those jobs have stopped before discarding it.

`workspace-library.spec.ts` tests actual scene/camera restoration, analysis settings, project save/update/open/ZIP export/import, library rename/archive/trash/restore, and visible failure handling. Its molecular uploads are disposable fixtures. Backend portability/security tests use pytest temporary directories independently of any live service.
