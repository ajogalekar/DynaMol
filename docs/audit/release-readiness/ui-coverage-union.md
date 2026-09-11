# UI coverage union

**82 distinct cases; latest outcomes: 82 passed.**

This is the latest-per-test union of separate runs, not one monolithic run of all tests on one source revision.
The preparation rejection case was renamed to its intended readiness replacement and is counted once; old failures remain in history.
Historical product hashes were not captured. Report hashes check integrity of the evidence bytes; current source hashes describe the final source snapshot only and are not assigned retroactively to earlier runs.

| Evidence run | Passed | Failed | Skipped |
|---|---:|---:|---:|
| [workspace-ui.json](workspace-ui.json) | 27 | 0 | 0 |
| [workspace-final-ui.json](workspace-final-ui.json) | 10 | 0 | 0 |
| [workspace-legacy-ui.json](workspace-legacy-ui.json) | 44 | 5 | 0 |
| [workspace-release-ui.json](workspace-release-ui.json) | 13 | 0 | 0 |
| [workspace-cancellation-ui.json](workspace-cancellation-ui.json) | 2 | 0 | 0 |
| [analysis-ui.json](analysis-ui.json) | 4 | 0 | 0 |
| [recovery-validation.json](recovery-validation.json) | 1 | 0 | 0 |

Each case and its earlier failures, exact commands/environments, evidence SHA-256 hashes and a final source snapshot are recorded in [ui-coverage-union.json](ui-coverage-union.json). The final source snapshot is not evidence that every older case ran on that revision.

## Final coverage by spec

| Spec | Unique cases | Latest failed/skipped |
|---|---:|---:|
| e2e.spec.ts | 7 | 0 |
| fullscreen-header.spec.ts | 5 | 0 |
| hydrogen-context.spec.ts | 1 | 0 |
| live-measurements.spec.ts | 3 | 0 |
| modified-residue-inspection.spec.ts | 1 | 0 |
| pinch-zoom.spec.ts | 4 | 0 |
| preparation-jobs.spec.ts | 1 | 0 |
| preparation-progress.spec.ts | 3 | 0 |
| preparation.e2e.spec.ts | 5 | 0 |
| recovery-controls.spec.ts | 1 | 0 |
| structural-analysis.spec.ts | 4 | 0 |
| ui-audit-controls.spec.ts | 11 | 0 |
| ui-audit-imports.spec.ts | 15 | 0 |
| ui-audit-jobs.spec.ts | 4 | 0 |
| ui-audit-recovery.spec.ts | 4 | 0 |
| ui-audit-sources.spec.ts | 8 | 0 |
| workspace-library.spec.ts | 5 | 0 |

## Resolved initial failures

- Missing-heavy-atom and unresolved-loop tests now verify the intended readiness blocker and recovery after correcting the setting.
- The native complex download fixture includes its completed preparation job and original parameters; the retry passed.
- Early Stop/Cancel worker-exit classification is verified by rerunning the original native assertions after the supervisor fix.

Tests used isolated ports 8767, 8768 and 8874. The initial recovery attempt with a shared Playwright artifact-directory collision is excluded; its successful isolated retry is retained. No cases are counted as passed merely because a native fixture was unavailable.
