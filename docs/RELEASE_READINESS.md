# DynaMol final-workbench additions

This pass implements the approved workspace, recovery, readiness, library and analysis additions. The app remains a bounded local alpha. Large-trajectory streaming was deliberately deferred; the existing file/atom/frame limits apply.

| Addition | Where to use it | Evidence |
| --- | --- | --- |
| Automatic workspace restoration and named projects | Saved locally indicator; Projects toolbar button | [Workspace contract](WORKSPACES.md), backend archive/identity tests and real-browser scene/camera round trips |
| Checkpoint recovery | Resume on a stopped/interrupted simulation card | [Recovery contract](CHECKPOINT_RECOVERY.md), actual OpenMM and GROMACS interruption/continuation checks |
| Readiness and resource estimates | Above preparation and simulation actions | Shared submission guards; browser blocked-to-ready transitions; no jobs created by readiness |
| Searchable managed molecule library | Dataset name or Projects → Molecule library | Rename/archive/trash/restore browser workflows; dependency and source preservation tests |
| Live native diagnostics | Simulation cards → Simulation diagnostics | Native temperature/energy while running; full recorded observations in Files |
| RMSD/RMSF and named selections | Scene → Named selections; Structural analysis below the plot | Analytic known-answer tests, independent numerical comparison and actual CSV/JSON/plot interactions |
| Broader chemistry checks | Preparation inspection and readiness details | [43 positive native fixtures and 6 boundary/pipeline cases](CHEMISTRY_SUPPORT.md) |

The new calculation and export behavior is specified in [ANALYSIS.md](ANALYSIS.md). Native chemistry tests verify finite execution and parameter continuity in selected fixtures. They do not establish force-field accuracy, ionization populations, metal coordination quality or convergence. Unregistered unnatural amino acids, specialized metal models and covalent ligands/cofactors still receive explicit blockers.

## Review findings addressed

This integration found and fixed NAD atom-map recovery, native ligand improper-torsion conversion, ion-only explicit-solvent requirements and misleading CYM/CYX state reporting. The numerical review caught a planar quaternion-fit failure, malformed diagnostic step handling and non-increasing time axes. Structural analysis uses a proper Kabsch rotation and explicit saved-frame fallback. The browser audit caught a compressed analysis pane and verified its corrected desktop/narrow layout. Exact scene restoration also required saving whether analysis was expanded, because that changes canvas dimensions.

Portable backups bind saved scene data to ordered atom identities and physical-coordinate hashes. Imports reject malformed/oversized arrays and unsafe parameter files before committing. Download responses preserve private snapshots across concurrent clients/resume and clean up temporary files even on rejected ranges or disconnected clients.

## Validation records

The current backend output records **251 passed** in [backend-tests.txt](audit/release-readiness/backend-tests.txt). The combined UI evidence records **82 distinct cases with passing latest results and zero skips**, as a union of isolated runs rather than a claim that every test ran in one command. See [UI coverage](audit/release-readiness/ui-coverage-union.md), [analysis UI](audit/release-readiness/analysis-ui.json), [recovery validation](audit/release-readiness/recovery-validation.json), [native recovery](audit/release-readiness/recovery-native.json), and [chemistry support](CHEMISTRY_SUPPORT.md).

The [aggregate review](audit/release-readiness/release-review.md) separates software acceptance from unsupported scientific claims. Its uncertainty gate remains blocking for accuracy/convergence promotion: deterministic unit checks and five-step molecular fixtures cannot supply scientifically meaningful uncertainty intervals. This is an automated, model-generated evidence review, not an independent human scientific certification.

The local Apple Silicon package includes its runtimes and these guides/evidence. Apple Developer signing/notarization, an independent clean-machine test and complete public redistribution materials remain separate release work. The package does not add chemistry support or production-simulation validation beyond the checks described here.
