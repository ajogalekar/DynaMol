"""Recompute an explicitly non-monolithic latest-per-case UI evidence union."""
from pathlib import Path
import hashlib, json, re, shlex
from collections import Counter
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ['workspace-ui', 'workspace-final-ui', 'workspace-legacy-ui', 'workspace-release-ui', 'workspace-cancellation-ui', 'analysis-ui']
ASSETS = {'workspace-final-ui': 'index-CG7UGSbg.js', 'workspace-legacy-ui': 'index-CG7UGSbg.js', 'workspace-release-ui': 'index-YBadVVti.js', 'workspace-cancellation-ui': 'index-YBadVVti.js', 'analysis-ui': 'index-YBadVVti.js'}
OLD = 'a rejected preparation exposes the real error beside the header and restores the prep button'
NEW = 'preparation readiness blocks incomplete atoms before submission and clears when repair is restored'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def specs(suite):
    yield from suite.get('specs', [])
    for child in suite.get('suites', []):
        yield from specs(child)

runs, cases = [], {}
for name in RUNS:
    path = HERE / f'{name}.json'
    if not path.exists():
        continue
    report = json.loads(path.read_text())
    env = {'DYNAMOL_TEST_ISOLATED': '1', 'DYNAMOL_BASE_URL': f'http://127.0.0.1:{8768 if name == "analysis-ui" else 8767}'}
    if name in {'workspace-legacy-ui', 'workspace-release-ui'}:
        env.update(DYNAMOL_AUDIT_COMPLEX_DATASET='4cf5bcedf38847c1', DYNAMOL_MODIFIED_DATASET='f491cb784d7f4fbf')
    env['PLAYWRIGHT_JSON_OUTPUT_FILE'] = str(path)
    argv = report['config']['argv'][2:]
    runs.append({'report': path.name, 'sha256': sha(path), 'stats': report['stats'], 'environment': env, 'command': 'npx playwright ' + shlex.join(argv), 'cwd': 'frontend', 'reported_tested_bundle': ASSETS.get(name), 'historical_product_sha256': None, 'top_level_errors': report.get('errors', [])})
    for suite in report['suites']:
        for spec in specs(suite):
            title = NEW if spec['title'] == OLD else spec['title']
            for test in spec['tests']:
                key = f"{spec['file']} :: {title}"
                final = test['results'][-1]
                attempt = {'report': path.name, 'original_title': spec['title'], 'started_at': final['startTime'], 'status': final['status'], 'classification': test['status'], 'duration_ms': final['duration'], 'line': spec['line']}
                cases.setdefault(key, {'file': spec['file'], 'title': title, 'history': []})['history'].append(attempt)

recovery_path = HERE / 'recovery-validation.json'
if recovery_path.exists():
    recovery = json.loads(recovery_path.read_text())
    browser = recovery['browser']
    title = 'Stop saves a checkpoint, Resume continues it, and terminal jobs retain Files'
    cases[f'recovery-controls.spec.ts :: {title}'] = {'file': 'recovery-controls.spec.ts', 'title': title, 'history': [{'report': recovery_path.name, 'original_title': title, 'started_at': recovery['recorded_at'], 'status': 'passed' if browser['passed'] == 1 else 'failed', 'classification': 'expected' if browser['passed'] == 1 else 'unexpected', 'duration_ms': int(browser['seconds'] * 1000), 'evidence_kind': 'agent-recorded successful real-browser run; list reporter, .last-run.json and native job artifacts', 'job_id': browser['job_id']}]}
    runs.append({'report': recovery_path.name, 'sha256': sha(recovery_path), 'stats': {'expected': browser['passed'], 'unexpected': 0 if browser['passed'] == 1 else 1, 'skipped': 0, 'duration': int(browser['seconds'] * 1000)}, 'environment': {'DYNAMOL_TEST_ISOLATED': '1', 'DYNAMOL_BASE_URL': 'http://127.0.0.1:8874', 'DYNAMOL_RECOVERY_TEST_INPUT': 'recovery-ui-protein'}, 'command': 'npx playwright test recovery-controls.spec.ts --output=/Users/ashujo/Documents/Science/DynaMol/data/integration-checks/recovery-browser-artifacts --reporter=list', 'cwd': 'frontend', 'reported_tested_bundle': 'index-YBadVVti.js', 'historical_product_sha256': None, 'checks': browser['checks']})

for case in cases.values():
    case['history'].sort(key=lambda v: v['started_at'])
    case['latest'] = case['history'][-1]
ordered = sorted(cases.values(), key=lambda v: (v['file'], v['title']))
counts = Counter(c['latest']['status'] for c in ordered)
current_files = []
for folder, patterns in [(ROOT/'backend', ['*.py']), (ROOT/'frontend/src', ['*.ts','*.tsx','*.css']), (ROOT/'frontend/tests', ['*.ts']), (ROOT/'frontend/dist/assets', ['*.js','*.css'])]:
    for pattern in patterns:
        current_files.extend(p for p in folder.rglob(pattern) if '__pycache__' not in p.parts)
source_hashes = {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(set(current_files))}
notes = [
    'This is the latest-per-test union of separate runs, not one monolithic run of all tests on one source revision.',
    'The preparation rejection case was renamed to its intended readiness replacement and is counted once; old failures remain in history.',
    'Historical product hashes were not captured. Report hashes authenticate the evidence bytes; current source hashes describe the final source snapshot only and are not assigned retroactively to earlier runs.',
    'Final layout/workspace/readiness/native-complex and structural-analysis checks were rerun on index-YBadVVti.js; early cancellation cases were rerun after the supervisor fix.',
    'Native recovery continuity is a software check on short trajectories, not proof of equilibration, convergence or bitwise reproducibility.',
]
output = {'generated_at': datetime.now(timezone.utc).isoformat(), 'scope': 'latest-per-test UI coverage union', 'summary': {'unique_cases': len(ordered), 'statuses': dict(counts), 'all_latest_passed': counts.get('passed', 0) == len(ordered)}, 'notes': notes, 'runs': runs, 'cases': ordered, 'current_source_snapshot': {'sha256': hashlib.sha256(json.dumps(source_hashes, sort_keys=True).encode()).hexdigest(), 'files': source_hashes}}
(HERE/'ui-coverage-union.json').write_text(json.dumps(output, indent=2)+'\n')
lines = ['# UI coverage union', '', f"**{len(ordered)} distinct cases; latest outcomes: " + ', '.join(f'{n} {s}' for s,n in sorted(counts.items())) + '.**', '', *notes[:3], '', '| Evidence run | Passed | Failed | Skipped |', '|---|---:|---:|---:|']
for run in runs:
    st = run['stats']; lines.append(f"| [{run['report']}]({run['report']}) | {st.get('expected', 0)} | {st.get('unexpected', 0)} | {st.get('skipped', 0)} |")
lines += ['', 'Each case and its earlier failures, exact commands/environments, evidence SHA-256 hashes and a final source snapshot are recorded in [ui-coverage-union.json](ui-coverage-union.json). The final source snapshot is not evidence that every older case ran on that revision.', '', '## Final coverage by spec', '', '| Spec | Unique cases | Latest failed/skipped |', '|---|---:|---:|']
for filename in sorted({c['file'] for c in ordered}):
    matching = [c for c in ordered if c['file'] == filename]
    lines.append(f"| {filename} | {len(matching)} | {sum(c['latest']['status'] != 'passed' for c in matching)} |")
lines += ['', '## Resolved initial failures', '', '- Missing-heavy-atom and unresolved-loop tests now verify the intended readiness blocker and recovery after correcting the setting.', '- The native complex download fixture includes its completed preparation job and original parameters; the retry passed.', '- Early Stop/Cancel worker-exit classification is verified by rerunning the original native assertions after the supervisor fix.', '', 'Tests used isolated ports 8767, 8768 and 8874. The initial recovery attempt with a shared Playwright artifact-directory collision is excluded; its successful isolated retry is retained. No cases are counted as passed merely because a native fixture was unavailable.', '']
(HERE/'ui-coverage-union.md').write_text('\n'.join(lines))
print(json.dumps(output['summary']))
