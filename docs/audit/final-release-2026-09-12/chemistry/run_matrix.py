"""Fresh bounded native matrix; preserves the prior release evidence files."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
DATA = ROOT / 'data/integration-checks/final-release-2026-09-12/chemistry'
os.environ['DYNAMOL_DATA_DIR'] = str(DATA)
os.environ['DYNAMOL_CPU_THREADS'] = '2'
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT))
import check_chemistry_matrix as matrix
from backend import storage

matrix.OUT = OUT
matrix.WORK = DATA / 'chemistry-matrix'
matrix.REPORT = OUT / 'chemistry-matrix.json'

if __name__ == '__main__':
    matrix.main()
    report = json.loads(matrix.REPORT.read_text())
    try:
        zinc = matrix.run_ion('ZN', 'Zn', 2)
    except Exception as exc:
        import traceback
        zinc = {'id': 'ZN', 'passed': False, 'error': str(exc), 'traceback': traceback.format_exc()}
    report['rows'].append(zinc)
    report['passed'] = all(row['passed'] for row in report['rows'])
    report['finished_at'] = dt.datetime.now(dt.timezone.utc).isoformat()
    report['audit_driver_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    storage.atomic_json(matrix.REPORT, report)
    print(json.dumps({'passed': report['passed'], 'cases': len(report['rows'])}), flush=True)
    if not report['passed']:
        raise SystemExit(1)
    import check_chemistry_boundaries as boundaries
    boundaries.OUT = OUT / 'chemistry-boundaries.json'
    boundaries.main()
    if not json.loads(boundaries.OUT.read_text())['passed']:
        raise SystemExit(1)
