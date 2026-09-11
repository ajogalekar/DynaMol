"""Bounded, read-only native energy/temperature diagnostics for running jobs."""
import csv
import math
from fastapi import APIRouter, Query
from . import config, jobs, storage

router = APIRouter(prefix='/api')
FIELDS = ('potential_kj_mol', 'kinetic_kj_mol', 'temperature_k')


def read_diagnostics(job_id, max_points=600):
    job = jobs.get_job(job_id)
    path = config.JOBS_DIR / storage.safe_id(job_id) / 'energies.csv'
    rows, invalid = [], 0
    if path.is_file():
        if path.stat().st_size > 16 * 1024**2:
            raise ValueError('Diagnostic file exceeds the bounded reader limit.')
        # Ignore a final partially written line until the next poll.
        text = path.read_text(errors='replace')
        text = text[:text.rfind('\n') + 1]
        for raw in csv.DictReader(text.splitlines()):
            try:
                step = float(raw['step'])
                if not math.isfinite(step) or step < 0 or not step.is_integer():
                    raise ValueError('step')
                row = {'step': int(step), 'time_ps': float(raw['time_ps'])}
                if not math.isfinite(row['time_ps']) or row['time_ps'] < 0:
                    raise ValueError('time')
                for field in FIELDS:
                    value = raw.get(field)
                    row[field] = None if value in (None, '') else float(value)
                    if row[field] is not None and not math.isfinite(row[field]):
                        raise ValueError(field)
                if rows and (row['step'] < rows[-1]['step'] or row['time_ps'] < rows[-1]['time_ps']):
                    raise ValueError('Nonmonotonic diagnostic step')
                if rows and row['step'] == rows[-1]['step']:
                    rows[-1] = row
                else:
                    rows.append(row)
            except (ValueError, TypeError, KeyError, OverflowError):
                invalid += 1
    indices = set(range(len(rows)))
    if len(rows) > max_points:
        indices = {0, len(rows) - 1}
        size = math.ceil(len(rows) / max(1, (max_points - 2) // 6))
        for start in range(0, len(rows), size):
            bucket = range(start, min(start + size, len(rows)))
            for key in FIELDS:
                valid = [i for i in bucket if rows[i][key] is not None]
                if valid:
                    indices.add(min(valid, key=lambda i: rows[i][key]))
                    indices.add(max(valid, key=lambda i: rows[i][key]))
    rate, remaining = None, None
    if len(rows) >= 2:
        elapsed = job.get('elapsed_seconds', 0)
        completed = job.get('completed_steps', 0)
        # Includes setup time; explicitly a measured overall estimate, not a benchmark.
        if elapsed > 0 and completed > 0:
            rate = completed / elapsed
            if job['status'] in {'running', 'queued', 'cancelling'}:
                remaining = max(0, (job['total_steps'] - completed) / rate)
    return {'job_id': job_id, 'status': job['status'], 'stage': job['stage'],
            'points': [rows[i] for i in sorted(indices)], 'recorded_points': len(rows),
            'latest': rows[-1] if rows else None, 'invalid_rows': invalid,
            'available': [key for key in FIELDS if any(row[key] is not None for row in rows)],
            'target_temperature_k': job.get('config', {}).get('temperature_k'),
            'elapsed_seconds': job.get('elapsed_seconds', 0), 'estimated_remaining_seconds': remaining,
            'steps_per_second': rate, 'estimate_note': 'Overall measured throughput includes setup; remaining time is approximate.',
            'note': 'Native saved energy and kinetic-temperature observations. A stable-looking trace does not establish equilibration or convergence.'}


@router.get('/jobs/{job_id}/diagnostics')
def diagnostics(job_id: str, max_points: int = Query(default=600, ge=20, le=2000)):
    return read_diagnostics(job_id, max_points)
