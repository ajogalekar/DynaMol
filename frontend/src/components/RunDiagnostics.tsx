import { useEffect, useState } from 'react';
import { Activity, Download } from 'lucide-react';
import type { Job } from '../types';
import { scientificRequest, ScientificChart, saveCsv } from './ScientificCharts';
interface Point {
  step: number;
  time_ps: number;
  potential_kj_mol: number | null;
  kinetic_kj_mol: number | null;
  temperature_k: number | null;
}
interface Diagnostics {
  job_id: string;
  points: Point[];
  latest: Point | null;
  available: string[];
  recorded_points: number;
  invalid_rows: number;
  note: string;
  elapsed_seconds: number;
  estimated_remaining_seconds: number | null;
}
const duration = (s: number) =>
  s < 60
    ? `${Math.round(s)} s`
    : s < 3600
      ? `${(s / 60).toFixed(1)} min`
      : `${(s / 3600).toFixed(1)} h`;
export default function RunDiagnostics({ job }: { job: Job }) {
  const [result, setResult] = useState<Diagnostics | null>(null),
    [error, setError] = useState('');
  const [metric, setMetric] = useState<'potential_kj_mol' | 'kinetic_kj_mol' | 'temperature_k'>(
    'potential_kj_mol',
  );
  useEffect(() => {
    const controller = new AbortController();
    let timer: ReturnType<typeof setTimeout>;
    setResult(null);
    setError('');
    async function poll() {
      try {
        const data = await scientificRequest<Diagnostics>(
          `/api/jobs/${job.id}/diagnostics`,
          undefined,
          controller.signal,
        );
        if (!controller.signal.aborted) {
          setResult(data);
          setError('');
        }
      } catch (e) {
        if (!controller.signal.aborted) setError((e as Error).message);
      }
      if (!controller.signal.aborted && ['queued', 'running', 'cancelling'].includes(job.status))
        timer = setTimeout(poll, 2500);
    }
    void poll();
    return () => {
      controller.abort();
      clearTimeout(timer);
    };
  }, [job.id, job.status]);
  if (!['openmm', 'gromacs'].includes(job.engine)) return null;
  const current = result?.job_id === job.id ? result : null;
  const title =
    metric === 'temperature_k'
      ? 'Temperature'
      : metric === 'potential_kj_mol'
        ? 'Potential energy'
        : 'Kinetic energy';
  return (
    <details className="run-diagnostics" open>
      <summary>
        <Activity size={14} /> Simulation diagnostics{' '}
        <span>{current?.recorded_points ?? 0} observations</span>
      </summary>
      <div className="diagnostic-body">
        <div className="science-control-row">
          <label>
            Quantity
            <select
              aria-label="Diagnostic quantity"
              value={metric}
              onChange={(e) => setMetric(e.target.value as typeof metric)}
            >
              <option value="potential_kj_mol">Potential energy</option>
              <option value="kinetic_kj_mol">Kinetic energy</option>
              <option value="temperature_k">Temperature</option>
            </select>
          </label>
          <button
            type="button"
            className="text-button"
            disabled={!current?.points.length}
            onClick={() =>
              current &&
              saveCsv(`DynaMol-${job.id}-diagnostics.csv`, [
                ['step', 'time_ps', 'potential_kj_mol', 'kinetic_kj_mol', 'temperature_k'],
                ...current.points.map((p) => [
                  p.step,
                  p.time_ps,
                  p.potential_kj_mol,
                  p.kinetic_kj_mol,
                  p.temperature_k,
                ]),
              ])
            }
          >
            <Download size={13} /> CSV
          </button>
        </div>
        {error ? (
          <p className="inline-warning">{error}</p>
        ) : !current?.points.length ? (
          <p className="science-empty">
            Energy and temperature appear when dynamics records its first observation.
          </p>
        ) : !current.available.includes(metric) ? (
          <p className="science-empty">
            This run did not record {title.toLowerCase()}. New runs include native energy and
            temperature observations.
          </p>
        ) : (
          <ScientificChart
            x={current.points.map((p) => p.time_ps)}
            values={current.points.map((p) => p[metric])}
            title={title}
            unit={metric === 'temperature_k' ? 'K' : 'kJ/mol'}
            xLabel="Time (ps)"
          />
        )}
        {current?.latest && (
          <div className="diagnostic-status">
            <span>Elapsed {duration(current.elapsed_seconds)}</span>
            {current.estimated_remaining_seconds !== null && (
              <span>
                ~{duration(current.estimated_remaining_seconds)} remaining · includes setup
              </span>
            )}
            <span>Latest step {current.latest.step.toLocaleString()}</span>
          </div>
        )}
        {!!current?.invalid_rows && (
          <p className="inline-warning">
            {current.invalid_rows} invalid observations could not be displayed. Inspect the run log.
          </p>
        )}
        <p className="form-note">
          {current?.note ?? 'Native energy and temperature observations.'}
          {current && current.points.length < current.recorded_points
            ? ' Plot/CSV retain extrema from a reduced set; full observations are in the run Files.'
            : ''}
        </p>
      </div>
    </details>
  );
}
