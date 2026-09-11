import { useEffect, useRef, useState } from 'react';
import { CheckCircle2, AlertCircle, LoaderCircle, RotateCcw } from 'lucide-react';
import { scientificRequest } from './ScientificCharts';
interface Readiness {
  ready: boolean;
  blockers: string[];
  warnings: string[];
  model: string;
  resources: null | {
    estimated_atoms: number;
    estimate_kind: string;
    saved_frames: number;
    coordinate_bytes: number;
    disk_estimate_bytes: number;
    free_disk_bytes: number;
    note: string;
  };
  chemical_state: {
    ph: number;
    modifications: { residue: string; formal_charge?: number }[];
    ligands: { component_id?: string; residue?: string; formal_charge?: number }[];
    ions: { element: string; formal_charge: number }[];
  };
}
const bytes = (n: number) =>
  n >= 1024 ** 3 ? `${(n / 1024 ** 3).toFixed(1)} GiB` : `${Math.ceil(n / 1024 ** 2)} MiB`;
export default function ReadinessPanel({
  mode,
  datasetId,
  settings,
  onReadyChange,
  suspended = false,
}: {
  mode: 'preparation' | 'simulation';
  datasetId: string;
  settings: object;
  onReadyChange?: (ready: boolean) => void;
  suspended?: boolean;
}) {
  const [result, setResult] = useState<{ key: string; data?: Readiness; error?: string }>({
      key: '',
    }),
    [revision, setRevision] = useState(0);
  const callback = useRef(onReadyChange);
  callback.current = onReadyChange;
  const key = JSON.stringify([datasetId, mode, settings]);
  useEffect(() => {
    callback.current?.(false);
    if (suspended) return;
    const controller = new AbortController();
    const [id, kind, values] = JSON.parse(key);
    setResult({ key });
    const timer = setTimeout(() => {
      scientificRequest<Readiness>(
        `/api/datasets/${id}/readiness`,
        { mode: kind, settings: values },
        controller.signal,
      ).then(
        (data) => {
          if (!controller.signal.aborted) {
            setResult({ key, data });
            callback.current?.(data.ready);
          }
        },
        (e) => {
          if (!controller.signal.aborted) setResult({ key, error: (e as Error).message });
        },
      );
    }, 400);
    return () => {
      controller.abort();
      clearTimeout(timer);
    };
  }, [key, revision, suspended]);
  const value = result.key === key ? result : { key };
  const data = value.data;
  return (
    <div
      className={`readiness-panel ${data ? (data.ready ? 'readiness-ready' : 'readiness-blocked') : ''}`}
      aria-label={`${mode === 'preparation' ? 'Preparation' : 'Simulation'} readiness`}
      aria-live="polite"
    >
      <div className="readiness-heading">
        {data ? (
          data.ready ? (
            <CheckCircle2 size={16} />
          ) : (
            <AlertCircle size={16} />
          )
        ) : (
          <LoaderCircle size={16} className="spin" />
        )}
        <b>
          {suspended
            ? 'Computation in progress'
            : value.error
              ? 'Readiness check unavailable'
              : data
                ? data.ready
                  ? `Ready to ${mode === 'preparation' ? 'prepare' : 'run'}`
                  : 'Before you continue'
                : 'Checking structure and settings…'}
        </b>
        <button
          type="button"
          className="icon-button"
          aria-label="Recheck readiness"
          disabled={suspended}
          onClick={() => setRevision((v) => v + 1)}
        >
          <RotateCcw size={12} />
        </button>
      </div>
      {value.error && <p className="inline-warning">{value.error}</p>}
      {!!data?.blockers.length && (
        <ul className="readiness-blockers">
          {data.blockers.map((b) => (
            <li key={b}>{b}</li>
          ))}
        </ul>
      )}
      {data?.resources && (
        <div className="readiness-estimates">
          <span>
            <strong>{data.resources.estimated_atoms.toLocaleString()}</strong>{' '}
            {data.resources.estimate_kind === 'saved topology' ? 'atoms' : 'estimated atoms'}
          </span>
          <span>
            <strong>{data.resources.saved_frames.toLocaleString()}</strong> saved frames
          </span>
          <span>
            <strong>{bytes(data.resources.disk_estimate_bytes)}</strong> estimated disk
          </span>
          <span>
            <strong>{bytes(data.resources.free_disk_bytes)}</strong> free
          </span>
        </div>
      )}
      {data && (
        <details>
          <summary>Model and preparation details</summary>
          <p>{data.model}</p>
          <p>Requested/prepared pH {data.chemical_state.ph}; fixed-state assumptions apply.</p>
          {data.chemical_state.modifications.length > 0 && (
            <p>
              Modified residues:{' '}
              {data.chemical_state.modifications
                .map(
                  (m) =>
                    `${m.residue}${m.formal_charge !== undefined ? ` (${m.formal_charge} e)` : ''}`,
                )
                .join(', ')}
            </p>
          )}
          {data.chemical_state.ligands.length > 0 && (
            <p>
              Ligands:{' '}
              {data.chemical_state.ligands
                .map(
                  (m) =>
                    `${m.component_id ?? m.residue}${m.formal_charge !== undefined ? ` (${m.formal_charge} e)` : ''}`,
                )
                .join(', ')}
            </p>
          )}
          {data.chemical_state.ions.length > 0 && (
            <p>
              Ions:{' '}
              {data.chemical_state.ions
                .map((m) => `${m.element} (${m.formal_charge} e)`)
                .join(', ')}
            </p>
          )}
          {data.resources && <p>{data.resources.note}</p>}
          {data.warnings.map((w) => (
            <p key={w}>{w}</p>
          ))}
        </details>
      )}
    </div>
  );
}
