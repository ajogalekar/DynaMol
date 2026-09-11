import { useEffect, useState } from 'react';
import { LoaderCircle } from 'lucide-react';
import { api } from '../api';
import type { Dataset, MeasureKind } from '../types';
import './LiveMeasurement.css';

export interface MeasurementPreview {
  kind: MeasureKind;
  atoms: number[];
  unit: string;
  values: (number | null)[];
  times_ps: number[];
  frame_errors: (string | null)[];
  angle_values?: (number | null)[];
  geometry_passes?: (boolean | null)[];
  warnings: string[];
}

export interface LiveMeasurementState {
  status: 'idle' | 'loading' | 'ready' | 'error';
  kind: MeasureKind;
  atoms: number[];
  frame: number;
  value?: number;
  unit?: string;
  angle?: number;
  geometryPasses?: boolean;
  error?: string;
  retry: () => void;
}

/** Fetch physical geometry once per selection; playback never sends frame-by-frame requests. */
export function useLiveMeasurement(
  dataset: Dataset | null,
  kind: MeasureKind,
  atoms: number[],
  frame: number,
): LiveMeasurementState {
  const count = { distance: 2, angle: 3, dihedral: 4, hbond: 3 }[kind];
  const complete = !!dataset && atoms.length === count;
  const selectionKey = complete ? JSON.stringify([dataset.id, kind, atoms]) : '';
  const [revision, setRevision] = useState(0);
  const [response, setResponse] = useState<{
    key: string;
    preview?: MeasurementPreview;
    error?: string;
  }>({ key: '' });
  useEffect(() => {
    if (!selectionKey) {
      setResponse({ key: '' });
      return;
    }
    const [datasetId, selectedKind, selectedAtoms] = JSON.parse(selectionKey) as [
      string,
      MeasureKind,
      number[],
    ];
    const controller = new AbortController();
    setResponse({ key: selectionKey });
    api.previewMeasurement(datasetId, selectedKind, selectedAtoms, controller.signal).then(
      (preview) => {
        if (!controller.signal.aborted) setResponse({ key: selectionKey, preview });
      },
      (error: Error) => {
        if (!controller.signal.aborted) setResponse({ key: selectionKey, error: error.message });
      },
    );
    return () => controller.abort();
  }, [selectionKey, revision]);
  const savedFrame = Math.max(0, Math.min((dataset?.n_frames ?? 1) - 1, Math.round(frame)));
  const base = {
    kind,
    atoms,
    frame: savedFrame,
    retry: () => setRevision((value) => value + 1),
  };
  if (!complete) return { ...base, status: 'idle' };
  // A changed selection must never briefly inherit the preceding selection's numeric readout.
  if (response.key !== selectionKey) return { ...base, status: 'loading' };
  if (response.error) return { ...base, status: 'error', error: response.error };
  if (!response.preview) return { ...base, status: 'loading' };
  const preview = response.preview;
  const value = preview.values[savedFrame];
  const error = preview.frame_errors[savedFrame];
  if (error || value === null || value === undefined || !Number.isFinite(value))
    return { ...base, status: 'error', error: error || 'Geometry is unavailable for this frame.' };
  return {
    ...base,
    status: 'ready',
    value,
    unit: preview.unit,
    angle: preview.angle_values?.[savedFrame] ?? undefined,
    geometryPasses: preview.geometry_passes?.[savedFrame] ?? undefined,
  };
}

export function liveMeasurementLabel(preview: LiveMeasurementState): string {
  if (preview.status !== 'ready' || preview.value === undefined) return '';
  const value = `${preview.value.toFixed(preview.unit === '°' ? 1 : 2)} ${preview.unit}`;
  return preview.kind === 'hbond' && preview.angle !== undefined
    ? `D–A ${value} · D–H–A ${preview.angle.toFixed(1)}°`
    : value;
}

export default function LiveMeasurement({ preview }: { preview: LiveMeasurementState }) {
  if (preview.status === 'idle') return null;
  return (
    <div
      className={`live-measurement live-measurement--${preview.status}`}
      aria-label="Current measurement"
    >
      <div className="live-measurement__heading">
        <span>{preview.kind === 'hbond' ? 'Hydrogen-bond geometry' : 'Current measurement'}</span>
        <small>Frame {preview.frame + 1}</small>
      </div>
      {preview.status === 'loading' ? (
        <p>
          <LoaderCircle size={14} className="spin" /> Measuring selected atoms…
        </p>
      ) : preview.status === 'error' ? (
        <>
          <p className="live-measurement__error">{preview.error}</p>
          <button className="text-button" onClick={preview.retry}>
            Retry measurement
          </button>
        </>
      ) : (
        <>
          <output className="live-measurement__value" aria-live="off">
            {liveMeasurementLabel(preview)}
          </output>
          {preview.kind === 'hbond' && (
            <p className="live-measurement__criterion">
              {preview.geometryPasses ? 'Meets' : 'Outside'} geometric cutoffs: ≤3.5 Å, ≥150°.{' '}
              Chemical donor/acceptor eligibility still needs review.
            </p>
          )}
          <small className="live-measurement__note">
            Saved-frame geometry · updates as you move through the trajectory.
          </small>
        </>
      )}
    </div>
  );
}
