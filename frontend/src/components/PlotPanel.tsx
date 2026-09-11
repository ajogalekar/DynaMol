import { useMemo, useState } from 'react';
import { Activity, Download, Plus, X } from 'lucide-react';
import type { Dataset, Measurement } from '../types';

interface Props {
  dataset: Dataset | null;
  measurements: Measurement[];
  frame: number;
  activeId: string | null;
  onActive: (id: string) => void;
  onRemove: (id: string) => void;
  onSeek: (frame: number) => void;
  onAdd: () => void;
}

/** Use physical timestamps for spacing; angular summaries respect the torsion branch cut. */
export function createPlotData(current: Measurement, timeUnit = 'ps') {
  const finite = current.values.filter(Number.isFinite);
  if (!finite.length) return null;
  const lo = Math.min(...finite),
    hi = Math.max(...finite),
    margin = Math.max((hi - lo) * 0.18, current.unit === 'Å' ? 0.04 : 1);
  const min = lo - margin,
    max = hi + margin,
    w = 1000,
    h = 120,
    left = 56,
    right = 16,
    top = 12,
    bottom = 24;
  const timestampsUsable =
    current.times_ps.length === current.values.length &&
    current.times_ps.every(
      (time, i, all) => Number.isFinite(time) && (i === 0 || time >= all[i - 1]),
    ) &&
    (current.values.length === 1 || current.times_ps.at(-1)! > current.times_ps[0]);
  const axisValues = timestampsUsable ? current.times_ps : current.values.map((_, i) => i);
  const start = axisValues[0],
    end = axisValues.at(-1)!,
    span = end - start;
  const atAxis = (value: number) =>
    left + (span > 0 ? (value - start) / span : 0) * (w - left - right);
  const x = (i: number) => atAxis(axisValues[Math.max(0, Math.min(axisValues.length - 1, i))]);
  const y = (value: number) => top + ((max - value) / (max - min)) * (h - top - bottom);
  const segments: string[] = [];
  let segment = '';
  let previous: number | undefined;
  current.values.forEach((value, i) => {
    const crossesWrap =
      current.kind === 'dihedral' && previous !== undefined && Math.abs(value - previous) > 180;
    if (!Number.isFinite(value) || crossesWrap) {
      if (segment) segments.push(segment);
      segment = '';
    }
    if (Number.isFinite(value)) {
      segment += `${x(i)},${y(value)} `;
      previous = value;
    } else previous = undefined;
  });
  if (segment) segments.push(segment);
  let mean: number | null = finite.reduce((a, b) => a + b, 0) / finite.length;
  if (current.kind === 'dihedral') {
    const sine =
      finite.reduce((sum, value) => sum + Math.sin((value * Math.PI) / 180), 0) / finite.length;
    const cosine =
      finite.reduce((sum, value) => sum + Math.cos((value * Math.PI) / 180), 0) / finite.length;
    mean = Math.hypot(sine, cosine) < 1e-6 ? null : (Math.atan2(sine, cosine) * 180) / Math.PI;
  }
  const nearestIndex = (fraction: number) => {
    const target = start + Math.max(0, Math.min(1, fraction)) * span;
    let lower = 0,
      upper = axisValues.length - 1;
    while (lower < upper) {
      const middle = Math.floor((lower + upper) / 2);
      if (axisValues[middle] < target) lower = middle + 1;
      else upper = middle;
    }
    return lower > 0 &&
      Math.abs(axisValues[lower - 1] - target) <= Math.abs(axisValues[lower] - target)
      ? lower - 1
      : lower;
  };
  const axisCaption = !timestampsUsable
    ? 'Frame index (timestamps are not ordered)'
    : timeUnit === 'frame'
      ? 'Original frame index (time not supplied)'
      : 'Time (ps)';
  const ticks = (axisValues.length === 1 ? [0] : [0, 0.25, 0.5, 0.75, 1]).map((fraction) => ({
    value: start + fraction * span,
    x: atAxis(start + fraction * span),
    fraction,
  }));
  return {
    min,
    max,
    w,
    h,
    left,
    right,
    top,
    bottom,
    x,
    y,
    segments,
    mean,
    lo,
    hi,
    nearestIndex,
    axisCaption,
    ticks,
    timestampsUsable,
  };
}

export default function PlotPanel({
  dataset,
  measurements,
  frame,
  activeId,
  onActive,
  onRemove,
  onSeek,
  onAdd,
}: Props) {
  const [hover, setHover] = useState<number | null>(null);
  const current = measurements.find((m) => m.id === activeId) ?? measurements[0];
  const plot = useMemo(
    () => (current ? createPlotData(current, dataset?.time_unit) : null),
    [current, dataset?.time_unit],
  );
  const selectedIndex = current
    ? Math.min(current.values.length - 1, Math.max(0, hover ?? Math.round(frame)))
    : 0;
  function download() {
    if (!current) return;
    const rows = [
      `${dataset?.time_unit === 'frame' ? 'frame' : 'time_ps'},${current.kind}_${current.unit === 'Å' ? 'angstrom' : 'degrees'}${current.angle_values ? ',DHA_angle_degrees' : ''}`,
      ...current.values.map(
        (v, i) =>
          `${current.times_ps[i]},${v}${current.angle_values ? ',' + current.angle_values[i] : ''}`,
      ),
    ];
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([rows.join('\n')], { type: 'text/csv' }));
    a.download = `dynamol-${current.kind}-${current.atoms.join('-')}.csv`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }
  return (
    <section className="plot-panel" aria-label="Trajectory measurements">
      <div className="plot-heading">
        <div className="section-label">
          <Activity size={14} /> MOTION, MEASURED
        </div>
        <div className="plot-heading-right">
          <span className="muted desktop-label">Click the plot to jump to a frame</span>
          <button
            className="icon-button compact"
            title="Export measurement CSV"
            aria-label="Export measurement CSV"
            disabled={!current}
            onClick={download}
          >
            <Download size={15} />
          </button>
          <button className="text-button" onClick={onAdd}>
            <Plus size={14} /> Measure
          </button>
        </div>
      </div>
      {current && plot ? (
        <>
          <div className="measurement-tabs">
            {measurements.map((m) => (
              <div
                key={m.id}
                className={`measurement-tab ${m.id === current.id ? 'active' : ''}`}
                style={{ '--plot-color': m.color } as React.CSSProperties}
              >
                <button
                  onClick={() => {
                    setHover(null);
                    onActive(m.id);
                  }}
                >
                  <i />
                  {m.label}
                </button>
                <button aria-label={`Remove ${m.label}`} onClick={() => onRemove(m.id)}>
                  <X size={11} />
                </button>
              </div>
            ))}
          </div>
          <div className="plot-content">
            <div className="plot-readout">
              <strong style={{ color: current.color }}>
                {Number.isFinite(current.values[selectedIndex])
                  ? current.values[selectedIndex].toFixed(2)
                  : '—'}
                <small>{current.unit}</small>
              </strong>
              <span>{hover !== null ? 'At cursor' : 'Current frame'}</span>
              {current.occupancy !== undefined && (
                <span className="occupancy">
                  {(current.occupancy * 100).toFixed(1)}% geometric occupancy
                </span>
              )}
            </div>
            <div className="chart-wrap">
              <svg
                viewBox={`0 0 ${plot.w} ${plot.h}`}
                preserveAspectRatio="none"
                role="img"
                aria-label={`${current.label}, plotted across ${current.values.length} saved frames`}
                onMouseLeave={() => setHover(null)}
                onMouseMove={(e) => {
                  const r = e.currentTarget.getBoundingClientRect();
                  setHover(
                    plot.nearestIndex(
                      (((e.clientX - r.left) / r.width) * plot.w - plot.left) /
                        (plot.w - plot.left - plot.right),
                    ),
                  );
                }}
                onClick={(e) => {
                  const r = e.currentTarget.getBoundingClientRect();
                  onSeek(
                    plot.nearestIndex(
                      (((e.clientX - r.left) / r.width) * plot.w - plot.left) /
                        (plot.w - plot.left - plot.right),
                    ),
                  );
                }}
              >
                {[0, 0.5, 1].map((p) => {
                  const v = plot.min + (plot.max - plot.min) * p;
                  return (
                    <g key={p}>
                      <line
                        x1={plot.left}
                        x2={plot.w - plot.right}
                        y1={plot.y(v)}
                        y2={plot.y(v)}
                        stroke="#27333f"
                        strokeDasharray="3 5"
                      />
                      <text
                        x={plot.left - 12}
                        y={plot.y(v) + 3}
                        textAnchor="end"
                        className="axis-text"
                      >
                        {v.toFixed(current.unit === 'Å' ? 2 : 0)}
                      </text>
                    </g>
                  );
                })}
                {plot.segments.map((points, i) => (
                  <polyline
                    key={i}
                    points={points}
                    fill="none"
                    stroke={current.color}
                    strokeWidth="2"
                    vectorEffect="non-scaling-stroke"
                    strokeLinejoin="round"
                  />
                ))}
                <line
                  x1={plot.x(Math.round(frame))}
                  x2={plot.x(Math.round(frame))}
                  y1={plot.top - 2}
                  y2={plot.h - plot.bottom}
                  stroke="#d7e5ee"
                  strokeWidth="1"
                  strokeDasharray="3 3"
                />
                {Number.isFinite(current.values[selectedIndex]) && (
                  <circle
                    cx={plot.x(selectedIndex)}
                    cy={plot.y(current.values[selectedIndex])}
                    r="3.5"
                    fill={current.color}
                    stroke="#0b121b"
                    strokeWidth="2"
                  />
                )}
                {plot.ticks.map(({ value, x, fraction }, i) => (
                  <text
                    key={i}
                    x={x}
                    y={plot.h - 5}
                    textAnchor={fraction === 0 ? 'start' : fraction === 1 ? 'end' : 'middle'}
                    className="axis-text"
                  >
                    {value.toFixed(2)}
                  </text>
                ))}
              </svg>
              <span className="axis-caption">{plot.axisCaption}</span>
            </div>
          </div>
          <div className="plot-footer">
            <span>
              {current.kind === 'dihedral' ? 'Circular mean' : 'Mean'}{' '}
              <b>{plot.mean === null ? 'undefined' : `${plot.mean.toFixed(2)} ${current.unit}`}</b>
              <i />
              {current.kind === 'dihedral' ? 'Wrapped range' : 'Range'}{' '}
              <b>
                {plot.lo.toFixed(2)}–{plot.hi.toFixed(2)} {current.unit}
              </b>
            </span>
            <span>Saved-frame analysis{dataset?.has_unitcell ? ' · periodic boundaries' : ''}</span>
          </div>
          {!!current.warnings?.filter((w) => !w.startsWith('No periodic box')).length && (
            <div className="inline-warning">
              {current.warnings.filter((w) => !w.startsWith('No periodic box')).join(' ')}
            </div>
          )}
        </>
      ) : (
        <div className="plot-empty">
          <Activity size={25} />
          <div>
            <strong>Follow a molecular interaction.</strong>
            <p>
              Pick atoms in the scene to plot a distance, angle, dihedral, or hydrogen bond over
              time.
            </p>
          </div>
          <button className="secondary-button" onClick={onAdd}>
            Add a measurement <Plus size={14} />
          </button>
        </div>
      )}
    </section>
  );
}
