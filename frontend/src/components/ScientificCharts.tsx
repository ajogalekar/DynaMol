import { useMemo, useState } from 'react';
import './ScientificAnalysis.css';
export function saveJson(name: string, value: unknown) {
  const url = URL.createObjectURL(
    new Blob([JSON.stringify(value, null, 2)], { type: 'application/json' }),
  );
  const a = document.createElement('a');
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
export function saveCsv(name: string, rows: (string | number | null)[][]) {
  const csv = rows
    .map((row) => row.map((v) => `"${String(v ?? '').replaceAll('"', '""')}"`).join(','))
    .join('\n');
  const url = URL.createObjectURL(new Blob([csv], { type: 'text/csv;charset=utf-8' }));
  const a = document.createElement('a');
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
export async function scientificRequest<T>(
  path: string,
  body?: unknown,
  signal?: AbortSignal,
): Promise<T> {
  const response = await fetch(
    path,
    body === undefined
      ? { signal }
      : {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(body),
          signal,
        },
  );
  const json = await response.json();
  if (!response.ok)
    throw new Error(
      typeof json.detail === 'string'
        ? json.detail
        : 'The request could not be completed. Check the selected settings.',
    );
  return json;
}
export function ScientificChart({
  x,
  values,
  title,
  unit,
  xLabel,
  labels,
  currentIndex,
  onPick,
}: {
  x: number[];
  values: (number | null)[];
  title: string;
  unit: string;
  xLabel: string;
  labels?: string[];
  currentIndex?: number;
  onPick?: (index: number) => void;
}) {
  const [hover, setHover] = useState<number | null>(null);
  const plot = useMemo(() => {
    const valid = values.filter((v): v is number => v !== null && Number.isFinite(v));
    if (!valid.length || !x.length) return null;
    const lo = valid.reduce((a, b) => Math.min(a, b)),
      hi = valid.reduce((a, b) => Math.max(a, b)),
      pad = Math.max((hi - lo) * 0.1, Math.abs(hi) * 0.0001, 0.01);
    const min = lo - pad,
      max = hi + pad,
      start = x[0],
      end = x.at(-1)!;
    const px = (i: number) => 65 + ((x[i] - start) / (end - start || 1)) * 715;
    const py = (v: number) => 15 + ((max - v) / (max - min)) * 130;
    let d = '',
      drawing = false;
    values.forEach((v, i) => {
      if (v === null || !Number.isFinite(v)) {
        drawing = false;
        return;
      }
      d += `${drawing ? ' L' : ' M'}${px(i)},${py(v)}`;
      drawing = true;
    });
    return { lo, hi, min, max, start, end, px, py, d };
  }, [x, values]);
  if (!plot) return <p className="science-empty">Waiting for recorded values…</p>;
  const index = hover !== null && hover < values.length ? hover : currentIndex;
  const nearest = (e: React.MouseEvent<SVGSVGElement>) => {
    const rect = e.currentTarget.getBoundingClientRect(),
      point = (((e.clientX - rect.left) / rect.width) * 800 - 65) / 715;
    const target = plot.start + Math.max(0, Math.min(1, point)) * (plot.end - plot.start);
    let result = 0;
    x.forEach((value, i) => {
      if (Math.abs(value - target) < Math.abs(x[result] - target)) result = i;
    });
    return result;
  };
  return (
    <div className="scientific-chart">
      <div className="scientific-chart-heading">
        <b>{title}</b>
        <output>
          {index !== undefined &&
          index !== null &&
          values[index] !== undefined &&
          values[index] !== null
            ? `${labels?.[index] ?? x[index].toFixed(3)} · ${values[index]!.toFixed(3)} ${unit}`
            : `${plot.lo.toFixed(2)}–${plot.hi.toFixed(2)} ${unit}`}
        </output>
      </div>
      <svg
        viewBox="0 0 800 180"
        role="img"
        aria-label={`${title}, ${xLabel}, ${unit}`}
        onMouseMove={(e) => setHover(nearest(e))}
        onMouseLeave={() => setHover(null)}
        onClick={(e) => onPick?.(nearest(e))}
        className={onPick ? 'is-interactive' : ''}
      >
        {[0, 0.5, 1].map((f) => (
          <g key={f}>
            <line x1="65" x2="780" y1={15 + f * 130} y2={15 + f * 130} />
            <text x="59" y={19 + f * 130} textAnchor="end">
              {(plot.max - f * (plot.max - plot.min)).toLocaleString(undefined, {
                maximumFractionDigits: 2,
              })}
            </text>
          </g>
        ))}
        <path className="scientific-chart-line" d={plot.d} />
        {values.length === 1 && values[0] !== null && (
          <circle cx={plot.px(0)} cy={plot.py(values[0])} r="3" className="scientific-chart-dot" />
        )}
        {index !== undefined && index !== null && index < x.length && (
          <line
            className="scientific-chart-cursor"
            x1={plot.px(index)}
            x2={plot.px(index)}
            y1="15"
            y2="145"
          />
        )}
        {[0, 0.5, 1].map((f) => (
          <text
            key={f}
            x={65 + f * 715}
            y="162"
            textAnchor={f === 0 ? 'start' : f === 1 ? 'end' : 'middle'}
          >
            {(plot.start + f * (plot.end - plot.start)).toLocaleString(undefined, {
              maximumFractionDigits: 3,
            })}
          </text>
        ))}
        <text x="422" y="178" textAnchor="middle">
          {xLabel}
        </text>
      </svg>
    </div>
  );
}
