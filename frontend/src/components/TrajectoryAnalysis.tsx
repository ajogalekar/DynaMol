import { useEffect, useRef, useState } from 'react';
import { Activity, Download, LoaderCircle } from 'lucide-react';
import type { Dataset } from '../types';
import { scientificRequest, ScientificChart, saveCsv, saveJson } from './ScientificCharts';
export interface AnalysisSettings {
  kind: 'rmsd' | 'rmsf';
  selection: string;
  alignment: 'ca' | 'backbone' | 'selection' | 'none';
  reference: number;
  start: number;
  end: number;
  stride: number;
  periodic: 'whole' | 'cartesian';
  residue_average: boolean;
}
interface Result {
  dataset_id: string;
  kind: 'rmsd' | 'rmsf';
  values: number[];
  x: number[];
  x_label: string;
  labels: string[];
  frame_indices: number[];
  atom_groups: number[][] | null;
  method: string;
  warnings: string[];
  summary: { min: number; max: number; mean: number; frames: number; atoms: number };
}
interface NamedSelection {
  id: string;
  name: string;
  atoms: number[];
}
export default function TrajectoryAnalysis({
  dataset,
  frame,
  selectedAtoms,
  namedSelections = [],
  onSeek,
  onSelectAtoms,
  savedState,
  onStateChange,
  expanded,
  onExpandedChange,
}: {
  dataset: Dataset | null;
  frame: number;
  selectedAtoms: number[];
  namedSelections?: NamedSelection[];
  onSeek: (frame: number) => void;
  onSelectAtoms?: (atoms: number[]) => void;
  savedState?: AnalysisSettings | null;
  onStateChange?: (state: AnalysisSettings) => void;
  expanded?: boolean;
  onExpandedChange?: (open: boolean) => void;
}) {
  const defaults = (): AnalysisSettings => ({
    kind: 'rmsd',
    selection: 'ca',
    alignment: 'ca',
    reference: 1,
    start: 1,
    end: dataset?.n_frames ?? 1,
    stride: 1,
    periodic: 'whole',
    residue_average: true,
  });
  const [settings, setSettings] = useState<AnalysisSettings>(savedState ?? defaults),
    [result, setResult] = useState<Result | null>(null),
    [busy, setBusy] = useState(false),
    [error, setError] = useState('');
  const controller = useRef<AbortController | null>(null),
    currentDataset = useRef(dataset?.id);
  const onChange = useRef(onStateChange);
  onChange.current = onStateChange;
  useEffect(() => {
    controller.current?.abort();
    currentDataset.current = dataset?.id;
    setSettings(savedState ?? defaults());
    setResult(null);
    setBusy(false);
    setError('');
  }, [dataset?.id]);
  useEffect(() => () => controller.current?.abort(), []);
  const update = (values: Partial<AnalysisSettings>) => {
    controller.current?.abort();
    setBusy(false);
    const next = { ...settings, ...values };
    setSettings(next);
    setResult(null);
    onChange.current?.(next);
  };
  useEffect(() => {
    if (
      settings.selection.startsWith('named:') &&
      !namedSelections.some((s) => `named:${s.id}` === settings.selection)
    )
      update({ selection: 'ca' });
  }, [namedSelections, settings.selection]);
  async function calculate() {
    if (!dataset) return;
    controller.current?.abort();
    const abort = new AbortController();
    controller.current = abort;
    setBusy(true);
    setError('');
    setResult(null);
    const custom =
      settings.selection === 'current'
        ? selectedAtoms
        : namedSelections.find((s) => `named:${s.id}` === settings.selection)?.atoms;
    try {
      const data = await scientificRequest<Result>(
        `/api/datasets/${dataset.id}/structural-analysis`,
        {
          kind: settings.kind,
          selection: custom ? 'custom' : settings.selection,
          atoms: custom ?? [],
          alignment: settings.alignment,
          reference_frame: settings.reference - 1,
          start_frame: settings.start - 1,
          end_frame: settings.end - 1,
          stride: settings.stride,
          periodic: settings.periodic,
          residue_average: settings.residue_average,
        },
        abort.signal,
      );
      if (!abort.signal.aborted && currentDataset.current === dataset.id) setResult(data);
    } catch (e) {
      if (!abort.signal.aborted) setError((e as Error).message);
    } finally {
      if (!abort.signal.aborted) setBusy(false);
    }
  }
  if (!dataset) return null;
  const currentIndex =
    result?.kind === 'rmsd' ? result.frame_indices.indexOf(Math.round(frame)) : undefined;
  return (
    <details
      className="trajectory-analysis"
      open={expanded}
      onToggle={(e) => onExpandedChange?.(e.currentTarget.open)}
    >
      <summary>
        <Activity size={16} /> Structural analysis <small>RMSD · RMSF · saved selections</small>
      </summary>
      <div className="trajectory-analysis-body">
        <div className="analysis-controls">
          <label>
            Analysis
            <select
              aria-label="Structural analysis kind"
              value={settings.kind}
              onChange={(e) => update({ kind: e.target.value as AnalysisSettings['kind'] })}
            >
              <option value="rmsd">RMSD · change over time</option>
              <option value="rmsf">RMSF · local fluctuations</option>
            </select>
          </label>
          <label>
            Measure atoms
            <select
              aria-label="Analysis atom group"
              value={settings.selection}
              onChange={(e) => update({ selection: e.target.value })}
            >
              <option value="ca">Protein Cα</option>
              <option value="backbone">Protein backbone</option>
              <option value="protein-heavy">Protein heavy atoms</option>
              <option value="solute-heavy">All solute heavy atoms</option>
              <option value="current" disabled={!selectedAtoms.length}>
                Current selected atoms ({selectedAtoms.length})
              </option>
              {namedSelections.map((s) => (
                <option key={s.id} value={`named:${s.id}`}>
                  {s.name} ({s.atoms.length} atoms)
                </option>
              ))}
            </select>
          </label>
          <label>
            Fit before measuring
            <select
              aria-label="Analysis alignment"
              value={settings.alignment}
              onChange={(e) =>
                update({ alignment: e.target.value as AnalysisSettings['alignment'] })
              }
            >
              <option value="ca">Fit protein Cα</option>
              <option value="backbone">Fit protein backbone</option>
              <option value="selection">Fit measured atoms</option>
              <option value="none">No fitting</option>
            </select>
          </label>
          <label>
            Periodic coordinates
            <select
              aria-label="Analysis periodic handling"
              value={settings.periodic}
              onChange={(e) => update({ periodic: e.target.value as AnalysisSettings['periodic'] })}
            >
              <option value="whole">Make molecules whole</option>
              <option value="cartesian">Saved Cartesian coordinates</option>
            </select>
          </label>
          <label>
            Reference frame
            <input
              aria-label="Analysis reference frame"
              type="number"
              min="1"
              max={dataset.n_frames}
              value={settings.reference}
              onChange={(e) => update({ reference: Number(e.target.value) })}
            />
          </label>
          <label>
            First frame
            <input
              aria-label="Analysis first frame"
              type="number"
              min="1"
              max={dataset.n_frames}
              value={settings.start}
              onChange={(e) => update({ start: Number(e.target.value) })}
            />
          </label>
          <label>
            Last frame
            <input
              aria-label="Analysis last frame"
              type="number"
              min="1"
              max={dataset.n_frames}
              value={settings.end}
              onChange={(e) => update({ end: Number(e.target.value) })}
            />
          </label>
          <label>
            Analyze every N frames
            <input
              aria-label="Analysis stride"
              type="number"
              min="1"
              max={dataset.n_frames}
              value={settings.stride}
              onChange={(e) => update({ stride: Number(e.target.value) })}
            />
          </label>
        </div>
        <div className="analysis-actions">
          <button
            type="button"
            className="primary-button"
            disabled={busy || (settings.kind === 'rmsf' && dataset.n_frames < 2)}
            onClick={() => void calculate()}
          >
            {busy ? <LoaderCircle size={14} className="spin" /> : <Activity size={14} />}
            {busy ? 'Analyzing saved frames…' : `Calculate ${settings.kind.toUpperCase()}`}
          </button>
          {settings.kind === 'rmsf' && (
            <label className="analysis-note">
              <input
                type="checkbox"
                checked={settings.residue_average}
                onChange={(e) => update({ residue_average: e.target.checked })}
              />{' '}
              Group by residue
            </label>
          )}
          {result && (
            <button
              type="button"
              className="text-button"
              onClick={() =>
                saveCsv(`DynaMol-${dataset.name}-${result.kind}.csv`, [
                  [result.x_label, 'label', `${result.kind}_angstrom`],
                  ...result.values.map((v, i) => [result.x[i], result.labels[i], v]),
                ])
              }
            >
              <Download size={14} /> Export CSV
            </button>
          )}
          {result && (
            <button
              type="button"
              className="text-button"
              onClick={() =>
                saveJson(`DynaMol-${dataset.name}-${result.kind}-record.json`, {
                  format: 'DynaMol structural analysis',
                  version: 1,
                  dataset: { id: dataset.id, name: dataset.name },
                  analysis: result,
                })
              }
            >
              <Download size={14} /> Analysis record
            </button>
          )}
        </div>
        {error && (
          <p className="analysis-warning" role="alert">
            {error}
          </p>
        )}
        {result && (
          <>
            <div className="analysis-summary">
              <span>{result.summary.atoms.toLocaleString()} measured atoms</span>
              <span>{result.summary.frames.toLocaleString()} frames</span>
              <span>Mean {result.summary.mean.toFixed(3)} Å</span>
              <span>Maximum {result.summary.max.toFixed(3)} Å</span>
            </div>
            <ScientificChart
              x={result.x}
              values={result.values}
              title={result.kind.toUpperCase()}
              unit="Å"
              xLabel={result.x_label}
              labels={result.labels}
              currentIndex={
                currentIndex !== undefined && currentIndex >= 0 ? currentIndex : undefined
              }
              onPick={(i) =>
                result.kind === 'rmsd'
                  ? onSeek(result.frame_indices[i])
                  : result.atom_groups && onSelectAtoms?.(result.atom_groups[i])
              }
            />
            <p className="analysis-note">
              {result.method}{' '}
              {result.kind === 'rmsd'
                ? 'Click the plot to seek the trajectory.'
                : 'Click a residue to select its measured atoms.'}
            </p>
            {result.warnings.map((w) => (
              <p className="analysis-note" key={w}>
                {w}
              </p>
            ))}
          </>
        )}
        {!result && (
          <p className="analysis-note">
            Uses original saved coordinates, independent of the camera and display animation. RMSF
            is measured about the mean structure over the chosen frame window. Values describe this
            trajectory; they do not establish convergence.
          </p>
        )}
      </div>
    </details>
  );
}
