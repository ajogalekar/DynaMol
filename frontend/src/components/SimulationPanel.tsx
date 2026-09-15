import { useEffect, useRef, useState, type ReactNode } from 'react';
import {
  ArrowDownToLine,
  ArrowRight,
  Check,
  ChevronDown,
  ChevronRight,
  Clock3,
  ChartNoAxesCombined,
  Cpu,
  Droplets,
  FlaskConical,
  LoaderCircle,
  Play,
  Plus,
  RotateCcw,
  Square,
  Terminal,
  X,
  Zap,
} from 'lucide-react';
import type { Dataset, Health, Job, Measurement, SimulationConfig } from '../types';
import { api } from '../api';
import StructureWorkbench from './StructureWorkbench';
import StructureJobMonitor from './StructureJobMonitor';
import RunDiagnostics from './RunDiagnostics';
import ReadinessPanel from './ReadinessPanel';
interface RecoveryInfo {
  available: boolean;
  checkpoint_saved?: boolean;
  step?: number;
  reason: string;
}
const recovery = (job: Job) => (job as Job & { recovery?: RecoveryInfo }).recovery;
const active = (j: Job) => !['completed', 'failed', 'cancelled', 'interrupted'].includes(j.status);
const structureJob = (j: Job) => ['preparation', 'solvation'].includes(j.engine);
export default function SimulationPanel({
  dataset,
  newSetupAt,
  continuingJobs,
  engine,
  onEngineChange,
  health,
  jobs,
  viewerReady,
  viewerError,
  onClose,
  onStarted,
  onLoad,
  onRefresh,
  onDatasetLoaded,
  onWaterVisibility,
  measurements,
  onToggleTracking,
  onAddTrackedMeasurement,
  measurementEditor,
}: {
  dataset: Dataset | null;
  newSetupAt: number;
  continuingJobs: string[];
  engine: 'openmm' | 'gromacs';
  onEngineChange: (engine: 'openmm' | 'gromacs') => void;
  health: Health | null;
  jobs: Job[];
  viewerReady: boolean;
  viewerError: string;
  onClose: () => void;
  onStarted: (job: Job) => void;
  onLoad: (id: string, showWater?: boolean) => Promise<void>;
  onRefresh: () => void;
  onDatasetLoaded: (
    d: Dataset,
    options?: { showWater?: boolean; showHydrogens?: boolean },
  ) => Promise<void>;
  onWaterVisibility: (show: boolean) => void;
  measurements: Measurement[];
  onToggleTracking: (id: string, tracked: boolean) => void;
  onAddTrackedMeasurement: () => void;
  measurementEditor?: ReactNode;
}) {
  const [name, setName] = useState('My molecular journey'),
    [duration, setDuration] = useState(10),
    [temp, setTemp] = useState(300),
    [solvent, setSolvent] = useState<'implicit' | 'explicit'>('implicit'),
    [step, setStep] = useState(2),
    [interval, setInterval] = useState(50),
    [seed, setSeed] = useState(42),
    [friction, setFriction] = useState(1),
    [equil, setEquil] = useState(500),
    [padding, setPadding] = useState(1),
    [minimize, setMinimize] = useState(true),
    [advanced, setAdvanced] = useState(false),
    [busy, setBusy] = useState(false),
    [error, setError] = useState(''),
    [openLog, setOpenLog] = useState<string | null>(null);
  const [runReady, setRunReady] = useState(false);
  const [sourceBusy, setSourceBusy] = useState(false);
  const [trackingOpen, setTrackingOpen] = useState(false);
  const trackedMeasurements = measurements.filter((measurement) => measurement.trackDuringRun);
  const trackingExpanded = trackingOpen;
  const hasMeasurementEditor = !!measurementEditor;
  useEffect(() => {
    if (hasMeasurementEditor) setTrackingOpen(true);
  }, [hasMeasurementEditor]);
  const [openingHistory, setOpeningHistory] = useState<string | null>(null);
  const [resuming, setResuming] = useState<string | null>(null);
  const [recoveryDetails, setRecoveryDetails] = useState<Record<string, RecoveryInfo>>({});
  const checkedRecovery = useRef(new Set<string>());
  const recoveryFor = (job: Job) =>
    recoveryDetails[`${job.id}:${job.status}:${recovery(job)?.step ?? ''}`] ?? recovery(job);
  const [ph, setPh] = useState(dataset?.preparation?.ph ?? 7);
  const [submitting, setSubmitting] = useState<'preparation' | 'solvation' | null>(null);
  const [monitorOperation, setMonitorOperation] = useState<'preparation' | 'solvation'>('preparation');
  const [monitorId, setMonitorId] = useState<string | null>(null);
  const [dismissed, setDismissed] = useState<string[]>([]);
  const [monitorError, setMonitorError] = useState('');
  const [loadingResult, setLoadingResult] = useState(false);
  const [pending, setPending] = useState<{
    id: string;
    sourceId: string;
    operation: 'preparation' | 'solvation';
  } | null>(null);
  const datasetRef = useRef(dataset);
  datasetRef.current = dataset;
  const loadingJob = useRef<string | null>(null);
  const jobCards = useRef(new Map<string, HTMLElement>());
  const continuedJobs = useRef(new Set(continuingJobs));
  for (const job of jobs) if (structureJob(job) && active(job)) continuedJobs.current.add(job.id);
  const pendingJob = jobs.find((j) => j.id === pending?.id);
  const monitoredJob =
    jobs.find((j) => j.id === monitorId) ??
    jobs.find(
      (j) =>
        structureJob(j) &&
        !dismissed.includes(j.id) &&
        (active(j) ||
          continuedJobs.current.has(j.id) ||
          !newSetupAt ||
          Date.parse(j.created_at) >= newSetupAt) &&
        (active(j) || j.config.dataset_id === dataset?.id || j.dataset_id === dataset?.id),
    );
  useEffect(() => {
    for (const job of jobs) {
      if (
        active(job) ||
        structureJob(job) ||
        job.status === 'completed' ||
        !recovery(job)?.available
      )
        continue;
      const key = `${job.id}:${job.status}:${recovery(job)?.step ?? ''}`;
      if (checkedRecovery.current.has(key)) continue;
      checkedRecovery.current.add(key);
      void fetch(`/api/jobs/${job.id}/recovery`)
        .then(async (response) => {
          if (!response.ok) throw new Error('Could not verify this checkpoint.');
          const info = (await response.json()) as RecoveryInfo;
          setRecoveryDetails((previous) => ({ ...previous, [key]: info }));
        })
        .catch(() => checkedRecovery.current.delete(key));
    }
  }, [jobs]);
  const anyActive = jobs.some(active);
  const engineLocked =
    busy || sourceBusy || anyActive || !!pending || !!submitting || loadingResult;
  const incompatibleGromacsInput =
    engine === 'gromacs' && !!(dataset?.preparation || dataset?.solvation);
  const resultSelected = !!monitoredJob?.dataset_id && dataset?.id === monitoredJob.dataset_id;
  const complexPrepared = !!dataset?.preparation?.ligand_parameters;
  const modifiedPrepared = !!dataset?.preparation?.modified_residues?.length;
  const requiresExplicit = complexPrepared || !!dataset?.preparation?.requires_explicit_solvent;
  const paddingError =
    !Number.isFinite(padding) || padding < 1 || padding > 3
      ? 'Enter a box padding between 1 and 3 nm.'
      : '';
  const resultViewerError =
    resultSelected && viewerError ? `Could not display the prepared structure: ${viewerError}` : '';
  useEffect(() => {
    setError('');
    setMonitorError('');
    setPh(dataset?.preparation?.ph ?? 7);
    if (dataset?.monomer_selection) {
      setMonitorId(null);
      setPending(null);
      setSolvent('implicit');
    } else if (dataset?.solvation) {
      setSolvent('explicit');
      setPadding(dataset.solvation.padding_nm);
    } else if (requiresExplicit) {
      setSolvent('explicit');
    }
  }, [dataset?.id]);
  // Reconnect to a background preparation when the studio is reopened.
  useEffect(() => {
    if (pending || !dataset) return;
    const job = jobs.find(
      (j) => structureJob(j) && active(j) && j.config.dataset_id === dataset.id,
    );
    if (!job) return;
    setMonitorId(job.id);
    setPending({
      id: job.id,
      sourceId: dataset.id,
      operation: job.engine as 'preparation' | 'solvation',
    });
  }, [jobs, pending, dataset?.id]);
  useEffect(() => {
    if (!pending || !pendingJob || active(pendingJob)) return;
    const action = pending;
    if (pendingJob.status === 'completed' && pendingJob.dataset_id) {
      if (loadingJob.current === pendingJob.id) return;
      if (datasetRef.current?.id !== action.sourceId) {
        setPending(null);
        return;
      }
      loadingJob.current = pendingJob.id;
      setLoadingResult(true);
      api
        .dataset(pendingJob.dataset_id)
        .then(async (next) => {
          if (datasetRef.current?.id !== action.sourceId) return;
          await onDatasetLoaded(next, {
            showWater: action.operation === 'solvation',
            showHydrogens: action.operation === 'preparation',
          });
          if (action.operation === 'solvation') {
            setSolvent('explicit');
            onWaterVisibility(true);
          } else
            setSolvent(
              next.preparation?.ligand_parameters || next.preparation?.requires_explicit_solvent
                ? 'explicit'
                : 'implicit',
            );
        })
        .catch((e) =>
          setMonitorError(
            `The structure is prepared, but could not be opened: ${(e as Error).message}`,
          ),
        )
        .finally(() => {
          setPending(null);
          setLoadingResult(false);
          loadingJob.current = null;
        });
    } else setPending(null);
  }, [pending, pendingJob, onDatasetLoaded, onWaterVisibility]);
  function prepared(job: Job) {
    if (!dataset) return;
    setError('');
    setMonitorError('');
    setMonitorOperation('preparation');
    setMonitorId(job.id);
    setPending({ id: job.id, sourceId: dataset.id, operation: 'preparation' });
    setOpenLog(job.id);
    onStarted(job);
  }
  async function openResult() {
    if (!monitoredJob?.dataset_id || loadingResult) return;
    const sourceId = datasetRef.current?.id;
    setLoadingResult(true);
    setMonitorError('');
    try {
      const result = await api.dataset(monitoredJob.dataset_id);
      if (datasetRef.current?.id !== sourceId) return;
      await onDatasetLoaded(result, {
        showWater: monitoredJob.engine === 'solvation',
        showHydrogens: monitoredJob.engine === 'preparation',
      });
    } catch (e) {
      setMonitorError(`Could not open the prepared structure: ${(e as Error).message}`);
    } finally {
      setLoadingResult(false);
    }
  }
  async function previewWater(force = false) {
    if (!dataset || busy || engine !== 'openmm') return;
    setSolvent('explicit');
    onWaterVisibility(true);
    setError('');
    if (!dataset.preparation) {
      setError(
        'Prepare the structure above first, then build the explicit-water box. The preview will retain the selected protonation states.',
      );
      return;
    }
    if (dataset.solvation && !force) return;
    if (paddingError) return;
    setBusy(true);
    setSubmitting('solvation');
    setMonitorOperation('solvation');
    setMonitorId(null);
    setMonitorError('');
    try {
      const source = force && dataset.solvation ? dataset.solvation.parent_dataset_id : dataset.id;
      const job = await api.solvate(source, padding, ph, seed);
      setMonitorId(job.id);
      setPending({ id: job.id, sourceId: dataset.id, operation: 'solvation' });
      setOpenLog(job.id);
      onStarted(job);
    } catch (e) {
      setMonitorError((e as Error).message);
    } finally {
      setBusy(false);
      setSubmitting(null);
    }
  }
  async function changeSolvent(value: 'implicit' | 'explicit') {
    if (value === 'implicit' && requiresExplicit) return;
    if (value === 'explicit') {
      await previewWater();
      return;
    }
    setSolvent('implicit');
    onWaterVisibility(false);
    if (dataset?.solvation?.parent_dataset_id) {
      setBusy(true);
      try {
        await onDatasetLoaded(await api.dataset(dataset.solvation.parent_dataset_id), {
          showWater: false,
        });
      } catch (e) {
        setError((e as Error).message);
      } finally {
        setBusy(false);
      }
    }
  }
  const selected = health?.engines.find((e) => e.id === engine),
    steps = Math.round((duration * 1000) / step),
    frames = Math.ceil(steps / interval) + 1;
  const simulationSettings = {
    dataset_id: dataset?.id ?? '',
    engine,
    name,
    duration_ps: duration,
    temperature_k: temp,
    timestep_fs: step,
    report_interval: interval,
    friction_ps: friction,
    seed,
    solvent: engine === 'gromacs' ? 'explicit' : solvent,
    minimize,
    equilibration_steps: equil,
    padding_nm:
      engine === 'gromacs' || solvent === 'explicit'
        ? (dataset?.solvation?.padding_nm ?? padding)
        : 1,
    measurements: trackedMeasurements.map(({ id, kind, atoms, label, color }) => ({
      id,
      kind,
      atoms,
      label,
      color,
    })),
  } as SimulationConfig;
  const canStart =
    !busy &&
    !sourceBusy &&
    viewerReady &&
    !viewerError &&
    runReady &&
    !submitting &&
    !pending &&
    !anyActive &&
    !hasMeasurementEditor &&
    !incompatibleGromacsInput &&
    !((engine === 'gromacs' || solvent === 'explicit') && paddingError) &&
    !(engine === 'openmm' && solvent === 'explicit' && !dataset?.solvation) &&
    !!dataset &&
    !!selected?.available &&
    frames >= 1;
  async function start(e: React.FormEvent) {
    e.preventDefault();
    if (!dataset || !canStart) return;
    setBusy(true);
    setError('');
    try {
      const job = await api.start(simulationSettings);
      onStarted(job);
      setOpenLog(job.id);
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }
  return (
    <div className="studio-dock">
      <section
        className="simulation-drawer"
        role="dialog"
        aria-modal="false"
        aria-labelledby="simulation-title"
      >
        <header className="drawer-header">
          <div>
            <span className="eyebrow">
              <FlaskConical size={13} /> SIMULATION STUDIO
            </span>
            <h2 id="simulation-title">Set molecules in motion.</h2>
          </div>
          <button className="icon-button" onClick={onClose} aria-label="Close simulation studio">
            <X size={20} />
          </button>
        </header>
        <StructureJobMonitor
          complexPreparation={!!dataset?.atoms.some((atom) => atom.category === 'ligands')}
          job={submitting || (monitorError && !monitorId) ? undefined : monitoredJob}
          submitting={submitting}
          requestOperation={monitorOperation}
          loadingResult={loadingResult || (resultSelected && !viewerReady && !viewerError)}
          loadError={monitorError || resultViewerError}
          resultInView={resultSelected && viewerReady && !viewerError}
          onCancel={async () => {
            if (!monitoredJob) return;
            try {
              await api.cancel(monitoredJob.id);
              onRefresh();
            } catch (e) {
              setMonitorError(`Could not cancel preparation: ${(e as Error).message}`);
            }
          }}
          onOpenResult={() => void openResult()}
          onViewLog={() => {
            if (!monitoredJob) return;
            setOpenLog(monitoredJob.id);
            jobCards.current
              .get(monitoredJob.id)
              ?.scrollIntoView({ behavior: 'smooth', block: 'start' });
          }}
          onDismiss={() => {
            if (monitoredJob) setDismissed((ids) => [...ids, monitoredJob.id]);
            setMonitorId(null);
            setMonitorError('');
          }}
        />
        <div className="drawer-scroll">
          <form onSubmit={start}>
            <p className="modal-subtitle">
              A few good defaults. All the controls when you need them.
            </p>
            <div className="field-heading">
              <span>01</span>
              <h3>Choose your engine</h3>
              <span className="local-note">
                <Cpu size={11} /> Runs locally
              </span>
            </div>
            <div className="engine-options">
              {(['openmm', 'gromacs'] as const).map((id) => {
                const info = health?.engines.find((e) => e.id === id);
                return (
                  <button
                    type="button"
                    key={id}
                    className={`engine-card ${engine === id ? 'selected' : ''}`}
                    aria-label={id === 'openmm' ? 'OpenMM' : 'GROMACS'}
                    aria-pressed={engine === id}
                    onClick={() => {
                      if (id === engine) return;
                      onEngineChange(id);
                      setError('');
                      setRunReady(false);
                      if (id === 'gromacs') {
                        setSolvent('explicit');
                        onWaterVisibility(true);
                      }
                    }}
                    disabled={engineLocked}
                  >
                    <div className="engine-title">
                      {id === 'openmm' ? <Zap size={18} /> : <Cpu size={18} />}
                      <strong>{id === 'openmm' ? 'OpenMM' : 'GROMACS'}</strong>
                      <i>{engine === id && <Check size={11} />}</i>
                    </div>
                    <span>
                      {id === 'openmm'
                        ? 'Protein & complex preparation'
                        : 'Native setup for standard proteins'}
                    </span>
                    <small className={info?.available ? 'available' : ''}>
                      <b />
                      {info?.available ? `Installed · ${info.version ?? 'ready'}` : 'Not installed'}
                    </small>
                  </button>
                );
              })}
            </div>
            {!selected?.available && (
              <div className="inline-warning">
                {selected?.message ?? 'Checking the local engine…'}
              </div>
            )}
            <p className="form-note engine-preparation-note">
              {engine === 'openmm'
                ? 'Prepare with Amber ff14SB below. Recorded protonation states, supported ligand parameters and the water preview are retained by OpenMM.'
                : 'GROMACS builds its own Amber99SB-ILDN topology and TIP3P solvent when you start. Its defaults differ from the OpenMM preparation.'}
            </p>
            {incompatibleGromacsInput && (
              <div className="inline-warning engine-compatibility" role="status">
                This structure has OpenMM preparation or a saved solvent box. DynaMol cannot yet
                transfer that state to GROMACS. Use OpenMM, or load the original structure for a
                separate GROMACS setup.
                <button
                  type="button"
                  className="text-button"
                  disabled={engineLocked}
                  onClick={() => {
                    onEngineChange('openmm');
                    setRunReady(false);
                  }}
                >
                  Use OpenMM for this structure
                </button>
              </div>
            )}
            <StructureWorkbench
              dataset={dataset}
              engine={engine}
              ph={ph}
              onPh={setPh}
              seed={seed}
              locked={engineLocked}
              preparing={!!pending && pending.operation === 'preparation'}
              onDatasetLoaded={onDatasetLoaded}
              onSourceRequest={setSourceBusy}
              onPreparationStarted={prepared}
              onPreparationRequest={(starting, requestError) => {
                setSubmitting(starting ? 'preparation' : null);
                setMonitorOperation('preparation');
                setMonitorError(requestError ?? '');
                if (starting) setMonitorId(null);
              }}
            />
            <section className="water-environment" aria-labelledby="water-environment-title">
              <div className="field-heading">
                <span>03</span>
                <h3 id="water-environment-title">Water &amp; environment</h3>
                <Droplets size={16} className="water-heading-icon" />
              </div>
              <label className="full-label">
                Solvent environment
                <select
                  value={engine === 'gromacs' ? 'explicit' : solvent}
                  onChange={(e) => void changeSolvent(e.target.value as 'implicit' | 'explicit')}
                  disabled={engine === 'gromacs' || engineLocked}
                >
                  <option value="implicit" disabled={requiresExplicit}>
                    Implicit water · faster exploration
                  </option>
                  <option value="explicit">Explicit water · periodic box</option>
                </select>
              </label>
              {engine === 'openmm' && (
                <div className="solvent-preview-card">
                  <Droplets size={17} />
                  <div>
                    <b>
                      {dataset?.solvation ? 'Explicit water box is ready' : 'Build a water box'}
                    </b>
                    <span>
                      {dataset?.solvation
                        ? `${dataset.atoms.filter((a) => a.category === 'water').length.toLocaleString()} water atoms · ${dataset.solvation.water_model ?? 'TIP3P'} · ${dataset.solvation.padding_nm} nm padding`
                        : dataset?.preparation
                          ? 'Build a real solvent box and inspect it before you run.'
                          : 'Prepare the structure above to create a visible solvent box.'}
                    </span>
                  </div>
                  <button
                    type="button"
                    className="primary-button water-build-button"
                    disabled={engineLocked || !dataset?.preparation || !!paddingError}
                    onClick={() => void previewWater(true)}
                  >
                    {submitting === 'solvation' || pending?.operation === 'solvation' ? (
                      <LoaderCircle size={14} className="spin" />
                    ) : (
                      <Droplets size={14} />
                    )}
                    {dataset?.solvation ? 'Update box' : 'Build water'}
                  </button>
                </div>
              )}
              {engine === 'gromacs' && (
                <div className="solvent-preview-card" aria-label="GROMACS solvent setup">
                  <Droplets size={17} />
                  <div>
                    <b>TIP3P water · built during native setup</b>
                    <span>
                      GROMACS creates its periodic box and adds neutralizing ions after you start.
                      Follow these stages in Background activity. A separate GROMACS water preview
                      is not available.
                    </span>
                  </div>
                </div>
              )}
              <label className="water-padding">
                Box padding
                <div className="input-unit">
                  <input
                    type="number"
                    min="1"
                    max="3"
                    step="0.1"
                    value={padding}
                    onChange={(e) => setPadding(Number(e.target.value))}
                    aria-invalid={!!paddingError}
                    aria-describedby="water-padding-help"
                    required
                    disabled={engineLocked || (engine === 'openmm' && solvent === 'implicit')}
                  />
                  <span>nm</span>
                </div>
              </label>
              <p
                id="water-padding-help"
                className={paddingError ? 'error-box' : 'form-note'}
                role={paddingError ? 'alert' : undefined}
              >
                {paddingError || 'Choose 1–3 nm of water padding around the structure.'}
              </p>
            </section>
            <div className="field-heading">
              <span>04</span>
              <h3>Make it your simulation</h3>
            </div>
            <label className="full-label">
              Simulation name
              <input
                value={name}
                onChange={(e) => setName(e.target.value)}
                required
                maxLength={100}
              />
            </label>
            <div className="form-row">
              <label>
                Duration
                <div className="input-unit">
                  <input
                    type="number"
                    min="0.002"
                    max="1000"
                    step="any"
                    value={duration}
                    onChange={(e) => setDuration(Number(e.target.value))}
                    required
                  />
                  <span>ps</span>
                </div>
              </label>
              <label>
                Temperature
                <div className="input-unit">
                  <input
                    type="number"
                    min="50"
                    max="500"
                    step="any"
                    value={temp}
                    onChange={(e) => setTemp(Number(e.target.value))}
                    required
                  />
                  <span>K</span>
                </div>
              </label>
            </div>
            <p className="form-note">
              {engine === 'openmm'
                ? modifiedPrepared
                  ? `Amber ff14SB with recorded modified-residue templates${complexPrepared ? ' + GAFF2 / AM1-BCC ligands' : ''} · TIP3P water · Langevin dynamics.`
                  : complexPrepared
                    ? 'Amber ff14SB protein + GAFF2 / AM1-BCC ligands · TIP3P water · Langevin dynamics. Ligand states and parameters are retained in the saved box.'
                    : 'Amber ff14SB · Langevin dynamics. Implicit uses GBn2; explicit uses TIP3P.'
                : 'Amber99SB-ILDN · TIP3P water · stochastic dynamics.'}{' '}
              NVT exploration.{' '}
              {engine === 'openmm'
                ? 'Prepare protein–ligand complexes above, then build water and minimize before dynamics.'
                : 'Start with a complete standard protein. This preset does not support prepared ligand complexes or modified-residue parameter bundles.'}
            </p>
            <section className="tracking-section" aria-label="Track during simulation">
              <button
                type="button"
                className="tracking-disclosure"
                aria-expanded={trackingExpanded}
                aria-controls="simulation-tracking-content"
                onClick={() => setTrackingOpen(!trackingOpen)}
              >
                <ChartNoAxesCombined size={17} />
                <span>Track during simulation</span>
                <small>
                  {trackedMeasurements.length
                    ? `${trackedMeasurements.length} selected`
                    : 'Optional'}
                </small>
                {trackingExpanded ? <ChevronDown size={15} /> : <ChevronRight size={15} />}
              </button>
              {trackingExpanded && (
                <div id="simulation-tracking-content" className="tracking-content">
                  <p className="form-note">
                    Choose a distance, angle, dihedral or hydrogen bond to watch as frames are
                    saved.
                    {engine === 'gromacs' &&
                      ' Track heavy-atom geometry here. Hydrogen-bond tracking requires OpenMM’s retained hydrogen states.'}
                  </p>
                  {measurements.length > 0 && (
                    <div
                      className="tracking-measurements"
                      role="group"
                      aria-label="Measurements to track"
                    >
                      {measurements.map((measurement) => (
                        <label className="tracking-measurement" key={measurement.id}>
                          <input
                            type="checkbox"
                            checked={!!measurement.trackDuringRun}
                            disabled={
                              engineLocked ||
                              (!measurement.trackDuringRun && trackedMeasurements.length >= 12)
                            }
                            onChange={(event) =>
                              onToggleTracking(measurement.id, event.target.checked)
                            }
                            aria-label={`Track ${measurement.label}`}
                          />
                          <i style={{ background: measurement.color }} aria-hidden="true" />
                          <span>{measurement.label}</span>
                          <small>
                            {measurement.kind === 'hbond' ? 'H-bond' : measurement.kind}
                          </small>
                        </label>
                      ))}
                    </div>
                  )}
                  <button
                    type="button"
                    className="secondary-button tracking-add"
                    onClick={onAddTrackedMeasurement}
                    disabled={
                      !dataset || !viewerReady || engineLocked || trackedMeasurements.length >= 12
                    }
                  >
                    <Plus size={14} /> Add measurement
                  </button>
                  {trackedMeasurements.length >= 12 && (
                    <p className="form-note">Up to 12 measurements per simulation.</p>
                  )}
                  {measurementEditor && <div className="tracking-editor">{measurementEditor}</div>}
                </div>
              )}
            </section>
            <button
              type="button"
              className="advanced-button"
              onClick={() => setAdvanced(!advanced)}
            >
              {advanced ? <ChevronDown size={15} /> : <ChevronRight size={15} />} Advanced controls{' '}
              <span>Optional</span>
            </button>
            {advanced && (
              <div className="advanced-fields">
                <div className="form-row">
                  <label>
                    Integration step
                    <div className="input-unit">
                      <input
                        type="number"
                        min="0.1"
                        max="2"
                        step="0.1"
                        value={step}
                        onChange={(e) => setStep(Number(e.target.value))}
                        required
                      />
                      <span>fs</span>
                    </div>
                  </label>
                  <label>
                    Save every
                    <div className="input-unit">
                      <input
                        type="number"
                        min="1"
                        max="100000"
                        value={interval}
                        onChange={(e) => setInterval(Number(e.target.value))}
                        required
                      />
                      <span>steps</span>
                    </div>
                  </label>
                </div>
                <div className="form-row">
                  <label>
                    Random seed
                    <input
                      type="number"
                      min="1"
                      max="2147483646"
                      value={seed}
                      onChange={(e) => setSeed(Number(e.target.value))}
                      required
                    />
                  </label>
                  <label>
                    Friction
                    <div className="input-unit">
                      <input
                        type="number"
                        min="0.01"
                        max="100"
                        step="any"
                        value={friction}
                        onChange={(e) => setFriction(Number(e.target.value))}
                        required
                      />
                      <span>ps⁻¹</span>
                    </div>
                  </label>
                </div>
                <div className="form-row">
                  <label>
                    Equilibration
                    <div className="input-unit">
                      <input
                        type="number"
                        min="0"
                        max="100000"
                        value={equil}
                        onChange={(e) => setEquil(Number(e.target.value))}
                        required
                      />
                      <span>steps</span>
                    </div>
                  </label>
                </div>
                <label className="checkbox-label">
                  <input
                    type="checkbox"
                    checked={minimize}
                    onChange={(e) => setMinimize(e.target.checked)}
                  />{' '}
                  Minimize energy before equilibration
                </label>
              </div>
            )}
            <div className="run-summary">
              <span>
                <Clock3 size={14} />
                {steps.toLocaleString()} production steps
              </span>
              <span>~{frames.toLocaleString()} saved frames</span>
            </div>
            <p className="form-note">
              Starts from the first frame. Initial relaxation precedes production. Native
              checkpoints are saved automatically during production; Resume retains the original
              configuration. Short runs do not establish convergence.
            </p>
            {error && (
              <div className="error-box" role="alert">
                {error}
              </div>
            )}
            {dataset && (
              <ReadinessPanel
                mode="simulation"
                datasetId={dataset.id}
                settings={simulationSettings}
                onReadyChange={setRunReady}
                suspended={anyActive}
              />
            )}
            <button className="primary-button full-width run-button" disabled={!canStart}>
              {busy ? (
                <LoaderCircle size={17} className="spin" />
              ) : (
                <Play size={16} fill="currentColor" />
              )}
              {busy ? 'Starting simulation…' : 'Start simulation'}
              <ArrowRight size={17} />
            </button>
            <p className="under-button">
              {hasMeasurementEditor
                ? 'Finish or cancel the measurement selection before starting.'
                : engine === 'openmm' && solvent === 'explicit' && !dataset?.solvation
                  ? 'Prepare the structure and build its water preview before starting.'
                  : 'Keep exploring while your simulation runs in the background.'}
            </p>
          </form>
          <section className="jobs-section">
            <div className="field-heading">
              <span>
                <ActivityIcon />
              </span>
              <h3>Background activity</h3>
              <b>{jobs.length}</b>
            </div>
            {jobs.length === 0 ? (
              <div className="jobs-empty">
                Your runs will appear here, with live progress and logs.
              </div>
            ) : (
              jobs.map((job) => (
                <article
                  className="job-card"
                  key={job.id}
                  ref={(node) => {
                    if (node) jobCards.current.set(job.id, node);
                    else jobCards.current.delete(job.id);
                  }}
                >
                  <div className="job-top">
                    <div>
                      <strong>{job.name}</strong>
                      <span>
                        {job.engine} · {job.stage ?? job.status}
                      </span>
                    </div>
                    <span className={`job-status status-${job.status}`}>
                      {active(job) && <LoaderCircle size={11} className="spin" />}
                      {job.status}
                    </span>
                  </div>
                  <div className="progress-track">
                    <div style={{ width: `${Math.max(0, Math.min(100, job.progress))}%` }} />
                  </div>
                  <div className="job-progress">
                    <span>
                      {job.completed_steps.toLocaleString()} / {job.total_steps.toLocaleString()}{' '}
                      {job.engine === 'preparation' || job.engine === 'solvation'
                        ? 'stages'
                        : 'steps'}
                    </span>
                    <b>{Math.round(job.progress)}%</b>
                  </div>
                  {job.error && <div className="error-box">{job.error}</div>}
                  <div className="job-actions">
                    <button
                      type="button"
                      className="text-button"
                      onClick={() => setOpenLog(openLog === job.id ? null : job.id)}
                    >
                      <Terminal size={13} /> {openLog === job.id ? 'Hide' : 'View'} log
                    </button>
                    {!active(job) && (
                      <a className="text-button" href={`/api/jobs/${job.id}/download`}>
                        <ArrowDownToLine size={13} /> Files
                      </a>
                    )}
                    {!active(job) && !structureJob(job) && job.status !== 'completed' && (
                      <button
                        type="button"
                        className="text-button accent"
                        disabled={anyActive || resuming === job.id || !recoveryFor(job)?.available}
                        title={
                          recoveryFor(job)?.reason ??
                          'No native checkpoint is available for this job.'
                        }
                        onClick={async () => {
                          setResuming(job.id);
                          setError('');
                          try {
                            const response = await fetch(`/api/jobs/${job.id}/resume`, {
                              method: 'POST',
                            });
                            const result = await response.json();
                            if (!response.ok)
                              throw new Error(
                                typeof result.detail === 'string'
                                  ? result.detail
                                  : 'Could not resume this checkpoint.',
                              );
                            onStarted(result as Job);
                            setOpenLog(job.id);
                          } catch (e) {
                            setError((e as Error).message);
                          } finally {
                            setResuming(null);
                          }
                        }}
                      >
                        {resuming === job.id ? (
                          <LoaderCircle size={13} className="spin" />
                        ) : (
                          <RotateCcw size={13} />
                        )}
                        {resuming === job.id ? 'Resuming…' : 'Resume'}
                      </button>
                    )}
                    {active(job) ? (
                      <button
                        type="button"
                        className="text-button danger"
                        disabled={job.status === 'cancelling'}
                        onClick={async () => {
                          try {
                            await api.cancel(job.id);
                            onRefresh();
                          } catch (e) {
                            setError((e as Error).message);
                          }
                        }}
                      >
                        <Square size={11} /> {job.status === 'cancelling' ? 'Stopping…' : 'Stop'}
                      </button>
                    ) : job.status === 'completed' && job.dataset_id ? (
                      <>
                        <button
                          type="button"
                          className="text-button accent"
                          disabled={openingHistory !== null}
                          onClick={async () => {
                            setOpeningHistory(job.id);
                            try {
                              await onLoad(
                                job.dataset_id!,
                                job.config.solvent === 'explicit' || job.engine === 'solvation',
                              );
                            } finally {
                              setOpeningHistory(null);
                            }
                          }}
                        >
                          {openingHistory === job.id
                            ? 'Opening…'
                            : structureJob(job)
                              ? 'Open structure'
                              : 'Open trajectory'}{' '}
                          {openingHistory === job.id ? (
                            <LoaderCircle size={13} className="spin" />
                          ) : (
                            <ArrowRight size={13} />
                          )}
                        </button>
                      </>
                    ) : null}
                  </div>
                  {!active(job) && !structureJob(job) && job.status !== 'completed' && (
                    <p className="form-note">
                      {recoveryFor(job)?.reason ??
                        'No production checkpoint is available. Files still contains the saved inputs and logs.'}
                    </p>
                  )}
                  {!structureJob(job) && (active(job) || openLog === job.id) && (
                    <RunDiagnostics job={job} />
                  )}
                  {openLog === job.id && (
                    <pre className="job-log">
                      {job.logs.length ? job.logs.join('\n') : 'Waiting for engine output…'}
                    </pre>
                  )}
                </article>
              ))
            )}
          </section>
        </div>
      </section>
    </div>
  );
}
function ActivityIcon() {
  return <Cpu size={14} />;
}
