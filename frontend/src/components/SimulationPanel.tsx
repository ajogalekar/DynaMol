import { useState } from 'react';
import {
  ArrowDownToLine,
  ArrowRight,
  Check,
  ChevronDown,
  ChevronRight,
  Clock3,
  Cpu,
  FlaskConical,
  LoaderCircle,
  Play,
  Square,
  Terminal,
  X,
  Zap,
} from 'lucide-react';
import type { Dataset, Health, Job, SimulationConfig } from '../types';
import { api } from '../api';
const active = (j: Job) => !['completed', 'failed', 'cancelled', 'interrupted'].includes(j.status);
export default function SimulationPanel({
  dataset,
  health,
  jobs,
  onClose,
  onStarted,
  onLoad,
  onRefresh,
}: {
  dataset: Dataset | null;
  health: Health | null;
  jobs: Job[];
  onClose: () => void;
  onStarted: (job: Job) => void;
  onLoad: (id: string) => void;
  onRefresh: () => void;
}) {
  const [engine, setEngine] = useState<'openmm' | 'gromacs'>('openmm'),
    [name, setName] = useState('My molecular journey'),
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
  const selected = health?.engines.find((e) => e.id === engine),
    steps = Math.round((duration * 1000) / step),
    frames = Math.ceil(steps / interval) + 1;
  async function start(e: React.FormEvent) {
    e.preventDefault();
    if (!dataset) return;
    setBusy(true);
    setError('');
    try {
      const job = await api.start({
        dataset_id: dataset.id,
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
        padding_nm: padding,
      } as SimulationConfig);
      onStarted(job);
      setOpenLog(job.id);
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }
  return (
    <div
      className="modal-backdrop drawer-backdrop"
      onMouseDown={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <section
        className="simulation-drawer"
        role="dialog"
        aria-modal="true"
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
        <div className="drawer-scroll">
          <form onSubmit={start}>
            <p className="modal-subtitle">
              A few good defaults. All the controls when you need them.
            </p>
            <div className="source-summary">
              <FileMolecule />
              <div>
                <span>Starting structure · first saved frame</span>
                <strong>{dataset?.name ?? 'Open a structure to begin'}</strong>
                <small>
                  {dataset
                    ? `${dataset.n_atoms.toLocaleString()} atoms · ${dataset.n_residues} residues`
                    : ''}
                </small>
              </div>
              <Check size={15} />
            </div>
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
                    onClick={() => {
                      setEngine(id);
                      if (id === 'gromacs') setSolvent('explicit');
                    }}
                  >
                    <div className="engine-title">
                      {id === 'openmm' ? <Zap size={18} /> : <Cpu size={18} />}
                      <strong>{id === 'openmm' ? 'OpenMM' : 'GROMACS'}</strong>
                      <i>{engine === id && <Check size={11} />}</i>
                    </div>
                    <span>
                      {id === 'openmm'
                        ? 'Flexible. Fast to get started.'
                        : 'The established MD workhorse.'}
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
            <div className="field-heading">
              <span>02</span>
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
            <label className="full-label">
              Solvent environment
              <select
                value={engine === 'gromacs' ? 'explicit' : solvent}
                onChange={(e) => setSolvent(e.target.value as 'implicit' | 'explicit')}
                disabled={engine === 'gromacs'}
              >
                <option value="implicit">Implicit water · faster exploration</option>
                <option value="explicit">Explicit water · periodic box</option>
              </select>
            </label>
            <p className="form-note">
              {engine === 'openmm'
                ? 'Amber ff14SB · Langevin dynamics. Implicit uses GBn2; explicit uses TIP3P.'
                : 'Amber99SB-ILDN · TIP3P water · stochastic dynamics.'}{' '}
              NVT exploration with standard protein residues. Ligands and unsupported chemistry
              require separate parameterization.
            </p>
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
                  <label>
                    Box padding
                    <div className="input-unit">
                      <input
                        type="number"
                        min="1"
                        max="3"
                        step="0.1"
                        value={padding}
                        onChange={(e) => setPadding(Number(e.target.value))}
                        required
                        disabled={engine === 'openmm' && solvent === 'implicit'}
                      />
                      <span>nm</span>
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
              Starts from the first frame. Equilibration precedes production; short runs are for
              exploration and do not establish convergence. Configuration, seed, logs, and outputs
              are saved with each job.
            </p>
            {error && (
              <div className="error-box" role="alert">
                {error}
              </div>
            )}
            <button
              className="primary-button full-width run-button"
              disabled={busy || !dataset || !selected?.available || frames < 1}
            >
              {busy ? (
                <LoaderCircle size={17} className="spin" />
              ) : (
                <Play size={16} fill="currentColor" />
              )}
              {busy ? 'Starting simulation…' : 'Start simulation'}
              <ArrowRight size={17} />
            </button>
            <p className="under-button">
              Keep exploring while your simulation runs in the background.
            </p>
          </form>
          <section className="jobs-section">
            <div className="field-heading">
              <span>
                <ActivityIcon />
              </span>
              <h3>Simulation activity</h3>
              <b>{jobs.length}</b>
            </div>
            {jobs.length === 0 ? (
              <div className="jobs-empty">
                Your runs will appear here, with live progress and logs.
              </div>
            ) : (
              jobs.map((job) => (
                <article className="job-card" key={job.id}>
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
                      steps
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
                    {active(job) ? (
                      <button
                        type="button"
                        className="text-button danger"
                        onClick={async () => {
                          try {
                            await api.cancel(job.id);
                            onRefresh();
                          } catch (e) {
                            setError((e as Error).message);
                          }
                        }}
                      >
                        <Square size={11} /> Stop
                      </button>
                    ) : job.status === 'completed' && job.dataset_id ? (
                      <>
                        <a className="text-button" href={`/api/jobs/${job.id}/download`}>
                          <ArrowDownToLine size={13} /> Files
                        </a>
                        <button
                          type="button"
                          className="text-button accent"
                          onClick={() => onLoad(job.dataset_id!)}
                        >
                          Open trajectory <ArrowRight size={13} />
                        </button>
                      </>
                    ) : null}
                  </div>
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
function FileMolecule() {
  return (
    <div className="source-icon">
      <FlaskConical size={22} />
    </div>
  );
}
function ActivityIcon() {
  return <Cpu size={14} />;
}
