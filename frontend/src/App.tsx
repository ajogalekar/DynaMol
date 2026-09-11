import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import {
  Activity,
  ArrowLeft,
  ArrowRight,
  Atom,
  BookOpen,
  Box,
  Camera,
  Check,
  ChevronDown,
  ChevronRight,
  CircleHelp,
  Crosshair,
  Droplets,
  Expand,
  Eye,
  EyeOff,
  FlaskConical,
  Focus,
  FolderOpen,
  Info,
  Layers3,
  LoaderCircle,
  Minus,
  MousePointer2,
  MoveUpRight,
  Pause,
  Play,
  Plus,
  Repeat2,
  RotateCcw,
  Search,
  Settings2,
  SkipBack,
  SkipForward,
  Sparkles,
  Terminal,
  X,
  Zap,
} from 'lucide-react';
import MolecularViewer from './components/MolecularViewer';
import type { ViewerHandle } from './components/MolecularViewer';
import PlotPanel from './components/PlotPanel';
import ImportDialog from './components/ImportDialog';
import SimulationPanel from './components/SimulationPanel';
import { api } from './api';
import type {
  AtomInfo,
  Dataset,
  Health,
  Job,
  MeasureKind,
  Measurement,
  Representation,
  Visibility,
} from './types';

const palette = ['#74dfc5', '#e9b970', '#ad9dff', '#7fc5ff', '#f493a1'];
const atomLabel = (a: AtomInfo) => `${a.residue}${a.resid} · ${a.name}`;
const measureConfig: Record<
  MeasureKind,
  { name: string; count: number; hint: string; short: string }
> = {
  distance: {
    name: 'Distance',
    count: 2,
    hint: 'Pick any two atoms to follow their separation.',
    short: 'A — B',
  },
  angle: {
    name: 'Angle',
    count: 3,
    hint: 'Pick three atoms in order. The second is the vertex.',
    short: 'A — B — C',
  },
  dihedral: {
    name: 'Dihedral',
    count: 4,
    hint: 'Pick four atoms in bond order to track their torsion.',
    short: 'A — B — C — D',
  },
  hbond: {
    name: 'Hydrogen bond',
    count: 3,
    hint: 'Pick donor, its bonded hydrogen, then acceptor. Shows D–A distance and geometric occupancy (≤3.5 Å, ≥150°).',
    short: 'D — H ··· A',
  },
};
const isActive = (j: Job) =>
  !['completed', 'failed', 'cancelled', 'interrupted'].includes(j.status);

export default function App() {
  const [dataset, setDataset] = useState<Dataset | null>(null),
    [coordinates, setCoordinates] = useState<Float32Array | null>(null),
    [datasets, setDatasets] = useState<Dataset[]>([]),
    [health, setHealth] = useState<Health | null>(null),
    [jobs, setJobs] = useState<Job[]>([]);
  const [loading, setLoading] = useState(true),
    [loadingMessage, setLoadingMessage] = useState('Opening your molecular workspace…'),
    [error, setError] = useState(''),
    [toast, setToast] = useState(''),
    [ready, setReady] = useState(false);
  const [frame, setFrame] = useState(0),
    [playing, setPlaying] = useState(false),
    [speed, setSpeed] = useState(1),
    [loop, setLoop] = useState(true),
    [spin, setSpin] = useState(false);
  const [visibility, setVisibility] = useState<Visibility>({
      protein: true,
      water: false,
      ligands: true,
      ions: true,
      hydrogens: 'none',
    }),
    [representation, setRepresentation] = useState<Representation>('cartoon'),
    [colorScheme, setColorScheme] = useState<'chain' | 'residue' | 'element'>('residue');
  const [measurements, setMeasurements] = useState<Measurement[]>([]),
    [activeMeasurement, setActiveMeasurement] = useState<string | null>(null),
    [kind, setKind] = useState<MeasureKind>('distance'),
    [picking, setPicking] = useState(false),
    [selectedAtoms, setSelectedAtoms] = useState<number[]>([]),
    [measureBusy, setMeasureBusy] = useState(false),
    [atomSearch, setAtomSearch] = useState('');
  const [modal, setModal] = useState<'import' | 'simulation' | 'help' | null>(null),
    [library, setLibrary] = useState(false),
    [details, setDetails] = useState(false),
    [inspector, setInspector] = useState(true),
    [focusQuery, setFocusQuery] = useState('');
  const viewer = useRef<ViewerHandle>(null),
    frameRef = useRef(0),
    loadToken = useRef(0),
    workspace = useRef<HTMLDivElement>(null);
  const visibleMeasurements = useMemo(
    () => measurements.filter((m) => m.visible !== false),
    [measurements],
  );
  const time = dataset?.times_ps[Math.min((dataset?.n_frames ?? 1) - 1, Math.round(frame))] ?? 0;
  const refreshJobs = useCallback(() => {
    api
      .jobs()
      .then(setJobs)
      .catch(() => {});
  }, []);
  const refreshLibrary = useCallback(() => {
    api
      .datasets()
      .then(setDatasets)
      .catch(() => {});
  }, []);

  const loadDataset = useCallback(
    async (
      d: Dataset,
      options: { keepStudio?: boolean; showWater?: boolean; showHydrogens?: boolean } = {},
    ) => {
      const token = ++loadToken.current;
      setLoading(true);
      setReady(false);
      setLoadingMessage(`Loading ${d.name}…`);
      setPlaying(false);
      setError('');
      setMeasurements([]);
      setMeasureBusy(false);
      setActiveMeasurement(null);
      setSelectedAtoms([]);
      setPicking(false);
      setFrame(0);
      frameRef.current = 0;
      try {
        const coords = await api.coordinates(d);
        if (token !== loadToken.current) return;
        setDataset(d);
        setCoordinates(coords);
        setVisibility({
          protein: true,
          water: options.showWater ?? !!d.solvation,
          ligands: true,
          ions: true,
          hydrogens: options.showHydrogens ? 'polar' : 'none',
        });
        setRepresentation('cartoon');
        if (!options.keepStudio) setModal(null);
        setLibrary(false);
        refreshLibrary();
        const ca = d.atoms.filter((a) => a.name === 'CA' && a.category === 'protein');
        if (ca.length > 5 && d.n_frames > 1) {
          const atoms = [
            ca[Math.floor(ca.length * 0.27)].index,
            ca[Math.floor(ca.length * 0.67)].index,
          ];
          api
            .measure(d.id, 'distance', atoms)
            .then((m) => {
              if (token !== loadToken.current) return;
              const id = crypto.randomUUID();
              setMeasurements([
                {
                  ...m,
                  id,
                  color: palette[0],
                  visible: false,
                  label: `${d.atoms[atoms[0]].residue}${d.atoms[atoms[0]].resid} ↔ ${d.atoms[atoms[1]].residue}${d.atoms[atoms[1]].resid}`,
                },
              ]);
              setActiveMeasurement(id);
            })
            .catch(() => {});
        }
      } catch (e) {
        if (token === loadToken.current) setError((e as Error).message);
      } finally {
        if (token === loadToken.current) setLoading(false);
      }
    },
    [refreshLibrary],
  );
  useEffect(() => {
    let live = true;
    api
      .health()
      .then((h) => {
        if (live) setHealth(h);
      })
      .catch(() => {
        if (live)
          setError(
            'The local molecular service is offline. Start DynaMol with ./start.sh, then reload this page.',
          );
      });
    api
      .demo()
      .then((d) => {
        if (live) void loadDataset(d);
      })
      .catch((e) => {
        if (live) {
          setError((e as Error).message);
          setLoading(false);
        }
      });
    refreshJobs();
    refreshLibrary();
    const timer = window.setInterval(() => {
      refreshJobs();
      api
        .health()
        .then((h) => {
          if (live) setHealth(h);
        })
        .catch(() => {});
    }, 2500);
    return () => {
      live = false;
      window.clearInterval(timer);
      ++loadToken.current;
    };
  }, [loadDataset, refreshJobs, refreshLibrary]);
  useEffect(() => {
    if (!toast) return;
    const t = window.setTimeout(() => setToast(''), 4500);
    return () => clearTimeout(t);
  }, [toast]);
  useEffect(() => {
    frameRef.current = frame;
  }, [frame]);
  useEffect(() => {
    if (!playing || !dataset || dataset.n_frames < 2) return;
    let request = 0,
      last = 0;
    const animate = (now: number) => {
      if (!last) last = now;
      const elapsed = Math.min(0.1, (now - last) / 1000);
      last = now;
      let next = frameRef.current + elapsed * 12 * speed;
      if (next >= dataset.n_frames - 1) {
        if (loop) {
          next = 0;
        } else {
          next = dataset.n_frames - 1;
          setPlaying(false);
        }
      }
      frameRef.current = next;
      setFrame(next);
      request = requestAnimationFrame(animate);
    };
    request = requestAnimationFrame(animate);
    return () => cancelAnimationFrame(request);
  }, [playing, dataset, speed, loop]);
  const seek = useCallback(
    (f: number) => {
      setPlaying(false);
      setFrame(Math.max(0, Math.min((dataset?.n_frames ?? 1) - 1, f)));
    },
    [dataset],
  );
  const togglePlayback = useCallback(() => {
    if (!dataset || dataset.n_frames < 2 || !ready) return;
    if (frameRef.current >= dataset.n_frames - 1) setFrame(0);
    setPlaying((p) => !p);
  }, [dataset, ready]);
  useEffect(() => {
    const key = (e: KeyboardEvent) => {
      const tag = (e.target as HTMLElement).tagName;
      if (
        ['INPUT', 'SELECT', 'TEXTAREA', 'BUTTON'].includes(tag) ||
        e.metaKey ||
        e.ctrlKey ||
        e.altKey
      )
        return;
      if (e.key === 'Escape') {
        setModal(null);
        setPicking(false);
        setSelectedAtoms([]);
        setLibrary(false);
      }
      if (modal && modal !== 'simulation') return;
      if (e.code === 'Space') {
        e.preventDefault();
        togglePlayback();
      }
      if (e.key === 'ArrowRight') {
        e.preventDefault();
        seek(Math.round(frameRef.current) + 1);
      }
      if (e.key === 'ArrowLeft') {
        e.preventDefault();
        seek(Math.round(frameRef.current) - 1);
      }
      if (e.key.toLowerCase() === 'f') viewer.current?.fit();
      if (e.key.toLowerCase() === 'm') setPicking((p) => !p);
      if (e.key === '?') setModal('help');
    };
    window.addEventListener('keydown', key);
    return () => window.removeEventListener('keydown', key);
  }, [modal, seek, togglePlayback]);
  useEffect(() => {
    if (!modal || modal === 'simulation') return;
    const previous = document.activeElement as HTMLElement | null;
    const dialog = document.querySelector<HTMLElement>('[role="dialog"]');
    const getFocusable = () =>
      Array.from(
        dialog?.querySelectorAll<HTMLElement>(
          'button:not(:disabled), input:not([hidden]):not(:disabled), select:not(:disabled), a[href], [tabindex="0"]',
        ) ?? [],
      ).filter((e) => e.offsetParent !== null);
    getFocusable()[0]?.focus();
    const trap = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        setModal(null);
        return;
      }
      if (e.key !== 'Tab') return;
      const elements = getFocusable();
      const first = elements[0],
        last = elements[elements.length - 1];
      if (e.shiftKey && document.activeElement === first) {
        e.preventDefault();
        last?.focus();
      } else if (!e.shiftKey && document.activeElement === last) {
        e.preventDefault();
        first?.focus();
      }
    };
    document.addEventListener('keydown', trap);
    return () => {
      document.removeEventListener('keydown', trap);
      previous?.focus();
    };
  }, [modal]);
  const onAtomPick = useCallback(
    (index: number) => {
      if (!picking) return;
      setSelectedAtoms((prev) => {
        if (prev.includes(index)) return prev.filter((i) => i !== index);
        return prev.length >= measureConfig[kind].count ? [index] : [...prev, index];
      });
    },
    [picking, kind],
  );
  async function measure() {
    if (!dataset || selectedAtoms.length !== measureConfig[kind].count) return;
    setMeasureBusy(true);
    setError('');
    const currentLoadToken = loadToken.current;
    try {
      const m = await api.measure(dataset.id, kind, selectedAtoms);
      if (currentLoadToken !== loadToken.current) return;
      const id = crypto.randomUUID();
      const label = selectedAtoms
        .map((i) => `${dataset.atoms[i].residue}${dataset.atoms[i].resid}:${dataset.atoms[i].name}`)
        .join(kind === 'distance' ? ' ↔ ' : ' · ');
      setMeasurements((prev) => [
        ...prev,
        { ...m, id, label, color: palette[prev.length % palette.length] },
      ]);
      setActiveMeasurement(id);
      setPicking(false);
      setSelectedAtoms([]);
      setToast('Measurement added. Click the plot to explore it.');
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setMeasureBusy(false);
    }
  }
  function focusSelection() {
    const q = focusQuery.trim().toLowerCase();
    if (!dataset || !q) return;
    const matches = dataset.atoms.filter(
      (a) =>
        String(a.resid) === q ||
        `${a.residue}${a.resid}`.toLowerCase() === q ||
        a.residue.toLowerCase() === q ||
        `${a.chain}:${a.resid}`.toLowerCase() === q,
    );
    if (matches.length) {
      viewer.current?.focus(matches.map((a) => a.index));
      setToast(`Focused on ${q.toUpperCase()}`);
    } else setToast('No matching residue. Try a residue number, GLY35, or A:35.');
  }
  async function snapshot() {
    try {
      const blob = await viewer.current?.snapshot();
      if (!blob) throw new Error('Snapshot unavailable until the scene is ready.');
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = `DynaMol-${dataset?.name ?? 'scene'}-frame-${Math.round(frame) + 1}.png`;
      a.click();
      setTimeout(() => URL.revokeObjectURL(a.href), 1000);
      setToast('Snapshot saved.');
    } catch (e) {
      setError((e as Error).message);
    }
  }
  async function openById(
    id: string,
    options: { keepStudio?: boolean; showWater?: boolean; showHydrogens?: boolean } = {},
  ) {
    const requestToken = ++loadToken.current;
    setPlaying(false);
    try {
      const next = await api.dataset(id);
      if (requestToken === loadToken.current) await loadDataset(next, options);
    } catch (e) {
      if (requestToken === loadToken.current) setError((e as Error).message);
    }
  }
  const counts =
    dataset?.atoms.reduce<Record<string, number>>((acc, a) => {
      acc[a.category] = (acc[a.category] ?? 0) + 1;
      return acc;
    }, {}) ?? {};
  const searchResults =
    dataset && atomSearch.trim()
      ? dataset.atoms
          .filter((a) =>
            `${atomLabel(a)} ${a.chain} ${a.index + 1}`
              .toLowerCase()
              .includes(atomSearch.toLowerCase()),
          )
          .slice(0, 30)
      : [];
  const currentConfig = measureConfig[kind],
    running = jobs.filter(isActive),
    lastTime = dataset?.times_ps[dataset.n_frames - 1] ?? 0;
  const handleReady = useCallback(() => setReady(true), []),
    handleViewerError = useCallback((message: string) => {
      setError(message);
      setReady(false);
    }, []);
  function openStudio() {
    setModal('simulation');
  }

  return (
    <div className={`app-shell ${modal === 'simulation' ? 'studio-open' : ''}`}>
      <header className="app-header">
        <a
          className="brand"
          href="#"
          onClick={(e) => {
            e.preventDefault();
            viewer.current?.fit();
          }}
          aria-label="DynaMol home"
        >
          <span className="brand-mark">
            <Atom size={27} strokeWidth={1.8} />
          </span>
          <span>
            Dyna<span>Mol</span>
            <small>MOLECULES IN MOTION</small>
          </span>
          <b>α</b>
        </a>
        <nav className="main-nav" aria-label="Workspace">
          <button className={modal !== 'simulation' ? 'active' : ''} onClick={() => setModal(null)}>
            <Box size={16} /> Explore
          </button>
          <button className={modal === 'simulation' ? 'active' : ''} onClick={openStudio}>
            <FlaskConical size={16} /> Simulate
            {running.length > 0 && <span className="nav-count">{running.length}</span>}
          </button>
        </nav>
        <div className="header-right">
          <span className="local-status">
            <i className={health ? 'online' : ''} />
            {health ? 'Local workspace' : 'Connecting…'}
          </span>
          <button
            className="icon-button"
            title="Quick guide and keyboard shortcuts"
            aria-label="Open quick guide"
            onClick={() => setModal('help')}
          >
            <CircleHelp size={18} />
          </button>
          <span className="version">v0.1</span>
        </div>
      </header>
      <div className="workspace-bar">
        <div className="breadcrumb">
          <span>Workspace</span>
          <ChevronRight size={13} />
          <button
            onClick={() => setLibrary(!library)}
            className="dataset-switch"
            title="Switch dataset"
          >
            {dataset?.name ?? 'Your molecules'}
            <ChevronDown size={14} />
          </button>
          {dataset && (
            <span className="tag">{dataset.n_frames > 1 ? 'TRAJECTORY' : 'STRUCTURE'}</span>
          )}
        </div>
        <div className="workspace-actions">
          <button className="secondary-button" onClick={() => setModal('import')}>
            <FolderOpen size={15} /> Open files
          </button>
          <button className="primary-button" onClick={openStudio}>
            <Plus size={16} /> New simulation
          </button>
        </div>
        {library && (
          <div className="library-popover">
            <span className="section-label">IN THIS WORKSPACE</span>
            {datasets.length ? (
              datasets.map((d) => (
                <button key={d.id} onClick={() => void openById(d.id)}>
                  <Box size={15} />
                  <span>
                    <b>{d.name}</b>
                    <small>
                      {d.n_atoms.toLocaleString()} atoms · {d.n_frames} frames
                    </small>
                  </span>
                  {dataset?.id === d.id && <Check size={14} />}
                </button>
              ))
            ) : (
              <p>No datasets yet.</p>
            )}
          </div>
        )}
      </div>
      {error && (
        <div className="global-error" role="alert">
          <Info size={16} />
          <span>{error}</span>
          <button onClick={() => setError('')} aria-label="Dismiss error">
            <X size={16} />
          </button>
        </div>
      )}
      <main
        className={`workspace-grid ${inspector && modal !== 'simulation' ? '' : 'inspector-hidden'}`}
        ref={workspace}
      >
        <aside className="scene-panel">
          <div className="panel-heading">
            <span className="section-label">
              <Layers3 size={14} /> SCENE
            </span>
            <span className="small-tag">01</span>
          </div>
          <div className="structure-card">
            <div className="structure-card-top">
              <div className="structure-icon">
                <Atom size={24} />
              </div>
              <span className="dot online" />
            </div>
            <h2>{dataset?.name ?? 'A space for discovery'}</h2>
            <p>
              {dataset
                ? `${dataset.n_residues} residues · ${dataset.n_atoms.toLocaleString()} atoms`
                : 'Open a structure to get started'}
            </p>
            <button className="text-button" onClick={() => setDetails(!details)}>
              <Info size={12} /> Structure details{' '}
              <ChevronRight size={12} className={details ? 'rotate-90' : ''} />
            </button>
            {details && (
              <div className="dataset-details">
                <p>{dataset?.description}</p>
                <p>Source: {dataset?.source}</p>
                <p>
                  {dataset?.n_frames} saved frames ·{' '}
                  {dataset?.has_unitcell ? 'periodic cell' : 'nonperiodic coordinates'}
                </p>
                {dataset?.warnings.map((w, i) => (
                  <p className="inline-warning" key={i}>
                    {w}
                  </p>
                ))}
              </div>
            )}
          </div>
          <div className="scene-section">
            <div className="section-title">
              <h3>Representation</h3>
              <Settings2 size={13} />
            </div>
            <div className="representation-grid">
              {(
                [
                  { id: 'cartoon', label: 'Ribbons', icon: <RibbonIcon /> },
                  { id: 'ball+stick', label: 'Ball & stick', icon: <Atom size={21} /> },
                  { id: 'licorice', label: 'Sticks', icon: <StickIcon /> },
                  { id: 'surface', label: 'Surface', icon: <Box size={21} /> },
                ] as const
              ).map((rep) => (
                <button
                  key={rep.id}
                  className={representation === rep.id ? 'selected' : ''}
                  onClick={() => setRepresentation(rep.id)}
                >
                  {rep.icon}
                  <span>{rep.label}</span>
                </button>
              ))}
            </div>
            <label className="color-select">
              <span>Color by</span>
              <select
                aria-label="Color molecules by"
                value={colorScheme}
                onChange={(e) => setColorScheme(e.target.value as typeof colorScheme)}
              >
                <option value="residue">Sequence</option>
                <option value="chain">Chain</option>
                <option value="element">Element</option>
              </select>
              <span className="color-dots">
                <i />
                <i />
                <i />
                <i />
              </span>
            </label>
          </div>
          <div className="scene-section">
            <div className="section-title">
              <h3>Visibility</h3>
              <Eye size={13} />
            </div>
            <div className="visibility-list">
              {(
                [
                  {
                    id: 'protein',
                    name: 'Protein & nucleic acids',
                    count: (counts.protein ?? 0) + (counts.nucleic ?? 0),
                    icon: <Activity size={16} />,
                    color: '#9da8f4',
                  },
                  {
                    id: 'ligands',
                    name: 'Ligands',
                    count: (counts.ligands ?? 0) + (counts.ligand ?? 0),
                    icon: <Atom size={16} />,
                    color: '#e5b578',
                  },
                  {
                    id: 'water',
                    name: 'Water',
                    count: counts.water ?? 0,
                    icon: <Droplets size={16} />,
                    color: '#80bfeb',
                  },
                  {
                    id: 'ions',
                    name: 'Ions',
                    count: (counts.ions ?? 0) + (counts.ion ?? 0),
                    icon: <Plus size={16} />,
                    color: '#daa3d6',
                  },
                ] as const
              ).map((group) => (
                <button
                  key={group.id}
                  className={`visibility-row ${visibility[group.id] ? 'visible' : ''}`}
                  onClick={() => setVisibility((v) => ({ ...v, [group.id]: !v[group.id] }))}
                  aria-pressed={visibility[group.id]}
                  aria-label={`${visibility[group.id] ? 'Hide' : 'Show'} ${group.name.toLowerCase()}`}
                >
                  <span className="group-icon" style={{ color: group.color }}>
                    {group.icon}
                  </span>
                  <span>
                    {group.name}
                    <small>{group.count.toLocaleString()} atoms</small>
                  </span>
                  {visibility[group.id] ? <Eye size={15} /> : <EyeOff size={15} />}
                </button>
              ))}
            </div>
            <div className="hydrogen-heading">
              <span className="hydrogen-symbol">H</span>
              <span>Hydrogens</span>
            </div>
            <div className="segmented small">
              {(
                [
                  { id: 'none', label: 'Hidden' },
                  { id: 'polar', label: 'Polar only' },
                  { id: 'all', label: 'All' },
                ] as const
              ).map((h) => (
                <button
                  key={h.id}
                  className={visibility.hydrogens === h.id ? 'active' : ''}
                  onClick={() => setVisibility((v) => ({ ...v, hydrogens: h.id }))}
                >
                  {h.label}
                </button>
              ))}
            </div>
          </div>
          <div className="scene-section focus-section">
            <div className="section-title">
              <h3>Go to a residue</h3>
              <Focus size={14} />
            </div>
            <form
              className="search-input"
              onSubmit={(e) => {
                e.preventDefault();
                focusSelection();
              }}
            >
              <Search size={14} />
              <input
                aria-label="Find a residue"
                placeholder="e.g. 35, GLY35, A:35"
                value={focusQuery}
                onChange={(e) => setFocusQuery(e.target.value)}
              />
              <button aria-label="Focus residue" title="Focus residue">
                <ArrowRight size={13} />
              </button>
            </form>
          </div>
          <div className="scene-bottom">
            <span>
              <i />
              PRIVATE BY DESIGN
            </span>
            <p>
              Your files. Your computer.
              <br />
              An open source workspace.
            </p>
            <button className="text-button" onClick={() => setModal('help')}>
              Meet DynaMol <MoveUpRight size={12} />
            </button>
          </div>
        </aside>
        <section className="viewer-column">
          <div className="molecule-viewport">
            <div className="viewport-top">
              <div>
                <span className="eyebrow">THE MOLECULAR CANVAS</span>
                <div className="canvas-title">
                  {dataset?.name.replace(/\s*[·–—]\s*in motion$/i, '') ?? 'Welcome to DynaMol'}{' '}
                  <span className="canvas-subtitle">
                    {dataset?.n_frames && dataset.n_frames > 1 ? 'in motion' : 'in focus'}
                  </span>
                </div>
              </div>
              <button
                className={`icon-button ${inspector ? 'selected' : ''}`}
                title="Toggle measurement inspector"
                aria-label="Toggle measurement inspector"
                onClick={() => setInspector(!inspector)}
              >
                <Settings2 size={17} />
              </button>
            </div>
            <div className="ngl-host">
              <MolecularViewer
                ref={viewer}
                dataset={dataset}
                coordinates={coordinates}
                frame={frame}
                visibility={visibility}
                representation={representation}
                colorScheme={colorScheme}
                selectedAtoms={selectedAtoms}
                measurements={visibleMeasurements}
                picking={picking}
                spin={spin}
                onAtomPick={onAtomPick}
                onReady={handleReady}
                onError={handleViewerError}
              />
            </div>
            {(loading || (!ready && dataset && !error)) && (
              <div className="canvas-loading">
                <div className="loading-orbit">
                  <Atom size={37} />
                </div>
                <strong>{loading ? loadingMessage : 'Building your molecular scene…'}</strong>
                <span>Preparing something worth exploring.</span>
              </div>
            )}
            {!dataset && !loading && (
              <div className="canvas-empty">
                <Atom size={54} />
                <h2>Every molecule has a story.</h2>
                <p>Open a structure or a trajectory to start exploring.</p>
                <button className="primary-button" onClick={() => setModal('import')}>
                  <FolderOpen size={16} /> Open molecular files
                </button>
                <button
                  className="text-button"
                  onClick={() => {
                    setLoading(true);
                    api
                      .demo()
                      .then(loadDataset)
                      .catch((e) => {
                        setError((e as Error).message);
                        setLoading(false);
                      });
                  }}
                >
                  Try the ubiquitin example <ArrowRight size={13} />
                </button>
              </div>
            )}
            <div className="canvas-tools">
              <button
                className={`icon-button ${visibleMeasurements.length ? 'selected' : ''}`}
                title={
                  visibleMeasurements.length
                    ? 'Hide measurement lines and labels (keep plots)'
                    : 'Show measurement lines and labels'
                }
                aria-label={visibleMeasurements.length ? 'Hide measurements' : 'Show measurements'}
                aria-pressed={visibleMeasurements.length > 0}
                disabled={!measurements.length}
                onClick={() =>
                  setMeasurements((ms) =>
                    ms.map((m) => ({ ...m, visible: visibleMeasurements.length === 0 })),
                  )
                }
              >
                {visibleMeasurements.length ? <Eye size={17} /> : <EyeOff size={17} />}
              </button>
              <span />
              <button
                className="icon-button"
                title="Zoom in"
                aria-label="Zoom in"
                onClick={() => viewer.current?.zoom(1.2)}
              >
                <Plus size={18} />
              </button>
              <button
                className="icon-button"
                title="Zoom out"
                aria-label="Zoom out"
                onClick={() => viewer.current?.zoom(0.8)}
              >
                <Minus size={18} />
              </button>
              <span />
              <button
                className="icon-button"
                title="Fit molecule (F)"
                aria-label="Fit molecule"
                onClick={() => viewer.current?.fit()}
              >
                <Focus size={18} />
              </button>
              <button
                className={`icon-button ${spin ? 'selected' : ''}`}
                title="Auto rotate"
                aria-label="Auto rotate"
                aria-pressed={spin}
                onClick={() => setSpin(!spin)}
              >
                <RotateCcw size={17} />
              </button>
              <span />
              <button
                className="icon-button"
                title="Save snapshot"
                aria-label="Save snapshot"
                onClick={() => void snapshot()}
                disabled={!ready}
              >
                <Camera size={17} />
              </button>
              <button
                className="icon-button"
                title="Full screen"
                aria-label="Full screen"
                onClick={() => {
                  if (document.fullscreenElement) void document.exitFullscreen();
                  else
                    void workspace.current
                      ?.requestFullscreen()
                      .catch(() => setToast('Fullscreen is unavailable in this browser.'));
                }}
              >
                <Expand size={16} />
              </button>
            </div>
            {picking && (
              <div className="picking-banner">
                <Crosshair size={14} />
                <span>
                  Pick {currentConfig.count} atoms · {selectedAtoms.length} selected
                </span>
                <button onClick={() => setPicking(false)} aria-label="Stop picking">
                  <X size={13} />
                </button>
              </div>
            )}
            <div className="canvas-footer">
              <span>
                <MousePointer2 size={12} /> Drag to rotate <i />
                Scroll to zoom <i />
                Right drag to pan
              </span>
              <span className="orientation">
                <b>Y</b>
                <i />
                <strong>X</strong>
                <em>Z</em>
              </span>
            </div>
          </div>
          <div className="timeline">
            <div className="timeline-top">
              <div className="playback-buttons">
                <button
                  className="icon-button"
                  aria-label="First frame"
                  title="First frame"
                  onClick={() => seek(0)}
                  disabled={!dataset}
                >
                  <SkipBack size={16} />
                </button>
                <button
                  className="icon-button"
                  aria-label="Previous frame"
                  title="Previous frame (←)"
                  onClick={() => seek(Math.round(frame) - 1)}
                  disabled={!dataset}
                >
                  <ArrowLeft size={15} />
                </button>
                <button
                  className="play-button"
                  onClick={togglePlayback}
                  disabled={!ready || (dataset?.n_frames ?? 0) < 2}
                  aria-label={playing ? 'Pause trajectory' : 'Play trajectory'}
                  title="Play / pause (Space)"
                >
                  {playing ? (
                    <Pause size={19} fill="currentColor" />
                  ) : (
                    <Play size={19} fill="currentColor" />
                  )}
                </button>
                <button
                  className="icon-button"
                  aria-label="Next frame"
                  title="Next frame (→)"
                  onClick={() => seek(Math.round(frame) + 1)}
                  disabled={!dataset}
                >
                  <ArrowRight size={15} />
                </button>
                <button
                  className="icon-button"
                  aria-label="Last frame"
                  title="Last frame"
                  onClick={() => seek((dataset?.n_frames ?? 1) - 1)}
                  disabled={!dataset}
                >
                  <SkipForward size={16} />
                </button>
              </div>
              <div className="time-display">
                <b>{time.toFixed(2)}</b>
                <span>
                  {dataset?.time_unit === 'frame' ? 'frame index' : 'ps'} <i>/</i>{' '}
                  {lastTime.toFixed(2)} {dataset?.time_unit === 'frame' ? '' : 'ps'}
                </span>
              </div>
              <div className="playback-options">
                <select
                  aria-label="Playback speed"
                  value={speed}
                  onChange={(e) => setSpeed(Number(e.target.value))}
                >
                  <option value="0.25">0.25×</option>
                  <option value="0.5">0.5×</option>
                  <option value="1">1×</option>
                  <option value="2">2×</option>
                  <option value="4">4×</option>
                </select>
                <button
                  className={`icon-button ${loop ? 'selected' : ''}`}
                  aria-label="Loop playback"
                  aria-pressed={loop}
                  title="Loop playback"
                  onClick={() => setLoop(!loop)}
                >
                  <Repeat2 size={17} />
                </button>
              </div>
            </div>
            <input
              className="timeline-slider"
              type="range"
              aria-label="Trajectory frame"
              min="0"
              max={Math.max(1, (dataset?.n_frames ?? 1) - 1)}
              step="1"
              value={Math.round(frame)}
              disabled={!dataset || dataset.n_frames < 2}
              onChange={(e) => seek(Number(e.target.value))}
              style={
                {
                  '--progress': `${(frame / Math.max(1, (dataset?.n_frames ?? 1) - 1)) * 100}%`,
                } as React.CSSProperties
              }
            />
            <div className="timeline-meta">
              <span>
                FRAME <b>{dataset ? Math.round(frame) + 1 : 0}</b> / {dataset?.n_frames ?? 0}
              </span>
              <span>
                {dataset?.has_unitcell ? 'Saved-frame playback' : 'Smooth interpolation'}
                <i />
                {dataset?.n_frames ?? 0} saved frames
              </span>
              <span>
                {dataset?.time_unit === 'frame'
                  ? 'Frame index axis · time not supplied'
                  : 'Time in picoseconds'}
              </span>
            </div>
          </div>
          <PlotPanel
            dataset={dataset}
            measurements={measurements}
            frame={frame}
            activeId={activeMeasurement}
            onActive={setActiveMeasurement}
            onRemove={(id) => setMeasurements((ms) => ms.filter((m) => m.id !== id))}
            onToggleVisibility={(id) =>
              setMeasurements((ms) =>
                ms.map((m) => (m.id === id ? { ...m, visible: m.visible === false } : m)),
              )
            }
            onSeek={seek}
            onAdd={() => {
              setInspector(true);
              setPicking(true);
              setPlaying(false);
            }}
          />
          {dataset?.warnings.some((w) =>
            /time.*(unavailable|unknown|not|index)|frame indices|timestamps|physical time/i.test(w),
          ) && (
            <div className="timing-warning">
              <Info size={12} />
              {dataset.warnings
                .filter((w) =>
                  /time.*(unavailable|unknown|not|index)|frame indices|timestamps|physical time/i.test(
                    w,
                  ),
                )
                .join(' ')}
            </div>
          )}
        </section>
        {inspector && modal !== 'simulation' && (
          <aside className="measurement-panel">
            <div className="panel-heading">
              <span className="section-label">
                <Crosshair size={14} /> INSPECT & MEASURE
              </span>
              <button
                className="icon-button compact"
                onClick={() => setInspector(false)}
                aria-label="Close measurement inspector"
              >
                <X size={14} />
              </button>
            </div>
            <div className="measure-intro">
              <div className="measure-graphic">
                <span />
                <span />
                <span />
                <span />
                <i />
                <b>θ</b>
              </div>
              <h2>A closer look.</h2>
              <p>Turn a moment of molecular motion into a measurable insight.</p>
            </div>
            <div className="measure-section">
              <h3>What would you like to follow?</h3>
              <div className="measure-types">
                {(Object.keys(measureConfig) as MeasureKind[]).map((k) => (
                  <button
                    key={k}
                    className={kind === k ? 'selected' : ''}
                    onClick={() => {
                      setKind(k);
                      setSelectedAtoms([]);
                      if (k === 'hbond') setVisibility((v) => ({ ...v, hydrogens: 'polar' }));
                    }}
                  >
                    <span>{measureConfig[k].short}</span>
                    {measureConfig[k].name}
                  </button>
                ))}
              </div>
              <p className="measurement-hint">{currentConfig.hint}</p>
              <button
                className={`pick-button ${picking ? 'picking' : ''}`}
                onClick={() => {
                  setPicking(!picking);
                  setPlaying(false);
                  if (kind === 'hbond') setVisibility((v) => ({ ...v, hydrogens: 'polar' }));
                }}
                disabled={!dataset || !ready}
              >
                <MousePointer2 size={16} />
                {picking ? 'Picking atoms in the scene…' : 'Pick atoms in the scene'}
                {picking && <span className="pulse-dot" />}
              </button>
              <div className="selected-atoms">
                {Array.from({ length: currentConfig.count }, (_, i) => {
                  const idx = selectedAtoms[i],
                    a = idx === undefined ? null : dataset?.atoms[idx];
                  return (
                    <div key={i} className={`atom-slot ${a ? 'filled' : ''}`}>
                      <span>
                        {kind === 'hbond' ? ['D', 'H', 'A'][i] : String.fromCharCode(65 + i)}
                      </span>
                      {a ? (
                        <>
                          <div>
                            <b>
                              {a.residue} {a.resid}
                              <small>{a.name}</small>
                            </b>
                            <span>
                              Chain {a.chain || '—'} · atom {a.index + 1}
                            </span>
                          </div>
                          <button
                            className="icon-button compact"
                            onClick={() =>
                              setSelectedAtoms((prev) => prev.filter((x) => x !== idx))
                            }
                            aria-label={`Remove atom ${a.index + 1}`}
                          >
                            <X size={12} />
                          </button>
                        </>
                      ) : (
                        <p>
                          {kind === 'hbond'
                            ? [
                                'Select donor atom',
                                'Select bonded hydrogen',
                                'Select acceptor atom',
                              ][i]
                            : `Select atom ${i + 1}`}
                        </p>
                      )}
                    </div>
                  );
                })}
              </div>
              <div className="search-input atom-search">
                <Search size={13} />
                <input
                  value={atomSearch}
                  onChange={(e) => setAtomSearch(e.target.value)}
                  placeholder="Or find an atom by name…"
                  aria-label="Find an atom"
                />
                {atomSearch && (
                  <button aria-label="Clear atom search" onClick={() => setAtomSearch('')}>
                    <X size={12} />
                  </button>
                )}
              </div>
              {atomSearch && (
                <div className="atom-search-results">
                  {searchResults.length ? (
                    searchResults.map((a) => (
                      <button
                        key={a.index}
                        onClick={() => {
                          setSelectedAtoms((prev) =>
                            prev.includes(a.index)
                              ? prev
                              : prev.length >= currentConfig.count
                                ? [a.index]
                                : [...prev, a.index],
                          );
                          setAtomSearch('');
                          viewer.current?.focus([a.index]);
                        }}
                      >
                        <b>{atomLabel(a)}</b>
                        <span>
                          #{a.index + 1} · {a.element}
                        </span>
                      </button>
                    ))
                  ) : (
                    <span>No atoms found.</span>
                  )}
                </div>
              )}
              <button
                className="primary-button full-width plot-action"
                disabled={selectedAtoms.length !== currentConfig.count || measureBusy}
                onClick={() => void measure()}
              >
                {measureBusy ? <LoaderCircle size={15} className="spin" /> : <Activity size={15} />}{' '}
                Plot over time <ArrowRight size={14} />
              </button>
              {selectedAtoms.length > 0 && (
                <button
                  className="text-button center"
                  onClick={() => {
                    setSelectedAtoms([]);
                    setPicking(false);
                  }}
                >
                  Clear selection
                </button>
              )}
            </div>
            <div className="inspector-tip">
              <Sparkles size={16} />
              <div>
                <b>Find your focus</b>
                <p>
                  Pick a residue from the scene panel to fly straight to it. Press <kbd>F</kbd> to
                  bring the whole molecule back.
                </p>
              </div>
            </div>
            <div className="analysis-note">
              <Check size={12} />
              <span>
                Analysis uses saved coordinates.
                <br />
                {dataset?.has_unitcell
                  ? 'Periodic boundaries are respected.'
                  : 'Display smoothing does not alter results.'}
              </span>
            </div>
          </aside>
        )}
      </main>
      <footer className="status-bar">
        <span>
          <i className={ready ? 'online' : ''} />
          {loading
            ? 'Loading molecular data'
            : ready
              ? 'Ready to explore'
              : 'Waiting for a structure'}
          {dataset && (
            <>
              <em /> {dataset.n_atoms.toLocaleString()} atoms
            </>
          )}
        </span>
        <span>
          {running.length ? (
            <button onClick={openStudio}>
              <LoaderCircle size={11} className="spin" /> {running.length} background job
              {running.length > 1 ? 's' : ''} running · view progress <ChevronRight size={11} />
            </button>
          ) : (
            <>
              <Zap size={11} />
              {health?.engines
                .filter((e) => e.available)
                .map((e) => e.name)
                .join(' + ') || 'Checking engines'}
            </>
          )}
        </span>
        <button onClick={() => setModal('help')}>
          <Terminal size={11} /> Keyboard shortcuts <kbd>?</kbd>
        </button>
      </footer>
      {modal === 'import' && (
        <ImportDialog onClose={() => setModal(null)} onLoaded={(d) => void loadDataset(d)} />
      )}{' '}
      {modal === 'simulation' && (
        <SimulationPanel
          dataset={dataset}
          health={health}
          jobs={jobs}
          onClose={() => setModal(null)}
          onStarted={(j) => {
            setJobs((prev) => [j, ...prev.filter((existing) => existing.id !== j.id)]);
            setToast(
              ['preparation', 'solvation'].includes(j.engine)
                ? 'Preparing your structure. Progress appears in the studio.'
                : 'Simulation started. You can keep exploring.',
            );
          }}
          onLoad={(id, showWater) => void openById(id, { keepStudio: true, showWater })}
          onDatasetLoaded={async (d, options) => {
            await loadDataset(d, { keepStudio: true, ...options });
          }}
          onWaterVisibility={(show) =>
            setVisibility((v) => ({ ...v, water: show, ions: show || v.ions }))
          }
          onRefresh={refreshJobs}
        />
      )}
      {modal === 'help' && (
        <div
          className="modal-backdrop"
          onMouseDown={(e) => {
            if (e.target === e.currentTarget) setModal(null);
          }}
        >
          <section
            className="modal help-modal"
            role="dialog"
            aria-modal="true"
            aria-labelledby="help-title"
          >
            <div className="modal-header">
              <span className="eyebrow">MADE FOR CURIOSITY</span>
              <button
                className="icon-button"
                aria-label="Close quick guide"
                onClick={() => setModal(null)}
              >
                <X size={18} />
              </button>
            </div>
            <h2 id="help-title">Meet DynaMol.</h2>
            <p className="modal-subtitle">An open source home for molecular motion.</p>
            <div className="help-steps">
              <div>
                <span>1</span>
                <section>
                  <h3>Bring a molecule</h3>
                  <p>
                    Open a PDB, mmCIF, or GRO structure. Pair a matching topology with XTC, DCD,
                    TRR, NetCDF, or H5 trajectory data. The starter scene is a short, real OpenMM
                    simulation of ubiquitin.
                  </p>
                </section>
              </div>
              <div>
                <span>2</span>
                <section>
                  <h3>Explore the motion</h3>
                  <p>
                    Drag to rotate, scroll to zoom, and right-drag to pan. Change the representation
                    or hide solvent in the Scene pane. The timeline and plots stay linked.
                  </p>
                </section>
              </div>
              <div>
                <span>3</span>
                <section>
                  <h3>Follow an interaction</h3>
                  <p>
                    Choose a measurement, then pick atoms in order. Atom detail appears while
                    picking. Measurements use saved frames, with minimum-image periodic geometry
                    when cell data is available.
                  </p>
                </section>
              </div>
              <div>
                <span>4</span>
                <section>
                  <h3>Start something new</h3>
                  <p>
                    Choose OpenMM or GROMACS, adjust defaults, and start a local background run.
                    Follow the live log, cancel a job, or load and download its completed results in
                    Simulation Studio.
                  </p>
                </section>
              </div>
            </div>
            <div className="shortcut-grid">
              <span>
                Play / pause <kbd>Space</kbd>
              </span>
              <span>
                Previous / next frame <kbd>← →</kbd>
              </span>
              <span>
                Fit molecule <kbd>F</kbd>
              </span>
              <span>
                Pick atoms <kbd>M</kbd>
              </span>
              <span>
                Clear / close <kbd>Esc</kbd>
              </span>
              <span>
                This guide <kbd>?</kbd>
              </span>
            </div>
            <div className="help-limits">
              <b>Version 0.1 · an exploratory workbench</b>
              <p>
                Automatic preparation currently supports standard proteins. Arbitrary ligands,
                missing heavy atoms, membranes, and production equilibration protocols need
                additional preparation. Imports are capped for memory; use a stride for larger
                trajectories. Periodic datasets play saved frames to avoid smoothing across box
                boundaries. Short trajectories do not establish convergence.
              </p>
            </div>
            <div className="help-links">
              <a href="https://nglviewer.org/ngl/api/manual/" target="_blank" rel="noreferrer">
                <BookOpen size={13} /> NGL
              </a>
              <a href="https://docs.openmm.org/latest/userguide/" target="_blank" rel="noreferrer">
                OpenMM <MoveUpRight size={11} />
              </a>
              <a href="https://manual.gromacs.org/current/" target="_blank" rel="noreferrer">
                GROMACS <MoveUpRight size={11} />
              </a>
              <a href="https://mdtraj.org/" target="_blank" rel="noreferrer">
                MDTraj <MoveUpRight size={11} />
              </a>
            </div>
          </section>
        </div>
      )}
      {toast && (
        <div className="toast" role="status">
          <Check size={15} />
          {toast}
          <button aria-label="Dismiss notification" onClick={() => setToast('')}>
            <X size={13} />
          </button>
        </div>
      )}
    </div>
  );
}
function RibbonIcon() {
  return (
    <svg width="23" height="23" viewBox="0 0 24 24" fill="none">
      <path
        d="M5 20C3 15 18 17 18 12S4 10 5 5M7 22C5 17 20 19 20 14S6 12 7 7"
        stroke="currentColor"
        strokeWidth="2"
      />
      <path d="M5 5L9 2L8 8Z" fill="currentColor" />
    </svg>
  );
}
function StickIcon() {
  return (
    <svg
      width="23"
      height="23"
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth="2"
      strokeLinecap="round"
    >
      <path d="M3 16L9 9L15 13L21 5M9 9L7 3M15 13L17 21" />
    </svg>
  );
}
