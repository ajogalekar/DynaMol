import { useEffect, useRef, useState } from 'react';
import {
  AlertTriangle,
  ArrowDownToLine,
  Check,
  ChevronDown,
  ChevronRight,
  FileUp,
  FlaskConical,
  Globe2,
  LoaderCircle,
  Search,
  Sparkles,
  Upload,
  WandSparkles,
} from 'lucide-react';
import type { Dataset, Inspection, Job } from '../types';
import { api } from '../api';

export default function StructureWorkbench({
  dataset,
  ph,
  onPh,
  seed,
  locked,
  preparing,
  onDatasetLoaded,
  onPreparationStarted,
  onPreparationRequest,
}: {
  dataset: Dataset | null;
  ph: number;
  onPh: (ph: number) => void;
  seed: number;
  locked: boolean;
  preparing: boolean;
  onDatasetLoaded: (
    d: Dataset,
    options?: { showWater?: boolean; showHydrogens?: boolean },
  ) => Promise<void>;
  onPreparationStarted: (job: Job) => void;
  onPreparationRequest: (busy: boolean, error?: string) => void;
}) {
  const [tab, setTab] = useState<'upload' | 'fetch' | 'smiles'>('upload');
  const [provider, setProvider] = useState<'pdb' | 'pubchem'>('pdb');
  const [identifier, setIdentifier] = useState(''),
    [smiles, setSmiles] = useState(''),
    [molName, setMolName] = useState('');
  const [busy, setBusy] = useState(false),
    [error, setError] = useState(''),
    [inspection, setInspection] = useState<Inspection | null>(null),
    [inspecting, setInspecting] = useState(false),
    [inspectionError, setInspectionError] = useState('');
  const [submittingPreparation, setSubmittingPreparation] = useState(false);
  const [options, setOptions] = useState(false),
    [buildMissing, setBuildMissing] = useState(false),
    [addAtoms, setAddAtoms] = useState(true),
    [optimize, setOptimize] = useState(true),
    [removeWaters, setRemoveWaters] = useState(true),
    [removeHeterogens, setRemoveHeterogens] = useState(false);
  const input = useRef<HTMLInputElement>(null);
  const [ligandOverrides, setLigandOverrides] = useState<Record<string, string>>({});
  const [inspectionRevision, setInspectionRevision] = useState(0);
  useEffect(() => {
    setInspection(null);
    setBuildMissing(false);
    setRemoveHeterogens(false);
    setLigandOverrides({});
    setError('');
  }, [dataset?.id]);
  useEffect(() => {
    let current = true;
    setInspectionError('');
    if (!dataset) {
      setInspection(null);
      setInspecting(false);
      return;
    }
    setInspecting(true);
    const timer = window.setTimeout(() => {
      api
        .inspect(
          dataset.id,
          ph,
          Object.fromEntries(Object.entries(ligandOverrides).filter(([, value]) => value.trim())),
        )
        .then((r) => {
          if (current) setInspection(r);
        })
        .catch((e) => {
          if (current) setInspectionError((e as Error).message);
        })
        .finally(() => {
          if (current) setInspecting(false);
        });
    }, 350);
    return () => {
      current = false;
      window.clearTimeout(timer);
    };
  }, [dataset?.id, ph, ligandOverrides, inspectionRevision]);
  async function load(operation: () => Promise<Dataset>) {
    setBusy(true);
    setError('');
    try {
      const d = await operation();
      await onDatasetLoaded(d, { showWater: false });
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }
  async function upload(file: File) {
    const form = new FormData();
    form.append('file', file);
    await load(() => api.importStructure(form));
  }
  async function prepare() {
    if (!dataset || locked || busy) return;
    setBusy(true);
    setSubmittingPreparation(true);
    onPreparationRequest(true);
    setError('');
    let requestError: string | undefined;
    try {
      const job = await api.prepare({
        dataset_id: dataset.id,
        name: `${dataset.name.slice(0, 68)} · prepared`,
        ph,
        add_missing_atoms: addAtoms,
        build_missing_residues: buildMissing,
        optimize_sidechains: optimize,
        remove_waters: removeWaters,
        remove_heterogens: removeHeterogens,
        ligand_overrides: Object.fromEntries(
          Object.entries(ligandOverrides).filter(([, value]) => value.trim()),
        ),
        seed,
      });
      onPreparationStarted(job);
    } catch (e) {
      requestError = (e as Error).message;
    } finally {
      setBusy(false);
      setSubmittingPreparation(false);
      onPreparationRequest(false, requestError);
    }
  }
  const missingCount = inspection?.missing_residues.reduce((n, r) => n + r.count, 0) ?? 0;
  const canBuild = inspection?.missing_residues.some((r) => r.buildable) ?? false;
  const hasProtein =
    (inspection?.protein_atoms ??
      dataset?.atoms.filter((a) => a.category === 'protein').length ??
      0) > 0;
  const warningCount = missingCount + (inspection?.gaps.length ?? 0);
  const disabled = locked || busy;
  const complex = hasProtein && !!inspection?.ligands?.length && !removeHeterogens;
  return (
    <section className="structure-workbench" aria-label="Structure loading and protein preparation">
      <div className="field-heading source-heading">
        <span>
          <Upload size={13} />
        </span>
        <h3>Start with a structure</h3>
        <small>Loads in the live view</small>
      </div>
      <div className="source-tabs" role="tablist" aria-label="Structure source">
        <button
          type="button"
          role="tab"
          aria-selected={tab === 'upload'}
          onClick={() => setTab('upload')}
          className={tab === 'upload' ? 'active' : ''}
        >
          <FileUp size={13} /> Upload
        </button>
        <button
          type="button"
          role="tab"
          aria-selected={tab === 'fetch'}
          onClick={() => setTab('fetch')}
          className={tab === 'fetch' ? 'active' : ''}
        >
          <Globe2 size={13} /> Fetch
        </button>
        <button
          type="button"
          role="tab"
          aria-selected={tab === 'smiles'}
          onClick={() => setTab('smiles')}
          className={tab === 'smiles' ? 'active' : ''}
        >
          <FlaskConical size={13} /> SMILES
        </button>
      </div>
      {tab === 'upload' && (
        <button
          type="button"
          className="studio-upload"
          onClick={() => input.current?.click()}
          disabled={disabled}
          onDragOver={(e) => e.preventDefault()}
          onDrop={(e) => {
            e.preventDefault();
            if (!disabled && e.dataTransfer.files[0]) void upload(e.dataTransfer.files[0]);
          }}
        >
          <Upload size={21} />
          <span>
            <b>Drop a structure, or choose a file</b>
            <small>PDB · mmCIF · MOL2 · SDF · SMILES</small>
          </span>
          <ChevronRight size={14} />
        </button>
      )}
      <input
        hidden
        ref={input}
        type="file"
        aria-label="Upload simulation structure"
        accept=".pdb,.cif,.mmcif,.pdbx,.mol2,.sdf,.mol,.smi,.smiles"
        onChange={(e) => {
          const f = e.target.files?.[0];
          e.target.value = '';
          if (f) void upload(f);
        }}
      />
      {tab === 'fetch' && (
        <div className="studio-source-form">
          <label>
            Database
            <select
              aria-label="Structure database"
              value={provider}
              onChange={(e) => setProvider(e.target.value as typeof provider)}
              disabled={disabled}
            >
              <option value="pdb">RCSB Protein Data Bank</option>
              <option value="pubchem">PubChem small molecules</option>
            </select>
          </label>
          <div className="source-fetch-row">
            <input
              aria-label={provider === 'pdb' ? 'PDB accession' : 'PubChem name or CID'}
              value={identifier}
              onChange={(e) => setIdentifier(e.target.value)}
              placeholder={provider === 'pdb' ? 'PDB ID, e.g. 1UBQ' : 'Name or CID, e.g. caffeine'}
              disabled={disabled}
            />
            <button
              type="button"
              className="secondary-button"
              disabled={disabled || !identifier.trim()}
              onClick={() => void load(() => api.fetchStructure(provider, identifier.trim()))}
            >
              <Search size={13} /> Fetch
            </button>
          </div>
          <p className="form-note">
            {provider === 'pdb'
              ? 'Fetches the deposited structure and sequence records used to identify missing residues.'
              : 'Fetches a molecular structure from PubChem. Small molecules can be viewed; protein MD presets do not parameterize them.'}
          </p>
        </div>
      )}
      {tab === 'smiles' && (
        <div className="studio-source-form">
          <label>
            Molecule name <small>optional</small>
            <input
              aria-label="SMILES molecule name"
              placeholder="My molecule"
              value={molName}
              onChange={(e) => setMolName(e.target.value)}
              disabled={disabled}
            />
          </label>
          <label>
            SMILES
            <textarea
              aria-label="SMILES string"
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              placeholder="e.g. CCO"
              disabled={disabled}
              rows={2}
            />
          </label>
          <button
            type="button"
            className="secondary-button full-width"
            disabled={disabled || !smiles.trim()}
            onClick={() => void load(() => api.smiles(smiles.trim(), molName.trim() || undefined))}
          >
            <Sparkles size={13} /> Build 3D structure & view
          </button>
          <p className="form-note">
            Generates a local 3D conformer with RDKit. This does not assign MD force-field
            parameters.
          </p>
        </div>
      )}
      {busy && !submittingPreparation && (
        <div className="studio-working" role="status">
          <LoaderCircle size={13} className="spin" /> Preparing your request…
        </div>
      )}
      {error && (
        <div className="error-box" role="alert">
          {error}
        </div>
      )}
      {dataset && (
        <>
          <div className="source-summary">
            <div className="source-icon">
              <FlaskConical size={22} />
            </div>
            <div>
              <span>
                In view ·{' '}
                {dataset.n_frames > 1 ? 'first saved frame for preparation' : 'starting structure'}
              </span>
              <strong>{dataset.name}</strong>
              <small>
                {dataset.n_atoms.toLocaleString()} atoms · {dataset.n_residues} residues
              </small>
            </div>
            <Check size={15} />
          </div>
          <div className="prep-section">
            <div className="prep-heading">
              <div>
                <WandSparkles size={15} />
                <h3>{complex ? 'Protein + ligand preparation' : 'Protein preparation'}</h3>
              </div>
              <span>{dataset.preparation ? 'Prepared structure' : 'Review & repair'}</span>
            </div>
            {inspecting ? (
              <div className="studio-working">
                <LoaderCircle size={12} className="spin" /> Checking residues and ligand chemistry…
              </div>
            ) : inspectionError ? (
              <div className="inline-warning">
                {inspectionError}{' '}
                <button
                  type="button"
                  className="text-button"
                  onClick={() => setInspectionRevision((n) => n + 1)}
                >
                  Retry inspection
                </button>
              </div>
            ) : (
              inspection && (
                <>
                  <div className={`inspection-summary ${warningCount ? 'has-issues' : ''}`}>
                    {warningCount ? <AlertTriangle size={14} /> : <Check size={14} />}
                    <span>
                      {!hasProtein
                        ? 'Small-molecule structure'
                        : warningCount
                          ? `${missingCount} missing residues · ${inspection.gaps.length} chain gaps`
                          : `${inspection.missing_atoms.length} residues with missing atoms · ${inspection.hydrogen_atoms.toLocaleString()} hydrogens`}
                    </span>
                  </div>
                  {!!inspection.missing_residues.length && (
                    <div className="missing-residues">
                      <b>Missing sequence regions</b>
                      {inspection.missing_residues.map((r, i) => (
                        <div key={i}>
                          <span>
                            Chain {r.chain || '—'} · {r.count} residues{' '}
                            {r.terminal ? 'at terminus' : ''}
                          </span>
                          <small>
                            {r.residues.join('–')} ·{' '}
                            {r.buildable
                              ? 'can build a starting model'
                              : r.terminal ? 'reported · prepare the observed terminus' : 'requires additional modeling'}
                          </small>
                        </div>
                      ))}
                    </div>
                  )}
                  {!!inspection.gaps.length && (
                    <div className="inline-warning">
                      {inspection.gaps.map((g) => g.message).join(' ')}
                    </div>
                  )}
                  {!inspection.has_sequence && hasProtein && (
                    <p className="form-note prep-sequence-note">
                      Full sequence records are absent. Unknown loop identities cannot be inferred;
                      fetch the deposited PDB or upload a file with its sequence records.
                    </p>
                  )}
                  {!hasProtein && (
                    <p className="form-note">
                      Protein preparation is for standard amino-acid proteins. This molecule is
                      available in the viewer; ligand parameterization is a separate workflow.
                    </p>
                  )}
                </>
              )
            )}
            {!!inspection?.ligands?.length && (
              <div className="ligand-preparation-list" aria-label="Ligand preparation choices">
                <div className="ligand-preparation-intro">
                  <FlaskConical size={15} />
                  <div>
                    <b>
                      {inspection.ligands.length} ligands ·{' '}
                      {removeHeterogens ? 'removal selected' : 'kept in the complex'}
                    </b>
                    <span>
                      Bound poses stay in place. Review the fixed charge states before preparing.
                    </span>
                  </div>
                </div>
                {inspection.ligands.map((ligand) => (
                  <details
                    className={`ligand-preparation-card${ligand.error ? ' has-issues' : ''}`}
                    key={ligand.key}
                  >
                    <summary>
                      <span>
                        <b>{ligand.component_id}</b> · {ligand.chain}:{ligand.resid}
                      </span>
                      <strong>
                        {ligand.error
                          ? 'Needs chemistry'
                          : `${(ligand.formal_charge ?? 0) > 0 ? '+' : ''}${ligand.formal_charge ?? 0} charge`}
                      </strong>
                    </summary>
                    {ligand.error && <p className="inline-warning">{ligand.error}</p>}
                    <p className="form-note">{ligand.protonation_method}</p>
                    {ligand.warnings?.map((warning, i) => (
                      <p className="form-note" key={i}>
                        {warning}
                      </p>
                    ))}
                    {ligand.selected_smiles && (
                      <p className="ligand-smiles" aria-label={`Selected SMILES ${ligand.key}`}>
                        {ligand.selected_smiles}
                      </p>
                    )}
                    <label className="full-label">
                      Explicit-state SMILES · optional
                      <textarea
                        aria-label={`Ligand state override ${ligand.key}`}
                        value={ligandOverrides[ligand.key] ?? ''}
                        onChange={(e) =>
                          setLigandOverrides((values) => ({
                            ...values,
                            [ligand.key]: e.target.value,
                          }))
                        }
                        placeholder="Provide the same molecule with your chosen charge and stereochemistry"
                        disabled={disabled || removeHeterogens}
                        rows={3}
                        maxLength={10000}
                      />
                    </label>
                  </details>
                ))}
                {!!inspection.metal_environment?.retained_coordinating_waters.length && (
                  <p className="form-note">
                    {inspection.metal_environment.retained_coordinating_waters.length}{' '}
                    metal-coordinating waters will be retained with the ions. Nearby protein donor
                    sidechains are protected.
                  </p>
                )}
                {!inspection.ligand_runtime?.available && (
                  <div className="inline-warning">
                    {inspection.ligand_runtime?.message ??
                      'Ligand preparation tools are unavailable.'}
                  </div>
                )}
              </div>
            )}
            {!!inspection?.blockers.length && (
              <div className="inline-warning">{inspection.blockers.join(' ')}</div>
            )}
            <div className="prep-action-row">
              <label className="prep-ph">
                Target pH
                <input
                  aria-label="Preparation pH"
                  type="number"
                  min="0"
                  max="14"
                  step="0.1"
                  value={ph}
                  onChange={(e) => onPh(Number(e.target.value))}
                  disabled={disabled}
                />
              </label>
              <button
                type="button"
                className="primary-button"
                disabled={disabled || !hasProtein || inspecting}
                aria-busy={submittingPreparation || preparing}
                onClick={() => void prepare()}
              >
                {submittingPreparation || preparing ? (
                  <LoaderCircle size={14} className="spin" />
                ) : (
                  <WandSparkles size={14} />
                )}
                {submittingPreparation
                  ? 'Starting preparation…'
                  : preparing
                    ? complex
                      ? 'Preparing complex…'
                      : 'Preparing protein…'
                    : complex
                      ? 'Prep complex'
                      : 'Prep protein'}
              </button>
            </div>
            <p className="form-note">
              Rebuild protein hydrogens for the selected pH, repair missing heavy atoms, and sample
              side-chain clashes.
              {complex
                ? ' Ligands receive GAFF2 / AM1-BCC parameters for explicit-water OpenMM runs.'
                : ' Protein relaxation uses backbone restraints.'}{' '}
              Ionization uses recorded defaults, not a site-specific pKₐ calculation.
            </p>
            <button
              type="button"
              className="prep-options-toggle"
              aria-expanded={options}
              onClick={() => setOptions(!options)}
            >
              {options ? <ChevronDown size={13} /> : <ChevronRight size={13} />} Preparation options
            </button>
            {(options || canBuild) && (
              <div className="prep-options">
                <label className="checkbox-label">
                  <input
                    type="checkbox"
                    checked={addAtoms}
                    onChange={(e) => setAddAtoms(e.target.checked)}
                    disabled={disabled}
                  />{' '}
                  Add missing heavy atoms
                </label>
                <label className="checkbox-label">
                  <input
                    type="checkbox"
                    checked={optimize}
                    onChange={(e) => setOptimize(e.target.checked)}
                    disabled={disabled}
                  />{' '}
                  Refine side-chain rotamers and clashes
                </label>
                <label className="checkbox-label">
                  <input
                    type="checkbox"
                    checked={removeWaters}
                    onChange={(e) => setRemoveWaters(e.target.checked)}
                    disabled={disabled}
                  />{' '}
                  Remove existing waters before preparation (keep metal-coordinating waters)
                </label>
                <label className="checkbox-label">
                  <input
                    type="checkbox"
                    checked={removeHeterogens}
                    onChange={(e) => setRemoveHeterogens(e.target.checked)}
                    disabled={disabled}
                  />{' '}
                  Remove ligands and other non-protein residues
                </label>
                <label className={`checkbox-label loop-option ${canBuild ? 'available' : ''}`}>
                  <input
                    type="checkbox"
                    checked={buildMissing}
                    onChange={(e) => setBuildMissing(e.target.checked)}
                    disabled={disabled || !canBuild}
                  />{' '}
                  Build supported missing loops / residues
                </label>
                <p className="form-note">
                  {canBuild
                    ? 'Optional PDBFixer starting models for short internal gaps. Rebuilt loops are uncertain and need inspection; long gaps and terminal regions require additional modeling.'
                    : 'Missing-region construction needs known sequence identities and a supported short internal gap.'}
                </p>
              </div>
            )}
            {dataset.preparation && (
              <div className="preparation-result">
                <b>
                  <Check size={13} /> Preparation recorded · pH {dataset.preparation.ph}
                </b>
                {dataset.preparation.summary?.map((line, i) => (
                  <p key={i}>{line}</p>
                ))}
                {dataset.preparation.warnings?.map((line, i) => (
                  <p className="prep-result-warning" key={i}>
                    {line}
                  </p>
                ))}
                <a
                  href={`/api/datasets/${dataset.id}/prepared`}
                  download={`${dataset.name}.pdb`}
                  className="text-button"
                >
                  <ArrowDownToLine size={12} /> Prepared PDB
                </a>
                {dataset.preparation.job_id && dataset.preparation.ligand_parameters && (
                  <a
                    className="text-button"
                    href={`/api/jobs/${dataset.preparation.job_id}/download`}
                  >
                    <ArrowDownToLine size={12} /> Parameters & preparation files
                  </a>
                )}
              </div>
            )}
          </div>
        </>
      )}
    </section>
  );
}
