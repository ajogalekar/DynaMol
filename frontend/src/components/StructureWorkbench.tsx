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
  Layers3,
  Wrench,
  Search,
  Sparkles,
  Upload,
  WandSparkles,
} from 'lucide-react';
import type { Dataset, Inspection, Job } from '../types';
import { api } from '../api';
import './ModifiedResidues.css';
import ReadinessPanel from './ReadinessPanel';

export default function StructureWorkbench({
  dataset,
  engine,
  ph,
  onPh,
  seed,
  locked,
  preparing,
  onDatasetLoaded,
  onPreparationStarted,
  onPreparationRequest,
  onSourceRequest,
}: {
  dataset: Dataset | null;
  engine: 'openmm' | 'gromacs';
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
  onSourceRequest?: (busy: boolean) => void;
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
  const [monomers, setMonomers] = useState<Awaited<ReturnType<typeof api.monomers>> | null>(null);
  const [monomerChain, setMonomerChain] = useState<number | null>(null);
  const [monomerLoading, setMonomerLoading] = useState(false);
  const [monomerError, setMonomerError] = useState('');
  const [monomerRevision, setMonomerRevision] = useState(0);
  const [selectingMonomer, setSelectingMonomer] = useState(false);
  const datasetIdRef = useRef(dataset?.id);
  datasetIdRef.current = dataset?.id;
  const operationRevision = useRef(0);
  const sourceRequest = useRef(onSourceRequest);
  sourceRequest.current = onSourceRequest;
  const mounted = useRef(false);
  useEffect(() => {
    mounted.current = true;
    return () => {
      mounted.current = false;
      operationRevision.current++;
      sourceRequest.current?.(false);
    };
  }, []);
  const [prepReady, setPrepReady] = useState(false);
  const [options, setOptions] = useState(false),
    [buildMissing, setBuildMissing] = useState(false),
    [addAtoms, setAddAtoms] = useState(true),
    [optimize, setOptimize] = useState(true),
    [removeWaters, setRemoveWaters] = useState(true),
    [removeHeterogens, setRemoveHeterogens] = useState(false);
  const input = useRef<HTMLInputElement>(null);
  const [ligandOverrides, setLigandOverrides] = useState<Record<string, string>>({});
  const [ligandActions, setLigandActions] = useState<Record<string, 'repair' | 'remove'>>({});
  const [inspectionRevision, setInspectionRevision] = useState(0);
  useEffect(() => {
    setInspection(null);
    setBuildMissing(false);
    setRemoveHeterogens(false);
    setLigandOverrides({});
    setLigandActions({});
    setError('');
    setBusy(false);
    setSelectingMonomer(false);
    sourceRequest.current?.(false);
  }, [dataset?.id]);
  useEffect(() => {
    let current = true;
    setInspectionError('');
    if (!dataset || engine !== 'openmm') {
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
          ligandActions,
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
  }, [dataset?.id, engine, ph, ligandOverrides, ligandActions, inspectionRevision]);
  useEffect(() => {
    let current = true;
    setMonomers(null);
    setMonomerChain(null);
    setMonomerError('');
    if (!dataset || !dataset.atoms.some((atom) => atom.category === 'protein')) {
      setMonomerLoading(false);
      return;
    }
    setMonomerLoading(true);
    api
      .monomers(dataset.id)
      .then(
        (choices) => {
          if (!current) return;
          setMonomers(choices);
          setMonomerChain(choices.recommended_chain_index ?? choices.chains[0]?.index ?? null);
        },
        (reason) => {
          if (current) setMonomerError((reason as Error).message);
        },
      )
      .finally(() => {
        if (current) setMonomerLoading(false);
      });
    return () => {
      current = false;
    };
  }, [dataset?.id, monomerRevision]);
  async function load(operation: () => Promise<Dataset>, selectingChain = false) {
    const revision = ++operationRevision.current;
    const sourceId = datasetIdRef.current;
    const current = () =>
      mounted.current &&
      operationRevision.current === revision &&
      datasetIdRef.current === sourceId;
    setBusy(true);
    setSelectingMonomer(selectingChain);
    sourceRequest.current?.(true);
    setError('');
    try {
      const d = await operation();
      if (!current()) return;
      await onDatasetLoaded(d, { showWater: false });
    } catch (e) {
      if (current()) setError((e as Error).message);
    } finally {
      if (mounted.current && operationRevision.current === revision) {
        setBusy(false);
        setSelectingMonomer(false);
        sourceRequest.current?.(false);
      }
    }
  }
  async function upload(file: File) {
    const form = new FormData();
    form.append('file', file);
    await load(() => api.importStructure(form));
  }
  async function prepare() {
    if (!dataset || engine !== 'openmm' || locked || busy || !prepReady) return;
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
        ligand_actions: ligandActions,
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
  const loopPolicy = inspection?.loop_policy;
  const internalMissing = inspection?.missing_residues.filter((region) => !region.terminal) ?? [];
  const internalMissingCount = internalMissing.reduce((count, region) => count + region.count, 0);
  const longestMissingGap = internalMissing.reduce(
    (longest, region) => Math.max(longest, region.count),
    0,
  );
  const gapLimitExceeded = !!loopPolicy && longestMissingGap > loopPolicy.max_gap_residues;
  const totalLimitExceeded = !!loopPolicy && internalMissingCount > loopPolicy.max_total_residues;
  const isExtendedLoop = (region: Inspection['missing_residues'][number]) =>
    region.modeling === 'extended' ||
    (!region.modeling && !!loopPolicy && region.count > loopPolicy.short_gap_residues);
  const hasExtendedBuildableLoop = internalMissing.some(
    (region) => region.buildable && isExtendedLoop(region),
  );
  const hasProtein =
    (inspection?.protein_atoms ??
      dataset?.atoms.filter((a) => a.category === 'protein').length ??
      0) > 0;
  const structuralGaps = inspection?.gaps.filter((gap) => gap.structural_break !== false) ?? [];
  const numberingGaps = inspection?.gaps.filter((gap) => gap.structural_break === false) ?? [];
  const warningCount = missingCount + structuralGaps.length;
  const disabled = locked || busy;
  const retainedLigands = inspection?.ligands?.filter((ligand) => !ligand.removed) ?? [];
  const repairableLigands =
    inspection?.ligands?.filter(
      (ligand) =>
        ligand.can_repair &&
        ligand.missing_heavy_atoms?.length &&
        ligandActions[ligand.key] !== 'remove',
    ) ?? [];
  const complex = hasProtein && !!retainedLigands.length && !removeHeterogens;
  return (
    <section className="structure-workbench" aria-label="Structure loading and protein preparation">
      <div className="field-heading source-heading">
        <span>02</span>
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
          <LoaderCircle size={13} className="spin" />{' '}
          {selectingMonomer
            ? 'Selecting one protein chain and associated molecules…'
            : 'Preparing your request…'}
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
          {monomerLoading ? (
            <div className="monomer-checking" role="status">
              <LoaderCircle size={12} className="spin" /> Checking protein chains…
            </div>
          ) : monomerError ? (
            <div className="inline-warning monomer-error">
              Could not check protein chains: {monomerError}{' '}
              <button
                type="button"
                className="text-button"
                disabled={disabled}
                onClick={() => setMonomerRevision((revision) => revision + 1)}
              >
                Retry chain check
              </button>
            </div>
          ) : (
            !!monomers &&
            monomers.chains.length > 1 && (
              <section className="monomer-setup" aria-labelledby="monomer-setup-title">
                <div className="monomer-heading">
                  <Layers3 size={16} />
                  <h4 id="monomer-setup-title">Choose one protein chain</h4>
                  <span>{monomers.chains.length} chains</span>
                </div>
                <p>
                  Use one chain as your starting structure, with its nearby ligands, cofactors and
                  ions.
                </p>
                <div className="monomer-actions">
                  <label>
                    Protein chain
                    <select
                      aria-label="Protein chain for monomer"
                      value={monomerChain ?? ''}
                      disabled={disabled}
                      onChange={(event) => setMonomerChain(Number(event.target.value))}
                    >
                      {monomers.chains.map((chain) => (
                        <option key={chain.index} value={chain.index}>
                          {chain.label} · {chain.n_residues} residues
                          {chain.index === monomers.recommended_chain_index ? ' · suggested' : ''}
                        </option>
                      ))}
                    </select>
                  </label>
                  <button
                    type="button"
                    className="secondary-button"
                    disabled={disabled || monomerChain === null}
                    onClick={() => {
                      if (monomerChain !== null)
                        void load(() => api.useMonomer(dataset.id, monomerChain, true), true);
                    }}
                  >
                    {selectingMonomer ? (
                      <LoaderCircle size={14} className="spin" />
                    ) : (
                      <Layers3 size={14} />
                    )}
                    {selectingMonomer ? 'Selecting monomer…' : 'Use one monomer'}
                  </button>
                </div>
                <p className="monomer-note">
                  Selects one protein chain; chain count alone does not establish the biological
                  assembly. Creates a new structure to prepare, with bulk solvent removed.
                </p>
                {!!monomers.warnings.length && (
                  <details className="monomer-notes">
                    <summary>Selection notes</summary>
                    {monomers.warnings.map((warning, index) => (
                      <p key={index}>{warning}</p>
                    ))}
                  </details>
                )}
              </section>
            )
          )}
          {dataset.monomer_selection && (
            <div className="monomer-result" role="status">
              <Check size={14} />
              <div>
                <b>{dataset.monomer_selection.chain_label} selected</b>
                <span>
                  {dataset.preparation
                    ? 'This selected chain is prepared for OpenMM.'
                    : 'One protein chain is loaded. Review missing regions and prepare this structure.'}
                </span>
                {!!(
                  dataset.monomer_selection.retained_associated_molecules.length ||
                  dataset.monomer_selection.excluded_molecules.length
                ) && (
                  <span>
                    {dataset.monomer_selection.retained_associated_molecules.length} associated
                    molecules retained · {dataset.monomer_selection.excluded_molecules.length} other
                    molecules omitted
                  </span>
                )}
                {!!dataset.monomer_selection.notes.length && (
                  <details className="monomer-notes">
                    <summary>Selection details</summary>
                    {dataset.monomer_selection.notes.map((note, index) => (
                      <p key={index}>{note}</p>
                    ))}
                  </details>
                )}
                <button
                  type="button"
                  className="text-button"
                  disabled={disabled}
                  onClick={() =>
                    void load(() => api.dataset(dataset.monomer_selection!.parent_dataset_id))
                  }
                >
                  Restore full structure
                </button>
              </div>
            </div>
          )}
          {engine === 'gromacs' ? (
            <div className="prep-section engine-preparation" aria-label="GROMACS preparation">
              <div className="prep-heading">
                <div>
                  <WandSparkles size={15} />
                  <h3>GROMACS preparation</h3>
                </div>
                <span>Native setup at Start</span>
              </div>
              <p className="form-note">
                Load a complete standard protein, review the simulation settings, then start.
                GROMACS generates its topology and rebuilds hydrogens using its default protonation
                and terminal states. Existing hydrogen choices are not retained.
              </p>
              <p className="form-note">
                This path does not offer pH selection, missing-loop repair or
                ligand/modified-residue parameterization. Choose OpenMM above for DynaMol’s
                supported protein and complex preparation tools.
              </p>
            </div>
          ) : (
            <div className="prep-section">
              <div className="prep-heading">
                <div>
                  <WandSparkles size={15} />
                  <h3>{complex ? 'Protein + ligand preparation' : 'Protein preparation'}</h3>
                </div>
                <span>
                  {dataset.preparation ? 'Prepared for OpenMM' : 'OpenMM · review & repair'}
                </span>
              </div>
              <section
                className="missing-structure-setup"
                aria-labelledby="missing-structure-title"
              >
                <div className="missing-structure-heading">
                  <Wrench size={14} />
                  <h4 id="missing-structure-title">Missing atoms &amp; residues</h4>
                </div>
                <div className="missing-structure-options">
                  <label className="checkbox-label">
                    <input
                      type="checkbox"
                      checked={addAtoms}
                      onChange={(event) => setAddAtoms(event.target.checked)}
                      disabled={disabled || !hasProtein}
                    />
                    Add missing heavy atoms
                  </label>
                  <label className={`checkbox-label loop-option ${canBuild ? 'available' : ''}`}>
                    <input
                      type="checkbox"
                      checked={buildMissing}
                      onChange={(event) => setBuildMissing(event.target.checked)}
                      disabled={disabled || !hasProtein}
                      aria-describedby="missing-loop-help"
                    />
                    Build supported missing loops / residues
                  </label>
                </div>
                {inspecting ? (
                  <div className="studio-working" role="status">
                    <LoaderCircle size={12} className="spin" /> Checking residues and ligand
                    chemistry…
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
                              ? `${missingCount} missing residues · ${structuralGaps.length} backbone breaks`
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
                                  ? isExtendedLoop(r)
                                    ? 'longer loop · provisional starting model'
                                    : 'can build a starting model'
                                  : r.terminal
                                    ? 'reported · prepare the observed terminus'
                                    : loopPolicy && r.count > loopPolicy.max_gap_residues
                                      ? `exceeds ${loopPolicy.max_gap_residues}-residue gap limit`
                                      : 'requires additional modeling'}
                              </small>
                            </div>
                          ))}
                        </div>
                      )}
                      {!!structuralGaps.length && (
                        <div className="inline-warning">
                          {structuralGaps.map((g) => g.message).join(' ')}
                        </div>
                      )}
                      {!!numberingGaps.length && (
                        <p className="form-note">
                          Residue numbering skips were found. Numbering skips alone are not missing
                          loops.
                        </p>
                      )}
                      {!inspection.has_sequence && hasProtein && (
                        <p className="form-note prep-sequence-note">
                          Full sequence records are absent. Unknown loop identities cannot be
                          inferred; fetch the deposited PDB or upload a file with its sequence
                          records.
                        </p>
                      )}
                      {!hasProtein && (
                        <p className="form-note">
                          Protein preparation is for amino-acid polymers. This molecule is available
                          in the viewer; ligand parameterization is a separate workflow.
                        </p>
                      )}
                    </>
                  )
                )}
                {!inspecting && loopPolicy && (gapLimitExceeded || totalLimitExceeded) && (
                  <div className="inline-warning missing-build-limit" role="status">
                    {gapLimitExceeded &&
                      `The longest missing internal region has ${longestMissingGap} residues; the computational limit is ${loopPolicy.max_gap_residues} per gap. `}
                    {totalLimitExceeded &&
                      `This structure has ${internalMissingCount} missing internal residues; the computational limit is ${loopPolicy.max_total_residues} total. `}
                    Supply a more complete model for regions beyond these limits.
                  </div>
                )}
                {!inspecting && loopPolicy && internalMissingCount > 0 && (
                  <p className="form-note">
                    Current computational limits: {loopPolicy.max_gap_residues} residues per
                    internal gap and {loopPolicy.max_total_residues} residues total.
                  </p>
                )}
                {!inspecting && inspection && hasProtein && (
                  <p className="form-note missing-build-note" id="missing-loop-help">
                    {canBuild
                      ? `Use the original sequence to build a provisional starting model. ${hasExtendedBuildableLoop ? 'Longer loops take more work. ' : ''}Geometry is checked; loop conformations remain uncertain.`
                      : missingCount
                        ? internalMissingCount
                          ? 'Unsupported internal gaps require additional modeling before preparation.'
                          : 'Terminal regions remain omitted; the observed termini are prepared.'
                        : inspection.has_sequence
                          ? 'No missing protein sequence regions were detected. This option has nothing to add for the current structure; ligand atom repair is separate below.'
                          : 'Loop building needs known sequence records and a supported internal gap.'}
                  </p>
                )}
              </section>
              {!!inspection?.modified_residues?.length && (
                <div className="modified-residue-list" aria-label="Modified protein residues">
                  <div className="modified-residue-intro">
                    <WandSparkles size={15} />
                    <div>
                      <b>{inspection.modified_residues.length} modified protein residues</b>
                      <span>Kept in the protein with their chemical identity.</span>
                    </div>
                  </div>
                  {inspection.modified_residues.map((residue, index) => (
                    <details
                      className={`modified-residue-card${residue.supported ? '' : ' has-issues'}`}
                      key={`${residue.chain}:${residue.resid}:${residue.insertion_code ?? ''}:${residue.residue}:${index}`}
                    >
                      <summary>
                        <span>
                          <b>{residue.residue}</b> · {residue.chain}:{residue.resid}
                          {residue.insertion_code ?? ''}
                        </span>
                        <strong>
                          {residue.supported ? 'Protein template' : 'Needs a template'}
                        </strong>
                      </summary>
                      {residue.parent_residue && (
                        <p>
                          {residue.name ? `${residue.name} · ` : ''}
                          Modified {residue.parent_residue} · original {residue.residue} identity
                          retained.
                        </p>
                      )}
                      {residue.forcefield && (
                        <p>
                          {residue.forcefield} · template {residue.template}
                        </p>
                      )}
                      {residue.formal_charge !== undefined && (
                        <p>
                          Fixed residue charge: {residue.formal_charge > 0 ? '+' : ''}
                          {residue.formal_charge}.
                        </p>
                      )}
                      {residue.protonation_note && <p>{residue.protonation_note}</p>}
                      {residue.error && <p className="inline-warning">{residue.error}</p>}
                      {!residue.supported && (
                        <p>
                          Provide parameters for this covalently connected amino acid and its chosen
                          chemical state to enable preparation.
                        </p>
                      )}
                    </details>
                  ))}
                  <p className="form-note">
                    Ligand removal preserves modified protein residues. Unsupported residues require
                    an exact amino-acid template.
                  </p>
                </div>
              )}
              {!!inspection?.ligands?.length && (
                <div className="ligand-preparation-list" aria-label="Ligand preparation choices">
                  <div className="ligand-preparation-intro">
                    <FlaskConical size={15} />
                    <div>
                      <b>
                        {inspection.ligands.length} molecules ·{' '}
                        {removeHeterogens
                          ? 'removal selected'
                          : `${retainedLigands.length} retained`}
                      </b>
                      <span>
                        Observed atoms stay in place. Choose how to handle incomplete molecules
                        before preparing.
                      </span>
                    </div>
                  </div>
                  {!!repairableLigands.some((ligand) => ligandActions[ligand.key] !== 'repair') && (
                    <button
                      type="button"
                      className="secondary-button repair-ligands-button"
                      disabled={disabled || removeHeterogens || inspecting}
                      onClick={() =>
                        setLigandActions((actions) => ({
                          ...actions,
                          ...Object.fromEntries(
                            repairableLigands.map((ligand) => [ligand.key, 'repair' as const]),
                          ),
                        }))
                      }
                    >
                      <Wrench size={14} /> Repair missing ligand atoms
                    </button>
                  )}
                  {inspection.ligands.map((ligand) => (
                    <details
                      className={`ligand-preparation-card${ligand.error ? ' has-issues' : ''}`}
                      key={ligand.key}
                      open={!!ligand.error || !!ligandActions[ligand.key]}
                    >
                      <summary>
                        <span>
                          <b>{ligand.component_id}</b> · {ligand.chain}:{ligand.resid}
                        </span>
                        <strong>
                          {ligand.removed
                            ? 'Removal selected'
                            : ligandActions[ligand.key] === 'repair' && !ligand.error
                              ? 'Repair selected'
                              : ligand.error
                                ? 'Needs chemistry'
                                : `${(ligand.formal_charge ?? 0) > 0 ? '+' : ''}${ligand.formal_charge ?? 0} charge`}
                        </strong>
                      </summary>
                      {!!ligand.missing_heavy_atoms?.length && (
                        <p className="ligand-missing-atoms">
                          <b>{ligand.missing_heavy_atoms.length} missing heavy atoms</b> ·{' '}
                          {ligand.missing_heavy_atoms.join(', ')}
                        </p>
                      )}
                      <label className="full-label ligand-action-label">
                        Preparation choice
                        <select
                          aria-label={`Preparation choice ${ligand.key}`}
                          value={ligandActions[ligand.key] ?? ''}
                          disabled={disabled || removeHeterogens}
                          onChange={(event) =>
                            setLigandActions((actions) => {
                              const next = { ...actions };
                              if (event.target.value)
                                next[ligand.key] = event.target.value as 'repair' | 'remove';
                              else delete next[ligand.key];
                              return next;
                            })
                          }
                        >
                          <option value="">Keep observed molecule</option>
                          <option value="repair" disabled={!ligand.can_repair}>
                            Build missing atoms from chemical reference
                          </option>
                          <option value="remove" disabled={ligand.can_remove === false}>
                            Remove this molecule from prepared structure
                          </option>
                        </select>
                      </label>
                      {ligandActions[ligand.key] === 'repair' && (
                        <p className="form-note">
                          Missing coordinates will be modeled during Prep while observed atoms stay
                          fixed. Inspect the rebuilt part afterward; its conformation is uncertain.
                        </p>
                      )}
                      {ligandActions[ligand.key] === 'remove' && (
                        <p className="form-note">
                          Only this molecule will be omitted from the prepared structure. The
                          original dataset remains available.
                        </p>
                      )}
                      {!ligand.can_repair && ligand.repair_reason && (
                        <p className="form-note">{ligand.repair_reason}</p>
                      )}
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
                          disabled={
                            disabled || removeHeterogens || ligandActions[ligand.key] === 'remove'
                          }
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
                  {!!retainedLigands.length && !inspection.ligand_runtime?.available && (
                    <div className="inline-warning">
                      {inspection.ligand_runtime?.message ??
                        'Ligand preparation tools are unavailable.'}
                    </div>
                  )}
                </div>
              )}
              {!!inspection?.ions?.some((ion) => ion.supported) && (
                <details className="ion-preparation-summary">
                  <summary>
                    Ions ·{' '}
                    {inspection.ions
                      .filter((ion) => ion.supported)
                      .map(
                        (ion) =>
                          `${ion.element}${ion.formal_charge > 0 ? '+' : ''}${ion.formal_charge}`,
                      )
                      .join(', ')}
                    {removeHeterogens ? ' · removal selected' : ' · retained'}
                  </summary>
                  {inspection.ions
                    .filter((ion) => ion.supported)
                    .map((ion) => (
                      <p className="form-note" key={ion.key}>
                        {ion.key} · {ion.model}
                      </p>
                    ))}
                </details>
              )}
              {!!inspection?.blockers.length && (
                <div className="inline-warning">{inspection.blockers.join(' ')}</div>
              )}
              <ReadinessPanel
                mode="preparation"
                datasetId={dataset.id}
                settings={{
                  ph,
                  add_missing_atoms: addAtoms,
                  build_missing_residues: buildMissing,
                  optimize_sidechains: optimize,
                  remove_waters: removeWaters,
                  remove_heterogens: removeHeterogens,
                  ligand_overrides: Object.fromEntries(
                    Object.entries(ligandOverrides).filter(([, value]) => value.trim()),
                  ),
                  ligand_actions: ligandActions,
                  seed,
                }}
                onReadyChange={setPrepReady}
                suspended={locked || preparing}
              />
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
                  disabled={disabled || !hasProtein || inspecting || !prepReady}
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
                Rebuild protein hydrogens for the selected pH, repair missing heavy atoms, and
                sample side-chain clashes.
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
                {options ? <ChevronDown size={13} /> : <ChevronRight size={13} />} Preparation
                options
              </button>
              {options && (
                <div className="prep-options">
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
                  {dataset.preparation.job_id &&
                    (dataset.preparation.ligand_parameters ||
                      !!dataset.preparation.modified_residues) && (
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
          )}
        </>
      )}
    </section>
  );
}
