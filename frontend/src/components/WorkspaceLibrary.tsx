import { useCallback, useEffect, useRef, useState } from 'react';
import {
  Archive,
  ArchiveRestore,
  Check,
  Download,
  FolderOpen,
  HardDrive,
  LoaderCircle,
  Pencil,
  Plus,
  Save,
  Search,
  Trash2,
  Upload,
  X,
} from 'lucide-react';
import type { Dataset } from '../types';
import {
  formatBytes,
  workspaceApi,
  type LibraryIndex,
  type SavedProject,
  type WorkspaceState,
} from '../workspace';
import './WorkspaceLibrary.css';

export default function WorkspaceLibrary({
  dataset,
  state,
  onClose,
  onOpenDataset,
  onOpenProject,
  onDatasetChanged,
  saveStatus,
  initialTab = 'projects',
}: {
  dataset: Dataset | null;
  state: () => WorkspaceState | null;
  onClose: () => void;
  onOpenDataset: (id: string) => Promise<void>;
  onOpenProject: (state: WorkspaceState) => Promise<void>;
  onDatasetChanged: (dataset?: Dataset) => void;
  saveStatus: string;
  initialTab?: 'projects' | 'library';
}) {
  const [tab, setTab] = useState<'projects' | 'library'>(initialTab);
  const [projects, setProjects] = useState<SavedProject[]>([]);
  const [index, setIndex] = useState<LibraryIndex | null>(null);
  const [query, setQuery] = useState('');
  const [filter, setFilter] = useState<'active' | 'archived' | 'trash'>('active');
  const [name, setName] = useState(dataset?.name ? `${dataset.name} study` : '');
  const [busy, setBusy] = useState('loading');
  const [error, setError] = useState('');
  const [notice, setNotice] = useState('');
  const [edit, setEdit] = useState<{
    kind: 'rename' | 'trash' | 'remove-project' | 'update-project';
    id: string;
    name: string;
  } | null>(null);
  const [editName, setEditName] = useState('');
  const fileInput = useRef<HTMLInputElement>(null);
  const refresh = useCallback(async () => {
    const [projects, library] = await Promise.all([
      workspaceApi.projects(),
      workspaceApi.library(),
    ]);
    setProjects(projects);
    setIndex(library);
  }, []);
  useEffect(() => {
    void refresh()
      .catch((e: Error) => setError(e.message))
      .finally(() => setBusy(''));
  }, [refresh]);
  async function act(key: string, work: () => Promise<void>) {
    setBusy(key);
    setError('');
    setNotice('');
    try {
      await work();
      await refresh();
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy('');
    }
  }
  function openEdit(
    kind: 'rename' | 'trash' | 'remove-project' | 'update-project',
    id: string,
    name: string,
  ) {
    setEdit({ kind, id, name });
    setEditName(kind === 'rename' ? name : '');
    setError('');
  }
  async function confirmEdit() {
    if (!edit) return;
    await act('edit', async () => {
      if (edit.kind === 'rename') {
        const changed = await workspaceApi.rename(edit.id, editName);
        onDatasetChanged(changed);
      } else if (edit.kind === 'trash') {
        await workspaceApi.trash(edit.id, editName);
        onDatasetChanged();
        setNotice('Moved to recoverable trash. Open Trash to restore it.');
      } else if (edit.kind === 'remove-project') {
        await workspaceApi.removeProject(edit.id);
        setNotice('Named project removed. All molecular datasets remain in the library.');
      } else {
        const current = state();
        if (!current) throw new Error('Load a structure first.');
        await workspaceApi.saveProject(edit.name, current, edit.id);
        setNotice('Project updated with the current scene.');
      }
      setEdit(null);
    });
  }
  const filtered =
    index?.datasets.filter(
      (d) =>
        (filter === 'archived' ? d.archived : !d.archived) &&
        `${d.name} ${d.source} ${d.id}`.toLowerCase().includes(query.toLowerCase()),
    ) ?? [];
  return (
    <div
      className="modal-backdrop"
      onMouseDown={(e) => {
        if (e.target === e.currentTarget && (!busy || busy.startsWith('open-'))) onClose();
      }}
    >
      <section
        className="modal workspace-library"
        role="dialog"
        aria-modal="true"
        aria-labelledby="workspace-library-title"
      >
        <div className="workspace-library__heading">
          <div>
            <span className="eyebrow">YOUR MOLECULAR WORKSPACE</span>
            <h2 id="workspace-library-title">Pick up where you left off.</h2>
            <p>
              Scenes save automatically on this computer. Named projects keep a deliberate snapshot.
            </p>
          </div>
          <button
            className="icon-button"
            aria-label="Close workspace library"
            disabled={!!busy && !busy.startsWith('open-')}
            onClick={onClose}
          >
            <X size={19} />
          </button>
        </div>
        <div className="workspace-library__tabs">
          <button
            aria-pressed={tab === 'projects'}
            onClick={() => {
              setTab('projects');
              setEdit(null);
            }}
          >
            Projects
          </button>
          <button
            aria-pressed={tab === 'library'}
            onClick={() => {
              setTab('library');
              setEdit(null);
            }}
          >
            Molecule library
          </button>
          <span role="status">{saveStatus}</span>
        </div>
        {error && (
          <p className="workspace-library__error" role="alert">
            {error}
          </p>
        )}
        {notice && (
          <p className="workspace-library__notice" role="status">
            <Check size={14} />
            {notice}
          </p>
        )}
        {tab === 'projects' ? (
          <>
            <form
              className="workspace-library__save"
              onSubmit={(e) => {
                e.preventDefault();
                void act('save', async () => {
                  const current = state();
                  if (!current) throw new Error('Load a structure first.');
                  await workspaceApi.saveProject(name, current);
                  setNotice('Project saved: scene, frame, plots and named selections.');
                });
              }}
            >
              <label>
                New project name
                <input
                  aria-label="New project name"
                  value={name}
                  maxLength={120}
                  onChange={(e) => setName(e.target.value)}
                  placeholder="e.g. Binding-site distances"
                />
              </label>
              <button className="primary-button" disabled={!!busy || !dataset || !name.trim()}>
                <Save size={15} />
                {busy === 'save' ? 'Saving…' : 'Save new project'}
              </button>
            </form>
            <div className="workspace-library__tools">
              <span>
                {projects.length} saved {projects.length === 1 ? 'project' : 'projects'}
              </span>
              <button
                className="secondary-button"
                disabled={!!busy}
                onClick={() => fileInput.current?.click()}
              >
                <Upload size={14} />
                Import project backup
              </button>
            </div>
            <input
              type="file"
              hidden
              ref={fileInput}
              aria-label="Import project backup file"
              accept=".zip"
              onChange={(e) => {
                const file = e.target.files?.[0];
                e.target.value = '';
                if (file)
                  void act('import', async () => {
                    const result = await workspaceApi.importProject(file);
                    onDatasetChanged();
                    setNotice(
                      `Imported “${result.project.name}” with ${result.imported_datasets.length} new and ${result.reused_datasets.length} existing datasets. Open it below.`,
                    );
                  });
              }}
            />
            <div className="workspace-library__rows">
              {projects.map((project) => (
                <div className="workspace-library__row" key={project.id}>
                  <FolderOpen size={19} />
                  <div className="workspace-library__description">
                    <b>{project.name}</b>
                    <small>
                      {project.unavailable
                        ? 'This project record is damaged. Its molecular files remain in the library.'
                        : `Saved ${new Date(project.updated_at).toLocaleString()}`}
                    </small>
                  </div>
                  <button
                    className="secondary-button"
                    disabled={!!busy || project.unavailable}
                    onClick={() =>
                      void act(`open-${project.id}`, async () => {
                        const full = await workspaceApi.project(project.id);
                        await onOpenProject(full.state);
                        onClose();
                      })
                    }
                  >
                    Open
                  </button>
                  <button
                    className="icon-button"
                    title="Update this snapshot from the current scene"
                    aria-label={`Update project ${project.name}`}
                    disabled={!!busy || !dataset || project.unavailable}
                    onClick={() => openEdit('update-project', project.id, project.name)}
                  >
                    <Save size={15} />
                  </button>
                  <button
                    className="icon-button"
                    title="Download project and all required molecular data"
                    aria-label={`Export project ${project.name}`}
                    disabled={!!busy || project.unavailable}
                    onClick={() =>
                      void act(`export-${project.id}`, async () => {
                        await workspaceApi.exportProject(project.id);
                        setNotice(
                          'Backup downloaded with molecular files, source ancestry, parameters and historical run records.',
                        );
                      })
                    }
                  >
                    <Download size={16} />
                  </button>
                  <button
                    className="icon-button"
                    aria-label={`Remove project ${project.name}`}
                    disabled={!!busy}
                    onClick={() => openEdit('remove-project', project.id, project.name)}
                  >
                    <Trash2 size={15} />
                  </button>
                </div>
              ))}
              {!projects.length && (
                <p className="workspace-library__empty">
                  Save your first project to keep this exact view, measurements and atom selections.
                </p>
              )}
            </div>
            <p className="workspace-library__footnote">
              Backups include related source structures and parameter files. Historical runs are
              retained as provenance. Import does not start a simulation. Limit: 500 MiB ZIP / 2 GiB
              expanded.
            </p>
          </>
        ) : (
          <>
            <div className="workspace-library__tools">
              <label className="workspace-library__search">
                <Search size={15} />
                <input
                  aria-label="Search molecule library"
                  value={query}
                  onChange={(e) => setQuery(e.target.value)}
                  placeholder="Search structures, sources…"
                />
              </label>
              <select
                aria-label="Library filter"
                value={filter}
                onChange={(e) => setFilter(e.target.value as typeof filter)}
              >
                <option value="active">Active</option>
                <option value="archived">Archived</option>
                <option value="trash">Trash</option>
              </select>
            </div>
            <div className="workspace-library__rows">
              {filter !== 'trash'
                ? filtered.map((d) => (
                    <div className="workspace-library__row" key={d.id}>
                      <FolderOpen size={19} />
                      <div className="workspace-library__description">
                        <b>
                          {d.name}
                          {d.id === dataset?.id && (
                            <span className="workspace-library__open">OPEN</span>
                          )}
                        </b>
                        <small>
                          {d.n_atoms.toLocaleString()} atoms · {d.n_frames.toLocaleString()} frames
                          · {formatBytes(d.bytes)}
                        </small>
                      </div>
                      <button
                        className="secondary-button"
                        disabled={!!busy && !busy.startsWith('open-')}
                        onClick={() =>
                          void act(`open-${d.id}`, async () => {
                            await onOpenDataset(d.id);
                            onClose();
                          })
                        }
                      >
                        Open
                      </button>
                      <button
                        className="icon-button"
                        aria-label={`Rename ${d.name}`}
                        disabled={!!busy}
                        onClick={() => openEdit('rename', d.id, d.name)}
                      >
                        <Pencil size={15} />
                      </button>
                      <button
                        className="icon-button"
                        aria-label={`${d.archived ? 'Unarchive' : 'Archive'} ${d.name}`}
                        disabled={!!busy}
                        onClick={() =>
                          void act(`archive-${d.id}`, async () => {
                            await workspaceApi.archive(d.id, !d.archived);
                            onDatasetChanged();
                          })
                        }
                      >
                        {d.archived ? <ArchiveRestore size={16} /> : <Archive size={16} />}
                      </button>
                      <button
                        className="icon-button"
                        aria-label={`Move ${d.name} to trash`}
                        disabled={!!busy}
                        onClick={() => openEdit('trash', d.id, d.name)}
                      >
                        <Trash2 size={15} />
                      </button>
                    </div>
                  ))
                : index?.trash
                    .filter((d) => d.name.toLowerCase().includes(query.toLowerCase()))
                    .map((d) => (
                      <div className="workspace-library__row" key={d.id}>
                        <Trash2 size={18} />
                        <div className="workspace-library__description">
                          <b>{d.name}</b>
                          <small>
                            {formatBytes(d.bytes)} · moved{' '}
                            {new Date(d.trashed_at).toLocaleDateString()}
                          </small>
                        </div>
                        <button
                          className="secondary-button"
                          disabled={!!busy}
                          onClick={() =>
                            void act(`restore-${d.id}`, async () => {
                              await workspaceApi.restore(d.id);
                              onDatasetChanged();
                              setNotice('Structure restored to the molecule library.');
                            })
                          }
                        >
                          <ArchiveRestore size={15} />
                          Restore
                        </button>
                      </div>
                    ))}
              {((filter === 'trash' && !index?.trash.length) ||
                (filter !== 'trash' && !filtered.length)) && (
                <p className="workspace-library__empty">
                  No matching {filter === 'trash' ? 'trashed' : filter} molecules.
                </p>
              )}
            </div>
            {index && (
              <div className="workspace-library__disk">
                <HardDrive size={17} />
                <div>
                  <b>Storage on this computer</b>
                  <p>
                    Structures {formatBytes(index.disk.datasets_bytes)} · Runs{' '}
                    {formatBytes(index.disk.jobs_bytes)} · Trash{' '}
                    {formatBytes(index.disk.trash_bytes)}
                    <br />
                    {formatBytes(index.disk.free_bytes)} free
                  </p>
                </div>
              </div>
            )}
            <p className="workspace-library__footnote">
              Archive hides a structure from the active list. Trash is recoverable and keeps its
              disk space; sources used by projects or runs are protected.
            </p>
          </>
        )}
        {edit && (
          <form
            className="workspace-library__confirm"
            onSubmit={(e) => {
              e.preventDefault();
              void confirmEdit();
            }}
          >
            <b>
              {edit.kind === 'rename'
                ? `Rename “${edit.name}”`
                : edit.kind === 'trash'
                  ? `Move “${edit.name}” to recoverable trash?`
                  : edit.kind === 'remove-project'
                    ? `Remove the “${edit.name}” snapshot?`
                    : `Replace “${edit.name}” with the current scene?`}
            </b>
            {edit.kind === 'rename' || edit.kind === 'trash' ? (
              <label>
                {edit.kind === 'trash'
                  ? `Type the exact name to confirm: ${edit.name}`
                  : 'New molecule name'}
                <input
                  autoFocus
                  aria-label={edit.kind === 'trash' ? 'Confirm dataset name' : 'Rename molecule'}
                  value={editName}
                  maxLength={120}
                  onChange={(e) => setEditName(e.target.value)}
                />
              </label>
            ) : (
              <p>
                {edit.kind === 'remove-project'
                  ? 'Molecular datasets remain in your library.'
                  : 'The previous camera, plots and selections in this snapshot will be replaced.'}
              </p>
            )}
            <div>
              <button
                type="button"
                className="secondary-button"
                disabled={!!busy}
                onClick={() => setEdit(null)}
              >
                Cancel
              </button>
              <button
                className="primary-button"
                disabled={
                  !!busy ||
                  (edit.kind === 'trash'
                    ? editName !== edit.name
                    : edit.kind === 'rename'
                      ? !editName.trim()
                      : false)
                }
              >
                {busy === 'edit' ? 'Working…' : 'Confirm'}
              </button>
            </div>
          </form>
        )}
        {busy && (
          <div className="workspace-library__busy" role="status">
            <LoaderCircle className="spin" size={16} />
            {busy === 'import'
              ? 'Validating and importing the project…'
              : busy.startsWith('export')
                ? 'Collecting sources, checking files and building your backup…'
                : 'Working…'}
          </div>
        )}
      </section>
    </div>
  );
}

export function NamedSelections({
  selections,
  selectedAtoms,
  onSave,
  onSelect,
  onRemove,
}: {
  selections: { id: string; name: string; atoms: number[] }[];
  selectedAtoms: number[];
  onSave: (name: string) => void;
  onSelect: (atoms: number[]) => void;
  onRemove: (id: string) => void;
}) {
  const [name, setName] = useState('');
  return (
    <div className="scene-section named-selections">
      <div className="section-title">
        <h3>Named selections</h3>
        <span>{selections.length}</span>
      </div>
      <form
        onSubmit={(e) => {
          e.preventDefault();
          if (selectedAtoms.length && name.trim()) {
            onSave(name.trim());
            setName('');
          }
        }}
      >
        <input
          aria-label="Selection name"
          placeholder="Name these atoms…"
          maxLength={120}
          value={name}
          onChange={(e) => setName(e.target.value)}
        />
        <button
          className="icon-button"
          title="Save current atom selection"
          aria-label="Save named selection"
          disabled={!selectedAtoms.length || !name.trim() || selections.length >= 100}
        >
          <Plus size={16} />
        </button>
      </form>
      <p>
        {selectedAtoms.length
          ? `${selectedAtoms.length.toLocaleString()} selected atoms`
          : 'Pick atoms or focus a residue, then save a selection.'}
      </p>
      {selections.map((selection) => (
        <div className="named-selections__row" key={selection.id}>
          <button title="Select and focus these atoms" onClick={() => onSelect(selection.atoms)}>
            <b>{selection.name}</b>
            <span>{selection.atoms.length} atoms</span>
          </button>
          <button
            className="icon-button"
            aria-label={`Remove selection ${selection.name}`}
            onClick={() => onRemove(selection.id)}
          >
            <X size={12} />
          </button>
        </div>
      ))}
    </div>
  );
}
