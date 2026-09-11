import type { Dataset, Measurement, MeasureKind, Representation, Visibility } from './types';
import type { AnalysisSettings } from './components/TrajectoryAnalysis';

export interface NamedSelection {
  id: string;
  name: string;
  atoms: number[];
}
export interface WorkspaceState {
  version: 1;
  dataset_id: string;
  atom_signature?: string;
  trajectory_signature?: string;
  frame: number;
  camera: number[] | null;
  visibility: Visibility;
  representation: Representation;
  color_scheme: 'chain' | 'residue' | 'element';
  measurements: Measurement[];
  active_measurement: string | null;
  selected_atoms: number[];
  named_selections: NamedSelection[];
  measure_kind: MeasureKind;
  inspector: boolean;
  speed: number;
  loop: boolean;
  analysis_settings?: AnalysisSettings | null;
  analysis_expanded?: boolean;
}
export interface SavedProject {
  id: string;
  name: string;
  created_at: string;
  updated_at: string;
  dataset_id: string;
  state: WorkspaceState;
  unavailable?: boolean;
}
export interface LibraryDataset extends Dataset {
  archived: boolean;
  bytes: number;
}
export interface TrashEntry {
  id: string;
  dataset_id: string;
  name: string;
  trashed_at: string;
  bytes: number;
}
export interface LibraryIndex {
  datasets: LibraryDataset[];
  trash: TrashEntry[];
  disk: { datasets_bytes: number; jobs_bytes: number; trash_bytes: number; free_bytes: number };
}
async function request<T>(url: string, body?: unknown): Promise<T> {
  const response = await fetch(
    url,
    body === undefined
      ? undefined
      : {
          method: 'POST',
          headers: body instanceof FormData ? undefined : { 'Content-Type': 'application/json' },
          body: body instanceof FormData ? body : JSON.stringify(body),
        },
  );
  if (!response.ok) {
    const error = await response.json().catch(() => null);
    throw new Error(error?.detail || `${response.status} ${response.statusText}`);
  }
  return response.json();
}
export const workspaceApi = {
  current: () => request<{ state: WorkspaceState | null; warning?: string }>('/api/workspace'),
  saveCurrent: (state: WorkspaceState) =>
    request<{ updated_at: string; atom_signature: string; trajectory_signature: string }>(
      '/api/workspace',
      { state },
    ),
  projects: () => request<SavedProject[]>('/api/projects'),
  project: (id: string) => request<SavedProject>(`/api/projects/${id}`),
  saveProject: (name: string, state: WorkspaceState, id?: string) =>
    request<SavedProject>('/api/projects', { name, state, id }),
  removeProject: (id: string) => request(`/api/projects/${id}/delete`, { confirm: true }),
  importProject: (file: File) => {
    const form = new FormData();
    form.append('file', file);
    return request<{
      project: SavedProject;
      imported_datasets: string[];
      reused_datasets: string[];
    }>('/api/projects/import', form);
  },
  async exportProject(id: string) {
    const response = await fetch(`/api/projects/${id}/export`);
    if (!response.ok)
      throw new Error(
        (await response.json().catch(() => null))?.detail || 'Project export failed.',
      );
    const blob = await response.blob();
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = `DynaMol-project-${id}.zip`;
    link.click();
    window.setTimeout(() => URL.revokeObjectURL(link.href), 30_000);
  },
  library: () => request<LibraryIndex>('/api/library'),
  rename: (id: string, name: string) => request<Dataset>(`/api/library/${id}/rename`, { name }),
  archive: (id: string, archived: boolean) => request(`/api/library/${id}/archive`, { archived }),
  trash: (id: string, confirm_name: string) =>
    request(`/api/library/${id}/trash`, { confirm_name }),
  restore: (id: string) => request<Dataset>(`/api/library/trash/${id}/restore`, {}),
};
export function formatBytes(bytes: number) {
  if (bytes >= 1024 ** 3) return `${(bytes / 1024 ** 3).toFixed(2)} GiB`;
  if (bytes >= 1024 ** 2) return `${(bytes / 1024 ** 2).toFixed(1)} MiB`;
  return `${Math.max(0.1, bytes / 1024).toFixed(1)} KiB`;
}
