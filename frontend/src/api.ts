import type {
  Dataset,
  Health,
  Job,
  MeasureKind,
  Measurement,
  SimulationConfig,
  Inspection,
  PreparationConfig,
  RunMeasurement,
  RunMeasurementSnapshot,
  MonomerOptions,
} from './types';
import type { MeasurementPreview } from './components/LiveMeasurement';
async function request<T>(url: string, options?: RequestInit): Promise<T> {
  const response = await fetch(url, options);
  if (!response.ok) {
    let message = `${response.status} ${response.statusText}`;
    try {
      const body = await response.json();
      message = typeof body.detail === 'string' ? body.detail : JSON.stringify(body.detail ?? body);
    } catch {
      /* response not JSON */
    }
    throw new Error(message);
  }
  return response.json();
}
export const api = {
  health: () => request<Health>('/api/health'),
  demo: () => request<Dataset>('/api/datasets/demo'),
  datasets: () => request<Dataset[]>('/api/datasets'),
  dataset: (id: string) => request<Dataset>(`/api/datasets/${id}`),
  monomers: (id: string) => request<MonomerOptions>(`/api/datasets/${id}/monomers`),
  useMonomer: (id: string, chain_index: number, keep_associated_molecules = true) =>
    request<Dataset>(`/api/datasets/${id}/monomer`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ chain_index, keep_associated_molecules }),
    }),
  async coordinates(dataset: Dataset) {
    const r = await fetch(dataset.coordinates_url);
    if (!r.ok) throw new Error('Could not load trajectory coordinates.');
    const data = new Float32Array(await r.arrayBuffer());
    if (data.length !== dataset.n_frames * dataset.n_atoms * 3)
      throw new Error('Coordinate count does not match topology.');
    return data;
  },
  upload: (data: FormData) =>
    request<Dataset>('/api/datasets/upload', { method: 'POST', body: data }),
  importStructure: (data: FormData) =>
    request<Dataset>('/api/structures/upload', { method: 'POST', body: data }),
  fetchStructure: (provider: 'pdb' | 'pubchem', identifier: string) =>
    request<Dataset>('/api/structures/fetch', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ provider, identifier }),
    }),
  smiles: (smiles: string, name?: string) =>
    request<Dataset>('/api/structures/smiles', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles, name }),
    }),
  inspect: (id: string, ph = 7, overrides: Record<string, string> = {}) =>
    Object.keys(overrides).length
      ? request<Inspection>(`/api/datasets/${id}/inspection`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ ph, ligand_overrides: overrides }),
        })
      : request<Inspection>(`/api/datasets/${id}/inspection?ph=${ph}`),
  prepare: (config: PreparationConfig) =>
    request<Job>('/api/preparations', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(config),
    }),
  solvate: (id: string, padding_nm: number, ph: number, seed: number) =>
    request<Job>(`/api/datasets/${id}/solvate`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ padding_nm, ph, seed }),
    }),
  measure: (id: string, kind: MeasureKind, atoms: number[], signal?: AbortSignal) =>
    request<Omit<Measurement, 'id' | 'label' | 'color'>>(`/api/datasets/${id}/measurements`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ kind, atoms }),
      signal,
    }),
  previewMeasurement: (id: string, kind: MeasureKind, atoms: number[], signal: AbortSignal) =>
    request<MeasurementPreview>(`/api/datasets/${id}/measurement-preview`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ kind, atoms }),
      signal,
    }),
  jobs: () => request<Job[]>('/api/jobs'),
  jobMeasurements: (id: string, signal?: AbortSignal) =>
    request<RunMeasurementSnapshot>(`/api/jobs/${id}/measurements`, { signal }),
  remapMeasurements: (target: string, source: string, measurements: RunMeasurement[]) =>
    request<{ measurements: RunMeasurement[]; warnings: string[]; errors: string[] }>(
      `/api/datasets/${target}/remap-measurements`,
      {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ source_dataset_id: source, measurements }),
      },
    ),
  start: (config: SimulationConfig) =>
    request<Job>('/api/jobs', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(config),
    }),
  cancel: (id: string) => request<Job>(`/api/jobs/${id}/cancel`, { method: 'POST' }),
};
