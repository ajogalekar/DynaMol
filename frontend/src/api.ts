import type { Dataset, Health, Job, MeasureKind, Measurement, SimulationConfig } from './types';
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
  measure: (id: string, kind: MeasureKind, atoms: number[]) =>
    request<Omit<Measurement, 'id' | 'label' | 'color'>>(`/api/datasets/${id}/measurements`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ kind, atoms }),
    }),
  jobs: () => request<Job[]>('/api/jobs'),
  start: (config: SimulationConfig) =>
    request<Job>('/api/jobs', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(config),
    }),
  cancel: (id: string) => request<Job>(`/api/jobs/${id}/cancel`, { method: 'POST' }),
};
