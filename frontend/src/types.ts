export interface AtomInfo {
  index: number;
  name: string;
  element: string;
  residue: string;
  resid: number;
  chain: string;
  category: string;
  nonpolar_hydrogen: boolean;
}
export interface Dataset {
  id: string;
  name: string;
  n_atoms: number;
  n_residues: number;
  n_frames: number;
  times_ps: number[];
  time_unit: string;
  source: string;
  description: string;
  topology_url: string;
  coordinates_url: string;
  atoms: AtomInfo[];
  bonds: [number, number][];
  has_unitcell: boolean;
  warnings: string[];
}
export interface Visibility {
  protein: boolean;
  water: boolean;
  ligands: boolean;
  ions: boolean;
  hydrogens: 'all' | 'polar' | 'none';
}
export type Representation = 'cartoon' | 'ball+stick' | 'licorice' | 'surface';
export type MeasureKind = 'distance' | 'angle' | 'dihedral' | 'hbond';
export interface Measurement {
  id: string;
  kind: MeasureKind;
  atoms: number[];
  label: string;
  unit: string;
  values: number[];
  times_ps: number[];
  color: string;
  occupancy?: number;
  warnings?: string[];
  angle_values?: number[];
}
export interface Job {
  id: string;
  name: string;
  engine: string;
  status: string;
  stage: string;
  progress: number;
  completed_steps: number;
  total_steps: number;
  elapsed_seconds: number;
  logs: string[];
  dataset_id?: string;
  error?: string;
  config: Record<string, unknown>;
  created_at: string;
}
export interface Engine {
  id: string;
  name: string;
  available: boolean;
  version: string | null;
  message: string;
}
export interface Health {
  status: string;
  engines: Engine[];
}
export interface SimulationConfig {
  dataset_id: string;
  engine: 'openmm' | 'gromacs';
  name: string;
  duration_ps: number;
  temperature_k: number;
  timestep_fs: number;
  report_interval: number;
  friction_ps: number;
  seed: number;
  solvent: 'implicit' | 'explicit';
  minimize: boolean;
  equilibration_steps: number;
  padding_nm: number;
}
