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
  parent_dataset_id?: string;
  preparation?: {
    ph: number;
    method?: string;
    summary?: string[];
    warnings?: string[];
    job_id?: string;
    ligand_parameters?: {
      forcefield: string;
      charge_method: string;
      requires_explicit_solvent: boolean;
      ligands: LigandInspection[];
    };
    [key: string]: unknown;
  };
  solvation?: {
    parent_dataset_id: string;
    padding_nm: number;
    water_model: string;
    added_water_atoms?: number;
    [key: string]: unknown;
  };
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
  visible?: boolean;
  occupancy?: number;
  warnings?: string[];
  angle_values?: number[];
}

export interface LigandInspection {
  key: string;
  component_id: string;
  chain: string;
  resid: string;
  residue: string;
  formal_charge?: number;
  selected_smiles?: string;
  protonation_method?: string;
  warnings?: string[];
  error?: string | null;
}

export interface Inspection {
  dataset_id: string;
  protein_atoms: number;
  hydrogen_atoms: number;
  water_atoms: number;
  heterogen_residues: string[];
  can_prepare: boolean;
  has_sequence: boolean;
  missing_atoms: { chain: string; resid: string; residue: string; atoms: string[] }[];
  missing_residues: {
    chain: string;
    position: number;
    residues: string[];
    count: number;
    terminal: boolean;
    buildable: boolean;
  }[];
  gaps: { chain: string; after: string; before: string; message: string }[];
  warnings: string[];
  blockers: string[];
  ligands: LigandInspection[];
  ligand_errors: string[];
  ligand_runtime?: { available: boolean; message?: string } | null;
  metal_environment?: { retained_coordinating_waters: string[][]; contacts: unknown[] };
}

export interface PreparationConfig {
  dataset_id: string;
  name?: string;
  ph: number;
  add_missing_atoms: boolean;
  build_missing_residues: boolean;
  optimize_sidechains: boolean;
  remove_waters: boolean;
  remove_heterogens: boolean;
  ligand_overrides?: Record<string, string>;
  seed: number;
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
