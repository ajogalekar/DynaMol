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
  monomer_selection?: {
    parent_dataset_id: string;
    chain_index: number;
    chain_id: string;
    chain_label: string;
    retained_atom_indices: number[];
    excluded_protein_chains: unknown[];
    retained_associated_molecules: unknown[];
    excluded_molecules: unknown[];
    association_cutoff_angstrom: number;
    requires_preparation: boolean;
    notes: string[];
  };
  preparation?: {
    ph: number;
    method?: string;
    summary?: string[];
    warnings?: string[];
    job_id?: string;
    modified_residues?: ModifiedResidueInspection[];
    requires_explicit_solvent?: boolean;
    recommend_explicit_solvent?: boolean;
    solvent_recommendation?: string;
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
export interface MonomerOptions {
  chains: {
    index: number;
    chain_id: string;
    label: string;
    n_residues: number;
    n_atoms: number;
    sequence_group?: string;
    recommended?: boolean;
  }[];
  recommended_chain_index: number | null;
  warnings: string[];
}
export type Representation = 'cartoon' | 'ball+stick' | 'licorice' | 'surface';
export type MeasureKind = 'distance' | 'angle' | 'dihedral' | 'hbond';
export interface Measurement {
  id: string;
  kind: MeasureKind;
  atoms: number[];
  label: string;
  unit: string;
  values: (number | null)[];
  times_ps: number[];
  color: string;
  visible?: boolean;
  occupancy?: number | null;
  warnings?: string[];
  angle_values?: (number | null)[];
  frame_errors?: (string | null)[];
  trackDuringRun?: boolean;
}

export type RunMeasurement = Pick<Measurement, 'id' | 'kind' | 'atoms' | 'label' | 'color'>;
export interface RunMeasurementSnapshot {
  job_id: string;
  status: string;
  source_dataset_id: string;
  output_dataset_id?: string | null;
  measurements: (RunMeasurement & {
    output_atoms?: number[] | null;
    values: (number | null)[];
    times_ps: number[];
    unit: string;
    frame_errors: (string | null)[];
    warnings: string[];
    angle_values?: (number | null)[];
    occupancy?: number | null;
  })[];
  warnings: string[];
  errors: string[];
  last_time_ps: number | null;
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
  missing_heavy_atoms?: string[];
  extra_heavy_atoms?: string[];
  can_repair?: boolean;
  can_remove?: boolean;
  repair_reason?: string;
  selected_action?: 'repair' | 'remove' | 'keep';
  removed?: boolean;
}

export interface Inspection {
  dataset_id: string;
  protein_atoms: number;
  hydrogen_atoms: number;
  water_atoms: number;
  heterogen_residues: string[];
  can_prepare: boolean;
  has_sequence: boolean;
  loop_policy?: {
    max_gap_residues: number;
    max_total_residues: number;
    short_gap_residues: number;
  };
  missing_atoms: { chain: string; resid: string; residue: string; atoms: string[] }[];
  missing_residues: {
    chain: string;
    position: number;
    residues: string[];
    count: number;
    terminal: boolean;
    buildable: boolean;
    modeling?: 'short' | 'extended' | 'unsupported' | 'terminal';
  }[];
  gaps: {
    chain: string;
    after: string;
    before: string;
    message: string;
    structural_break?: boolean;
  }[];
  warnings: string[];
  blockers: string[];
  ligands: LigandInspection[];
  modified_residues?: ModifiedResidueInspection[];
  ligand_errors: string[];
  ions?: {
    key: string;
    residue: string;
    element: string;
    formal_charge: number;
    supported: boolean;
    error?: string | null;
    model?: string;
  }[];
  ligand_runtime?: { available: boolean; message?: string } | null;
  metal_environment?: { retained_coordinating_waters: string[][]; contacts: unknown[] };
}

export interface ModifiedResidueInspection {
  name?: string;
  chain: string;
  resid: string;
  insertion_code?: string;
  residue: string;
  supported: boolean;
  parent_residue?: string;
  template?: string;
  formal_charge?: number;
  forcefield?: string;
  protonation_note?: string;
  error?: string;
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
  ligand_actions?: Record<string, 'repair' | 'remove'>;
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
  measurements?: RunMeasurement[];
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
