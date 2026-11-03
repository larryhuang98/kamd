export interface Region {
  id: number;
  name: string;
  color: string;
  region_type: number; // 1 solution, 2 transmembrane
  refer_link: number; // region id, 0 = own direction
}

export interface Settings {
  pseudobd_cycles: number;
  nzero_crossings: number;
  shell_cut: number;
  surface_cut: number;
  membrane_cut: number;
  kamd_shell_mode: number | null;
  temperature: number;
  seed: number | null;
  membrane: { axis: number; lo: number; hi: number } | null;
  small_atoms: number;
  large_atoms: number;
}

export interface ProjectDoc {
  kamd_regions_project: number;
  inputs: Record<string, string>;
  n_residues: number;
  residue_names: string[];
  residue_table?: Record<string, unknown>; // sequence of the system the project was made on (migration)
  regions: Region[];
  assignment: [number, number, number][];
  settings: Settings;
  method: Record<string, unknown>;
  edited: boolean;
}

export interface Residues {
  index: number[];
  name: string[];
  one: string[];
  chain: string[];
  number: number[];
  icode: string[];
  kind: string[];
  ss: string[];
  n_atoms: number[];
  n_heavy: number[];
  in_slab: boolean[];
  tm: boolean[];
  msf: number[];
}

export interface Membrane {
  axis: number;
  axis_name: string;
  lo: number;
  hi: number;
  source: string;
  thickness: number;
  n_lipids: number;
}

export interface SystemInfo {
  version: string;
  files: Record<string, string>;
  n_atoms: number;
  n_residues: number;
  has_topology: boolean;
  topology_format: string;
  kinds: Record<string, number>;
  chains: { chain: string; start: number; stop: number }[];
  residues: Residues;
  membrane: Membrane | null;
  bath: {
    solvent_atoms: number;
    membrane_atoms: number;
    solvent_contiguous: boolean;
    membrane_contiguous: boolean;
    n_lipids: number;
  };
  hinges: number[];
  methods: Record<string, { label: string; params: Record<string, number | boolean> } | Record<string, number | boolean>>;
  notes: string[];
}

export interface Issue {
  code: string;
  message: string;
  engines: string[];
  region: number | null;
}

export interface RegionStats {
  n_residues: number;
  n_atoms: number;
  n_heavy: number;
  summary: string;
  frac_in_slab: number | null;
  issues: { level: string; message: string }[];
}

export interface Validation {
  errors: Issue[];
  warnings: Issue[];
  regions: Record<string, RegionStats>;
  ok: Record<string, boolean>;
  n_unassigned_protein: number;
}

export interface StructureData {
  solute_pdb: string;
  solute_res: number[];
  lipid_pdb: string | null;
  ion_pdb: string | null;
}

export interface VerificationItem {
  item: string;
  ok: boolean;
  n_atoms?: number;
  message?: string;
}

export interface ExportResult {
  ok: boolean;
  blocked: string[];
  validation: Validation;
  files: Record<string, string>;
  verification: {
    ok: boolean;
    kamd_config: { ok: boolean | null; message: string };
    amber: VerificationItem[];
    gromacs: VerificationItem[];
    openmm: VerificationItem[];
  } | null;
  written?: string[];
}
