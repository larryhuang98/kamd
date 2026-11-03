// Thin adapter so the 3D engine (3Dmol.js today, Mol* later) can be swapped.
// Residues are identified by the backend's global topology residue index.

export type Representation = "cartoon" | "surface" | "sticks";

export interface PickEvent {
  residue: number;
  shift: boolean;
  ctrl: boolean;
}

export interface SlabSpec {
  axis: number; // 1 x, 2 y, 3 z
  lo: number;
  hi: number;
}

export interface ViewerData {
  solutePdb: string;
  soluteRes: number[]; // atom (in file order) -> residue index
  ss: Map<number, string>; // residue -> H/E/C
  lipidPdb: string | null;
  ionPdb: string | null;
}

export interface MolViewer {
  load(data: ViewerData): void;
  setColors(color: (residue: number) => string): void;
  setRepresentation(rep: Representation): void;
  setSelection(residues: Set<number>, highlight: string): void;
  setMembrane(slab: SlabSpec | null, visible: boolean, color: string): void;
  setLipidsVisible(v: boolean): void;
  setIonsVisible(v: boolean): void;
  setBackground(color: string): void;
  onPick(cb: (e: PickEvent) => void): void;
  onHover(cb: (residue: number | null) => void): void;
  zoomTo(residues?: number[]): void;
  orientNormalUp(axis: number): void;
  resize(): void;
  dispose(): void;
}
