// Editable state with undo/redo.  Residues are addressed by their position in the list
// of assignable residues ("pos"); the backend's global topology index is residues.index[pos].
import { colorFor } from "./colors";
import type { ProjectDoc, Region, Residues, Settings } from "./types";

export interface Doc {
  regions: Region[];
  assign: Int32Array; // per pos: region id, 0 = unassigned
  settings: Settings;
  method: Record<string, unknown>;
  edited: boolean;
}

export interface History {
  past: Doc[];
  present: Doc;
  future: Doc[];
}

export function decode(p: ProjectDoc, res: Residues): Doc {
  const posOf = new Map<number, number>();
  res.index.forEach((g, i) => posOf.set(g, i));
  const assign = new Int32Array(res.index.length);
  for (const [a, b, v] of p.assignment) {
    for (let g = a; g <= b; g++) {
      const pos = posOf.get(g);
      if (pos !== undefined) assign[pos] = v;
    }
  }
  return {
    regions: p.regions.map((r) => ({ ...r })),
    assign,
    settings: { ...p.settings },
    method: p.method ?? {},
    edited: p.edited,
  };
}

export function encode(d: Doc, res: Residues, base: { inputs: Record<string, string>; n_residues: number; residue_names: string[] }): ProjectDoc {
  const runs: [number, number, number][] = [];
  let cur: [number, number, number] | null = null;
  for (let i = 0; i < d.assign.length; i++) {
    const v = d.assign[i];
    const g = res.index[i];
    if (v === 0) {
      cur = null;
      continue;
    }
    if (cur && cur[2] === v && cur[1] === g - 1) cur[1] = g;
    else {
      cur = [g, g, v];
      runs.push(cur);
    }
  }
  return {
    kamd_regions_project: 1,
    inputs: base.inputs,
    n_residues: base.n_residues,
    residue_names: base.residue_names,
    regions: d.regions,
    assignment: runs,
    settings: d.settings,
    method: d.method,
    edited: d.edited,
  };
}

// --------------------------------------------------------------------------- actions
export type Action =
  | { type: "replace"; doc: Doc; keepHistory?: boolean }
  | { type: "assign"; positions: number[]; region: number }
  | { type: "setAssign"; assign: Int32Array }
  | { type: "newRegion"; positions: number[] }
  | { type: "deleteRegions"; ids: number[] }
  | { type: "mergeRegions"; ids: number[] }
  | { type: "splitRegion"; id: number; positions: number[] }
  | { type: "updateRegion"; id: number; patch: Partial<Region> }
  | { type: "moveRegion"; id: number; delta: number }
  | { type: "settings"; patch: Partial<Settings> }
  | { type: "undo" }
  | { type: "redo" };

function nextId(d: Doc): number {
  return d.regions.reduce((m, r) => Math.max(m, r.id), 0) + 1;
}

function defaultName(d: Doc): string {
  const used = new Set(d.regions.map((r) => r.name));
  for (let k = d.regions.length + 1; ; k++) if (!used.has(`region ${k}`)) return `region ${k}`;
}

function usedColors(d: Doc): Set<string> {
  return new Set(d.regions.map((r) => r.color.toLowerCase()));
}

function freshColor(d: Doc): string {
  const used = usedColors(d);
  for (let i = 0; i < 64; i++) {
    const c = colorFor(i);
    if (!used.has(c.toLowerCase())) return c;
  }
  return colorFor(d.regions.length);
}

export function apply(d: Doc, a: Action): Doc {
  switch (a.type) {
    case "assign": {
      const assign = d.assign.slice();
      for (const p of a.positions) assign[p] = a.region;
      return { ...d, assign, edited: true };
    }
    case "setAssign":
      return { ...d, assign: a.assign, edited: true };
    case "newRegion": {
      const id = nextId(d);
      const r: Region = { id, name: defaultName(d), color: freshColor(d), region_type: 1, refer_link: 0 };
      const assign = d.assign.slice();
      for (const p of a.positions) assign[p] = id;
      return { ...d, regions: [...d.regions, r], assign, edited: true };
    }
    case "deleteRegions": {
      const del = new Set(a.ids);
      const assign = d.assign.map((v) => (del.has(v) ? 0 : v));
      const regions = d.regions
        .filter((r) => !del.has(r.id))
        .map((r) => (del.has(r.refer_link) ? { ...r, refer_link: 0 } : r));
      return { ...d, regions, assign, edited: true };
    }
    case "mergeRegions": {
      if (a.ids.length < 2) return d;
      const order = d.regions.filter((r) => a.ids.includes(r.id)).map((r) => r.id);
      const keep = order[0];
      const drop = new Set(order.slice(1));
      const assign = d.assign.map((v) => (drop.has(v) ? keep : v));
      const regions = d.regions
        .filter((r) => !drop.has(r.id))
        .map((r) => (drop.has(r.refer_link) ? { ...r, refer_link: r.id === keep ? 0 : keep } : r));
      return { ...d, regions, assign, edited: true };
    }
    case "splitRegion": {
      const pos = a.positions.filter((p) => d.assign[p] === a.id);
      if (!pos.length) return d;
      const src = d.regions.find((r) => r.id === a.id);
      const id = nextId(d);
      const r: Region = {
        id,
        name: `${src?.name ?? "region"} (split)`,
        color: freshColor(d),
        region_type: src?.region_type ?? 1,
        refer_link: 0,
      };
      const assign = d.assign.slice();
      for (const p of pos) assign[p] = id;
      const idx = d.regions.findIndex((x) => x.id === a.id);
      const regions = [...d.regions];
      regions.splice(idx + 1, 0, r);
      return { ...d, regions, assign, edited: true };
    }
    case "updateRegion":
      return { ...d, regions: d.regions.map((r) => (r.id === a.id ? { ...r, ...a.patch } : r)), edited: true };
    case "moveRegion": {
      const i = d.regions.findIndex((r) => r.id === a.id);
      const j = i + a.delta;
      if (i < 0 || j < 0 || j >= d.regions.length) return d;
      const regions = [...d.regions];
      [regions[i], regions[j]] = [regions[j], regions[i]];
      return { ...d, regions, edited: true };
    }
    case "settings":
      return { ...d, settings: { ...d.settings, ...a.patch } };
    default:
      return d;
  }
}

export function historyReducer(h: History, a: Action): History {
  if (a.type === "undo") {
    if (!h.past.length) return h;
    return { past: h.past.slice(0, -1), present: h.past[h.past.length - 1], future: [h.present, ...h.future] };
  }
  if (a.type === "redo") {
    if (!h.future.length) return h;
    return { past: [...h.past, h.present], present: h.future[0], future: h.future.slice(1) };
  }
  if (a.type === "replace") {
    if (a.keepHistory) return { past: [...h.past, h.present].slice(-200), present: a.doc, future: [] };
    return { past: [], present: a.doc, future: [] };
  }
  const next = apply(h.present, a);
  if (next === h.present) return h;
  return { past: [...h.past, h.present].slice(-200), present: next, future: [] };
}

// --------------------------------------------------------------------------- helpers
export function summary(positions: number[], res: Residues): string {
  if (!positions.length) return "—";
  const sorted = [...positions].sort((a, b) => a - b);
  const parts: string[] = [];
  let s = sorted[0];
  let p = s;
  const flush = () => {
    const c = res.chain[s] || "-";
    const a = `${res.number[s]}${res.icode[s] || ""}`;
    const b = `${res.number[p]}${res.icode[p] || ""}`;
    parts.push(p > s ? `${c}:${a}-${b}` : `${c}:${a}`);
  };
  for (let k = 1; k < sorted.length; k++) {
    const x = sorted[k];
    if (x === p + 1 && res.chain[x] === res.chain[p] && res.index[x] === res.index[p] + 1) {
      p = x;
      continue;
    }
    flush();
    s = p = x;
  }
  flush();
  return parts.join(", ");
}

export function positionsOf(d: Doc, id: number): number[] {
  const out: number[] = [];
  d.assign.forEach((v, i) => {
    if (v === id) out.push(i);
  });
  return out;
}

export function residueLabel(res: Residues, pos: number, d?: Doc, regions?: Region[]): string {
  const base = `${res.chain[pos] || "-"}:${res.name[pos]} ${res.number[pos]}${res.icode[pos] || ""}`;
  if (!d) return base;
  const v = d.assign[pos];
  if (!v) return `${base} (unassigned)`;
  const k = (regions ?? d.regions).findIndex((r) => r.id === v);
  const r = (regions ?? d.regions)[k];
  return `${base} (region ${k + 1}${r ? ": " + r.name : ""})`;
}
