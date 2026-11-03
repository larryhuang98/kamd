// Offline snapshot mode: the page carries precomputed backend results in window.__KAMD_DEMO__
// (written by `kamd-regions snapshot`). Validation of edits runs here with the same rules as
// kamd_regions/validate.py (sizes, membrane placement, links, bath); export needs the local server.
import type { ProjectDoc, StructureData, SystemInfo, Validation } from "./types";

export interface DemoData {
  session: { system: SystemInfo; project: ProjectDoc; validation: Validation; out_dir: string };
  structure: StructureData;
  suggest: Record<string, { project: ProjectDoc; notes: string[]; validation: Validation }>;
  title?: string;
}

export function demoData(): DemoData | null {
  return ((window as unknown as { __KAMD_DEMO__?: DemoData }).__KAMD_DEMO__) ?? null;
}

type Issue = Validation["errors"][number];
const ALL = ["amber", "gromacs", "openmm"];

export function validateLocal(sys: SystemInfo, p: ProjectDoc): Validation {
  const res = sys.residues;
  const pos = new Map<number, number>();
  res.index.forEach((g, i) => pos.set(g, i));
  const errors: Issue[] = [];
  const warnings: Issue[] = [];
  const stats: Validation["regions"] = {};
  const members = new Map<number, number[]>();
  for (const [a, b, v] of p.assignment) for (let g = a; g <= b; g++) {
    const q = pos.get(g);
    if (q !== undefined) (members.get(v) ?? members.set(v, []).get(v)!).push(q);
  }
  const noH = res.n_atoms.every((n, i) => n === res.n_heavy[i]);
  const st = p.settings;
  const err = (code: string, message: string, region: number | null = null, engines = ALL) => {
    errors.push({ code, message, engines, region });
    if (region !== null) stats[String(region)].issues.push({ level: "error", message });
  };
  const warn = (code: string, message: string, region: number | null = null, engines = ALL) => {
    warnings.push({ code, message, engines, region });
    if (region !== null) stats[String(region)].issues.push({ level: "warning", message });
  };
  for (const r of p.regions) {
    const m = members.get(r.id) ?? [];
    const inSlab = m.filter((q) => res.in_slab[q]).length;
    stats[String(r.id)] = {
      n_residues: m.length,
      n_atoms: m.reduce((s, q) => s + res.n_atoms[q], 0),
      n_heavy: m.reduce((s, q) => s + res.n_heavy[q], 0),
      summary: "",
      frac_in_slab: sys.membrane && m.length ? inSlab / m.length : null,
      issues: [],
    };
  }
  const nonempty = p.regions.filter((r) => stats[String(r.id)].n_residues > 0);
  for (const r of p.regions) if (!stats[String(r.id)].n_residues) err("empty_group", `${r.name} is empty (no residues)`, r.id);
  if (nonempty.length < 2) err("too_few_groups", `KAMD needs at least 2 reference groups; ${nonempty.length} defined`);
  if (p.regions.length > 32) err("too_many_groups", `Amber allows at most 32 groups; ${p.regions.length} defined`, null, ["amber"]);
  for (const r of nonempty) {
    const s = stats[String(r.id)];
    const n = noH ? 2 * s.n_heavy : s.n_atoms;
    const est = noH ? " (estimated with hydrogens)" : "";
    if (n < st.small_atoms) warn("small_group", `${r.name} has only ${n} atoms${est}; resets of small groups are noisy (a few hundred to a few thousand atoms is typical)`, r.id);
    if (n > st.large_atoms) warn("large_group", `${r.name} has ${n} atoms${est}; very large groups average out relative motion`, r.id);
  }
  const b = sys.bath;
  if (!b.solvent_atoms) {
    err("no_solvent", "no water or ions found: the solvent bath is required. Load the solvated system (prmtop + coordinates, or the solvated PDB/gro).", null, ["gromacs", "openmm"]);
    warn("no_solvent_amber", 'no water or ions in this structure: the Amber input uses the placeholder solvent_mask = ":WAT,Na+,Cl-"; check it against your prmtop', null, ["amber"]);
  }
  const hasTm = nonempty.some((r) => r.region_type === 2);
  const mem = st.membrane ?? sys.membrane;
  if (hasTm) {
    if (!b.membrane_atoms) err("no_membrane_bath", "transmembrane groups (region_type 2) need a membrane bath, but the system has no lipid residues" + (mem ? " (the OPM slab is used only for the suggestion)" : ""));
    if (!mem) err("no_membrane_axis", "transmembrane groups need membrane_axis, but no membrane was found; set the slab manually");
  }
  if (!b.solvent_contiguous) err("solvent_not_contiguous", "The solvent bath is not one contiguous atom range (Amber); reorder the system.", null, ["amber"]);
  if (hasTm && !b.membrane_contiguous) err("membrane_not_contiguous", "The membrane bath is not one contiguous atom range (Amber); reorder the system.", null, ["amber"]);
  if (mem) for (const r of nonempty) {
    const f = stats[String(r.id)].frac_in_slab;
    if (f === null) continue;
    if (r.region_type === 2 && f < 0.5) warn("tm_outside", `${r.name} is transmembrane but only ${Math.round(f * 100)}% of its residues lie in the membrane slab`, r.id);
    if (r.region_type === 1 && f > 0.6 && b.membrane_atoms) warn("solution_in_membrane", `${r.name} is a solution group but ${Math.round(f * 100)}% of its residues are buried in the membrane; consider region_type 2`, r.id);
  }
  const ids = new Set(p.regions.map((r) => r.id));
  const linked = p.regions.filter((r) => r.refer_link);
  for (const r of linked) {
    if (!ids.has(r.refer_link)) err("bad_link", `${r.name}: refer_link points to a deleted region`, r.id);
    else if (r.refer_link === r.id) err("self_link", `${r.name} cannot link to itself`, r.id);
  }
  if (linked.length) {
    const targets = new Set(linked.map((r) => r.refer_link).filter((k) => ids.has(k)));
    const free = nonempty.filter((r) => !r.refer_link);
    if (targets.size === 1 && free.length <= 1) warn("all_linked", "every group shares one direction through refer_link: the reset then pushes all groups the same way and no relative motion is excited");
    for (const r of linked) {
      const t = p.regions.find((q) => q.id === r.refer_link);
      if (t?.refer_link) warn("link_chain", `${r.name} links to ${t.name}, which itself links to another group; KAMD applies links in group order, so the result depends on order`, r.id);
    }
  }
  if (!sys.has_topology) warn("no_topology", "No topology loaded: residue numbers come from the structure file order. The Amber masks may not match your prmtop (and atom numbers may not match a GROMACS .tpr). Load the prmtop/top/tpr for exact export.");
  const ok: Record<string, boolean> = {};
  for (const e of ALL) ok[e] = !errors.some((x) => x.engines.includes(e));
  let unassigned = res.index.length;
  for (const m of members.values()) unassigned -= m.length;
  return { errors, warnings, regions: stats, ok, n_unassigned_protein: unassigned };
}
