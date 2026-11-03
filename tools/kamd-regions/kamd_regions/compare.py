"""Compare automatic suggestions with a manual (imported) grouping."""

from __future__ import annotations

import numpy as np

from .masks import residue_summary


def adjusted_rand(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a)
    b = np.asarray(b)
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    n = a.size
    if n < 2:
        return 1.0
    table = np.zeros((ai.max() + 1, bi.max() + 1))
    np.add.at(table, (ai, bi), 1)

    def c2(x):
        return x * (x - 1) / 2

    sum_ij = c2(table).sum()
    sa = c2(table.sum(1)).sum()
    sb = c2(table.sum(0)).sum()
    exp = sa * sb / c2(n)
    mx = 0.5 * (sa + sb)
    return float((sum_ij - exp) / (mx - exp)) if mx != exp else 1.0


def compare(session, manual, auto) -> dict:
    s = session.system
    an = session.an
    poly = an.polymer
    m = manual.assignment
    x = auto.assignment
    rows = []
    for r in manual.regions:
        res = manual.residues_of(r.id)
        atoms = int(s.res_natoms[res].sum())
        best = None
        for q in auto.regions:
            rq = auto.residues_of(q.id)
            inter = int(s.res_natoms[np.intersect1d(res, rq)].sum())
            union = int(s.res_natoms[np.union1d(res, rq)].sum())
            if inter and (best is None or inter > best[1]):
                best = (q, inter, union)
        cov = best[1] / atoms if best else 0.0
        jac = best[1] / best[2] if best else 0.0
        rows.append({
            "manual": r.name, "summary": residue_summary(s, res), "atoms": atoms,
            "best": best[0].name if best else "(unassigned)", "covered": cov, "jaccard": jac,
            "unassigned": int(s.res_natoms[res[x[res] == 0]].sum()),
        })
    auto_rows = []
    for q in auto.regions:
        rq = auto.residues_of(q.id)
        auto_rows.append({"name": q.name, "summary": residue_summary(s, rq), "atoms": int(s.res_natoms[rq].sum()),
                          "type": q.region_type})
    both = (m[poly] > 0) | (x[poly] > 0)
    return {
        "rows": rows,
        "auto": auto_rows,
        "ari_all": adjusted_rand(m[poly], x[poly]),
        "ari_assigned": adjusted_rand(m[poly][both], x[poly][both]) if both.any() else 1.0,
        "manual_atoms": int(s.res_natoms[poly][m[poly] > 0].sum()),
        "auto_atoms": int(s.res_natoms[poly][x[poly] > 0].sum()),
        "shared_atoms": int(s.res_natoms[poly][(m[poly] > 0) & (x[poly] > 0)].sum()),
    }


def _bar(p, s, poly) -> str:
    sym = ".123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return "".join(sym[min(int(p.assignment[i]), len(sym) - 1)] for i in poly)


def compare_report(session, methods, params=None) -> str:
    manual = session.project
    s = session.system
    poly = session.an.polymer
    lines = [
        "# kamd-regions: automatic vs. manual reference groups",
        "",
        f"System: {', '.join(v.split('/')[-1] for v in s.files.values())} — {s.n_atoms} atoms, "
        f"{len(poly)} polymer residues, {int(s.res_nheavy[poly].sum())} heavy atoms in the polymer.",
        "",
        "Manual groups (imported):",
        "",
        "| group | residues | atoms |",
        "|---|---|---|",
    ]
    for r in manual.regions:
        res = manual.residues_of(r.id)
        lines.append(f"| {r.name} | {residue_summary(s, res)} | {int(s.res_natoms[res].sum())} |")
    lines += ["", "Sequence view (one character per residue; `.` = unassigned):", "", "```",
              f"{'manual':>9s}  {_bar(manual, s, poly)}"]
    results = []
    for meth in methods:
        auto, notes = session.suggest(meth, params or {})
        results.append((meth, auto, notes, compare(session, manual, auto)))
        lines.append(f"{meth:>9s}  {_bar(auto, s, poly)}")
    lines += ["```", ""]
    for meth, auto, notes, c in results:
        lines += [
            f"## {meth}",
            "",
            f"{len(auto.regions)} groups. Adjusted Rand index vs. manual: {c['ari_all']:.2f} over all polymer "
            f"residues (unassigned counted as a class), {c['ari_assigned']:.2f} over residues assigned in either. "
            f"Atoms in groups: manual {c['manual_atoms']}, {meth} {c['auto_atoms']}, in both {c['shared_atoms']}.",
            "",
            "| auto group | type | residues | atoms |",
            "|---|---|---|---|",
        ]
        for a in c["auto"]:
            lines.append(f"| {a['name']} | {'TM' if a['type'] == 2 else 'solution'} | {a['summary']} | {a['atoms']} |")
        lines += [
            "",
            "| manual group | atoms | best-overlapping auto group | manual atoms covered | Jaccard | atoms left "
            "unassigned |",
            "|---|---|---|---|---|---|",
        ]
        for r in c["rows"]:
            lines.append(f"| {r['manual']} | {r['atoms']} | {r['best']} | {r['covered']:.0%} | {r['jaccard']:.2f} | "
                         f"{r['unassigned']} |")
        if notes:
            lines += ["", "Notes: " + " ".join(notes)]
        lines.append("")
    return "\n".join(lines)
