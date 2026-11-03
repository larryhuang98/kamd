"""Move regions from the system a project was made on to another system by sequence alignment.

Typical use: groups are planned on a PDB/AlphaFold model, then the MD system is built
(tleap, CHARMM-GUI, pdb2gmx). The builder may renumber residues, drop or add residues
(missing loops, caps, termini), rename histidines or chains, and puts water, ions and lipids
around the protein, so residue indices change.

Procedure
---------
1. Polymer chains are taken from both systems (runs of protein / nucleic residues).
2. Every old chain is aligned to every new chain (global Needleman-Wunsch on one-letter
   codes, linear gap penalty); chains are paired greedily by score, preferring equal chain
   labels and equal chain order for identical copies (homo-oligomers).
3. Aligned residue pairs carry the region id over; mutations (aligned, different residue) are
   carried over and reported; residues in gaps are reported as lost or new (unassigned).
4. Ligands and other non-polymer residues in regions are matched by name and residue number.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .project import Project, Region, Settings
from .residues import ASSIGNABLE_KINDS

MATCH, MISMATCH, GAP = 2.0, -1.0, -2.0


def residue_table(system) -> dict:
    """Compact description of the assignable residues (stored in project.json)."""
    idx = [r.index for r in system.residues if r.kind in ASSIGNABLE_KINDS]
    res = system.residues
    return {
        "index": idx,
        "chain": [res[i].chain for i in idx],
        "number": [res[i].number for i in idx],
        "icode": [res[i].icode for i in idx],
        "name": [res[i].name for i in idx],
        "one": [res[i].one for i in idx],
        "kind": [res[i].kind for i in idx],
    }


def align(a: str, b: str) -> list[tuple[int | None, int | None]]:
    """Global alignment of two sequences; returns (i, j) pairs (None = gap)."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return [(i, None) for i in range(n)] + [(None, j) for j in range(m)]
    A = np.frombuffer(a.encode(), dtype=np.uint8)
    B = np.frombuffer(b.encode(), dtype=np.uint8)
    H = np.zeros((n + 1, m + 1))
    P = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up (gap in b), 2 left (gap in a)
    H[0, :] = GAP * np.arange(m + 1)
    H[:, 0] = GAP * np.arange(n + 1)
    P[0, 1:] = 2
    P[1:, 0] = 1
    ks = np.arange(m + 1)
    for i in range(1, n + 1):
        s = np.where(B == A[i - 1], MATCH, MISMATCH)
        diag = H[i - 1, :-1] + s
        up = H[i - 1, 1:] + GAP
        best = np.maximum(diag, up)
        src = np.where(diag >= up, 0, 1).astype(np.int8)
        # left moves: H[j] = max_k<=j (best[k] - GAP*k) + GAP*j, vectorised with a running max
        full = np.concatenate([[H[i, 0]], best])
        run = np.maximum.accumulate(full - GAP * ks) + GAP * ks
        H[i, 1:] = run[1:]
        P[i, 1:] = np.where(run[1:] > full[1:] + 1e-9, 2, src)
    out = []
    i, j = n, m
    while i > 0 or j > 0:
        d = P[i, j]
        if i > 0 and j > 0 and d == 0:
            out.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif i > 0 and (j == 0 or d == 1):
            out.append((i - 1, None))
            i -= 1
        else:
            out.append((None, j - 1))
            j -= 1
    return out[::-1]


def _chains(table: dict) -> list[dict]:
    """Polymer chains: runs of protein/nucleic residues with one chain label and consecutive index."""
    out: list[dict] = []
    for k, (g, c, kind) in enumerate(zip(table["index"], table["chain"], table["kind"], strict=True)):
        if kind not in ("protein", "nucleic"):
            continue
        if out and out[-1]["chain"] == c and table["index"][out[-1]["rows"][-1]] == g - 1:
            out[-1]["rows"].append(k)
        else:
            out.append({"chain": c, "rows": [k]})
    for n, ch in enumerate(out):
        ch["order"] = n
        ch["seq"] = "".join(table["one"][k] for k in ch["rows"])
    return out


@dataclass
class MigrationReport:
    chain_pairs: list[dict] = field(default_factory=list)
    carried: int = 0  # assigned residues mapped onto the new system
    lost: list[str] = field(default_factory=list)  # assigned residues with no counterpart
    mutated: list[str] = field(default_factory=list)  # aligned but different residue
    per_region: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def summary(self) -> list[str]:
        out = []
        for c in self.chain_pairs:
            out.append(
                f"chain {c['old'] or '-'} ({c['old_len']} res) -> chain {c['new'] or '-'} ({c['new_len']} res): "
                f"{c['identity']:.0%} identical, {c['gaps_old']} residues missing in the new system, "
                f"{c['gaps_new']} new residues"
            )
        out.append(f"{self.carried} assigned residues carried over")
        if self.lost:
            out.append(
                f"{len(self.lost)} assigned residues have no counterpart and were dropped: "
                + ", ".join(self.lost[:12])
                + (" ..." if len(self.lost) > 12 else "")
            )
        if self.mutated:
            out.append(
                f"{len(self.mutated)} aligned residues differ (mutation or renaming), kept in their region: "
                + ", ".join(self.mutated[:12])
                + (" ..." if len(self.mutated) > 12 else "")
            )
        return out + self.notes


def migrate(project: Project, new_system, old_table: dict | None = None) -> tuple[Project, MigrationReport]:
    old = old_table or project.residue_table
    if not old or "index" not in old:
        raise ValueError(
            "this project has no residue table (saved by an older kamd-regions); load it once together "
            "with the structure it was made on and save it again, then migrate"
        )
    new = residue_table(new_system)
    rep = MigrationReport()
    oc, nc = _chains(old), _chains(new)
    # score every chain pair, then pair greedily
    cand = []
    for a in oc:
        for b in nc:
            aln = align(a["seq"], b["seq"])
            same = sum(1 for i, j in aln if i is not None and j is not None and a["seq"][i] == b["seq"][j])
            ident = same / max(len(a["seq"]), 1)
            tie = (a["chain"] == b["chain"]) * 0.01 + (a["order"] == b["order"]) * 0.005
            cand.append((ident + tie, ident, a["order"], b["order"], aln))
    cand.sort(key=lambda c: -c[0])
    used_a, used_b = set(), set()
    mapping: dict[int, int] = {}  # old global residue index -> new global residue index
    for _, ident, ia, ib, aln in cand:
        if ia in used_a or ib in used_b or ident < 0.5:
            continue
        used_a.add(ia)
        used_b.add(ib)
        a, b = oc[ia], nc[ib]
        gaps_old = gaps_new = 0
        for i, j in aln:
            if i is not None and j is not None:
                ka, kb = a["rows"][i], b["rows"][j]
                mapping[old["index"][ka]] = new["index"][kb]
                if old["name"][ka] != new["name"][kb] and a["seq"][i] != b["seq"][j]:
                    rep.mutated.append(_label(old, ka) + "->" + new["name"][kb])
            elif i is not None:
                gaps_old += 1
            else:
                gaps_new += 1
        rep.chain_pairs.append(
            {
                "old": a["chain"],
                "new": b["chain"],
                "old_len": len(a["seq"]),
                "new_len": len(b["seq"]),
                "identity": ident,
                "gaps_old": gaps_old,
                "gaps_new": gaps_new,
            }
        )
    for ch in oc:
        if ch["order"] not in used_a:
            rep.notes.append(f"old chain {ch['chain'] or '-'} ({len(ch['seq'])} residues) matches no new chain")
    # non-polymer residues: by name + number
    new_other = {
        (n, num): g
        for g, n, num, k in zip(new["index"], new["name"], new["number"], new["kind"], strict=True)
        if k == "other"
    }
    for k, (g, kind) in enumerate(zip(old["index"], old["kind"], strict=True)):
        if kind == "other" and (old["name"][k], old["number"][k]) in new_other:
            mapping[g] = new_other[(old["name"][k], old["number"][k])]
    # carry the assignment
    assign = np.zeros(new_system.n_residues, dtype=np.int64)
    pos_old = {g: k for k, g in enumerate(old["index"])}
    for r in project.regions:
        rep.per_region[r.id] = {"old": 0, "carried": 0}
    for g in np.flatnonzero(project.assignment > 0):
        rid = int(project.assignment[g])
        if rid in rep.per_region:
            rep.per_region[rid]["old"] += 1
        tgt = mapping.get(int(g))
        if tgt is None:
            k = pos_old.get(int(g))
            rep.lost.append(_label(old, k) if k is not None else f"residue {g + 1}")
            continue
        assign[tgt] = rid
        rep.carried += 1
        if rid in rep.per_region:
            rep.per_region[rid]["carried"] += 1
    regions = [Region(r.id, r.name, r.color, r.region_type, r.refer_link) for r in project.regions]
    for r in regions:
        st = rep.per_region.get(r.id, {"old": 0, "carried": 0})
        if st["old"] and st["carried"] == 0:
            rep.notes.append(f"{r.name} has no residues in the new system")
    out = Project(
        n_residues=new_system.n_residues,
        regions=regions,
        assignment=assign,
        settings=Settings(**vars(project.settings)),
        inputs=dict(new_system.files),
        method={**project.method, "migrated_from": project.inputs},
        residue_names=[r.name for r in new_system.residues],
        residue_table=new,
        edited=project.edited,
    )
    return out, rep


def _label(table: dict, k: int) -> str:
    return f"{table['chain'][k] or '-'}:{table['name'][k]} {table['number'][k]}{table['icode'][k]}"
