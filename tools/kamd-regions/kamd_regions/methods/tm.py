"""Transmembrane segments as reference groups (region_type = 2).

DSSP helices and strands that cross the membrane slab are transmembrane segments.
A segment counts as crossing when enough of its residues sit inside the slab and its
extent along the membrane normal covers a good part of the slab thickness (this keeps
amphipathic interface helices out).  Each segment, or each run of ``per_group``
neighbouring segments (``0`` = as many as needed to reach the size threshold),
becomes one group.  Short turns between segments of one group are included.

Parts of the protein outside the membrane that are large enough (for example a
soluble domain) are decomposed with GNM and become solution groups.
"""

from __future__ import annotations

import numpy as np

from ..analysis import Analysis
from ..secstruct import segments
from . import gnm

DEFAULTS = {
    "per_group": 0,
    "min_heavy_tm": 150,
    "turn_max": 6,
    "soluble": True,
    "slab_margin": 2.0,
}


def tm_segments(an: Analysis, margin: float = 2.0) -> list[dict]:
    if an.membrane is None:
        return []
    mem = an.membrane
    ax = mem.axis - 1
    inside = an.in_slab(margin)
    x = an.rep_xyz
    out = []
    for c in an.s.chains():
        idx = np.arange(c["start"], c["stop"])
        codes = ["H" if an.ss[i] in ("H", "G", "I") else "E" if an.ss[i] in ("E", "B") else "C" for i in idx]
        # close beta bulges / 1-2 residue DSSP breaks inside a strand
        for k in range(1, len(codes) - 1):
            for w in (1,):
                if (
                    k + w < len(codes)
                    and codes[k - 1] == "E"
                    and codes[k + w] == "E"
                    and all(c == "C" for c in codes[k : k + w])
                ):
                    codes[k : k + w] = ["E"] * w
        for kind, min_len, min_in in (("H", 10, 8), ("E", 5, 4)):
            for a, b in segments(codes, kind, min_len):
                r = idx[a:b]
                n_in = int(inside[r].sum())
                span = float(np.ptp(x[r, ax]))
                if n_in >= min_in and span >= 0.45 * mem.thickness:
                    out.append({"start": int(r[0]), "stop": int(r[-1]) + 1, "kind": kind, "chain": c["chain"]})
    out.sort(key=lambda d: d["start"])
    return out


def suggest(an: Analysis, params: dict) -> tuple[np.ndarray, dict]:
    p = {**DEFAULTS, **{k: v for k, v in params.items() if k in DEFAULTS}}
    labels = np.zeros(an.s.n_residues, dtype=np.int64)
    info: dict = {"notes": [], "region_type": {}}
    if an.membrane is None:
        info["notes"].append("no membrane found (no lipids and no OPM DUM atoms); the tm method needs a slab")
        return labels, info
    segs = tm_segments(an, float(p["slab_margin"]))
    info["tm_segments"] = segs
    if not segs:
        info["notes"].append("no transmembrane helix or strand crosses the membrane slab")
        return labels, info
    nheavy = an.s.res_nheavy
    groups: list[list[dict]] = []
    per = int(p["per_group"])
    for seg in segs:
        if groups:
            last = groups[-1]
            same_chain = last[-1]["chain"] == seg["chain"]
            if per > 0:
                full = len(last) >= per
            else:
                full = sum(int(nheavy[s["start"] : s["stop"]].sum()) for s in last) >= int(p["min_heavy_tm"])
            if same_chain and not full:
                last.append(seg)
                continue
        groups.append([seg])
    # a small trailing group joins its predecessor (auto mode)
    if per <= 0 and len(groups) > 2:
        tail = groups[-1]
        if sum(int(nheavy[s["start"] : s["stop"]].sum()) for s in tail) < int(p["min_heavy_tm"]) // 2:
            if groups[-2][-1]["chain"] == tail[0]["chain"]:
                groups[-2].extend(tail)
                groups.pop()
    rid = 0
    rtype = {}
    for g in groups:
        rid += 1
        rtype[rid] = 2
        for j, s in enumerate(g):
            labels[s["start"] : s["stop"]] = rid
            if j + 1 < len(g):
                gap = g[j + 1]["start"] - s["stop"]
                if 0 < gap <= int(p["turn_max"]):
                    labels[s["stop"] : g[j + 1]["start"]] = rid
    # soluble parts outside the membrane
    if p["soluble"]:
        outside = (~an.in_slab(6.0)) & (labels == 0) & np.isin(an.s.res_kind, ["protein", "nucleic"])
        heavy_out = int(nheavy[outside].sum())
        if heavy_out >= 2 * 300:
            g_lab, g_info = gnm.suggest(an, {"n_regions": max(2, heavy_out // 1500)})
            sol = np.where(outside, g_lab, 0)
            for v in np.unique(sol[sol > 0]):
                members = sol == v
                if int(nheavy[members].sum()) >= 300:
                    rid += 1
                    rtype[rid] = 1
                    labels[members] = rid
            info["notes"].append("extramembrane parts were decomposed with GNM into solution groups")
    info["region_type"] = rtype
    info["n_tm_segments"] = len(segs)
    return labels, info
