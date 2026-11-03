"""Post-processing shared by every suggestion method.

1. sequence smoothing (runs shorter than ``min_run`` residues join a neighbour);
2. flexible termini and long, mobile coil loops become unassigned;
3. regions below the heavy-atom threshold are merged into the neighbour they touch most;
4. at most ``max_regions`` (Amber: 32);
5. region ids renumbered 1..K in sequence order.
"""

from __future__ import annotations

import numpy as np

from ..analysis import Analysis
from ..secstruct import segments

DEFAULTS = {
    "min_run": 3,  # residues
    "min_heavy": 150,  # heavy atoms per region (~300 atoms with hydrogens)
    "unassign_flexible": True,
    "loop_min_len": 8,  # interior coil runs at least this long can be dropped
    "flex_factor": 1.5,  # ... if their GNM fluctuation is this many times the median
    "max_regions": 32,
}


def chain_slices(an: Analysis) -> list[np.ndarray]:
    return [np.arange(c["start"], c["stop"]) for c in an.s.chains()]


def smooth(labels: np.ndarray, an: Analysis, min_run: int) -> np.ndarray:
    lab = labels.copy()
    if min_run <= 1:
        return lab
    for ch in chain_slices(an):
        for _ in range(3):
            seq = lab[ch]
            changed = False
            i = 0
            n = len(seq)
            while i < n:
                j = i
                while j < n and seq[j] == seq[i]:
                    j += 1
                if seq[i] != 0 and j - i < min_run and n > min_run:
                    left = seq[i - 1] if i > 0 else None
                    right = seq[j] if j < n else None
                    cand = [x for x in (left, right) if x is not None and x != 0]
                    if cand:
                        # prefer the longer neighbouring run
                        target = cand[0] if len(cand) == 1 or left == right else _longer_side(seq, i, j)
                        seq[i:j] = target
                        changed = True
                i = j
            lab[ch] = seq
            if not changed:
                break
    return lab


def _longer_side(seq, i, j):
    left_len = 0
    k = i - 1
    while k >= 0 and seq[k] == seq[i - 1]:
        left_len += 1
        k -= 1
    right_len = 0
    k = j
    while k < len(seq) and seq[k] == seq[j]:
        right_len += 1
        k += 1
    return seq[i - 1] if left_len >= right_len else seq[j]


def unassign_flexible(labels: np.ndarray, an: Analysis, loop_min_len: int, flex_factor: float) -> np.ndarray:
    lab = labels.copy()
    ss = an.ss
    msf = an.msf
    for ch in chain_slices(an):
        codes = ["C" if (ss[i] in ("C", "") and an.s.res_kind[i] == "protein") else "S" for i in ch]
        for a, b in segments(codes, "C"):
            idx = ch[a:b]
            terminal = a == 0 or b == len(ch)
            mobile = float(np.mean(msf[idx])) >= flex_factor
            if terminal and (b - a) >= 2 and mobile:
                # trim a terminus to its first structured residue
                lab[idx] = 0
            elif not terminal and (b - a) >= loop_min_len and mobile:
                lab[idx] = 0
    return lab


def heavy_per_label(labels: np.ndarray, an: Analysis) -> dict[int, int]:
    out: dict[int, int] = {}
    for g in np.unique(labels):
        if g > 0:
            out[int(g)] = int(an.s.res_nheavy[labels == g].sum())
    return out


def label_contacts(labels: np.ndarray, an: Analysis) -> dict[tuple[int, int], float]:
    m = an.contact_matrix().tocoo()
    la, lb = labels[m.row], labels[m.col]
    keep = (la > 0) & (lb > 0) & (la != lb)
    out: dict[tuple[int, int], float] = {}
    for a, b, w in zip(la[keep], lb[keep], m.data[keep], strict=True):
        out[(int(a), int(b))] = out.get((int(a), int(b)), 0.0) + float(w)
    # sequence neighbours also count as touching
    for ch in chain_slices(an):
        seq = labels[ch]
        for x, y in zip(seq[:-1], seq[1:], strict=True):
            if x > 0 and y > 0 and x != y:
                out[(int(x), int(y))] = out.get((int(x), int(y)), 0.0) + 5.0
                out[(int(y), int(x))] = out.get((int(y), int(x)), 0.0) + 5.0
    return out


def merge_small(
    labels: np.ndarray, an: Analysis, min_heavy: int, keep_at_least: int = 2, same_type: dict[int, int] | None = None
) -> tuple[np.ndarray, list[str]]:
    lab = labels.copy()
    notes: list[str] = []
    while True:
        heavy = heavy_per_label(lab, an)
        if len(heavy) <= keep_at_least:
            small = [g for g, h in heavy.items() if h < min_heavy]
            if small:
                notes.append(
                    f"{len(small)} region(s) are below {min_heavy} heavy atoms but were kept so that "
                    f"at least {keep_at_least} regions remain."
                )
            break
        g, h = min(heavy.items(), key=lambda kv: kv[1])
        if h >= min_heavy:
            break
        cont = label_contacts(lab, an)
        nb = {b: w for (a, b), w in cont.items() if a == g}
        if same_type:
            # never merge across region types (a soluble scrap must not grow a TM group); leave it unassigned
            nb = {b: w for b, w in nb.items() if same_type.get(b) == same_type.get(g)}
        if nb:
            target = max(nb.items(), key=lambda kv: (kv[1], -heavy[kv[0]]))[0]
            lab[lab == g] = target
        else:
            lab[lab == g] = 0
            notes.append(f"a region of {h} heavy atoms with no neighbour of its own type was left unassigned")
    return lab, notes


def cap_regions(labels: np.ndarray, an: Analysis, max_regions: int) -> np.ndarray:
    lab = labels.copy()
    while len(heavy_per_label(lab, an)) > max_regions:
        heavy = heavy_per_label(lab, an)
        g = min(heavy, key=heavy.get)
        cont = label_contacts(lab, an)
        nb = {b: w for (a, b), w in cont.items() if a == g}
        target = max(nb, key=nb.get) if nb else min((x for x in heavy if x != g), key=heavy.get)
        lab[lab == g] = target
    return lab


def renumber(labels: np.ndarray, extra: dict[int, dict] | None = None):
    mapping: dict[int, int] = {}
    out = np.zeros_like(labels)
    for i, g in enumerate(labels):
        if g > 0:
            out[i] = mapping.setdefault(int(g), len(mapping) + 1)
    if extra is not None:
        return out, {mapping[k]: v for k, v in extra.items() if k in mapping}
    return out


def postprocess(labels: np.ndarray, an: Analysis, params: dict, region_type: dict[int, int] | None = None):
    p = {**DEFAULTS, **{k: v for k, v in params.items() if k in DEFAULTS}}
    notes: list[str] = []
    lab = np.asarray(labels, dtype=np.int64).copy()
    lab[~an.assignable] = 0
    lab = smooth(lab, an, int(p["min_run"]))
    if p["unassign_flexible"]:
        lab = unassign_flexible(lab, an, int(p["loop_min_len"]), float(p["flex_factor"]))
    total = int(an.s.res_nheavy[lab > 0].sum())
    thr = int(p["min_heavy"])
    eff = min(thr, max(total // 3, 1))
    if eff < thr and total > 0:
        notes.append(
            f"The protein has only {total} heavy atoms in structured parts, so the minimum region size was "
            f"lowered from {thr} to {eff} heavy atoms to keep at least 2-3 regions."
        )
    lab, n2 = merge_small(lab, an, eff, same_type=region_type)
    notes += n2
    lab = cap_regions(lab, an, int(p["max_regions"]))
    if region_type is not None:
        lab, rt = renumber(lab, region_type)
        return lab, notes, rt
    return renumber(lab), notes, None
