"""Per-residue secondary structure: 'H' helix, 'E' strand, 'C' coil, '' not protein.

Uses mdtraj's DSSP (simplified alphabet).  If mdtraj is missing or fails, falls back
to a coarse C-alpha geometry rule (i -> i+3 distance) that is good enough to find
helices and strands for the region heuristics.
"""

from __future__ import annotations

import warnings

import numpy as np

from .residues import STANDARD_NAME
from .system import System


def secondary_structure(s: System) -> np.ndarray:
    ss = np.array([""] * s.n_residues, dtype=object)
    prot = [r for r in s.residues if r.kind == "protein"]
    if not prot:
        return ss
    try:
        sub = _dssp_mdtraj(s, prot)
    except Exception as exc:  # pragma: no cover - fallback path
        warnings.warn(f"DSSP via mdtraj failed ({exc}); using a C-alpha geometry rule", stacklevel=2)
        sub = _ca_rule(s, prot)
    for r, c in zip(prot, sub, strict=True):
        ss[r.index] = c
    return ss


def _dssp_mdtraj(s: System, prot) -> list[str]:
    import mdtraj as md
    from mdtraj.core import element as mel

    top = md.Topology()
    keep = []
    chain = None
    last_chain = None
    for r in prot:
        if r.chain != last_chain:
            chain = top.add_chain()
            last_chain = r.chain
        res = top.add_residue(STANDARD_NAME.get(r.name.upper(), r.name.upper()), chain, resSeq=r.number)
        for a in range(r.start, r.stop):
            name = str(s.atom_name[a])
            if name in ("OT1", "OC1", "O1"):
                name = "O"
            if name == "HN":
                name = "H"
            el = str(s.element[a]) or "C"
            try:
                e = mel.get_by_symbol(el)
            except KeyError:
                e = mel.virtual
            top.add_atom(name, e, res)
            keep.append(a)
    # bonds are needed for nothing in DSSP; geometry only
    traj = md.Trajectory(s.xyz[keep][None, :, :] / 10.0, top)
    codes = md.compute_dssp(traj, simplified=True)[0]
    return ["C" if c in ("C", "NA", " ") else c for c in codes]


def _ca_rule(s: System, prot) -> list[str]:
    ca = []
    for r in prot:
        x = None
        for a in range(r.start, r.stop):
            if s.atom_name[a] == "CA":
                x = s.xyz[a]
                break
        ca.append(x if x is not None else s.xyz[r.start : r.stop].mean(axis=0))
    ca = np.asarray(ca)
    n = len(ca)
    out = ["C"] * n
    for i in range(n - 3):
        d3 = np.linalg.norm(ca[i + 3] - ca[i])
        if d3 < 6.0:
            for j in range(i, i + 4):
                out[j] = "H"
    for i in range(n - 2):
        if out[i] == "C" and np.linalg.norm(ca[i + 2] - ca[i]) > 6.4:
            out[i] = "E"
    return out


def segments(codes, value: str, min_len: int = 1) -> list[tuple[int, int]]:
    """Runs [start, stop) of `value` in a sequence of codes."""
    out = []
    i = 0
    n = len(codes)
    while i < n:
        if codes[i] == value:
            j = i
            while j < n and codes[j] == value:
                j += 1
            if j - i >= min_len:
                out.append((i, j))
            i = j
        else:
            i += 1
    return out
