"""Membrane and bath detection.

* Bath: water + ions (solvent) and lipids (membrane), identified by residue name.
* Membrane normal: z is tried first; the choice is confirmed from the distribution of
  the lipid phosphorus atoms (two well separated layers).  Without lipids, the OPM
  ``DUM`` pseudo-atoms are used.
* Slab: the mean positions of the two phosphate layers (or the OPM hydrophobic
  boundaries).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .system import System

AXIS_NAMES = {1: "x", 2: "y", 3: "z"}


@dataclass
class Membrane:
    axis: int  # 1 = x, 2 = y, 3 = z (Amber membrane_axis)
    lo: float  # Angstrom, along axis
    hi: float
    source: str  # "lipids", "opm", "manual"
    score: float = 0.0
    n_lipids: int = 0

    @property
    def thickness(self) -> float:
        return self.hi - self.lo

    @property
    def center(self) -> float:
        return 0.5 * (self.lo + self.hi)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["axis_name"] = AXIS_NAMES[self.axis]
        d["thickness"] = self.thickness
        return d

    def inside(self, x: np.ndarray, margin: float = 0.0) -> np.ndarray:
        c = np.asarray(x)[..., self.axis - 1]
        return (c >= self.lo - margin) & (c <= self.hi + margin)


def _two_layers(c: np.ndarray) -> tuple[float, float, float]:
    """Best 1-D split into two clusters: (mean_lo, mean_hi, separation / pooled std)."""
    c = np.sort(np.asarray(c, dtype=np.float64))
    n = c.size
    if n < 4:
        return 0.0, 0.0, 0.0
    cs = np.cumsum(c)
    cs2 = np.cumsum(c * c)
    best = (np.inf, 0)
    for k in range(2, n - 1):
        m1, m2 = cs[k - 1] / k, (cs[-1] - cs[k - 1]) / (n - k)
        v1 = cs2[k - 1] / k - m1 * m1
        v2 = (cs2[-1] - cs2[k - 1]) / (n - k) - m2 * m2
        w = k * v1 + (n - k) * v2
        if w < best[0]:
            best = (w, k)
    k = best[1]
    lo, hi = c[:k], c[k:]
    pooled = np.sqrt(max(best[0] / n, 1e-6))
    return float(lo.mean()), float(hi.mean()), float((hi.mean() - lo.mean()) / pooled)


def detect_membrane(s: System) -> Membrane | None:
    lip_res = [r for r in s.residues if r.kind == "lipid"]
    if lip_res:
        p_atoms = [
            a
            for r in lip_res
            for a in range(r.start, r.stop)
            if str(s.element[a]) == "P" or str(s.atom_name[a]).upper() in ("P", "P8", "P31")
        ]
        if len(p_atoms) < 8:
            # cholesterol-only or no phosphate: use all lipid heavy atoms' extent
            p_atoms = [a for r in lip_res for a in range(r.start, r.stop) if s.heavy[a]]
        x = s.xyz[p_atoms].copy()
        prot = s.atoms_of_kind(["protein"])
        ref = s.xyz[prot].mean(axis=0) if prot.size else x.mean(axis=0)
        n_lip = len({r.molecule for r in lip_res})
        best = None
        for axis in (3, 1, 2):  # try z first
            d = x[:, axis - 1] - ref[axis - 1]
            if s.box is not None and s.box[axis - 1] > 0:
                L = s.box[axis - 1]
                d = d - L * np.round(d / L)
            lo, hi, score = _two_layers(d)
            cand = Membrane(axis, float(lo + ref[axis - 1]), float(hi + ref[axis - 1]), "lipids", float(score), n_lip)
            if axis == 3 and score > 8.0 and 25.0 < hi - lo < 60.0:
                return cand
            if best is None or score > best.score:
                best = cand
        if best is not None and best.score > 6.0 and 20.0 < best.thickness < 70.0:
            return best
        return None
    dum = [a for r in s.residues if r.kind == "dummy" for a in range(r.start, r.stop)]
    if dum:
        x = s.xyz[dum]
        for axis in (3, 1, 2):
            vals = np.unique(np.round(x[:, axis - 1], 1))
            if vals.size == 2 or (vals.size <= 4 and np.ptp(vals) > 10):
                return Membrane(axis, float(vals.min()), float(vals.max()), "opm", 99.0, 0)
        lo, hi, score = _two_layers(x[:, 2])
        return Membrane(3, lo, hi, "opm", score, 0)
    return None


# --------------------------------------------------------------------------------------
# bath
# --------------------------------------------------------------------------------------
def runs(sorted_idx: np.ndarray) -> list[tuple[int, int]]:
    """Contiguous runs [first, last] (inclusive) of a sorted integer array."""
    a = np.asarray(sorted_idx, dtype=np.int64)
    if a.size == 0:
        return []
    br = np.flatnonzero(np.diff(a) != 1)
    starts = np.concatenate([[0], br + 1])
    ends = np.concatenate([br, [a.size - 1]])
    return [(int(a[i]), int(a[j])) for i, j in zip(starts, ends, strict=True)]


@dataclass
class Bath:
    solvent_atoms: np.ndarray
    membrane_atoms: np.ndarray
    solvent_residues: np.ndarray
    membrane_residues: np.ndarray

    def contiguity(self, s: System) -> list[dict]:
        """Problems with the Amber 'one contiguous atom range' rule, with advice."""
        out = []
        for name, atoms, _res in (
            ("solvent", self.solvent_atoms, self.solvent_residues),
            ("membrane", self.membrane_atoms, self.membrane_residues),
        ):
            rr = runs(atoms)
            if len(rr) > 1:
                gaps = []
                for (_, e), (b, _) in zip(rr[:-1], rr[1:], strict=True):
                    between = sorted({s.residues[s.atom_res[k]].name for k in range(e + 1, b)})
                    gaps.append(f"atoms {e + 2}-{b} ({', '.join(between[:6])})")
                out.append(
                    {
                        "bath": name,
                        "n_ranges": len(rr),
                        "ranges": [(a + 1, b + 1) for a, b in rr[:10]],
                        "message": (
                            f"The {name} bath is split into {len(rr)} atom ranges; Amber needs one contiguous "
                            f"range. Interrupted by: {'; '.join(gaps[:3])}"
                            + (" ..." if len(gaps) > 3 else "")
                            + ". Reorder the system so that all "
                            + ("water and ion" if name == "solvent" else "lipid")
                            + " residues follow each other (e.g. in tleap: load protein, then "
                            + ("lipids, then water and ions" if name == "solvent" else "all lipids, then water/ions")
                            + "; with CPPTRAJ/ParmEd you can reorder residues and write a new prmtop)."
                        ),
                    }
                )
        return out


def detect_bath(s: System) -> Bath:
    kinds = s.res_kind
    sol_res = np.flatnonzero(np.isin(kinds, ["water", "ion"]))
    mem_res = np.flatnonzero(kinds == "lipid")
    return Bath(s.atoms_of(sol_res), s.atoms_of(mem_res), sol_res, mem_res)
