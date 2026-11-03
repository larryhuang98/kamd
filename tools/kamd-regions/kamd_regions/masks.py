"""Range compression and the residue / atom numbering conversions.

Conventions
-----------
* internal residue index: 0-based, sequential in the topology
* AmberMask ``:N``: 1-based sequential residue number (= internal + 1), no chains
* AmberMask ``@N`` and GROMACS index files: 1-based atom number (= atom index + 1)
* OpenMM / kamd_openmm: 0-based atom index
"""

from __future__ import annotations

import numpy as np


def compress(values) -> str:
    """Sorted integers -> '1-77,120-131,140' (input numbers are used as given)."""
    v = np.unique(np.asarray(list(values), dtype=np.int64))
    if v.size == 0:
        return ""
    parts = []
    start = prev = int(v[0])
    for x in v[1:]:
        x = int(x)
        if x == prev + 1:
            prev = x
            continue
        parts.append(f"{start}-{prev}" if prev > start else f"{start}")
        start = prev = x
    parts.append(f"{start}-{prev}" if prev > start else f"{start}")
    return ",".join(parts)


def expand(text: str) -> np.ndarray:
    """'1-77,120-131' -> array of ints (inverse of :func:`compress`)."""
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return np.asarray(out, dtype=np.int64)


def residue_mask(residue_indices) -> str:
    """0-based residue indices -> AmberMask ':1-77,120-131'."""
    r = np.asarray(list(residue_indices), dtype=np.int64)
    return ":" + compress(r + 1) if r.size else ""


def atom_mask(atom_indices) -> str:
    """0-based atom indices -> AmberMask '@1-100,200-300'."""
    a = np.asarray(list(atom_indices), dtype=np.int64)
    return "@" + compress(a + 1) if a.size else ""


def parse_simple_mask(mask: str) -> tuple[str, np.ndarray]:
    """Parse the masks this package writes (':ranges' or '@ranges'), returning 0-based numbers."""
    m = mask.strip()
    if m.startswith(":") and not any(c.isalpha() for c in m[1:]):
        return "residue", expand(m[1:]) - 1
    if m.startswith("@") and not any(c.isalpha() for c in m[1:]):
        return "atom", expand(m[1:]) - 1
    raise ValueError(f"not a simple numeric mask: {mask!r}")


def residue_summary(system, residue_indices) -> str:
    """Display form 'A:1-77, A:120-131' with chain IDs and PDB residue numbers."""
    idx = np.sort(np.asarray(list(residue_indices), dtype=np.int64))
    if idx.size == 0:
        return ""
    parts = []
    res = system.residues
    start = prev = int(idx[0])
    for x in list(idx[1:]) + [None]:
        if x is not None:
            x = int(x)
        if (
            x is not None
            and x == prev + 1
            and res[x].chain == res[prev].chain
            and res[x].number in (res[prev].number, res[prev].number + 1)  # a numbering gap starts a new range
        ):
            prev = x
            continue
        a, b = res[start], res[prev]
        c = a.chain or "-"
        parts.append(f"{c}:{a.number}{a.icode}" + (f"-{b.number}{b.icode}" if prev > start else ""))
        if x is not None:
            start = prev = x
    return ", ".join(parts)
