"""GNM slow-mode decomposition.

The Kirchhoff matrix of the GNM is the Laplacian of the residue contact network
(C-alpha within ``cutoff``).  Its slowest non-zero modes describe the largest-scale
relative motions; residues whose slow-mode components agree move together, and the
nodes (sign changes) of the slowest mode are the hinges.  We cluster the residues on
the ``k-1`` slowest modes (k-means, i.e. spectral clustering of the contact network)
and choose ``k`` by the eigengap when it is not given.

The matrix is small (one node per residue), so it is built and diagonalised with
numpy/scipy directly; ``prody.GNM`` gives the same modes (checked in the tests when
ProDy is installed).
"""

from __future__ import annotations

import numpy as np

from ..analysis import Analysis

DEFAULTS = {"n_regions": 0, "cutoff": 10.0, "max_auto": 8}


def kmeans(x: np.ndarray, k: int, n_init: int = 10, seed: int = 0, iters: int = 100) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = len(x)
    best = (np.inf, None)
    for _ in range(n_init):
        # k-means++ initialisation
        c = [x[rng.integers(n)]]
        for _ in range(1, k):
            d = np.min(((x[:, None, :] - np.asarray(c)[None]) ** 2).sum(-1), axis=1)
            p = d / d.sum() if d.sum() > 0 else np.full(n, 1 / n)
            c.append(x[rng.choice(n, p=p)])
        c = np.asarray(c)
        lab = np.zeros(n, dtype=np.int64)
        for _ in range(iters):
            d = ((x[:, None, :] - c[None]) ** 2).sum(-1)
            new = np.argmin(d, axis=1)
            if np.array_equal(new, lab) and _ > 0:
                break
            lab = new
            for j in range(k):
                if np.any(lab == j):
                    c[j] = x[lab == j].mean(axis=0)
        inertia = float(((x - c[lab]) ** 2).sum())
        if inertia < best[0]:
            best = (inertia, lab.copy())
    return best[1]


def hinges(an: Analysis, cutoff: float = 10.0) -> list[int]:
    """Residue indices where the slowest GNM mode changes sign (within a chain)."""
    vals, vecs, nz = an.gnm(cutoff)
    if vecs.shape[1] <= nz:
        return []
    v = vecs[:, nz]
    poly = an.gnm_nodes
    out = []
    chain_of = {}
    for c in an.s.chains():
        for i in range(c["start"], c["stop"]):
            chain_of[i] = c["start"]
    for j in range(len(poly) - 1):
        a, b = poly[j], poly[j + 1]
        if chain_of.get(a) == chain_of.get(b) and v[j] * v[j + 1] < 0:
            h = int(a if abs(v[j]) < abs(v[j + 1]) else b)
            if not out or out[-1] != h:
                out.append(h)
    return out


def choose_k(vals: np.ndarray, nz: int, kmax: int) -> tuple[int, list[float]]:
    lam = vals[nz:]
    gaps = []
    for k in range(2, kmax + 1):
        if k < lam.size and lam[k - 2] > 0:
            gaps.append(float(lam[k - 1] / lam[k - 2]))
        else:
            gaps.append(0.0)
    if not gaps:
        return 2, gaps
    return int(np.argmax(gaps)) + 2, gaps


def suggest(an: Analysis, params: dict) -> tuple[np.ndarray, dict]:
    p = {**DEFAULTS, **{k: v for k, v in params.items() if k in DEFAULTS}}
    poly = an.gnm_nodes
    labels = np.zeros(an.s.n_residues, dtype=np.int64)
    if poly.size < 4:
        return labels, {"notes": ["fewer than 4 polymer residues; nothing to decompose"]}
    vals, vecs, nz = an.gnm(float(p["cutoff"]), n_modes=max(20, int(p["max_auto"]) + 2))
    total_heavy = int(an.s.res_nheavy[poly].sum())
    kmax = int(max(2, min(int(p["max_auto"]), total_heavy // 150, 32)))
    k = int(p["n_regions"])
    gaps: list[float] = []
    if k <= 0:
        k, gaps = choose_k(vals, nz, kmax)
    k = max(2, min(k, len(poly)))
    # components: zero modes separate disconnected parts; slow modes split each part
    emb = vecs[:, nz : min(nz + k - 1, vecs.shape[1])]
    if nz > 1:
        emb = np.hstack([vecs[:, :nz] * 10.0, emb])
    # no row normalisation: residues with small slow-mode amplitude (near a hinge) should not be
    # pushed onto the unit sphere; benchmark in docs/gnm_multidomain.md
    lab = kmeans(emb, k)
    labels[poly] = lab + 1
    info = {
        "k": k,
        "eigengap": gaps,
        "eigenvalues": [float(x) for x in vals[nz : nz + 10]],
        "hinges": hinges(an, float(p["cutoff"])),
        "notes": [],
    }
    return labels, info
