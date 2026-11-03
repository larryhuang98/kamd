"""Contact-map communities (Leiden).

Nodes are polymer residues, edge weights are heavy-atom contact counts within
``cutoff`` plus a bonus for sequence neighbours.  Communities are found with the
Leiden algorithm (``leidenalg``, RB-configuration model) at the given resolution;
higher resolution gives more, smaller blocks.  Falls back to networkx's Louvain
when ``leidenalg``/``igraph`` are not installed.
"""

from __future__ import annotations

import numpy as np

from ..analysis import Analysis

DEFAULTS = {"resolution": 0.3, "cutoff": 4.5, "seq_weight": 2.0, "seed": 0}


def suggest(an: Analysis, params: dict) -> tuple[np.ndarray, dict]:
    p = {**DEFAULTS, **{k: v for k, v in params.items() if k in DEFAULTS}}
    poly = an.polymer
    labels = np.zeros(an.s.n_residues, dtype=np.int64)
    if poly.size < 2:
        return labels, {"notes": ["fewer than 2 polymer residues"]}
    m = an.contact_matrix(float(p["cutoff"]))
    sub = m[poly][:, poly].tocoo()
    w: dict[tuple[int, int], float] = {}
    for a, b, x in zip(sub.row, sub.col, sub.data, strict=True):
        if a < b:
            w[(int(a), int(b))] = w.get((int(a), int(b)), 0.0) + float(x)
    chains = {i: c["start"] for c in an.s.chains() for i in range(c["start"], c["stop"])}
    for j in range(len(poly) - 1):
        if chains.get(poly[j]) == chains.get(poly[j + 1]) and poly[j + 1] == poly[j] + 1:
            w[(j, j + 1)] = w.get((j, j + 1), 0.0) + float(p["seq_weight"])
    edges = list(w.keys())
    weights = [w[e] for e in edges]
    backend = "leiden"
    try:
        import igraph as ig
        import leidenalg

        g = ig.Graph(n=len(poly), edges=edges)
        g.es["weight"] = weights
        part = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=float(p["resolution"]),
            seed=int(p["seed"]),
            n_iterations=-1,
        )
        memb = np.asarray(part.membership, dtype=np.int64)
    except ImportError:  # pragma: no cover
        import networkx as nx

        backend = "louvain (networkx)"
        g = nx.Graph()
        g.add_nodes_from(range(len(poly)))
        g.add_weighted_edges_from((a, b, x) for (a, b), x in w.items())
        comms = nx.community.louvain_communities(
            g, weight="weight", resolution=float(p["resolution"]), seed=int(p["seed"])
        )
        memb = np.zeros(len(poly), dtype=np.int64)
        for i, c in enumerate(comms):
            memb[list(c)] = i
    labels[poly] = memb + 1
    return labels, {"backend": backend, "n_communities": int(memb.max() + 1), "notes": []}
