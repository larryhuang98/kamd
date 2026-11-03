"""Per-system data shared by the suggestion methods (computed once, cached)."""

from __future__ import annotations

from functools import cached_property

import numpy as np
from scipy.spatial import cKDTree

from .membrane import Bath, Membrane, detect_bath, detect_membrane
from .secstruct import secondary_structure
from .system import System


class Analysis:
    def __init__(self, system: System, membrane: Membrane | None | str = "auto"):
        self.s = system
        self.membrane = detect_membrane(system) if membrane == "auto" else membrane
        self.bath: Bath = detect_bath(system)

    # ---------------------------------------------------------------- residues
    @cached_property
    def ss(self) -> np.ndarray:
        return secondary_structure(self.s)

    @cached_property
    def assignable(self) -> np.ndarray:
        return self.s.assignable()

    @cached_property
    def polymer(self) -> np.ndarray:
        """Residue indices used by the network methods (protein + nucleic)."""
        return np.flatnonzero(np.isin(self.s.res_kind, ["protein", "nucleic"]))

    @cached_property
    def gnm_nodes(self) -> np.ndarray:
        """Residues of the GNM network: protein residues (nucleic acids only if there is no protein).

        Nucleic-acid strands bound to a protein are left out: they are loosely connected to the
        protein network and their slow modes would dominate the decomposition.
        """
        prot = np.flatnonzero(self.s.res_kind == "protein")
        return prot if prot.size >= 4 else self.polymer

    @cached_property
    def rep_xyz(self) -> np.ndarray:
        """One representative point per residue (CA, C4', or heavy-atom centroid)."""
        s = self.s
        out = np.zeros((s.n_residues, 3))
        for r in s.residues:
            names = s.atom_name[r.start : r.stop]
            pick = None
            for want in ("CA", "C4'", "C4*", "P"):
                hit = np.flatnonzero(names == want)
                if hit.size and not (want == "CA" and r.n_atoms == 1):
                    pick = r.start + int(hit[0])
                    break
            if pick is not None:
                out[r.index] = s.xyz[pick]
            else:
                h = np.arange(r.start, r.stop)[s.heavy[r.start : r.stop]]
                out[r.index] = s.xyz[h if h.size else np.arange(r.start, r.stop)].mean(axis=0)
        return out

    # ---------------------------------------------------------------- contacts
    def contact_matrix(self, cutoff: float = 4.5):
        """Sparse symmetric (n_res x n_res) heavy-atom contact counts between assignable residues."""
        key = round(cutoff, 3)
        cache = self.__dict__.setdefault("_contacts", {})
        if key in cache:
            return cache[key]
        from scipy.sparse import coo_matrix

        s = self.s
        res = np.flatnonzero(self.assignable)
        atoms = s.atoms_of(res)
        atoms = atoms[s.heavy[atoms]]
        tree = cKDTree(s.xyz[atoms])
        pairs = tree.query_pairs(cutoff, output_type="ndarray")
        ra = s.atom_res[atoms[pairs[:, 0]]]
        rb = s.atom_res[atoms[pairs[:, 1]]]
        keep = ra != rb
        ra, rb = ra[keep], rb[keep]
        n = s.n_residues
        m = coo_matrix((np.ones(ra.size), (ra, rb)), shape=(n, n)).tocsr()
        m = (m + m.T).tocsr()
        cache[key] = m
        return m

    # ---------------------------------------------------------------- GNM
    def gnm(self, cutoff: float = 10.0, n_modes: int = 20):
        """GNM of :attr:`gnm_nodes`: (eigvals, eigvecs, n_zero) of the Kirchhoff matrix.

        Eigenvector rows follow the order of :attr:`gnm_nodes`.
        """
        key = (round(cutoff, 3), n_modes)
        cache = self.__dict__.setdefault("_gnm", {})
        if key in cache:
            return cache[key]
        x = self.rep_xyz[self.gnm_nodes]
        n = len(x)
        tree = cKDTree(x)
        pairs = tree.query_pairs(cutoff, output_type="ndarray")
        from scipy.sparse import coo_matrix

        a = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)).tocsr()
        a = a + a.T
        deg = np.asarray(a.sum(axis=1)).ravel()
        from scipy.sparse import diags

        kirchhoff = (diags(deg) - a).tocsr()
        from scipy.sparse.csgraph import connected_components

        ncomp, _ = connected_components(a, directed=False)
        k = min(n - 1, n_modes + ncomp) if n > 1 else 1
        if n <= 4000:
            from scipy.linalg import eigh

            vals, vecs = eigh(kirchhoff.toarray(), subset_by_index=[0, max(k - 1, 0)])
        else:  # pragma: no cover - large systems
            from scipy.sparse.linalg import eigsh

            vals, vecs = eigsh(kirchhoff.astype(float), k=k, sigma=-1e-3, which="LM")
            o = np.argsort(vals)
            vals, vecs = vals[o], vecs[:, o]
        vals = np.where(np.abs(vals) < 1e-8, 0.0, vals)
        nzero = max(int(np.sum(vals < 1e-6)), ncomp)
        cache[key] = (vals, vecs, nzero)
        return cache[key]

    @cached_property
    def msf(self) -> np.ndarray:
        """GNM mean-square fluctuation per residue, normalised to median 1 (0 for non-polymer)."""
        out = np.zeros(self.s.n_residues)
        if self.gnm_nodes.size < 4:
            return out
        vals, vecs, nz = self.gnm()
        v = vecs[:, nz:]
        lam = vals[nz:]
        if lam.size == 0:
            return out
        f = (v * v / lam[None, :]).sum(axis=1)
        med = np.median(f)
        out[self.gnm_nodes] = f / med if med > 0 else f
        return out

    # ---------------------------------------------------------------- membrane
    def in_slab(self, margin: float = 0.0) -> np.ndarray:
        """Per residue: representative point inside the membrane slab."""
        if self.membrane is None:
            return np.zeros(self.s.n_residues, dtype=bool)
        return self.membrane.inside(self.rep_xyz, margin)
