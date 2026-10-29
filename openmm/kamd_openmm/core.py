"""Engine-independent KAMD algorithm (numpy only).

This module is a direct port of the pmemd.cuda KAMD implementation
(ambers/kamd, branch kamd_membrane_shell_balance, src/pmemd/src/runmd.F90 and
src/pmemd/src/cuda/kKAMD.h).  It knows nothing about OpenMM; the engine layer
(``kamd_openmm.reporter``) feeds it numpy arrays.

Units used throughout this module (OpenMM's native units):
    length nm, time ps, mass amu (g/mol), energy kJ/mol,
    velocity nm/ps, force kJ/mol/nm.
In these units F/m is directly in nm/ps^2 and kB*T/m in nm^2/ps^2.

User-facing cutoffs in :class:`KAMDConfig` are in Angstrom, the same numbers
you would put into the Amber ``&kamd`` namelist.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

KB_KJ_PER_MOL_K = 0.00831446261815324  # Boltzmann constant, kJ/(mol K)

REGION_SOLUTION = 1
REGION_TRANSMEMBRANE = 2


class KAMDInputError(ValueError):
    """Invalid KAMD setup (mirrors the &kamd validation errors in pmemd)."""


# --------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------
@dataclass
class KAMDConfig:
    """KAMD parameters.  Names follow the Amber ``&kamd`` namelist.

    Index conventions
    -----------------
    * Atom indices are **0-based** (OpenMM convention).
    * ``refer_link`` and ``membrane_axis`` keep the Amber **1-based** meaning so an
      Amber input can be transcribed literally: ``refer_link[i] = 0`` means group
      ``i+1`` draws its own direction, ``refer_link[i] = k`` means it shares the
      direction of group ``k``.  ``membrane_axis`` is 1=x, 2=y, 3=z.
    * ``region_type`` per group: 1 = solution, 2 = transmembrane.
    """

    refer_groups: Sequence[Sequence[int]]
    solvent_atoms: Sequence[int]
    temperature: float  # K, target temperature for the Maxwell-Boltzmann resample
    pseudobd_cycles: int = 250
    nzero_crossings: int = 2
    shell_cut: float = 5.0  # Angstrom
    surface_cut: float = 3.0  # Angstrom
    refer_link: Sequence[int] | None = None
    region_type: Sequence[int] | None = None
    membrane_atoms: Sequence[int] | None = None
    membrane_axis: int = 0
    membrane_cut: float = 3.0  # Angstrom
    seed: int | None = None

    # filled by validate()
    groups: list = field(default_factory=list, init=False, repr=False)

    def validate(self, n_atoms: int | None = None) -> KAMDConfig:
        if self.pseudobd_cycles <= 0:
            raise KAMDInputError("pseudobd_cycles must be > 0 (0 disables KAMD; just don't attach it)")
        if self.nzero_crossings <= 0:
            raise KAMDInputError("nzero_crossings must be > 0")
        if self.shell_cut <= 0 or self.surface_cut <= 0:
            raise KAMDInputError("shell_cut and surface_cut must be > 0")
        if self.surface_cut > self.shell_cut:
            raise KAMDInputError("surface_cut must be <= shell_cut")
        if self.temperature < 0:
            raise KAMDInputError("temperature must be >= 0")

        groups = [np.unique(np.asarray(g, dtype=np.int64)) for g in self.refer_groups]
        if len(groups) == 0:
            raise KAMDInputError("at least one reference group must be given")
        for i, g in enumerate(groups):
            if g.size == 0:
                raise KAMDInputError(f"reference group {i + 1} is empty")
        ng = len(groups)

        link = np.zeros(ng, dtype=np.int64) if self.refer_link is None else np.asarray(self.refer_link, dtype=np.int64)
        if link.shape != (ng,):
            raise KAMDInputError(f"refer_link must have one entry per group ({ng})")
        for g in range(ng):
            if link[g] < 0 or link[g] > ng:
                raise KAMDInputError(f"refer_link({g + 1}) = {link[g]} is out of range [0, {ng}]")
            if link[g] == g + 1:
                raise KAMDInputError(f"refer_link({g + 1}) cannot reference its own group")

        rtype = (
            np.full(ng, REGION_SOLUTION, dtype=np.int64)
            if self.region_type is None
            else np.asarray(self.region_type, dtype=np.int64)
        )
        if rtype.shape != (ng,):
            raise KAMDInputError(f"region_type must have one entry per group ({ng})")
        for g in range(ng):
            if rtype[g] not in (REGION_SOLUTION, REGION_TRANSMEMBRANE):
                raise KAMDInputError(f"region_type({g + 1}) = {rtype[g]} is invalid (1=solution, 2=transmembrane)")

        solvent = np.unique(np.asarray(self.solvent_atoms, dtype=np.int64))
        if solvent.size == 0:
            raise KAMDInputError("solvent_atoms is required")
        membrane = (
            np.zeros(0, dtype=np.int64)
            if self.membrane_atoms is None
            else np.unique(np.asarray(self.membrane_atoms, dtype=np.int64))
        )

        has_tm = bool(np.any(rtype == REGION_TRANSMEMBRANE))
        if has_tm:
            if membrane.size == 0:
                raise KAMDInputError("membrane_atoms is required when any region_type == 2")
            if self.membrane_axis not in (1, 2, 3):
                raise KAMDInputError(f"membrane_axis = {self.membrane_axis} is invalid (1=x, 2=y, 3=z)")
            if self.membrane_cut <= 0:
                raise KAMDInputError("membrane_cut must be > 0")

        # overlap checks (pmemd rejects all of these)
        allref = np.concatenate(groups)
        uniq, counts = np.unique(allref, return_counts=True)
        if np.any(counts > 1):
            a = int(uniq[counts > 1][0])
            owners = [i + 1 for i, g in enumerate(groups) if a in set(g.tolist())]
            raise KAMDInputError(f"atom {a} is in reference groups {owners}; reference groups must not overlap")
        if np.intersect1d(allref, solvent).size:
            raise KAMDInputError("reference groups and solvent bath overlap")
        if membrane.size:
            if np.intersect1d(allref, membrane).size:
                raise KAMDInputError("reference groups and membrane bath overlap")
            if np.intersect1d(solvent, membrane).size:
                raise KAMDInputError("solvent bath and membrane bath overlap")
        if n_atoms is not None:
            for name, arr in (("reference", allref), ("solvent", solvent), ("membrane", membrane)):
                if arr.size and (arr.min() < 0 or arr.max() >= n_atoms):
                    raise KAMDInputError(f"{name} atom index out of range for a {n_atoms}-atom system")

        self.groups = groups
        self._link = link
        self._rtype = rtype
        self._solvent = solvent
        self._membrane = membrane
        return self

    # convenience accessors (valid after validate())
    @property
    def ngroups(self) -> int:
        return len(self.groups)

    @property
    def link(self) -> np.ndarray:
        return self._link

    @property
    def rtype(self) -> np.ndarray:
        return self._rtype

    @property
    def solvent(self) -> np.ndarray:
        return self._solvent

    @property
    def membrane(self) -> np.ndarray:
        return self._membrane


# --------------------------------------------------------------------------------------
# Momentum monitor  (runmd.F90: "Per-cycle momentum monitoring")
# --------------------------------------------------------------------------------------
class MomentumMonitor:
    """Zero-crossing detector for the group momenta.

    Exactly reproduces the pmemd logic, including its quirks:
      * the reference velocity is the COM velocity of *all reference groups*;
      * for transmembrane groups the membrane-axis component of the relative
        momentum is zeroed before projecting;
      * after a reset ``last`` is zero, so the first check after a reset can never
        register a crossing (product == 0 is not < 0);
      * with a single group the relative momentum is zero up to round-off, so
        the crossings of a one-group setup are floating-point noise (same formula
        as pmemd); use two or more groups.
    """

    def __init__(self, group_masses: np.ndarray, region_type: np.ndarray, membrane_axis: int, nzero_crossings: int):
        self.group_masses = np.asarray(group_masses, dtype=np.float64)
        self.region_type = np.asarray(region_type)
        self.membrane_axis = int(membrane_axis)
        self.nzero = int(nzero_crossings)
        self.reset()

    def reset(self) -> None:
        ng = self.group_masses.size
        self.last = np.zeros(ng)
        self.count = np.zeros(ng, dtype=np.int64)

    def projections(self, group_momenta: np.ndarray) -> np.ndarray:
        p = np.asarray(group_momenta, dtype=np.float64)  # (ng, 3)
        vcom = p.sum(axis=0) / self.group_masses.sum()
        rel = p - self.group_masses[:, None] * vcom[None, :]
        tm = self.region_type == REGION_TRANSMEMBRANE
        if np.any(tm):
            rel[tm, self.membrane_axis - 1] = 0.0
        return rel @ vcom

    def update(self, group_momenta: np.ndarray) -> bool:
        """Register one monitoring check.  Returns True if a reset must fire now."""
        s = self.projections(group_momenta)
        crossed = (s * self.last < 0.0) & (self.count < self.nzero)
        self.count[crossed] += 1
        self.last = s
        return bool(np.all(self.count >= self.nzero))

    def state_dict(self) -> dict:
        return {"last": self.last.tolist(), "count": self.count.tolist()}

    def load_state_dict(self, d: dict) -> None:
        self.last = np.asarray(d["last"], dtype=np.float64)
        self.count = np.asarray(d["count"], dtype=np.int64)


def group_momenta(velocities: np.ndarray, masses: np.ndarray, groups: Sequence[np.ndarray]) -> np.ndarray:
    """Translational momentum sum(m v) of each group, shape (ngroups, 3)."""
    return np.stack([(masses[g, None] * velocities[g]).sum(axis=0) for g in groups])


# --------------------------------------------------------------------------------------
# Random directions  (runmd.F90 "Step 2" + TM projection)
# --------------------------------------------------------------------------------------
def group_directions(
    rng: np.random.Generator, link: np.ndarray, region_type: np.ndarray, membrane_axis: int
) -> np.ndarray:
    """One unit vector per group, with refer_link sharing and TM in-plane projection.

    Same order of operations as pmemd: draw isotropic vectors, apply refer_link in
    group order (so a link to a later group copies that group's *own* draw), then
    project transmembrane groups onto the membrane plane and renormalise.
    """
    ng = link.size
    d = rng.standard_normal((ng, 3))
    for g in range(ng):
        if link[g] > 0:
            d[g] = d[link[g] - 1]
    norms = np.linalg.norm(d, axis=1)
    d = d / np.where(norms > 0, norms, 1.0)[:, None]
    for g in range(ng):
        if region_type[g] == REGION_TRANSMEMBRANE:
            ax = membrane_axis - 1
            d[g, ax] = 0.0
            w = np.linalg.norm(d[g])
            if w > 1e-10:
                d[g] /= w
            else:  # pathological: draw was along the membrane normal
                d[g] = 0.0
                d[g, 1 if ax == 0 else 0] = 1.0
    return d


# --------------------------------------------------------------------------------------
# First-shell detection  (kKAMD_FindFirstShell_kernel)
# --------------------------------------------------------------------------------------
def _is_orthorhombic(box: np.ndarray) -> bool:
    off = box - np.diag(np.diag(box))
    return bool(np.all(np.abs(off) < 1e-12 * max(1.0, np.abs(box).max())))


def _min_image(d: np.ndarray, box: np.ndarray, box_inv: np.ndarray) -> np.ndarray:
    """Fractional-coordinate minimum image, same convention as kKAMD.h (rows = vectors)."""
    f = d @ box_inv
    f -= np.round(f)
    return f @ box


def _nearest_bruteforce(query: np.ndarray, points: np.ndarray, box: np.ndarray, cut: float, chunk: int = 2048):
    """For each query point: (index of nearest point, squared distance) or (-1, inf) if none < cut."""
    box_inv = np.linalg.inv(box)
    best_i = np.full(len(query), -1, dtype=np.int64)
    best_r2 = np.full(len(query), np.inf)
    cut2 = cut * cut
    for s in range(0, len(query), chunk):
        q = query[s : s + chunk]
        for t in range(0, len(points), chunk):
            p = points[t : t + chunk]
            d = _min_image(q[:, None, :] - p[None, :, :], box, box_inv)
            r2 = np.einsum("ijk,ijk->ij", d, d)
            j = np.argmin(r2, axis=1)
            r2min = r2[np.arange(len(q)), j]
            better = (r2min < best_r2[s : s + chunk]) & (r2min < cut2)
            best_i[s : s + chunk][better] = t + j[better]
            best_r2[s : s + chunk][better] = r2min[better]
    return best_i, best_r2


def _nearest(query: np.ndarray, points: np.ndarray, box: np.ndarray, cut: float, use_kdtree: bool):
    if len(points) == 0 or len(query) == 0:
        return np.full(len(query), -1, dtype=np.int64), np.full(len(query), np.inf)
    if use_kdtree:
        from scipy.spatial import cKDTree

        L = np.diag(box).copy()
        tree = cKDTree(np.mod(points, L) % L, boxsize=L)
        # small pad on the bound, exact strict-< test applied afterwards
        r, idx = tree.query(np.mod(query, L) % L, k=1, distance_upper_bound=cut * (1 + 1e-9))
        ok = np.isfinite(r) & (r * r < cut * cut)
        idx = np.where(ok, idx, -1).astype(np.int64)
        r2 = np.where(ok, r * r, np.inf)
        return idx, r2
    return _nearest_bruteforce(query, points, box, cut)


@dataclass
class ShellResult:
    """Output of the first-shell search.

    ``bath_atoms[k]`` is assigned to 0-based group ``bath_group[k]``;
    ``surface_atoms[k]`` is a surface reference atom of group ``surface_group[k]``.
    """

    bath_atoms: np.ndarray
    bath_group: np.ndarray
    surface_atoms: np.ndarray
    surface_group: np.ndarray

    def counts(self, ngroups: int):
        return (np.bincount(self.bath_group, minlength=ngroups), np.bincount(self.surface_group, minlength=ngroups))


def find_first_shell(positions: np.ndarray, box: np.ndarray, cfg: KAMDConfig, method: str = "auto") -> ShellResult:
    """First-shell bath atoms and surface reference atoms.

    positions : (N, 3) nm.  box : (3, 3) nm, rows are the periodic box vectors.

    Semantics of kKAMD_FindFirstShell_kernel:
      * solution groups search the solvent bath with ``shell_cut`` (assignment) and
        ``surface_cut`` (surface marking); transmembrane groups search the membrane
        bath with ``membrane_cut`` for both;
      * a bath atom goes to the group owning its nearest reference atom (r < cut);
      * a reference atom is "surface" if any bath atom of its region is at r < cut;
      * all tests are strict ``<`` with the minimum-image convention.
    """
    if method == "auto":
        method = "kdtree" if _is_orthorhombic(box) else "bruteforce"
        if method == "kdtree":
            try:
                import scipy.spatial  # noqa: F401
            except ImportError:
                method = "bruteforce"
    use_kd = method == "kdtree"
    if use_kd and not _is_orthorhombic(box):
        raise ValueError("kdtree search requires an orthorhombic box")

    to_nm = 0.1
    bath_atoms, bath_group, surf_atoms, surf_group = [], [], [], []
    for region, bath, acut, scut in (
        (REGION_SOLUTION, cfg.solvent, cfg.shell_cut * to_nm, cfg.surface_cut * to_nm),
        (REGION_TRANSMEMBRANE, cfg.membrane, cfg.membrane_cut * to_nm, cfg.membrane_cut * to_nm),
    ):
        gids = [g for g in range(cfg.ngroups) if cfg.rtype[g] == region]
        if not gids or bath.size == 0:
            continue
        ref = np.concatenate([cfg.groups[g] for g in gids])
        ref_gid = np.concatenate([np.full(cfg.groups[g].size, g) for g in gids])
        ref_x = positions[ref]
        bath_x = positions[bath]

        # bath atom -> nearest reference atom within the assignment cutoff
        j, _ = _nearest(bath_x, ref_x, box, acut, use_kd)
        hit = j >= 0
        bath_atoms.append(bath[hit])
        bath_group.append(ref_gid[j[hit]])

        # reference atom with any bath atom within the surface cutoff
        k, _ = _nearest(ref_x, bath_x, box, scut, use_kd)
        s = k >= 0
        surf_atoms.append(ref[s])
        surf_group.append(ref_gid[s])

    cat = lambda xs: np.concatenate(xs).astype(np.int64) if xs else np.zeros(0, dtype=np.int64)  # noqa: E731
    return ShellResult(cat(bath_atoms), cat(bath_group), cat(surf_atoms), cat(surf_group))


# --------------------------------------------------------------------------------------
# Velocity reset + redirect  (kRV.h + kKAMD_RedirectVelocities_kernel)
# --------------------------------------------------------------------------------------
def maxwell_boltzmann(rng: np.random.Generator, masses: np.ndarray, temperature: float) -> np.ndarray:
    """Full-step velocities drawn from Maxwell-Boltzmann (nm/ps).  Massless particles get 0."""
    v = np.zeros((masses.size, 3))
    if temperature < 1e-6:
        return v
    m = masses > 0
    sd = np.sqrt(KB_KJ_PER_MOL_K * temperature / masses[m])
    v[m] = sd[:, None] * rng.standard_normal((int(m.sum()), 3))
    return v


def redirect(
    v_full: np.ndarray,
    atoms: np.ndarray,
    atom_group: np.ndarray,
    directions: np.ndarray,
    region_type: np.ndarray,
    membrane_axis: int,
) -> None:
    """Speed-preserving redirect of full-step velocities, in place.

    Solution group: v <- d_g |v|.
    Transmembrane group: the component along the membrane axis is kept, the in-plane
    part becomes d_g |v_perp| (d_g already lies in the plane).
    """
    if atoms.size == 0:
        return
    v = v_full[atoms]
    d = directions[atom_group]
    tm = region_type[atom_group] == REGION_TRANSMEMBRANE
    out = np.empty_like(v)
    sol = ~tm
    if np.any(sol):
        out[sol] = d[sol] * np.linalg.norm(v[sol], axis=1)[:, None]
    if np.any(tm):
        ax = membrane_axis - 1
        vpar = v[tm, ax]
        vperp = np.sqrt(np.maximum((v[tm] ** 2).sum(axis=1) - vpar**2, 0.0))
        out[tm] = d[tm] * vperp[:, None]
        out[tm, ax] += vpar
    v_full[atoms] = out


def pseudobd_reset(
    rng: np.random.Generator,
    positions: np.ndarray,
    forces: np.ndarray,
    masses: np.ndarray,
    box: np.ndarray,
    cfg: KAMDConfig,
    dt: float,
    half_step_velocities: bool = True,
    method: str = "auto",
):
    """One complete pseudo-BD reset.  Returns (new_velocities, ShellResult, directions).

    pmemd sequence (runmd.F90 "Pseudo-Brownian dynamics reset"):
      1. resample *all* velocities from Maxwell-Boltzmann at T0 (as full-step values);
      2. draw one direction per group (refer_link, TM projection);
      3. find first-shell bath atoms and surface reference atoms;
      4. redirect those atoms' full-step velocities, preserving speed;
      5. for a leapfrog integrator, convert every velocity to the half-step value
         v(t - dt/2) = v(t) - F/m dt/2  (kRV.h does this for the resample; the
         redirect kernel does +F/m dt/2, redirect, -F/m dt/2, which composes to the same).
    """
    v_full = maxwell_boltzmann(rng, masses, cfg.temperature)
    directions = group_directions(rng, cfg.link, cfg.rtype, cfg.membrane_axis)
    shell = find_first_shell(positions, box, cfg, method=method)
    atoms = np.concatenate([shell.bath_atoms, shell.surface_atoms])
    grp = np.concatenate([shell.bath_group, shell.surface_group])
    redirect(v_full, atoms, grp, directions, cfg.rtype, cfg.membrane_axis)
    if half_step_velocities:
        m = masses > 0
        v_full[m] -= forces[m] / masses[m, None] * (0.5 * dt)
    return v_full, shell, directions
