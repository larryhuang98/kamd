"""Toy system for testing KAMD canonical sampling.

Argon-like LJ liquid (sigma 3.4 A, eps 0.996 kJ/mol, m 39.9, T 120 K, rho* 0.8) in a
cubic periodic box, plus two 13-particle icosahedral LJ clusters held together by
stiff harmonic bonds (12 spokes + 30 edges). A quartic double well acts on the
COM-COM distance:   U(r) = h [((r - r0)/a)^2 - 1]^2,   wells at r0 +- a.

Construction: an FCC lattice of 4 n^3 sites; two 13-site cuboctahedra (a site +
its 12 nearest neighbours) are replaced by icosahedra, so the particle count and
density are exactly those of the lattice.
"""

from __future__ import annotations

import numpy as np
import openmm as mm
from openmm import unit

KB = 0.00831446261815324  # kJ/mol/K

SIGMA = 0.34  # nm
EPS = 0.996  # kJ/mol
MASS = 39.9
TEMP = 120.0
RHO_STAR = 0.8
CUTOFF = 2.5 * SIGMA

ICO_R = 0.37  # nm, centre-vertex distance (edge = 1.0515 R = 0.389 nm ~ LJ minimum)
K_BOND = 20000.0  # kJ/mol/nm^2

# double well on COM-COM distance (nm)
DW_R0 = 1.50
DW_A = 0.25
DW_H_KT = 2.0  # barrier of the bare double well, in kT


def icosahedron(R: float) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """13 positions (centre first) and the 42 bonds (12 spokes + 30 edges)."""
    phi = (1 + 5**0.5) / 2
    v = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            v += [(0, s1, s2 * phi), (s1, s2 * phi, 0), (s2 * phi, 0, s1)]
    v = np.array(v, float)
    v *= R / np.linalg.norm(v[0])
    pos = np.vstack([np.zeros(3), v])
    d = np.linalg.norm(v[:, None] - v[None], axis=-1)
    edge = d[d > 1e-9].min()
    bonds = [(0, i + 1) for i in range(12)]
    bonds += [(i + 1, j + 1) for i in range(12) for j in range(i + 1, 12) if abs(d[i, j] - edge) < 1e-6]
    assert len(bonds) == 42
    return pos, bonds


def build(n_cells: int = 8, h_kT: float = DW_H_KT, r0: float = DW_R0, a_dw: float = DW_A, start: str = "outer", seed: int = 0):
    """Return (system, positions[nm], info dict)."""
    rng = np.random.default_rng(seed)
    N = 4 * n_cells**3
    rho = RHO_STAR / (SIGMA**3)  # nm^-3
    L = (N / rho) ** (1 / 3)
    a = L / n_cells
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    cells = np.array([[i, j, k] for i in range(n_cells) for j in range(n_cells) for k in range(n_cells)])
    lat = ((cells[:, None, :] + basis[None]) * a).reshape(-1, 3)

    def cubo(center_site):
        d = lat - lat[center_site]
        d -= L * np.round(d / L)
        r = np.linalg.norm(d, axis=1)
        return np.argsort(r)[:13]

    # two cluster centres separated by ~ the outer (or inner) well distance along x
    c1 = np.argmin(np.linalg.norm(lat - np.array([0.25 * L, 0.5 * L, 0.5 * L]), axis=1))
    c2 = np.argmin(np.linalg.norm(lat - (lat[c1] + np.array([r0 + (a_dw if start == "outer" else -a_dw), 0, 0])), axis=1))
    s1, s2 = cubo(c1), cubo(c2)
    assert not set(s1) & set(s2)
    removed = np.concatenate([s1, s2])
    solvent_pos = np.delete(lat, removed, axis=0)

    ico, bonds = icosahedron(ICO_R)
    pos = np.vstack([lat[c1] + ico, lat[c2] + ico, solvent_pos])
    g1 = np.arange(13)
    g2 = np.arange(13, 26)
    solvent = np.arange(26, N)

    system = mm.System()
    for _ in range(N):
        system.addParticle(MASS)
    system.setDefaultPeriodicBoxVectors(mm.Vec3(L, 0, 0), mm.Vec3(0, L, 0), mm.Vec3(0, 0, L))

    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(CUTOFF)
    nb.setUseDispersionCorrection(True)
    nb.setUseSwitchingFunction(False)
    for _ in range(N):
        nb.addParticle(0.0, SIGMA, EPS)
    # all intra-cluster pairs excluded: clusters are held by the bonds only
    for g in (g1, g2):
        for i in range(13):
            for j in range(i + 1, 13):
                nb.addException(int(g[i]), int(g[j]), 0.0, SIGMA, 0.0)
    system.addForce(nb)

    hb = mm.HarmonicBondForce()
    hb.setUsesPeriodicBoundaryConditions(True)
    for g in (g1, g2):
        for i, j in bonds:
            b0 = float(np.linalg.norm(ico[i] - ico[j]))
            hb.addBond(int(g[i]), int(g[j]), b0, K_BOND)
    system.addForce(hb)

    kT = KB * TEMP
    dw = mm.CustomCentroidBondForce(2, "h*(((distance(g1,g2)-r0)/a)^2-1)^2")
    dw.addGlobalParameter("h", h_kT * kT)
    dw.addGlobalParameter("r0", r0)
    dw.addGlobalParameter("a", a_dw)
    dw.addGroup([int(i) for i in g1])
    dw.addGroup([int(i) for i in g2])
    dw.addBond([0, 1], [])
    dw.setUsesPeriodicBoundaryConditions(True)
    system.addForce(dw)

    info = dict(N=N, L=L, g1=g1, g2=g2, solvent=solvent, kT=kT, bonds=bonds)
    return system, pos, info


def make_simulation(system, pos, friction=1.0, dt_fs=5.0, seed=1, threads=1):
    integ = mm.LangevinMiddleIntegrator(TEMP * unit.kelvin, friction / unit.picosecond, dt_fs * unit.femtosecond)
    integ.setRandomNumberSeed(int(seed))
    plat = mm.Platform.getPlatformByName("CPU")
    from openmm import app

    top = app.Topology()
    ch = top.addChain()
    el = app.Element.getBySymbol("Ar")
    for i in range(system.getNumParticles()):
        top.addAtom("AR", el, top.addResidue("AR", ch))
    sim = app.Simulation(top, system, integ, plat, {"Threads": str(threads)})
    sim.context.setPositions(pos)
    return sim


def com_distance(x, g1, g2, L):
    """Minimum-image COM-COM distance; clusters are made whole relative to atom 0 of each."""
    def com(g):
        y = x[g]
        y = y[0] + (y - y[0]) - L * np.round((y - y[0]) / L)
        return y.mean(axis=0), y
    c1, y1 = com(g1)
    c2, y2 = com(g2)
    d = c2 - c1
    d -= L * np.round(d / L)
    return np.linalg.norm(d), c1, c2, y1, y2
