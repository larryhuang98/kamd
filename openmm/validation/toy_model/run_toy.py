"""Run one replica of the KAMD toy model and store compact observables (npz).

    python run_toy.py --mode plain|kamd|mbonly|sepmon --seed 1 --ns 4 --out runs/x.npz
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from openmm import unit
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import toy_system as ts  # noqa: E402

GR_RMAX = 1.2  # nm
GR_NB = 240
SG_RMAX = 1.5
SG_NB = 300
U_MAX = 6.0
U_NB = 240
TSR_EDGES_PS = np.array([0.0, 0.25, 1.0, 3.0, 1e9])  # time-since-reset classes


class Sampler:
    """Collects observables. Must be appended AFTER the KAMD reporter."""

    def __init__(self, info, dt_ps, kamd=None, every=50, vel_every=100, gr_every=1000, sg_every=200):
        self.info, self.dt, self.kamd = info, dt_ps, kamd
        self.every, self.vel_every, self.gr_every, self.sg_every = every, vel_every, gr_every, sg_every
        self.g1, self.g2, self.solv = info["g1"], info["g2"], info["solvent"]
        self.L = info["L"]
        self.r, self.rg, self.step = [], [], []
        self.tlog = []  # step, tsr_ps, T_all(sync), T_solv, T_solute, T_all(half), KEcom1/kT, KEcom2/kT
        nc = len(TSR_EDGES_PS) - 1
        self.uhist = np.zeros((nc, 2, U_NB))  # class, (solvent, solute), bins
        self.gr_ss = np.zeros(GR_NB)
        self.n_gr = 0
        self.gr_us = np.zeros((2, SG_NB))  # (solute atom-solvent, COM-solvent)
        self.n_sg = 0
        self.start = None

    def describeNextReport(self, sim):
        s = sim.currentStep
        steps = self.every - s % self.every
        need_v = (s + steps) % self.vel_every == 0
        return (steps, True, need_v, need_v, False)

    def _tsr(self, step):
        if self.kamd is None or not self.kamd.events:
            return np.inf if self.kamd is None else (step - self.start) * self.dt
        return (step - self.kamd.events[-1]["step"]) * self.dt

    def report(self, sim, state):
        step = sim.currentStep
        if self.start is None:
            self.start = step
        x = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        r, c1, c2, y1, y2 = ts.com_distance(x, self.g1, self.g2, self.L)
        rg1 = np.sqrt(((y1 - c1) ** 2).sum(1).mean())
        rg2 = np.sqrt(((y2 - c2) ** 2).sum(1).mean())
        self.r.append(r)
        self.rg.append((rg1, rg2))
        self.step.append(step)
        if step % self.vel_every == 0:
            vh = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
            f = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            v = vh + 0.5 * self.dt * f / ts.MASS
            kT = self.info["kT"]
            u = v / np.sqrt(kT / ts.MASS)
            tsr = self._tsr(step)
            c = min(int(np.searchsorted(TSR_EDGES_PS, tsr, side="right") - 1), len(TSR_EDGES_PS) - 2)
            sol = np.concatenate([self.g1, self.g2])
            for k, idx in enumerate((self.solv, sol)):
                h, _ = np.histogram(u[idx].ravel(), bins=U_NB, range=(-U_MAX, U_MAX))
                self.uhist[c, k] += h
            T = lambda w, idx: ts.MASS * (w[idx] ** 2).sum() / (3 * len(idx) * ts.KB)  # noqa: E731
            allidx = np.arange(len(v))
            vc1 = v[self.g1].mean(0)
            vc2 = v[self.g2].mean(0)
            kec = lambda vc: 0.5 * 13 * ts.MASS * (vc**2).sum() / kT  # noqa: E731
            self.tlog.append(
                (step, tsr, T(v, allidx), T(v, self.solv), T(v, sol), T(vh, allidx), kec(vc1), kec(vc2))
            )
        if step % self.sg_every == 0:
            xs = x[self.solv]
            sol = np.concatenate([self.g1, self.g2])
            d = x[sol][:, None, :] - xs[None]
            d -= self.L * np.round(d / self.L)
            h, _ = np.histogram(np.linalg.norm(d, axis=-1), bins=SG_NB, range=(0, SG_RMAX))
            self.gr_us[0] += h
            d = np.stack([c1, c2])[:, None, :] - xs[None]
            d -= self.L * np.round(d / self.L)
            h, _ = np.histogram(np.linalg.norm(d, axis=-1), bins=SG_NB, range=(0, SG_RMAX))
            self.gr_us[1] += h
            self.n_sg += 1
        if step % self.gr_every == 0:
            xs = np.mod(x[self.solv], self.L)
            tree = cKDTree(xs, boxsize=self.L)
            pairs = tree.query_pairs(GR_RMAX, output_type="ndarray")
            d = xs[pairs[:, 0]] - xs[pairs[:, 1]]
            d -= self.L * np.round(d / self.L)
            h, _ = np.histogram(np.linalg.norm(d, axis=1), bins=GR_NB, range=(0, GR_RMAX))
            self.gr_ss += 2 * h
            self.n_gr += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="plain", choices=["plain", "kamd", "mbonly", "sepmon"])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ns", type=float, default=1.0)
    ap.add_argument("--equil_ps", type=float, default=100.0)
    ap.add_argument("--friction", type=float, default=1.0)
    ap.add_argument("--dt_fs", type=float, default=5.0)
    ap.add_argument("--pbd", type=int, default=250)
    ap.add_argument("--nz", type=int, default=2)
    ap.add_argument("--shell_cut", type=float, default=5.3)
    ap.add_argument("--surface_cut", type=float, default=4.0)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--dw_h", type=float, default=ts.DW_H_KT)
    ap.add_argument("--dw_r0", type=float, default=ts.DW_R0)
    ap.add_argument("--dw_a", type=float, default=ts.DW_A)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    start = "outer" if a.seed % 2 else "inner"  # half of the replicas start in each well
    system, pos, info = ts.build(h_kT=a.dw_h, r0=a.dw_r0, a_dw=a.dw_a, start=start)
    sim = ts.make_simulation(system, pos, friction=a.friction, dt_fs=a.dt_fs, seed=1000 + a.seed, threads=a.threads)
    sim.minimizeEnergy(maxIterations=1000)
    sim.context.setVelocitiesToTemperature(ts.TEMP * unit.kelvin, 7919 * a.seed + 13)
    dt = a.dt_fs * 1e-3
    t0 = time.time()
    # equilibration is ALWAYS plain Langevin (gamma = 1/ps), then the production condition
    sim.integrator.setFriction(1.0)
    sim.step(int(a.equil_ps / dt))
    sim.integrator.setFriction(a.friction)
    sim.currentStep = 0

    kamd = None
    if a.mode != "plain":
        from kamd_openmm import KAMDConfig, KAMDReporter

        from kamd_variants import SeparationMonitorReporter

        sc, sf = (a.shell_cut, a.surface_cut) if a.mode != "mbonly" else (1e-3, 1e-3)
        cfg = KAMDConfig(
            refer_groups=[info["g1"], info["g2"]],
            solvent_atoms=info["solvent"],
            temperature=ts.TEMP,
            pseudobd_cycles=a.pbd,
            nzero_crossings=a.nz,
            shell_cut=sc,
            surface_cut=sf,
            seed=99991 * a.seed + 7,
        )
        cls = SeparationMonitorReporter if a.mode == "sepmon" else KAMDReporter
        kamd = cls(cfg)
        sim.reporters.append(kamd)
    smp = Sampler(info, dt, kamd)
    sim.reporters.append(smp)
    nsteps = int(round(a.ns * 1e3 / dt))
    t1 = time.time()
    sim.step(nsteps)
    t2 = time.time()

    ev = kamd.events if kamd else []
    meta = dict(vars(a), N=info["N"], L=info["L"], kT=info["kT"], equil_wall=t1 - t0, prod_wall=t2 - t1,
                steps_per_s=nsteps / (t2 - t1), n_resets=len(ev), n_checks=kamd.n_checks if kamd else 0,
                reset_time=kamd.reset_time if kamd else 0.0)
    np.savez_compressed(
        a.out,
        r=np.asarray(smp.r, np.float32),
        rg=np.asarray(smp.rg, np.float32),
        step=np.asarray(smp.step, np.int64),
        tlog=np.asarray(smp.tlog, np.float64),
        uhist=smp.uhist,
        gr_ss=smp.gr_ss,
        n_gr=smp.n_gr,
        gr_us=smp.gr_us,
        n_sg=smp.n_sg,
        reset_steps=np.array([e["step"] for e in ev], np.int64),
        reset_shell=np.array([e["shell"] for e in ev], np.int64).reshape(-1, 2),
        reset_surface=np.array([e["surface"] for e in ev], np.int64).reshape(-1, 2),
        tsr_edges=TSR_EDGES_PS,
        meta=json.dumps(meta),
    )
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
