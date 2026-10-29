"""Run an Amber system (prmtop + restart + mdin) in OpenMM with or without KAMD.

The &cntrl settings used by the KAMD test (ntt=3, ntp=1, ntc=2, cut, dt, temp0,
gamma_ln) are mapped to OpenMM; the &kamd namelist is read by kamd_openmm.
Writes <out>.json with timing, reset steps and thermodynamic averages.
"""

import argparse
import json
import time

import numpy as np
import openmm as mm
from openmm import app, unit

from kamd_openmm import KAMDReporter
from kamd_openmm.amber_input import config_from_amber, parse_namelist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--restart", required=True)
    ap.add_argument("--mdin", required=True)
    ap.add_argument("--nsteps", type=int, default=250000)
    ap.add_argument("--platform", default="CUDA")
    ap.add_argument("--precision", default="mixed")
    ap.add_argument("--no-kamd", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pme-alpha", type=float, help="Ewald coefficient in 1/Angstrom (Amber mdout 'Ewald Coefficient')")
    ap.add_argument("--pme-grid", type=int, nargs=3, help="PME grid (Amber NFFT1-3)")
    ap.add_argument("--cmm-freq", type=int, help="CMMotionRemover frequency (Amber nscm)")
    a = ap.parse_args()

    text = open(a.mdin).read()
    c = parse_namelist(text, "cntrl")
    T = float(c.get("temp0", 300.0))
    dt = float(c.get("dt", 0.002))
    gamma = float(c.get("gamma_ln", 1.0))
    cut = float(c.get("cut", 8.0))

    prm = app.AmberPrmtopFile(a.prmtop)
    crd = app.AmberInpcrdFile(a.restart)
    system = prm.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=cut * unit.angstrom,
        constraints=app.HBonds if int(c.get("ntc", 1)) == 2 else None,
        rigidWater=True,
    )
    if a.pme_alpha or a.pme_grid:
        nb = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
        nb.setPMEParameters(a.pme_alpha * 10.0, *a.pme_grid)  # 1/A -> 1/nm
    if a.cmm_freq:
        for f in system.getForces():
            if isinstance(f, mm.CMMotionRemover):
                f.setFrequency(a.cmm_freq)
    if int(c.get("ntp", 0)) > 0:
        if int(c.get("barostat", 1)) != 2:
            print("note: Amber barostat=1 (Berendsen) has no OpenMM equivalent; using Monte Carlo")
        system.addForce(
            mm.MonteCarloBarostat(float(c.get("pres0", 1.0)) * unit.bar, T * unit.kelvin, int(c.get("mcbarint", 100)))
        )
    integ = mm.LangevinMiddleIntegrator(T * unit.kelvin, gamma / unit.picosecond, dt * unit.picosecond)
    integ.setRandomNumberSeed(a.seed)
    plat = mm.Platform.getPlatformByName(a.platform)
    props = {"Precision": a.precision} if a.platform in ("CUDA", "OpenCL") else {}
    sim = app.Simulation(prm.topology, system, integ, plat, props)
    sim.context.setPeriodicBoxVectors(*crd.getBoxVectors())
    sim.context.setPositions(crd.positions)
    if crd.velocities is not None:
        sim.context.setVelocities(crd.velocities)
    else:
        sim.context.setVelocitiesToTemperature(T * unit.kelvin, a.seed)

    sim.reporters.append(
        app.StateDataReporter(
            a.out + ".thermo",
            max(1, min(2500, a.nsteps // 10)),
            step=True,
            temperature=True,
            potentialEnergy=True,
            density=True,
            speed=True,
        )
    )
    rep = None
    if not a.no_kamd:
        cfg = config_from_amber(text, a.prmtop, seed=a.seed)
        rep = KAMDReporter(cfg, log=a.out + ".kamd.log")
        sim.reporters.append(rep)

    sim.step(10)  # warm-up / kernel compilation, excluded from timing
    t0 = time.perf_counter()
    sim.step(a.nsteps)
    wall = time.perf_counter() - t0

    th = np.genfromtxt(a.out + ".thermo", delimiter=",", names=True)
    summary = {
        "platform": a.platform,
        "precision": a.precision,
        "kamd": rep is not None,
        "natoms": system.getNumParticles(),
        "nsteps": a.nsteps,
        "wall_s": wall,
        "ns_per_day": a.nsteps * dt / 1000.0 / (wall / 86400.0),
        "T_mean": float(np.mean(th[th.dtype.names[2]])),
        "density_mean": float(np.mean(th[th.dtype.names[3]])),
    }
    if rep is not None:
        # KAMD starts counting at step 0 (during warm-up); +1 gives pmemd's "step S" numbering
        steps = [e["step"] + 1 for e in rep.events]
        summary.update(
            {
                "n_checks": rep.n_checks,
                "n_resets": rep.n_resets,
                "reset_steps": steps,
                "reset_intervals": np.diff([0] + steps).tolist(),
                "reset_wall_s": rep.reset_time,
                "mean_shell_atoms": float(np.mean([sum(e["shell"]) for e in rep.events])) if rep.events else 0,
                "mean_surface_atoms": float(np.mean([sum(e["surface"]) for e in rep.events])) if rep.events else 0,
            }
        )
        rep.close()
    json.dump(summary, open(a.out + ".json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("reset_steps",)}, indent=1))


if __name__ == "__main__":
    main()
