"""dt stability check: NVE (Verlet) energy drift after Langevin equilibration, dt = 2, 4, 5 fs."""
import sys
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit

sys.path.insert(0, str(Path(__file__).resolve().parent))
import toy_system as ts

system, pos, info = ts.build()
sim = ts.make_simulation(system, pos, dt_fs=4.0, seed=5)
sim.minimizeEnergy(maxIterations=1000)
sim.context.setVelocitiesToTemperature(ts.TEMP * unit.kelvin, 3)
sim.step(10000)
st = sim.context.getState(getPositions=True, getVelocities=True)
x0, v0 = st.getPositions(), st.getVelocities()
for dt in (2.0, 4.0, 5.0):
    integ = mm.VerletIntegrator(dt * unit.femtosecond)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("CPU"), {"Threads": "4"})
    ctx.setPositions(x0)
    ctx.setVelocities(v0)
    nstep = int(20000 / dt)  # 20 ps
    E = []
    for k in range(40):
        integ.step(nstep // 40)
        s = ctx.getState(getEnergy=True)
        E.append((s.getPotentialEnergy() + s.getKineticEnergy()).value_in_unit(unit.kilojoule_per_mole))
    E = np.array(E)
    t = np.arange(1, 41) * 0.5
    slope = np.polyfit(t, E, 1)[0]
    print(f"dt {dt} fs: NVE drift {slope / info['N'] / info['kT'] * 1000:.3e} kT/particle/ns, "
          f"E fluct sd {E.std() / info['N'] / info['kT']:.2e} kT/particle")
