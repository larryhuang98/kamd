"""OpenMM integration: KAMD as a Simulation reporter.

Usage::

    from kamd_openmm import KAMDConfig, KAMDReporter
    cfg = KAMDConfig(refer_groups=[g1, g2], solvent_atoms=water, temperature=300)
    simulation.reporters.append(KAMDReporter(cfg, log="kamd.log"))
    simulation.step(1_000_000)

The reporter is called every ``pseudobd_cycles`` steps.  It reads velocities,
updates the zero-crossing monitor, and when a reset fires it reads positions and
forces, builds the new velocities on the host and writes them back with
``Context.setVelocities``.  Host transfers happen only at those points.
"""

from __future__ import annotations

import json
import time

import numpy as np
from openmm import unit

from .core import KAMDConfig, MomentumMonitor, group_momenta, pseudobd_reset

# Integrators whose stored velocities lag the positions by dt/2 (leapfrog family).
_HALF_STEP = {
    "VerletIntegrator",
    "LangevinMiddleIntegrator",
    "LangevinIntegrator",
    "VariableVerletIntegrator",
    "VariableLangevinIntegrator",
    "DrudeLangevinIntegrator",
}
# Integrators whose stored velocities are synchronous with positions (velocity Verlet).
_FULL_STEP = {"NoseHooverIntegrator", "DrudeNoseHooverIntegrator"}

_NM_PS = unit.nanometer / unit.picosecond
_KJ_NM = unit.kilojoule_per_mole / unit.nanometer


def velocity_offset_for(integrator) -> str:
    name = type(integrator).__name__
    if name in _HALF_STEP:
        return "half"
    if name in _FULL_STEP:
        return "full"
    if name == "BrownianIntegrator":
        raise ValueError("KAMD needs velocities; BrownianIntegrator has none")
    raise ValueError(
        f"Don't know the velocity convention of {name}. Pass velocity_offset='half' "
        "(velocities stored at t-dt/2, leapfrog) or 'full' (velocities at t)."
    )


class KAMDReporter:
    """Kinetically Accelerated MD for OpenMM.

    Parameters
    ----------
    config : KAMDConfig
    log : path or file-like, optional
        Event log (one line per reset; per-check lines too if ``log_checks``).
    velocity_offset : 'auto' | 'half' | 'full'
        Whether the integrator stores half-step (leapfrog) or full-step velocities.
    search : 'auto' | 'kdtree' | 'bruteforce'
        First-shell search backend.  'auto' uses a periodic k-d tree for
        rectangular boxes and a brute-force minimum-image search otherwise.
    apply_velocity_constraints : bool
        Call ``Context.applyVelocityConstraints`` after a reset.  pmemd does not
        (SHAKE removes the constraint components on the next step); default False.
    log_checks : bool
        Also log every monitoring check (projections and crossing counters).
    """

    def __init__(
        self,
        config: KAMDConfig,
        log=None,
        velocity_offset: str = "auto",
        search: str = "auto",
        apply_velocity_constraints: bool = False,
        log_checks: bool = False,
    ):
        self.cfg = config
        self.velocity_offset = velocity_offset
        self.search = search
        self.apply_velocity_constraints = apply_velocity_constraints
        self.log_checks = log_checks
        self.rng = np.random.default_rng(config.seed)
        self.events: list = []
        self.n_checks = 0
        self.n_resets = 0
        self.reset_time = 0.0
        self._start: int | None = None
        self._masses: np.ndarray | None = None
        self._monitor: MomentumMonitor | None = None
        self._half: bool | None = None
        self._own_log = isinstance(log, str)
        self._log = open(log, "w") if isinstance(log, str) else log

    # ---------------------------------------------------------------- setup
    def _setup(self, simulation) -> None:
        system = simulation.system
        n = system.getNumParticles()
        self.cfg.validate(n_atoms=n)
        self._masses = np.array([system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(n)])
        gm = np.array([self._masses[g].sum() for g in self.cfg.groups])
        self._monitor = MomentumMonitor(gm, self.cfg.rtype, self.cfg.membrane_axis, self.cfg.nzero_crossings)
        off = velocity_offset_for(simulation.integrator) if self.velocity_offset == "auto" else self.velocity_offset
        if off not in ("half", "full"):
            raise ValueError("velocity_offset must be 'auto', 'half' or 'full'")
        self._half = off == "half"
        self._start = simulation.currentStep
        if self.cfg.ngroups == 1:
            self._write(
                "# WARNING: one reference group -> its momentum relative to the reference COM is "
                "zero up to round-off, so zero crossings (and resets) are driven by "
                "floating-point noise.  pmemd uses the same formula.  Use >= 2 groups."
            )
        self._write(
            f"# KAMD for OpenMM: {self.cfg.ngroups} groups, sizes "
            f"{[int(g.size) for g in self.cfg.groups]}, solvent {self.cfg.solvent.size} atoms, "
            f"membrane {self.cfg.membrane.size} atoms, pseudobd_cycles {self.cfg.pseudobd_cycles}, "
            f"velocities {'half-step' if self._half else 'full-step'}"
        )
        self._write("# step\tevent\tshell_atoms_per_group\tsurface_atoms_per_group\tdirections")

    def _write(self, line: str) -> None:
        if self._log is not None:
            self._log.write(line + "\n")
            self._log.flush()

    # ------------------------------------------------------- reporter API
    def describeNextReport(self, simulation):
        if self._start is None:
            self._setup(simulation)
        n = simulation.currentStep - self._start  # completed steps since KAMD started
        steps = self.cfg.pseudobd_cycles - ((n + 1) % self.cfg.pseudobd_cycles)
        if steps == 0:
            steps = self.cfg.pseudobd_cycles
        # (steps, positions, velocities, forces, energies)
        return (steps, False, True, False, False)

    def report(self, simulation, state) -> None:
        v = state.getVelocities(asNumpy=True).value_in_unit(_NM_PS)
        p = group_momenta(v, self._masses, self.cfg.groups)
        fire = self._monitor.update(p)
        self.n_checks += 1
        step = simulation.currentStep
        if self.log_checks:
            self._write(
                f"{step}\tcheck\tproj={np.array2string(self._monitor.last, precision=4)}\t"
                f"count={self._monitor.count.tolist()}"
            )
        if fire:
            self._do_reset(simulation, step)

    # ------------------------------------------------------------ reset
    def _do_reset(self, simulation, step: int) -> None:
        t0 = time.perf_counter()
        ctx = simulation.context
        st = ctx.getState(getPositions=True, getForces=True)
        pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        frc = st.getForces(asNumpy=True).value_in_unit(_KJ_NM)
        box = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        dt = simulation.integrator.getStepSize().value_in_unit(unit.picosecond)
        newv, shell, dirs = pseudobd_reset(
            self.rng,
            np.asarray(pos),
            np.asarray(frc),
            self._masses,
            np.asarray(box),
            self.cfg,
            dt,
            half_step_velocities=self._half,
            method=self.search,
        )
        ctx.setVelocities(newv * _NM_PS)
        if self.apply_velocity_constraints:
            ctx.applyVelocityConstraints(1e-5)
        self._monitor.reset()
        self.n_resets += 1
        self.reset_time += time.perf_counter() - t0
        nb, ns = shell.counts(self.cfg.ngroups)
        ev = {"step": int(step), "shell": nb.tolist(), "surface": ns.tolist(), "directions": dirs.round(6).tolist()}
        self.events.append(ev)
        self._write(f"{step}\treset\t{ev['shell']}\t{ev['surface']}\t{json.dumps(ev['directions'])}")
        self.last_shell = shell

    # ---------------------------------------------------- checkpointing
    def state_dict(self) -> dict:
        """Monitor counters + RNG state, for continuing a run from a checkpoint."""
        return {
            "monitor": self._monitor.state_dict() if self._monitor else None,
            "rng": self.rng.bit_generator.state,
            "n_resets": self.n_resets,
            "n_checks": self.n_checks,
        }

    def load_state_dict(self, d: dict, simulation) -> None:
        """Call after loading the OpenMM checkpoint (sets the phase to the current step)."""
        if self._start is None:
            self._setup(simulation)
        if d.get("monitor"):
            self._monitor.load_state_dict(d["monitor"])
        self.rng.bit_generator.state = d["rng"]
        self.n_resets = d.get("n_resets", 0)
        self.n_checks = d.get("n_checks", 0)

    def close(self) -> None:
        if self._own_log and self._log is not None:
            self._log.close()
            self._log = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
