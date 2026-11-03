"""Experimental KAMD variants for the toy-model test. core.py is NOT modified.

SeparationMonitorReporter: same reset as KAMDReporter, but the zero-crossing
monitor projects each group's momentum (relative to the reference COM) onto the
unit vector from the reference COM to the group COM, i.e. onto the inter-group
separation for two groups, instead of onto the reference COM velocity.
"""

from __future__ import annotations

import numpy as np
from openmm import unit

from kamd_openmm.core import group_momenta
from kamd_openmm.reporter import _NM_PS, KAMDReporter


class SeparationMonitorReporter(KAMDReporter):
    def describeNextReport(self, simulation):
        steps, _, v, f, e = super().describeNextReport(simulation)
        return (steps, True, True, False, False)  # also positions

    def report(self, simulation, state) -> None:
        v = state.getVelocities(asNumpy=True).value_in_unit(_NM_PS)
        x = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        L = np.diag(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
        m = self._masses
        coms = []
        for g in self.cfg.groups:
            y = x[g]
            y = y[0] + (y - y[0]) - L * np.round((y - y[0]) / L)
            coms.append((m[g, None] * y).sum(0) / m[g].sum())
        coms = np.array(coms)
        # reference COM, with all groups imaged next to group 0
        d = coms - coms[0]
        coms = coms[0] + d - L * np.round(d / L)
        gm = self._monitor.group_masses
        rcom = (gm[:, None] * coms).sum(0) / gm.sum()
        u = coms - rcom
        u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-12)
        p = group_momenta(v, m, self.cfg.groups)
        vcom = p.sum(0) / gm.sum()
        rel = p - gm[:, None] * vcom
        s = np.einsum("ij,ij->i", rel, u)
        mon = self._monitor
        crossed = (s * mon.last < 0.0) & (mon.count < mon.nzero)
        mon.count[crossed] += 1
        mon.last = s
        self.n_checks += 1
        if np.all(mon.count >= mon.nzero):
            self._do_reset(simulation, simulation.currentStep)
