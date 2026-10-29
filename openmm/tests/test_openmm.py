import numpy as np
import pytest

mm = pytest.importorskip("openmm")
from openmm import app, unit  # noqa: E402

from kamd_openmm import KAMDConfig, KAMDReporter  # noqa: E402


def _water_box(seed=0):
    ff = app.ForceField("amber14/tip3p.xml")
    mod = app.Modeller(app.Topology(), [])
    mod.addSolvent(ff, boxSize=mm.Vec3(2.4, 2.4, 2.4) * unit.nanometer)
    system = ff.createSystem(
        mod.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.8 * unit.nanometer, constraints=app.HBonds
    )
    integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 5 / unit.picosecond, 0.002 * unit.picosecond)
    integ.setRandomNumberSeed(seed)
    sim = app.Simulation(mod.topology, system, integ, mm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(mod.positions)
    sim.minimizeEnergy(maxIterations=200)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)
    st = sim.context.getState(getPositions=True)
    L = np.diag(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    pos = np.mod(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer), L)
    # "solutes": the waters whose oxygen is in two opposite corners
    res = list(mod.topology.residues())
    g1, g2, solvent = [], [], []
    for r in res:
        idx = [a.index for a in r.atoms()]
        o = pos[idx[0]]
        if np.linalg.norm(o - 0.25 * L) < 0.45:
            g1 += idx
        elif np.linalg.norm(o - 0.75 * L) < 0.45:
            g2 += idx
        else:
            solvent += idx
    assert len(g1) >= 15 and len(g2) >= 15
    return sim, g1, g2, solvent


def test_reporter_fires_and_redirects():
    sim, g1, g2, solvent = _water_box()
    cfg = KAMDConfig(
        refer_groups=[g1, g2], solvent_atoms=solvent, temperature=300.0, pseudobd_cycles=10, nzero_crossings=1, seed=1
    )
    rep = KAMDReporter(cfg, log_checks=True)
    sim.reporters.append(rep)
    sim.step(1999)
    assert rep.n_checks == 200  # checks at completed steps 9, 19, ..., 1999
    assert rep.n_resets > 0
    for ev in rep.events:
        assert (ev["step"] + 1) % 10 == 0
    T = sim.context.getState(getEnergy=True).getKineticEnergy()
    assert np.isfinite(T.value_in_unit(unit.kilojoule_per_mole))


def test_reset_velocities_on_device_are_redirected():
    sim, g1, g2, solvent = _water_box(seed=3)
    sim.step(50)
    cfg = KAMDConfig(refer_groups=[g1, g2], solvent_atoms=solvent, temperature=300.0, seed=5)
    rep = KAMDReporter(cfg)
    rep.describeNextReport(sim)  # performs setup
    rep._do_reset(sim, sim.currentStep)
    st = sim.context.getState(getVelocities=True, getForces=True)
    v = st.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
    F = st.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    m = rep._masses
    vfull = v + F / m[:, None] * 0.001
    sh = rep.last_shell
    dirs = np.array(rep.events[-1]["directions"])
    atoms = np.concatenate([sh.bath_atoms, sh.surface_atoms])
    grp = np.concatenate([sh.bath_group, sh.surface_group])
    assert atoms.size > 0
    u = vfull[atoms] / np.linalg.norm(vfull[atoms], axis=1)[:, None]
    cos = (u * dirs[grp]).sum(1)
    assert cos.min() > 1 - 1e-5  # CPU platform stores velocities in single precision
    # continuing the run is stable
    sim.step(200)
    ke = sim.context.getState(getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    ndof = 3 * sim.system.getNumParticles() - sim.system.getNumConstraints() - 3
    T = 2 * ke / (ndof * 0.00831446261815324)
    assert 250 < T < 350


def test_amber_style_input(tmp_path):
    from kamd_openmm.amber_input import parse_namelist

    txt = (
        "&cntrl temp0=310, /\n&kamd pseudobd_cycles=100, refer_mask(1)=':1', refer_mask(2)=':2',"
        " solvent_mask=':WAT', /\n"
    )
    assert parse_namelist(txt, "kamd")["refer_mask"] == {1: ":1", 2: ":2"}
