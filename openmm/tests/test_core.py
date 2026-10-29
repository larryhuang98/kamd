import numpy as np
import pytest

from kamd_openmm.amber_input import parse_namelist
from kamd_openmm.core import (
    KB_KJ_PER_MOL_K,
    KAMDConfig,
    KAMDInputError,
    MomentumMonitor,
    find_first_shell,
    group_directions,
    maxwell_boltzmann,
    pseudobd_reset,
    redirect,
)

# ----------------------------------------------------------------------------- namelist
MDIN = """
&cntrl
   imin=0, irest=1, ntx=5, dt=0.002, temp0=355, ig=1060,
/
&kamd
   pseudobd_cycles = 250,     ! comment, with a / slash inside
   shell_cut       = 5.0d0,
   refer_mask(1) = "@20-96",
   refer_mask(2) = ':285-339 & !@H=',
   refer_link = 0, 1,
   region_type(2) = 2,
   solvent_mask  = ":WAT",
/
"""


def test_parse_namelist():
    k = parse_namelist(MDIN, "kamd")
    assert k["pseudobd_cycles"] == 250
    assert k["shell_cut"] == 5.0
    assert k["refer_mask"] == {1: "@20-96", 2: ":285-339 & !@H="}
    assert k["refer_link"] == {1: 0, 2: 1}
    assert k["region_type"] == {2: 2}
    assert k["solvent_mask"] == ":WAT"
    c = parse_namelist(MDIN, "cntrl")
    assert c["temp0"] == 355 and c["dt"] == 0.002


# ----------------------------------------------------------------------------- config
def _cfg(**kw):
    base = dict(refer_groups=[[0, 1], [2, 3]], solvent_atoms=range(4, 20), temperature=300.0)
    base.update(kw)
    return KAMDConfig(**base)


@pytest.mark.parametrize(
    "kw, msg",
    [
        (dict(refer_groups=[[0, 1], [1, 2]]), "must not overlap"),
        (dict(surface_cut=6.0), "surface_cut must be <= shell_cut"),
        (dict(refer_link=[0, 2]), "cannot reference its own group"),
        (dict(refer_link=[0, 3]), "out of range"),
        (dict(region_type=[1, 2]), "membrane_atoms is required"),
        (dict(region_type=[1, 2], membrane_atoms=[30, 31]), "membrane_axis"),
        (dict(solvent_atoms=[1, 5]), "solvent bath overlap"),
        (dict(region_type=[1, 3]), "invalid"),
    ],
)
def test_validation_errors(kw, msg):
    with pytest.raises(KAMDInputError, match=msg):
        _cfg(**kw).validate()


# ----------------------------------------------------------------------------- monitor
def test_monitor_crossings_and_fire():
    m = MomentumMonitor(np.array([1.0, 1.0]), np.array([1, 1]), 0, nzero_crossings=2)

    # COM velocity along +x; group 1 ahead/behind alternately -> projection sign flips
    def p(sign):
        return np.array([[1.0 + sign, 0, 0], [1.0 - sign, 0, 0]])

    assert not m.update(p(+1))  # first check: last == 0, no crossing
    assert m.count.tolist() == [0, 0]
    assert not m.update(p(-1))  # crossing 1 for both groups
    assert m.count.tolist() == [1, 1]
    assert not m.update(p(-1))  # same sign: nothing
    assert m.update(p(+1))  # crossing 2 -> fire
    m.reset()
    assert m.count.tolist() == [0, 0] and np.all(m.last == 0)


def test_monitor_counter_capped_and_needs_all_groups():
    m = MomentumMonitor(np.array([1.0, 1.0, 2.0]), np.array([1, 1, 1]), 0, nzero_crossings=2)
    rng = np.random.default_rng(0)
    fired = False
    for _ in range(200):
        fired = m.update(rng.normal(size=(3, 3)) + [1.0, 0, 0])
        assert np.all(m.count <= 2)
        if fired:
            break
    assert fired and np.all(m.count == 2)


def test_single_group_projection_is_roundoff():
    m = MomentumMonitor(np.array([3.0]), np.array([1]), 0, 2)
    rng = np.random.default_rng(1)
    for _ in range(100):
        p = rng.normal(size=(1, 3))
        assert abs(m.projections(p)[0]) < 1e-12 * (p**2).sum()


def test_monitor_tm_axis_zeroed():
    m = MomentumMonitor(np.array([1.0, 1.0]), np.array([2, 2]), 3, 2)
    p = np.array([[1.0, 0.0, 5.0], [1.0, 0.0, -3.0]])  # relative motion only along z
    assert np.allclose(m.projections(p), 0.0)


# ----------------------------------------------------------------------------- directions
def test_directions_link_and_tm():
    rng = np.random.default_rng(3)
    d = group_directions(rng, np.array([0, 1, 0, 3]), np.array([1, 1, 2, 2]), membrane_axis=3)
    assert np.allclose(np.linalg.norm(d, axis=1), 1.0)
    assert np.allclose(d[1], d[0])  # linked to group 1
    assert d[2, 2] == 0.0 and d[3, 2] == 0.0  # TM: in-plane
    assert np.allclose(d[3], d[2])  # linked TM group


# ----------------------------------------------------------------------------- first shell
def _random_system(seed=0, L=3.0, n=3000):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0, L, size=(n, 3))
    x[:40] = rng.normal([0.1, 1.5, 1.5], 0.15, size=(40, 3))  # group 1 straddles x=0 boundary
    x[40:80] = rng.normal([1.5, 1.5, 2.9], 0.15, size=(40, 3))  # group 2 straddles z=L boundary
    box = np.diag([L, L, L])
    cfg = KAMDConfig(
        refer_groups=[range(0, 40), range(40, 80)], solvent_atoms=range(80, n), temperature=300.0
    ).validate()
    return x, box, cfg


def test_first_shell_kdtree_matches_bruteforce():
    x, box, cfg = _random_system()
    a = find_first_shell(x, box, cfg, method="kdtree")
    b = find_first_shell(x, box, cfg, method="bruteforce")
    assert dict(zip(a.bath_atoms, a.bath_group, strict=False)) == dict(zip(b.bath_atoms, b.bath_group, strict=False))
    assert dict(zip(a.surface_atoms, a.surface_group, strict=False)) == dict(
        zip(b.surface_atoms, b.surface_group, strict=False)
    )
    assert a.bath_atoms.size > 0 and a.surface_atoms.size > 0


def test_first_shell_reference_semantics():
    """Independent O(N^2) oracle written exactly like the pmemd CPU diagnostic loop."""
    x, box, cfg = _random_system(seed=7)
    L = np.diag(box)
    shell2, surf2 = (cfg.shell_cut * 0.1) ** 2, (cfg.surface_cut * 0.1) ** 2
    rec = {}
    dist = {}
    surf = {}
    for g, atoms in enumerate(cfg.groups):
        for j in atoms:
            d = np.abs(x[cfg.solvent] - x[j])
            d = np.minimum(d, L - d)
            r2 = (d**2).sum(1)
            if np.any(r2 < surf2):
                surf[j] = g
            for i, r in zip(cfg.solvent[r2 < shell2], r2[r2 < shell2], strict=False):
                if r < dist.get(i, shell2):
                    dist[i], rec[i] = r, g
    res = find_first_shell(x, box, cfg)
    assert dict(zip(res.bath_atoms.tolist(), res.bath_group.tolist(), strict=False)) == rec
    assert dict(zip(res.surface_atoms.tolist(), res.surface_group.tolist(), strict=False)) == surf


def test_first_shell_triclinic_bruteforce_runs():
    x, box, cfg = _random_system(seed=2)
    tric = box.copy()
    tric[1, 0] = 0.3
    res = find_first_shell(x, tric, cfg)  # auto -> bruteforce
    assert res.bath_atoms.size > 0


def test_tm_uses_membrane_bath_and_cut():
    x = np.array(
        [
            [0.0, 0, 0],
            [1.0, 0, 0],  # group 1 (solution), group 2 (TM)
            [0.04, 0, 0],
            [1.04, 0, 0],  # solvent near g1, solvent near g2
            [1.02, 0, 0],
            [0.02, 0, 0.0],
        ]
    )  # membrane near g2, membrane near g1
    x[5, 1] = 0.02
    cfg = KAMDConfig(
        refer_groups=[[0], [1]],
        solvent_atoms=[2, 3],
        membrane_atoms=[4, 5],
        region_type=[1, 2],
        membrane_axis=3,
        temperature=300.0,
    ).validate()
    res = find_first_shell(x, np.diag([5.0, 5.0, 5.0]), cfg)
    got = dict(zip(res.bath_atoms.tolist(), res.bath_group.tolist(), strict=False))
    assert got == {2: 0, 4: 1}  # solvent only for solution group, membrane only for TM group


# ----------------------------------------------------------------------------- redirect / reset
def test_redirect_preserves_speed():
    rng = np.random.default_rng(4)
    v = rng.normal(size=(10, 3))
    v0 = v.copy()
    d = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    redirect(v, np.arange(10), np.array([0] * 5 + [1] * 5), d, np.array([1, 2]), membrane_axis=3)
    assert np.allclose(np.linalg.norm(v[:5], axis=1), np.linalg.norm(v0[:5], axis=1))
    assert np.allclose(v[:5, :2], 0.0)
    assert np.allclose(v[5:, 2], v0[5:, 2])  # axial kept
    assert np.allclose(np.abs(v[5:, 0]), np.linalg.norm(v0[5:, :2], axis=1))
    assert np.allclose(v[5:, 1], 0.0)


def test_maxwell_boltzmann_temperature():
    rng = np.random.default_rng(5)
    m = np.tile([16.0, 1.008, 1.008], 20000)
    v = maxwell_boltzmann(rng, m, 300.0)
    ke = 0.5 * (m[:, None] * v**2).sum()
    T = 2 * ke / (3 * m.size * KB_KJ_PER_MOL_K)
    assert abs(T - 300.0) < 3.0
    assert np.all(maxwell_boltzmann(rng, np.array([0.0, 1.0]), 300)[0] == 0)


def test_pseudobd_reset_half_step_algebra():
    x, box, cfg = _random_system(seed=9)
    rng = np.random.default_rng(11)
    m = rng.uniform(1, 16, size=len(x))
    F = rng.normal(scale=500, size=x.shape)
    dt = 0.002
    v, shell, dirs = pseudobd_reset(np.random.default_rng(12), x, F, m, box, cfg, dt)
    vfull = v + F / m[:, None] * dt / 2
    atoms = np.concatenate([shell.bath_atoms, shell.surface_atoms])
    grp = np.concatenate([shell.bath_group, shell.surface_group])
    u = vfull[atoms] / np.linalg.norm(vfull[atoms], axis=1)[:, None]
    assert np.allclose(u, dirs[grp], atol=1e-12)
    # full-step speeds of redirected atoms are the MB draws (same RNG stream)
    vmb = maxwell_boltzmann(np.random.default_rng(12), m, cfg.temperature)
    assert np.allclose(np.linalg.norm(vfull, axis=1), np.linalg.norm(vmb, axis=1))
    # full-step convention: no force correction
    v2, _, _ = pseudobd_reset(np.random.default_rng(12), x, F, m, box, cfg, dt, half_step_velocities=False)
    assert np.allclose(v2, vfull)
