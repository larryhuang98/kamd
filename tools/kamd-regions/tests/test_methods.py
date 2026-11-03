import numpy as np
import pytest

from kamd_regions.methods import postprocess, suggest
from kamd_regions.methods.tm import tm_segments


def check_shape(an, res):
    lab = res["labels"]
    assert lab.shape == (an.s.n_residues,)
    ids = sorted(set(lab[lab > 0].tolist()))
    assert ids == list(range(1, len(ids) + 1))  # renumbered 1..K
    assert len(ids) >= 2
    assert np.all(lab[~an.assignable] == 0)  # bath / dummies never assigned
    assert set(res["region_type"]) == set(ids)
    for g in ids:  # no region below the (possibly lowered) size threshold unless needed
        assert an.s.res_nheavy[lab == g].sum() > 0
    return ids


@pytest.mark.parametrize("method", ["gnm", "contacts"])
def test_network_methods_on_ompa(ompa, method):
    res = suggest(ompa, method)
    ids = check_shape(ompa, res)
    assert len(ids) <= 32
    for g in ids:
        assert ompa.s.res_nheavy[res["labels"] == g].sum() >= postprocess.DEFAULTS["min_heavy"]
    # detergent (C8E) residues are assignable but not used by the auto methods
    c8e = [r.index for r in ompa.s.residues if r.name == "C8E"]
    assert np.all(res["labels"][c8e] == 0)


def test_gnm_matches_prody(ompa):
    prody = pytest.importorskip("prody")
    vals, vecs, nz = ompa.gnm(7.3)
    g = prody.GNM()
    g.buildKirchhoff(ompa.rep_xyz[ompa.gnm_nodes], cutoff=7.3)
    g.calcModes(10, zeros=False)
    assert np.allclose(vals[nz : nz + 10], g.getEigvals(), rtol=1e-6)


def test_gnm_k_parameter(ompa):
    res = suggest(ompa, "gnm", {"n_regions": 3, "min_heavy": 50})
    assert res["info"]["k"] == 3
    assert len(set(res["labels"][res["labels"] > 0].tolist())) == 3


def test_contacts_resolution_changes_granularity(ompa):
    lo = suggest(ompa, "contacts", {"resolution": 0.1, "min_heavy": 20})
    hi = suggest(ompa, "contacts", {"resolution": 1.5, "min_heavy": 20})
    assert hi["labels"].max() >= lo["labels"].max()


def test_tm_on_vdac_finds_19_strands(vdac):
    assert vdac.membrane is not None and vdac.membrane.source == "opm"
    assert vdac.membrane.axis == 3
    assert vdac.membrane.lo == pytest.approx(-11.5) and vdac.membrane.hi == pytest.approx(11.5)
    segs = tm_segments(vdac)
    assert len(segs) == 19 and all(s["kind"] == "E" for s in segs)
    res = suggest(vdac, "tm")
    ids = check_shape(vdac, res)
    assert all(res["region_type"][g] == 2 for g in ids)
    # every TM segment ends up (mostly) in a group, the N-terminal helix inside the barrel does not
    for s in segs:
        assert (res["labels"][s["start"] : s["stop"]] > 0).mean() > 0.8
    assert np.all(res["labels"][:20] == 0)


def test_tm_per_group(ompa):
    res = suggest(ompa, "tm", {"per_group": 1, "min_heavy": 1})
    assert len(tm_segments(ompa)) == 8
    assert res["labels"].max() == 8


def test_auto_method(vdac, ompa, amber_session):
    from kamd_regions.methods import auto_method

    assert auto_method(vdac) == "tm"
    assert auto_method(amber_session.an) == "gnm"


def test_methods_on_small_amber_peptide(amber_session):
    an = amber_session.an
    for m in ("gnm", "contacts"):
        res = suggest(an, m)
        ids = check_shape(an, res)
        assert 2 <= len(ids) <= 4
    res = suggest(an, "tm")
    assert res["labels"].max() == 0 and any("no membrane" in n for n in res["notes"])


# --------------------------------------------------------------------------------------------
# multi-domain proteins with literature domains (benchmark: docs/gnm_multidomain.md)
ADK = {"CORE": [(1, 29), (60, 121), (160, 214)], "NMP": [(30, 59)], "LID": [(122, 159)]}
T4L = {"C": [(1, 12), (60, 164)], "N": [(13, 59)]}


def _truth(an, ref):
    num = np.array([r.number for r in an.s.residues])
    t = np.zeros(an.s.n_residues, dtype=int)
    for k, rr in enumerate(ref.values()):
        for a, b in rr:
            t[(num >= a) & (num <= b) & (an.s.res_kind == "protein")] = k + 1
    return t


@pytest.mark.parametrize(
    "name,ref,params,min_ari",
    [
        ("adk_open.pdb", ADK, {"n_regions": 3}, 0.6),
        ("t4l.pdb", T4L, {}, 0.55),
    ],
)
def test_gnm_recovers_domains(name, ref, params, min_ari):
    from kamd_regions.analysis import Analysis
    from kamd_regions.compare import adjusted_rand
    from kamd_regions.system import load_system

    from conftest import DATA

    an = Analysis(load_system([str(DATA / name)]))
    res = suggest(an, "gnm", params)
    lab, t = res["labels"], _truth(an, ref)
    assert lab.max() == len(ref)
    both = (lab > 0) & (t > 0)
    assert adjusted_rand(t[both], lab[both]) >= min_ari
    # every literature domain is the majority of exactly one region
    owners = {int(np.bincount(t[(lab == g) & (t > 0)]).argmax()) for g in range(1, lab.max() + 1)}
    assert owners == set(range(1, len(ref) + 1))


def test_tm_on_gpcr_b2ar_t4l():
    """beta2-adrenergic receptor with T4 lysozyme in ICL3 (OPM 2rh1, one copy)."""
    from kamd_regions.analysis import Analysis
    from kamd_regions.masks import residue_summary
    from kamd_regions.system import load_system

    from conftest import DATA

    an = Analysis(load_system([str(DATA / "2rh1_monomer.pdb")]))
    assert an.membrane is not None and an.membrane.thickness == pytest.approx(31.8)
    res = suggest(an, "tm")
    lab, rt = res["labels"], res["region_type"]
    tm = [g for g in rt if rt[g] == 2]
    sol = [g for g in rt if rt[g] == 1]
    assert len(tm) == 7  # one group per helix
    assert len(sol) == 1
    t4l = [r.index for r in an.s.residues if r.number >= 1002 and r.kind == "protein"]
    assert (lab[t4l] == sol[0]).mean() == 1.0  # the T4L insert is the solution group
    # summaries break at numbering gaps (230 -> 1002 -> 263)
    assert "C:1002-1161" in residue_summary(an.s, np.flatnonzero(lab == sol[0]))
