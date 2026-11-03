import numpy as np
import pytest

from kamd_regions.methods import suggest
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
        assert ompa.s.res_nheavy[res["labels"] == g].sum() >= 300
    # detergent (C8E) residues are assignable but not used by the auto methods
    c8e = [r.index for r in ompa.s.residues if r.name == "C8E"]
    assert np.all(res["labels"][c8e] == 0)


def test_gnm_matches_prody(ompa):
    prody = pytest.importorskip("prody")
    vals, vecs, nz = ompa.gnm(7.3)
    g = prody.GNM()
    g.buildKirchhoff(ompa.rep_xyz[ompa.polymer], cutoff=7.3)
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
