import numpy as np
import pytest

from kamd_regions.masks import atom_mask, compress, expand, parse_simple_mask, residue_mask
from kamd_regions.project import Project, rle_decode, rle_encode


def test_compress_examples():
    assert compress([1, 2, 3, 77, 120, 121, 131]) == "1-3,77,120-121,131"
    assert compress([5]) == "5"
    assert compress([]) == ""
    assert compress([3, 1, 2, 2]) == "1-3"


@pytest.mark.parametrize("seed", range(20))
def test_compress_expand_roundtrip(seed):
    rng = np.random.default_rng(seed)
    v = np.unique(rng.integers(1, 3000, size=rng.integers(1, 400)))
    assert np.array_equal(expand(compress(v)), v)


def test_residue_and_atom_mask_are_one_based():
    assert residue_mask([0, 1, 2, 3, 76, 119, 120]) == ":1-4,77,120-121"
    assert atom_mask([19, 20, 21]) == "@20-22"
    kind, idx = parse_simple_mask(":1-4,77")
    assert kind == "residue" and idx.tolist() == [0, 1, 2, 3, 76]
    kind, idx = parse_simple_mask("@623-625")
    assert kind == "atom" and idx.tolist() == [622, 623, 624]
    with pytest.raises(ValueError):
        parse_simple_mask(":WAT")


@pytest.mark.parametrize("seed", range(10))
def test_rle_roundtrip(seed):
    rng = np.random.default_rng(seed)
    a = np.repeat(rng.integers(0, 5, 30), rng.integers(1, 9, 30))
    assert np.array_equal(rle_decode(rle_encode(a), a.size), a)


def test_project_json_roundtrip(tmp_path):
    lab = np.array([0, 1, 1, 2, 2, 2, 0, 3, 1])
    p = Project.from_labels(lab, {1: 1, 2: 2, 3: 1}, names=["ALA"] * 9)
    p.regions[0].refer_link = p.regions[2].id
    p.settings.membrane = {"axis": 3, "lo": -18.0, "hi": 18.0}
    f = tmp_path / "project.json"
    p.save(f)
    q = Project.load(f)
    assert np.array_equal(q.assignment, lab)
    assert [r.region_type for r in q.regions] == [1, 2, 1]
    assert q.regions[0].refer_link == q.regions[2].id
    assert q.settings.membrane["hi"] == 18.0
