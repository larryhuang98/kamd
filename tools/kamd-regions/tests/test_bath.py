from kamd_regions.analysis import Analysis
from kamd_regions.membrane import runs
from kamd_regions.system import load_system

from conftest import synthetic_membrane


def test_runs():
    assert runs([1, 2, 3, 7, 8, 10]) == [(1, 3), (7, 8), (10, 10)]
    assert runs([]) == []


def test_contiguous_bath(tmp_path):
    s = load_system([str(synthetic_membrane(tmp_path, "lipid21"))])
    an = Analysis(s)
    assert an.bath.contiguity(s) == []
    assert an.bath.solvent_atoms.size == 18
    assert an.bath.membrane_atoms.size == 24 * 12


def test_split_bath_is_reported_with_advice(tmp_path):
    s = load_system([str(synthetic_membrane(tmp_path, "lipid21", split_water=True))])
    an = Analysis(s)
    probs = an.bath.contiguity(s)
    assert len(probs) == 1 and probs[0]["bath"] == "solvent"
    assert probs[0]["n_ranges"] == 2
    assert "LIG" in probs[0]["message"] and "Reorder" in probs[0]["message"]
