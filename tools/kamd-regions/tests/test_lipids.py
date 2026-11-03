import pytest

from conftest import synthetic_membrane
from kamd_regions.analysis import Analysis
from kamd_regions.residues import classify
from kamd_regions.system import load_system


def test_lipid21_fragments_are_merged_into_molecules(tmp_path):
    s = load_system([str(synthetic_membrane(tmp_path, "lipid21"))])
    lip = [r for r in s.residues if r.kind == "lipid"]
    assert len(lip) == 3 * 24
    mols = {}
    for r in lip:
        mols.setdefault(r.molecule, []).append(r.name)
    assert len(mols) == 24
    assert all(v == ["PA", "PC", "OL"] for v in mols.values())


def test_lipid21_merge_from_prmtop_bonds_not_needed_for_charmm(tmp_path):
    s = load_system([str(synthetic_membrane(tmp_path, "charmm"))])
    lip = [r for r in s.residues if r.kind == "lipid"]
    assert len(lip) == 24 and len({r.molecule for r in lip}) == 24


@pytest.mark.parametrize("style,axis", [("lipid21", 3), ("charmm", 3), ("charmm", 1), ("lipid21", 2)])
def test_membrane_normal_and_slab(tmp_path, style, axis):
    s = load_system([str(synthetic_membrane(tmp_path, style, axis))])
    m = Analysis(s).membrane
    assert m is not None and m.source == "lipids"
    assert m.axis == axis
    expect = 19.0 if style == "lipid21" else 17.0
    assert m.lo == pytest.approx(-expect, abs=0.5) and m.hi == pytest.approx(expect, abs=0.5)
    assert m.n_lipids == 24


def test_classify_names():
    assert classify("POPC", ["P"] * 30, ["P"] * 30) == "lipid"
    assert classify("PA", ["C"] * 10, ["C"] * 10) == "lipid"
    assert classify("CL", ["CL"], ["Cl"]) == "ion"
    assert classify("CA", ["CA"], ["Ca"]) == "ion"
    assert classify("WAT", ["O", "H1", "H2"], ["O", "H", "H"]) == "water"
    assert classify("SOD", ["SOD"], ["Na"]) == "ion"
    assert classify("HSD", ["N", "CA", "C", "O", "CB"], ["N", "C", "C", "O", "C"]) == "protein"
    assert classify("DUM", ["O"], ["O"]) == "dummy"
