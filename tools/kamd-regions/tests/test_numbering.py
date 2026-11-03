"""Topology residue index vs. PDB chain/resid vs. Amber masks vs. atom numbers."""

import numpy as np

from kamd_regions.analysis import Analysis
from kamd_regions.export import amber_namelist, gromacs_files, read_ndx
from kamd_regions.masks import residue_summary
from kamd_regions.project import Project
from kamd_regions.system import load_system

from conftest import AMBER_TEST, pdb_line


def two_chain_pdb(tmp_path):
    lines = []
    serial = 1
    # chain B numbered 101.., chain A numbered 5.. but listed second
    for chain, start in (("B", 101), ("A", 5)):
        for k in range(4):
            for j, (an, el) in enumerate((("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"))):
                lines.append(
                    pdb_line(serial, an, "GLY", chain, start + k, 3.8 * k + 0.4 * j, 10.0 * (chain == "A"), 0.0, el)
                )
                serial += 1
        lines.append("TER")
    for k in range(3):
        lines.append(pdb_line(serial, "O", "HOH", "W", 1 + k, 30 + 3 * k, 0, 0, "O", het=True))
        serial += 1
    p = tmp_path / "two.pdb"
    p.write_text("\n".join(lines) + "\nEND\n")
    return p


def test_internal_index_is_file_order(tmp_path):
    s = load_system([str(two_chain_pdb(tmp_path))])
    assert [r.chain for r in s.residues[:8]] == ["B"] * 4 + ["A"] * 4
    assert [r.number for r in s.residues[:8]] == [101, 102, 103, 104, 5, 6, 7, 8]
    assert not s.has_topology
    an = Analysis(s)
    lab = np.zeros(s.n_residues, dtype=int)
    lab[[0, 1]] = 1  # B:101-102
    lab[[4, 5, 6]] = 2  # A:5-7
    p = Project.from_labels(lab)
    text = amber_namelist(an, p)
    # Amber masks are sequential 1-based residue numbers, not PDB resids, and have no chains
    assert 'refer_mask(1) = ":1-2"' in text
    assert 'refer_mask(2) = ":5-7"' in text
    assert "WARNING: written without a topology" in text
    assert residue_summary(s, [4, 5, 6]) == "A:5-7"
    assert residue_summary(s, [0, 1, 4]) == "B:101-102, A:5"
    # GROMACS: 1-based atom numbers
    nd = read_ndx(gromacs_files(an, p)["index.ndx"])
    assert nd["region_1"].tolist() == list(range(0, 8))
    assert nd["region_2"].tolist() == list(range(16, 28))
    assert nd["KAMD_SOLVENT"].tolist() == [32, 33, 34]


def test_prmtop_numbers_are_sequential():
    s = load_system([str(AMBER_TEST / "prmtop"), str(AMBER_TEST / "restrt")])
    assert s.has_topology and s.topology_format == "amber"
    assert [r.number for r in s.residues[:3]] == [1, 2, 3]
    assert s.residues[0].chain == "A"  # labelled from the molecule table
    # residue 2 (LYS) starts at atom 20 in the prmtop
    assert s.residues[1].start == 19 and s.residues[1].name == "LYS"
    assert s.n_atoms == 29467
