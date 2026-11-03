from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).parent
DATA = HERE / "data"
REPO = HERE.parents[2]
AMBER_TEST = REPO / "amber" / "test"


def pdb_line(serial, name, resn, chain, resi, x, y, z, el, het=False):
    rec = "HETATM" if het else "ATOM  "
    nm = f" {name:<3s}" if len(name) < 4 else name
    return f"{rec}{serial:5d} {nm} {resn:<4s}{chain}{resi:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}"


def synthetic_membrane(tmp_path: Path, style: str = "lipid21", axis: int = 3, split_water: bool = False) -> Path:
    """Small helix bundle-free test system: a 3-residue peptide, a lipid bilayer, water and ions.

    style "lipid21": every lipid is 3 residues (PA PC OL) without bonds in the file;
    style "charmm": one POPC residue per lipid.
    The membrane normal is `axis`; phosphates at +-19 A.
    """
    lines = ["CRYST1   60.000   60.000   80.000  90.00  90.00  90.00 P 1           1"]
    serial = 1
    resi = 1

    def place(u, v, w):
        # (in-plane u, in-plane v, normal w) -> xyz
        p = [0.0, 0.0, 0.0]
        others = [k for k in range(3) if k != axis - 1]
        p[others[0]], p[others[1]], p[axis - 1] = u, v, w
        return p

    # peptide (chain A), roughly along the normal
    for k, rn in enumerate(["ALA", "GLY", "LEU"]):
        for j, (an, el) in enumerate((("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"))):
            x, y, z = place(0.5 * j, 0.3 * j, -4.0 + 3.8 * k + 0.9 * j)
            lines.append(pdb_line(serial, an, rn, "A", resi, x, y, z, el))
            serial += 1
        resi += 1
    lines.append("TER")
    n_lip = 0
    for s in (-1, 1):
        for i in range(4):
            for j in range(4):
                u, v = -20 + 13 * i, -20 + 13 * j
                if abs(u) < 8 and abs(v) < 8:
                    continue
                n_lip += 1
                # head at |w| = 19, tails pointing to the centre, consecutive atoms 1.5 A apart
                if style == "lipid21":
                    # tail - head - tail, bonded in a chain (1.5 A steps); the second tail runs back
                    # down next to the head so the molecule stays compact
                    frags = [
                        ("PA", [(0.0, w) for w in (10.0, 11.5, 13.0, 14.5)]),
                        ("PC", [(0.0, w) for w in (16.0, 17.5, 19.0, 20.5)]),
                        ("OL", [(1.5, w) for w in (20.5, 19.0, 17.5, 16.0)]),
                    ]
                    names = ["C12", "C13", "C14", "C15"]
                    for rn, coords in frags:
                        for k, (du, w) in enumerate(coords):
                            an = ["C11", "O11", "P", "O12"][k] if rn == "PC" else names[k]
                            x, y, z = place(u + du, v, s * w)
                            lines.append(pdb_line(serial, an, rn, " ", resi, x, y, z, an[0], het=True))
                            serial += 1
                        resi += 1
                else:
                    for k, (an, el) in enumerate((("C1", "C"), ("C2", "C"), ("P", "P"), ("O11", "O"), ("C21", "C"))):
                        x, y, z = place(u, v, s * (14.0 + 1.5 * k))
                        lines.append(pdb_line(serial, an, "POPC", " ", resi, x, y, z, el, het=True))
                        serial += 1
                    resi += 1
                lines.append("TER")
    waters = 6

    def water(k):
        nonlocal serial, resi
        x, y, z = place(-25 + 5 * k, 25, 32.0)
        for an, el, d in (("O", "O", 0.0), ("H1", "H", 0.96), ("H2", "H", -0.96)):
            lines.append(pdb_line(serial, an, "WAT", " ", resi, x + d, y, z, el, het=True))
            serial += 1
        resi += 1

    for k in range(waters):
        if split_water and k == 3:
            x, y, z = place(25, -25, -32.0)
            lines.append(pdb_line(serial, "NA", "NA", " ", resi, x, y, z, "NA", het=True))
            serial += 1
            resi += 1
            # a ligand between the water blocks breaks the contiguity
            x, y, z = place(25, -20, -32.0)
            lines.append(pdb_line(serial, "C1", "LIG", " ", resi, x, y, z, "C", het=True))
            serial += 1
            resi += 1
        water(k)
    lines.append("END")
    p = tmp_path / f"memb_{style}_{axis}{'_split' if split_water else ''}.pdb"
    p.write_text("\n".join(lines) + "\n")
    synthetic_membrane.n_lipids = n_lip
    return p


@pytest.fixture(scope="session")
def amber_session():
    from kamd_regions.session import Session

    return Session([str(AMBER_TEST / "prmtop"), str(AMBER_TEST / "restrt"), str(AMBER_TEST / "md_test.in")])


@pytest.fixture(scope="session")
def vdac():
    from kamd_regions.analysis import Analysis
    from kamd_regions.system import load_system

    return Analysis(load_system([str(DATA / "3emn.pdb")]))


@pytest.fixture(scope="session")
def ompa():
    from kamd_regions.analysis import Analysis
    from kamd_regions.system import load_system

    return Analysis(load_system([str(DATA / "1qjp.pdb")]))


def rng():
    return np.random.default_rng(0)
