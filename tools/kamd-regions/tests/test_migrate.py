"""Migrating regions between systems by sequence alignment."""

import numpy as np
import pytest

from kamd_regions.migrate import align, migrate
from kamd_regions.project import Project
from kamd_regions.system import load_system

from conftest import DATA

T4L = DATA / "t4l.pdb"


def rewrite(src, dst, chain_map=None, offset=0, drop=(), mutate=None, heavy_only=False, copies=1):
    """Write a modified copy of a PDB: rename chains, shift numbers, drop residues, rename residues."""
    lines = [ln for ln in open(src) if ln.startswith(("ATOM", "HETATM"))]
    out = []
    for c in range(copies):
        for ln in lines:
            l = ln  # noqa: E741
            n = int(l[22:26])
            if n in drop:
                continue
            name = l[12:16].strip()
            if heavy_only and (name.startswith("H") or l[76:78].strip() == "H"):
                continue
            chain = l[21] if copies == 1 else "AB"[c]
            chain = (chain_map or {}).get(chain, chain)
            resn = (mutate or {}).get(n, l[17:20])
            x = float(l[30:38]) + 40.0 * c
            out.append(l[:17] + f"{resn:>3s}" + " " + chain + f"{n + offset:4d}" + l[26:30] + f"{x:8.3f}" + l[38:])
        out.append("TER\n")
    dst.write_text("".join(out) + "END\n")
    return dst


def project_on(system, rule):
    lab = np.zeros(system.n_residues, dtype=int)
    for r in system.residues:
        if r.kind == "protein":
            lab[r.index] = rule(r)
    p = Project.from_labels(lab)
    from kamd_regions.migrate import residue_table

    p.residue_table = residue_table(system)
    p.residue_names = [r.name for r in system.residues]
    return p


def test_align_basic():
    pairs = align("ABCDEF", "ABXDF")
    got = {(i, j) for i, j in pairs if i is not None and j is not None}
    assert {(0, 0), (1, 1), (3, 3), (5, 4)} <= got
    assert (4, None) in pairs  # E deleted
    assert align("", "AB") == [(None, 0), (None, 1)]


def test_renumbered_truncated_renamed(tmp_path):
    # planning structure: chain X, numbers +10, first 3 residues missing, no hydrogens
    old = load_system([str(rewrite(T4L, tmp_path / "plan.pdb", {"A": "X"}, 10, drop=(1, 2, 3), heavy_only=True))])
    new = load_system([str(T4L)])
    p = project_on(old, lambda r: 1 if r.number - 10 < 60 else 2)
    q, rep = migrate(p, new)
    assert rep.chain_pairs[0]["identity"] == 1.0 and rep.chain_pairs[0]["gaps_new"] == 3
    for r in new.residues:
        if r.kind != "protein":
            continue
        want = 0 if r.number <= 3 else (1 if r.number < 60 else 2)
        assert q.assignment[r.index] == want, r.label()
    assert not rep.lost and rep.carried == old.n_residues - 0


def test_missing_loop_and_mutation_are_reported(tmp_path):
    old = load_system([str(T4L)])
    new = load_system([str(rewrite(T4L, tmp_path / "new.pdb", drop=(60, 61, 62), mutate={100: "TRP"}))])
    p = project_on(old, lambda r: 1 if r.number < 70 else 2)
    q, rep = migrate(p, new)
    assert sorted(x.split()[-1] for x in rep.lost) == ["60", "61", "62"]
    assert any(x.startswith("A:") and "100" in x and x.endswith("->TRP") for x in rep.mutated)
    num = {r.index: r.number for r in new.residues}
    for g in np.flatnonzero(q.assignment):
        assert q.assignment[g] == (1 if num[g] < 70 else 2)


def test_homodimer_pairs_chains_by_order_and_label(tmp_path):
    old = load_system([str(rewrite(T4L, tmp_path / "dimer.pdb", copies=2))])
    # region per chain and per half, so a chain swap would be detected
    p = project_on(old, lambda r: (1 if r.number < 80 else 2) + (2 if r.chain == "B" else 0))
    new = load_system([str(rewrite(T4L, tmp_path / "dimer2.pdb", {"A": "C", "B": "D"}, offset=5, copies=2))])
    q, rep = migrate(p, new)
    assert [(c["old"], c["new"]) for c in rep.chain_pairs] == [("A", "C"), ("B", "D")]
    for r in new.residues:
        if r.kind == "protein":
            want = (1 if r.number - 5 < 80 else 2) + (2 if r.chain == "D" else 0)
            assert q.assignment[r.index] == want


def test_session_adopts_project_from_another_system(tmp_path):
    from kamd_regions.session import Session

    plan = rewrite(T4L, tmp_path / "plan.pdb", {"A": "X"}, 10, drop=(1,), heavy_only=True)
    s_old = Session([str(plan)], "gnm")
    s_old.project.regions[0].name = "first"
    path = tmp_path / "plan.json"
    s_old.save(s_old.project.to_dict(), str(path))
    s_new = Session([str(T4L), str(path)])
    assert any("migrated" in n for n in s_new.notes)
    assert s_new.project.regions[0].name == "first"
    assert s_new.project.n_residues == s_new.system.n_residues
    old_sets = {
        r.name: {s_old.system.residues[i].number - 10 for i in s_old.project.residues_of(r.id)}
        for r in s_old.project.regions
    }
    new_sets = {
        r.name: {s_new.system.residues[i].number for i in s_new.project.residues_of(r.id)}
        for r in s_new.project.regions
    }
    assert old_sets == new_sets
    # same residues in the same order: used as it is, no migration needed
    s_same = Session([str(rewrite(T4L, tmp_path / "same.pdb", {"A": "Y"}, 20, drop=(1,))), str(path)])
    assert not any("migrated" in n for n in s_same.notes)


def test_migrate_onto_built_md_system(tmp_path):
    """Plan on a crystal-like PDB, build a solvated amber14 system with OpenMM, write a prmtop, migrate."""
    openmm = pytest.importorskip("openmm")
    import parmed
    from openmm import app, unit

    from kamd_regions.analysis import Analysis
    from kamd_regions.export import export_files, verify
    from kamd_regions.session import Session

    pdb = app.PDBFile(str(T4L))
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    m = app.Modeller(pdb.topology, pdb.positions)
    m.addSolvent(ff, padding=0.6 * unit.nanometer, ionicStrength=0.15 * unit.molar)
    system = ff.createSystem(m.topology, nonbondedMethod=app.PME, rigidWater=False)
    st = parmed.openmm.load_topology(m.topology, system, xyz=m.positions)
    st.save(str(tmp_path / "sys.parm7"))
    st.save(str(tmp_path / "sys.rst7"))
    assert openmm.__version__

    plan = rewrite(T4L, tmp_path / "plan.pdb", {"A": "X"}, 10, drop=(1, 2, 3), heavy_only=True)
    s_old = Session([str(plan)], "gnm")
    s_old.save(s_old.project.to_dict(), str(tmp_path / "plan.json"))
    s_new = Session([str(tmp_path / "sys.parm7"), str(tmp_path / "sys.rst7"), str(tmp_path / "plan.json")])
    assert s_new.system.has_topology and s_new.an.bath.solvent_atoms.size > 1000
    q = s_new.project
    for r in s_old.project.regions:
        old_nums = {s_old.system.residues[i].number - 10 for i in s_old.project.residues_of(r.id)}
        new_nums = {s_new.system.residues[i].number for i in q.residues_of(r.id)}
        assert old_nums == new_nums
    files = export_files(s_new.an, q, ("amber", "gromacs", "openmm"))
    rep = verify(s_new.an, q, files)
    assert rep["ok"], rep
    assert isinstance(s_new.an, Analysis)


def test_gui_load_project_migrates(tmp_path):
    from fastapi.testclient import TestClient

    from kamd_regions.server import create_app
    from kamd_regions.session import Session

    s_old = Session([str(rewrite(T4L, tmp_path / "plan.pdb", {"A": "X"}, 10, drop=(1, 2), heavy_only=True))], "gnm")
    c = TestClient(create_app(Session([str(T4L)], "gnm"), str(tmp_path)))
    r = c.post("/api/project/load", json={"project": s_old.project.to_dict()}).json()
    assert any("migrated" in n for n in r["notes"])
    assert r["project"]["n_residues"] == 162
    assert len(r["project"]["regions"]) == len(s_old.project.regions)
