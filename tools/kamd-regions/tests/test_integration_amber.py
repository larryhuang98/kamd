"""Integration test on amber/test (prmtop + restrt + md_test.in, four manual groups)."""

import json

import numpy as np

from conftest import AMBER_TEST, REPO
from kamd_openmm.amber_input import amber_mask_indices, parse_namelist
from kamd_regions.compare import compare_report
from kamd_regions.export import export_files, read_ndx, verify

MASKS = ["@20-96", "@285-339", "@366-482", "@568-611"]


def test_import_matches_original_masks(amber_session):
    s = amber_session.system
    p = amber_session.project
    assert len(p.regions) == 4
    for r, m in zip(p.regions, MASKS, strict=True):
        want = amber_mask_indices(s.structure, m)
        got = s.atoms_of(p.residues_of(r.id))
        assert np.array_equal(np.sort(got), np.sort(want)), m
    assert not any("covers" in n for n in amber_session.notes)  # masks are whole residues
    assert p.settings.temperature == 355.0 and p.settings.seed == 1060


def test_export_roundtrip_atom_by_atom(amber_session, tmp_path):
    an, p = amber_session.an, amber_session.project
    files = export_files(an, p)
    rep = verify(an, p, files)
    assert rep["ok"], rep
    assert rep["kamd_config"]["ok"] is True
    nl = parse_namelist(files["kamd_amber.in"], "kamd")
    for i, m in enumerate(MASKS):
        a = amber_mask_indices(an.s.structure, nl["refer_mask"][i + 1])
        b = amber_mask_indices(an.s.structure, m)
        assert np.array_equal(np.sort(a), np.sort(b))
    sol = amber_mask_indices(an.s.structure, nl["solvent_mask"])
    assert np.array_equal(np.sort(sol), np.arange(622, 29467))  # == @623-29467
    # identical atoms to the hand-written GROMACS index in gromacs/validation
    ref = read_ndx((REPO / "gromacs" / "validation" / "index.ndx").read_text())
    ours = read_ndx(files["index.ndx"])
    for k, name in enumerate(["dom1", "dom2", "dom3", "dom4"]):
        assert np.array_equal(ours[f"group_{k + 1}"], ref[name])
    assert np.array_equal(ours["KAMD_SOLVENT"], ref["bath"])
    d = json.loads(files["kamd_openmm.json"])
    assert d["refer_groups"][0] == list(range(19, 96))
    # OpenMM: the same config as kamd_openmm.config_from_amber on the original input
    from kamd_openmm import config_from_amber

    cfg = config_from_amber(str(AMBER_TEST / "md_test.in"), an.s.structure)
    for g, h in zip(cfg.refer_groups, d["refer_groups"], strict=True):
        assert list(g) == h


def test_openmm_example_runs(amber_session, tmp_path, monkeypatch):
    files = export_files(amber_session.an, amber_session.project, ("openmm",))
    (tmp_path / "kamd_openmm.json").write_text(files["kamd_openmm.json"])
    code = files["kamd_openmm_example.py"].split("# `simulation`")[0]
    code = code.replace("KAMDConfig, KAMDReporter", "KAMDConfig")  # the reporter needs OpenMM
    monkeypatch.chdir(tmp_path)
    ns: dict = {}
    exec(code, ns)
    ns["cfg"].validate(29467)
    assert ns["cfg"].ngroups == 4


def test_compare_report(amber_session):
    text = compare_report(amber_session, ["gnm", "contacts"])
    assert "## gnm" in text and "## contacts" in text and "Jaccard" in text
