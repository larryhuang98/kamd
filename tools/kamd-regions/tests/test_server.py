from fastapi.testclient import TestClient

from kamd_regions.server import create_app


def test_api_flow(amber_session, tmp_path):
    c = TestClient(create_app(amber_session, str(tmp_path)))
    s = c.get("/api/session").json()
    assert s["system"]["n_atoms"] == 29467
    assert len(s["system"]["residues"]["index"]) == 39
    assert len(s["project"]["regions"]) == 4
    st = c.get("/api/structure").json()
    assert len(st["solute_res"]) == st["solute_pdb"].count("\nATOM") + 1
    r = c.post("/api/suggest", json={"method": "gnm", "params": {}}).json()
    assert len(r["project"]["regions"]) >= 2
    v = c.post("/api/validate", json={"project": r["project"]}).json()
    assert v["ok"]["amber"]
    e = c.post("/api/export", json={"project": s["project"], "engines": ["amber", "gromacs", "openmm"]}).json()
    assert e["ok"] and e["verification"]["ok"]
    assert set(e["files"]) >= {"kamd_amber.in", "index.ndx", "kamd.dat", "kamd_openmm.json", "project.json"}
    z = c.post("/api/export.zip", json={"project": s["project"]})
    assert z.status_code == 200 and z.content[:2] == b"PK"
    saved = c.post("/api/project/save", json={"project": s["project"]}).json()
    assert (tmp_path / "project.json").exists() and saved["path"].endswith("project.json")
    imp = c.post("/api/import", json={"kind": "amber", "text": "&kamd\n refer_mask(1)=':2-5', refer_mask(2)=':10-12',\n/"})
    assert imp.status_code == 200 and len(imp.json()["project"]["regions"]) == 2
    bad = c.post("/api/import", json={"kind": "amber", "text": "nothing"})
    assert bad.status_code == 400


def test_cli_suggest_and_export(tmp_path):
    from conftest import AMBER_TEST, DATA
    from kamd_regions.cli import main

    out = tmp_path / "p.json"
    rc = main(["suggest", str(AMBER_TEST / "prmtop"), str(AMBER_TEST / "restrt"), "--method", "gnm", "--out",
               str(out), "--export-dir", str(tmp_path / "ex")])
    assert rc == 0 and out.exists() and (tmp_path / "ex" / "kamd_amber.in").exists()
    rc = main(["export", str(AMBER_TEST / "prmtop"), str(AMBER_TEST / "restrt"), str(out), "--out-dir",
               str(tmp_path / "ex2")])
    assert rc == 0
    rc = main(["suggest", str(DATA / "3emn.pdb"), "--out", str(tmp_path / "v.json")])
    assert rc == 0
