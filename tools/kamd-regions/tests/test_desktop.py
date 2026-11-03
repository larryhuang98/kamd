"""Start screen / open-files flow used by the desktop app and `serve` without files."""

from fastapi.testclient import TestClient

from kamd_regions.server import create_app

from conftest import AMBER_TEST, DATA


def test_open_flow_and_recent(tmp_path, monkeypatch):
    monkeypatch.setenv("KAMD_REGIONS_CONFIG", str(tmp_path / "cfg"))
    c = TestClient(create_app(None, None, desktop=True))
    info = c.get("/api/app").json()
    assert info == {**info, "loaded": False, "desktop": True, "recent": []}
    assert c.get("/api/session").status_code == 409  # nothing loaded yet
    bad = c.post("/api/open", json={"paths": [str(tmp_path / "missing.pdb")]})
    assert bad.status_code == 400 and "not found" in bad.json()["detail"]
    files = [str(AMBER_TEST / "prmtop"), str(AMBER_TEST / "restrt"), str(AMBER_TEST / "md_test.in")]
    assert c.post("/api/open", json={"paths": files}).json()["ok"]
    s = c.get("/api/session").json()
    assert len(s["project"]["regions"]) == 4 and s["system"]["n_atoms"] == 29467
    info = c.get("/api/app").json()
    assert info["loaded"] and info["recent"][0] == files
    assert info["out_dir"].endswith("kamd_regions_out")
    # opening another system replaces the session and moves it to the top of the recent list
    assert c.post("/api/open", json={"paths": [str(DATA / "t4l.pdb")]}).json()["ok"]
    assert c.get("/api/session").json()["system"]["n_residues"] == 162
    assert [len(r) for r in c.get("/api/app").json()["recent"]] == [1, 3]


def test_bridge_save_text(tmp_path):
    from kamd_regions.desktop import Bridge

    class FakeWindow:
        def create_file_dialog(self, kind, **kw):
            return str(tmp_path / kw.get("save_filename", "x"))

    b = Bridge()
    b.window = FakeWindow()
    path = b.save_text("kamd_amber.in", "&kamd\n/\n")
    assert path.endswith("kamd_amber.in") and (tmp_path / "kamd_amber.in").read_text() == "&kamd\n/\n"
