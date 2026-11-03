"""Local FastAPI server for the GUI (``kamd-regions serve``)."""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

from fastapi import Body, FastAPI, HTTPException
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles

from .project import Settings
from .session import Session

STATIC = Path(__file__).parent / "static"


def create_app(session: Session | None, default_out: str | None = None, desktop: bool = False) -> FastAPI:
    """`session` may be None: the GUI then shows its start screen and opens files via /api/open."""
    app = FastAPI(title="kamd-regions", docs_url="/api/docs")
    state = {"session": session, "out": default_out or str(Path.cwd() / "kamd_regions_out")}

    def ses() -> Session:
        if state["session"] is None:
            raise HTTPException(409, "no system is loaded; open a structure first")
        return state["session"]

    @app.get("/api/app")
    def get_app():
        from . import __version__
        from .recent import recent_files

        return {
            "version": __version__,
            "desktop": desktop,
            "loaded": state["session"] is not None,
            "recent": recent_files(),
            "out_dir": state["out"],
        }

    @app.post("/api/open")
    def post_open(body: dict = Body(...)):
        """Load a new system (and optional KAMD input / project) from local paths."""
        from .recent import add_recent
        from .system import LoadError

        paths = [str(p) for p in body.get("paths") or [] if str(p).strip()]
        if not paths:
            raise HTTPException(400, "no files given")
        try:
            s = Session(paths, body.get("method", "auto"))
        except (LoadError, ValueError, FileNotFoundError, OSError) as exc:
            raise HTTPException(400, str(exc)) from exc
        state["session"] = s
        if body.get("out_dir"):
            state["out"] = str(body["out_dir"])
        else:  # write next to the input by default
            state["out"] = str(Path(paths[0]).resolve().parent / "kamd_regions_out")
        add_recent(paths)
        return {"ok": True, "notes": s.notes}

    @app.get("/api/session")
    def get_session():
        s = ses()
        return {
            "system": s.system_summary(),
            "project": s.project.to_dict(),
            "validation": s.validate(),
            "out_dir": state["out"],
        }

    @app.get("/api/structure")
    def get_structure():
        return ses().viewer_structure()

    @app.post("/api/suggest")
    def post_suggest(body: dict = Body(...)):
        s = ses()
        st = None
        if body.get("settings"):
            st = Settings(**{k: v for k, v in body["settings"].items() if k in Settings.__dataclass_fields__})
        try:
            p, notes = s.suggest(body.get("method", "auto"), body.get("params") or {}, st)
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc
        return {"project": p.to_dict(), "notes": notes, "validation": s.validate(p.to_dict())}

    @app.post("/api/validate")
    def post_validate(body: dict = Body(...)):
        try:
            return ses().validate(body["project"])
        except (KeyError, ValueError) as exc:
            raise HTTPException(400, str(exc)) from exc

    @app.post("/api/export")
    def post_export(body: dict = Body(...)):
        engines = tuple(body.get("engines") or ("amber", "gromacs", "openmm"))
        out_dir = body.get("out_dir") if body.get("write", False) else None
        res = ses().export(body["project"], engines, out_dir=out_dir)
        return res

    @app.post("/api/export.zip")
    def post_export_zip(body: dict = Body(...)):
        engines = tuple(body.get("engines") or ("amber", "gromacs", "openmm"))
        res = ses().export(body["project"], engines)
        if not res["ok"]:
            raise HTTPException(
                400, "validation errors: " + "; ".join(e["message"] for e in res["validation"]["errors"])
            )
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for name, text in res["files"].items():
                z.writestr(name, text)
        return Response(
            buf.getvalue(),
            media_type="application/zip",
            headers={"Content-Disposition": "attachment; filename=kamd_inputs.zip"},
        )

    @app.post("/api/project/save")
    def post_save(body: dict = Body(...)):
        path = body.get("path") or str(Path(state["out"]) / "project.json")
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        return {"path": ses().save(body["project"], path)}

    @app.post("/api/import")
    def post_import(body: dict = Body(...)):
        try:
            p, notes = ses().import_text(body["kind"], body["text"], body.get("ndx_text"))
        except Exception as exc:  # user input: report, don't crash
            raise HTTPException(400, f"import failed: {exc}") from exc
        return {"project": p.to_dict(), "notes": notes, "validation": ses().validate(p.to_dict())}

    @app.post("/api/project/load")
    def post_load(body: dict = Body(...)):
        from .project import Project

        try:
            p = Project.from_dict(body["project"])
        except Exception as exc:
            raise HTTPException(400, f"not a project file: {exc}") from exc
        try:
            p, notes = ses().adopt(p)  # migrates by sequence alignment when made on another system
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc
        return {"project": p.to_dict(), "notes": notes, "validation": ses().validate(p.to_dict())}

    if (STATIC / "index.html").exists():
        app.mount("/assets", StaticFiles(directory=STATIC / "assets"), name="assets")

        @app.get("/")
        def index():
            return FileResponse(STATIC / "index.html")
    else:

        @app.get("/")
        def no_frontend():
            return Response(
                "<h1>kamd-regions</h1><p>The frontend is not built. Run <code>npm run build</code> in "
                "tools/kamd-regions/frontend, or use the API at <a href='/api/docs'>/api/docs</a>.</p>",
                media_type="text/html",
            )

    return app
