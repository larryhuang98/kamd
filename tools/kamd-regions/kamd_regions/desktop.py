"""Desktop application: the GUI in a native window (``kamd-regions app``).

The backend server runs in a background thread on 127.0.0.1 (random free port); the window is a
native web view (pywebview: WebKit on macOS, WebView2 on Windows, Qt WebEngine on Linux), so there
is no browser and no address bar. A small bridge gives the page native file dialogs and the menu.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

STRUCTURE_TYPES = (
    "Structures and topologies (*.pdb;*.ent;*.cif;*.mmcif;*.gro;*.prmtop;*.parm7;*.top;*.tpr;*.inpcrd;*.rst7;"
    "*.ncrst;*.crd;*.in;*.mdin;*.dat;*.ndx;*.json)",
    "All files (*.*)",
)


def _say(msg: str) -> None:
    if sys.stdout is not None:  # None in a windowed (no console) Windows build
        print(msg, flush=True)


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _dialog_kind(webview, name: str):
    fd = getattr(webview, "FileDialog", None)  # pywebview >= 5
    if fd is not None:
        return getattr(fd, name)
    return getattr(webview, {"OPEN": "OPEN_DIALOG", "SAVE": "SAVE_DIALOG", "FOLDER": "FOLDER_DIALOG"}[name])


class Bridge:
    """Methods callable from the page as ``window.pywebview.api.<name>(...)``."""

    def __init__(self):
        self.window = None

    def _webview(self):
        import webview

        return webview

    def open_dialog(self, multiple: bool = True):
        wv = self._webview()
        r = self.window.create_file_dialog(
            _dialog_kind(wv, "OPEN"), allow_multiple=multiple, file_types=STRUCTURE_TYPES
        )
        return list(r) if r else []

    def folder_dialog(self):
        wv = self._webview()
        r = self.window.create_file_dialog(_dialog_kind(wv, "FOLDER"))
        return (r[0] if isinstance(r, (list, tuple)) else r) if r else None

    def save_dialog(self, name: str):
        wv = self._webview()
        r = self.window.create_file_dialog(_dialog_kind(wv, "SAVE"), save_filename=name)
        return (r[0] if isinstance(r, (list, tuple)) else r) if r else None

    def save_text(self, name: str, text: str):
        """Ask where to save `name`, write `text` there, return the path (None if cancelled)."""
        path = self.save_dialog(name)
        if not path:
            return None
        Path(path).write_text(text)
        return str(path)

    def reveal(self, path: str):
        """Show a file or folder in the system file manager."""
        p = Path(path)
        target = p if p.is_dir() else p.parent
        if sys.platform == "darwin":
            subprocess.Popen(["open", str(target)])
        elif os.name == "nt":
            os.startfile(str(target))  # noqa: S606 - local folder chosen by the user
        else:
            subprocess.Popen(["xdg-open", str(target)])
        return True

    def quit(self):
        self.window.destroy()


def _menu(webview, call):
    try:
        from webview.menu import Menu, MenuAction, MenuSeparator
    except ImportError:  # pragma: no cover - old pywebview
        return None
    return [
        Menu(
            "File",
            [
                MenuAction("Open…", lambda: call("open")),
                MenuSeparator(),
                MenuAction("Save project…", lambda: call("save")),
                MenuAction("Load project…", lambda: call("load")),
                MenuAction("Import KAMD input…", lambda: call("import")),
                MenuAction("Export to folder…", lambda: call("export")),
            ],
        ),
        Menu(
            "Edit",
            [
                MenuAction("Undo", lambda: call("undo")),
                MenuAction("Redo", lambda: call("redo")),
            ],
        ),
        Menu(
            "View",
            [
                MenuAction("Toggle dark theme", lambda: call("theme")),
                MenuAction("Zoom to selection", lambda: call("zoom")),
            ],
        ),
    ]


def run(
    files: list[str] | None = None,
    method: str = "auto",
    out_dir: str | None = None,
    selftest: str | None = None,
    debug: bool = False,
) -> int:
    try:
        import webview
    except ImportError:
        print("the desktop app needs pywebview: pip install 'kamd-regions[desktop]'", file=sys.stderr)
        return 2
    import uvicorn

    from .recent import add_recent
    from .server import create_app
    from .session import Session

    session = None
    if files:
        session = Session(files, method)
        add_recent(files)
        if out_dir is None:
            out_dir = str(Path(files[0]).resolve().parent / "kamd_regions_out")
    app = create_app(session, out_dir, desktop=True)
    port = _free_port()
    server = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=port, log_level="warning"))
    th = threading.Thread(target=server.run, daemon=True)
    th.start()
    for _ in range(200):
        if server.started:
            break
        time.sleep(0.05)

    bridge = Bridge()
    holder: dict = {}

    def call(cmd: str):
        holder["w"].evaluate_js(f"window.kamdMenu && window.kamdMenu({json.dumps(cmd)})")

    menu = _menu(webview, call)
    extra = {"menu": menu} if menu else {}
    window = webview.create_window(
        "KAMD Regions",
        f"http://127.0.0.1:{port}/",
        js_api=bridge,
        width=1440,
        height=920,
        min_size=(900, 600),
        text_select=True,
        **extra,
    )
    holder["w"] = window
    bridge.window = window

    result = {"ok": False}

    def check():
        # wait for the page to report itself ready, print what it shows, close.
        # evaluate_js before the web view has loaded can block forever (WebView2), so wait first
        loaded = getattr(getattr(window, "events", None), "loaded", None)
        if loaded is not None:
            loaded.wait(90)
        deadline = time.time() + 90
        state = None
        while time.time() < deadline:
            try:
                state = window.evaluate_js("JSON.stringify(window.__kamdState || null)")
            except Exception:
                state = None
            if state and state != "null" and json.loads(state).get("screen") in ("main", "start"):
                break
            time.sleep(0.5)
        result.update(ok=bool(state and state != "null"), state=json.loads(state) if state else None)
        _say("SELFTEST-READY")
        time.sleep(float(os.environ.get("KAMD_SELFTEST_HOLD", "1")))  # time for an external screenshot
        _say("SELFTEST " + json.dumps(result))
        if selftest and selftest != "-":  # windowed apps on Windows have no stdout: also write a file
            Path(selftest).write_text(json.dumps(result))
        window.destroy()

    kwargs = {"debug": debug}
    if sys.platform.startswith("linux") and not os.environ.get("PYWEBVIEW_GUI"):
        try:  # Qt WebEngine if available (bundled in the Linux build), else pywebview picks GTK
            import qtpy  # noqa: F401

            kwargs["gui"] = "qt"
        except ImportError:
            pass
    if selftest:
        kwargs["func"] = check
    webview.start(**kwargs)
    server.should_exit = True
    return 0 if (not selftest or result["ok"]) else 1


def _log_path() -> Path:
    from .recent import _path

    return _path().parent / "kamd-regions.log"


def main() -> int:  # entry point of the frozen app (PyInstaller)
    import argparse
    import faulthandler
    import logging

    # a windowed app has no console: send output and errors to a log file
    log = _log_path()
    try:
        log.parent.mkdir(parents=True, exist_ok=True)
        fh = open(log, "w", buffering=1, encoding="utf-8")  # noqa: SIM115 - lives for the whole run
        if sys.stdout is None:
            sys.stdout = fh
        if sys.stderr is None:
            sys.stderr = fh
        faulthandler.enable(fh)
        logging.basicConfig(stream=fh, level=logging.INFO)
    except OSError:
        pass

    ap = argparse.ArgumentParser(prog="KAMD Regions")
    ap.add_argument("files", nargs="*")
    ap.add_argument("--selftest", nargs="?", const="-", default=None)
    ap.add_argument("--debug", action="store_true")
    a, _ = ap.parse_known_args()
    # macOS passes -psn_... when started from Finder
    files = [f for f in a.files if not f.startswith("-psn")]
    return run(files or None, selftest=a.selftest, debug=a.debug)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
