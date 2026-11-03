"""Write a self-contained HTML snapshot of the GUI (``kamd-regions snapshot``).

The page embeds the built frontend (JS + CSS inline), the 3D structure and the results of the
suggestion methods, so it opens from a file or a web host without the Python server. Editing,
undo and validation run in the browser; export needs ``kamd-regions serve``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from .server import STATIC
from .session import Session


def _inline_assets(html: str) -> str:
    def css(m):
        return "<style>\n" + (STATIC / m.group(1)).read_text() + "\n</style>"

    def js(m):
        code = (STATIC / m.group(1)).read_text().replace("</script", "<\\/script")
        return '<script type="module">\n' + code + "\n</script>"

    html = re.sub(r'<link rel="stylesheet"[^>]*href="\./([^"]+)"[^>]*>', css, html)
    html = re.sub(r'<script type="module"[^>]*src="\./([^"]+)"[^>]*></script>', js, html)
    return html


def snapshot(
    session: Session, methods=("gnm", "contacts", "tm"), title: str | None = None, body_only: bool = False
) -> str:
    index = STATIC / "index.html"
    if not index.exists():
        raise FileNotFoundError("the frontend is not built (npm run build in tools/kamd-regions/frontend)")
    suggestions = {}
    for m in methods:
        if m == "tm" and session.an.membrane is None:
            continue
        p, notes = session.suggest(m, {})
        suggestions[m] = {"project": p.to_dict(), "notes": notes, "validation": session.validate(p.to_dict())}
    # file paths of the machine that made the snapshot mean nothing to the reader
    system = session.system_summary()
    system["files"] = {k: Path(v).name for k, v in system["files"].items()}
    project = session.project.to_dict()
    project["inputs"] = system["files"]
    for s in suggestions.values():
        s["project"]["inputs"] = system["files"]
    data = {
        "session": {"system": system, "project": project, "validation": session.validate(), "out_dir": "-"},
        "structure": session.viewer_structure(),
        "suggest": suggestions,
    }
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    html = _inline_assets(index.read_text())
    name = title or "KAMD Regions"
    data_tag = f"<script>window.__KAMD_DEMO__={payload};</script>"
    if body_only:
        # for hosts that wrap the page in their own document: title first, no html/head/body tags
        style = "".join(re.findall(r"<style>.*?</style>", html, re.S))
        script = "".join(re.findall(r'<script type="module">.*?</script>', html, re.S))
        return f'<title>{name}</title>\n{style}\n{data_tag}\n<div id="root"></div>\n{script}\n'
    html = re.sub(r"<title>.*?</title>", f"<title>{name}</title>", html)
    # keep the title near the top of the file; the data block follows it
    return html.replace("</title>", "</title>\n" + data_tag, 1)
