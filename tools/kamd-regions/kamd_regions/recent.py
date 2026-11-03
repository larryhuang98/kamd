"""Recently opened file sets (desktop start screen), stored in the user's config directory."""

from __future__ import annotations

import json
import os
from pathlib import Path

MAX = 10


def _path() -> Path:
    base = os.environ.get("KAMD_REGIONS_CONFIG")
    if base:
        return Path(base) / "recent.json"
    if os.name == "nt":
        root = Path(os.environ.get("APPDATA", Path.home()))
    elif os.uname().sysname == "Darwin":
        root = Path.home() / "Library" / "Application Support"
    else:
        root = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return root / "kamd-regions" / "recent.json"


def recent_files() -> list[list[str]]:
    try:
        items = json.loads(_path().read_text())
    except (OSError, ValueError):
        return []
    return [x for x in items if isinstance(x, list) and all(os.path.exists(p) for p in x)][:MAX]


def add_recent(paths: list[str]) -> None:
    entry = [str(Path(p).resolve()) for p in paths]
    items = [x for x in recent_files() if x != entry]
    items.insert(0, entry)
    try:
        f = _path()
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps(items[:MAX], indent=1))
    except OSError:
        pass  # recent files are a convenience only
