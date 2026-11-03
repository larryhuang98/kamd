"""Command line: ``kamd-regions serve|suggest|validate|export|compare``."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _params(items: list[str] | None) -> dict:
    out: dict = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"--param expects key=value, got {it!r}")
        k, v = it.split("=", 1)
        try:
            out[k] = json.loads(v)
        except json.JSONDecodeError:
            out[k] = v
    return out


def cmd_serve(a) -> int:
    import threading
    import webbrowser

    import uvicorn

    from .server import create_app
    from .session import Session

    print("loading", " ".join(a.files), "...", flush=True)
    s = Session(a.files, a.method, _params(a.param))
    for n in s.notes:
        print("note:", n)
    app = create_app(s, a.out_dir)
    url = f"http://{a.host}:{a.port}/"
    print(f"kamd-regions GUI at {url}  (Ctrl+C to stop)", flush=True)
    if not a.no_browser:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    uvicorn.run(app, host=a.host, port=a.port, log_level="warning")
    return 0


def _print_validation(v: dict) -> None:
    for e in v["errors"]:
        print(f"ERROR   [{','.join(e['engines'])}] {e['message']}")
    for w in v["warnings"]:
        print(f"warning [{','.join(w['engines'])}] {w['message']}")


def cmd_suggest(a) -> int:
    from .session import Session

    s = Session(a.files, a.method, _params(a.param))
    if any(k in s.roles for k in ("project", "mdin", "kamd_dat", "openmm_json")):
        # an input file was imported; still run the requested method
        s.project, notes = s.suggest(a.method, _params(a.param))
    else:
        notes = []
    p = s.project
    for n in s.notes + notes:
        print("note:", n)
    v = s.validate(p.to_dict())
    print(f"method {p.method.get('method')}: {len(p.regions)} regions")
    for r in p.regions:
        st = v["regions"][str(r.id)]
        kind = "TM " if r.region_type == 2 else "sol"
        print(f"  {r.id:3d} {kind} {st['n_atoms']:7d} atoms  {st['summary']}")
    _print_validation(v)
    if a.out:
        s.save(p.to_dict(), a.out)
        print("wrote", a.out)
    if a.export_dir:
        res = s.export(p.to_dict(), tuple(a.engines.split(",")), out_dir=a.export_dir, force=a.force)
        if not res["ok"]:
            print("export blocked by validation errors for:", ", ".join(res["blocked"]))
            return 2
        _print_verification(res["verification"])
        print("wrote", ", ".join(res["written"]))
    return 0


def _print_verification(rep: dict) -> None:
    print("verification:", rep["kamd_config"]["message"])
    for eng in ("amber", "gromacs", "openmm"):
        bad = [x for x in rep[eng] if not x["ok"]]
        if rep[eng]:
            print(
                f"  {eng}: {len(rep[eng]) - len(bad)}/{len(rep[eng])} items match"
                + ("" if not bad else " -- MISMATCH: " + ", ".join(x["item"] for x in bad))
            )


def cmd_validate(a) -> int:
    from .session import Session

    s = Session(a.files)
    v = s.validate()
    _print_validation(v)
    return 1 if v["errors"] else 0


def cmd_export(a) -> int:
    from .session import Session

    s = Session(a.files)
    res = s.export(None, tuple(a.engines.split(",")), out_dir=a.out_dir, force=a.force)
    _print_validation(res["validation"])
    if not res["ok"]:
        print("export blocked by validation errors for:", ", ".join(res["blocked"]))
        return 2
    _print_verification(res["verification"])
    print("wrote", ", ".join(res["written"]))
    return 0 if res["verification"]["ok"] else 1


def cmd_compare(a) -> int:
    from .compare import compare_report
    from .session import Session

    s = Session(a.files)
    text = compare_report(s, a.methods.split(","), _params(a.param))
    if a.out:
        Path(a.out).write_text(text)
        print("wrote", a.out)
    else:
        print(text)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="kamd-regions", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    files_help = (
        "structure / topology / coordinate files (PDB, mmCIF, gro, prmtop + inpcrd/rst7, top/tpr), "
        "optionally an existing &kamd mdin, kamd.dat + index.ndx, or a project.json"
    )

    p = sub.add_parser("serve", help="start the GUI")
    p.add_argument("files", nargs="+", help=files_help)
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--no-browser", action="store_true")
    p.add_argument("--method", default="auto")
    p.add_argument("--param", action="append", help="method parameter key=value (repeatable)")
    p.add_argument("--out-dir", default=None, help="where 'save'/'export to disk' write (default ./kamd_regions_out)")
    p.set_defaults(fn=cmd_serve)

    p = sub.add_parser("suggest", help="suggest regions without the GUI")
    p.add_argument("files", nargs="+", help=files_help)
    p.add_argument("--method", default="auto", help="auto, gnm, contacts or tm")
    p.add_argument("--param", action="append", help="method parameter key=value (repeatable)")
    p.add_argument("--out", help="write project.json here")
    p.add_argument("--export-dir", help="also export engine inputs here")
    p.add_argument("--engines", default="amber,gromacs,openmm")
    p.add_argument("--force", action="store_true", help="export even with validation errors")
    p.set_defaults(fn=cmd_suggest)

    p = sub.add_parser("validate", help="validate a project or imported input")
    p.add_argument("files", nargs="+", help=files_help)
    p.set_defaults(fn=cmd_validate)

    p = sub.add_parser("export", help="export a project.json / imported input to engine inputs")
    p.add_argument("files", nargs="+", help=files_help)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--engines", default="amber,gromacs,openmm")
    p.add_argument("--force", action="store_true")
    p.set_defaults(fn=cmd_export)

    p = sub.add_parser("compare", help="compare automatic methods with the imported (manual) groups")
    p.add_argument("files", nargs="+", help=files_help + "; must include the manual groups (mdin, kamd.dat, ...)")
    p.add_argument("--methods", default="gnm,contacts")
    p.add_argument("--param", action="append")
    p.add_argument("--out")
    p.set_defaults(fn=cmd_compare)

    a = ap.parse_args(argv)
    try:
        return a.fn(a)
    except (ValueError, FileNotFoundError) as exc:
        print("error:", exc, file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
