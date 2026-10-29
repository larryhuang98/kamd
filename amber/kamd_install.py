#!/usr/bin/env python3
"""Install KAMD into an Amber 26 pmemd source tree (and remove it again).

    python3 kamd_install.py --check     /path/to/amber26_src   # dry run, changes nothing
    python3 kamd_install.py             /path/to/amber26_src   # install
    python3 kamd_install.py --uninstall /path/to/amber26_src   # restore the original files

The Amber source tree is the directory that contains src/pmemd/src (the unpacked
Amber 26 tarball). If omitted, $AMBERHOME is used. After installing, rebuild Amber as
usual (cmake ... -DCUDA=TRUE; make install); KAMD runs in pmemd.cuda_SPFP/_DPFP.

How it works: kamd_src/manifest.json lists 25 edits. Each edit names a pmemd file, a few
lines of unchanged Amber code that identify one unique place ("anchor"), the Amber lines
to replace there (often none), and a KAMD code block from kamd_src/blocks/. One new file
(cuda/kKAMD.h) is added. Every edit is located and checked before any file is written;
if any anchor is missing or ambiguous nothing is changed. Originals are kept in
src/pmemd/src/.kamd_backup/ and the result is recorded in src/pmemd/src/.kamd_install.json.
Requires only the Python 3.6+ standard library.
"""
import sys
if sys.version_info < (3, 6):
    sys.exit("kamd_install: ERROR: Python 3.6 or newer is required (run it with python3)")
import argparse, datetime, hashlib, json, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "kamd_src")
ENC = "latin-1"                      # byte-exact round trip

def split_lines(text):
    """Split on '\n' only (str.splitlines also splits on \x85 etc., which occur inside
    UTF-8 characters such as the Angstrom sign when the file is read as latin-1)."""
    parts = text.split("\n")
    out = [p + "\n" for p in parts[:-1]]
    if parts[-1]:
        out.append(parts[-1])
    return out

PMEMD = os.path.join("src", "pmemd", "src")
STAMP = ".kamd_install.json"
BACKUP = ".kamd_backup"


def die(msg):
    sys.stdout.flush()
    sys.stderr.write("kamd_install: ERROR: %s\n" % msg)
    sys.exit(1)


def read(path):
    with open(path, encoding=ENC, newline="") as f:
        return f.read()


def write_atomic(path, text):
    tmp = path + ".kamd_tmp"
    with open(tmp, "w", encoding=ENC, newline="") as f:
        f.write(text)
    shutil.copymode(path, tmp) if os.path.exists(path) else None
    os.replace(tmp, path)


def sha(text):
    return hashlib.sha256(text.encode(ENC)).hexdigest()


def norm(line):
    return " ".join(line.split())


def find(lines, seq):
    """Start indices where seq occurs; exact first, then ignoring whitespace differences."""
    k = len(seq)
    hits = [i for i in range(len(lines) - k + 1) if lines[i:i + k] == seq]
    if hits:
        return hits, False
    ns = [norm(l) for l in seq]
    nl = [norm(l) for l in lines]
    return [i for i in range(len(lines) - k + 1) if nl[i:i + k] == ns], True


def plan(root, manifest):
    """Locate every edit; return {file: new_text}, report lines, number of problems."""
    out, report, bad = {}, [], 0
    byfile = {}
    for e in manifest["edits"]:
        byfile.setdefault(e["file"], []).append(e)
    for rel, edits in byfile.items():
        path = os.path.join(root, rel)
        if not os.path.isfile(path):
            report.append("  [FAIL] %s: file not found" % rel); bad += 1; continue
        lines = split_lines(read(path))
        located = []
        for e in edits:
            seq = e["before"] + e["old"] + e["after"]
            hits, fuzzy = find(lines, seq)
            tag = "edit %2d  %s" % (e["id"], rel.replace(PMEMD + "/", "").replace("src/pmemd/src/", ""))
            if len(hits) != 1:
                why = "anchor not found" if not hits else "anchor found %d times" % len(hits)
                report.append("  [FAIL] %s: %s" % (tag, why)); bad += 1; continue
            start = hits[0] + len(e["before"])
            block = read(os.path.join(DATA, e["block"]))
            located.append((start, start + len(e["old"]), block, e))
            report.append("  [ ok ] %s%s" % (tag, "  (matched ignoring whitespace)" if fuzzy else ""))
        located.sort(key=lambda t: t[0])
        for a, b in zip(located, located[1:]):
            if a[1] > b[0]:
                report.append("  [FAIL] edits %d and %d overlap" % (a[3]["id"], b[3]["id"])); bad += 1
        for s, t, block, _ in reversed(located):
            lines[s:t] = split_lines(block)
        out[rel] = "".join(lines)
    return out, report, bad


def main():
    ap = argparse.ArgumentParser(description="Install KAMD into an Amber 26 pmemd source tree.")
    ap.add_argument("amber_src", nargs="?", default=os.environ.get("AMBERHOME"),
                    help="Amber source tree containing src/pmemd/src (default: $AMBERHOME)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true", help="dry run: report what would be done")
    g.add_argument("--uninstall", action="store_true", help="restore the original pmemd files")
    ap.add_argument("--force", action="store_true",
                    help="with --uninstall: restore even if patched files were edited since")
    a = ap.parse_args()
    if not a.amber_src:
        die("give the Amber source directory (or set AMBERHOME)")
    root = os.path.abspath(a.amber_src)
    if not os.path.isdir(os.path.join(root, PMEMD)):
        die("%s does not contain %s; point to the unpacked Amber 26 source tree" % (root, PMEMD))
    manifest = json.load(open(os.path.join(DATA, "manifest.json")))
    pm = os.path.join(root, PMEMD)
    stamp_path = os.path.join(pm, STAMP)
    bdir = os.path.join(pm, BACKUP)
    print("KAMD %s installer for Amber %s pmemd" % (manifest["kamd_version"], manifest["amber"]))
    print("Amber source tree: %s" % root)

    if a.uninstall:
        if not os.path.isfile(stamp_path):
            die("KAMD is not installed here (no %s)" % os.path.join(PMEMD, STAMP))
        st = json.load(open(stamp_path))
        changed = [f for f, h in st["patched_sha256"].items()
                   if not os.path.isfile(os.path.join(root, f)) or sha(read(os.path.join(root, f))) != h]
        if changed and not a.force:
            die("these files changed after KAMD was installed (use --force to restore anyway, "
                "discarding those changes):\n  " + "\n  ".join(changed))
        for f in st["modified"]:
            shutil.copy2(os.path.join(bdir, f), os.path.join(root, f))
        for f in st["added"]:
            if os.path.exists(os.path.join(root, f)):
                os.remove(os.path.join(root, f))
        shutil.rmtree(bdir)
        os.remove(stamp_path)
        print("KAMD removed: %d files restored, %d removed. Rebuild Amber." % (
            len(st["modified"]), len(st["added"])))
        return

    if os.path.isfile(stamp_path):
        st = json.load(open(stamp_path))
        die("KAMD %s is already installed here (%s). Use --uninstall first." % (
            st.get("kamd_version"), st.get("date")))
    for nf in manifest["new_files"]:
        if os.path.exists(os.path.join(root, nf["file"])):
            die("%s already exists: this tree seems to contain KAMD already (e.g. from the old "
                "diff). Start from a clean Amber 26 source tree." % nf["file"])

    new_text, report, bad = plan(root, manifest)
    print("\n".join(report))
    for nf in manifest["new_files"]:
        print("  [ ok ] new file %s" % nf["file"].replace("src/pmemd/src/", ""))
    if bad:
        die("%d edit(s) could not be placed; nothing was changed. This installer supports "
            "Amber %s (pmemd from the Amber%s release tarball, with or without updates)."
            % (bad, manifest["amber"], manifest["amber"]))
    if a.check:
        print("Check passed: all %d edits can be applied. Nothing was changed." % len(manifest["edits"]))
        return

    os.makedirs(bdir)
    rec = {"kamd_version": manifest["kamd_version"], "amber": manifest["amber"],
           "date": datetime.datetime.now().isoformat(timespec="seconds"),
           "modified": sorted(new_text), "added": [], "patched_sha256": {}}
    for rel, text in sorted(new_text.items()):
        dst = os.path.join(bdir, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(os.path.join(root, rel), dst)
    for rel, text in sorted(new_text.items()):
        write_atomic(os.path.join(root, rel), text)
        rec["patched_sha256"][rel] = sha(text)
    for nf in manifest["new_files"]:
        text = read(os.path.join(DATA, nf["source"]))
        with open(os.path.join(root, nf["file"]), "w", encoding=ENC, newline="") as f:
            f.write(text)
        rec["added"].append(nf["file"])
        rec["patched_sha256"][nf["file"]] = sha(text)
    with open(stamp_path, "w") as f:
        json.dump(rec, f, indent=1)
    print("KAMD installed: %d files modified, %d added (originals in %s)." % (
        len(rec["modified"]), len(rec["added"]), os.path.join(PMEMD, BACKUP)))
    print("Next: build Amber with CUDA (see amber/README.md); KAMD runs in pmemd.cuda_SPFP/_DPFP.")


if __name__ == "__main__":
    main()
