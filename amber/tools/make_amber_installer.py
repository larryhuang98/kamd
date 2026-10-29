#!/usr/bin/env python3
"""Maintainer tool: convert the tested KAMD pmemd diff into the data used by kamd_install.py.

  make_amber_installer.py PATCH OUTDIR TREE [TREE ...]

PATCH  unified diff (git diff) of src/pmemd/src, KAMD branch vs its Amber base
OUTDIR receives manifest.json, new files under files/, inserted code under blocks/
TREE   pristine Amber source roots; every anchor must be unique in every tree

Each change (a run of -/+ lines) becomes one edit: {file, before, old, after, block}.
'before'/'after' are the shortest unchanged context that makes before+old+after unique.
"""
import json, os, re, sys

ENC = "latin-1"   # byte-exact round trip for any source file

def split_lines(text):
    """Split on '\n' only (str.splitlines also splits on \x85 etc., which occur inside
    UTF-8 characters such as the Angstrom sign when the file is read as latin-1)."""
    parts = text.split("\n")
    out = [p + "\n" for p in parts[:-1]]
    if parts[-1]:
        out.append(parts[-1])
    return out



def read_lines(path):
    with open(path, encoding=ENC, newline="") as f:
        return split_lines(f.read())


def parse(patch):
    files, cur, hunk = [], None, None
    for line in split_lines(open(patch, encoding=ENC, newline="").read()):
        if line.startswith("diff --git"):
            cur = {"path": line.split()[3][2:], "new": False, "hunks": []}
            files.append(cur)
        elif line.startswith("new file mode"):
            cur["new"] = True
        elif line.startswith("@@"):
            m = re.match(r"@@ -(\d+)(?:,(\d+))? \+(\d+)", line)
            hunk = {"old_start": int(m.group(1)), "lines": []}
            cur["hunks"].append(hunk)
        elif hunk is not None and line[:1] in " +-" and not line.startswith(("--- ", "+++ ")):
            hunk["lines"].append(line)
        elif line.startswith("\\"):
            # "\ No newline at end of file": strip the newline from the previous line
            prev = hunk["lines"][-1]
            hunk["lines"][-1] = prev.rstrip("\r\n")
    return files


def count(lines, seq):
    n, k = 0, len(seq)
    for i in range(len(lines) - k + 1):
        if lines[i:i + k] == seq:
            n += 1
    return n


def main():
    patch, out, trees = sys.argv[1], sys.argv[2], sys.argv[3:]
    os.makedirs(os.path.join(out, "files"), exist_ok=True)
    os.makedirs(os.path.join(out, "blocks"), exist_ok=True)
    manifest = {"edits": [], "new_files": []}
    nblock = 0
    for f in parse(patch):
        rel = f["path"]
        if f["new"]:
            body = "".join(l[1:] for h in f["hunks"] for l in h["lines"])
            dst = os.path.join(out, "files", os.path.basename(rel))
            open(dst, "w", encoding=ENC, newline="").write(body)
            manifest["new_files"].append({"file": rel, "source": "files/" + os.path.basename(rel)})
            continue
        orig = read_lines(os.path.join(trees[0], rel))
        others = [read_lines(os.path.join(t, rel)) for t in trees]
        # rebuild the base file positions of every change from the hunks
        changes = []
        for h in f["hunks"]:
            # locate hunk in the first tree by its old-side text (offsets differ between trees)
            oldside = [l[1:] for l in h["lines"] if l[0] in " -"]
            starts = [i for i in range(len(orig) - len(oldside) + 1) if orig[i:i + len(oldside)] == oldside]
            if len(starts) != 1:
                sys.exit("hunk of %s not located uniquely (%d)" % (rel, len(starts)))
            pos, i, L = starts[0], 0, h["lines"]
            while i < len(L):
                if L[i][0] == " ":
                    pos += 1; i += 1; continue
                old, new = [], []
                while i < len(L) and L[i][0] == "-":
                    old.append(L[i][1:]); i += 1
                while i < len(L) and L[i][0] == "+":
                    new.append(L[i][1:]); i += 1
                changes.append((pos, old, new))
                pos += len(old)
        for pos, old, new in changes:
            nb = na = 0
            while True:
                before = orig[max(0, pos - nb):pos]
                after = orig[pos + len(old):pos + len(old) + na]
                seq = before + old + after
                if seq and all(count(t, seq) == 1 for t in others):
                    break
                if nb < 12 and pos - nb > 0:
                    nb += 1
                elif na < 6:
                    na += 1
                else:
                    sys.exit("no unique anchor for change at %s:%d" % (rel, pos + 1))
            nblock += 1
            ext = os.path.basename(rel)
            bname = "blocks/%02d_%s.txt" % (nblock, ext)
            open(os.path.join(out, bname), "w", encoding=ENC, newline="").write("".join(new))
            manifest["edits"].append({"id": nblock, "file": rel, "before": before, "old": old,
                                      "after": after, "block": bname})
    json.dump(manifest, open(os.path.join(out, "manifest.json"), "w"), indent=1)
    print("edits: %d, new files: %d" % (len(manifest["edits"]), len(manifest["new_files"])))
    for e in manifest["edits"]:
        print("  %2d %-40s before=%d old=%d after=%d new=%d" % (
            e["id"], e["file"].replace("src/pmemd/src/", ""), len(e["before"]), len(e["old"]),
            len(e["after"]), len(open(os.path.join(out, e["block"]), encoding=ENC).read().splitlines())))


if __name__ == "__main__":
    main()
