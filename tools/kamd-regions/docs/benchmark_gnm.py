"""GNM domain recovery on multi-domain proteins with literature domains.

    python docs/benchmark_gnm.py            # writes the tables used in docs/gnm_multidomain.md

adk_open.pdb and t4l.pdb are in tests/data; the closed ADK and HhaI methyltransferase files are
downloaded from the MDAnalysis and ProDy test suites on first use.
"""

from __future__ import annotations

import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np

from kamd_regions.analysis import Analysis
from kamd_regions.compare import adjusted_rand
from kamd_regions.masks import residue_summary
from kamd_regions.methods import suggest
from kamd_regions.system import load_system

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "tests" / "data"
CACHE = HERE / ".benchmark_cache"
URLS = {
    "adk_closed.pdb": "https://raw.githubusercontent.com/MDAnalysis/mdanalysis/develop/testsuite/MDAnalysisTests/data/adk_closed.pdb",
    "3mht.pdb": "https://raw.githubusercontent.com/prody/ProDy/main/prody/tests/datafiles/pdb3mht.pdb",
}
# literature domains (residue numbers)
CASES = {
    "adk_open.pdb": (
        "ADK open (4AKE)",
        None,
        {"CORE": [(1, 29), (60, 121), (160, 214)], "NMP": [(30, 59)], "LID": [(122, 159)]},
    ),
    "adk_closed.pdb": (
        "ADK closed (1AKE)",
        None,
        {"CORE": [(1, 29), (60, 121), (160, 214)], "NMP": [(30, 59)], "LID": [(122, 159)]},
    ),
    "t4l.pdb": ("T4 lysozyme L99A", None, {"C-domain": [(1, 12), (60, 164)], "N-domain": [(13, 59)]}),
    "3mht.pdb": ("HhaI methyltransferase + DNA (3MHT)", "A", {"large": [(1, 193), (304, 327)], "small": [(194, 303)]}),
}
SYM = ".123456789"


def path_of(name: str) -> Path:
    p = DATA / name
    if p.exists():
        return p
    CACHE.mkdir(exist_ok=True)
    p = CACHE / name
    if not p.exists():
        urllib.request.urlretrieve(URLS[name], p)
    return p


def evaluate(an: Analysis, chain, ref, params):
    s = an.s
    prot = np.array([r.index for r in s.residues if r.kind == "protein" and (chain is None or r.chain == chain)])
    num = np.array([s.residues[i].number for i in prot])
    t = np.zeros(s.n_residues, dtype=int)
    for k, rr in enumerate(ref.values()):
        for a, b in rr:
            t[prot[(num >= a) & (num <= b)]] = k + 1
    res = suggest(an, "gnm", params)
    lab = res["labels"]
    both = (lab[prot] > 0) & (t[prot] > 0)
    ari = adjusted_rand(t[prot][both], lab[prot][both])
    rows = []
    for g in range(1, lab.max() + 1):
        m = lab == g
        k = int(np.bincount(t[m], minlength=len(ref) + 1)[1:].argmax()) + 1
        name = list(ref)[k - 1]
        rows.append(
            (
                g,
                residue_summary(s, np.flatnonzero(m & np.isin(np.arange(s.n_residues), prot))),
                int(s.res_natoms[m].sum()),
                name,
                float((t[m] == k).mean()),
                float((m & (t == k)).sum() / (t == k).sum()),
            )
        )
    bar_ref = "".join(SYM[t[i]] for i in prot)
    bar = "".join(SYM[min(lab[i], 9)] for i in prot)
    return res, ari, rows, bar_ref, bar, int((lab[prot] == 0).sum())


def main() -> int:
    warnings.filterwarnings("ignore")
    out = []
    for name, (title, chain, ref) in CASES.items():
        an = Analysis(load_system([str(path_of(name))]))
        out.append(f"### {title}\n")
        bars = []
        for label, params in (("defaults", {}), (f"n_regions={len(ref)}", {"n_regions": len(ref)})):
            res, ari, rows, bar_ref, bar, nun = evaluate(an, chain, ref, params)
            if not bars:
                bars.append(f"{'reference':>12s}  {bar_ref}")
            bars.append(f"{label:>12s}  {bar}")
            out.append(
                f"**{label}**: {len(rows)} regions (k = {res['info']['k']}), ARI {ari:.2f}, "
                f"{nun} residues unassigned.\n"
            )
            out.append("| region | residues | atoms | best domain | precision | recall |\n|---|---|---|---|---|---|")
            for g, summ, n, dom, prec, rec in rows:
                out.append(f"| {g} | {summ} | {n} | {dom} | {prec:.0%} | {rec:.0%} |")
            out.append("")
        out.append("```\n" + "\n".join(bars) + "\n```\n")
    text = "\n".join(out)
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
