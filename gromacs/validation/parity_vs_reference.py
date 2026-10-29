"""Compare the GROMACS KAMD first-shell result (kamd_diag.txt) with kamd_openmm's search on the same
coordinates.  kamd_openmm itself matches pmemd.cuda atom by atom (openmm/validation/RESULTS.md)."""
import json
import sys

import numpy as np

from kamd_openmm.core import KAMDConfig, find_first_shell


def read_ndx(path):
    groups, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith("["):
            cur = line.strip("[] ")
            groups[cur] = []
        elif line and cur:
            groups[cur] += [int(t) - 1 for t in line.split()]
    return groups


diag, ndx = sys.argv[1], sys.argv[2]
lines = open(diag).read().splitlines()
step = lines[0].split()[-1]
box = np.array(list(map(float, lines[1].split()[1:]))).reshape(3, 3)
data = np.loadtxt(lines[3:])
pos = data[:, 1:4]
rec = data[:, 4].astype(int)

g = read_ndx(ndx)
names = ["dom1", "dom2", "dom3", "dom4"]
cfg = KAMDConfig(refer_groups=[g[n] for n in names], solvent_atoms=g["bath"], temperature=355.0).validate()
res = find_first_shell(pos, box, cfg)
ref_shell = dict(zip(res.bath_atoms.tolist(), res.bath_group.tolist()))
ref_surf = dict(zip(res.surface_atoms.tolist(), res.surface_group.tolist()))
gmx_shell = {i: r - 1 for i, r in enumerate(rec) if r > 0}
gmx_surf = {i: -r - 1 for i, r in enumerate(rec) if r < 0}

print(f"first reset at step {step}")
print(f"{'group':>5} {'shell GROMACS':>14} {'shell ref':>10} {'surf GROMACS':>13} {'surf ref':>9}")
for k in range(4):
    print(f"{k + 1:>5} {sum(v == k for v in gmx_shell.values()):>14} {sum(v == k for v in ref_shell.values()):>10}"
          f" {sum(v == k for v in gmx_surf.values()):>13} {sum(v == k for v in ref_surf.values()):>9}")
ms = set(gmx_shell.items()) ^ set(ref_shell.items())
mf = set(gmx_surf.items()) ^ set(ref_surf.items())
print(f"total {len(gmx_shell):>14} {len(ref_shell):>10} {len(gmx_surf):>13} {len(ref_surf):>9}")
print(f"mismatched shell entries: {len(ms)}, mismatched surface entries: {len(mf)}")
json.dump({"step": step, "shell": len(gmx_shell), "surface": len(gmx_surf),
           "shell_mismatch": len(ms), "surface_mismatch": len(mf)}, open("parity.json", "w"))
