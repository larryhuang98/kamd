"""Atom-by-atom comparison of the first-shell search: kamd_openmm vs pmemd.cuda.

Inputs come from two pmemd.cuda runs started from the same restart:
  A) KAMD_DIAG=1 run: stops at the first reset step S and writes kamd_diag_gpu.txt
     (the GPU kernel's per-atom result for the state at completed step S-1);
  B) the same run with nstlim = S-1 and an ASCII restart: the coordinates and box
     that kernel saw.
"""

import argparse
import json

import numpy as np
from openmm import app, unit

from kamd_openmm.amber_input import config_from_amber
from kamd_openmm.core import find_first_shell


def load_gpu_diag(path):
    surface, assign = set(), {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        atom, rec, grp, _ = line.split()
        a = int(atom) - 1
        if int(rec) == -1:
            surface.add(a)
        elif int(grp) > 0:
            assign[a] = int(grp) - 1
    return surface, assign


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--mdin", required=True)
    ap.add_argument("--restart", required=True, help="ASCII restart at completed step S-1")
    ap.add_argument("--gpu-diag", required=True)
    ap.add_argument("--json", default="parity_first_shell.json")
    a = ap.parse_args()

    cfg = config_from_amber(a.mdin, a.prmtop).validate()
    crd = app.AmberInpcrdFile(a.restart)
    pos = np.asarray(crd.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    bv = crd.getBoxVectors()
    box = np.array(
        [
            [float(c) for c in (v.value_in_unit(unit.nanometer) if unit.is_quantity(v) else v)]
            for v in (bv.value_in_unit(unit.nanometer) if unit.is_quantity(bv) else bv)
        ]
    )

    res = find_first_shell(pos, box, cfg)
    res_bf = find_first_shell(pos, box, cfg, method="bruteforce")
    ours_assign = dict(zip(res.bath_atoms.tolist(), res.bath_group.tolist(), strict=False))
    ours_surf = dict(zip(res.surface_atoms.tolist(), res.surface_group.tolist(), strict=False))
    assert ours_assign == dict(zip(res_bf.bath_atoms.tolist(), res_bf.bath_group.tolist(), strict=False))
    assert ours_surf == dict(zip(res_bf.surface_atoms.tolist(), res_bf.surface_group.tolist(), strict=False))

    amb_surf, amb_assign = load_gpu_diag(a.gpu_diag)
    ref_group = {int(x): g for g, grp in enumerate(cfg.groups) for x in grp}

    ng = cfg.ngroups
    rows = []
    for g in range(ng):
        rows.append(
            {
                "group": g + 1,
                "shell_amber": sum(1 for v in amb_assign.values() if v == g),
                "shell_openmm": sum(1 for v in ours_assign.values() if v == g),
                "surface_amber": sum(1 for x in amb_surf if ref_group.get(x) == g),
                "surface_openmm": sum(1 for v in ours_surf.values() if v == g),
            }
        )
    shell_mis = sorted(set(amb_assign.items()) ^ set(ours_assign.items()))
    surf_mis = sorted(amb_surf ^ set(ours_surf))

    # distance of mismatched atoms to the nearest reference atom (to show whether they sit on the cutoff)
    allref = np.concatenate(cfg.groups)
    L = np.diag(box)

    def nearest(i, atoms):
        d = np.abs(pos[atoms] - pos[i])
        d = np.minimum(d, L - d)
        r = np.sqrt((d**2).sum(1)) * 10
        k = int(np.argmin(r))
        return float(r[k]), int(atoms[k])

    detail = []
    for atom, _ in shell_mis[:50]:
        r, j = nearest(atom, allref)
        detail.append(
            {
                "atom": atom + 1,
                "kind": "shell",
                "amber": amb_assign.get(atom, 0) + (atom in amb_assign),
                "openmm": ours_assign.get(atom, -1) + 1,
                "r_nearest_ref_A": round(r, 6),
            }
        )
    for atom in surf_mis[:50]:
        r, j = nearest(atom, cfg.solvent)
        detail.append(
            {
                "atom": atom + 1,
                "kind": "surface",
                "amber": atom in amb_surf,
                "openmm": atom in ours_surf,
                "r_nearest_solvent_A": round(r, 6),
            }
        )

    out = {
        "groups": rows,
        "shell_total_amber": len(amb_assign),
        "shell_total_openmm": len(ours_assign),
        "surface_total_amber": len(amb_surf),
        "surface_total_openmm": len(ours_surf),
        "shell_mismatches": len(set(a for a, _ in shell_mis)),
        "surface_mismatches": len(surf_mis),
        "mismatch_detail": detail,
    }
    print(f"{'group':>5} {'shell Amber':>12} {'shell OpenMM':>13} {'surf Amber':>11} {'surf OpenMM':>12}")
    for r in rows:
        print(
            f"{r['group']:>5} {r['shell_amber']:>12} {r['shell_openmm']:>13} "
            f"{r['surface_amber']:>11} {r['surface_openmm']:>12}"
        )
    print(f"total {len(amb_assign):>12} {len(ours_assign):>13} {len(amb_surf):>11} {len(ours_surf):>12}")
    print(f"mismatched shell atoms: {out['shell_mismatches']}, mismatched surface atoms: {out['surface_mismatches']}")
    for d in detail:
        print("  ", d)
    json.dump(out, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
