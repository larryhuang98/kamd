"""Import existing KAMD inputs (Amber &kamd, GROMACS kamd.dat + index.ndx, OpenMM JSON)."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from .analysis import Analysis
from .project import Project, Region, Settings, default_color


def _groups_to_project(
    an: Analysis,
    groups: list[np.ndarray],
    names: list[str],
    rtype: list[int],
    links: list[int],
    settings: Settings,
    source: str,
) -> tuple[Project, list[str]]:
    s = an.s
    notes: list[str] = []
    assignment = np.zeros(s.n_residues, dtype=np.int64)
    regions = []
    for i, g in enumerate(groups):
        rid = i + 1
        g = np.unique(np.asarray(g, dtype=np.int64))
        if g.size and (g.min() < 0 or g.max() >= s.n_atoms):
            raise ValueError(f"group {rid} has atom indices outside the {s.n_atoms}-atom system")
        res, counts = np.unique(s.atom_res[g], return_counts=True)
        for r, c in zip(res, counts, strict=True):
            full = s.res_natoms[r]
            if c < full:
                rr = s.residues[int(r)]
                keep = c * 2 >= full
                notes.append(
                    f"group {rid}: the mask covers {c} of {full} atoms of {rr.label()}; regions are whole residues, "
                    f"so the residue is {'included' if keep else 'left out'}"
                )
                if not keep:
                    continue
            if assignment[r] and assignment[r] != rid:
                notes.append(
                    f"{s.residues[int(r)].label()} is in groups {assignment[r]} and {rid}; kept in {assignment[r]}"
                )
                continue
            assignment[r] = rid
        regions.append(
            Region(
                rid,
                names[i] if i < len(names) else f"group {rid}",
                default_color(i),
                int(rtype[i]) if i < len(rtype) else 1,
                int(links[i]) if i < len(links) else 0,
            )
        )
    bad = np.flatnonzero((assignment > 0) & ~an.assignable)
    if bad.size:
        notes.append(f"{bad.size} residues of the imported groups are water/ion/lipid by name; they stay in the group")
    p = Project(
        n_residues=s.n_residues,
        regions=regions,
        assignment=assignment,
        settings=settings,
        method={"method": "import", "source": source},
        residue_names=[r.name for r in s.residues],
    )
    return p, notes


def _bath_note(an: Analysis, name: str, atoms: np.ndarray | None, notes: list[str]) -> None:
    if atoms is None:
        return
    want = an.bath.solvent_atoms if name == "solvent" else an.bath.membrane_atoms
    if not np.array_equal(np.unique(atoms), want):
        notes.append(
            f"the imported {name} bath ({atoms.size} atoms) differs from the detected one ({want.size} atoms); "
            "the export writes the detected bath"
        )


def import_amber(an: Analysis, mdin_text: str) -> tuple[Project, list[str]]:
    from kamd_openmm.amber_input import amber_mask_indices, parse_namelist

    st = an.s.structure
    if st is None:
        raise ValueError("Amber masks need a structure loaded with ParmEd")
    kamd = parse_namelist(mdin_text, "kamd")
    try:
        cntrl = parse_namelist(mdin_text, "cntrl")
    except KeyError:
        cntrl = {}
    raw = kamd.get("refer_mask", {})
    raw = raw if isinstance(raw, dict) else {1: raw}
    if sorted(raw) != list(range(1, len(raw) + 1)):
        raise ValueError("refer_mask entries must be packed (1..N without gaps)")
    masks = [str(raw[i]) for i in range(1, len(raw) + 1)]
    groups = [amber_mask_indices(st, m) for m in masks]

    def per(key, default):
        d = kamd.get(key, {})
        d = d if isinstance(d, dict) else {1: d}
        return [int(d.get(i, default)) for i in range(1, len(masks) + 1)]

    settings = Settings(
        pseudobd_cycles=int(kamd.get("pseudobd_cycles", 250)) or 250,
        nzero_crossings=int(kamd.get("nzero_crossings", 2)),
        shell_cut=float(kamd.get("shell_cut", 5.0)),
        surface_cut=float(kamd.get("surface_cut", 3.0)),
        membrane_cut=float(kamd.get("membrane_cut", 3.0)),
        kamd_shell_mode=int(kamd["kamd_shell_mode"]) if "kamd_shell_mode" in kamd else None,
        temperature=float(cntrl.get("temp0", 300.0)),
        seed=int(cntrl["ig"]) if "ig" in cntrl and int(cntrl["ig"]) > 0 else None,
    )
    p, notes = _groups_to_project(
        an,
        groups,
        [f"group {i + 1}" for i in range(len(groups))],
        per("region_type", 1),
        per("refer_link", 0),
        settings,
        "amber",
    )
    p.method["masks"] = masks
    if kamd.get("solvent_mask"):
        _bath_note(an, "solvent", amber_mask_indices(st, str(kamd["solvent_mask"])), notes)
    if kamd.get("membrane_mask"):
        _bath_note(an, "membrane", amber_mask_indices(st, str(kamd["membrane_mask"])), notes)
    if kamd.get("membrane_axis") and an.membrane is not None and int(kamd["membrane_axis"]) != an.membrane.axis:
        notes.append(f"membrane_axis = {kamd['membrane_axis']} differs from the detected axis {an.membrane.axis}")
    return p, notes


def _read_dat(text: str) -> dict[str, str]:
    out = {}
    for line in text.splitlines():
        line = re.split(r"[;#]", line, maxsplit=1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip().lower().replace("_", "-")] = v.strip()
    return out


def import_gromacs(an: Analysis, dat_text: str, ndx_text: str) -> tuple[Project, list[str]]:
    from .export import read_ndx

    d = _read_dat(dat_text)
    nd = read_ndx(ndx_text)
    names = d.get("refer-groups", "").split()
    if not names:
        raise ValueError("kamd.dat has no refer-groups")
    missing = [n for n in names if n not in nd]
    if missing:
        raise ValueError(f"groups {missing} are not in the index file")
    groups = [nd[n] for n in names]
    n = len(names)

    def ints(key, default):
        v = d.get(key)
        return [int(x) for x in v.split()] if v else [default] * n

    settings = Settings(
        pseudobd_cycles=int(d.get("pseudobd-cycles", 250)),
        nzero_crossings=int(d.get("nzero-crossings", 2)),
        shell_cut=float(d.get("shell-cut", 0.5)) * 10,
        surface_cut=float(d.get("surface-cut", 0.3)) * 10,
        membrane_cut=float(d.get("membrane-cut", 0.3)) * 10,
        temperature=float(d.get("temperature", 300.0)),
        seed=int(d["seed"]) if d.get("seed") and int(d["seed"]) >= 0 else None,
    )
    p, notes = _groups_to_project(an, groups, names, ints("region-type", 1), ints("refer-link", 0), settings, "gromacs")
    if d.get("solvent-group") in nd:
        _bath_note(an, "solvent", nd[d["solvent-group"]], notes)
    if d.get("membrane-group") in nd:
        _bath_note(an, "membrane", nd[d["membrane-group"]], notes)
    return p, notes


def import_openmm(an: Analysis, text: str) -> tuple[Project, list[str]]:
    d = json.loads(text)
    groups = [np.asarray(g) for g in d["refer_groups"]]
    n = len(groups)
    settings = Settings(
        pseudobd_cycles=int(d.get("pseudobd_cycles", 250)),
        nzero_crossings=int(d.get("nzero_crossings", 2)),
        shell_cut=float(d.get("shell_cut", 5.0)),
        surface_cut=float(d.get("surface_cut", 3.0)),
        membrane_cut=float(d.get("membrane_cut", 3.0)),
        temperature=float(d.get("temperature", 300.0)),
        seed=d.get("seed"),
    )
    names = d.get("group_names") or [f"group {i + 1}" for i in range(n)]
    p, notes = _groups_to_project(
        an, groups, names, d.get("region_type") or [1] * n, d.get("refer_link") or [0] * n, settings, "openmm"
    )
    if "solvent_atoms" in d:
        _bath_note(an, "solvent", np.asarray(d["solvent_atoms"]), notes)
    return p, notes


def import_file(an: Analysis, path: str, ndx: str | None = None) -> tuple[Project, list[str]]:
    text = Path(path).read_text()
    if "&kamd" in text.lower():
        return import_amber(an, text)
    if path.endswith(".json"):
        return import_openmm(an, text)
    if ndx is None:
        cand = Path(path).with_name("index.ndx")
        d = _read_dat(text)
        if d.get("index"):
            cand = Path(path).parent / d["index"]
        ndx = str(cand)
    return import_gromacs(an, text, Path(ndx).read_text())
