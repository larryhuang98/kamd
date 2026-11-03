"""Live validation of a region set against the KAMD rules of each engine.

Overlap cannot occur: a residue belongs to at most one region by construction.
"""

from __future__ import annotations

import numpy as np

from .analysis import Analysis
from .masks import residue_summary
from .membrane import Membrane, runs
from .project import Project

ENGINES = ("amber", "gromacs", "openmm")
AMBER_MAX_GROUPS = 32


def effective_membrane(an: Analysis, project: Project) -> Membrane | None:
    m = project.settings.membrane
    if m:
        return Membrane(int(m["axis"]), float(m["lo"]), float(m["hi"]), "manual")
    return an.membrane


def region_stats(an: Analysis, project: Project) -> dict[int, dict]:
    s = an.s
    mem = effective_membrane(an, project)
    inside = mem.inside(an.rep_xyz) if mem is not None else np.zeros(s.n_residues, dtype=bool)
    out = {}
    for r in project.regions:
        res = project.residues_of(r.id)
        out[r.id] = {
            "n_residues": int(res.size),
            "n_atoms": int(s.res_natoms[res].sum()),
            "n_heavy": int(s.res_nheavy[res].sum()),
            "summary": residue_summary(s, res),
            "frac_in_slab": float(inside[res].mean()) if res.size and mem is not None else None,
            "issues": [],
        }
    return out


def validate(an: Analysis, project: Project, engines=ENGINES) -> dict:
    s = an.s
    errors: list[dict] = []
    warnings: list[dict] = []
    stats = region_stats(an, project)

    def err(code, msg, eng=ENGINES, region=None):
        errors.append({"code": code, "message": msg, "engines": list(eng), "region": region})
        if region is not None:
            stats[region]["issues"].append({"level": "error", "message": msg})

    def warn(code, msg, eng=ENGINES, region=None):
        warnings.append({"code": code, "message": msg, "engines": list(eng), "region": region})
        if region is not None:
            stats[region]["issues"].append({"level": "warning", "message": msg})

    st = project.settings
    regions = project.regions
    nonempty = [r for r in regions if stats[r.id]["n_residues"] > 0]
    for r in regions:
        if stats[r.id]["n_residues"] == 0:
            err("empty_group", f"{r.name} is empty (no residues)", region=r.id)
    if len(nonempty) < 2:
        err("too_few_groups", f"KAMD needs at least 2 reference groups; {len(nonempty)} defined")
    if len(regions) > AMBER_MAX_GROUPS:
        err("too_many_groups", f"Amber allows at most {AMBER_MAX_GROUPS} groups; {len(regions)} defined",
            eng=("amber",))

    bad = np.flatnonzero((project.assignment > 0) & ~an.assignable)
    if bad.size:
        err("bath_in_group", f"{bad.size} water/ion/lipid residues are assigned to a group (e.g. "
            f"{s.residues[int(bad[0])].label()}); groups and bath must not overlap")
    unknown = set(np.unique(project.assignment[project.assignment > 0]).tolist()) - {r.id for r in regions}
    if unknown:
        err("unknown_region", f"residues assigned to unknown region ids {sorted(unknown)}")

    # sizes
    for r in nonempty:
        n = stats[r.id]["n_atoms"]
        if n < st.small_atoms:
            warn("small_group", f"{r.name} has only {n} atoms; resets of small groups are noisy "
                 f"(a few hundred to a few thousand atoms is typical)", region=r.id)
        if n > st.large_atoms:
            warn("large_group", f"{r.name} has {n} atoms; very large groups average out relative motion",
                 region=r.id)

    # bath
    bath = an.bath
    if bath.solvent_atoms.size == 0:
        err("no_solvent", "no water or ions found: the solvent bath is required. Load the solvated system "
            "(prmtop + coordinates, or the solvated PDB/gro).", eng=("gromacs", "openmm"))
        warn("no_solvent_amber", "no water or ions in this structure: the Amber input uses the placeholder "
             "solvent_mask = \":WAT,Na+,Cl-\"; check it against your prmtop", eng=("amber",))
    has_tm = any(r.region_type == 2 for r in nonempty)
    mem = effective_membrane(an, project)
    if has_tm:
        if bath.membrane_atoms.size == 0:
            err("no_membrane_bath", "transmembrane groups (region_type 2) need a membrane bath, but the system "
                "has no lipid residues" + (" (the OPM slab is used only for the suggestion)" if mem else ""))
        if mem is None:
            err("no_membrane_axis", "transmembrane groups need membrane_axis, but no membrane was found; set the "
                "slab manually")
    for problem in bath.contiguity(s):
        if problem["bath"] == "membrane" and not has_tm:
            continue
        err(f"{problem['bath']}_not_contiguous", problem["message"], eng=("amber",))
    if bath.membrane_atoms.size and bath.solvent_atoms.size:
        pass  # disjoint by construction (residue kinds)

    # membrane placement
    if mem is not None:
        for r in nonempty:
            f = stats[r.id]["frac_in_slab"]
            if r.region_type == 2 and f is not None and f < 0.5:
                warn("tm_outside", f"{r.name} is transmembrane but only {f:.0%} of its residues lie in the "
                     "membrane slab", region=r.id)
            if r.region_type == 1 and f is not None and f > 0.6 and bath.membrane_atoms.size:
                warn("solution_in_membrane", f"{r.name} is a solution group but {f:.0%} of its residues are buried "
                     "in the membrane; consider region_type 2", region=r.id)

    # refer_link
    pos = project.group_position()
    ids = {r.id for r in regions}
    links = {r.id: r.refer_link for r in regions if r.refer_link}
    for rid, k in links.items():
        if k not in ids:
            err("bad_link", f"{project.region(rid).name}: refer_link points to a deleted region", region=rid)
        elif k == rid:
            err("self_link", f"{project.region(rid).name} cannot link to itself", region=rid)
    if links:
        targets = {k for k in links.values() if k in ids}
        free = [r.id for r in nonempty if r.id not in links]
        if len(targets) == 1 and len(free) <= 1:
            warn("all_linked", "every group shares one direction through refer_link: the reset then pushes "
                 "all groups the same way and no relative motion is excited")
        for rid, k in links.items():
            if k in links and k in ids:
                warn("link_chain", f"{project.region(rid).name} links to {project.region(k).name}, which itself "
                     "links to another group; KAMD applies links in group order, so the result depends on order",
                     region=rid)
        if len(targets) == 1 and len(links) == 1 and len(nonempty) > 2:
            pass

    # numbering
    if not s.has_topology:
        warn("no_topology", "No topology loaded: residue numbers come from the structure file order. The Amber "
             "masks may not match your prmtop (and atom numbers may not match a GROMACS .tpr). Load the "
             "prmtop/top/tpr for exact export.", eng=("amber", "gromacs", "openmm"))

    ok = {e: not any(e in x["engines"] for x in errors) for e in engines}
    return {
        "errors": errors,
        "warnings": warnings,
        "regions": {str(k): v for k, v in stats.items()},
        "ok": ok,
        "group_position": {str(k): v for k, v in pos.items()},
        "n_unassigned_protein": int(((project.assignment == 0) & an.assignable).sum()),
    }


def bath_ranges(atoms: np.ndarray) -> list[tuple[int, int]]:
    return runs(atoms)
