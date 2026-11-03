"""Export to Amber (&kamd), GROMACS (index.ndx + kamd.dat) and OpenMM (JSON + example).

After writing, :func:`verify` re-reads the outputs: every Amber mask is evaluated with
ParmEd's AmberMask and compared atom by atom with the GUI groups, the index file is
parsed back, and the groups are run through ``kamd_openmm.KAMDConfig.validate``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from .analysis import Analysis
from .masks import atom_mask, residue_mask, residue_summary
from .membrane import AXIS_NAMES, runs
from .project import Project
from .validate import effective_membrane

PLACEHOLDER_SOLVENT = ":WAT,Na+,Cl-"


# --------------------------------------------------------------------------------------
# groups
# --------------------------------------------------------------------------------------
def group_atoms(an: Analysis, project: Project) -> list[np.ndarray]:
    return [an.s.atoms_of(project.residues_of(r.id)) for r in project.regions]


def bath_mask(an: Analysis, residues: np.ndarray, atoms: np.ndarray) -> str:
    """Residue range mask when the bath is whole contiguous residues, else atom ranges."""
    if residues.size and len(runs(residues)) == 1:
        return residue_mask(residues)
    return atom_mask(atoms)


def links_1based(project: Project) -> list[int]:
    pos = project.group_position()
    return [pos.get(r.refer_link, 0) if r.refer_link else 0 for r in project.regions]


def has_tm(project: Project) -> bool:
    return any(r.region_type == 2 for r in project.regions)


def gmx_names(project: Project) -> list[str]:
    seen: set[str] = set()
    out = []
    for i, r in enumerate(project.regions):
        n = re.sub(r"[^A-Za-z0-9_]+", "_", r.name).strip("_") or f"group{i + 1}"
        if n[0].isdigit():
            n = "g" + n
        base, k = n, 2
        while n in seen or n in ("KAMD_SOLVENT", "KAMD_MEMBRANE"):
            n = f"{base}_{k}"
            k += 1
        seen.add(n)
        out.append(n)
    return out


# --------------------------------------------------------------------------------------
# Amber
# --------------------------------------------------------------------------------------
def amber_namelist(an: Analysis, project: Project) -> str:
    s = an.s
    st = project.settings
    lines = []
    if not s.has_topology:
        lines += [
            "! WARNING: written without a topology. Residue numbers are the order of residues in",
            f"! {Path(s.files.get('structure', '?')).name}; they may NOT match your prmtop. Check every mask",
            "! (e.g. with cpptraj 'mask' or parmed) before running.",
        ]
    lines.append("&kamd")
    lines.append(f"   pseudobd_cycles = {st.pseudobd_cycles},")
    lines.append(f"   nzero_crossings = {st.nzero_crossings},")
    lines.append(f"   shell_cut       = {st.shell_cut},")
    lines.append(f"   surface_cut     = {st.surface_cut},")
    links = links_1based(project)
    for i, r in enumerate(project.regions):
        res = project.residues_of(r.id)
        lines.append(f'   refer_mask({i + 1}) = "{residue_mask(res)}",   ! {r.name}: {residue_summary(s, res)}')
    for i, r in enumerate(project.regions):
        if r.region_type != 1:
            lines.append(f"   region_type({i + 1}) = {r.region_type},")
    for i, k in enumerate(links):
        if k:
            lines.append(f"   refer_link({i + 1}) = {k},")
    b = an.bath
    if b.solvent_atoms.size:
        names = sorted({s.residues[i].name for i in b.solvent_residues})
        lines.append(f'   solvent_mask    = "{bath_mask(an, b.solvent_residues, b.solvent_atoms)}",'
                     f"   ! {', '.join(names[:8])}")
    else:
        lines.append(f'   solvent_mask    = "{PLACEHOLDER_SOLVENT}",   ! PLACEHOLDER: no solvent in the structure')
    if has_tm(project):
        mem = effective_membrane(an, project)
        if b.membrane_atoms.size:
            names = sorted({s.residues[i].name for i in b.membrane_residues})
            lines.append(f'   membrane_mask   = "{bath_mask(an, b.membrane_residues, b.membrane_atoms)}",'
                         f"   ! {', '.join(names[:8])}")
        if mem is not None:
            lines.append(f"   membrane_axis   = {mem.axis},")
        lines.append(f"   membrane_cut    = {st.membrane_cut},")
        mode = st.kamd_shell_mode if st.kamd_shell_mode is not None else 2
        lines.append(f"   kamd_shell_mode = {mode},")
    elif st.kamd_shell_mode is not None:
        lines.append(f"   kamd_shell_mode = {st.kamd_shell_mode},")
    lines.append("/")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------------------
# GROMACS
# --------------------------------------------------------------------------------------
def ndx_group(name: str, atoms0: np.ndarray) -> str:
    nums = (np.asarray(atoms0) + 1).tolist()
    rows = [" ".join(str(x) for x in nums[i : i + 15]) for i in range(0, len(nums), 15)]
    return f"[ {name} ]\n" + "\n".join(rows) + "\n"


def gromacs_files(an: Analysis, project: Project) -> dict[str, str]:
    st = project.settings
    names = gmx_names(project)
    groups = group_atoms(an, project)
    b = an.bath
    ndx = []
    for n, g in zip(names, groups, strict=True):
        ndx.append(ndx_group(n, g))
    if b.solvent_atoms.size:
        ndx.append(ndx_group("KAMD_SOLVENT", b.solvent_atoms))
    tm = has_tm(project)
    if tm and b.membrane_atoms.size:
        ndx.append(ndx_group("KAMD_MEMBRANE", b.membrane_atoms))
    links = links_1based(project)
    d = ["; KAMD input for gmx mdrun -kamd kamd.dat (written by kamd-regions)"]
    if not an.s.has_topology:
        d.append("; WARNING: written without a topology; atom numbers follow the structure file and may not "
                 "match your .tpr")
    d += [
        "index           = index.ndx",
        f"refer-groups    = {' '.join(names)}",
        "solvent-group   = KAMD_SOLVENT",
        f"pseudobd-cycles = {st.pseudobd_cycles}",
        f"nzero-crossings = {st.nzero_crossings}",
        f"shell-cut       = {st.shell_cut / 10:g}     ; nm",
        f"surface-cut     = {st.surface_cut / 10:g}     ; nm",
    ]
    if any(links):
        d.append(f"refer-link      = {' '.join(str(k) for k in links)}")
    if tm:
        mem = effective_membrane(an, project)
        d.append(f"region-type     = {' '.join(str(r.region_type) for r in project.regions)}")
        d.append("membrane-group  = KAMD_MEMBRANE")
        if mem is not None:
            d.append(f"membrane-axis   = {AXIS_NAMES[mem.axis]}")
        d.append(f"membrane-cut    = {st.membrane_cut / 10:g}     ; nm")
    d.append(f"temperature     = {st.temperature:g}    ; must equal ref-t")
    if st.seed is not None:
        d.append(f"seed            = {st.seed}")
    d.append("log             = kamd.log")
    return {"index.ndx": "".join(ndx), "kamd.dat": "\n".join(d) + "\n"}


def read_ndx(text: str) -> dict[str, np.ndarray]:
    out: dict[str, list[int]] = {}
    cur = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("["):
            cur = line.strip("[] ").strip()
            out[cur] = []
        elif line and cur is not None:
            out[cur].extend(int(x) for x in line.split())
    return {k: np.asarray(v, dtype=np.int64) - 1 for k, v in out.items()}


# --------------------------------------------------------------------------------------
# OpenMM
# --------------------------------------------------------------------------------------
def openmm_dict(an: Analysis, project: Project) -> dict:
    st = project.settings
    b = an.bath
    tm = has_tm(project)
    mem = effective_membrane(an, project)
    d = {
        "comment": "KAMDConfig input written by kamd-regions; atom indices are 0-based",
        "group_names": [r.name for r in project.regions],
        "refer_groups": [g.tolist() for g in group_atoms(an, project)],
        "solvent_atoms": b.solvent_atoms.tolist(),
        "temperature": st.temperature,
        "pseudobd_cycles": st.pseudobd_cycles,
        "nzero_crossings": st.nzero_crossings,
        "shell_cut": st.shell_cut,
        "surface_cut": st.surface_cut,
        "refer_link": links_1based(project),
        "region_type": [r.region_type for r in project.regions],
        "seed": st.seed,
        "n_atoms": an.s.n_atoms,
    }
    if tm:
        d["membrane_atoms"] = b.membrane_atoms.tolist()
        d["membrane_axis"] = mem.axis if mem is not None else 0
        d["membrane_cut"] = st.membrane_cut
    return d


OPENMM_EXAMPLE = '''"""Attach KAMD to an OpenMM simulation using the groups chosen in kamd-regions."""
import json

from kamd_openmm import KAMDConfig, KAMDReporter

with open("kamd_openmm.json") as fh:
    d = json.load(fh)

cfg = KAMDConfig(
    refer_groups=d["refer_groups"],          # 0-based atom indices, one list per group
    solvent_atoms=d["solvent_atoms"],
    temperature=d["temperature"],            # must equal the thermostat temperature
    pseudobd_cycles=d["pseudobd_cycles"],
    nzero_crossings=d["nzero_crossings"],
    shell_cut=d["shell_cut"],                # Angstrom
    surface_cut=d["surface_cut"],
    refer_link=d["refer_link"],
    region_type=d["region_type"],
    membrane_atoms=d.get("membrane_atoms"),
    membrane_axis=d.get("membrane_axis", 0),
    membrane_cut=d.get("membrane_cut", 3.0),
    seed=d.get("seed"),
)

# `simulation` is your openmm.app.Simulation, built from the same topology
cfg.validate(simulation.system.getNumParticles())
simulation.reporters.append(KAMDReporter(cfg, log="kamd.log"))
'''


def openmm_files(an: Analysis, project: Project) -> dict[str, str]:
    return {
        "kamd_openmm.json": json.dumps(openmm_dict(an, project), indent=None, separators=(",", ":")),
        "kamd_openmm_example.py": OPENMM_EXAMPLE,
    }


# --------------------------------------------------------------------------------------
# all + verification
# --------------------------------------------------------------------------------------
def export_files(an: Analysis, project: Project, engines=("amber", "gromacs", "openmm")) -> dict[str, str]:
    files: dict[str, str] = {}
    if "amber" in engines:
        files["kamd_amber.in"] = amber_namelist(an, project)
    if "gromacs" in engines:
        files.update(gromacs_files(an, project))
    if "openmm" in engines:
        files.update(openmm_files(an, project))
    files["project.json"] = json.dumps(project.to_dict(), indent=1)
    return files


def kamd_config_check(an: Analysis, project: Project) -> dict:
    try:
        from kamd_openmm.core import KAMDConfig, KAMDInputError
    except ImportError:
        return {"ok": None, "message": "kamd_openmm is not installed; KAMDConfig check skipped"}
    d = openmm_dict(an, project)
    try:
        KAMDConfig(
            refer_groups=d["refer_groups"], solvent_atoms=d["solvent_atoms"], temperature=d["temperature"],
            pseudobd_cycles=d["pseudobd_cycles"], nzero_crossings=d["nzero_crossings"],
            shell_cut=d["shell_cut"], surface_cut=d["surface_cut"], refer_link=d["refer_link"],
            region_type=d["region_type"], membrane_atoms=d.get("membrane_atoms"),
            membrane_axis=d.get("membrane_axis", 0), membrane_cut=d.get("membrane_cut", 3.0),
        ).validate(an.s.n_atoms)
    except KAMDInputError as exc:
        return {"ok": False, "message": f"KAMDConfig.validate: {exc}"}
    return {"ok": True, "message": f"KAMDConfig.validate({an.s.n_atoms}) passed"}


def _mask_atoms(structure, mask: str) -> np.ndarray:
    from parmed.amber import AmberMask

    return np.fromiter(AmberMask(structure, mask).Selected(), dtype=np.int64)


def verify(an: Analysis, project: Project, files: dict[str, str]) -> dict:
    rep: dict = {"kamd_config": kamd_config_check(an, project), "amber": [], "gromacs": [], "openmm": []}
    groups = group_atoms(an, project)
    b = an.bath
    if "kamd_amber.in" in files:
        from kamd_openmm.amber_input import parse_namelist

        nl = parse_namelist(files["kamd_amber.in"], "kamd")
        masks = nl.get("refer_mask", {})
        masks = masks if isinstance(masks, dict) else {1: masks}
        struct = an.s.structure
        for i, g in enumerate(groups):
            m = masks.get(i + 1)
            if struct is None or m is None:
                rep["amber"].append({"item": f"refer_mask({i + 1})", "ok": m is not None, "message": "not checked"})
                continue
            sel = _mask_atoms(struct, m)
            ok = np.array_equal(np.sort(sel), np.sort(g))
            rep["amber"].append({"item": f"refer_mask({i + 1})", "ok": bool(ok), "n_atoms": int(sel.size),
                                 "message": "ParmEd AmberMask selects exactly the group atoms" if ok else
                                 f"mismatch: mask selects {sel.size} atoms, group has {g.size}"})
        for key, want in (("solvent_mask", b.solvent_atoms), ("membrane_mask", b.membrane_atoms)):
            m = nl.get(key)
            if m is None or struct is None or (key == "solvent_mask" and m == PLACEHOLDER_SOLVENT and not want.size):
                continue
            sel = _mask_atoms(struct, str(m))
            ok = np.array_equal(np.sort(sel), np.sort(want))
            contiguous = len(runs(np.sort(sel))) == 1
            rep["amber"].append({"item": key, "ok": bool(ok and contiguous), "n_atoms": int(sel.size),
                                 "message": ("matches the bath, one contiguous range" if ok and contiguous else
                                             "matches the bath but is NOT one contiguous range" if ok else
                                             "does not match the bath atoms")})
    if "index.ndx" in files:
        nd = read_ndx(files["index.ndx"])
        for n, g in zip(gmx_names(project), groups, strict=True):
            ok = np.array_equal(nd.get(n, np.zeros(0)), g)
            rep["gromacs"].append({"item": n, "ok": bool(ok), "n_atoms": int(g.size)})
        if b.solvent_atoms.size:
            rep["gromacs"].append({"item": "KAMD_SOLVENT",
                                   "ok": bool(np.array_equal(nd.get("KAMD_SOLVENT"), b.solvent_atoms)),
                                   "n_atoms": int(b.solvent_atoms.size)})
    if "kamd_openmm.json" in files:
        d = json.loads(files["kamd_openmm.json"])
        for i, g in enumerate(groups):
            rep["openmm"].append({"item": f"refer_groups[{i}]", "ok": d["refer_groups"][i] == g.tolist(),
                                  "n_atoms": int(g.size)})
    rep["ok"] = all(x["ok"] for k in ("amber", "gromacs", "openmm") for x in rep[k]) and rep["kamd_config"]["ok"] in (
        True, None)
    return rep


def write_files(files: dict[str, str], out_dir: str | Path) -> list[str]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, text in files.items():
        p = out / name
        p.write_text(text)
        paths.append(str(p))
    return paths
