"""One loaded system + its analysis + the current project (used by the server and the CLI)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import __version__
from .analysis import Analysis
from .export import export_files, verify, write_files
from .importers import import_amber, import_file, import_gromacs, import_openmm
from .membrane import Membrane
from .methods import method_info
from .methods import suggest as run_suggest
from .methods.gnm import hinges
from .methods.tm import tm_segments
from .project import Project, Settings
from .system import classify_files, load_system
from .validate import validate


class Session:
    def __init__(self, files: list[str], method: str = "auto", params: dict | None = None):
        roles = classify_files(files)
        sys_roles = {k: v for k, v in roles.items() if k in ("structure", "topology", "coordinates")}
        self.roles = roles
        self.system = load_system(sys_roles)
        self.an = Analysis(self.system)
        self.notes: list[str] = list(self.system.notes)
        if "structure2" in roles:
            self.notes.append("A second structure was given; DynDom-style analysis is on the roadmap and is not used.")
        if "pae" in roles:
            self.notes.append("A PAE file was given; the PAE method is on the roadmap and is not used yet.")
        self.project: Project
        if "project" in roles:
            self.project = Project.load(roles["project"])
            self.notes += self.project.check_against(self.system)
        elif "mdin" in roles:
            self.project, n = import_file(self.an, roles["mdin"])
            self.notes += n
        elif "kamd_dat" in roles:
            self.project, n = import_file(self.an, roles["kamd_dat"], roles.get("ndx"))
            self.notes += n
        elif "openmm_json" in roles:
            self.project, n = import_file(self.an, roles["openmm_json"])
            self.notes += n
        else:
            self.project, _ = self.suggest(method, params or {})
        self.project.inputs = {k: v for k, v in self.system.files.items()}

    # ---------------------------------------------------------------- suggestion
    def suggest(self, method: str = "auto", params: dict | None = None, settings: Settings | None = None):
        params = dict(params or {})
        st = settings or (self.project.settings if hasattr(self, "project") else Settings())
        saved = self.an.membrane
        if st.membrane:
            self.an.membrane = Membrane(
                int(st.membrane["axis"]), float(st.membrane["lo"]), float(st.membrane["hi"]), "manual"
            )
        try:
            res = run_suggest(self.an, method, params)
        finally:
            self.an.membrane = saved
        p = Project.from_labels(
            res["labels"], res["region_type"], names=[r.name for r in self.system.residues], settings=st
        )
        if self.an.membrane is not None and any(r.region_type == 2 for r in p.regions):
            count = {1: 0, 2: 0}
            for r in p.regions:  # TM1, TM2, ... and sol 1, sol 2, ... in sequence order
                count[r.region_type] += 1
                r.name = f"TM{count[2]}" if r.region_type == 2 else f"sol {count[1]}"
        p.method = {"method": res["method"], "params": params, "info": _jsonable(res["info"])}
        p.inputs = dict(self.system.files)
        return p, res["notes"]

    # ---------------------------------------------------------------- summaries
    def system_summary(self) -> dict:
        s = self.system
        an = self.an
        idx = np.flatnonzero(an.assignable)
        inside = an.in_slab() if an.membrane is not None else np.zeros(s.n_residues, dtype=bool)
        tm = np.zeros(s.n_residues, dtype=bool)
        for seg in tm_segments(an):
            tm[seg["start"] : seg["stop"]] = True
        msf = an.msf
        kinds: dict[str, int] = {}
        for k in s.res_kind:
            kinds[str(k)] = kinds.get(str(k), 0) + 1
        b = an.bath
        return {
            "version": __version__,
            "files": s.files,
            "n_atoms": s.n_atoms,
            "n_residues": s.n_residues,
            "has_topology": s.has_topology,
            "topology_format": s.topology_format,
            "kinds": kinds,
            "chains": s.chains(),
            "residues": {
                "index": idx.tolist(),
                "name": [s.residues[i].name for i in idx],
                "one": [s.residues[i].one for i in idx],
                "chain": [s.residues[i].chain for i in idx],
                "number": [s.residues[i].number for i in idx],
                "icode": [s.residues[i].icode for i in idx],
                "kind": [s.residues[i].kind for i in idx],
                "ss": [str(an.ss[i]) for i in idx],
                "n_atoms": s.res_natoms[idx].tolist(),
                "n_heavy": s.res_nheavy[idx].tolist(),
                "in_slab": inside[idx].tolist(),
                "tm": tm[idx].tolist(),
                "msf": np.round(msf[idx], 3).tolist(),
            },
            "membrane": an.membrane.to_dict() if an.membrane is not None else None,
            "bath": {
                "solvent_atoms": int(b.solvent_atoms.size),
                "membrane_atoms": int(b.membrane_atoms.size),
                "solvent_contiguous": not any(p["bath"] == "solvent" for p in b.contiguity(s)),
                "membrane_contiguous": not any(p["bath"] == "membrane" for p in b.contiguity(s)),
                "n_lipids": an.membrane.n_lipids if an.membrane is not None else 0,
            },
            "hinges": hinges(an) if an.polymer.size > 4 else [],
            "methods": method_info(),
            "notes": self.notes,
        }

    def viewer_structure(self) -> dict:
        """PDB text for the 3D view: solute heavy atoms (+ lipids, ions) and atom -> residue maps."""
        s = self.system
        solute = s.atoms_of(np.flatnonzero(self.an.assignable))
        solute = solute[s.heavy[solute]]
        lip = s.atoms_of(self.an.bath.membrane_residues)
        lip = lip[s.heavy[lip]]
        ions = s.atoms_of(np.flatnonzero(s.res_kind == "ion"))
        return {
            "solute_pdb": _pdb(s, solute),
            "solute_res": s.atom_res[solute].tolist(),
            "lipid_pdb": _pdb(s, lip) if lip.size else None,
            "ion_pdb": _pdb(s, ions) if ions.size else None,
        }

    # ---------------------------------------------------------------- actions
    def project_from(self, d: dict | None) -> Project:
        return Project.from_dict(d) if d else self.project

    def validate(self, d: dict | None = None) -> dict:
        p = self.project_from(d)
        return validate(self.an, p)

    def export(
        self, d: dict | None, engines=("amber", "gromacs", "openmm"), out_dir: str | None = None, force: bool = False
    ) -> dict:
        p = self.project_from(d)
        p.inputs = dict(self.system.files)
        v = validate(self.an, p, engines)
        blocked = [e for e in engines if not v["ok"].get(e, True)]
        if blocked and not force:
            return {"ok": False, "blocked": blocked, "validation": v, "files": {}, "verification": None}
        files = export_files(self.an, p, engines)
        rep = verify(self.an, p, files)
        written = write_files(files, out_dir) if out_dir else []
        return {"ok": True, "blocked": [], "validation": v, "files": files, "verification": rep, "written": written}

    def import_text(self, kind: str, text: str, ndx_text: str | None = None):
        if kind == "amber":
            p, n = import_amber(self.an, text)
        elif kind == "gromacs":
            p, n = import_gromacs(self.an, text, ndx_text or "")
        elif kind == "openmm":
            p, n = import_openmm(self.an, text)
        else:
            raise ValueError(kind)
        p.inputs = dict(self.system.files)
        return p, n

    def save(self, d: dict | None, path: str) -> str:
        p = self.project_from(d)
        p.inputs = dict(self.system.files)
        p.save(path)
        return str(Path(path).resolve())


def _pdb(s, atoms: np.ndarray) -> str:
    out = []
    for k, a in enumerate(atoms):
        r = s.residues[int(s.atom_res[a])]
        name = str(s.atom_name[a])[:4]
        el = str(s.element[a])[:2]
        nm = f" {name:<3s}" if len(name) < 4 and len(el) == 1 else f"{name:<4s}"
        x, y, z = s.xyz[a]
        rec = "ATOM  " if r.kind in ("protein", "nucleic") else "HETATM"
        out.append(
            f"{rec}{(k + 1) % 100000:5d} {nm} {r.name[:4]:<4s}{(r.chain or 'A')[:1]}{r.number % 10000:4d}"
            f"{(r.icode or ' ')[:1]}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}"
        )
    out.append("END")
    return "\n".join(out) + "\n"


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x
