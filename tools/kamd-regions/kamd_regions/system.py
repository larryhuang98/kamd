"""Load a structure (and optionally a topology) into one engine-neutral model.

Primary key everywhere in this package is the **topology's sequential residue index**
(0-based in Python, 1-based in AmberMask ``:N``).  PDB chain IDs and residue numbers are
kept only for display.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .residues import ASSIGNABLE_KINDS, classify, is_lipid21_fragment, one_letter

STRUCTURE_EXT = {".pdb", ".ent", ".cif", ".mmcif", ".pdbx", ".gro", ".pqr"}
COORD_EXT = {".inpcrd", ".rst7", ".rst", ".crd", ".restrt", ".ncrst", ".nc", ".netcdf"}
AMBER_TOP_EXT = {".prmtop", ".parm7", ".top", ".prm"}
GMX_TOP_EXT = {".top", ".tpr"}


class LoadError(ValueError):
    pass


@dataclass
class Residue:
    index: int  # sequential topology residue index, 0-based
    name: str
    chain: str
    number: int  # residue number shown to the user (PDB resid / GROMACS resnr)
    icode: str
    start: int  # first atom (0-based)
    stop: int  # one past the last atom
    kind: str
    molecule: int = -1
    one: str = "X"

    @property
    def n_atoms(self) -> int:
        return self.stop - self.start

    def label(self) -> str:
        c = self.chain or "-"
        return f"{c}:{self.name} {self.number}{self.icode.strip()}"


@dataclass
class System:
    atom_name: np.ndarray
    element: np.ndarray
    mass: np.ndarray
    xyz: np.ndarray  # (n, 3) Angstrom
    residues: list[Residue]
    box: np.ndarray | None = None  # (a, b, c) Angstrom, orthorhombic part only
    has_topology: bool = False
    topology_format: str = ""  # "amber", "gromacs", ""
    files: dict = field(default_factory=dict)
    structure: object = None  # parmed.Structure when available (AmberMask parsing)
    notes: list = field(default_factory=list)

    def __post_init__(self):
        self.n_atoms = int(self.xyz.shape[0])
        self.atom_res = np.empty(self.n_atoms, dtype=np.int64)
        for r in self.residues:
            self.atom_res[r.start : r.stop] = r.index
        self.heavy = np.array([e not in ("H", "D", "EP", "LP", "") for e in self.element])
        self.res_start = np.array([r.start for r in self.residues], dtype=np.int64)
        self.res_stop = np.array([r.stop for r in self.residues], dtype=np.int64)
        self.res_kind = np.array([r.kind for r in self.residues])
        self.res_natoms = self.res_stop - self.res_start
        self.res_nheavy = np.add.reduceat(self.heavy.astype(np.int64), self.res_start) if self.residues else np.zeros(0)

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    def assignable(self) -> np.ndarray:
        return np.isin(self.res_kind, ASSIGNABLE_KINDS)

    def atoms_of(self, residue_indices) -> np.ndarray:
        idx = np.asarray(residue_indices, dtype=np.int64)
        if idx.size == 0:
            return np.zeros(0, dtype=np.int64)
        return np.concatenate([np.arange(self.res_start[i], self.res_stop[i]) for i in np.sort(idx)])

    def atoms_of_kind(self, kinds) -> np.ndarray:
        return self.atoms_of(np.flatnonzero(np.isin(self.res_kind, list(kinds))))

    def chains(self) -> list[dict]:
        """Polymer chains (consecutive assignable residues sharing a chain label)."""
        out: list[dict] = []
        for r in self.residues:
            if r.kind not in ASSIGNABLE_KINDS:
                continue
            if out and out[-1]["chain"] == r.chain and out[-1]["stop"] == r.index:
                out[-1]["stop"] = r.index + 1
            else:
                out.append({"chain": r.chain, "start": r.index, "stop": r.index + 1})
        return out


# --------------------------------------------------------------------------------------
# file roles
# --------------------------------------------------------------------------------------
def _head(path: str, n: int = 4096) -> str:
    try:
        with open(path, "rb") as fh:
            return fh.read(n).decode("latin-1", "replace")
    except OSError:
        return ""


def classify_files(paths: list[str]) -> dict[str, str]:
    """Assign a role to every input file by extension and content."""
    roles: dict[str, str] = {}

    def put(role, p):
        if role in roles:
            raise LoadError(f"two files for role '{role}': {roles[role]} and {p}")
        roles[role] = p

    for p in paths:
        if not os.path.exists(p):
            raise LoadError(f"file not found: {p}")
        ext = Path(p).suffix.lower()
        head = _head(p)
        name = Path(p).name.lower()
        if "%FLAG" in head and "%VERSION" in head or ext in (".prmtop", ".parm7"):
            put("topology", p)
        elif ext == ".tpr":
            put("topology", p)
        elif ext in (".top", ".itp"):
            put("topology", p)
        elif ext == ".ndx":
            put("ndx", p)
        elif ext == ".json":
            text = _head(p, 2_000_000)
            if '"kamd_regions_project"' in text:
                put("project", p)
            elif "predicted_aligned_error" in text or '"pae"' in text:
                put("pae", p)
            elif '"refer_groups"' in text:
                put("openmm_json", p)
            else:
                raise LoadError(f"unrecognised JSON file: {p}")
        elif "&kamd" in head.lower() or ext in (".in", ".mdin"):
            put("mdin", p)
        elif ext == ".dat" or "refer-groups" in head or "refer_groups" in head:
            put("kamd_dat", p)
        elif ext in STRUCTURE_EXT:
            if "structure" in roles:
                put("structure2", p)
            else:
                put("structure", p)
        elif ext in COORD_EXT or name in ("restrt", "inpcrd", "rst7"):
            put("coordinates", p)
        else:
            # Amber ASCII restart without extension: title line, then natom (+time)
            lines = head.splitlines()
            if len(lines) > 2 and lines[1].split() and lines[1].split()[0].isdigit():
                put("coordinates", p)
            else:
                raise LoadError(f"cannot tell what kind of file this is: {p}")
    if "structure" not in roles and "coordinates" not in roles and "topology" not in roles:
        raise LoadError("no structure file given (PDB/mmCIF/gro, or a topology with coordinates)")
    return roles


# --------------------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------------------
def load_system(paths: list[str] | dict[str, str]) -> System:
    roles = classify_files(list(paths)) if not isinstance(paths, dict) else dict(paths)
    top = roles.get("topology")
    crd = roles.get("coordinates") or roles.get("structure")
    if top and Path(top).suffix.lower() == ".tpr" or (top and _is_gromacs_top(top)):
        system = _load_mdanalysis(top, crd)
    elif top:
        system = _load_parmed_amber(top, crd)
    else:
        if crd is None:
            raise LoadError("no structure or coordinate file")
        system = _load_parmed_structure(crd)
    system.files = {k: str(Path(v).resolve()) for k, v in roles.items()}
    return system


def _is_gromacs_top(path: str) -> bool:
    h = _head(path)
    return Path(path).suffix.lower() in (".top", ".itp") and "%FLAG" not in h


def _elements_parmed(struct) -> list[str]:
    from parmed.periodic_table import Element

    out = []
    for a in struct.atoms:
        z = a.element
        out.append(Element[z] if z and z > 0 else _guess_element(a.name, a.mass))
    return out


def _guess_element(name: str, mass: float | None = None) -> str:
    if mass:
        for el, m in (("H", 1.008), ("C", 12.011), ("N", 14.007), ("O", 15.999), ("P", 30.974), ("S", 32.06)):
            if abs(mass - m) < 0.6:
                return el
    n = name.strip().lstrip("0123456789").upper()
    if n[:2] in ("CL", "NA", "MG", "ZN", "FE", "CA") and len(n) <= 3 and mass and mass > 20:
        return n[:2].capitalize()
    return n[:1] if n else ""


def _from_parmed(struct, has_topology: bool, fmt: str) -> System:
    n = len(struct.atoms)
    if struct.coordinates is None:
        raise LoadError("the structure has no coordinates; give a coordinate file (inpcrd/rst7/pdb)")
    xyz = np.asarray(struct.coordinates, dtype=np.float64).reshape(n, 3)
    names = np.array([a.name for a in struct.atoms], dtype=object)
    elements = np.array(_elements_parmed(struct), dtype=object)
    mass = np.array([a.mass or 0.0 for a in struct.atoms], dtype=np.float64)
    residues = []
    for i, r in enumerate(struct.residues):
        idx = [a.idx for a in r.atoms]
        start, stop = min(idx), max(idx) + 1
        if stop - start != len(idx):
            raise LoadError(f"residue {i + 1} ({r.name}) is not a contiguous atom range")
        chain = (r.chain or "").strip() or (getattr(r, "segid", "") or "").strip()
        kind = classify(r.name, [a.name for a in r.atoms], list(elements[start:stop]))
        residues.append(
            Residue(i, r.name, chain, int(r.number), (r.insertion_code or "").strip(), start, stop, kind)
        )
    # ParmEd numbers prmtop residues from 0 when there is no RESIDUE_NUMBER flag
    if residues and all(r.number == r.index for r in residues):
        for r in residues:
            r.number = r.index + 1
    box = None
    if struct.box is not None:
        box = np.asarray(struct.box[:3], dtype=np.float64)
    bonds = None
    if has_topology and struct.bonds:
        bonds = np.array([(b.atom1.idx, b.atom2.idx) for b in struct.bonds], dtype=np.int64)
    sysm = System(names, elements, mass, xyz, residues, box, has_topology, fmt, structure=struct)
    _finish(sysm, bonds)
    return sysm


def _load_parmed_amber(top: str, crd: str | None) -> System:
    import parmed

    if crd is None:
        raise LoadError("a topology needs a coordinate file (inpcrd/rst7/restrt/pdb)")
    try:
        struct = parmed.load_file(top, xyz=crd)
    except Exception:
        struct = parmed.load_file(top)
        c = parmed.load_file(crd)
        coords = c.coordinates if hasattr(c, "coordinates") else None
        if coords is None or np.asarray(coords).reshape(-1, 3).shape[0] != len(struct.atoms):
            raise LoadError(f"{crd}: atom count does not match {top}") from None
        struct.coordinates = np.asarray(coords).reshape(-1, 3)
        if getattr(c, "box", None) is not None:
            struct.box = c.box
    if struct.coordinates is None:
        raise LoadError(f"could not read coordinates from {crd}")
    fmt = "amber" if isinstance(struct, parmed.amber.AmberParm) else "other"
    return _from_parmed(struct, True, fmt)


def _load_parmed_structure(path: str) -> System:
    import parmed

    try:
        struct = parmed.load_file(path)
    except Exception as exc:  # pragma: no cover - depends on file
        raise LoadError(f"could not read {path}: {exc}") from exc
    if isinstance(struct, list):
        struct = struct[0]
    s = _from_parmed(struct, False, "")
    s.notes.append(
        "No topology given: residues are numbered in file order. AmberMask residue numbers "
        "may not match a prmtop built from this structure."
    )
    return s


def _load_mdanalysis(top: str, crd: str | None) -> System:
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise LoadError("GROMACS topologies need MDAnalysis: pip install 'kamd-regions[gromacs]'") from exc
    u = mda.Universe(top, crd) if crd else mda.Universe(top)
    n = len(u.atoms)
    names = np.array(list(u.atoms.names), dtype=object)
    try:
        elements = np.array([e.capitalize() if e else "" for e in u.atoms.elements], dtype=object)
    except Exception:
        elements = np.array([_guess_element(a.name, a.mass) for a in u.atoms], dtype=object)
    if any(e == "" for e in elements):
        elements = np.array([_guess_element(a.name, a.mass) for a in u.atoms], dtype=object)
    mass = np.asarray(u.atoms.masses, dtype=np.float64)
    xyz = np.asarray(u.atoms.positions, dtype=np.float64).reshape(n, 3)
    residues = []
    for i, r in enumerate(u.residues):
        idx = r.atoms.indices
        start, stop = int(idx.min()), int(idx.max()) + 1
        if stop - start != len(idx):
            raise LoadError(f"residue {i + 1} ({r.resname}) is not a contiguous atom range")
        chain = ""
        try:
            chain = str(r.atoms[0].chainID).strip()
        except Exception:
            pass
        if not chain:
            try:
                seg = str(r.segment.segid)
                chain = "" if seg.startswith("seg_") or seg == "SYSTEM" else seg
            except Exception:
                pass
        kind = classify(r.resname, list(r.atoms.names), list(elements[start:stop]))
        residues.append(Residue(i, r.resname, chain, int(r.resid), "", start, stop, kind))
    box = np.asarray(u.dimensions[:3], dtype=np.float64) if u.dimensions is not None else None
    bonds = None
    try:
        if len(u.bonds):
            bonds = np.asarray(u.bonds.indices, dtype=np.int64)
    except Exception:
        pass
    struct = _parmed_from_arrays(names, elements, mass, xyz, residues)
    s = System(names, elements, mass, xyz, residues, box, True, "gromacs", structure=struct)
    _finish(s, bonds)
    return s


def _parmed_from_arrays(names, elements, mass, xyz, residues):
    """Minimal ParmEd structure (names + residues) so AmberMask can be evaluated."""
    import parmed
    from parmed.periodic_table import AtomicNum

    st = parmed.Structure()
    for r in residues:
        for a in range(r.start, r.stop):
            z = AtomicNum.get(str(elements[a]).capitalize(), 0)
            st.add_atom(parmed.Atom(name=str(names[a]), type=str(names[a]), atomic_number=z, mass=mass[a]),
                        r.name, r.number, r.chain)
    st.coordinates = xyz
    return st


# --------------------------------------------------------------------------------------
# molecules and chains
# --------------------------------------------------------------------------------------
def _finish(s: System, bonds: np.ndarray | None) -> None:
    s.mol = molecule_ids(s, bonds)
    for r in s.residues:
        r.molecule = int(s.mol[r.index])
        r.one = one_letter(r.name) if r.kind in ("protein", "nucleic") else "X"
    # prmtop has no chain IDs: label every polymer molecule A, B, C, ...
    if not any(r.chain for r in s.residues if r.kind in ASSIGNABLE_KINDS):
        letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
        mol_chain: dict[int, str] = {}
        for r in s.residues:
            if r.kind in ("protein", "nucleic"):
                if r.molecule not in mol_chain:
                    mol_chain[r.molecule] = letters[len(mol_chain) % len(letters)]
                r.chain = mol_chain[r.molecule]
        s.notes.append("The topology has no chain IDs; polymer molecules were labelled A, B, ... in order.")


def molecule_ids(s: System, bonds: np.ndarray | None) -> np.ndarray:
    """Molecule index per residue.

    With bonds (topology) molecules are the connected components of the bond graph.
    Without bonds, consecutive polymer residues of one chain form one molecule and
    consecutive Amber Lipid21 fragments that are covalently bonded (closest atoms
    < 2 A) are merged into one lipid.
    """
    nres = s.n_residues
    if bonds is not None and len(bonds):
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components

        rb = s.atom_res[bonds]
        rb = rb[rb[:, 0] != rb[:, 1]]
        g = coo_matrix((np.ones(len(rb)), (rb[:, 0], rb[:, 1])), shape=(nres, nres))
        _, lab = connected_components(g, directed=False)
        # renumber in order of first residue
        order = {}
        out = np.empty(nres, dtype=np.int64)
        for i, x in enumerate(lab):
            out[i] = order.setdefault(int(x), len(order))
        return out
    out = np.empty(nres, dtype=np.int64)
    m = -1
    prev = None
    for r in s.residues:
        join = False
        if prev is not None:
            if r.kind in ("protein", "nucleic") and prev.kind == r.kind and prev.chain == r.chain:
                join = _bonded(s, prev, r, 2.2) or _ca_close(s, prev, r)
            elif r.kind == "lipid" and prev.kind == "lipid" and (
                is_lipid21_fragment(r.name) or is_lipid21_fragment(prev.name)
            ):
                join = _bonded(s, prev, r, 2.0)
        if not join:
            m += 1
        out[r.index] = m
        prev = r
    return out


def _mic(d: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    if box is None or np.any(box <= 0):
        return d
    return d - box * np.round(d / box)


def _bonded(s: System, a: Residue, b: Residue, cut: float) -> bool:
    xa = s.xyz[a.start : a.stop]
    xb = s.xyz[b.start : b.stop]
    d = _mic(xa[:, None, :] - xb[None, :, :], s.box)
    return bool(np.min(np.einsum("ijk,ijk->ij", d, d)) < cut * cut)


def _ca_close(s: System, a: Residue, b: Residue) -> bool:
    def ca(r):
        for i in range(r.start, r.stop):
            if s.atom_name[i] == "CA":
                return s.xyz[i]
        return None

    x, y = ca(a), ca(b)
    return x is not None and y is not None and np.linalg.norm(x - y) < 4.3
