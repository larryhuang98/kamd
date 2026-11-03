"""Project state: regions, residue -> region assignment, KAMD settings.

The assignment is indexed by the topology's sequential residue index (0-based).
Region order in :attr:`Project.regions` is the KAMD group order (group 1, 2, ...).
``refer_link`` is stored as a region *id* and converted to a 1-based group position
on export.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

PROJECT_VERSION = 1

# Okabe-Ito without black (black disappears on the dark theme and under text),
# then tab20 for groups 8+.
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
TAB20 = [
    "#1f77b4",
    "#aec7e8",
    "#ff7f0e",
    "#ffbb78",
    "#2ca02c",
    "#98df8a",
    "#d62728",
    "#ff9896",
    "#9467bd",
    "#c5b0d5",
    "#8c564b",
    "#c49c94",
    "#e377c2",
    "#f7b6d2",
    "#7f7f7f",
    "#c7c7c7",
    "#bcbd22",
    "#dbdb8d",
    "#17becf",
    "#9edae5",
]
PALETTE = OKABE_ITO + TAB20
UNASSIGNED_COLOR = "#C8C8C8"


def default_color(i: int) -> str:
    return PALETTE[i % len(PALETTE)]


@dataclass
class Region:
    id: int
    name: str
    color: str
    region_type: int = 1  # 1 solution, 2 transmembrane
    refer_link: int = 0  # region id whose direction is shared, 0 = own


@dataclass
class Settings:
    pseudobd_cycles: int = 250
    nzero_crossings: int = 2
    shell_cut: float = 5.0  # Angstrom
    surface_cut: float = 3.0
    membrane_cut: float = 3.0
    kamd_shell_mode: int | None = None  # None: 2 for membrane systems, else not written
    temperature: float = 300.0
    seed: int | None = None
    membrane: dict | None = None  # manual slab {"axis": 3, "lo": -20, "hi": 20}; None = auto
    small_atoms: int = 200  # warn below (all atoms)
    large_atoms: int = 15000  # warn above (all atoms)


@dataclass
class Project:
    n_residues: int
    regions: list[Region] = field(default_factory=list)
    assignment: np.ndarray | None = None
    settings: Settings = field(default_factory=Settings)
    inputs: dict = field(default_factory=dict)
    method: dict = field(default_factory=dict)
    residue_names: list[str] = field(default_factory=list)  # consistency check on reload
    # chain / number / name of the assignable residues of the system this project was made on;
    # used to migrate the regions to another system by sequence alignment (migrate.py)
    residue_table: dict = field(default_factory=dict)
    edited: bool = False

    def __post_init__(self):
        if self.assignment is None:
            self.assignment = np.zeros(self.n_residues, dtype=np.int64)
        self.assignment = np.asarray(self.assignment, dtype=np.int64)

    # ------------------------------------------------------------------ helpers
    def region(self, rid: int) -> Region | None:
        for r in self.regions:
            if r.id == rid:
                return r
        return None

    def residues_of(self, rid: int) -> np.ndarray:
        return np.flatnonzero(self.assignment == rid)

    def next_id(self) -> int:
        return max([r.id for r in self.regions], default=0) + 1

    def group_position(self) -> dict[int, int]:
        """region id -> 1-based KAMD group number."""
        return {r.id: i + 1 for i, r in enumerate(self.regions)}

    @classmethod
    def from_labels(
        cls, labels: np.ndarray, region_type: dict | None = None, names: list[str] | None = None, **kw
    ) -> Project:
        labels = np.asarray(labels, dtype=np.int64)
        ids = [int(x) for x in np.unique(labels[labels > 0])]
        # order regions by their first residue
        ids.sort(key=lambda g: int(np.argmax(labels == g)))
        regions = [
            Region(g, f"region {k + 1}", default_color(k), int((region_type or {}).get(g, 1)), 0)
            for k, g in enumerate(ids)
        ]
        p = cls(n_residues=len(labels), regions=regions, assignment=labels.copy(), **kw)
        if names is not None:
            p.residue_names = list(names)
        return p

    # ------------------------------------------------------------------ JSON
    def to_dict(self) -> dict:
        return {
            "kamd_regions_project": PROJECT_VERSION,
            "inputs": self.inputs,
            "n_residues": self.n_residues,
            "residue_names": self.residue_names,
            "residue_table": self.residue_table,
            "regions": [asdict(r) for r in self.regions],
            "assignment": rle_encode(self.assignment),
            "settings": asdict(self.settings),
            "method": self.method,
            "edited": self.edited,
        }

    @classmethod
    def from_dict(cls, d: dict) -> Project:
        if "kamd_regions_project" not in d:
            raise ValueError("not a kamd-regions project file")
        n = int(d["n_residues"])
        a = d.get("assignment")
        assignment = rle_decode(a, n) if isinstance(a, list) and (not a or isinstance(a[0], list)) else np.asarray(a)
        st = Settings(**{k: v for k, v in d.get("settings", {}).items() if k in Settings.__dataclass_fields__})
        regions = [
            Region(**{k: v for k, v in r.items() if k in Region.__dataclass_fields__}) for r in d.get("regions", [])
        ]
        return cls(
            n_residues=n,
            regions=regions,
            assignment=assignment,
            settings=st,
            inputs=d.get("inputs", {}),
            method=d.get("method", {}),
            residue_names=d.get("residue_names", []),
            residue_table=d.get("residue_table") or {},
            edited=bool(d.get("edited", False)),
        )

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=1))

    @classmethod
    def load(cls, path: str | Path) -> Project:
        return cls.from_dict(json.loads(Path(path).read_text()))

    def fits(self, system) -> bool:
        """True if the residue assignment can be used on `system` as it is."""
        if self.n_residues != system.n_residues:
            return False
        names = [r.name for r in system.residues]
        return not self.residue_names or self.residue_names == names

    def check_against(self, system) -> list[str]:
        """Warnings if this project was made for a different topology."""
        out = []
        if self.n_residues != system.n_residues:
            out.append(
                f"project has {self.n_residues} residues but the loaded system has {system.n_residues}; "
                "the residue assignment does not fit this system"
            )
        elif self.residue_names:
            names = [r.name for r in system.residues]
            diff = sum(a != b for a, b in zip(self.residue_names, names, strict=False))
            if diff:
                out.append(f"{diff} residue names differ between the project and the loaded system")
        return out


def rle_encode(a: np.ndarray) -> list[list[int]]:
    """[[first, last, value], ...] (0-based, inclusive), zeros omitted."""
    a = np.asarray(a, dtype=np.int64)
    out = []
    i = 0
    n = a.size
    while i < n:
        j = i
        while j + 1 < n and a[j + 1] == a[i]:
            j += 1
        if a[i] != 0:
            out.append([i, j, int(a[i])])
        i = j + 1
    return out


def rle_decode(runs: list, n: int) -> np.ndarray:
    a = np.zeros(n, dtype=np.int64)
    for first, last, v in runs:
        a[int(first) : int(last) + 1] = int(v)
    return a
