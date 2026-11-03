"""Residue-name tables: protein, nucleic acid, water, ions and lipids.

Lipids cover CHARMM36 (one residue per lipid, e.g. ``POPC``) and Amber Lipid21 /
Lipid17 / Lipid14 (one lipid split into head and tail residues, e.g. ``PA PC OL``).
The Lipid21 fragments are merged back into whole molecules in
:func:`kamd_regions.system.molecule_ids`.
"""

from __future__ import annotations

import re

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # protonation / force-field variants (Amber, CHARMM, GROMACS)
    "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "HISD": "H", "HISE": "H",
    "HISH": "H", "HISA": "H", "HISB": "H", "CYX": "C", "CYM": "C", "CYS2": "C", "ASH": "D", "GLH": "E",
    "LYN": "K", "ARN": "R", "ASPP": "D", "GLUP": "E", "LSN": "K", "TYM": "Y", "MSE": "M", "SEP": "S",
    "TPO": "T", "PTR": "Y", "HYP": "P", "NLE": "L", "ACE": "X", "NME": "X", "NMA": "X", "NH2": "X",
}
# names mdtraj's DSSP understands
STANDARD_NAME = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HISD": "HIS",
    "HISE": "HIS", "HISH": "HIS", "HISA": "HIS", "HISB": "HIS", "CYX": "CYS", "CYM": "CYS", "CYS2": "CYS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "ARN": "ARG", "ASPP": "ASP", "GLUP": "GLU", "LSN": "LYS",
    "TYM": "TYR", "MSE": "MET", "SEP": "SER", "TPO": "THR", "PTR": "TYR", "HYP": "PRO", "NLE": "LEU",
}

NUCLEIC = {
    "DA", "DC", "DG", "DT", "DU", "A", "C", "G", "U", "T", "RA", "RC", "RG", "RU",
    "DA5", "DA3", "DC5", "DC3", "DG5", "DG3", "DT5", "DT3", "RA5", "RA3", "RC5", "RC3", "RG5", "RG3",
    "RU5", "RU3", "A5", "A3", "C5", "C3", "G5", "G3", "U5", "U3", "ADE", "CYT", "GUA", "THY", "URA",
}

WATER = {
    "WAT", "HOH", "H2O", "SOL", "TIP3", "TIP3P", "TIP4", "TIP4P", "TIP5", "TP3", "TP4", "T3P", "T4P",
    "T4E", "SPC", "SPCE", "OPC", "OPC3", "SWM4", "DOD", "TIP", "WATR", "W",
}

IONS = {
    "NA", "NA+", "Na+", "SOD", "K", "K+", "POT", "CL", "CL-", "Cl-", "CLA", "MG", "MG2", "Mg+", "Mg2+",
    "CA", "CAL", "Ca2+", "ZN", "ZN2", "Zn2+", "LI", "Li+", "LIT", "RB", "Rb+", "CS", "Cs+", "CES", "BR",
    "Br-", "I", "I-", "IOD", "F", "F-", "MN", "Mn2+", "FE", "FE2", "Fe2+", "CU", "Cu2+", "CD", "Cd2+",
    "BA", "Ba2+", "SR", "Sr2+", "NI", "Ni2+", "CO", "Co2+", "IB", "IB+", "NIO", "CIO", "KIO",
}

# CHARMM36 lipids (CHARMM-GUI names). A regular expression below catches most other
# phospholipids of the form <2 letters><headgroup>.
CHARMM_LIPIDS = {
    "POPC", "POPE", "POPS", "POPG", "POPA", "POPI", "DPPC", "DPPE", "DPPS", "DPPG", "DMPC", "DMPE", "DMPG",
    "DOPC", "DOPE", "DOPS", "DOPG", "DOPA", "DLPC", "DLPE", "DSPC", "DSPE", "DYPC", "DYPE", "SAPC", "SAPE",
    "SOPC", "SDPC", "PLPC", "PSPC", "PAPC", "PAPE", "PAPS", "SAPI", "SAPI14", "SAPI15", "SAPI24", "SAPI25",
    "POPI14", "POPI15", "POPI24", "POPI25", "POPI2A", "POPI2B", "POPI33", "PIP2", "PI2", "PI3",
    "CHL1", "CHOL", "CHL", "ERG", "SITO", "STIG", "PSM", "SSM", "NSM", "LSM", "ASM", "OSM", "CER160",
    "CER180", "TOCL", "TOCL1", "TOCL2", "TLCL", "TLCL2", "PMCL", "PVCL", "LPPC", "LPPE", "DHPC",
    "DLIPC", "DLIPE", "POPP", "PVPC", "PVPE", "PVPG", "OPPC", "YOPC", "YOPE", "YOPS", "DIPC", "DIPE",
    "DUPC", "DEPC", "DGPC", "DNPC", "DXPC", "QMPE", "PYPE", "PYPG", "DBPC", "PRPC", "PGPC",
}
_CHARMM_LIPID_RE = re.compile(r"^[A-Z]{2}(PC|PE|PS|PG|PA|PI|SM|CL)\d{0,2}$")

# Amber Lipid21 / Lipid17 / Lipid14 fragments
LIPID21_HEADS = {"PC", "PE", "PS", "PGR", "PGS", "PH-", "PI", "SPM", "CER", "PA-", "H2-", "CL"}
LIPID21_TAILS = {"LA", "MY", "PA", "ST", "OL", "LEO", "LEN", "AR", "DHA", "SA", "PO", "EI", "BE", "DO"}
LIPID21_WHOLE = {"CHL", "CHOL", "CHL1"}
LIPID21_FRAGMENTS = LIPID21_HEADS | LIPID21_TAILS


def is_lipid_name(name: str) -> bool:
    n = name.strip().upper()
    return (
        n in CHARMM_LIPIDS
        or n in LIPID21_FRAGMENTS
        or n in LIPID21_WHOLE
        or bool(_CHARMM_LIPID_RE.match(n))
    )


def is_lipid21_fragment(name: str) -> bool:
    return name.strip().upper() in LIPID21_FRAGMENTS


def one_letter(name: str) -> str:
    n = name.strip().upper()
    if n in THREE_TO_ONE:
        return THREE_TO_ONE[n]
    if n in NUCLEIC:
        base = n.rstrip("35")
        return base[-1] if base else "X"
    return "X"


def classify(name: str, atom_names: list[str], elements: list[str]) -> str:
    """Kind of a residue: protein, nucleic, water, ion, lipid or other."""
    n = name.strip()
    u = n.upper()
    if u == "DUM":  # OPM membrane boundary pseudo-atoms
        return "dummy"
    # 'CA' is both C-alpha and calcium: a one-atom residue is never an amino acid
    if u in THREE_TO_ONE and len(atom_names) > 1:
        return "protein"
    if {"N", "CA", "C"} <= set(atom_names) and len(atom_names) > 3 and u not in WATER:
        return "protein"
    if u in WATER or n in WATER:
        return "water"
    if len(atom_names) <= 4 and elements and set(elements) <= {"O", "H", "EP", ""} and "O" in elements:
        return "water"
    if (u in IONS or n in IONS) and len(atom_names) <= 2:
        return "ion"
    if u in NUCLEIC and ("P" in atom_names or "C1'" in atom_names or "C4'" in atom_names):
        return "nucleic"
    # Lipid21 tail "PA" (palmitoyl) clashes with nothing in the tables above; "CL"
    # (cardiolipin head) clashes with chloride, which is a single atom.
    if is_lipid_name(u) and len(atom_names) > 2:
        return "lipid"
    if len(atom_names) == 1 and elements and elements[0] not in ("C", "N", "O", "S", "P", "H"):
        return "ion"
    return "other"


ASSIGNABLE_KINDS = ("protein", "nucleic", "other")
BATH_KINDS = ("water", "ion", "lipid")
