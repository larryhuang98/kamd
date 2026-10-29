"""Build a KAMDConfig from an Amber mdin (&kamd namelist) + prmtop.

This lets the same input drive pmemd.cuda and OpenMM, which is how the port is
validated.  Masks are resolved with ParmEd's AmberMask (same syntax as pmemd).
"""

from __future__ import annotations

import re

import numpy as np

from .core import KAMDConfig

_KEY = re.compile(r"([A-Za-z_]\w*)\s*(?:\(\s*(\d+)\s*\))?\s*=")


def _strip_comments(text: str) -> str:
    out, q = [], None
    for line in text.splitlines():
        buf = []
        for ch in line:
            if q:
                buf.append(ch)
                if ch == q:
                    q = None
            elif ch in "'\"":
                q = ch
                buf.append(ch)
            elif ch == "!":
                break
            else:
                buf.append(ch)
        out.append("".join(buf))
    return "\n".join(out)


def _block(text: str, name: str) -> str:
    text = _strip_comments(text)
    m = re.search(r"&" + name + r"\b", text, re.IGNORECASE)
    if not m:
        raise KeyError(f"&{name} namelist not found")
    i, q = m.end(), None
    while i < len(text):
        ch = text[i]
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "/" or (ch == "&" and text[i : i + 4].lower() == "&end"):
            return text[m.end() : i]
        i += 1
    raise ValueError(f"&{name} namelist is not terminated")


def _items(s: str):
    """Values of one assignment: quoted strings stay strings, bare tokens are converted."""
    out = []
    for m in re.finditer(r'"([^"]*)"|\'([^\']*)\'|([^,\s]+)', s):
        if m.group(1) is not None:
            out.append(m.group(1))
        elif m.group(2) is not None:
            out.append(m.group(2))
        else:
            out.append(_value(m.group(3)))
    return out


def _value(tok: str):
    t = tok.strip()
    try:
        return int(t)
    except ValueError:
        pass
    try:
        return float(t.lower().replace("d", "e"))
    except ValueError:
        return t


def parse_namelist(text: str, name: str) -> dict[str, object]:
    """Parse one Fortran namelist into {key: scalar | {index: value}} (keys lower-case)."""
    body = _block(text, name)
    # look for "key =" only outside quoted strings (masks contain '=' e.g. "@H=")
    masked = re.sub(r'"[^"]*"|\'[^\']*\'', lambda m: "#" * len(m.group(0)), body)
    keys = list(_KEY.finditer(masked))
    res: dict[str, object] = {}
    for k, m in enumerate(keys):
        end = keys[k + 1].start() if k + 1 < len(keys) else len(body)
        raw = body[m.end() : end]
        vals = _items(raw)
        key = m.group(1).lower()
        if m.group(2) is not None or len(vals) > 1:
            start = int(m.group(2)) if m.group(2) else 1
            d = res.setdefault(key, {})
            if not isinstance(d, dict):
                d = res[key] = {}
            for j, v in enumerate(vals):
                d[start + j] = v
        elif vals:
            res[key] = vals[0]
    return res


def amber_mask_indices(structure, mask: str) -> np.ndarray:
    """0-based atom indices selected by an Amber mask (ParmEd AmberMask)."""
    from parmed.amber import AmberMask

    return np.fromiter(AmberMask(structure, mask).Selected(), dtype=np.int64)


def _packed(d, n_max=32):
    if not isinstance(d, dict):
        return [d] if d not in (None, "") else []
    return [d[i] for i in range(1, n_max + 1) if i in d]


def config_from_amber(mdin: str, prmtop, temperature=None, seed=None) -> KAMDConfig:
    """KAMDConfig from an Amber mdin file (or text) and a prmtop path / ParmEd structure.

    The temperature defaults to ``temp0`` of &cntrl.  Only the pseudo-BD part of
    &kamd is used; ion_traj_* and kamd_shell_* (experimental membrane shell) are
    ignored with a warning.
    """
    import warnings

    try:
        text = open(mdin).read()
    except (OSError, TypeError, ValueError):
        text = mdin
    kamd = parse_namelist(text, "kamd")
    if temperature is None:
        cntrl = parse_namelist(text, "cntrl")
        temperature = float(cntrl.get("temp0", 300.0))

    if isinstance(prmtop, str):
        import parmed

        structure = parmed.load_file(prmtop)
    else:
        structure = prmtop

    masks = _packed(kamd.get("refer_mask", {}))
    if not masks:
        raise ValueError("no refer_mask(i) in &kamd")
    raw = kamd.get("refer_mask", {})
    if isinstance(raw, dict) and sorted(raw) != list(range(1, len(raw) + 1)):
        raise ValueError("refer_mask entries must be packed (no gaps), as in pmemd")
    ng = len(masks)
    groups = [amber_mask_indices(structure, m) for m in masks]
    for i, (m, g) in enumerate(zip(masks, groups, strict=False)):
        if g.size == 0:
            raise ValueError(f'refer_mask({i + 1}) = "{m}" matched no atoms')

    def per_group(key, default):
        d = kamd.get(key, {})
        d = d if isinstance(d, dict) else {1: d}
        return [int(d.get(i, default)) for i in range(1, ng + 1)]

    if "solvent_mask" not in kamd:
        raise ValueError("solvent_mask is required in &kamd")
    solvent = amber_mask_indices(structure, str(kamd["solvent_mask"]))
    membrane = amber_mask_indices(structure, str(kamd["membrane_mask"])) if kamd.get("membrane_mask") else None

    for k in kamd:
        if k.startswith(("ion_traj", "ion_mask", "kamd_shell")):
            warnings.warn(f"&kamd {k} is not supported by kamd_openmm and is ignored", stacklevel=2)

    return KAMDConfig(
        refer_groups=groups,
        solvent_atoms=solvent,
        temperature=float(temperature),
        pseudobd_cycles=int(kamd.get("pseudobd_cycles", 0)),
        nzero_crossings=int(kamd.get("nzero_crossings", 2)),
        shell_cut=float(kamd.get("shell_cut", 5.0)),
        surface_cut=float(kamd.get("surface_cut", 3.0)),
        refer_link=per_group("refer_link", 0),
        region_type=per_group("region_type", 1),
        membrane_atoms=membrane,
        membrane_axis=int(kamd.get("membrane_axis", 0)),
        membrane_cut=float(kamd.get("membrane_cut", 3.0)),
        seed=seed,
    )
