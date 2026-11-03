"""Automatic region suggestion.  Every method returns residue -> region id (0 = unassigned)."""

from __future__ import annotations

import numpy as np

from ..analysis import Analysis
from . import contacts, gnm, postprocess, tm

METHODS = {"gnm": gnm, "contacts": contacts, "tm": tm}


def method_info() -> dict:
    return {
        "gnm": {"label": "GNM slow modes", "params": gnm.DEFAULTS},
        "contacts": {"label": "Contact map + Leiden", "params": contacts.DEFAULTS},
        "tm": {"label": "Transmembrane segments", "params": tm.DEFAULTS},
        "postprocess": postprocess.DEFAULTS,
    }


def auto_method(an: Analysis) -> str:
    if an.membrane is not None and tm.tm_segments(an):
        return "tm"
    return "gnm"


def suggest(an: Analysis, method: str = "auto", params: dict | None = None) -> dict:
    """Run one method plus post-processing.

    Returns ``{"method", "labels" (np.ndarray per residue), "region_type" {id: 1|2},
    "info", "notes"}``.
    """
    params = dict(params or {})
    if method == "auto":
        method = auto_method(an)
    if method not in METHODS:
        raise ValueError(f"unknown method '{method}' (choose from {', '.join(METHODS)})")
    if method == "tm" and "min_heavy" not in params:
        params["min_heavy"] = params.get("min_heavy_tm", tm.DEFAULTS["min_heavy_tm"])
    raw, info = METHODS[method].suggest(an, params)
    rtype_raw = info.pop("region_type", None) if isinstance(info.get("region_type"), dict) else None
    labels, notes, rtype = postprocess.postprocess(raw, an, params, rtype_raw or None)
    if rtype is None:
        rtype = {int(g): 1 for g in np.unique(labels[labels > 0])}
    notes = list(info.pop("notes", [])) + notes
    return {"method": method, "labels": labels, "region_type": rtype, "info": info, "notes": notes}
