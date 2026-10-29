"""KAMD (Kinetically Accelerated Molecular Dynamics) for OpenMM.

Port of the pmemd.cuda implementation by Z. Huang, H. Wei and R. Luo (UC Irvine).
Method: H. Wei et al., bioRxiv 2025, doi:10.1101/2025.06.14.659665.
"""

from .core import (
    KAMDConfig,
    KAMDInputError,
    MomentumMonitor,
    ShellResult,
    find_first_shell,
    group_directions,
    group_momenta,
    maxwell_boltzmann,
    pseudobd_reset,
    redirect,
)

__version__ = "0.1.0"

__all__ = [
    "KAMDConfig",
    "KAMDInputError",
    "MomentumMonitor",
    "ShellResult",
    "find_first_shell",
    "group_directions",
    "group_momenta",
    "maxwell_boltzmann",
    "pseudobd_reset",
    "redirect",
    "KAMDReporter",
    "config_from_amber",
    "parse_namelist",
    "__version__",
]


def __getattr__(name):  # lazy: core works without OpenMM / ParmEd installed
    if name == "KAMDReporter":
        from .reporter import KAMDReporter

        return KAMDReporter
    if name in ("config_from_amber", "parse_namelist", "amber_mask_indices"):
        from . import amber_input

        return getattr(amber_input, name)
    raise AttributeError(name)
