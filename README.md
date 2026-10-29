# KAMD — Kinetic Accelerated Molecular Dynamics

KAMD speeds up conformational sampling of biomolecules. It does not modify the
potential energy surface and needs no collective variables. Every
`pseudobd_cycles` MD steps, KAMD checks the translational momentum of user-defined
reference groups (protein domains, transmembrane strands, …) relative to their
common centre of mass. Once every group has reversed direction `nzero_crossings`
times, KAMD performs a pseudo-Brownian reset:

* all velocities are redrawn from Maxwell-Boltzmann at the target temperature;
* the first-shell solvent (or membrane) atoms and the surface atoms of each group
  are turned onto one random direction per group, each atom keeping its speed.

Method: H. Wei, *J. Chem. Inf. Model.* **66**, 4116–4129 (2026), [doi:10.1021/acs.jcim.5c03093](https://doi.org/10.1021/acs.jcim.5c03093).

This repository contains three implementations that share the same algorithm,
input names and defaults:

| folder | engine | form | GPU | status |
|---|---|---|---|---|
| [`amber/`](amber/) | Amber 26 `pmemd.cuda` | installer program that adds KAMD to the pmemd source | all KAMD work on the GPU, ~0.4 % overhead | reference implementation; production use |
| [`gromacs/`](gromacs/) | GROMACS 2026.4 `gmx mdrun` | patch to the GROMACS source | CPU + GPU update supported, ~0.2 % overhead | validated against Amber |
| [`openmm/`](openmm/) | OpenMM ≥ 8.2 | Python package (`pip install`) | resets on the host, ~4 % overhead | validated against Amber |

## How the ports were validated

All three engines were run on the same 29,467-atom solvated protein (Amber KAMD test system) at 355 K, NPT, 2 fs:

* **First-shell search, atom by atom.** pmemd.cuda dumped its GPU first-shell assignment at a reset. The OpenMM and
  GROMACS searches were then run on the same coordinates. Both matched pmemd for all 1108 shell and 211 surface atoms,
  with **0 mismatches** (GROMACS was checked in three independent runs).
* **Reset statistics.** Over 250,000 steps from the same restart, the reset counts and mean intervals were:

  | engine | resets | mean interval |
  |---|---|---|
  | pmemd | 145 | 1724 steps |
  | OpenMM | 147 | 1689 steps |
  | GROMACS | 144 | 1731 steps |

  The interval distributions are indistinguishable (KS p = 0.81 OpenMM vs pmemd; 1.00 GROMACS vs pmemd).
* **Overhead.** KAMD on vs off on one RTX 5090 under matched conditions:

  | engine | overhead |
  |---|---|
  | pmemd | 0.4 % |
  | GROMACS | 0.2 % (CPU or GPU update) |
  | OpenMM | 4.3 % |

Details and scripts are in `openmm/validation/` and `gromacs/validation/`.

## Requirements at a glance

| | Amber | GROMACS | OpenMM |
|---|---|---|---|
| engine version | Amber 26 (`Amber26_release`) | 2026.4 | ≥ 8.2 (tested 8.2, 8.3.1, 8.6.1) |
| build tools | as for Amber 26 (tested CMake 3.31, GCC/gfortran 11.3) | CMake ≥ 3.28, GCC ≥ 11 | Python 3.10–3.13 |
| CUDA | required for KAMD (tested 12.8) | optional (≥ 12.1; tested 12.8) | optional |
| tested on | RHEL 8, GCC 11.3, CUDA 12.8, RTX 4090/5090 | RHEL 8, GCC 11.3, CUDA 12.8, RTX 4090/5090 | Linux x86-64 |

## Citing

If you use KAMD, please cite the method paper and this software (see [`CITATION.cff`](CITATION.cff); the
"Cite this repository" button on GitHub gives APA and BibTeX):

```bibtex
@article{Wei2026KAMD,
  author  = {Wei, Haixin},
  title   = {Enhanced Sampling on Domain/Motif Level with Kinetic Accelerated Molecular Dynamics},
  journal = {Journal of Chemical Information and Modeling},
  volume  = {66},
  number  = {7},
  pages   = {4116--4129},
  year    = {2026},
  doi     = {10.1021/acs.jcim.5c03093}
}
```

## Licenses

Each folder carries the license of the engine it modifies or builds on:

* `gromacs/`: **LGPL-2.1-or-later**, the license of GROMACS.
* `openmm/`: **MIT**.
* `amber/`: the KAMD code and installer are made available under the same terms as the Amber pmemd source they
  modify. They only work with a licensed Amber source distribution, which this repository does not contain.

See [`LICENSE.md`](LICENSE.md).
