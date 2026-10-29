# KAMD — Kinetic Accelerated Molecular Dynamics

[![DOI](https://zenodo.org/badge/1438016791.svg)](https://doi.org/10.5281/zenodo.23755858)

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
See [How to cite](#how-to-cite).

This repository contains three implementations that share the same algorithm,
input names and defaults:

| folder | engine | form | GPU | status |
|---|---|---|---|---|
| [`amber/`](amber/) | Amber 26 `pmemd.cuda` | installer program that adds KAMD to the pmemd source | all KAMD work on the GPU, ~0.4 % overhead | reference implementation; production use |
| [`gromacs/`](gromacs/) | GROMACS 2026.4 `gmx mdrun` | patch to the GROMACS source | CPU + GPU update supported, ~0.2 % overhead | validated against Amber |
| [`openmm/`](openmm/) | OpenMM ≥ 8.2 | Python package (`pip install`) | resets on the host, ~4 % overhead | validated against Amber |

## Using KAMD

### Workflow

1. **Prepare and equilibrate the system as for ordinary MD.** KAMD starts from an equilibrated NVT/NPT state
   and uses your normal thermostat and barostat. It changes no force-field term.
2. **Choose reference groups.** These are the structural units whose *relative* motion you want to sample:
   protein domains, motifs, subunits, or transmembrane helices/strands. Use **two or more** groups; the groups
   must not overlap.
3. **Choose the bath.** The solvent around the groups (water + ions) and, for membrane proteins, the
   membrane. KAMD redirects bath atoms that touch a group (the "first shell") together with the group's own
   surface atoms.
4. **Add the KAMD input**: a `&kamd` namelist in the Amber `mdin`, a `kamd.dat` file for `gmx mdrun -kamd`,
   or a `KAMDConfig` for OpenMM (examples below). The defaults are a good start; usually only the groups and
   the bath have to be set.
5. **Run, and check that resets happen** (see [Checking that it works](#checking-that-it-works)).
6. **Analyse the conformations.** The potential is unchanged, so the conformations visited are those of the
   unbiased force field. The time axis, however, is accelerated by the resets: do not read rates, diffusion
   constants or ion currents directly from a KAMD trajectory. To measure such quantities, start plain MD from
   KAMD snapshots.

### Minimal input in each engine

Amber (`mdin`, together with the usual `&cntrl`):

```
&kamd
   pseudobd_cycles = 250,          ! 0 = KAMD off
   refer_mask(1)   = ":1-77",      ! reference group 1 (AmberMask)
   refer_mask(2)   = ":78-132",    ! reference group 2
   solvent_mask    = ":WAT,Na+,Cl-",
/
```

GROMACS (`gmx mdrun -s topol.tpr -kamd kamd.dat -ntmpi 1`):

```
index           = index.ndx       ; groups below are read from this index file
refer-groups    = dom1 dom2
solvent-group   = SOL_ION
pseudobd-cycles = 250
```

OpenMM:

```python
from kamd_openmm import KAMDConfig, KAMDReporter

cfg = KAMDConfig(refer_groups=[dom1_atoms, dom2_atoms],   # 0-based atom indices
                 solvent_atoms=solvent_atoms, temperature=300.0, pseudobd_cycles=250)
simulation.reporters.append(KAMDReporter(cfg, log="kamd.log"))
# or reuse an Amber input: cfg = kamd_openmm.config_from_amber("md.in", "system.prmtop")
```

### Parameters

The same names are used in all three engines (GROMACS writes `-` for `_` and uses nm instead of Å).

| Parameter (Amber / OpenMM) | GROMACS `kamd.dat` | Default | Meaning |
|---|---|---|---|
| `pseudobd_cycles` | `pseudobd-cycles` | Amber 0 (off); GROMACS, OpenMM 250 | MD steps between momentum checks. 250 = 0.5 ps at 2 fs. In Amber, 0 switches KAMD off. |
| `nzero_crossings` | `nzero-crossings` | 2 | How many times every group's relative momentum must change sign before a reset fires. Larger = fewer, more widely spaced resets. |
| `refer_mask(i)` / `refer_groups` | `refer-groups` | required | The reference groups: AmberMask per group in Amber (slots 1..N without gaps), group names from the index file in GROMACS, lists of atom indices in OpenMM. Amber allows up to 32 groups. |
| `solvent_mask` / `solvent_atoms` | `solvent-group` | required | Solvent bath (water and ions). In Amber it must be one contiguous atom range. |
| `shell_cut` | `shell-cut` | 5.0 Å (0.5 nm) | A bath atom within this distance of a group belongs to that group's first shell (it goes to the nearest group). |
| `surface_cut` | `surface-cut` | 3.0 Å (0.3 nm) | Group atoms within this distance of the bath are the group's surface atoms and are redirected too. Must be ≤ `shell_cut`. |
| `refer_link(i)` | `refer-link` | 0 | If set to *k* (1-based), group *i* gets the same random direction as group *k*, so linked groups are pushed together. 0 = own direction. |
| `region_type(i)` | `region-type` | 1 | 1 = solution group (in contact with solvent); 2 = transmembrane group (in contact with the membrane). |
| `membrane_mask` / `membrane_atoms` | `membrane-group` | – | Membrane bath. Required if any group has `region_type = 2`. Contiguous in Amber. |
| `membrane_axis` | `membrane-axis` | – | Membrane normal: 1/2/3 (Amber, OpenMM) or x/y/z (GROMACS). Required for transmembrane groups. |
| `membrane_cut` | `membrane-cut` | 3.0 Å (0.3 nm) | The single cutoff for transmembrane groups, used both for the membrane first shell and for the surface atoms. |
| temperature | `temperature` | Amber: `temp0`; GROMACS: `ref-t`; OpenMM: required | Temperature of the Maxwell-Boltzmann redraw. Must equal the thermostat temperature. |
| seed | `seed` | Amber: `ig`; GROMACS: `ld-seed`; OpenMM: `None` | Random seed for the redraw and the directions. |
| `ion_traj_outfreq`, `ion_mask`, `ion_traj_file` | – | 0 (off) | Amber only: write ion positions relative to the membrane every N steps, for offline crossing counts. Works without KAMD. |
| `kamd_shell_mode` and `kamd_shell_*` | – | 0 (off) | Amber only, **experimental**: molecular membrane shell (whole lipids are turned in-plane, with an optional area gate). See the [tutorial](amber/doc/kamd_tutorial.md). |

How to set them:

* **Reference groups matter most.** Choose units that move relative to each other in the transition you want
  to sample, for example the two lobes of a kinase, the domains of a multidomain protein, or the strands of a
  β-barrel. Rigid units of a few hundred to a few thousand atoms work well. Units that should move together
  can share a direction with `refer_link`. Do not use a single group: its relative momentum is zero up to
  round-off, so resets would be driven by numerical noise.
* **Solution vs. transmembrane.** Use `region_type = 2` for groups that sit in the bilayer. Their bath is the
  membrane, and the velocity component along `membrane_axis` is kept, so they are pushed only within the
  membrane plane. Solution and transmembrane groups can be mixed in one system.
* **`pseudobd_cycles` and `nzero_crossings`** set how often resets happen. The defaults (250, 2) gave about one
  reset every 1,700 steps (3.4 ps) for the 29,467-atom test protein. Increase either one for gentler sampling.
  If resets are much rarer than this, check that the groups really move relative to each other.
* **Cutoffs.** The defaults correspond to the first hydration shell of water. Change them only for unusual
  solvents or coarse-grained models.
* **Thermostat and barostat.** Keep your usual settings (Langevin/`sd`, Monte Carlo or c-rescale barostat). The
  redraw uses the thermostat temperature, so the temperature stays at its target.

### Checking that it works

* **Amber:** `mdout` lists the atoms matched by each mask at startup, then one line per group at each reset:
  `KAMD: velocity reset, group 1, step 12500`.
* **GROMACS:** `md.log` records the KAMD setup (groups, bath, cutoffs). `kamd.log` (the `log` key) has one
  line per reset and ends with the total number of resets and the time spent in them.
* **OpenMM:** the `log` file has one line per reset; `reporter.n_resets` counts them.
* Inputs are checked before the run starts. Overlapping groups, an empty mask, a missing bath, or an unsupported
  setup stop the run with an error that names the problem.

### Limits

* **Amber:** single-GPU `pmemd.cuda` only (CPU `pmemd` and `pmemd.cuda.MPI` reject `&kamd`). Solvent and
  membrane must each be one contiguous atom range. At most 32 groups.
* **GROMACS:** `integrator = md` or `sd`, one PP rank (`-ntmpi 1`, or `-ntmpi 2 -npme 1`). CPU and GPU update.
* **OpenMM:** any platform. Leap-frog and velocity-Verlet integrators are handled automatically; a
  `CustomIntegrator` needs `velocity_offset` set.
* **Energies printed at a reset step (Amber, GROMACS):** the kinetic energy and temperature reported for the
  step at which a reset fires mix velocities from before and after the reset (in Amber the printed value can
  be far off). Ignore those single values; the trajectory and the other steps are unaffected.

Full references: [Amber tutorial](amber/doc/kamd_tutorial.md) (all options, examples, error messages),
[GROMACS](gromacs/README.md#run), [OpenMM](openmm/README.md#use).

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

## How to cite

If you use KAMD in published work, please cite **both** the method paper and the software:

1. **Method:** H. Wei, Enhanced Sampling on Domain/Motif Level with Kinetic Accelerated Molecular Dynamics.
   *J. Chem. Inf. Model.* **2026**, *66* (7), 4116–4129.
   [doi:10.1021/acs.jcim.5c03093](https://doi.org/10.1021/acs.jcim.5c03093)
2. **Software:** Z. Huang, Y. Duan, R. Luo, H. Wei, KAMD: Kinetic Accelerated Molecular Dynamics for Amber,
   GROMACS and OpenMM. Zenodo, **2026**.
   [doi:10.5281/zenodo.23755858](https://doi.org/10.5281/zenodo.23755858)

Which software DOI to use:
* [10.5281/zenodo.23755858](https://doi.org/10.5281/zenodo.23755858) is the *concept* DOI. It covers all
  versions and always resolves to the latest release. Use it for general citations.
* Each release also has its own DOI, for example
  [10.5281/zenodo.23755860](https://doi.org/10.5281/zenodo.23755860) for v0.1.0. Use the version DOI (and
  state the version in your methods section) when you need to identify the exact code that produced your results.

Please also cite the MD engine you ran KAMD in (Amber, GROMACS or OpenMM), following that engine's own
citation instructions.

BibTeX:

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

@software{Huang2026KAMDsoftware,
  author    = {Huang, Zhen and Duan, Yong and Luo, Ray and Wei, Haixin},
  title     = {{KAMD: Kinetic Accelerated Molecular Dynamics for Amber, GROMACS and OpenMM}},
  publisher = {Zenodo},
  year      = {2026},
  doi       = {10.5281/zenodo.23755858},
  url       = {https://github.com/larryhuang98/kamd}
}
```

The same information is in [`CITATION.cff`](CITATION.cff); the "Cite this repository" button on the GitHub page
exports it as APA or BibTeX.

## Licenses

Each folder carries the license of the engine it modifies or builds on:

* `gromacs/`: **LGPL-2.1-or-later**, the license of GROMACS.
* `openmm/`: **MIT**.
* `amber/`: the KAMD code and installer are made available under the same terms as the Amber pmemd source they
  modify. They only work with a licensed Amber source distribution, which this repository does not contain.

See [`LICENSE.md`](LICENSE.md).
