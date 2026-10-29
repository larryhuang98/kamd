# kamd-openmm

KAMD (Kinetic Accelerated Molecular Dynamics) for [OpenMM](https://openmm.org).

This is a port of the pmemd.cuda KAMD implementation by Zhen Huang, Yong Duan,
Ray Luo and Haixin Wei. Method paper: H. Wei *et al.*, bioRxiv 2025,
doi:10.1101/2025.06.14.659665.

KAMD speeds up conformational sampling without changing the potential and without
collective variables. The user defines reference groups (domains of the
macromolecule). Every `pseudobd_cycles` steps the code checks the momentum of
each group relative to the reference COM, projected on the COM velocity. Once
every group has crossed zero `nzero_crossings` times, it does a pseudo-Brownian
reset:

1. All velocities are redrawn from Maxwell-Boltzmann at the target temperature.
2. Each group draws a random unit vector. `refer_link` groups share one vector.
   Transmembrane groups get an in-plane vector.
3. The search finds first-shell bath atoms and surface reference atoms.
   * Solution groups search the solvent bath, using `shell_cut` for assignment
     and `surface_cut` for surface marking.
   * Transmembrane groups search the membrane bath, using `membrane_cut` for both.
   * Each bath atom goes to the group of its nearest reference atom.
4. The velocities of those atoms are turned onto their group's direction, keeping
   each atom's speed. Transmembrane groups keep the component along the membrane
   normal.

## Requirements

* Python 3.10–3.13
* OpenMM ≥ 8.2 (the oldest release with PyPI wheels). Conda-forge builds work as well.
* numpy and scipy. ParmEd is optional; it is needed only to read Amber inputs and masks.

Tested in fresh virtual environments:

| Python | OpenMM |
|---|---|
| 3.10 | 8.3.1 |
| 3.11 | 8.3.1 |
| 3.12 | 8.2.0 |
| 3.12 | 8.3.1 |
| 3.13 | 8.3.1 |
| 3.12 | 8.6.1 (conda-forge) |

All 24 tests pass in each.

## Install

From a clone of this repository:

```bash
pip install ./openmm             # numpy, scipy, openmm
pip install "./openmm[amber]"    # + parmed, to read Amber &kamd inputs and masks
```

Or directly from GitHub:

```bash
pip install "kamd-openmm[amber] @ git+https://github.com/larryhuang98/kamd#subdirectory=openmm"
```

For CUDA, install OpenMM with the CUDA wheels (`pip install "openmm[cuda12]"`) or install it from conda-forge.

## Use

```python
from kamd_openmm import KAMDConfig, KAMDReporter

cfg = KAMDConfig(
    refer_groups=[domain1_atoms, domain2_atoms],  # 0-based atom indices
    solvent_atoms=water_atoms,
    temperature=300.0,  # K, same as the thermostat
    pseudobd_cycles=250,
    seed=1,
)
simulation.reporters.append(KAMDReporter(cfg, log="kamd.log"))
simulation.step(5_000_000)
```

You can also reuse an Amber input as is:

```python
from kamd_openmm import config_from_amber

cfg = config_from_amber("md.in", "system.prmtop")  # reads &kamd and temp0 from &cntrl
```

### Parameters (same names as the Amber `&kamd` namelist)

| parameter | default | meaning |
|---|---|---|
| `refer_groups` | required | list of atom-index lists, one per reference group (max: no limit) |
| `solvent_atoms` | required | solvent bath atoms |
| `temperature` | required | Maxwell-Boltzmann temperature, K |
| `pseudobd_cycles` | 250 | steps between momentum checks |
| `nzero_crossings` | 2 | zero crossings per group before a reset |
| `shell_cut` / `surface_cut` | 5.0 / 3.0 Å | first-shell assignment / surface marking cutoffs |
| `refer_link` | all 0 | 1-based: group *i* shares group `refer_link[i]`'s direction |
| `region_type` | all 1 | 1 = solution, 2 = transmembrane |
| `membrane_atoms`, `membrane_axis` (1/2/3), `membrane_cut` | –, –, 3.0 Å | transmembrane setup |
| `seed` | None | RNG seed for the MB draw and directions |

`KAMDReporter` options:
* `velocity_offset` is `'auto'` by default. Leapfrog-family integrators (`LangevinMiddle`,
  `Langevin`, `Verlet`) store velocities at *t − dt/2*. The redirect is applied to
  the synchronous velocity *v + F dt/2m* and then converted back, as in pmemd.
  Velocity-Verlet integrators (`NoseHoover`) need no correction. A `CustomIntegrator`
  needs `'half'` or `'full'` set explicitly.
* `search` is `'auto'` by default: a periodic k-d tree for rectangular boxes, and a
  brute-force minimum-image search for triclinic boxes.
* `reporter.state_dict()` / `load_state_dict()` save and restore the monitor counters
  and RNG state when a run continues from a checkpoint.

### Differences from pmemd.cuda

* Only the pseudo-BD algorithm is ported. Not ported:
  * the experimental molecular membrane shell (`kamd_shell_mode`);
  * the ion-trajectory output (`ion_traj_*`). OpenMM reporters already cover that.
* Masks do not need to be contiguous, and there is no limit on the number of groups.
* The first-shell search runs on the host. It needs one host transfer of
  positions, velocities and forces per reset (about every 1000–2000 steps), plus a
  velocity download per check.
* With a single reference group the relative momentum is zero up to round-off, so
  resets are driven by floating-point noise. pmemd uses the same formula. Use two
  or more groups.

## Tests and validation

```bash
pip install "./openmm[amber,test]"
pytest openmm/tests     # algorithm unit tests + small OpenMM water-box tests (CPU, ~5 s)
sbatch validation/run_validation.sbatch   # parity with pmemd.cuda on the 29,467-atom KAMD test system
```
