# Validation against pmemd.cuda (2026-10-27)

System: `ambers/kamd/test/kamd` (29,467 atoms, 622 solute atoms + TIP3P), 4 reference
groups (`@20-96`, `@285-339`, `@366-482`, `@568-611`), solvent `@623-29467`.
Settings: 355 K Langevin (γ = 5 ps⁻¹), NPT, 2 fs, SHAKE/HBonds, 8 Å cutoff,
`pseudobd_cycles = 250`, `nzero_crossings = 2`.
Hardware: one RTX 5090 node. Amber binary: `build_shell_balance` pmemd.cuda_SPFP
(branch `kamd_membrane_shell_balance`). OpenMM 8.3.1, CUDA, mixed precision.
Job script: `run_validation.sbatch`. Outputs: `amber29k/`.

## 1. First-shell search, atom by atom

pmemd.cuda was run with `KAMD_DIAG=1` to its first reset (step 1500), and the GPU
kernel's per-atom result was dumped. The same run was then stopped one step
earlier, which gives the coordinates and box that kernel saw. `kamd_openmm`
ran its search on those coordinates.

| group | shell Amber | shell OpenMM | surface Amber | surface OpenMM |
|---|---|---|---|---|
| 1 | 271 | 271 | 48 | 48 |
| 2 | 229 | 229 | 42 | 42 |
| 3 | 409 | 409 | 88 | 88 |
| 4 | 199 | 199 | 33 | 33 |
| total | 1108 | 1108 | 211 | 211 |

There were 0 mismatched atoms. The total of 1319 equals pmemd's `KAMD-DIAG: flagged=1319`.
The k-d tree and brute-force backends agree as well.

## 2. Dynamics, 250,000 steps (500 ps) from the same restart

| | resets | mean interval (steps) | median | sd | ⟨T⟩ (K) | ns/day |
|---|---|---|---|---|---|---|
| pmemd.cuda KAMD | 145 | 1724 | 1750 | 450 | 354.9* | 483 |
| OpenMM KAMD | 147 | 1689 | 1500 | 498 | 356.5 | 706 |
| OpenMM, no KAMD | – | – | – | – | 356.0 | 743 |

The two engines' reset-interval distributions are indistinguishable
(Mann-Whitney p = 0.35, KS p = 0.81). On average a reset turns 1088 shell atoms
and 201 surface atoms.

The ns/day values in this table are **not** a fair speed comparison. pmemd ran
`ntp=1` with the default Berendsen barostat, which computes the virial on every
step, while OpenMM used a Monte Carlo barostat. Switching pmemd to `barostat=2`
alone takes it from 484 to 671 ns/day (see `bench/`).

## 3. Speed with matched conditions (`run_fair.sbatch`, outputs in `fair/`)

Both engines used: an MC barostat every 25 steps; Amber's Ewald coefficient
(0.34864 Å⁻¹) and 64³ PME grid forced in OpenMM; COM removal every 1000 steps;
the same thermostat, cutoff, constraints and time step. Each case ran 100,000
steps on one RTX 5090 (node B).

| | no KAMD | KAMD | KAMD overhead | resets |
|---|---|---|---|---|
| pmemd.cuda (SPFP) | 688 ns/day | 685 ns/day | 0.4% | 56 |
| OpenMM 8.3.1 (CUDA, mixed) | 765 ns/day | 734 ns/day | 4.3% | 53 |

Each case was run once. Node-to-node variation looks like about 2–3%: pmemd
with the MC barostat gave 671 ns/day on node A and 688 on node B.

The KAMD overhead in OpenMM comes from two sources:
* each reset costs about 16 ms on the host (getState, the k-d tree search, then
  setVelocities);
* every check downloads the velocities.

An OpenMM plugin with device kernels would be needed to reach pmemd's 0.4%.

\* pmemd prints a bogus kinetic energy on any `ntpr` step that coincides with a
KAMD reset. Example: step 1500 reports TEMP = −9294 K. This happens because
runmd.F90 sets `reset_velocities = .true.`, which takes the
`si_kin_ene = etot_save - si_pot_ene` path. Those prints also corrupt the
"A V E R A G E S" block (TEMP = −918 K). 354.9 K is the mean of the other 434
prints. This is only a reporting bug and does not affect the dynamics.
