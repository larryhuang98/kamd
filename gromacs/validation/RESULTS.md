# GROMACS KAMD validation (2026-10-27)

**System.** The pmemd KAMD test system (`ambers/kamd/test/kamd`), converted to
GROMACS with ParmEd (`prepare_system.py`). It has 29,467 atoms and 4 reference
groups (atoms 20–96, 285–339, 366–482, 568–611); the bath is atoms 623–29467.

**Software and hardware.** GROMACS 2026.4 + `kamd-gromacs-2026.4.patch`, CUDA
12.8, RTX 5090 nodes, `-ntmpi 1 -ntomp 8`.

## 1. First-shell search, atom by atom

`mdrun` was run with `diag = yes`, which dumps the coordinates, box and per-atom
result at the first reset. `parity_vs_reference.py` then reran the search on the
same coordinates with `kamd_openmm`. That implementation matches pmemd.cuda atom
by atom (see `../openmm/validation/RESULTS.md`).

| run | step | shell atoms | surface atoms | mismatches |
|---|---|---|---|---|
| CPU smoke test (`pseudobd-cycles 20`) | 99 | 1137 | 210 | 0 |
| GPU job 1107040 | 3499 | 1147 | 201 | 0 |
| GPU job 1107043 (OpenMP search) | 3499 | 1144 | 221 | 0 |

## 2. Reset statistics, 250,000 steps from the same restart

**Settings.** `sd` at 355 K with τ = 0.2 ps (the same friction as Amber
`gamma_ln = 5`), C-rescale barostat, h-bond constraints, 0.8 nm cutoff, PME.

| | resets | mean interval (steps) | median | sd |
|---|---|---|---|---|
| pmemd.cuda | 145 | 1724 | 1750 | 450 |
| GROMACS | 144 | 1731 | 1750 | 515 |

The interval distributions agree: KS p = 1.00, Mann-Whitney p = 0.84. ⟨T⟩ is
355.4 K. A second GROMACS run gave 141 resets. It differs from the first only
because GPU nonbonded accumulation is not bitwise deterministic.

## 3. Cost, `md` + v-rescale, 100,000 steps

Each case ran 3 times, alternating KAMD on and off (one RTX 5090 node):

| update | no KAMD | KAMD | overhead |
|---|---|---|---|
| GPU (`-update gpu`) | 1389 ns/day | 1386 ns/day | 0.2% |
| CPU (`-update cpu`) | 1050 ns/day | 1047 ns/day | 0.2% |

The run-to-run spread is about 1.5%, so the overhead is within noise.
* **Monitoring only** (a run that checks but never resets) costs nothing
  measurable. That check copies the velocities to the host every 250 steps.
* **Each reset** costs about 0.8 ms of host time with OpenMP: about 0.05 ms for
  the MB draw, 0.7 ms for the first-shell search and 0.03 ms for the redirect.
  For 51 resets that adds up to 39 ms per 100k steps.
* **Before OpenMP**, the MB draw and the search were single-threaded, which cost
  about 5 ms per reset and about 2% overall.

These ns/day values are not comparable with the pmemd/OpenMM numbers: the setup
differs (v-rescale, GPU update, GROMACS' buffered Verlet list).
