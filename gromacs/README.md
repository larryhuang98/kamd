# KAMD for GROMACS

A patch for GROMACS 2026.4 that adds KAMD to `gmx mdrun` through a new `-kamd` option.
It is a port of the reference pmemd.cuda implementation (`../amber`), and it uses the same input names and defaults,
with lengths in nm.

## Requirements

GROMACS 2026.4 itself requires:

* CMake ≥ 3.28;
* a C++17 compiler (GCC ≥ 11, or a recent Clang);
* for GPU builds, CUDA ≥ 12.1.

The KAMD code uses only standard C++17 and OpenMP and has no additional dependencies.

Tested:

| | |
|---|---|
| OS | RHEL 8 |
| compiler | GCC 11.3.0 |
| build tools | CMake 3.31.6 |
| CUDA | 12.8; both `GMX_GPU=CUDA` and `GMX_GPU=OFF` builds |
| GPUs | RTX 4090, RTX 5090 |
| CI | a GitHub Actions workflow (ubuntu-24.04, GCC 13, CPU build) is included but has not run yet |

## Install

```bash
curl -O https://ftp.gromacs.org/gromacs/gromacs-2026.4.tar.gz
tar xzf gromacs-2026.4.tar.gz
cd gromacs-2026.4
patch -p1 < /path/to/kamd/gromacs/kamd-gromacs-2026.4.patch     # 9 files, no fuzz expected
mkdir build && cd build
cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
         -DGMX_GPU=CUDA -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2026.4-kamd
make -j 16 && make install
make mdrun-kamd-test && bin/mdrun-kamd-test   # optional: KAMD unit tests (14 tests)
```

**Set the compilers explicitly.** On many clusters the default `cc`/`c++` is older than GCC 11, and CMake stops with
"GCC version 11 or later required". Add whatever other options you normally use (MPI, SIMD level, `GMX_GPU=OFF`, …).

The patch makes these changes:
* it adds `src/gromacs/mdrun/kamd.{h,cpp}` and their unit tests (`src/gromacs/mdrun/tests/`);
* it adds the KAMD hook to the leap-frog loop in `md.cpp` (about 60 lines);
* it registers the `-kamd` option;
* it adds checks that refuse `-kamd` where it would otherwise be ignored.

The `.tpr` format and `grompp` are unchanged.

## Run

```bash
gmx mdrun -s topol.tpr -kamd kamd.dat -ntmpi 1      # -update gpu works too
```

`kamd.dat` (`;` or `#` starts a comment; `-` and `_` are interchangeable in keys):

```
index           = index.ndx              ; required: groups below are read from here
refer-groups    = dom1 dom2 dom3 dom4    ; required: reference groups (>= 2)
solvent-group   = SOL_ION                ; required: solvent bath
pseudobd-cycles = 250                    ; steps between momentum checks
nzero-crossings = 2                      ; zero crossings per group before a reset
shell-cut       = 0.5                    ; first-shell assignment cutoff (nm)
surface-cut     = 0.3                    ; surface reference-atom cutoff (nm), <= shell-cut
refer-link      = 0 1 0 0                ; optional, 1-based: group i shares group k's direction
region-type     = 1 1 1 2                ; optional: 1 solution, 2 transmembrane
membrane-group  = MEMB                   ; required if any region-type is 2
membrane-axis   = z                      ; membrane normal (x, y or z)
membrane-cut    = 0.3                    ; transmembrane cutoff (nm)
temperature     = 310                    ; default: ref-t (all tc-groups must agree)
seed            = -1                     ; default: ld-seed from the .tpr
log             = kamd.log               ; one line per reset
diag            = no                     ; yes: dump coordinates + first shell at the first reset
```

`examples/kamd.dat` and `validation/` contain a complete setup for the 29,467-atom test system.

### Supported setups

* `integrator = md` or `sd` (leap-frog). With `md`, GPU update (`-update gpu`) works.
* One PP rank: `-ntmpi 1`, or `-ntmpi 2 -npme 1`. Multi-rank domain decomposition is not supported yet.
* Rectangular and triclinic boxes. Rectangular boxes use a cell search; triclinic boxes use a brute-force
  minimum-image search with the same convention as pmemd.

`mdrun` stops with an error, rather than ignoring `-kamd`, in these cases:
* any other integrator: minimization, `bd`, `md-vv`, TPI, rerun;
* multiple PP ranks;
* CUDA graphs (`GMX_CUDA_GRAPH`);
* the modular simulator (`GMX_USE_MODULAR_SIMULATOR`).

### Implementation notes

* **Where the hook runs.** After the forces at *t* and before the leap-frog update, which is the same point as pmemd
  and the place where GROMACS applies Andersen randomization. The stored velocities are v(t − dt/2). The redirect acts
  on v(t) = v(t − dt/2) + F dt/2m and converts back.
* **Random numbers.** ThreeFry with counters (seed, step, global atom index). Results are reproducible, restart-safe,
  and independent of thread count and atom ordering; a unit test checks this.
* **GPU update.** Velocities are copied to the host at each check (every `pseudobd-cycles` steps). Coordinates and
  forces are copied at a reset, and the new velocities are copied back. The measured cost is within noise.
* **Monitor counters.** They restart in every `mdrun` invocation, checkpoint continuations included, as in pmemd.
* **Energies at a reset step.** The leap-frog kinetic energy reported at that step averages a pre-reset and a
  post-reset half step.

## Validation (`validation/`, details in `validation/RESULTS.md`)

* **First-shell search vs pmemd.** At the first reset the coordinates and first shell were dumped, and the search
  was re-run with the OpenMM port, which matches pmemd atom by atom. 0 mismatches in three runs (1137–1147 shell
  atoms, 201–221 surface atoms).
* **Reset statistics.** Over 250,000 steps (Langevin, 355 K, NPT) from the pmemd restart: 144 resets, mean interval
  1731 steps. pmemd gave 145 resets and 1724 steps. KS p = 1.00.
* **Cost.** KAMD on vs off, 3 alternating repeats each:

  | update | overhead |
  |---|---|
  | GPU | 1389 → 1386 ns/day (0.2 %) |
  | CPU | 1050 → 1047 ns/day (0.2 %) |

## How to cite

Please cite the KAMD method paper and the KAMD software; see
[How to cite](../README.md#how-to-cite) in the top-level README.
