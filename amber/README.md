# KAMD for Amber (pmemd.cuda)

This is the reference implementation of KAMD. All KAMD work runs on the GPU: the momentum monitor, the first-shell
search (it reuses the PME neighbour-list cells) and the velocity redirect. The overhead is about 0.4 % relative to
plain pmemd.cuda.

## Requirements

* **Amber 26 source** (pmemd26 + AmberTools26). The installer was made and tested against the `Amber26_release`
  tag of the Amber git repository (2026-04-24). A licensed Amber 26 source distribution is required. This repository
  contains only the KAMD changes.
* **The usual Amber 26 build requirements** (CMake, a C/C++/Fortran compiler, CUDA).
  * Tested with CMake 3.31.6, GCC/gfortran 11.3.0 and CUDA 12.8 on RHEL 8, on RTX 4090 and RTX 5090 GPUs.
  * **CUDA is required:** KAMD exists only in `pmemd.cuda` (`pmemd.cuda_SPFP` / `pmemd.cuda_DPFP`). CPU `pmemd` and
    multi-GPU `pmemd.cuda.MPI` stop with an error if an input contains `&kamd`.

## Install

KAMD is added to an existing Amber 26 source tree by a small program, `kamd_install.py`, which edits the pmemd
source in place. It needs only Python 3.6 or newer (standard library) and changes nothing outside `src/pmemd/src/`.

```bash
# 1. unpack Amber 26 as usual (and apply official updates first, if you use them)
tar xjf AmberTools26.tar.bz2
tar xjf pmemd26.tar.bz2
cd pmemd26_src                       # the source directory that contains src/pmemd/

# 2. add KAMD to pmemd
python3 /path/to/kamd/amber/kamd_install.py --check .   # dry run: every edit must report [ ok ]
python3 /path/to/kamd/amber/kamd_install.py .           # install

# 3. build Amber with CUDA as usual. Either edit and run build/run_cmake, or:
mkdir build_kamd && cd build_kamd
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/amber26_kamd -DCOMPILER=GNU -DMPI=FALSE -DCUDA=TRUE \
         -DDOWNLOAD_MINICONDA=FALSE -DCHECK_UPDATES=FALSE -DAPPLY_UPDATES=FALSE
make -j 16 install
```

What the installer does:
* It applies 25 small edits to 11 pmemd files and adds one file (`src/pmemd/src/cuda/kKAMD.h`). Each edit is
  located by a few lines of unchanged Amber code that must occur exactly once (differences in indentation or
  trailing spaces are tolerated). The KAMD code it inserts is in [`kamd_src/blocks/`](kamd_src/blocks/), and the
  locations are listed in [`kamd_src/manifest.json`](kamd_src/manifest.json).
* All edits are checked before anything is written. If one cannot be placed (for example in a different Amber
  version, or after an Amber update that touched the same lines), it names the edit and changes nothing.
* The original files are kept in `src/pmemd/src/.kamd_backup/`. `src/pmemd/src/.kamd_install.json` records what
  was changed.
* `python3 kamd_install.py --uninstall .` restores the original files (then rebuild). It refuses if you have
  edited a patched file since; add `--force` to restore anyway.
* It refuses to run twice, or on a tree that already contains KAMD.
* With no directory argument it uses `$AMBERHOME`.

Verification: on a pristine `Amber26_release` source tree, the installer produces source files byte-identical to
the tree that was built and tested here (smoke test and validation below). Uninstalling restores the pristine
files exactly.

## Use

Add a `&kamd` namelist to an ordinary `mdin`:

```
&cntrl
   imin=0, irest=1, ntx=5, nstlim=500000, dt=0.002,
   ntt=3, gamma_ln=5, temp0=300, ntp=1, barostat=2, ntc=2, ntf=2,
/
&kamd
   pseudobd_cycles = 250,
   refer_mask(1) = ":1-77",
   refer_mask(2) = ":78-132",
   solvent_mask  = ":WAT",
/
```

The full namelist reference is in [`doc/kamd_tutorial.md`](doc/kamd_tutorial.md). It covers:
* reference groups and `refer_link`;
* transmembrane regions;
* ion-trajectory output;
* the experimental membrane shell mode.

## Test

`test/` holds a 29,467-atom solvated protein (`prmtop`, `restrt`) and an input with four reference groups.
`run_test.sh` runs 3,000 steps. The first reset typically comes after 1,000–2,500 steps (step 1,500 on an RTX 5090 with this input) (one `KAMD: velocity reset, group 1`
line per reset in `mdout`). The `KAMD-DIAG` mode (`export KAMD_DIAG=1`) dumps the first-shell assignment at the first
reset and exits. The OpenMM and GROMACS ports were validated against that dump.

## Known issues

* **Wrong energies at reset steps.** If an `ntpr` step coincides with a KAMD reset, the kinetic energy and temperature
  printed for that step are wrong (for example TEMP = −9294 K). The `A V E R A G E S` block is contaminated as well.
  The trajectory is not affected. Average over the other print steps, or compute temperatures from the trajectory.
* **One reference group.** With a single reference group the relative momentum is zero up to round-off, so resets
  are driven by floating-point noise. Use two or more groups.
* **Multi-GPU and CPU builds.** KAMD has not been validated in `pmemd.cuda.MPI` (multi-GPU). It is rejected there,
  and in CPU `pmemd`.
* **Membrane shell mode.** The membrane shell mode (`kamd_shell_mode = 1`) is experimental; see the tutorial.

## How to cite

Please cite the KAMD method paper and the KAMD software; see
[How to cite](../README.md#how-to-cite) in the top-level README.
