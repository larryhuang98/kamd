# Changelog

## 0.2.0 (unreleased)

* **amber/**: whole-lipid membrane shell modes, `kamd_shell_mode = 2` (recommended), 3 and 4. Lipids touching a
  transmembrane group are turned as whole molecules, as KAMD does for first-shell water. Lipid molecules come from
  the prmtop molecule table instead of residues, so split-residue lipids (Amber Lipid21) are handled correctly.
* **amber/**: fix the build with CUDA versions that still target compute capability < 6.0 (e.g. CUDA 12.3, which
  builds sm_50 to sm_90): double-precision atomic adds in the membrane kernels now use a fallback there. In
  0.1.0 the membrane shell kernel did not compile with such toolkits. Tested with CUDA 12.3 / GCC 10.2 and
  CUDA 12.8 / GCC 11.3.
* **amber/**: removed the `KAMD_DBG RESET FIRED` debug line from `mdout`; each reset is still reported by the
  `KAMD: velocity reset, group N, step` lines.
* Documentation: the momentum criterion is stated precisely (projection onto the common centre-of-mass velocity);
  OpenMM and GROMACS README files say that only the per-atom membrane shell (mode 0) is ported; removed references
  to files that are not in this repository.
* Documentation: nonbonded settings for CHARMM36 lipids (`cut=12, fswitch=10, vdwmeth=0`); NPT equilibration
  followed by NVT KAMD for membrane proteins.

## 0.1.0 (2026-10-29)

First public release.

* **amber/**: KAMD in pmemd.cuda for Amber 26, added to an Amber source tree by `kamd_install.py`
  (dry-run check, install, uninstall; no Amber source is distributed).
  * mask-based reference groups;
  * `refer_link`;
  * solution and transmembrane regions;
  * GPU first-shell search and redirect;
  * ion-trajectory output;
  * experimental molecular membrane shell with an area gate.

  `&kamd` is now rejected outside single-GPU pmemd.cuda; before this, CPU pmemd ignored it silently.
* **gromacs/**: `gmx mdrun -kamd` for GROMACS 2026.4.
  * leap-frog `md` and `sd`, CPU and GPU update, single PP rank;
  * index-group input;
  * GoogleTest unit tests (`mdrun-kamd-test`).
* **openmm/**: `kamd_openmm` 0.1.0, an OpenMM reporter.
  * reads Amber `&kamd` inputs directly (with ParmEd);
  * k-d tree first-shell search;
  * pytest suite.

The OpenMM and GROMACS ports were developed with AI assistance (Claude, Anthropic). The authors reviewed them and
validated them against the reference Amber implementation.
