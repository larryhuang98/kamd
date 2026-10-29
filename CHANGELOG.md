# Changelog

## 0.1.0 (2026-10-29)

First public release.

* **amber/**: KAMD in pmemd.cuda for Amber 26.
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
