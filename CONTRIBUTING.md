# Contributing

Bug reports and pull requests are welcome through the GitHub issue tracker.

**Reporting a problem.** Please include:
* which engine you used and its version;
* the KAMD input (`&kamd` namelist, `kamd.dat`, or the `KAMDConfig`);
* the log or error message;
* if you can, a small system that reproduces it.

**Changing the algorithm.** All three implementations must stay equivalent. A change to the momentum monitor, the
first-shell search or the redirect rule needs:
1. the same change in every engine, or a clear statement of why it applies to one engine only;
2. a passing first-shell parity check against pmemd (`openmm/validation/parity_first_shell.py`,
   `gromacs/validation/parity_vs_reference.py`);
3. updated unit tests.

**Before opening a pull request, run the tests:**

```bash
# OpenMM
cd openmm && pip install -e .[amber,test] && pytest tests && ruff check . && ruff format --check .
# GROMACS (after building with the patch)
make mdrun-kamd-test && bin/mdrun-kamd-test
```

GROMACS code follows the GROMACS style. Format it with clang-format 18 and the `.clang-format` from the GROMACS
source tree.
