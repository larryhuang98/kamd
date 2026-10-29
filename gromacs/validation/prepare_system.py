"""Convert the Amber KAMD test system to GROMACS and write the KAMD index/input files."""
import os

import numpy as np
import parmed
from openmm import app, unit

SRC = os.environ.get("AMBER_TEST", "../../amber/test")  # amber/test of this repository
s = parmed.load_file(f"{SRC}/prmtop")
crd = app.AmberInpcrdFile(f"{SRC}/restrt")
s.coordinates = crd.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
L = [v.value_in_unit(unit.angstrom) if unit.is_quantity(v) else v for v in crd.getBoxVectors()]
L = np.array([[float(c) for c in (v.value_in_unit(unit.angstrom) if unit.is_quantity(v) else v)] for v in
              (crd.getBoxVectors().value_in_unit(unit.angstrom) if unit.is_quantity(crd.getBoxVectors()) else crd.getBoxVectors())])
s.box = [L[0, 0], L[1, 1], L[2, 2], 90, 90, 90]
s.velocities = crd.getVelocities(asNumpy=True).value_in_unit(unit.angstrom / unit.picosecond)
s.save("system.top", format="gromacs", overwrite=True)
s.save("system.gro", overwrite=True)

# same groups as the Amber test: @20-96, @285-339, @366-482, @568-611, bath @623-29467 (1-based)
groups = {"dom1": (20, 96), "dom2": (285, 339), "dom3": (366, 482), "dom4": (568, 611), "bath": (623, len(s.atoms))}
with open("index.ndx", "w") as f:
    for name, (a, b) in groups.items():
        f.write(f"[ {name} ]\n")
        ids = list(range(a, b + 1))
        for k in range(0, len(ids), 15):
            f.write(" ".join(map(str, ids[k:k + 15])) + "\n")
print("atoms", len(s.atoms), "box (A)", s.box[:3])
