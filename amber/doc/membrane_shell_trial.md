# Molecular membrane shell trial

This is an opt-in KAMD experiment for membrane proteins. The default
`kamd_shell_mode=0` keeps the existing atom-wise membrane redirect.

The motivation is that atom-wise velocity redirection can shear a lipid and
create repeated, coherent lateral impulses. Lipids in MD rearrange slowly, so
these impulses may accumulate before the bilayer can relax. The trial keeps
KAMD's protein and local environment coupling while reducing the membrane
impulse applied at each reset.

With `kamd_shell_mode=1`, each complete membrane residue is treated as one
lipid. For each lipid, the masses of its atoms in each transmembrane KAMD
first shell are summed. The group with the largest contacted mass owns the
lipid. An existing owner is retained while its contacted mass is at least
80% of the largest score. The lipid's membrane-plane center-of-mass velocity
is turned toward that group's random direction, capped by
`kamd_shell_max_angle` and scaled by the square root of its contacted mass
fraction. The magnitude of this velocity and the normal component are
preserved. The molecular turn applies the same velocity change to every atom
in the lipid, preserving their relative velocities after KAMD's usual global
Maxwell-Boltzmann resampling. After all molecular turns, a uniform in-plane
velocity correction removes the net momentum those turns added to the
membrane. This leaves the nearby lipids moving with the protein while keeping
the membrane's total in-plane momentum at its pre-turn value. The correction
does not change any lipid's internal relative velocities or normal velocity.

Two optional reset gates can give the membrane more time between impulses:

- `kamd_shell_relax_steps` is the minimum number of MD steps between resets.
- `kamd_shell_area_floor` defers resets when membrane-plane box area falls
  below this fraction of the reference area.
- `kamd_shell_area_reference` sets a fixed reference area in Å². Use the
  original equilibrated box area for every segment of a long restart chain;
  zero (the default) uses each segment's own starting area.

Example for 2 fs steps and a z-normal bilayer:

```
&kamd
  pseudobd_cycles = 250,
  ... reference masks, region types, and membrane mask ...
  membrane_axis = 3,
  kamd_shell_mode = 1,
  kamd_shell_max_angle = 45.0,
  kamd_shell_relax_steps = 5000,
  kamd_shell_area_floor = 0.95,
  kamd_shell_area_reference = 10025.9872,
/
```

The area gate is a coarse safety measure, not a local stress measurement.
This mode requires the membrane mask to contain complete lipid residues, one
lipid per residue. A speed-preserving COM turn may also give less local clearance to
the moving protein than the original atom-wise redirect. Compare reset count,
protein motion, box area, bilayer thickness, height fluctuations, local lipid
packing, and crash rate against the legacy mode before using this in a long
production run.

The VDAC smoke tests compared both modes from the same initial structure and random seed.
