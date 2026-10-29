# KAMD Tutorial

This is a user-facing tutorial for KAMD (Kinetically Accelerated Molecular
Dynamics) in pmemd.cuda.  For the algorithm and the GPU implementation
internals, see `kamd_implementation.tex`.

KAMD is opt-in.  If you do not include a `&kamd` namelist (or you set
`pseudobd_cycles = 0`), pmemd runs as a normal MD with no KAMD overhead.

---

## 1. What KAMD does, in one paragraph

KAMD accelerates conformational sampling by periodically resetting the
velocities of first-shell solvent (and surface reference atoms) to a fresh
Maxwell-Boltzmann sample, then redirecting them along a random unit vector
per reference group.  The kicks are speed-preserving (temperature is
preserved) and trigger only when each reference group has accumulated a
configurable number of momentum-projection zero-crossings.  In practice this
shakes the system out of basins without modifying the potential, requires
no collective variables, and adds <0.5% overhead vs standard GPU MD.

---

## 2. Quick start

Put this `&kamd` block in your `mdin` file alongside your `&cntrl`:

```
&cntrl
   imin=0, irest=1, ntx=5,
   nstlim=50000, dt=0.002,
   ntt=3, gamma_ln=5, temp0=300,
   ntp=1, ntc=2, ntf=2, ntb=2,
/

&kamd
   pseudobd_cycles = 250,

   refer_mask(1) = ":1-77",
   refer_mask(2) = ":78-132",

   solvent_mask  = ":WAT",
/
```

This defines two reference domains and tells KAMD to monitor them against the
WAT bath, with a momentum check every 250 MD steps.  Reset fires when both
groups have hit two zero-crossings.

---

## 3. Full namelist reference

### General controls

| Variable | Default | Meaning |
|---|---|---|
| `pseudobd_cycles` | 0 | MD steps between momentum checks. **0 disables KAMD.** Typical: 250 (= 0.5 ps with 2 fs `dt`). |
| `nzero_crossings` | 2 | Per-group zero-crossing budget before a reset fires. |
| `shell_cut` | 5.0 Å | First-shell solvent assignment cutoff (solution groups). |
| `surface_cut` | 3.0 Å | Surface ref-atom marking cutoff (solution groups). |

### Reference groups

| Variable | Default | Meaning |
|---|---|---|
| `refer_mask(i)` | "" | Per-group AmberMask, **packed** in slots 1..N (no gaps). |
| `refer_link(i)` | 0 | If `>0`, group `i` shares group `refer_link(i)`'s random direction. 0 = own direction. |
| `region_type(i)` | 1 | Region type for group i: `1` = solution, `2` = transmembrane. |

Up to **`KAMD_MAX_GROUPS = 32`** groups.

### Solvent / membrane

| Variable | Default | Meaning |
|---|---|---|
| `solvent_mask` | "" | AmberMask for the solvent bath. **Required when KAMD is active.** Must resolve to a *contiguous* atom-index range. |
| `membrane_mask` | "" | AmberMask for the membrane bath. **Required iff any** `region_type(i) == 2`.  Also must resolve to a contiguous range. |
| `membrane_axis` | 0 | Membrane normal: `1` = x, `2` = y, `3` = z.  Single global axis.  **Required iff any TM region.** |
| `membrane_cut` | 3.0 Å | Single TM cutoff used for **both** first-shell membrane and surface ref-atom marking. |

### Ion-trajectory output (independent of pseudo-BD)

| Variable | Default | Meaning |
|---|---|---|
| `ion_traj_outfreq` | 0 | MD steps between ion-position records. **0 disables it.** Typical: 5000 (= 10 ps at 2 fs `dt`). Works in plain MD — no `pseudobd_cycles` needed. |
| `ion_mask` | "" | AmberMask selecting the atoms to record (the ions). **Required when `ion_traj_outfreq > 0`.** May be discontiguous. |
| `ion_traj_file` | `kamd_ions.dat` | Output filename. |
| `membrane_mask` | "" | AmberMask for the membrane, used as the COM reference. **Required when `ion_traj_outfreq > 0`**; must resolve to a contiguous range. |

---

## 4. AmberMask syntax cheat sheet

Same syntax cpptraj/sander/parmed use.  Most useful forms:

| Form | Selects |
|---|---|
| `:1-77` | Residues 1 through 77, all atoms |
| `:WAT` | All residues named `WAT` |
| `:WAT,Na+,Cl-` | Multiple residue names |
| `:1-100 & !@H=` | Heavy atoms of residues 1–100 |
| `:250@CA,CB` | Just `CA`, `CB` of residue 250 |
| `@20-96` | Atoms 20–96 (raw atom indices) |
| `:1-77 \| :150-220` | Union of two residue ranges (discontiguous group) |
| `:POPC,POPE` | Both POPC and POPE lipids (typical bilayer) |

Reference groups can be **discontiguous**.  Solvent and membrane masks must
collapse to a contiguous atom-index range — the GPU first-shell kernel uses
`[start, end]` filters.  KAMD will tell you (with the resolved span and the
gap count) if a bath mask is non-contiguous.

---

## 5. Solution vs. transmembrane regions

### Solution (default, `region_type = 1`)

A solution group is "exposed to the solvent bath."  Its **first shell** is
all solvent atoms within `shell_cut` (5 Å); its **surface ref atoms** are
group atoms within `surface_cut` (3 Å) of any solvent atom.  Each pseudo-BD
reset:

1. Solvent atom velocities are randomized per group along the group's random
   direction `ĝ`, speed-preserving: `v_new = ĝ * |v|`.
2. Surface ref atoms get the same redirect.

The two cutoffs are independent; you can tune them for non-water solvents or
low-density boxes.

### Transmembrane (`region_type = 2`)

A TM group is "embedded in the membrane bath."  Its first shell is
**membrane atoms** within `membrane_cut` of any TM ref atom; its surface ref
atoms are TM ref atoms within `membrane_cut` of any membrane atom.  Note:
TM uses **one cutoff for both** (3 Å by default).

Redirect rule for TM atoms:

- The component of velocity **along `membrane_axis`** is preserved (no
  cross-bilayer kick).
- The 2D perpendicular component is randomized to a random in-plane
  direction, preserving perpendicular speed.

This matches what membrane proteins want: lateral diffusion in the bilayer
plane is sampled stochastically, but the protein doesn't get punted through
the membrane.

#### Experimental molecular membrane shell (`kamd_shell_mode = 1`)

The default (`0`) retains the original atom-wise KAMD behavior. Mode `1` is
an opt-in trial for lipid membranes. A complete membrane residue is treated as
one molecule. The sum of contacted atom masses selects its reference group;
an existing owner is retained while its contacted mass remains within 80% of
the best group. Only the lipid's in-plane center-of-mass velocity is turned
toward the same group direction as the protein surface. This molecular turn
preserves the internal and membrane-normal velocities after the standard KAMD
Maxwell-Boltzmann resampling. The turn is capped per reset and scaled down
when only a small fraction of the lipid contacts the reference group.

```
&kamd
  ...
  kamd_shell_mode        = 1,
  kamd_shell_max_angle   = 45.0,   ! degrees; 180 allows full alignment
  kamd_shell_relax_steps = 50000,  ! minimum 100 ps at dt=0.002 ps
  kamd_shell_area_floor  = 0.95,   ! defer resets below 95% of starting xy area
/
```

The area gate is a coarse safety signal, not a measure of local lipid stress.
It uses the membrane-plane area at the start of each run as its reference. Resets resume
when the interval has passed and the area has recovered. Use an equilibrated,
healthy starting structure. The shell mode assumes each lipid is represented
by one complete topology residue in `membrane_mask`; a partial residue is an
error. Compare protein transitions and membrane structure against the default
mode before using this experimental mode for production.

You can mix solution and TM groups in the same simulation freely (e.g., a
soluble domain plus a TM helix on the same protein).

---

## 6. Examples

### 6a. Solution-only multi-domain protein

```
&kamd
   pseudobd_cycles = 250,

   refer_mask(1) = ":1-100",
   refer_mask(2) = ":120-220",
   refer_mask(3) = ":250-400",

   solvent_mask  = ":WAT",
/
```

### 6b. Linked groups (correlated motion)

Groups 2 and 3 share group 1's random direction, so all three domains kick
in the same direction at every reset:

```
&kamd
   pseudobd_cycles = 250,

   refer_mask(1) = ":1-100",
   refer_mask(2) = ":150-250",
   refer_mask(3) = ":300-400",
   refer_link(:) = 0, 1, 1,

   solvent_mask  = ":WAT",
/
```

### 6c. Discontiguous reference groups

A binding loop spread across the protein (residues 30-45 plus a beta-hairpin
at 180-195), restricted to backbone Cα/Cβ:

```
&kamd
   pseudobd_cycles = 250,

   refer_mask(1) = "(:30-45 | :180-195) & @CA,CB",

   solvent_mask  = ":WAT",
/
```

### 6d. Mixed solvent (water + ions)

Solvent mask covers WAT plus counterions, as long as they sit in a
contiguous atom-index block in the prmtop (typical for tleap-built systems
where ions follow water):

```
&kamd
   pseudobd_cycles = 250,
   refer_mask(1)   = ":1-300",
   solvent_mask    = ":WAT,Na+,Cl-",
/
```

### 6e. Transmembrane protein (GPCR-like)

Soluble extracellular domain (residues 1-77), TM helix bundle (residues
78-300), all within a POPC bilayer:

```
&kamd
   pseudobd_cycles = 250,
   shell_cut       = 5.0,
   surface_cut     = 3.0,
   membrane_cut    = 3.0,

   refer_mask(1)  = ":1-77",
   refer_mask(2)  = ":78-300 & !@H=",
   region_type(:) = 1, 2,

   solvent_mask   = ":WAT,Na+,Cl-",
   membrane_mask  = ":POPC,POPE",
   membrane_axis  = 3,
/
```

### 6f. Tuning cutoffs

```
&kamd
   pseudobd_cycles = 250,
   shell_cut       = 6.0,    ! looser first shell
   surface_cut     = 3.5,    ! looser surface marking
   nzero_crossings = 3,      ! wait longer between resets
   refer_mask(1)   = ":1-77",
   solvent_mask    = ":WAT",
/
```

### 6g. Many groups (up to 32)

Useful for highly modular assemblies — e.g. each subunit of a homo-octamer
as its own group:

```
&kamd
   pseudobd_cycles = 250,
   refer_mask(1)  = ":1-150",
   refer_mask(2)  = ":151-300",
   refer_mask(3)  = ":301-450",
   refer_mask(4)  = ":451-600",
   refer_mask(5)  = ":601-750",
   refer_mask(6)  = ":751-900",
   refer_mask(7)  = ":901-1050",
   refer_mask(8)  = ":1051-1200",
   solvent_mask   = ":WAT",
/
```

---

## 7. Diagnostics in `mdout`

When KAMD is active, the master output contains, after the &cntrl echo:

```
KAMD refer_mask(1) = "..." matches N atoms
KAMD refer_mask(2) = "..." matches M atoms
...
KAMD solvent_mask = "..." -> contiguous range [start, end]
KAMD membrane_mask = "..." -> contiguous range [start, end]    (only if TM)
```

At each reset event you will see one line per active group:

```
KAMD: velocity reset, group 1, step  N
KAMD: velocity reset, group 2, step  N
...
```

This confirms which step the kick fired and which groups participated.

---

## 8. Validation errors and what they mean

KAMD validates the entire input before any MD runs.  Errors print
`| ERROR:  &kamd: ...` and pmemd terminates.  Common ones:

| Error | Cause | Fix |
|---|---|---|
| `solvent_mask is required when pseudobd_cycles > 0` | Missing solvent mask | Set `solvent_mask = ":WAT"` (or whatever) |
| `at least one refer_mask(i) must be set` | No reference groups | Add at least one `refer_mask(i)` |
| `refer_mask entries must be packed; refer_mask(N) is empty` | Gap in indices (e.g., set `(1)` and `(3)` skipping `(2)`) | Renumber so masks fill 1..N consecutively |
| `refer_mask(i) = "..." matched no atoms` | Mask resolves to empty set | Check syntax against `parmed`/`cpptraj` |
| `solvent_mask = "..." must resolve to a contiguous atom range` | Solvent atoms have gaps in their indices | Reorder solvent in the prmtop, or split the run |
| `surface_cut must be <= shell_cut` | Cutoffs are inverted | Solution `surface_cut` must be inside `shell_cut` |
| `refer_link(i) = N is out of range` | Linking to a nonexistent group | Range is `[0, kamd_ngroups]` |
| `refer_link(i) cannot reference its own group` | Self-loop | Set `refer_link(i)` to 0 or another group |
| `region_type(i) = N is invalid` | Not 1 or 2 | Use `1` (solution) or `2` (transmembrane) |
| `membrane_mask is required when any region_type == 2` | TM region without membrane bath | Add `membrane_mask` |
| `membrane_axis = N is invalid` | Not 1, 2, or 3 | `1=x`, `2=y`, `3=z` |
| `atom N is matched by both refer_mask(i) and refer_mask(j)` | Two reference groups overlap | Make the masks disjoint |
| `atom N is matched by both refer_mask(i) and solvent_mask` | A reference atom is also in the solvent set | Exclude one or the other |
| `atom N is matched by both solvent_mask and membrane_mask` | Solvent and membrane baths overlap | Make the two baths disjoint |
| `ion_mask is required when ion_traj_outfreq > 0` | Ion trajectory on but no ion mask | Set `ion_mask`, or set `ion_traj_outfreq = 0` |
| `ion_mask = "..." matched no atoms` | Ion mask resolves to empty | Check the mask syntax against `parmed`/`cpptraj` |
| `membrane_mask is required for ion trajectory output (membrane COM)` | Ion trajectory on, no membrane mask | Set `membrane_mask` (the COM reference) |

The error message always names the offending atom by 1-based index, which
makes it easy to find with `parmed` or by grepping the prmtop.

---

## 9. Hard limits and current scope

- Up to **32 reference groups** (`KAMD_MAX_GROUPS`).  Easy to lift further if
  needed; just bump the parameter.
- Solvent and membrane each must resolve to a single **contiguous** atom-index
  range.  Reference groups have no such restriction.
- Reference groups, solvent, and membrane (when used) must be **mutually
  disjoint** — no atom may belong to more than one of these.
- Only one global `membrane_axis` (one bilayer plane per system).
- Binding/unbinding mode is on the roadmap, not yet shipped.

---

## 10. Files of interest

| File | What it contains |
|---|---|
| `doc/kamd/kamd_tutorial.md` | This document |
| `doc/kamd/kamd_implementation.tex` | Algorithm and GPU implementation details |
| `src/pmemd/src/mdin_ctrl_dat.F90` | `&kamd` namelist, validation, mask resolution |
| `src/pmemd/src/runmd.F90` | KAMD integration into the MD driver loop |
| `src/pmemd/src/cuda/kKAMD.h` | KAMD CUDA kernels |
| `src/pmemd/src/cuda/gpu.cpp` | Host-side `gpu_kamd_setup_` and `gpu_kamd_pseudobd_reset_` |
| `test/kamd/md_baseline.in` | Solution-only baseline test |
| `test/kamd/md_membrane_parse.in` | Mixed solution + TM example |
| `test/kamd/md_overlap_*.in`, `md_tm_*.in` | Validation trip-wires (each fails by design) |
| `test/kamd/md_iontraj.in`, `Run_Iontraj.sh`, `check_kamd_ions.py` | Ion-trajectory pipeline test |

---

## 11. Ion-trajectory output (`ion_traj_outfreq`)

Set `ion_traj_outfreq` (MD steps) to record the 3D positions of the atoms
selected by `ion_mask` every N steps, plus the per-frame membrane center-of-mass
(from `membrane_mask`) and the orthorhombic box. This is **independent of
`pseudobd_cycles`** — it works in a plain MD run, with negligible overhead. The
output is plain text; derive ionic current, permeation events, and occupancy
offline (see below).

```
&kamd
   ion_traj_outfreq = 5000,        ! 10 ps at dt=0.002; 0 disables
   ion_mask     = ":Na+,K+,Cl-",
   membrane_mask = ":POPC,POPE",   ! contiguous range; the COM reference
   ion_traj_file = "kamd_ions.dat",
/
```

`membrane_mask` must resolve to a contiguous atom-index range (it is iterated
as `[start, end]` for the COM reduction). `ion_mask` may be discontiguous.

### File format

```
# kamd ion trajectory
# membrane_axis = 3   dt = 0.00200 ps   ion_traj_outfreq = 5000 steps
# per frame: 1 FRAME line (step, time_ps, memCOM_x memCOM_y memCOM_z, Lx Ly Lz)
#          + N ion lines (atom_index atom_name x y z , Angstrom, wrapped)
# ion roster: index atom_name residue charge(e)
#   34201 NA  Na+  +1.0000
#   34202 CL  Cl-  -1.0000
FRAME 5000 10.000   1.234 -0.567 45.100   88.0 88.0 92.0
34201  NA   12.30  44.5  45.0
34202  CL   50.10   3.2  47.8
```

Each frame: one `FRAME` line (`step time_ps memCOM_x memCOM_y memCOM_z Lx Ly Lz`),
then one line per ion in roster order (`atom_index atom_name x y z`, Å, wrapped).
The file is flushed after every frame, so a walltime-killed job leaves a complete
file up to the last recorded frame. Each `pmemd` run starts a fresh file;
concatenate segments offline for a continuous series.

### Offline analysis (current, crossings, occupancy)

With the membrane COM and box in every frame, everything is derivable from the
ion coordinates alone:

- **Axial position relative to the membrane:** `z_rel = z_ion - memCOM_z`.
- **Pore crossing:** a sign change of `z_rel` (robust to box drift and `ntp`
  box-size changes).
- **Ionic current (charge-displacement):** `I = (1 / (Δt · L_z)) · Σ_i q_i · Δz_i`,
  block-averaged for the mean and for an I–V curve across voltages. `+efz`
  drives cations toward `+z`, which sets the current sign.
- **PBC unwrap:** use `L_z` (and `Lx, Ly` laterally) to reconstruct continuous
  trajectories.
