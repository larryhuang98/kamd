# kamd-regions: automatic vs. manual reference groups

System: amber/test (prmtop + restrt + md_test.in) — 29467 atoms, 39 polymer residues, 301 heavy atoms in the polymer.

Generated with

```bash
kamd-regions compare amber/test/prmtop amber/test/restrt amber/test/md_test.in --methods gnm,contacts
kamd-regions compare amber/test/prmtop amber/test/restrt amber/test/md_test.in --methods gnm,contacts \
    --param n_regions=4 --param min_heavy=40 --param resolution=1.0
```

**Summary.** The test protein is a 39-residue, mostly unstructured chain (DSSP: one short helix,
residues 9-13). The manual groups are four short motifs of 44-117 atoms each, chosen by hand to
exercise the code, not domains. The automatic methods are designed for domains of a few hundred
heavy atoms or more, so on this system they do not reproduce the manual choice:

* With the defaults, GNM and the contact method both return two halves (GNM: A:9-18, A:19-39;
  contacts: A:9-19, A:20-39). Manual groups 2-4 fall inside these halves. Manual group 1
  (A:2-5) is not covered at all: the post-processing marks the coil N-terminus (residues 1-8, high
  GNM fluctuation) as a flexible tail and leaves it unassigned.
* Asking for 4 regions with a 40-heavy-atom threshold gives three blocks (one candidate is merged
  as too small). Agreement improves (ARI 0.24 vs. 0.07 for GNM), but the blocks are contiguous
  stretches that fill the gaps the manual groups leave open.
* The automatic groups have 473 atoms in total against 293 in the manual groups: the manual
  setup leaves linkers between the motifs unassigned, while the methods assign every structured
  residue.

The test system therefore checks the import/export path exactly (see the integration tests) but is
not a meaningful benchmark for automatic domain detection. The VDAC test (3emn, `tm` method) is
a better check of the automatic path: all 19 β-strands are detected as transmembrane segments.
For multi-domain proteins with known domains see [gnm_multidomain.md](gnm_multidomain.md).

## Default parameters

Manual groups (imported):

| group | residues | atoms |
|---|---|---|
| group 1 | A:2-5 | 77 |
| group 2 | A:18-20 | 55 |
| group 3 | A:23-30 | 117 |
| group 4 | A:36-38 | 44 |

Sequence view (one character per residue; `.` = unassigned):

```
   manual  .1111............222..33333333.....444.
      gnm  ........1111111111222222222222222222222
 contacts  ........1111111111122222222222222222222
```

### gnm

2 groups. Adjusted Rand index vs. manual: 0.07 over all polymer residues (unassigned counted as a class), 0.18 over residues assigned in either. Atoms in groups: manual 293, gnm 473, in both 216.

| auto group | type | residues | atoms |
|---|---|---|---|
| region 1 | solution | A:9-18 | 154 |
| region 2 | solution | A:19-39 | 319 |

| manual group | atoms | best-overlapping auto group | manual atoms covered | Jaccard | atoms left unassigned |
|---|---|---|---|---|---|
| group 1 | 77 | (unassigned) | 0% | 0.00 | 77 |
| group 2 | 55 | region 2 | 65% | 0.11 | 0 |
| group 3 | 117 | region 2 | 100% | 0.37 | 0 |
| group 4 | 44 | region 2 | 100% | 0.14 | 0 |

Notes: The protein has only 233 heavy atoms in structured parts, so the minimum region size was lowered from 150 to 77 heavy atoms to keep at least 2-3 regions. 1 region(s) are below 77 heavy atoms but were kept so that at least 2 regions remain.

### contacts

2 groups. Adjusted Rand index vs. manual: 0.09 over all polymer residues (unassigned counted as a class), 0.20 over residues assigned in either. Atoms in groups: manual 293, contacts 473, in both 216.

| auto group | type | residues | atoms |
|---|---|---|---|
| region 1 | solution | A:9-19 | 176 |
| region 2 | solution | A:20-39 | 297 |

| manual group | atoms | best-overlapping auto group | manual atoms covered | Jaccard | atoms left unassigned |
|---|---|---|---|---|---|
| group 1 | 77 | (unassigned) | 0% | 0.00 | 77 |
| group 2 | 55 | region 1 | 75% | 0.22 | 0 |
| group 3 | 117 | region 2 | 100% | 0.39 | 0 |
| group 4 | 44 | region 2 | 100% | 0.15 | 0 |

Notes: The protein has only 233 heavy atoms in structured parts, so the minimum region size was lowered from 150 to 77 heavy atoms to keep at least 2-3 regions.

## n_regions=4, min_heavy=40, resolution=1.0

Sequence view (one character per residue; `.` = unassigned):

```
   manual  .1111............222..33333333.....444.
      gnm  ........1111111111222222222222333333333
 contacts  ........1111111222222223333333333333333
```

### gnm

3 groups. Adjusted Rand index vs. manual: 0.24 over all polymer residues (unassigned counted as a class), 0.38 over residues assigned in either. Atoms in groups: manual 293, gnm 473, in both 216.

| auto group | type | residues | atoms |
|---|---|---|---|
| region 1 | solution | A:9-18 | 154 |
| region 2 | solution | A:19-30 | 179 |
| region 3 | solution | A:31-39 | 140 |

| manual group | atoms | best-overlapping auto group | manual atoms covered | Jaccard | atoms left unassigned |
|---|---|---|---|---|---|
| group 1 | 77 | (unassigned) | 0% | 0.00 | 77 |
| group 2 | 55 | region 2 | 65% | 0.18 | 0 |
| group 3 | 117 | region 2 | 100% | 0.65 | 0 |
| group 4 | 44 | region 3 | 100% | 0.31 | 0 |

### contacts

3 groups. Adjusted Rand index vs. manual: 0.09 over all polymer residues (unassigned counted as a class), 0.19 over residues assigned in either. Atoms in groups: manual 293, contacts 473, in both 216.

| auto group | type | residues | atoms |
|---|---|---|---|
| region 1 | solution | A:9-15 | 113 |
| region 2 | solution | A:16-23 | 115 |
| region 3 | solution | A:24-39 | 245 |

| manual group | atoms | best-overlapping auto group | manual atoms covered | Jaccard | atoms left unassigned |
|---|---|---|---|---|---|
| group 1 | 77 | (unassigned) | 0% | 0.00 | 77 |
| group 2 | 55 | region 2 | 100% | 0.48 | 0 |
| group 3 | 117 | region 3 | 90% | 0.41 | 0 |
| group 4 | 44 | region 3 | 100% | 0.18 | 0 |
