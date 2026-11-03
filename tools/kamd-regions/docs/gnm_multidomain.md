# GNM on multi-domain proteins

The `gnm` method was tested on proteins whose domains are known from the literature:

| protein | structure | literature domains |
|---|---|---|
| adenylate kinase (E. coli), open | 4AKE (`tests/data/adk_open.pdb`, from the MDAnalysis test suite) | CORE 1–29, 60–121, 160–214; NMP 30–59; LID 122–159 |
| adenylate kinase, closed | 1AKE (MDAnalysis test suite) | same |
| T4 lysozyme L99A | `tests/data/t4l.pdb` (openmmtools) | N-domain 13–59; C-domain 1–12, 60–164 |
| HhaI DNA methyltransferase + DNA | 3MHT (ProDy test suite) | large 1–193, 304–327; small 194–303 |

Agreement is the adjusted Rand index (ARI) between the suggested regions and the literature domains,
over residues that are assigned in both (1 = identical partition, 0 = chance). Run
`python docs/benchmark_gnm.py` to reproduce the tables below.

## Summary

| protein | defaults | k = number of domains |
|---|---|---|
| ADK open | 0.37 (k=2: CORE vs. LID) | **0.71** (CORE / NMP / LID) |
| ADK closed | 0.47 (k=2: CORE vs. LID) | 0.36 |
| T4 lysozyme | **0.64** | 0.64 |
| 3MHT | **0.83** | 0.83 |

* **T4 lysozyme** and **3MHT**: the eigengap picks the right number of domains and the partition
  matches the literature closely. In T4 lysozyme the boundary falls at residue 71/72 instead of
  59/60: the long helix 60–80 that connects the domains is split, and GNM puts its first half with
  the N-domain. In 3MHT the loops around 172–191 and 283–295 are left unassigned as flexible.
* **ADK open**: with k = 3 the three domains are found (NMP = 31–69, LID = 118–164 vs. 30–59 and
  122–159 in the literature). The eigengap alone chooses k = 2 (CORE vs. LID, the largest motion),
  so set `n_regions = 3` (GUI: ⚙ next to Re-suggest).
* **ADK closed**: the LID is found, but in the compact closed state the NMP domain is fused with
  part of the CORE. A single closed structure does not show the NMP hinge clearly; this is a known
  limit of single-structure network models. Use the open structure, or the DynDom-style two-state
  analysis once it exists (roadmap).

## What changed after this test

The first implementation clustered row-normalised Kirchhoff eigenvectors with a 7.3 Å cutoff, over
all polymer residues (protein + nucleic acid). On these proteins it was the worst variant tried:

| variant (k = number of domains) | ADK open | ADK closed | T4L | 3MHT | mean |
|---|---|---|---|---|---|
| Kirchhoff modes, row-normalised, 7.3 Å, protein + DNA (first version) | 0.26 | 0.21 | 0.54 | 0.07 | 0.27 |
| Kirchhoff modes, not normalised, 10 Å, protein + DNA | 0.74 | 0.34 | 0.68 | 0.00 | 0.44 |
| normalised Laplacian, 10 Å | 0.42 | 0.30 | 0.62 | 0.72 | 0.52 |
| GNM cross-correlation (20 modes), k-means on rows, 10 Å | 0.38 | 0.34 | 0.68 | 0.73 | 0.53 |
| ANM cross-correlation (15 Å, 20 modes) | 0.20 | 0.15 | 0.08 | 0.05 | 0.12 |
| **Kirchhoff modes, not normalised, 10 Å, protein only (now)** | **0.74** | 0.34 | **0.68** | **0.77** | **0.63** |

(raw clustering, before post-processing; the full pipeline numbers are in the summary above.)

Four changes followed:

1. **No row normalisation.** Residues near a hinge have small slow-mode amplitudes; normalising
   them onto the unit sphere assigns them almost at random.
2. **Network of protein residues only.** In 3MHT the two short DNA strands are loosely attached to
   the protein, and the slowest modes describe the DNA moving against the protein, so the protein
   domains were not separated (ARI 0.00). Nucleic acids are now left out of the GNM network (and
   left unassigned by `gnm`); if there is no protein they are used.
3. **Cutoff 10 Å** (was 7.3 Å, the ProDy default), slightly better on ADK and T4L.
4. **Minimum region size 150 heavy atoms** (was 300). The ADK NMP domain has 212 heavy atoms
   (~420 atoms with hydrogens), a perfectly sized KAMD group, and was merged away at 300.

Two regression tests (`tests/test_methods.py::test_gnm_recovers_domains`) check ADK open (k = 3,
ARI ≥ 0.6, each domain the majority of exactly one region) and T4 lysozyme (defaults, ARI ≥ 0.55).

## Full results

### ADK open (4AKE)

**defaults**: 2 regions (k = 2), ARI 0.37, 10 residues unassigned.

| region | residues | atoms | best domain | precision | recall |
|---|---|---|---|---|---|
| 1 | 4AKE:1-113, 4AKE:167-214 | 2480 | CORE | 81% | 90% |
| 2 | 4AKE:114-134, 4AKE:145-166 | 701 | LID | 65% | 74% |

**n_regions=3**: 3 regions (k = 3), ARI 0.71, 10 residues unassigned.

| region | residues | atoms | best domain | precision | recall |
|---|---|---|---|---|---|
| 1 | 4AKE:1-30, 4AKE:70-117, 4AKE:164-214 | 2023 | CORE | 99% | 88% |
| 2 | 4AKE:31-69 | 588 | NMP | 74% | 97% |
| 3 | 4AKE:118-134, 4AKE:145-163 | 570 | LID | 78% | 74% |

```
   reference  1111111111111111111111111111122222222222222222222222222222211111111111111111111111111111111111111111111111111111111111111333333333333333333333333333333333333331111111111111111111111111111111111111111111111111111111
    defaults  11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111222222222222222222222..........2222222222222222222222111111111111111111111111111111111111111111111111
 n_regions=3  11111111111111111111111111111122222222222222222222222222222222222222211111111111111111111111111111111111111111111111133333333333333333..........3333333333333333333111111111111111111111111111111111111111111111111111
```

### ADK closed (1AKE)

**defaults**: 2 regions (k = 2), ARI 0.47, 0 residues unassigned.

| region | residues | atoms | best domain | precision | recall |
|---|---|---|---|---|---|
| 1 | 4AKE:1-115, 4AKE:165-214 | 2560 | CORE | 82% | 92% |
| 2 | 4AKE:116-164 | 781 | LID | 78% | 100% |

**n_regions=3**: 3 regions (k = 3), ARI 0.36, 0 residues unassigned.

| region | residues | atoms | best domain | precision | recall |
|---|---|---|---|---|---|
| 1 | 4AKE:1-28, 4AKE:75-85, 4AKE:102-113, 4AKE:177-214 | 1363 | CORE | 100% | 61% |
| 2 | 4AKE:29-74, 4AKE:86-101, 4AKE:162-176 | 1208 | CORE | 61% | 32% |
| 3 | 4AKE:114-161 | 770 | LID | 79% | 100% |

```
   reference  1111111111111111111111111111122222222222222222222222222222211111111111111111111111111111111111111111111111111111111111111333333333333333333333333333333333333331111111111111111111111111111111111111111111111111111111
    defaults  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111222222222222222222222222222222222222222222222222211111111111111111111111111111111111111111111111111
 n_regions=3  1111111111111111111111111111222222222222222222222222222222222222222222222211111111111222222222222222211111111111133333333333333333333333333333333333333333333333322222222222222211111111111111111111111111111111111111
```

### T4 lysozyme L99A

**defaults**: 2 regions (k = 2), ARI 0.64, 9 residues unassigned.

| region | residues | atoms | best domain | precision | recall |
|---|---|---|---|---|---|
| 1 | A:3-8, A:71-155 | 1487 | C-domain | 100% | 79% |
| 2 | A:9-70 | 972 | N-domain | 76% | 100% |

**n_regions=2**: 2 regions (k = 2), ARI 0.64, 9 residues unassigned.

| region | residues | atoms | best domain | precision | recall |
|---|---|---|---|---|---|
| 1 | A:3-8, A:71-155 | 1487 | C-domain | 100% | 79% |
| 2 | A:9-70 | 972 | N-domain | 76% | 100% |

```
   reference  111111111111222222222222222222222222222222222222222222222221111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
    defaults  ..111111222222222222222222222222222222222222222222222222222222222222221111111111111111111111111111111111111111111111111111111111111111111111111111111111111.......
 n_regions=2  ..111111222222222222222222222222222222222222222222222222222222222222221111111111111111111111111111111111111111111111111111111111111111111111111111111111111.......
```

### HhaI methyltransferase + DNA (3MHT)

**defaults**: 2 regions (k = 2), ARI 0.83, 49 residues unassigned.

| region | residues | atoms | best domain | precision | recall |
|---|---|---|---|---|---|
| 1 | A:1-153, A:164-172, A:304-323 | 1455 | large | 100% | 84% |
| 2 | A:154-163, A:192-214, A:228-282, A:295-303 | 747 | small | 88% | 77% |

**n_regions=2**: 2 regions (k = 2), ARI 0.83, 49 residues unassigned.

| region | residues | atoms | best domain | precision | recall |
|---|---|---|---|---|---|
| 1 | A:1-153, A:164-172, A:304-323 | 1455 | large | 100% | 84% |
| 2 | A:154-163, A:192-214, A:228-282, A:295-303 | 747 | small | 88% | 77% |

```
   reference  111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111122222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222111111111111111111111111.
    defaults  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111112222222222111111111...................22222222222222222222222.............2222222222222222222222222222222222222222222222222222222............22222222211111111111111111111.....
 n_regions=2  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111112222222222111111111...................22222222222222222222222.............2222222222222222222222222222222222222222222222222222222............22222222211111111111111111111.....
```
