# KAMD toy model: does KAMD sample the canonical configurational distribution? (2026-11-03)

## Conclusion

**KAMD biases the configurational distribution, and the bias is detectable even at the
default settings.** Plain Langevin MD and an "MB resample only" control agree with each
other within error. Every condition that includes the speed-preserving redirect of the
shell and surface atoms differs from the reference. The size of the bias grows with the reset rate:

| condition (γ = 1 ps⁻¹) | resets / ns | PMF barrier, outer side (kT) | Rg variance / ref | solute–solvent g(r) 1st peak, Δ | P(r) χ² p | transitions / ns |
|---|---|---|---|---|---|---|
| plain Langevin (reference) | – | **2.00 ± 0.20** | 1 | – | – | 1.72 ± 0.21 |
| MB resample only, trigger 20/1 | 2316 | 1.68 ± 0.19 | **1.001 ± 0.006** | +0.013 (χ² p = 1.0) | **0.99** | 0.62 ± 0.21 |
| **KAMD defaults 250/2** | 159 | **1.35 ± 0.15** | **0.943 ± 0.009** | **+0.034 (z = 8)** | **0.001** | 3.65 ± 0.23 (2.1×) |
| KAMD frequent 20/1 | 1891 | **0.26 ± 0.08** | **0.593 ± 0.004** | **+0.36 (z = 64)** | **10⁻⁴⁷** | 12.2 ± 0.8 (7.1×) |
| KAMD 20/1, separation-vector monitor | 1325 | 0.49 ± 0.16 | 0.663 ± 0.010 | +0.27 (z = 36) | 10⁻³¹ | 11.4 ± 1.1 (6.6×) |

* **Precision.** The bias test has these 1σ resolutions:
  * PMF barrier: ≈ 0.2–0.25 kT;
  * cluster Rg variance: ≈ 1 %;
  * solute–solvent g(r): ≈ 0.005 per 0.05 Å bin.
* **Default KAMD bias** (250/2, one reset every 6.3 ps):
  * the barrier seen from the outer well drops by **0.65 ± 0.25 kT** (2.6σ);
  * the cluster's internal fluctuations shrink: the Rg variance drops by **5.7 ± 0.9 %**, which corresponds to a cluster internal temperature of ≈ 113 K instead of 120 K;
  * the first solute–solvent g(r) peak rises by 0.034 (≈ 3 %, 8σ).
* **Frequent resets** (20/1): the barrier almost disappears (2.0 → 0.3 kT) and the PMF becomes flat between the wells.
* **ΔF between the wells does not detect the bias.** ΔF = F(inner) − F(outer) agrees with the reference in every condition. Two things cause this:
  * the reference ΔF error is large, ±0.38 kT (8 × 2.4 ns, about 4 transitions per replica);
  * both wells flatten roughly symmetrically, so their population ratio barely changes.

  On this observable we can only say "no detectable bias at ±0.4 kT". ΔF alone would have
  passed KAMD. The shape of P(r) and the internal and solvation structure do not.

## Mechanism (what the data show)

The redirect keeps every atom's speed. It still moves kinetic energy between modes: the
velocities of the 12 surface atoms of a cluster and of ~40 shell solvent atoms are
made parallel. Measured at the reset (figure `velocities.png`):

* **cluster COM kinetic energy: 1.5 kT → 12–14 kT** in the first 0.25 ps after a reset. It decays back to 1.5 kT in ≈ 1–2 ps at γ = 1 ps⁻¹. With the 20/1 trigger the time-averaged COM KE is 8.5 kT.
* **internal (relative) motion of the cluster is drained.** The solute kinetic temperature drops to **103 K** at 0.25–1 ps after a reset (20/1) and to **113.5 K** at 1–3 ps (250/2), then recovers.
* the whole-system kinetic temperature stays within 119.7–120.3 K. The solvent velocity distribution ≥ 1 ps after a reset matches Maxwell-Boltzmann to < 0.5 % for |u| < 3, with variance 1.000 and excess kurtosis 0.000.

So the system is not in equilibrium between resets. The inter-group coordinate is
**hot**, which flattens P(r), lowers the barrier and spreads the walls. The cluster's internal modes are
**cold**, which narrows the Rg distribution. The solute–solvent packing changes.
Langevin friction relaxes this, but the configurations sampled during relaxation are
part of the ensemble. A reset maps (x, p) to a momentum distribution that is not the
canonical conditional p | x. It is not even isotropic, and it depends on x through the shell
assignment. Nothing in the algorithm restores detailed balance, so there is no reason the stationary
distribution should be canonical. The bias scales with (relaxation time) / (reset interval):

| friction | condition | Rg variance / ref | barrier, outer side (kT) |
|---|---|---|---|
| γ = 0.1 ps⁻¹ | KAMD 20/1 | **0.27 ± 0.01** | 0.26 ± 0.04 (plain at γ = 0.1: 2.07 ± 0.24) |
| γ = 1 ps⁻¹ | KAMD 20/1 | 0.59 ± 0.004 | 0.26 ± 0.08 (plain: 2.00 ± 0.20) |
| γ = 10 ps⁻¹ | KAMD 20/1 | 0.89 ± 0.01 | –* |
| γ = 1 ps⁻¹ | KAMD 250/2 | 0.94 ± 0.01 | 1.35 ± 0.15 |

\* At γ = 10 ps⁻¹ plain MD made **0 transitions** in 4 × 1.1 ns, so P(r) and ΔF have no
reference. The Rg and g(r) biases are still significant: solute–solvent g(r) z = 21, P(r) χ² p = 0.002 against the unconverged reference.

**The trigger is not the problem.** The MB-resample-only control uses the same 20/1 trigger and fires
every 0.43 ps. It matches the reference in every observable: P(r) χ² p = 0.99, Rg
variance ratio 1.001 ± 0.006, both solute–solvent g(r) χ² p > 0.9. This also shows that the χ²
tests are not over-sensitive. Switching the monitor from the COM-velocity projection to the
inter-group separation vector (an experimental copy, `kamd_variants.py`) gives the same bias as the
default monitor, because the redirect is unchanged.

**Acceleration** (well-to-well transitions per ns, hysteresis cores r < 13.5 Å and r > 16.5 Å):
* 2.1× at the defaults;
* 7.1× at 20/1;
* 8.1× at 20/1 with γ = 0.1.

MB-only resampling every 0.4 ps is a strong Andersen-type thermostat and slows transitions to
0.36×. The acceleration therefore comes from the collective push, and the same push causes the bias.

## System and protocol

Code: `toy_system.py`, `run_toy.py`, `kamd_variants.py`, `make_jobs.sh`, `analyze.py`, `check_dt.py`.
`kamd_openmm/core.py` and `reporter.py` were not modified.

* **Solvent.** 2022 LJ particles: σ = 3.4 Å, ε = 0.996 kJ/mol, m = 39.9, T = 120 K (T* = 1.0), ρ* = 0.8. Cubic periodic box L = 46.5 Å, cutoff 2.5σ = 8.5 Å, long-range dispersion correction.
* **Solutes.** Two 13-particle icosahedral clusters (R = 3.7 Å), each held by 42 harmonic bonds (12 spokes + 30 edges, k = 20 000 kJ mol⁻¹ nm⁻², i.e. 200 kJ mol⁻¹ Å⁻²). All intra-cluster LJ pairs are excluded. The clusters interact with the solvent and with each other by the same LJ.
  * Construction: FCC lattice with 2048 sites; two 13-site cuboctahedra are swapped for icosahedra, so N = 2048 at exactly ρ* = 0.8.
* **Double well.** `CustomCentroidBondForce` on the periodic COM–COM distance: U = h[((r − r0)/a)² − 1]², with h = 2 kT, r0 = 15 Å, a = 2.5 Å (bare wells at 12.5 and 17.5 Å). With solvation, the reference PMF has its barrier at 15.1 Å and barrier heights of 2.0–2.3 kT.
* **Integrator.** `LangevinMiddleIntegrator`, γ = 1 ps⁻¹ unless stated, **dt = 5 fs**. Stability (`check_dt.out`): NVE drift is 4×10⁻⁴ (4 fs) and 3×10⁻² (5 fs) kT per particle per ns. The 2 fs run drifted 2×10⁻² kT, so the differences are single-precision noise. The γ = 1 ps⁻¹ thermostat dominates by orders of magnitude. A Langevin 4 fs vs 5 fs pilot gave the same Rg and T.
* **KAMD.** `refer_groups` = the two clusters, `solvent_atoms` = all LJ solvent.
  * `shell_cut` 5.3 Å: the measured solvent RDF first minimum is 5.33 Å.
  * `surface_cut` 4.0 Å.
  * T = 120 K.
  * Per reset: on average 40–43 shell atoms per group, and all 12 vertex atoms are surface atoms (the centre atom is buried).
  * MB-only control: `shell_cut` = `surface_cut` = 0.001 Å, which gives an empty shell.
* **Replicas.** 8 seeds × 2.5 ns for the four main conditions; 4 seeds × 1.2 ns for the γ scan and the separation-vector monitor. Each replica:
  * minimises, then runs 100 ps of plain Langevin (γ = 1 ps⁻¹) to equilibrate;
  * starts in the outer well for odd seeds and the inner well for even seeds;
  * discards the first 100 ps of production.

  Seed *s* uses the same starting structure in every condition (paired design).
* **Errors.** All error bars come from replica-to-replica scatter: SEM, or a bootstrap over replicas for ΔF and the PMF.
* **Compute.** Benchmark: 1 CPU thread, ≈ 470–490 steps/s (≈ 205 ns/day) per replica, 4 replicas in parallel on 4 cores. KAMD overhead is ≤ 4 % at 20/1. Production took 3.3 h wall on 4 cores (≈ 13 CPU-hours) for 52 runs and 104 ns of production in total. This is more than "a few CPU-hours", because the reference needs ≥ 2 ns per replica to make several transitions.

**Statistical tests** (all in `analyze.py`, numbers in `results.json` and `analysis_table.tsv`):

* **P(r):**
  * χ² over 1 Å bins, using replica variances;
  * two-sample KS on samples spaced 2τ apart (τ_ref = 143 ps for the well indicator; only 48–72 samples per condition, so KS has little power: p = 0.84 / 0.09 for 250/2 and 20/1).
* **PMF barrier:** F(r_b ± 0.5 Å) − min F in each well, bootstrap over replicas.
* **Rg:** KS test on samples every 5 ps, plus the ratio of per-replica variances.
* **g(r):** χ² over all bins, with replica variances.
* **Block averaging** (`block.png`): within-run block SEMs underestimate the replica SEM (0.04 vs 0.09 for p_inner of the reference), because the well-residence time is comparable to the block length. Replica SEMs are used everywhere.

## Further results

* **Solvent–solvent g(r).** No detectable change at the defaults: max |Δg| = 0.005, χ² p = 0.18. At 20/1 there is a small but significant change: max |Δg| = 0.010, < 0.4 % of the first peak. The solvent bulk is barely affected; the bias is localized at the groups.
* **Cluster COM–solvent g(r) first peak:**
  * +0.05 (z = 5) at the defaults;
  * +0.47 (z = 45) at 20/1;
  * +0.02 (n.s.) for MB-only.
* **Mean Rg** shifts only by −3×10⁻⁴ Å at 20/1. The Rg *width* is the sensitive observable: the bonds are stiff, so the mean barely moves.

## Limitations

* This is a toy: 13-atom groups, a ~40-atom shell each, and simple LJ liquid. The groups are small compared with their shell, so the collective push is large compared with the groups' own thermal motion. In proteins, each group's COM is a much smaller fraction of its degrees of freedom. The magnitude of the bias may be very different there. The mechanism (energy moved from internal and relative modes into collective translation, followed by finite-time relaxation) is not specific to the toy.
* The ΔF comparison is weak (±0.4 kT). Plain MD makes only ~1.7 transitions/ns, and the γ = 10 ps⁻¹ reference made none.
* The γ scan and the separation monitor have only 4 × 1.1 ns each.
* The "away from reset" velocity analysis for the 20/1 conditions has few frames ≥ 1 ps after a reset, because resets come every ≈ 0.5 ps.
* The kinetic temperatures use the synchronous estimate v(t − dt/2) + F dt / 2m, the same one KAMD uses.
* OpenMM CPU platform (single-precision forces); 5 fs time step. The pilot dt check showed no effect on Rg or T.
* Paired starting structures make conditions correlated at the start. The 100 ps burn-in plus 100 ps discard is shorter than τ for well residence, so well populations can carry a start memory. This is why half the seeds start in each well.

## Reproduce

```bash
pip install openmm scipy numpy matplotlib && pip install -e ./openmm
cd openmm/validation/toy_model
python check_dt.py > check_dt.out
./make_jobs.sh && cat jobs.txt | xargs -P 4 -I{} bash -c "{}"   # ~3.3 h on 4 cores
python analyze.py > analysis_table.tsv                           # figures/, results.json
```

The compact per-run outputs (`runs/*.npz`, 16 MB) are committed, so `analyze.py` can be re-run without simulating.

## Figures

* `figures/pmf_gamma1.png`: P(r), PMF and PMF difference vs reference (γ = 1 ps⁻¹); `pmf_gamma0.1.png` and `pmf_gamma10.png` for the friction scan.
* `figures/dF_rates.png`: ΔF between the wells and transition rates, all conditions.
* `figures/gr.png`: solvent–solvent, solute atom–solvent and cluster COM–solvent g(r), and differences (±2 SE bands).
* `figures/rg.png`: cluster Rg distributions.
* `figures/velocities.png`: solute kinetic temperature and cluster COM kinetic energy vs time since reset; solvent velocity distribution vs Maxwell-Boltzmann.
* `figures/block.png`: block-averaging convergence.
