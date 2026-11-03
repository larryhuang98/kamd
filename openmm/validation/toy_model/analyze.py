"""Analysis of the KAMD toy-model runs -> figures/*.png + results.json (+ printed tables).

    python analyze.py [runs_dir]
"""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
RUNS = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "runs"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

BURN_PS = 100.0  # discarded from each production run
R_EDGES = np.arange(9.0, 21.01, 0.25)  # Angstrom
R_BARRIER = 15.0  # dividing surface between the wells (refined from the reference PMF below)
CORE_IN, CORE_OUT = 13.5, 16.5  # hysteresis cores for counting transitions
COARSE = np.arange(9.0, 21.01, 1.0)
RNG = np.random.default_rng(12345)

ORDER = ["plain", "kamd_def", "kamd_freq", "mbonly", "sepmon_freq",
         "plain_g0.1", "kamd_freq_g0.1", "plain_g10", "kamd_freq_g10"]
LABEL = {
    "plain": "plain Langevin (ref)",
    "kamd_def": "KAMD 250/2",
    "kamd_freq": "KAMD 20/1",
    "mbonly": "MB resample only 20/1",
    "sepmon_freq": "KAMD 20/1, separation monitor",
    "plain_g0.1": "plain, γ=0.1",
    "kamd_freq_g0.1": "KAMD 20/1, γ=0.1",
    "plain_g10": "plain, γ=10",
    "kamd_freq_g10": "KAMD 20/1, γ=10",
}
REF_OF = {c: ("plain_g0.1" if "g0.1" in c else "plain_g10" if "g10" in c else "plain") for c in ORDER}
COL = {"plain": "#3a3a38", "kamd_def": "#2a78d6", "kamd_freq": "#eb6834", "mbonly": "#1baf7a",
       "sepmon_freq": "#eda100", "plain_g0.1": "#3a3a38", "kamd_freq_g0.1": "#eb6834",
       "plain_g10": "#3a3a38", "kamd_freq_g10": "#eb6834"}


def load():
    runs = defaultdict(list)
    for f in sorted(RUNS.glob("*.npz")):
        m = re.match(r"(.+)_s(\d+)$", f.stem)
        if not m:
            continue
        d = dict(np.load(f, allow_pickle=False))
        d["meta"] = json.loads(str(d["meta"]))
        d["seed"] = int(m.group(2))
        dt = d["meta"]["dt_fs"] * 1e-3
        t = d["step"] * dt
        keep = t >= BURN_PS
        d["t"] = t[keep]
        d["rA"] = d["r"][keep].astype(float) * 10.0
        d["rgA"] = d["rg"][keep].astype(float) * 10.0
        d["dt_sample"] = float(np.median(np.diff(t)))
        runs[m.group(1)].append(d)
    return {c: runs[c] for c in ORDER if c in runs}


# ----------------------------------------------------------------------------- helpers
def pmf_from_p(p, centers):
    with np.errstate(divide="ignore"):
        return -np.log(p / centers**2)


def boot_reps(vals_fn, reps, n=2000):
    """Bootstrap over replicas: vals_fn(list_of_reps) -> scalar/array."""
    out = []
    k = len(reps)
    for _ in range(n):
        idx = RNG.integers(0, k, k)
        out.append(vals_fn([reps[i] for i in idx]))
    return np.array(out)


def p_inner(reps, rb):
    return np.mean([np.mean(d["rA"] < rb) for d in reps])


def dF(p):
    return -np.log(p / (1 - p))  # F_inner - F_outer in kT


def count_transitions(r):
    s = np.where(r < CORE_IN, -1, np.where(r > CORE_OUT, 1, 0))
    s = s[s != 0]
    return int(np.sum(s[1:] != s[:-1]))


def acf_time(x, dt, maxlag_frac=0.5):
    """Integrated autocorrelation time (Sokal window, c=5)."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    if not np.any(x):
        return np.nan  # no transitions in this replica
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac /= ac[0]
    tau = 1.0
    for m in range(1, int(n * maxlag_frac)):
        tau = 1 + 2 * ac[1:m + 1].sum()
        if m >= 5 * tau:
            break
    return max(tau, 1.0) * dt


def block_sem(x, nblocks):
    n = len(x) // nblocks
    b = x[: n * nblocks].reshape(nblocks, n).mean(1)
    return b.std(ddof=1) / np.sqrt(nblocks)


# ----------------------------------------------------------------------------- main
def main():
    R = load()
    res = {"burn_ps": BURN_PS, "conditions": {}}
    centers = 0.5 * (R_EDGES[1:] + R_EDGES[:-1])

    # per replica histograms
    H = {}
    for c, reps in R.items():
        h = np.array([np.histogram(d["rA"], bins=R_EDGES, density=True)[0] for d in reps])
        H[c] = h
    # refine dividing surface: PMF maximum of the reference between the wells
    ref = H["plain"].mean(0)
    F = pmf_from_p(ref, centers)
    win = (centers > 13.0) & (centers < 17.0)
    rb = float(centers[win][np.argmax(F[win])])
    res["r_barrier"] = rb

    # decorrelation time of the reference (well indicator) -> subsampling for KS
    tau_ref = np.nanmean([acf_time((d["rA"] < rb).astype(float), d["dt_sample"]) for d in R["plain"]])
    res["tau_indicator_ref_ps"] = tau_ref

    rows = []
    for c, reps in R.items():
        meta = reps[0]["meta"]
        T_ns = np.array([(d["t"][-1] - d["t"][0] + d["dt_sample"]) / 1000 for d in reps])
        pin = np.array([np.mean(d["rA"] < rb) for d in reps])
        p = pin.mean()
        bs = boot_reps(lambda rr: p_inner(rr, rb), reps)
        dFc = dF(p)
        dF_err = np.std(dF(np.clip(bs, 1e-6, 1 - 1e-6)))
        ntr = np.array([count_transitions(d["rA"]) for d in reps])
        rate = ntr / T_ns
        tau_c = np.nanmean([acf_time((d["rA"] < rb).astype(float), d["dt_sample"]) for d in reps])
        # within-run block analysis of p_inner (10 blocks per run), pooled
        blk = np.sqrt(np.sum([block_sem((d["rA"] < rb).astype(float), 5) ** 2 for d in reps])) / len(reps)
        rg = np.concatenate([d["rgA"].ravel() for d in reps])
        rg_rep = np.array([d["rgA"].mean() for d in reps])
        tl = np.concatenate([d["tlog"] for d in reps])
        uh = np.sum([d["uhist"] for d in reps], axis=0)
        n_res = np.array([d["meta"]["n_resets"] for d in reps])
        sps = np.mean([d["meta"]["steps_per_s"] for d in reps])
        e = dict(
            n_reps=len(reps), ns_per_rep=float(T_ns.mean()), friction=meta["friction"], mode=meta["mode"],
            pbd=meta["pbd"], nz=meta["nz"],
            p_inner=float(p), p_inner_sem=float(pin.std(ddof=1) / np.sqrt(len(pin))), p_inner_reps=pin.tolist(),
            p_inner_block_sem=float(blk),
            dF=float(dFc), dF_err=float(dF_err),
            trans_per_ns=float(rate.mean()), trans_per_ns_sem=float(rate.std(ddof=1) / np.sqrt(len(rate))),
            n_trans_total=int(ntr.sum()), tau_indicator_ps=float(tau_c),
            rg_mean=float(rg.mean()), rg_sem=float(rg_rep.std(ddof=1) / np.sqrt(len(rg_rep))),
            T_mean=float(tl[:, 2].mean()), T_half_mean=float(tl[:, 5].mean()),
            T_solute=float(tl[:, 4].mean()), KEcom_kT=float(tl[:, 6:8].mean()),
            resets_per_ns=float((n_res / (np.array([d["meta"]["ns"] for d in reps]))).mean()),
            steps_per_s=float(sps),
        )
        if n_res.sum():
            e["shell_per_group"] = float(np.mean([d["reset_shell"].mean() for d in reps if len(d["reset_shell"])]))
            e["surface_per_group"] = float(np.mean([d["reset_surface"].mean() for d in reps if len(d["reset_surface"])]))
            e["reset_interval_ps"] = float(1000.0 / e["resets_per_ns"])
        # time-since-reset classes
        edges = reps[0]["tsr_edges"]
        cls = np.minimum(np.searchsorted(edges, tl[:, 1], side="right") - 1, len(edges) - 2)
        e["tsr_classes"] = []
        for k in range(len(edges) - 1):
            m = cls == k
            if not m.any():
                continue
            hs = uh[k, 0]
            x = np.linspace(-6, 6, hs.size + 1)
            xc = 0.5 * (x[1:] + x[:-1])
            var = (hs * xc**2).sum() / hs.sum()
            kurt = (hs * xc**4).sum() / hs.sum() / var**2 - 3
            e["tsr_classes"].append(dict(
                lo_ps=float(edges[k]), hi_ps=float(edges[k + 1]), n=int(m.sum()),
                T=float(tl[m, 2].mean()), T_solv=float(tl[m, 3].mean()), T_solute=float(tl[m, 4].mean()),
                KEcom_kT=float(tl[m, 6:8].mean()), solvent_u_var=float(var), solvent_u_exkurt=float(kurt)))
        res["conditions"][c] = e
        rows.append(c)

    # ---------------------------------------------------------- comparisons vs reference
    for c in rows:
        e = res["conditions"][c]
        rc = REF_OF[c]
        if rc not in res["conditions"] or rc == c:
            continue
        er = res["conditions"][rc]
        e["ref"] = rc
        e["ddF"] = e["dF"] - er["dF"]
        e["ddF_err"] = float(np.hypot(e["dF_err"], er["dF_err"]))
        e["ddF_z"] = e["ddF"] / e["ddF_err"]
        e["welch_p_inner_p"] = float(stats.ttest_ind(e["p_inner_reps"], er["p_inner_reps"], equal_var=False).pvalue)
        # chi2 on coarse P(r) with replica variances
        def coarse(reps):
            return np.array([np.histogram(d["rA"], bins=COARSE)[0] / len(d["rA"]) for d in reps])
        A, B = coarse(R[c]), coarse(R[rc])
        va = A.var(0, ddof=1) / len(A)
        vb = B.var(0, ddof=1) / len(B)
        ok = (va + vb) > 0
        chi2 = float((((A.mean(0) - B.mean(0)) ** 2)[ok] / (va + vb)[ok]).sum())
        dfree = int(ok.sum() - 1)
        e["chi2_Pr"] = chi2
        e["chi2_Pr_df"] = dfree
        e["chi2_Pr_p"] = float(stats.chi2.sf(chi2, dfree))
        # KS on decorrelated samples (stride = 2 * max(tau_c, tau_ref))
        tau = max(e["tau_indicator_ps"], er["tau_indicator_ps"])
        def sub(reps):
            st = max(1, int(round(2 * tau / reps[0]["dt_sample"])))
            return np.concatenate([d["rA"][::st] for d in reps])
        sa, sb = sub(R[c]), sub(R[rc])
        ks = stats.ks_2samp(sa, sb)
        e["ks_D"], e["ks_p"], e["ks_n"] = float(ks.statistic), float(ks.pvalue), [len(sa), len(sb)]
        # Rg: KS on samples every 5 ps (Rg decorrelates in << 1 ps)
        st = int(round(5.0 / R[c][0]["dt_sample"]))
        ga = np.concatenate([d["rgA"][::st].ravel() for d in R[c]])
        gb = np.concatenate([d["rgA"][::st].ravel() for d in R[rc]])
        k2 = stats.ks_2samp(ga, gb)
        e["rg_ks_D"], e["rg_ks_p"] = float(k2.statistic), float(k2.pvalue)
        e["rate_ratio"] = e["trans_per_ns"] / er["trans_per_ns"] if er["trans_per_ns"] > 0 else None

    # ---------------------------------------------------------- g(r)
    gres = {}
    for c, reps in R.items():
        L = reps[0]["meta"]["L"] * 10
        Ns = reps[0]["meta"]["N"] - 26
        rho = Ns / L**3
        e_ss = np.linspace(0, 12, reps[0]["gr_ss"].size + 1)
        e_us = np.linspace(0, 15, reps[0]["gr_us"].shape[1] + 1)
        sv = 4 / 3 * np.pi * (e_ss[1:] ** 3 - e_ss[:-1] ** 3)
        su = 4 / 3 * np.pi * (e_us[1:] ** 3 - e_us[:-1] ** 3)
        gss = np.array([d["gr_ss"] / (d["n_gr"] * Ns * rho * sv) for d in reps])
        gua = np.array([d["gr_us"][0] / (d["n_sg"] * 26 * rho * su) for d in reps])
        guc = np.array([d["gr_us"][1] / (d["n_sg"] * 2 * rho * su) for d in reps])
        gres[c] = dict(rss=0.5 * (e_ss[1:] + e_ss[:-1]), rus=0.5 * (e_us[1:] + e_us[:-1]), gss=gss, gua=gua, guc=guc)
    for c in rows:
        rc = REF_OF[c]
        if rc == c or rc not in gres:
            continue
        e = res["conditions"][c]
        for key, name in (("gss", "g_ss"), ("gua", "g_solute_atom_solvent"), ("guc", "g_com_solvent")):
            a, b = gres[c][key], gres[rc][key]
            d = a.mean(0) - b.mean(0)
            se = np.sqrt(a.var(0, ddof=1) / len(a) + b.var(0, ddof=1) / len(b))
            m = se > 0
            e[name + "_max_abs_diff"] = float(np.abs(d).max())
            e[name + "_max_abs_z"] = float(np.max(np.abs(d[m] / se[m])))
            chi2 = float(((d[m] / se[m]) ** 2).sum())
            e[name + "_chi2_p"] = float(stats.chi2.sf(chi2, int(m.sum())))
    # rdf minimum of the reference (for shell_cut)
    g = gres["plain"]["gss"].mean(0)
    rs = gres["plain"]["rss"]
    w = (rs > 4.0) & (rs < 7.0)
    res["solvent_rdf_first_min_A"] = float(rs[w][np.argmin(g[w])])

    (HERE / "results.json").write_text(json.dumps(res, indent=1))
    figures(R, H, centers, rb, res, gres)
    table(res)


# ----------------------------------------------------------------------------- figures
def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.6)


def figures(R, H, centers, rb, res, gres):
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "lines.linewidth": 1.6})
    groups = [("gamma1", ["plain", "kamd_def", "kamd_freq", "mbonly", "sepmon_freq"]),
              ("gamma0.1", ["plain_g0.1", "kamd_freq_g0.1"]),
              ("gamma10", ["plain_g10", "kamd_freq_g10"])]
    # --- P(r) / PMF
    for tag, cs in groups:
        cs = [c for c in cs if c in H]
        if len(cs) < 1:
            continue
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
        ref = REF_OF[cs[0]]
        Fr = None
        for c in cs:
            h = H[c]
            m = h.mean(0)
            se = h.std(0, ddof=1) / np.sqrt(len(h))
            ax[0].plot(centers, m, color=COL[c], label=LABEL[c])
            ax[0].fill_between(centers, m - se, m + se, color=COL[c], alpha=0.2, lw=0)
            # PMF and bootstrap error
            bsF = []
            for _ in range(500):
                idx = RNG.integers(0, len(h), len(h))
                F = pmf_from_p(h[idx].mean(0), centers)
                bsF.append(F)
            bsF = np.array(bsF)
            F = pmf_from_p(m, centers)
            # align on the outer well region (16.5-18.5 A) via weighted mean
            sel = (centers > 16.5) & (centers < 18.5)
            F = F - F[sel].mean()
            bsF = bsF - bsF[:, sel].mean(1, keepdims=True)
            eF = np.nanstd(np.where(np.isfinite(bsF), bsF, np.nan), 0)
            if c == ref:
                Fr, eFr = F, eF
            ax[1].plot(centers, F, color=COL[c], label=LABEL[c])
            ax[1].fill_between(centers, F - eF, F + eF, color=COL[c], alpha=0.2, lw=0)
            if Fr is not None and c != ref:
                d = F - Fr
                ax[2].errorbar(centers, d, yerr=np.hypot(eF, eFr), color=COL[c], label=LABEL[c], capsize=0,
                               lw=1.2, marker="o", ms=2.5)
        for a in ax:
            style(a)
            a.axvline(rb, color="#888", lw=0.8, ls=":")
        ax[0].set(xlabel="COM–COM distance r (Å)", ylabel="P(r) (1/Å)", title="Distance distribution")
        ax[1].set(xlabel="r (Å)", ylabel="PMF  −ln[P(r)/r²]  (kT)", title="PMF (zero = mean over 16.5–18.5 Å)",
                  ylim=(-1, 6))
        ax[2].axhline(0, color="#3a3a38", lw=0.8)
        ax[2].set(xlabel="r (Å)", ylabel="PMF − PMF(ref) (kT)", title=f"Difference vs {LABEL[ref]}", ylim=(-2, 2))
        ax[0].legend(frameon=False, fontsize=8)
        ax[2].legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG / f"pmf_{tag}.png", dpi=150)
        plt.close(fig)

    # --- dF and rates summary
    cs = [c for c in ORDER if c in res["conditions"]]
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    y = np.arange(len(cs))[::-1]
    for i, c in zip(y, cs):
        e = res["conditions"][c]
        ax[0].errorbar(e["dF"], i, xerr=e["dF_err"], fmt="o", color=COL[c], capsize=3)
        ax[1].errorbar(e["trans_per_ns"], i, xerr=e["trans_per_ns_sem"], fmt="o", color=COL[c], capsize=3)
    for a in ax:
        a.set_yticks(y)
        a.set_yticklabels([LABEL[c] for c in cs])
        style(a)
    for ref in ("plain", "plain_g0.1", "plain_g10"):
        if ref in res["conditions"]:
            ax[0].axvline(res["conditions"][ref]["dF"], color="#3a3a38", lw=0.6, ls=":")
    ax[0].set(xlabel="ΔF = F(inner) − F(outer)  (kT), ±1σ bootstrap over replicas", title="Free-energy difference between wells")
    ax[1].set(xlabel="well-to-well transitions per ns (±SEM over replicas)", title="Transition rate")
    fig.tight_layout()
    fig.savefig(FIG / "dF_rates.png", dpi=150)
    plt.close(fig)

    # --- g(r)
    cs1 = [c for c in ["plain", "kamd_def", "kamd_freq", "mbonly", "sepmon_freq"] if c in gres]
    fig, ax = plt.subplots(2, 3, figsize=(13, 6), sharex="col")
    for j, (key, title, xr) in enumerate((("gss", "solvent–solvent g(r)", "rss"),
                                          ("gua", "solute atom–solvent g(r)", "rus"),
                                          ("guc", "cluster COM–solvent g(r)", "rus"))):
        for c in cs1:
            g = gres[c][key]
            ax[0, j].plot(gres[c][xr], g.mean(0), color=COL[c], label=LABEL[c], lw=1.2)
            if c != "plain":
                d = g.mean(0) - gres["plain"][key].mean(0)
                se = np.sqrt(g.var(0, ddof=1) / len(g) + gres["plain"][key].var(0, ddof=1) / len(gres["plain"][key]))
                ax[1, j].plot(gres[c][xr], d, color=COL[c], lw=1.0)
                ax[1, j].fill_between(gres[c][xr], d - 2 * se, d + 2 * se, color=COL[c], alpha=0.15, lw=0)
        ax[0, j].set(title=title, ylabel="g(r)")
        ax[1, j].set(xlabel="r (Å)", ylabel="g − g(ref)  (band: ±2 SE)")
        ax[1, j].axhline(0, color="#3a3a38", lw=0.8)
        for a in ax[:, j]:
            style(a)
        ax[0, j].set_xlim(2, 12 if key == "gss" else 15)
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "gr.png", dpi=150)
    plt.close(fig)

    # --- Rg
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    bins = np.linspace(3.45, 3.65, 60)
    for c in cs1:
        x = np.concatenate([d["rgA"].ravel() for d in R[c]])
        ax.hist(x, bins=bins, density=True, histtype="step", color=COL[c], label=LABEL[c], lw=1.4)
    ax.set(xlabel="cluster Rg (Å)", ylabel="density", title="Cluster radius of gyration")
    style(ax)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "rg.png", dpi=150)
    plt.close(fig)

    # --- velocities vs time since reset
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    kc = [c for c in ["kamd_def", "kamd_freq", "mbonly", "sepmon_freq"] if c in res["conditions"]]
    for c in ["plain"] + kc:
        cl = res["conditions"][c]["tsr_classes"]
        if c == "plain":
            for a, key in ((ax[0], "T"), (ax[1], "KEcom_kT")):
                a.axhline(cl[0][key], color=COL[c], lw=1, ls="--", label=LABEL[c])
            continue
        x = [np.sqrt(max(k["lo_ps"], 0.05) * min(k["hi_ps"], 20)) for k in cl]
        ax[0].plot(x, [k["T"] for k in cl], "o-", color=COL[c], label=LABEL[c])
        ax[1].plot(x, [k["KEcom_kT"] for k in cl], "o-", color=COL[c], label=LABEL[c])
    ax[0].set(xscale="log", xlabel="time since last reset (ps; class midpoint)", ylabel="kinetic T (K)",
              title="Kinetic temperature (all atoms)")
    ax[1].axhline(1.5, color="#888", lw=0.8, ls=":")
    ax[1].set(xscale="log", xlabel="time since last reset (ps)", ylabel="cluster COM KE (kT)",
              title="Cluster COM kinetic energy (canonical: 1.5 kT)")
    # velocity distribution away from resets vs MB
    x = np.linspace(-6, 6, 241)
    xc = 0.5 * (x[1:] + x[:-1])
    gauss = stats.norm.pdf(xc)
    for c in ["plain"] + kc:
        uh = np.sum([d["uhist"] for d in R[c]], axis=0)
        k = uh.shape[0] - 1  # >= 3 ps since reset (or plain)
        for s, ls in ((0, "-"),):
            if uh[k, s].sum() == 0:
                continue
            h = uh[k, s] / uh[k, s].sum() / (x[1] - x[0])
            ax[2].plot(xc, h / gauss - 1, ls, color=COL[c], lw=1, label=LABEL[c])
    ax[2].set(xlim=(-4, 4), ylim=(-0.1, 0.1), xlabel="u = v·(m/kT)^½ (solvent, ≥3 ps after reset)",
              ylabel="P(u)/N(0,1) − 1", title="Velocity distribution vs Maxwell-Boltzmann")
    for a in ax:
        style(a)
    ax[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "velocities.png", dpi=150)
    plt.close(fig)

    # --- block analysis convergence for the reference
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for c in ["plain", "kamd_freq"]:
        if c not in R:
            continue
        nb = [2, 4, 5, 8, 10, 20, 40, 80]
        sem = [np.sqrt(np.sum([block_sem((d["rA"] < rb).astype(float), n) ** 2 for d in R[c]])) / len(R[c])
               for n in nb]
        L = R[c][0]["t"][-1] - R[c][0]["t"][0]
        ax.plot([L / n for n in nb], sem, "o-", color=COL[c], label=LABEL[c])
        ax.axhline(res["conditions"][c]["p_inner_sem"], color=COL[c], ls=":", lw=1)
    ax.set(xscale="log", xlabel="block length (ps)", ylabel="SEM of p(inner)",
           title="Block averaging (dotted: replica SEM)")
    style(ax)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "block.png", dpi=150)
    plt.close(fig)


def table(res):
    print(f"dividing surface r_b = {res['r_barrier']:.2f} A; solvent RDF 1st min = {res['solvent_rdf_first_min_A']:.2f} A;"
          f" tau(indicator, ref) = {res['tau_indicator_ref_ps']:.0f} ps")
    hdr = ("cond", "reps", "ns", "p_in", "dF", "±", "ddF", "z", "chi2p", "KSp", "trans/ns", "±", "ratio", "res/ns",
           "T", "Rg", "RgKSp")
    print("\t".join(hdr))
    for c, e in res["conditions"].items():
        f = lambda k, fmt="{:.3f}": fmt.format(e[k]) if e.get(k) is not None else "-"  # noqa: E731
        print("\t".join([c, str(e["n_reps"]), f"{e['ns_per_rep']:.2f}", f("p_inner"), f("dF"), f("dF_err"),
                         f("ddF"), f("ddF_z", "{:.2f}"), f("chi2_Pr_p"), f("ks_p"), f("trans_per_ns", "{:.2f}"),
                         f("trans_per_ns_sem", "{:.2f}"), f("rate_ratio", "{:.2f}"), f("resets_per_ns", "{:.0f}"),
                         f("T_mean", "{:.2f}"), f("rg_mean", "{:.4f}"), f("rg_ks_p")]))


if __name__ == "__main__":
    main()
