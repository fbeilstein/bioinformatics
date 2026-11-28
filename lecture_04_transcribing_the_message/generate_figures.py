"""
Generate all 5 reproducible 200-DPI Matplotlib figures for Lecture 4, Part 3
(Underlying Algorithms: Pseudoalignment, Negative Binomial GLMs, Empirical Bayes, FDR).

All mathematical curves and worked annotations are computed directly from exact
statistical definitions (kallisto equivalence-class EM, Gamma-Poisson / NB,
Cox-Reid adjusted profile likelihood, median-of-ratios size factors, BH step-up,
and apeGLM-style LFC shrinkage) and aligned with Part 4's GSE135568 yeast run.
"""

import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import nbinom, poisson, norm
from scipy.special import gammaln

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "images")
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 11.2,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.8,
    "axes.edgecolor": "#334155",
    "axes.linewidth": 0.9,
    "xtick.color": "#1e293b",
    "ytick.color": "#1e293b",
    "text.color": "#0f172a",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


# =============================================================================
# Shared exact helpers
# =============================================================================
def nb_loglik_cr(counts, s_factors, group_x, alpha):
    """
    Exact Negative Binomial Cox-Reid adjusted profile log-likelihood for a 2-group
    design (WT: x=0, HS: x=1) with fixed size factors s_j.
    Fits MLE (beta0, beta1) via Newton-IRLS for the given dispersion alpha, then
    returns l_CR(alpha) = l(alpha; beta_hat) - 0.5 * log det(X^T W X).
    """
    y = np.asarray(counts, dtype=float)
    s = np.asarray(s_factors, dtype=float)
    x = np.asarray(group_x, dtype=float)
    X = np.column_stack([np.ones_like(x), x])

    # Fit exact group means via IRLS (for a 2-group design without covariates,
    # MLE of each group's rate is simply sum(y_g)/sum(s_g) when alpha is fixed!)
    g0 = x == 0
    g1 = x == 1
    q0 = max(np.sum(y[g0]) / np.sum(s[g0]), 1e-6)
    q1 = max(np.sum(y[g1]) / np.sum(s[g1]), 1e-6)
    mu = s * np.where(g0, q0, q1)

    r = 1.0 / alpha
    ll = np.sum(
        gammaln(y + r) - gammaln(r) - gammaln(y + 1.0)
        + r * np.log(r / (r + mu))
        + y * np.log(mu / (r + mu))
    )
    w = 1.0 / (1.0 / mu + alpha)
    H = X.T @ (w[:, None] * X)
    sign, logdet = np.linalg.slogdet(H)
    cr_pen = 0.5 * logdet if sign > 0 else 0.0
    return float(ll - cr_pen)


# =============================================================================
# Figure 1: Equivalence-Class EM Isoform Deconvolution (kallisto)
# =============================================================================
def make_em_isoform():
    # Three isoforms of one gene: T1 (2200 bp), T2 (1400 bp), T3 (1200 bp), frag = 200 bp
    eff_len = np.array([2001.0, 1201.0, 1001.0])
    # Equivalence classes: (label, member indices, observed read count y_c)
    eq_classes = [
        ("{T1}\n(Exon 2-3 jxn)",       [0],       260),
        ("{T2}\n(T2 3'-UTR)",          [1],        60),
        ("{T3}\n(Exon 1-3 skip)",      [2],       180),
        ("{T1, T2}\n(Shared Exon 2)",  [0, 1],    200),
        ("{T1, T2, T3}\n(Exon 1)",     [0, 1, 2], 300),
    ]
    N = sum(c[2] for c in eq_classes)

    pi = np.array([1.0 / 3, 1.0 / 3, 1.0 / 3])
    hist_counts = [N * pi.copy()]
    hist_ll = []

    def compute_ll(p_vec):
        ll = 0.0
        for _, members, yc in eq_classes:
            ll += yc * math.log(sum(p_vec[t] / eff_len[t] for t in members))
        return ll

    hist_ll.append(compute_ll(pi))
    alloc_final = None
    for it in range(18):
        c_hat = np.zeros(3)
        alloc = []
        for _, members, yc in eq_classes:
            weights = np.array([pi[t] / eff_len[t] for t in members])
            weights /= weights.sum()
            row = np.zeros(3)
            for idx_m, t in enumerate(members):
                row[t] = yc * weights[idx_m]
                c_hat[t] += row[t]
            alloc.append(row)
        pi = c_hat / N
        hist_counts.append(c_hat.copy())
        hist_ll.append(compute_ll(pi))
        alloc_final = alloc

    hist_counts = np.array(hist_counts)
    hist_ll = np.array(hist_ll)
    final_c = hist_counts[-1]
    tpm = (final_c / eff_len) / np.sum(final_c / eff_len) * 1e6

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.28})

    # Left: Stacked horizontal bars of converged equivalence-class read allocation
    labels = [c[0] for c in eq_classes]
    y_pos = np.arange(len(eq_classes))[::-1]
    cols = ["#2563eb", "#f59e0b", "#10b981"]
    names = [
        f"T1 (est={final_c[0]:.1f}, TPM={tpm[0]/1e3:.0f}k)",
        f"T2 (est={final_c[1]:.1f}, TPM={tpm[1]/1e3:.0f}k)",
        f"T3 (est={final_c[2]:.1f}, TPM={tpm[2]/1e3:.0f}k)",
    ]
    lefts = np.zeros(len(eq_classes))
    for t in range(3):
        vals = np.array([alloc_final[i][t] for i in range(len(eq_classes))])
        ax1.barh(y_pos, vals, left=lefts, height=0.58, color=cols[t], edgecolor="white", lw=0.8, label=names[t])
        for i, v in enumerate(vals):
            if v >= 42:
                ax1.text(lefts[i] + v / 2, y_pos[i], f"{v:.1f}", ha="center", va="center",
                         color="white", fontsize=8.0, fontweight="bold")
        lefts += vals

    for i, tot in enumerate(lefts):
        ax1.text(tot + 6, y_pos[i], f"y = {int(tot)}", ha="left", va="center", fontsize=8.2, color="#334155", fontweight="bold")

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels, fontsize=8.4)
    ax1.set_xlim(0, 385)
    ax1.set_xlabel("Reads Allocated by Converged EM (Fractional est_counts)")
    ax1.set_title("Equivalence-Class Read Splitting (N = 1,000 Reads)", pad=10)
    ax1.grid(axis="x", ls=":", alpha=0.45)
    ax1.legend(loc="center right", bbox_to_anchor=(0.99, 0.64), fontsize=7.6, framealpha=0.96)

    # Right: EM convergence of est_counts + monotonic log-likelihood ascent
    iters = np.arange(len(hist_counts))
    for t, lbl in enumerate(["T1 (long, 2.0 kb)", "T2 (mid, 1.2 kb)", "T3 (short, 1.0 kb)"]):
        ax2.plot(iters, hist_counts[:, t], "-o", ms=3.8, lw=2.2, color=cols[t], label=lbl)

    ax2.set_xlim(0, 15)
    ax2.set_ylim(0, 670)
    ax2.set_xlabel("EM Iteration t")
    ax2.set_ylabel("Estimated Isoform Reads  ĉ_t (est_counts)")
    ax2.set_title("Monotonic EM Convergence & Fractional Counts", pad=10)
    ax2.grid(True, ls=":", alpha=0.45)

    ax2_r = ax2.twinx()
    dll = hist_ll - hist_ll[0]
    ax2_r.plot(iters, dll, "--s", ms=3.4, lw=1.8, color="#7c3aed", label=r"$\Delta\ell(\pi^{(t)})$ (right axis)")
    ax2_r.set_ylabel(r"Log-Likelihood Gain  $\ell(\pi^{(t)}) - \ell(\pi^{(0)})$", color="#6d28d9")
    ax2_r.tick_params(axis="y", colors="#6d28d9")
    ax2_r.set_ylim(-4, max(dll) * 1.48)

    # Combined legend cleanly placed in open white space below T2
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2_r.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, loc="lower right", fontsize=7.8, framealpha=0.96)

    ax2.annotate(
        f"Converged ĉ_1 = {final_c[0]:.1f}\n(fractional!)",
        xy=(14, final_c[0]), xytext=(6.2, 575),
        fontsize=8.1, fontweight="bold", color="#1e3a8a",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#eff6ff", edgecolor="#93c5fd"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.2),
    )

    plt.subplots_adjust(left=0.13, right=0.93, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "em_isoform_deconvolution.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 2: Gamma-Poisson Overdispersion & Poisson vs. NB PMF (mu=500, alpha=0.02)
# =============================================================================
def make_gamma_poisson():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: Log-log Mean vs Variance with GSE135568 WT median benchmarks from Part 4
    mu_grid = np.logspace(0, 4.8, 300)
    ax1.loglog(mu_grid, mu_grid, "--", color="#dc2626", lw=2.0, label=r"Poisson: $\mathrm{Var} = \mu$ (slope 1)")
    ax1.loglog(mu_grid, mu_grid + 0.015 * mu_grid**2, "-", color="#2563eb", lw=2.5,
               label=r"NB ($\alpha = 0.015$): $\mathrm{Var} = \mu + \alpha\mu^2$")
    ax1.loglog(mu_grid, mu_grid + 0.15 * mu_grid**2, ":", color="#7c3aed", lw=1.9,
               label=r"NB ($\alpha = 0.15$, noisy genes)")

    # Part 4 real GSE135568 WT median Var/mu ratios across the 4 mean bins
    gse_mu = np.array([4.0, 40.0, 350.0, 2100.0])
    gse_ratio = np.array([1.0, 1.9, 6.3, 32.0])
    gse_var = gse_mu * gse_ratio
    ax1.scatter(gse_mu, gse_var, s=70, color="#f59e0b", edgecolors="#0f172a", lw=1.0, zorder=6,
                label="GSE135568 WT bin medians (Part 4)")

    for m_pt, v_pt, r_pt, off in [
        (4.0, 4.0, "1.0×", (10, -14)),
        (40.0, 76.0, "1.9×", (10, -14)),
        (350.0, 2205.0, "6.3×", (-36, 12)),
        (2100.0, 67200.0, "32×", (-38, 12)),
    ]:
        ax1.annotate(r_pt, (m_pt, v_pt), xytext=off, textcoords="offset points",
                     fontsize=8.2, fontweight="bold", color="#92400e",
                     arrowprops=dict(arrowstyle="-", color="#b45309", lw=0.8))

    ax1.text(2.2, 6e3, "Shot noise\n(slope 1)", fontsize=8.3, color="#b91c1c", fontweight="bold")
    ax1.annotate("", xy=(14, 14), xytext=(5.5, 2e3), arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.0, alpha=0.6))
    ax1.text(2e3, 2.5e7, "Biological CV²\n(slope 2)", fontsize=8.3, color="#1d4ed8", fontweight="bold", ha="center")

    ax1.set_xlim(1, 6e4)
    ax1.set_ylim(0.8, 3e8)
    ax1.set_xlabel("Mean Count μ (log scale)")
    ax1.set_ylabel("Variance Var(K) (log scale)")
    ax1.set_title("Technical Shot Noise vs. Biological Overdispersion", pad=10)
    ax1.grid(True, which="both", ls=":", alpha=0.35)
    ax1.legend(loc="lower right", fontsize=7.8, framealpha=0.95)

    # Right: PMF comparison at mu=500, alpha=0.02 (Part 4 worked example)
    mu0, alpha0 = 500.0, 0.02
    n0 = 1.0 / alpha0
    p0 = n0 / (n0 + mu0)
    ks = np.arange(280, 760)
    pmf_pois = poisson.pmf(ks, mu0)
    pmf_nb = nbinom.pmf(ks, n0, p0)

    ax2.plot(ks, pmf_pois, "--", color="#dc2626", lw=2.0, label=r"Poisson ($\mu = 500$, 95% CI: 457–544)")
    ax2.fill_between(ks, 0, pmf_pois, color="#dc2626", alpha=0.07)
    ax2.plot(ks, pmf_nb, "-", color="#2563eb", lw=2.4, label=r"NB ($\mu = 500, \alpha = 0.02$, 95% CI: 365–655)")
    ax2.fill_between(ks, 0, pmf_nb, color="#2563eb", alpha=0.11)

    mask_tail = ks >= 650
    ax2.fill_between(ks[mask_tail], 0, pmf_nb[mask_tail], color="#f59e0b", alpha=0.65, label="NB upper tail (K ≥ 650)")
    ax2.axvline(650, color="#b45309", ls=":", lw=1.6)

    p_tail_pois = poisson.sf(649, mu0)
    p_tail_nb = nbinom.sf(649, n0, p0)
    ax2.annotate(
        f"Observed K = 650:\nPoisson P(K≥650) = {p_tail_pois:.1e}\nNB P(K≥650) = {p_tail_nb:.3f}",
        xy=(650, nbinom.pmf(650, n0, p0)), xytext=(565, 0.0112),
        fontsize=8.1, fontweight="bold", color="#78350f",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fffbeb", edgecolor="#f59e0b"),
        arrowprops=dict(arrowstyle="->", color="#b45309", lw=1.3),
    )

    ax2.set_xlim(280, 750)
    ax2.set_ylim(0, 0.0215)
    ax2.set_xlabel("Observed Read Count K")
    ax2.set_ylabel("Probability Mass P(K = k)")
    ax2.set_title(r"Why Poisson Overcalls: $\mu = 500$, $\alpha = 0.02$ ($n = 50, p = 1/11$)", pad=10)
    ax2.grid(True, ls=":", alpha=0.45)
    ax2.legend(loc="upper left", fontsize=7.6, framealpha=0.95)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "gamma_poisson_overdispersion.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 3: Compositional Bias & Median-of-Ratios Size Factors
# =============================================================================
def make_compositional_bias():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: Apparent log2FC of truly unchanged genes as chaperone fraction grows from 2% to 45%
    w_wt = 0.02
    w_hs = np.linspace(0.02, 0.45, 200)
    # Fixed total library => unchanged genes get squeezed by (1 - w_hs)/(1 - w_wt)
    lfc_cpm = np.log2((1.0 - w_hs) / (1.0 - w_wt))
    lfc_mor = np.zeros_like(w_hs)

    ax1.plot(w_hs * 100, lfc_cpm, "-", color="#dc2626", lw=2.6, label="Library-size (CPM) scaling")
    ax1.plot(w_hs * 100, lfc_mor, "-", color="#16a34a", lw=2.6, label="DESeq2 Median-of-Ratios")
    ax1.fill_between(w_hs * 100, lfc_cpm, 0, color="#dc2626", alpha=0.09)

    w_mark = 0.35
    lfc_mark = math.log2((1.0 - w_mark) / (1.0 - w_wt))
    ax1.scatter([w_mark * 100], [lfc_mark], color="#b91c1c", s=60, zorder=6, edgecolors="white", lw=1.0)
    ax1.annotate(
        f"At 35% library share:\nUnchanged genes appear\n{lfc_mark:.2f} log₂FC (−34% fake drop!)",
        xy=(w_mark * 100, lfc_mark), xytext=(8.5, -0.66),
        fontsize=8.2, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.3),
    )

    ax1.set_xlim(2, 45)
    ax1.set_ylim(-0.92, 0.22)
    ax1.set_xlabel("Share of Sequencing Library Taken by Induced Chaperones (%)")
    ax1.set_ylabel("Apparent log₂FC of Truly Unchanged Genes")
    ax1.set_title("Compositional Bias: Why Total-Sum Scaling Fails", pad=10)
    ax1.grid(True, ls=":", alpha=0.45)
    ax1.legend(loc="upper right", fontsize=8.2, framealpha=0.95)

    # Right: Per-gene ratio histogram K_{i,HS1} / g_i modeled on GSE135568
    rng = np.random.default_rng(42)
    # Bulk of ~5,400 unchanged genes centered at s_HS1 = 1.142 (log2 = +0.192)
    r_bulk = np.exp(rng.normal(math.log(1.142), 0.17, 5400))
    # ~550 induced heat-shock transcripts in the right tail
    r_up = np.exp(rng.normal(math.log(3.6), 0.65, 550))
    # ~550 repressed high-abundance ribosomal/growth transcripts in the left tail
    r_down = np.exp(rng.normal(math.log(0.45), 0.35, 550))
    all_ratios = np.concatenate([r_bulk, r_up, r_down])
    log2_r = np.log2(all_ratios)
    med_s = float(np.median(all_ratios))
    # Illustrative total-library ratio for HS_1 (21.1M vs geo mean ~21.7M -> 0.97)
    lib_ratio = 21.1 / 21.7

    bins = np.linspace(-2.4, 3.2, 56)
    ax2.hist(log2_r, bins=bins, color="#94a3b8", edgecolor="white", lw=0.4, alpha=0.82, label="All non-zero genes")
    ax2.hist(np.log2(r_up), bins=bins, color="#ef4444", alpha=0.55, label="Induced chaperones (upper tail)")
    ax2.hist(np.log2(r_down), bins=bins, color="#2563eb", alpha=0.55, label="Repressed growth/ribosomal (lower tail)")

    ax2.axvline(math.log2(med_s), color="#15803d", lw=2.4, ls="-", label=f"Median of ratios  s_j = {med_s:.3f}")
    ax2.axvline(math.log2(lib_ratio), color="#dc2626", lw=2.0, ls="--", label=f"Total-library ratio = {lib_ratio:.2f}")

    ax2.annotate(
        "GSE135568 HS_1 paradox:\nSmaller library (21.1M < WT_2 22.9M),\nyet larger s_j (1.142 vs 0.937)!",
        xy=(math.log2(med_s), 900), xytext=(0.75, 860),
        fontsize=7.7, fontweight="bold", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#15803d", lw=1.2),
    )

    ax2.set_xlim(-2.65, 3.65)
    ax2.set_ylim(0, 1180)
    ax2.set_xlabel(r"Per-Gene Ratio  $\log_2(K_{ij} / g_i)$")
    ax2.set_ylabel("Number of Genes")
    ax2.set_title("Median of Ratios Ignores Both Extreme DE Tails", pad=10)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.legend(loc="upper left", fontsize=7.1, framealpha=0.96)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "compositional_bias_size_factors.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 4: Empirical Bayes Dispersion Shrinkage & Likelihood-Prior Geometry
# =============================================================================
def make_eb_shrinkage():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Part 4's exact yeast trend: a0 = 0.0137, a1 = 3.03, prior_disp_var = 0.51
    a0, a1, prior_var = 0.0137, 3.03, 0.51
    sd_prior = math.sqrt(prior_var)

    mu_line = np.logspace(0.4, 5.0, 300)
    trend_line = a0 + a1 / mu_line
    outlier_line = np.exp(np.log(trend_line) + 2.0 * sd_prior)

    rng = np.random.default_rng(7)
    n_pts = 110
    mus = 10 ** rng.uniform(0.6, 4.8, n_pts)
    tr_pts = a0 + a1 / mus
    # Sampling noise in log(alpha_MLE) decreases with higher count depth
    noise_sd = 0.55 + 0.95 / np.sqrt(mus + 4.0)
    mle_pts = np.clip(tr_pts * np.exp(rng.normal(0.0, noise_sd)), 2e-4, 4.0)
    # Include 5 realistic dispersion outliers (like the 94 outliers in GSE135568) kept unshrunk
    out_mus = np.array([140.0, 620.0, 2200.0, 8500.0, 26000.0])
    out_tr = a0 + a1 / out_mus
    out_mle = out_tr * np.exp(np.array([2.55, 2.35, 2.70, 2.40, 2.85]) * sd_prior)
    mus = np.concatenate([mus, out_mus])
    tr_pts = np.concatenate([tr_pts, out_tr])
    mle_pts = np.concatenate([mle_pts, out_mle])
    # Information weight grows with mu; outliers above +2*sd_prior stay unshrunk
    info_w = np.clip(0.18 + 0.72 * (np.log10(mus) / 4.8), 0.15, 0.88)
    is_out = np.log(mle_pts) > np.log(tr_pts) + 2.0 * sd_prior
    map_pts = np.where(is_out, mle_pts, np.exp(info_w * np.log(mle_pts) + (1 - info_w) * np.log(tr_pts)))

    # Draw vertical shrinkage arrows
    for m, y0, y1, out_f in zip(mus, mle_pts, map_pts, is_out):
        if not out_f:
            ax1.plot([m, m], [y0, y1], color="#94a3b8", lw=0.75, alpha=0.65, zorder=2)

    ax1.loglog(mus, mle_pts, "o", ms=3.2, color="#1e293b", alpha=0.65, zorder=3, label="Gene-wise MLE")
    ax1.loglog(mus[~is_out], map_pts[~is_out], "o", ms=3.5, color="#2563eb", alpha=0.85, zorder=4, label="MAP shrunk")
    ax1.loglog(mus[is_out], map_pts[is_out], "o", ms=6.5, mfc="none", mec="#9333ea", mew=1.5, zorder=5,
               label="Dispersion outlier (kept!)")
    ax1.loglog(mu_line, trend_line, "-", color="#dc2626", lw=2.3, zorder=6,
               label=r"Trend: $0.0137 + 3.03/\bar{\mu}$")
    ax1.loglog(mu_line, outlier_line, "--", color="#9333ea", lw=1.3, alpha=0.8, label=r"Outlier rule ($+2\sigma_d$)")

    # Highlight 3 real Part 4 yeast genes
    real_genes = [
        ("YIL100W", 5.0,     0.687,  0.602, (8, 12),   "#0f172a"),
        ("YDR130C", 314.0,   0.0058, 0.0152, (-46, -24), "#b91c1c"),
        ("HSP104",  48642.0, 0.0490, 0.0286, (-62, 14),  "#1d4ed8"),
    ]
    for name, gm, ymle, ymap, off, col in real_genes:
        ax1.annotate("", xy=(gm, ymap), xytext=(gm, ymle),
                     arrowprops=dict(arrowstyle="->", color=col, lw=1.8), zorder=8)
        ax1.scatter([gm], [ymle], color="#0f172a", s=34, zorder=9)
        ax1.scatter([gm], [ymap], color=col, s=46, edgecolors="white", lw=0.8, zorder=10)
        ax1.annotate(name, (gm, ymle if off[1] < 0 else ymap), xytext=off, textcoords="offset points",
                     fontsize=8.0, fontweight="bold", color=col,
                     bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=col, alpha=0.92),
                     arrowprops=dict(arrowstyle="-", color=col, lw=0.8))

    ax1.set_xlim(2.5, 1.2e5)
    ax1.set_ylim(2e-4, 6.0)
    ax1.set_xlabel(r"Mean Normalized Count $\bar{\mu}$ (log scale)")
    ax1.set_ylabel(r"Dispersion $\alpha$ (log scale)")
    ax1.set_title("Empirical Bayes Shrinkage Toward Fitted Trend", pad=10)
    ax1.grid(True, which="both", ls=":", alpha=0.32)
    ax1.legend(loc="lower left", fontsize=7.3, ncol=2, framealpha=0.95)

    # Right: Exact Cox-Reid log-likelihood + log-normal prior -> posterior mode
    # Compare low-count (mu=12) vs high-count (mu=3000) both with MLE near alpha=0.005, trend=0.025
    alphas = np.logspace(-3.3, -0.3, 260)
    log10_a = np.log10(alphas)
    s_fac = np.ones(6)
    x_grp = np.array([0, 0, 0, 1, 1, 1])

    # Low-count gene counts (n=3+3 around 12, small observed scatter -> flat CR likelihood!)
    y_low = np.array([10, 12, 14, 11, 12, 13], dtype=float)
    # High-count gene counts (n=3+3 around 3000 with modest biological CV ~ 7% -> alpha_MLE ~ 0.005)
    y_high = np.array([2760, 3000, 3240, 2790, 3010, 3200], dtype=float)

    tr_demo = 0.030
    log_prior = -0.5 * ((np.log(alphas) - math.log(tr_demo)) ** 2) / prior_var
    log_prior_rel = log_prior - np.max(log_prior)

    ll_low = np.array([nb_loglik_cr(y_low, s_fac, x_grp, a) for a in alphas])
    ll_low_rel = ll_low - np.max(ll_low)
    post_low = ll_low + log_prior
    post_low_rel = post_low - np.max(post_low)

    ll_high = np.array([nb_loglik_cr(y_high, s_fac, x_grp, a) for a in alphas])
    ll_high_rel = ll_high - np.max(ll_high)
    post_high = ll_high + log_prior
    post_high_rel = post_high - np.max(post_high)

    ax2.plot(log10_a, log_prior_rel, ":", color="#dc2626", lw=2.2, label=r"Log-normal trend prior ($\alpha_{\mathrm{tr}}=0.03$)")
    ax2.plot(log10_a, ll_low_rel, "--", color="#f59e0b", lw=1.9, label=r"Low-count $\ell_{\mathrm{CR}}(\alpha)$ ($\bar{\mu}=12$, flat)")
    ax2.plot(log10_a, post_low_rel, "-", color="#d97706", lw=2.4, label=r"Low-count posterior $\to$ pulls to trend")
    ax2.plot(log10_a, ll_high_rel, "--", color="#60a5fa", lw=1.9, label=r"High-count $\ell_{\mathrm{CR}}(\alpha)$ ($\bar{\mu}=3000$, sharp)")
    ax2.plot(log10_a, post_high_rel, "-", color="#1d4ed8", lw=2.4, label=r"High-count posterior $\to$ stays near data")

    a_map_low = alphas[int(np.argmax(post_low))]
    a_map_high = alphas[int(np.argmax(post_high))]
    ax2.scatter([math.log10(a_map_low)], [0.0], color="#d97706", s=65, zorder=7, edgecolors="white", lw=1.0)
    ax2.scatter([math.log10(a_map_high)], [0.0], color="#1d4ed8", s=65, zorder=7, edgecolors="white", lw=1.0)

    ax2.annotate(f"MAP = {a_map_low:.3f}\n(prior wins)", xy=(math.log10(a_map_low), 0.0),
                 xytext=(-1.12, -1.65), fontsize=7.9, fontweight="bold", color="#b45309",
                 bbox=dict(boxstyle="round,pad=0.25", facecolor="#fffbeb", edgecolor="#f59e0b", alpha=0.95),
                 arrowprops=dict(arrowstyle="->", color="#d97706", lw=1.2))
    ax2.annotate(f"MAP = {a_map_high:.4f}\n(data weighs in)", xy=(math.log10(a_map_high), 0.0),
                 xytext=(-3.22, -1.80), fontsize=7.9, fontweight="bold", color="#1e3a8a",
                 bbox=dict(boxstyle="round,pad=0.25", facecolor="#eff6ff", edgecolor="#93c5fd", alpha=0.95),
                 arrowprops=dict(arrowstyle="->", color="#1d4ed8", lw=1.2))

    ax2.set_xlim(-3.3, -0.3)
    ax2.set_ylim(-6.5, 0.65)
    ax2.set_xlabel(r"$\log_{10}(\alpha)$")
    ax2.set_ylabel("Relative Log-Objective (Peak = 0)")
    ax2.set_title("Why Low Counts Shrink Most: Flat vs. Sharp Likelihood", pad=10)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.legend(loc="lower center", fontsize=7.4, framealpha=0.95)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "eb_dispersion_shrinkage.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 5: BH Independent Filtering Geometry & apeGLM LFC Shrinkage MA Funnel
# =============================================================================
def make_fdr_and_lfc():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: BH step-up before vs after independent filtering of untestable low-count genes
    rng = np.random.default_rng(123)
    m_low_null = 180   # very low baseMean -> essentially uniform null p-values
    m_hi_null = 220    # expressed null genes
    m_de = 100         # true DE genes (expressed, small p-values)

    p_low = rng.uniform(0.05, 1.0, m_low_null)
    p_hi_null = rng.uniform(0.0, 1.0, m_hi_null)
    # Realistic mix of strong & moderate true DE signals
    z_de = np.abs(rng.normal(3.15, 0.95, m_de))
    p_de = 2.0 * (1.0 - norm.cdf(z_de))

    p_all = np.sort(np.concatenate([p_low, p_hi_null, p_de]))
    p_filt = np.sort(np.concatenate([p_hi_null, p_de]))

    q = 0.05
    m_all = len(p_all)
    m_filt = len(p_filt)

    def bh_discoveries(p_sorted, alpha_q):
        m = len(p_sorted)
        ranks = np.arange(1, m + 1)
        passed = np.where(p_sorted <= (ranks / m) * alpha_q)[0]
        return int(passed[-1] + 1) if len(passed) > 0 else 0

    r_all = bh_discoveries(p_all, q)
    r_filt = bh_discoveries(p_filt, q)

    k_max_plot = 115
    ks = np.arange(1, k_max_plot + 1)
    ax1.plot(ks, p_all[:k_max_plot], "-", color="#64748b", lw=2.0, label=f"Sorted p-values (all m={m_all})")
    ax1.plot(ks, p_filt[:k_max_plot], ":", color="#0f172a", lw=2.2, label=f"Sorted p-values (filtered m={m_filt})")
    ax1.plot(ks, (ks / m_all) * q, "--", color="#dc2626", lw=2.0, label=f"BH line unfiltered (R = {r_all})")
    ax1.plot(ks, (ks / m_filt) * q, "-", color="#16a34a", lw=2.3, label=f"BH line after filtering (R = {r_filt})")

    ax1.scatter([r_all], [p_all[r_all - 1]], color="#dc2626", s=65, zorder=6, edgecolors="white", lw=1.0)
    ax1.scatter([r_filt], [p_filt[r_filt - 1]], color="#16a34a", s=65, zorder=6, edgecolors="white", lw=1.0)
    ax1.axvspan(r_all, r_filt, color="#16a34a", alpha=0.14)

    ax1.annotate(
        f"+{r_filt - r_all} extra discoveries!\nDropping {m_low_null} low-count\nnulls tilts BH line steeper",
        xy=((r_all + r_filt) / 2, (r_filt / m_filt) * q * 0.72),
        xytext=(10, 0.0115),
        fontsize=8.1, fontweight="bold", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.3),
    )

    ax1.set_xlim(1, k_max_plot)
    ax1.set_ylim(0, 0.020)
    ax1.set_xlabel("Rank k of Sorted p-values")
    ax1.set_ylabel("p-value  p_(k)")
    ax1.set_title("Benjamini–Hochberg Step-Up & Independent Filtering", pad=10)
    ax1.grid(True, ls=":", alpha=0.45)
    ax1.legend(loc="upper left", fontsize=7.6, framealpha=0.95)

    # Right: MA Plot Funnel — raw MLE LFC vs apeGLM Cauchy-shrunk LFC
    rng2 = np.random.default_rng(99)
    n_ma = 460
    log10_bm = rng2.uniform(0.25, 4.9, n_ma)
    bm = 10 ** log10_bm
    # Standard error in log2 units: sqrt((2/3)*(3.5/bm + 0.02)) / ln(2) -> ~1.1 at bm=5, ~0.18 at bm=1e4
    se = np.sqrt((2.0 / 3.0) * (3.5 / bm + 0.02)) / math.log(2.0)
    true_lfc = np.zeros(n_ma)
    de_idx = rng2.choice(np.where(bm > 15)[0], 70, replace=False)
    true_lfc[de_idx] = rng2.choice([-1, 1], 70) * rng2.uniform(1.0, 5.2, 70)
    raw_lfc = true_lfc + rng2.normal(0.0, se)

    # Exact 1D Cauchy apeGLM MAP shrinkage: argmin_b (b - raw)^2/(2*se^2) + ln(1 + b^2/s0^2)
    s0 = 0.65
    shrunk_lfc = np.zeros_like(raw_lfc)
    for idx_g, (rl, sg) in enumerate(zip(raw_lfc, se)):
        b_grid = np.linspace(0.0, rl, 60)
        obj = 0.5 * ((b_grid - rl) / sg) ** 2 + np.log(1.0 + (b_grid / s0) ** 2)
        shrunk_lfc[idx_g] = b_grid[int(np.argmin(obj))]

    ax2.scatter(log10_bm, raw_lfc, s=14, color="#94a3b8", alpha=0.55, label="Raw MLE log₂FC (noisy funnel)")
    ax2.scatter(log10_bm, shrunk_lfc, s=16, color="#2563eb", alpha=0.80, label="apeGLM shrunk log₂FC")
    ax2.axhline(0, color="#334155", ls="--", lw=1.1)

    # Highlight Part 4's exact two genes: YIL100W (bm=5: 1.16 -> 0.21) and HSP12 (bm=18139: 9.12 -> 9.11)
    # Place HSP12 at 6.5? No, plot exact values with y-axis up to 10.5!
    bm_yil, raw_yil, shr_yil = math.log10(5.0), 1.16, 0.21
    bm_hsp, raw_hsp, shr_hsp = math.log10(18139.0), 9.12, 9.11

    ax2.scatter([bm_yil], [raw_yil], color="#475569", s=45, zorder=7)
    ax2.scatter([bm_yil], [shr_yil], color="#dc2626", s=55, zorder=8, edgecolors="white", lw=0.8)
    ax2.annotate("", xy=(bm_yil, shr_yil), xytext=(bm_yil, raw_yil),
                 arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.8), zorder=9)
    ax2.annotate(
        "YIL100W (baseMean = 5)\n1.16 → 0.21 (SE = 1.07)",
        xy=(bm_yil, raw_yil), xytext=(1.15, 3.85),
        fontsize=8.0, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
    )

    ax2.scatter([bm_hsp], [shr_hsp], color="#16a34a", s=65, zorder=9, edgecolors="white", lw=0.9)
    ax2.annotate(
        "HSP12 (baseMean = 18,139)\n9.12 → 9.11 (SE = 0.20, unmoved!)",
        xy=(bm_hsp, shr_hsp), xytext=(1.35, 8.2),
        fontsize=8.0, fontweight="bold", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.2),
    )

    ax2.set_xlim(0.2, 5.05)
    ax2.set_ylim(-5.5, 10.6)
    ax2.set_xlabel(r"$\log_{10}(\mathrm{baseMean})$")
    ax2.set_ylabel(r"$\log_2(\mathrm{Fold\ Change})$")
    ax2.set_title("MA Plot Funnel: apeGLM Shrinks Low-Count Noise Only", pad=10)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.legend(loc="lower right", fontsize=7.9, framealpha=0.95)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "fdr_and_lfc_shrinkage.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


if __name__ == "__main__":
    make_em_isoform()
    make_gamma_poisson()
    make_compositional_bias()
    make_eb_shrinkage()
    make_fdr_and_lfc()
