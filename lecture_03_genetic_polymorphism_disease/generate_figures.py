#!/usr/bin/env python3
"""Generate 5 publication-quality figures for Lecture 3, Part 3 (Underlying Algorithms)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from scipy.stats import binom

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "images")
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.titlesize": 11.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "axes.edgecolor": "#475569",
    "axes.linewidth": 0.9,
    "xtick.color": "#334155",
    "ytick.color": "#334155",
    "text.color": "#1e293b",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 200,
})


# =============================================================================
# Figure 1: BGZF Independent Blocks + UCSC/BAI Hierarchical Binning Index
# =============================================================================
def make_bgzf_binning():
    fig, ax = plt.subplots(figsize=(12.8, 4.6))
    ax.set_xlim(0, 12.8)
    ax.set_ylim(0, 4.6)
    ax.axis("off")

    ax.text(6.4, 4.38, "Random-Access Streaming: BAI / Tabix Hierarchical Binning + BGZF 64 KB Blocks",
            ha="center", va="center", fontsize=12.0, fontweight="bold", color="#0f172a")

    # --- Top Half: Hierarchical Binning Tree over chr17 (83.3 Mb) ---
    levels = [
        ("Level 1 (64 Mb)",  3.62, 2,  [0],    ["Bin 1 (0–64 Mb)", "Bin 2 (64–128 Mb)"]),
        ("Level 2 (8 Mb)",   2.98, 10, [5],    None),
        ("Level 3 (1 Mb)",   2.34, 24, [12],   None),
        ("Level 5 (16 kb)",  1.70, 48, [24, 25], None),
    ]
    x0, x1 = 1.95, 12.45
    W = x1 - x0

    for lbl, y, nbins, hits, custom_names in levels:
        ax.text(x0 - 0.18, y + 0.18, lbl, ha="right", va="center",
                fontsize=9.2, fontweight="bold", color="#334155")
        bw = W / nbins
        for b in range(nbins):
            bx = x0 + b * bw
            is_hit = b in hits
            fc = "#2563eb" if is_hit else "#f1f5f9"
            ec = "#1d4ed8" if is_hit else "#cbd5e1"
            rect = Rectangle((bx, y), bw - 0.02, 0.36, facecolor=fc, edgecolor=ec, lw=1.0 if is_hit else 0.7)
            ax.add_patch(rect)
            if custom_names:
                ax.text(bx + bw / 2, y + 0.18, custom_names[b], ha="center", va="center",
                        fontsize=8.6, fontweight="bold" if is_hit else "normal",
                        color="white" if is_hit else "#475569")

    # Annotate exact SAMv1 bin IDs for BRCA1 (chr17:43,044,295-43,170,245)
    ax.text(x0 + 5.5 * (W / 10), 2.98 + 0.18, "Bin 14", ha="center", va="center",
            fontsize=8.3, fontweight="bold", color="white")
    ax.annotate("Bin 114\n(43–44 Mb)", xy=(x0 + 12.5 * (W / 24), 2.52), xytext=(8.6, 2.50),
                fontsize=8.4, fontweight="bold", color="#1d4ed8", va="center",
                arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.2))
    ax.annotate("Bins 7308–7315\n(8 × 16 kb = BRCA1)", xy=(x0 + 25.0 * (W / 48), 1.88), xytext=(8.6, 1.88),
                fontsize=8.4, fontweight="bold", color="#1d4ed8", va="center",
                arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.2))

    # --- Bottom Half: Compressed BGZF / CRAM Stream (64 KB blocks) ---
    ax.text(x0 - 0.18, 0.72, "Remote CRAM\n(BGZF Blocks)", ha="right", va="center",
            fontsize=9.2, fontweight="bold", color="#7c3aed")

    n_blocks = 28
    blw = W / n_blocks
    hit_blocks = {13, 14}
    for b in range(n_blocks):
        bx = x0 + b * blw
        is_hit = b in hit_blocks
        fc = "#7c3aed" if is_hit else "#e2e8f0"
        ec = "#5b21b6" if is_hit else "#94a3b8"
        rect = Rectangle((bx, 0.45), blw - 0.03, 0.52, facecolor=fc, edgecolor=ec, lw=1.3 if is_hit else 0.7)
        ax.add_patch(rect)

    # Virtual file offset arrow from Level 5 hit bins to BGZF blocks
    bin_cx = x0 + 24.9 * (W / 48)
    blk_cx = x0 + 14.0 * blw
    ax.annotate(
        "", xy=(blk_cx, 1.00), xytext=(bin_cx, 1.66),
        arrowprops=dict(arrowstyle="->", color="#7c3aed", lw=2.2)
    )
    ax.text((bin_cx + blk_cx) / 2 - 0.15, 1.33, "Virtual offset\n(block << 16 | intra)", ha="right", va="center",
            fontsize=8.4, fontweight="bold", color="#6d28d9")

    # Summary callout box on bottom right
    ax.text(
        9.95, 0.16,
        "HTTP Range Request fetches ONLY highlighted blocks (~2 MB)\nout of a multi-GB CRAM — skipping >99.9% of the file!",
        ha="center", va="center", fontsize=8.6, fontweight="bold", color="#4c1d95",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#f5f3ff", edgecolor="#c4b5fd", lw=1.2)
    )
    ax.text(
        3.85, 0.16,
        "Each BGZF block ≤ 64 KB is an independent gzip member\n→ seek directly without decompressing bytes 0..N",
        ha="center", va="center", fontsize=8.6, color="#334155",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8fafc", edgecolor="#cbd5e1", lw=1.0)
    )

    plt.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.03)
    out = os.path.join(OUT_DIR, "bgzf_binning_index.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 2: Prior vs. Data Interplay (Posterior P(1/1|D) + 2D Phase Diagram)
# =============================================================================
def bayes_posteriors(n_ref, n_alt, q=30.0, f=0.01, flat=False):
    eps = 10.0 ** (-q / 10.0)
    ll00 = n_ref * np.log(1.0 - eps) + n_alt * np.log(eps)
    ll01 = (n_ref + n_alt) * np.log(0.5)
    ll11 = n_ref * np.log(eps) + n_alt * np.log(1.0 - eps)
    if flat:
        lp = np.array([0.0, 0.0, 0.0])
    else:
        lp = np.log(np.array([(1.0 - f) ** 2, 2.0 * f * (1.0 - f), f ** 2]))
    lpost = np.array([ll00, ll01, ll11]) + lp
    m = np.max(lpost)
    p = np.exp(lpost - m)
    return p / np.sum(p)


def make_prior_vs_depth():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.7), gridspec_kw={"wspace": 0.25})

    # Left panel: P(G=1/1 | 0 REF, n ALT at Q30) across priors
    ns = np.arange(1, 17)
    configs = [
        ("Flat prior (uninformative)", None, True,  "#64748b", "--"),
        ("HWE f = 0.35 (common)",      0.35, False, "#16a34a", "-"),
        ("HWE f = 0.05",               0.05, False, "#2563eb", "-"),
        ("HWE f = 0.01 (rare)",        0.01, False, "#dc2626", "-"),
        ("HWE f = 0.001 (very rare)",  0.001,False, "#7c3aed", "-."),
    ]
    for lbl, f, flat, col, ls in configs:
        p11 = [bayes_posteriors(0, int(n), q=30.0, f=(f or 0.01), flat=flat)[2] for n in ns]
        ax1.plot(ns, p11, marker="o", ms=4.5, lw=2.2, color=col, ls=ls, label=lbl)

    ax1.axhline(0.5, color="#334155", ls=":", lw=1.2)
    ax1.text(15.8, 0.525, "MAP boundary (0/1 ↔ 1/1)", ha="right", fontsize=8.4, color="#334155", fontweight="bold")

    # Annotate DP=2 flip from Part 4
    ax1.axvline(2, color="#f59e0b", ls="--", lw=1.3, alpha=0.8)
    ax1.annotate(
        "DP = 2 (rs12516):\nFlat → 1/1 (80%)\nHWE f=0.01 → 0/1 (98%!)",
        xy=(2, 0.02), xytext=(4.3, 0.14),
        fontsize=8.4, fontweight="bold", color="#92400e",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fffbeb", edgecolor="#fcd34d"),
        arrowprops=dict(arrowstyle="->", color="#d97706", lw=1.3),
    )

    ax1.set_xlim(0.7, 16.3)
    ax1.set_ylim(-0.03, 1.05)
    ax1.set_xlabel("Pure ALT Read Depth n (0 REF, n ALT at Q30)")
    ax1.set_ylabel("Posterior Probability P(G = 1/1 | D)")
    ax1.set_title("Prior vs. Data: Depth Needed to Call Homozygous 1/1", pad=10)
    ax1.grid(True, ls=":", alpha=0.45)
    ax1.legend(loc="lower right", fontsize=8.1, framealpha=0.95)

    # Right panel: 2D MAP Phase Diagram over (n_alt, log10 f)
    logf_grid = np.linspace(-4.0, np.log10(0.5), 250)
    f_grid = 10.0 ** logf_grid
    eps = 10.0 ** (-3.0)
    n_star = np.log(2.0 * (1.0 - f_grid) / f_grid) / np.log(2.0 * (1.0 - eps))

    ax2.fill_betweenx(logf_grid, 0, np.minimum(n_star, 16.5), color="#dbeafe", alpha=0.85, label="MAP = 0/1 (Prior wins)")
    ax2.fill_betweenx(logf_grid, np.minimum(n_star, 16.5), 16.5, color="#fee2e2", alpha=0.85, label="MAP = 1/1 (Data win)")
    ax2.plot(n_star, logf_grid, color="#0f172a", lw=2.5, label=r"$n^* = \log_2\!\left(\frac{2(1-f)}{f}\right)$")

    # Mark Part 4 examples
    pts = [
        (2, np.log10(0.01), "DP=2, f=0.01\n→ 0/1", (12, -20), "#1d4ed8"),
        (2, np.log10(0.35), "DP=2, f=0.35\n→ 1/1", (14, -4),  "#b91c1c"),
        (14, np.log10(0.01), "DP≥15\nData dominate", (-62, -28), "#b91c1c"),
    ]
    for x, y, txt, off, col in pts:
        ax2.scatter([x], [y], color=col, s=55, zorder=6, edgecolors="white", lw=1.1)
        ax2.annotate(txt, (x, y), xytext=off, textcoords="offset points",
                     fontsize=8.4, fontweight="bold", color=col,
                     arrowprops=dict(arrowstyle="-", color=col, lw=0.9))

    ax2.set_xlim(0.5, 16.0)
    ax2.set_ylim(-4.0, np.log10(0.55))
    ax2.set_yticks([-4, -3, -2, -1, np.log10(0.5)])
    ax2.set_yticklabels(["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "0.5"])
    ax2.set_xlabel("Pure ALT Read Depth n (0 REF, n ALT at Q30)")
    ax2.set_ylabel("Population Allele Frequency f (log scale)")
    ax2.set_title("MAP Decision Boundary: Exact Crossover Depth n*", pad=10)
    ax2.grid(True, ls=":", alpha=0.45)
    ax2.legend(loc="lower right", fontsize=8.3, framealpha=0.95)

    plt.subplots_adjust(left=0.06, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "bayesian_prior_vs_depth.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 3: IEEE-754 Float64 Underflow vs. Log-Space Likelihoods
# =============================================================================
def make_log_underflow():
    fig, ax = plt.subplots(figsize=(7.2, 4.7))

    n = np.linspace(0, 1250, 500)
    log10_p = n * np.log10(0.5)

    ax.plot(n, log10_p, color="#2563eb", lw=2.8, label=r"True $\log_{10} P(D \mid 0/1) = n \log_{10}(0.5)$")
    ax.axhline(-308.25, color="#dc2626", ls="--", lw=1.8, label="IEEE-754 float64 min normal (2.2 × 10⁻³⁰⁸)")
    ax.axhline(-323.6, color="#991b1b", ls=":", lw=1.5, label="Float64 subnormal flush-to-zero (5 × 10⁻³²⁴)")
    ax.axhspan(-390, -308.25, color="#fee2e2", alpha=0.45)

    # Mark n = 1024 and BRCA1 peak depth = 1111
    ax.axvline(1111, color="#7c3aed", ls="-.", lw=1.8)
    ax.text(1100, -45, "BRCA1 max depth\n(1,111× reads)", ha="right", va="top",
            fontsize=8.8, fontweight="bold", color="#6d28d9",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f5f3ff", edgecolor="#c4b5fd"))

    ax.scatter([1024], [1024 * np.log10(0.5)], color="#dc2626", s=65, zorder=6, edgecolors="white", lw=1.2)
    ax.annotate(
        "n = 1,024: raw (0.5)ⁿ underflows!\nAt n = 1,075: 0.0 / 0.0 → NaN",
        xy=(1024, -308.25), xytext=(220, -245),
        fontsize=8.8, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.4),
    )

    ax.set_xlim(0, 1250)
    ax.set_ylim(-385, 15)
    ax.set_xlabel("Read Depth n at a Heterozygous Site")
    ax.set_ylabel(r"$\log_{10} P(D \mid G = 0/1)$")
    ax.set_title("Why Variant Callers Must Compute in Log Space", pad=10)
    ax.grid(True, ls=":", alpha=0.45)
    ax.legend(loc="lower left", fontsize=8.2, framealpha=0.95)

    plt.subplots_adjust(left=0.13, right=0.96, top=0.90, bottom=0.13)
    out = os.path.join(OUT_DIR, "log_underflow_lse.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 4: Site QUAL vs. Sample GQ Paradox (0 REF, n ALT at Q30, f = 0.01)
# =============================================================================
def make_qual_vs_gq():
    ns = np.arange(1, 16)
    quals = []
    gqs = []
    calls = []
    for n in ns:
        p = bayes_posteriors(0, int(n), q=30.0, f=0.01, flat=False)
        p00 = max(p[0], 1e-300)
        qual = -10.0 * np.log10(p00)
        best_p = np.max(p)
        err_p = max(1.0 - best_p, 1e-300)
        gq = -10.0 * np.log10(err_p)
        quals.append(qual)
        gqs.append(gq)
        calls.append(int(np.argmax(p)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.24})

    # Left: Site QUAL
    ax1.plot(ns, quals, "-o", color="#2563eb", lw=2.5, ms=5.5, label=r"$\mathrm{QUAL} = -10\log_{10} P(0/0 \mid D)$")
    ax1.axhline(20, color="#16a34a", ls="--", lw=1.4, label="Filter threshold QUAL ≥ 20")
    ax1.scatter([5], [quals[4]], color="#dc2626", s=75, zorder=6, edgecolors="white", lw=1.2)
    ax1.annotate(
        f"n = 5: QUAL = {quals[4]:.1f}\n(~10⁻¹² chance of 0/0!)",
        xy=(5, quals[4]), xytext=(6.8, 70),
        fontsize=8.7, fontweight="bold", color="#1e3a8a",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#eff6ff", edgecolor="#93c5fd"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.3),
    )
    ax1.set_xlim(0.5, 15.5)
    ax1.set_ylim(0, 440)
    ax1.set_xlabel("Pure ALT Read Depth n (0 REF, n ALT at Q30, f = 0.01)")
    ax1.set_ylabel("Site QUAL (Phred)")
    ax1.set_title("Site Confidence (QUAL): 'Is There Any Variant?'", pad=10)
    ax1.grid(True, ls=":", alpha=0.45)
    ax1.legend(loc="upper left", fontsize=8.4, framealpha=0.95)

    # Right: Genotype GQ (shows the dramatic dip at the 0/1 <-> 1/1 crossover!)
    colors_pt = ["#2563eb" if c == 1 else "#dc2626" for c in calls]
    ax2.plot(ns, gqs, "-", color="#475569", lw=2.0, zorder=2)
    ax2.scatter(ns, gqs, c=colors_pt, s=55, zorder=5, edgecolors="white", lw=1.0)
    ax2.axhline(20, color="#16a34a", ls="--", lw=1.4, label="High genotype confidence (GQ ≥ 20)")
    ax2.axvline(7.63, color="#f59e0b", ls=":", lw=1.8, label="Prior/data crossover n* = 7.6")

    # Legend proxies for call color
    ax2.scatter([], [], color="#2563eb", s=45, label="MAP call = 0/1")
    ax2.scatter([], [], color="#dc2626", s=45, label="MAP call = 1/1")

    ax2.annotate(
        f"n = 5–8: QUAL = 119–204,\nbut GQ plunges to {min(gqs):.1f}!\n(Sure it's variant, unsure 0/1 vs 1/1)",
        xy=(8, gqs[7]), xytext=(0.9, 23.2),
        fontsize=8.3, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.3),
    )

    ax2.set_xlim(0.5, 15.5)
    ax2.set_ylim(0, 32)
    ax2.set_xlabel("Pure ALT Read Depth n (0 REF, n ALT at Q30, f = 0.01)")
    ax2.set_ylabel("Genotype Quality GQ (Phred)")
    ax2.set_title("Genotype Confidence (GQ): '0/1 or 1/1?'", pad=10)
    ax2.grid(True, ls=":", alpha=0.45)
    ax2.legend(loc="upper right", fontsize=7.9, framealpha=0.95)

    plt.subplots_adjust(left=0.06, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "qual_vs_gq_paradox.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 5: Heterozygote Allelic Dropout & Binomial AAF Spread
# =============================================================================
def make_dropout_aaf():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.24})

    # Left: Allelic dropout P(k=0 or k=n) = 2*(0.5)^n on log scale
    ns = np.arange(1, 31)
    p_drop = 2.0 * (0.5 ** ns)
    p_miss_f20 = [binom.cdf(int(np.ceil(0.2 * n) - 1), n, 0.5) for n in ns]

    ax1.semilogy(ns, p_drop, "-o", ms=4.0, color="#dc2626", lw=2.4, label=r"Allelic dropout: $2(1/2)^n$ (0/1 looks 0/0 or 1/1)")
    ax1.semilogy(ns, p_miss_f20, "--", color="#f59e0b", lw=2.2, label=r"Missed by freebayes -F 0.2: $P(k/n < 0.2)$")

    benchmarks = [
        (4,  2 * 0.5**4,  "DP = 4 (BRCA1 median)\n12.5% dropout!", (18, 10), "#b91c1c"),
        (10, 2 * 0.5**10, "DP = 10\n0.20%", (-54, -26), "#334155"),
        (15, 2 * 0.5**15, "DP = 15 (Homework)\n0.006% dropout", (26, -6), "#15803d"),
    ]
    for x, y, txt, off, col in benchmarks:
        ax1.scatter([x], [y], color=col, s=55, zorder=6, edgecolors="white", lw=1.1)
        ax1.annotate(txt, (x, y), xytext=off, textcoords="offset points",
                     fontsize=8.3, fontweight="bold", color=col,
                     arrowprops=dict(arrowstyle="-", color=col, lw=0.9))

    ax1.axvline(15, color="#16a34a", ls=":", lw=1.5, alpha=0.8)
    ax1.set_xlim(1, 30)
    ax1.set_ylim(1e-9, 3.5)
    ax1.set_xlabel("Independent Read Depth n at a Heterozygous Site")
    ax1.set_ylabel("Failure Probability (log scale)")
    ax1.set_title("Heterozygote Sampling Failure vs. Read Depth", pad=10)
    ax1.grid(True, which="both", ls=":", alpha=0.38)
    ax1.legend(loc="lower left", fontsize=8.1, framealpha=0.95)

    # Right: Binomial AAF distribution for n = 5, 15, 40, 150 (with subtle ref bias p=0.48 shown at n=150)
    depths = [
        (5,   0.50, "#dc2626", "DP = 5   (SD = 0.224)"),
        (15,  0.50, "#f59e0b", "DP = 15  (SD = 0.129)"),
        (40,  0.50, "#2563eb", "DP = 40  (SD = 0.079)"),
        (150, 0.48, "#7c3aed", "DP = 150, ref bias p=0.48"),
    ]
    for n, p_alt, col, lbl in depths:
        ks = np.arange(0, n + 1)
        aafs = ks / n
        # Density-scaled PMF so curves at different n share comparable height
        dens = binom.pmf(ks, n, p_alt) * n
        ax2.plot(aafs, dens, lw=2.2, color=col, label=lbl)
        ax2.fill_between(aafs, 0, dens, color=col, alpha=0.10)

    ax2.axvline(0.5, color="#334155", ls="--", lw=1.3)
    ax2.axvline(0.2, color="#dc2626", ls=":", lw=1.4)
    ax2.text(0.21, 10.6, "freebayes\n-F 0.2", ha="left", va="top", fontsize=8.3, fontweight="bold", color="#b91c1c")
    ax2.text(0.515, 10.8, "Ideal het (0.50)", ha="left", va="top", fontsize=8.3, color="#334155", fontweight="bold")

    ax2.set_xlim(0.0, 1.0)
    ax2.set_ylim(0, 11.4)
    ax2.set_xlabel("Observed Alternate Allele Fraction  AAF = AO / (RO + AO)")
    ax2.set_ylabel("Probability Density")
    ax2.set_title(r"Heterozygote AAF Narrowing ($\mathrm{SD} = 1 / 2\sqrt{n}$) & Ref Bias", pad=10)
    ax2.grid(True, ls=":", alpha=0.45)
    ax2.legend(loc="center right", bbox_to_anchor=(0.99, 0.72), fontsize=7.8, framealpha=0.95)

    plt.subplots_adjust(left=0.06, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "dropout_and_aaf_spread.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


if __name__ == "__main__":
    make_bgzf_binning()
    make_prior_vs_depth()
    make_log_underflow()
    make_qual_vs_gq()
    make_dropout_aaf()
