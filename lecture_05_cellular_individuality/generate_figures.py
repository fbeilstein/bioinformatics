"""
Generate 5 reproducible 200-DPI Matplotlib figures for Lecture 5, Part 3
(Underlying Algorithms: Sparse Matrices, Cell Calling, Robust MAD, Doublets,
Dropout Kinetics & UMI Collisions).

Uses distinct filenames so Part 4's 5 existing PBMC3k images remain untouched.
"""

import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import norm

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
# Figure 1: Sparse Storage Scaling & CSC Array Layout
# =============================================================================
def make_sparse_storage():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.22, "width_ratios": [1.15, 1.0]})

    # Left: RAM vs Droplet Count at 32,738 genes & PBMC3k raw density (0.012039%)
    n_genes = 32738
    pbmc_obs = 737280
    pbmc_nnz = 2905861
    density = pbmc_nnz / (pbmc_obs * n_genes)

    obs_grid = np.logspace(3, 6.1, 200)
    nnz_grid = obs_grid * n_genes * density
    mb = 1024.0 ** 2

    ram_f64 = (8.0 * obs_grid * n_genes) / mb
    ram_f32 = (4.0 * obs_grid * n_genes) / mb
    ram_coo = (12.0 * nnz_grid) / mb
    ram_csc = (8.0 * nnz_grid + 4.0 * (n_genes + 1)) / mb

    ax1.loglog(obs_grid, ram_f64, ":", color="#991b1b", lw=2.0, label="Dense float64 (8 B/entry)")
    ax1.loglog(obs_grid, ram_f32, "-", color="#dc2626", lw=2.5, label="Dense float32 (4 B/entry)")
    ax1.loglog(obs_grid, ram_coo, "--", color="#d97706", lw=2.0, label="COO MatrixMarket (12 B/nnz)")
    ax1.loglog(obs_grid, ram_csc, "-", color="#2563eb", lw=2.6, label="CSC Sparse (8 B/nnz + indptr)")

    ax1.axhline(12 * 1024, color="#64748b", ls="-.", lw=1.4, label="Colab / Laptop RAM (12 GB)")
    ax1.axvline(pbmc_obs, color="#94a3b8", ls=":", lw=1.2)

    pbmc_f32_mb = (4.0 * pbmc_obs * n_genes) / mb
    pbmc_csc_mb = (8.0 * pbmc_nnz + 4.0 * (n_genes + 1)) / mb
    ax1.scatter([pbmc_obs, pbmc_obs], [pbmc_f32_mb, pbmc_csc_mb], color=["#dc2626", "#2563eb"],
                s=65, zorder=6, edgecolors="white", lw=1.1)

    ax1.annotate(
        f"PBMC3k Dense: {pbmc_f32_mb/1024:.1f} GB\n(Instant OOM!)",
        xy=(pbmc_obs, pbmc_f32_mb), xytext=(1.4e4, 1.4e5),
        fontsize=8.0, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.3),
    )
    ax1.annotate(
        f"PBMC3k CSC: {pbmc_csc_mb:.1f} MB\n(4,034× compression!)",
        xy=(pbmc_obs, pbmc_csc_mb), xytext=(1.2e4, 0.7),
        fontsize=8.0, fontweight="bold", color="#1e3a8a",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#eff6ff", edgecolor="#93c5fd"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.3),
    )

    ax1.set_xlim(1e3, 1.3e6)
    ax1.set_ylim(0.12, 8e5)
    ax1.set_xlabel("Number of Droplets n_obs (at 32,738 genes, 99.988% zeros)")
    ax1.set_ylabel("Memory Footprint (MB, log scale)")
    ax1.set_title("Dense vs. Sparse RAM Economics", pad=10)
    ax1.grid(True, which="both", ls=":", alpha=0.35)
    ax1.legend(loc="upper left", fontsize=7.6, framealpha=0.95)

    # Right: Crisp visual diagram of CSC (Compressed Sparse Column) arrays
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.axis("off")
    ax2.set_title("CSC Encoding: Storing Only Non-Zeros (nnz = 5)", pad=10)

    # 4x4 example count matrix on top-left
    M = [
        [3, 0, 0, 1],
        [0, 0, 4, 0],
        [2, 0, 0, 0],
        [0, 0, 5, 0],
    ]
    col_colors = ["#2563eb", "#94a3b8", "#16a34a", "#7c3aed"]
    ax2.text(2.4, 9.35, "Count Matrix (4 cells × 4 genes)", ha="center", va="center", fontsize=9.0, fontweight="bold")
    for j in range(4):
        ax2.text(1.0 + j * 0.92, 8.75, f"G{j}", ha="center", va="center", fontsize=8.4, fontweight="bold", color=col_colors[j])
    for i in range(4):
        ax2.text(0.35, 8.15 - i * 0.72, f"C{i}", ha="center", va="center", fontsize=8.4, fontweight="bold", color="#475569")
        for j in range(4):
            val = M[i][j]
            fc = col_colors[j] if val > 0 else "#f1f5f9"
            tc = "white" if val > 0 else "#94a3b8"
            ec = col_colors[j] if val > 0 else "#cbd5e1"
            box = FancyBboxPatch((0.62 + j * 0.92, 7.86 - i * 0.72), 0.76, 0.56,
                                 boxstyle="round,pad=0.04,rounding_size=0.1", facecolor=fc, edgecolor=ec, lw=1.0)
            ax2.add_patch(box)
            ax2.text(1.0 + j * 0.92, 8.14 - i * 0.72, str(val), ha="center", va="center",
                     fontsize=9.0, fontweight="bold", color=tc)

    # Summary note on top-right
    ax2.text(7.2, 8.35, "3 Contiguous Arrays\n(Zero entries omitted!)", ha="center", va="center",
             fontsize=8.6, fontweight="bold", color="#0f172a")
    ax2.text(7.2, 7.05, "Gene j slice:\ndata[indptr[j] : indptr[j+1]]\n\n• Column sum: O(nnz)\n• Empty gene G1:\n  indptr[1] == indptr[2] == 2",
             ha="center", va="center", fontsize=8.0, color="#334155",
             bbox=dict(boxstyle="round,pad=0.38", facecolor="#f8fafc", edgecolor="#cbd5e1"))

    # CSC Arrays at bottom
    data_arr = [3, 2, 4, 5, 1]
    idx_arr = [0, 2, 1, 3, 0]
    arr_cols = ["#2563eb", "#2563eb", "#16a34a", "#16a34a", "#7c3aed"]
    indptr_arr = [0, 2, 2, 4, 5]

    def draw_array(y0, label, vals, colors, note):
        ax2.text(0.2, y0 + 0.34, label, ha="left", va="center", fontsize=8.8, fontweight="bold", family="monospace")
        for k, v in enumerate(vals):
            fc = colors[k] if isinstance(colors, list) else colors
            box = FancyBboxPatch((2.15 + k * 0.92, y0 + 0.04), 0.78, 0.60,
                                 boxstyle="round,pad=0.03,rounding_size=0.08", facecolor=fc, edgecolor="#334155", lw=0.9)
            ax2.add_patch(box)
            tc = "white" if fc != "#f8fafc" else "#0f172a"
            ax2.text(2.54 + k * 0.92, y0 + 0.34, str(v), ha="center", va="center",
                     fontsize=9.0, fontweight="bold", color=tc, family="monospace")
        ax2.text(6.95, y0 + 0.34, note, ha="left", va="center", fontsize=7.9, color="#475569")

    draw_array(3.85, "data    =", data_arr, arr_cols, "non-zero UMI values")
    draw_array(2.45, "indices =", idx_arr,  arr_cols, "row (cell) index")
    draw_array(1.05, "indptr  =", indptr_arr, "#f8fafc", "col start/end offsets")

    # Index markers 0..4 above data
    for k in range(5):
        ax2.text(2.54 + k * 0.92, 4.72, f"[{k}]", ha="center", va="center", fontsize=7.5, color="#64748b", family="monospace")

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "sparse_storage_scaling.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 2: Barcode Rank ("Knee") Curve & EmptyDrops Ambient Rescue
# =============================================================================
def make_cell_calling_emptydrops():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: PBMC3k-calibrated log-log barcode rank curve matching Part 4's exact landmarks
    ranks_anchor = np.array([1, 10, 69, 500, 1500, 2000, 2500, 2700, 2850, 3000, 5000, 14397, 50000, 158890], dtype=float)
    umis_anchor = np.array([15844, 9200, 5000, 3200, 2195, 1775, 1102, 548, 220, 89, 28, 10, 3, 1], dtype=float)
    log_r = np.linspace(0, math.log10(158890), 500)
    log_u = np.interp(log_r, np.log10(ranks_anchor), np.log10(umis_anchor))
    ranks = 10 ** log_r
    umis = 10 ** log_u

    ax1.axvspan(1, 2000, color="#2563eb", alpha=0.08, label="1. Cellular plateau")
    ax1.axvspan(2000, 3000, color="#f59e0b", alpha=0.14, label="2. Cliff / knee zone")
    ax1.axvspan(3000, 2e5, color="#94a3b8", alpha=0.10, label="3. Ambient soup tail")

    ax1.loglog(ranks, umis, "-", color="#0f172a", lw=2.5, zorder=5)
    ax1.axhline(568.9, color="#7c3aed", ls="--", lw=1.6, label="Cell Ranger v2 (m/10 = 569)")
    ax1.axhline(500.0, color="#dc2626", ls=":", lw=1.7, label="Hard cutoff (UMI > 500)")

    ax1.scatter([2700], [548], color="#dc2626", s=55, zorder=7, edgecolors="white", lw=1.0)
    ax1.annotate(
        "Rank 2,700: 548 UMIs\n(90.4% of UMIs in <0.4%\nof droplets!)",
        xy=(2700, 548), xytext=(14, 14),
        fontsize=7.9, fontweight="bold", color="#1e293b",
        bbox=dict(boxstyle="round,pad=0.30", facecolor="#f8fafc", edgecolor="#64748b"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
    )
    ax1.annotate(
        "Low-RNA cells\nhidden here?",
        xy=(2900, 175), xytext=(7500, 1200),
        fontsize=8.0, fontweight="bold", color="#15803d",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.2),
    )

    ax1.set_xlim(1, 2e5)
    ax1.set_ylim(0.8, 3e4)
    ax1.set_xlabel("Barcode Rank (descending UMI, log scale)")
    ax1.set_ylabel("Total UMI Count (log scale)")
    ax1.set_title("Barcode Rank ('Knee') Curve on PBMC3k", pad=10)
    ax1.grid(True, which="both", ls=":", alpha=0.32)
    ax1.legend(loc="lower left", fontsize=7.3, framealpha=0.95)

    # Right: EmptyDrops Dirichlet-Multinomial rescue scatter
    rng = np.random.default_rng(501)
    # Ambient empty droplets (low-mid UMI, non-significant deviation from ambient soup)
    u_amb = 10 ** rng.uniform(1.3, 2.65, 220)
    sig_amb = rng.exponential(0.38, 220)
    # Normal high-UMI cells (well above knee AND distinct from ambient profile)
    u_cell = 10 ** rng.uniform(2.72, 3.85, 140)
    sig_cell = 3.2 + rng.uniform(1.5, 6.2, 140) + 1.1 * (np.log10(u_cell) - 2.7)
    # Genuine low-RNA cells (e.g., platelets / small lymphocytes, 150-500 UMI, distinct profile!)
    u_resc = 10 ** rng.uniform(2.15, 2.69, 32)
    sig_resc = rng.uniform(2.4, 6.2, 32)

    ax2.axvspan(100, 500, ymin=0.20, ymax=1.0, color="#16a34a", alpha=0.11)
    ax2.scatter(u_amb, sig_amb, s=16, color="#94a3b8", alpha=0.65, label="Empty soup droplets (6.6% MT-)")
    ax2.scatter(u_cell, sig_cell, s=18, color="#2563eb", alpha=0.75, label="Called cells (> 500 UMIs)")
    ax2.scatter(u_resc, sig_resc, s=28, color="#16a34a", edgecolors="white", lw=0.5, zorder=6,
                label="Rescued low-RNA cells (EmptyDrops)")

    ax2.axvline(500, color="#dc2626", ls=":", lw=1.6, label="Hard cutoff (500 UMIs)")
    ax2.axhline(2.0, color="#15803d", ls="--", lw=1.6, label="EmptyDrops FDR = 0.01")

    ax2.annotate(
        "Rescued by EmptyDrops!\nProfile deviates from\nambient soup p̂_amb",
        xy=(260, 4.4), xytext=(19, 7.6),
        fontsize=7.9, fontweight="bold", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.2),
    )

    ax2.set_xscale("log")
    ax2.set_xlim(15, 1.2e4)
    ax2.set_ylim(-0.2, 10.4)
    ax2.set_xlabel("Droplet Total UMI Count (log scale)")
    ax2.set_ylabel(r"EmptyDrops Evidence  $-\log_{10}(\mathrm{FDR})$")
    ax2.set_title("Dirichlet-Multinomial Ambient Test (EmptyDrops)", pad=10)
    ax2.grid(True, which="both", ls=":", alpha=0.35)
    ax2.legend(loc="lower right", fontsize=7.3, framealpha=0.95)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "cell_calling_emptydrops.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 3: Normal Consistency Factor (1.4826) & Breakdown Curve
# =============================================================================
def make_mad_breakdown():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: Standard Normal N(0,1) deriving k = 1 / Phi^-1(0.75) = 1.4826
    z = np.linspace(-3.6, 3.6, 400)
    pdf = norm.pdf(z)
    c = float(norm.ppf(0.75))  # 0.67448975...

    ax1.plot(z, pdf, "-", color="#1e293b", lw=2.3)
    mask_mid = np.abs(z) <= c
    ax1.fill_between(z[mask_mid], 0, pdf[mask_mid], color="#2563eb", alpha=0.30, label=f"Central 50% mass: |Z| ≤ {c:.4f}")
    ax1.fill_between(z[z < -c], 0, pdf[z < -c], color="#94a3b8", alpha=0.22, label="Outer tails: 25% + 25%")
    ax1.fill_between(z[z > c], 0, pdf[z > c], color="#94a3b8", alpha=0.22)
    ax1.axvline(-c, color="#2563eb", ls="--", lw=1.4)
    ax1.axvline(c, color="#2563eb", ls="--", lw=1.4)
    ax1.axvline(1.0, color="#dc2626", ls=":", lw=1.5, label="1 SD (σ = 1.00)")

    ax1.text(0, 0.16, "50%\nmass", ha="center", va="center", fontsize=9.2, fontweight="bold", color="#1e3a8a")
    ax1.text(-1.55, 0.05, "25%", ha="center", va="center", fontsize=8.4, color="#475569", fontweight="bold")
    ax1.text(1.55, 0.05, "25%", ha="center", va="center", fontsize=8.4, color="#475569", fontweight="bold")

    ax1.annotate(
        r"$\Phi^{-1}(0.75) = 0.6745$" + "\n" + r"$\Rightarrow k = \frac{1}{0.6745} \approx 1.4826$",
        xy=(c, norm.pdf(c)), xytext=(1.35, 0.30),
        fontsize=8.3, fontweight="bold", color="#1e3a8a",
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#eff6ff", edgecolor="#93c5fd"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.3),
    )

    ax1.set_xlim(-3.5, 3.5)
    ax1.set_ylim(0, 0.45)
    ax1.set_xlabel(r"Standard Normal $Z \sim \mathcal{N}(0, 1)$")
    ax1.set_ylabel("Probability Density")
    ax1.set_title("Deriving the Normal Consistency Factor k ≈ 1.4826", pad=10)
    ax1.grid(True, ls=":", alpha=0.4)
    ax1.legend(loc="upper left", fontsize=7.6, framealpha=0.95)

    # Right: Contamination breakdown curve (Mean + 3 SD vs Median + 3*1.4826*MAD)
    rng = np.random.default_rng(2026)
    n_tot = 2000
    clean_base = np.clip(rng.normal(2.03, 0.78, n_tot), 0.1, 4.5)
    apop_pool = np.clip(rng.normal(22.0, 5.5, n_tot), 10.0, 45.0)
    fracs = np.linspace(0.0, 0.56, 57)

    thr_sd = []
    thr_mad = []
    k_norm = 1.0 / norm.ppf(0.75)
    for f in fracs:
        n_bad = int(round(f * n_tot))
        vals = np.concatenate([clean_base[:n_tot - n_bad], apop_pool[:n_bad]])
        thr_sd.append(np.mean(vals) + 3.0 * np.std(vals, ddof=1))
        med = np.median(vals)
        mad_s = np.median(np.abs(vals - med)) * k_norm
        thr_mad.append(med + 3.0 * mad_s)

    ax2.plot(fracs * 100, thr_sd, "-", color="#dc2626", lw=2.5, label="Mean + 3 SD (0% breakdown!)")
    ax2.plot(fracs * 100, thr_mad, "-", color="#2563eb", lw=2.6, label="Median + 3 MAD (50% breakdown)")
    ax2.axvline(50.0, color="#7c3aed", ls="--", lw=1.6, label="Theoretical 50% breakdown cliff")

    # Mark PBMC3k observed point (~3.5% flagged apoptotic cells: Mean+3SD=5.94%, Med+3MAD=4.43%)
    ax2.scatter([3.5, 3.5], [5.94, 4.43], color=["#dc2626", "#2563eb"], s=58, zorder=7, edgecolors="white", lw=1.0)
    ax2.annotate(
        "PBMC3k (~3.5% dying):\nMean+3SD = 5.94% (masks!)\nMed+3MAD = 4.43% (clean)",
        xy=(3.5, 5.94), xytext=(15.5, 16.0),
        fontsize=7.9, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.30", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
    )

    ax2.set_xlim(0, 56)
    ax2.set_ylim(0, 55)
    ax2.set_xlabel("Fraction of Apoptotic Outlier Cells (MT ≈ 22%) (%)")
    ax2.set_ylabel("Computed Upper MT% Cutoff (%)")
    ax2.set_title("Finite-Sample Breakdown & Outlier Masking", pad=10)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.legend(loc="upper left", fontsize=7.6, framealpha=0.96)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "mad_breakdown_and_normal_scale.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 4: Poisson Microfluidic Loading & Scrublet 2D PCA Doublet Bridges
# =============================================================================
def make_poisson_scrublet():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: Poisson encapsulation curves vs loading rate lambda
    lam = np.linspace(0.002, 0.35, 300)
    p0 = np.exp(-lam)
    p1 = lam * np.exp(-lam)
    p_ge2 = 1.0 - p0 - p1
    cond_dbl = p_ge2 / (1.0 - p0)

    ax1.plot(lam, p0 * 100, "--", color="#64748b", lw=2.0, label=r"Empty droplets $P(X=0) = e^{-\lambda}$")
    ax1.plot(lam, p1 * 100, "-", color="#2563eb", lw=2.4, label=r"Singlets $P(X=1) = \lambda e^{-\lambda}$")
    ax1.plot(lam, cond_dbl * 100, "-", color="#dc2626", lw=2.6, label=r"Multiplet rate $P(X\geq 2 \mid X\geq 1)$")
    ax1.plot(lam, (lam / 2.0) * 100, ":", color="#f59e0b", lw=1.8, label=r"Small-$\lambda$ approx. $\lambda / 2$")

    lam0 = 0.05
    c0 = (1 - math.exp(-lam0) - lam0 * math.exp(-lam0)) / (1 - math.exp(-lam0)) * 100
    ax1.axvline(lam0, color="#16a34a", ls="-.", lw=1.5)
    ax1.scatter([lam0, lam0], [math.exp(-lam0) * 100, c0], color=["#64748b", "#dc2626"], s=58, zorder=7, edgecolors="white", lw=1.0)
    ax1.annotate(
        f"10x Sweet Spot (λ = 0.05):\n• 95.1% empty droplets\n• {c0:.2f}% doublets among cells",
        xy=(lam0, c0), xytext=(0.075, 29.0),
        fontsize=7.9, fontweight="bold", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.30", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.2),
    )

    ax1.set_xlim(0, 0.35)
    ax1.set_ylim(0, 104)
    ax1.set_xlabel(r"Mean Cells per Droplet  $\lambda$")
    ax1.set_ylabel("Percentage (%)")
    ax1.set_title("Poisson Loading Trade-Off: Empty vs. Doublets", pad=10)
    ax1.grid(True, ls=":", alpha=0.4)
    ax1.legend(loc="center right", bbox_to_anchor=(0.99, 0.62), fontsize=7.3, framealpha=0.96)

    # Right: Scrublet 2D PCA geometry (3 immune clusters + heterotypic bridges vs homotypic)
    rng = np.random.default_rng(77)
    centers = {
        "T-cells":   np.array([-3.0, -1.8]),
        "B-cells":   np.array([3.0,  -1.8]),
        "Monocytes": np.array([0.0,   3.0]),
    }
    c_cols = {"T-cells": "#2563eb", "B-cells": "#10b981", "Monocytes": "#7c3aed"}
    pts = []
    labs = []
    for k, mu in centers.items():
        cl = rng.normal(mu, 0.52, (95, 2))
        pts.append(cl)
        labs.extend([k] * len(cl))
        ax2.scatter(cl[:, 0], cl[:, 1], s=18, color=c_cols[k], alpha=0.72, label=k)
    pts = np.vstack(pts)
    labs = np.array(labs)

    # Simulate pairwise doublets (midpoint in normalized PCA space + small noise)
    n_sim = 130
    i1 = rng.integers(0, len(pts), n_sim)
    i2 = rng.integers(0, len(pts), n_sim)
    sim_pts = 0.5 * (pts[i1] + pts[i2]) + rng.normal(0, 0.25, (n_sim, 2))
    is_het = labs[i1] != labs[i2]

    ax2.scatter(sim_pts[~is_het, 0], sim_pts[~is_het, 1], s=22, marker="^", color="#94a3b8",
                alpha=0.55, label="Homotypic sim (hidden)")
    ax2.scatter(sim_pts[is_het, 0], sim_pts[is_het, 1], s=30, marker="D", color="#f97316",
                edgecolors="#9a3412", lw=0.5, alpha=0.85, zorder=6, label="Heterotypic sim (bridges!)")

    ax2.annotate(
        "Heterotypic bridges:\nhigh Scrublet k-NN score!",
        xy=(1.55, 0.55), xytext=(1.35, 3.25),
        fontsize=7.9, fontweight="bold", color="#9a3412",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#fff7ed", edgecolor="#fdba74"),
        arrowprops=dict(arrowstyle="->", color="#ea580c", lw=1.2),
    )

    ax2.set_xlim(-4.8, 4.9)
    ax2.set_ylim(-4.8, 4.6)
    ax2.set_xlabel("PC 1 (Joint Transcriptome PCA)")
    ax2.set_ylabel("PC 2 (Joint Transcriptome PCA)")
    ax2.set_title("Scrublet Geometry: Heterotypic Bridges vs. Homotypic", pad=10)
    ax2.grid(True, ls=":", alpha=0.35)
    ax2.legend(loc="lower center", ncol=2, fontsize=7.2, framealpha=0.96)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "poisson_loading_and_scrublet.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 5: Dropout Kinetics (PBMC3k vs Poisson/NB) & UMI Occupancy Saturation
# =============================================================================
def make_umi_and_dropout():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: Zero fraction vs log10(mean mu) with Part 4's exact 5 PBMC3k bins
    mu_g = np.logspace(-3, 1.6, 300)
    z_pois = np.exp(-mu_g)
    z_nb05 = (1.0 + 0.5 * mu_g) ** (-1.0 / 0.5)
    z_nb20 = (1.0 + 2.0 * mu_g) ** (-1.0 / 2.0)

    ax1.semilogx(mu_g, z_pois, "--", color="#dc2626", lw=2.2, label=r"Poisson capture: $e^{-\mu}$")
    ax1.semilogx(mu_g, z_nb05, "-", color="#2563eb", lw=2.2, label=r"NB ($\alpha = 0.5$): $(1+\alpha\mu)^{-1/\alpha}$")
    ax1.semilogx(mu_g, z_nb20, ":", color="#7c3aed", lw=2.2, label=r"NB ($\alpha = 2.0$, cell-type burst)")

    # Part 4's exact 5 PBMC3k observed zero fraction bins
    pbmc_mu = np.array([0.0032, 0.032, 0.32, 2.3, 16.0])
    pbmc_zero = np.array([0.9969, 0.9716, 0.8462, 0.2556, 0.0069])
    ax1.scatter(pbmc_mu, pbmc_zero, s=72, color="#f59e0b", edgecolors="#0f172a", lw=1.0, zorder=7,
                label="PBMC3k observed bins (Part 4)")

    ax1.annotate(
        "Low μ (< 0.1):\nMatches Poisson e^(−μ)\n(capture shot noise)",
        xy=(0.032, 0.9716), xytext=(0.0014, 0.58),
        fontsize=7.9, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
    )
    ax1.annotate(
        "High μ (1–10): 25.6% zeros\nvs. 10.1% Poisson\n(biological heterogeneity!)",
        xy=(2.3, 0.2556), xytext=(0.0014, 0.18),
        fontsize=7.9, fontweight="bold", color="#1e3a8a",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#eff6ff", edgecolor="#93c5fd"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.2),
    )

    ax1.set_xlim(1e-3, 35)
    ax1.set_ylim(-0.04, 1.07)
    ax1.set_xlabel(r"Mean UMI Count per Cell  $\mu$ (log scale)")
    ax1.set_ylabel("Fraction of Exact Zeros  P(X = 0)")
    ax1.set_title("Technical Dropout vs. Biological Overdispersion", pad=10)
    ax1.grid(True, which="both", ls=":", alpha=0.35)
    ax1.legend(loc="upper right", fontsize=7.4, framealpha=0.95)

    # Right: UMI Occupancy Saturation E(k) = N * (1 - exp(-k/N)) for L = 8, 10, 12
    # (No annotation at k=50,000, L=10 so Homework Task 4 stays unspoiled!)
    k_grid = np.logspace(2, 6.7, 350)
    ax2.loglog(k_grid, k_grid, "--", color="#64748b", lw=1.6, label="Ideal (no collisions): E = k")

    configs = [
        (8,  "#dc2626", "L = 8 nt  (N = 65,536)"),
        (10, "#2563eb", "L = 10 nt (N = 1.05M, 10x v2)"),
        (12, "#16a34a", "L = 12 nt (N = 16.8M, 10x v3)"),
    ]
    for L, col, lbl in configs:
        N = float(4 ** L)
        E = N * (-np.expm1(-k_grid / N))
        ax2.loglog(k_grid, E, "-", color=col, lw=2.4, label=lbl)
        # Mark 50% birthday collision onset k_bday = sqrt(2 ln 2 * N)
        k_bd = math.sqrt(2.0 * math.log(2.0) * N)
        E_bd = N * (-math.expm1(-k_bd / N))
        ax2.scatter([k_bd], [E_bd], color=col, s=52, marker="D", edgecolors="white", lw=0.8, zorder=7)

    # Dummy legend entry for birthday diamond
    ax2.scatter([], [], color="#334155", s=45, marker="D", label=r"50% birthday onset ($k \approx 1.18\sqrt{N}$)")

    ax2.annotate(
        "8-nt saturates early:\ncaps at 65,536 UMIs!",
        xy=(4e5, 65536), xytext=(1.5e4, 4.5e3),
        fontsize=7.9, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
    )
    ax2.annotate(
        "10x v2 → v3 (+2 nt):\n16× capacity expansion",
        xy=(1.5e6, 1.5e6 * (1 - 0.04)), xytext=(3e2, 4.5e5),
        fontsize=7.9, fontweight="bold", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.2),
    )

    ax2.set_xlim(1e2, 5e6)
    ax2.set_ylim(1e2, 8e6)
    ax2.set_xlabel("True Captured Molecules k (log scale)")
    ax2.set_ylabel("Expected Distinct UMIs Observed E (log scale)")
    ax2.set_title("UMI Occupancy Saturation & Barcode Capacity Scaling", pad=10)
    ax2.grid(True, which="both", ls=":", alpha=0.35)
    ax2.legend(loc="lower right", fontsize=7.3, framealpha=0.95)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "umi_occupancy_and_dropout.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


if __name__ == "__main__":
    make_sparse_storage()
    make_cell_calling_emptydrops()
    make_mad_breakdown()
    make_poisson_scrublet()
    make_umi_and_dropout()
