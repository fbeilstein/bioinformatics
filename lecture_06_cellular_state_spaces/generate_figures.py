"""
Generate 5 reproducible 200-DPI Matplotlib figures for Lecture 6, Part 3
(Underlying Algorithms: HVGs, Spectral PCA, Curse of Dimensionality,
Fuzzy k-NN, Leiden Refinement, UMAP Cross-Entropy, Diffusion Pseudotime & Waddington KDE).

Uses distinct filenames so Part 4's 6 existing PBMC3k images remain untouched.
"""

import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
from scipy.stats import gaussian_kde

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
# Figure 1: Mean-Binned HVG Selection & Z-Score Clipping (max_value=10)
# =============================================================================
def make_hvg_and_zscore():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: Mean log-expression vs Mean-binned normalized dispersion z_d
    rng = np.random.default_rng(601)
    n_hk = 230
    mu_hk = rng.uniform(0.08, 3.3, n_hk)
    # Housekeeping genes follow baseline mean-dispersion trend -> z_d ~ N(0, 0.45)
    zd_hk = rng.normal(-0.08, 0.38, n_hk)
    # High-expression abundant housekeeping genes have huge raw variance = mu * exp(...)
    raw_var_hk = (np.expm1(mu_hk) ** 1.35) * np.exp(zd_hk * 0.25)

    n_mk = 45
    mu_mk = rng.uniform(0.25, 2.75, n_mk)
    zd_mk = rng.uniform(0.75, 3.6, n_mk)
    raw_var_mk = (np.expm1(mu_mk) ** 1.35) * np.exp(zd_mk * 0.35)

    all_mu = np.concatenate([mu_hk, mu_mk])
    all_zd = np.concatenate([zd_hk, zd_mk])
    all_rv = np.concatenate([raw_var_hk, raw_var_mk])
    top_rv_idx = np.argsort(all_rv)[-18:]

    ax1.scatter(mu_hk, zd_hk, s=18, color="#94a3b8", alpha=0.68, label="Non-HVG housekeeping genes")
    ax1.scatter(mu_mk, zd_mk, s=28, color="#2563eb", edgecolors="white", lw=0.4, zorder=5,
                label="Mean-binned HVGs (z_d > 0.5)")
    ax1.scatter(all_mu[top_rv_idx], all_zd[top_rv_idx], s=72, facecolors="none", edgecolors="#dc2626",
                lw=1.5, zorder=6, label="Top raw-variance picks (bias!)")

    ax1.axhline(0.5, color="#2563eb", ls="--", lw=1.5)
    ax1.axvline(3.0, color="#64748b", ls=":", lw=1.2)

    ax1.annotate(
        "Raw variance picks\nabundant housekeeping!\n(high mean, low z_d)",
        xy=(3.10, -0.05), xytext=(1.35, -1.35),
        fontsize=7.9, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
    )
    ax1.annotate(
        "Binned dispersion z_d > 0.5:\nisolates true lineage markers!",
        xy=(1.2, 2.35), xytext=(0.12, 3.25),
        fontsize=7.9, fontweight="bold", color="#1e3a8a",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#eff6ff", edgecolor="#93c5fd"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.2),
    )

    ax1.set_xlim(0.0, 3.5)
    ax1.set_ylim(-1.65, 4.15)
    ax1.set_xlabel("Mean Log-Expression  μ_g")
    ax1.set_ylabel(r"Normalized Dispersion  $z_d = (d_g - \bar{d}_{\mathrm{bin}})/\sigma_{d,\mathrm{bin}}$")
    ax1.set_title("Why Raw Variance Fails vs. Mean-Binned HVGs", pad=10)
    ax1.grid(True, ls=":", alpha=0.4)
    ax1.legend(loc="lower left", fontsize=7.3, framealpha=0.95)

    # Right: Single-cell +40 sigma burst hijacking PC1 unclipped vs max_value=10 clipped
    rng2 = np.random.default_rng(602)
    n_c, n_g = 240, 30
    # True 2-lineage biological signal spread across first 12 marker genes
    lin = np.concatenate([-np.ones(n_c // 2), np.ones(n_c // 2)])
    X = rng2.normal(0.0, 1.0, (n_c, n_g))
    for g in range(12):
        X[:, g] += 0.95 * lin
    # Z-score genes cleanly
    X = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
    # Inject 1 extreme single-cell spike (+40 SD) on gene 29 (non-marker artifact)
    X_raw = X.copy()
    X_raw[0, 29] = 40.0
    X_clip = np.clip(X_raw, -10.0, 10.0)

    def pc1_loadings(M):
        Mc = M - M.mean(axis=0)
        _, _, vt = np.linalg.svd(Mc, full_matrices=False)
        return np.abs(vt[0])

    load_raw = pc1_loadings(X_raw)
    load_clip = pc1_loadings(X_clip)

    # Plot top 5 lineage markers + the 1 outlier burst gene
    idx_show = [0, 1, 2, 3, 4, 29]
    lbls = ["Marker 1", "Marker 2", "Marker 3", "Marker 4", "Marker 5", "1-Cell\n+40σ Spike"]
    x_pos = np.arange(len(idx_show))
    w = 0.36
    b1 = ax2.bar(x_pos - w / 2, load_raw[idx_show], width=w, color="#dc2626", alpha=0.88, label="Unclipped Z-scores")
    b2 = ax2.bar(x_pos + w / 2, load_clip[idx_show], width=w, color="#16a34a", alpha=0.90, label="Clipped (max_value=10)")

    ax2.annotate(
        f"1 outlier cell hijacks\nPC1 loading = {load_raw[29]:.2f}!",
        xy=(5 - w / 2, load_raw[29]), xytext=(2.05, 0.80),
        fontsize=8.0, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
    )
    ax2.annotate(
        f"Clipped to {load_clip[29]:.2f} ✓\n(lineage markers dominate PC1)",
        xy=(5 + w / 2, load_clip[29]), xytext=(2.25, 0.50),
        fontsize=7.9, fontweight="bold", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.2),
    )

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(lbls, fontsize=8.3)
    ax2.set_ylim(0, 1.08)
    ax2.set_ylabel("|PC1 Loading Weight|  (|v_1,g|)")
    ax2.set_title("Single-Cell Leverage: Why sc.pp.scale Clips at ±10σ", pad=10)
    ax2.grid(axis="y", ls=":", alpha=0.45)
    ax2.legend(loc="upper left", fontsize=7.8, framealpha=0.95)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "hvg_and_zscore_clipping.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 2: Spectral PCA Denoising (Marchenko-Pastur) & Curse of Dimensionality
# =============================================================================
def make_svd_and_curse():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: Exact SVD spectrum of 6 correlated biological modules + isotropic noise
    rng = np.random.default_rng(603)
    n, p = 400, 600
    Z_sig = rng.normal(0.0, 1.0, (n, 6))
    # Distinct module strengths
    W_sig = rng.normal(0.0, 1.0, (6, p)) * np.array([0.48, 0.40, 0.33, 0.28, 0.23, 0.19])[:, None]
    X_tot = Z_sig @ W_sig + rng.normal(0.0, 1.0, (n, p))
    X_tot -= X_tot.mean(axis=0)
    _, S, _ = np.linalg.svd(X_tot, full_matrices=False)
    eigs = (S ** 2) / (n - 1)
    # Theoretical Marchenko-Pastur upper edge for sigma^2 = 1: (1 + sqrt(p/n))^2
    lam_mp = (1.0 + math.sqrt(p / n)) ** 2

    k_show = 35
    pcs = np.arange(1, k_show + 1)
    ax1.semilogy(pcs[:6], eigs[:6], "o", ms=6.2, color="#2563eb", zorder=6, label="Correlated biological modules (signal)")
    ax1.semilogy(pcs[6:], eigs[6:k_show], "o", ms=4.5, color="#94a3b8", zorder=5, label="Uncorrelated noise bulk")
    ax1.semilogy(pcs, eigs[:k_show], "-", color="#475569", lw=1.2, alpha=0.6, zorder=4)
    ax1.axhline(lam_mp, color="#dc2626", ls="--", lw=1.8, label=rf"Marchenko–Pastur edge $\lambda_+ = {lam_mp:.1f}$")
    ax1.axhspan(0.8, lam_mp, color="#94a3b8", alpha=0.14)

    ax1.annotate(
        "Truncate above noise bulk:\nretains co-regulated programs,\ndiscards independent noise!",
        xy=(6, eigs[5]), xytext=(11.5, 42.0),
        fontsize=8.0, fontweight="bold", color="#1e3a8a",
        bbox=dict(boxstyle="round,pad=0.30", facecolor="#eff6ff", edgecolor="#93c5fd"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.2),
    )

    ax1.set_xlim(0.5, k_show + 0.5)
    ax1.set_ylim(1.2, 220)
    ax1.set_xlabel("Principal Component Rank i")
    ax1.set_ylabel(r"Eigenvalue  $\lambda_i = S_i^2 / (n - 1)$ (log scale)")
    ax1.set_title("PCA as a Spectral Denoiser (Signal vs. Noise Bulk)", pad=10)
    ax1.grid(True, which="both", ls=":", alpha=0.35)
    ax1.legend(loc="upper right", fontsize=7.4, framealpha=0.95)

    # Right: Curse of dimensionality distance concentration (p = 2 .. 2000)
    rng2 = np.random.default_rng(604)
    dims = np.array([2, 5, 10, 20, 40, 100, 250, 500, 1000, 1838, 4000])
    rel_contrast = []
    dist_cv = []
    n_pts = 220
    for d in dims:
        pts = rng2.normal(0.0, 1.0, (n_pts, d))
        # Distances from first 10 reference points to all others
        diffs = pts[:10, None, :] - pts[None, :, :]
        dists = np.sqrt(np.sum(diffs ** 2, axis=-1))
        # Exclude exact zero self-distance
        vals = dists[dists > 1e-9]
        rc = np.mean([(np.max(row[row > 1e-9]) - np.min(row[row > 1e-9])) / np.min(row[row > 1e-9]) for row in dists])
        rel_contrast.append(rc)
        dist_cv.append(np.std(vals) / np.mean(vals))

    ax2.loglog(dims, rel_contrast, "-o", ms=5.0, color="#dc2626", lw=2.3,
               label=r"Relative contrast $(\mathrm{dist}_{\max} - \mathrm{dist}_{\min})/\mathrm{dist}_{\min}$")
    ax2.loglog(dims, dist_cv, "--s", ms=4.5, color="#2563eb", lw=2.1,
               label=r"Pairwise distance CV  $\mathrm{SD}(d)/\bar{d} \approx 1/\sqrt{2p}$")

    ax2.axvline(40, color="#16a34a", ls="-.", lw=1.6)
    ax2.axvline(1838, color="#7c3aed", ls=":", lw=1.6)

    ax2.annotate(
        "40D PCA space\n(meaningful local k-NN!)",
        xy=(40, rel_contrast[4]), xytext=(3.2, 0.16),
        fontsize=7.9, fontweight="bold", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#f0fdf4", edgecolor="#86efac"),
        arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.2),
    )
    ax2.annotate(
        "1,838 HVGs:\ndistances concentrate\n(CV ≈ 1.6%!)",
        xy=(1838, dist_cv[9]), xytext=(160, 0.012),
        fontsize=7.9, fontweight="bold", color="#5b21b6",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#f5f3ff", edgecolor="#c4b5fd"),
        arrowprops=dict(arrowstyle="->", color="#7c3aed", lw=1.2),
    )

    ax2.set_xlim(1.6, 5500)
    ax2.set_ylim(0.008, 35)
    ax2.set_xlabel("Ambient Dimension p (log scale)")
    ax2.set_ylabel("Distance Spread Metric (log scale)")
    ax2.set_title("Curse of Dimensionality: Distance Concentration", pad=10)
    ax2.grid(True, which="both", ls=":", alpha=0.35)
    ax2.legend(loc="upper right", fontsize=7.3, framealpha=0.95)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "svd_denoising_and_curse.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 3: Curved Manifold Adaptive k-NN & Louvain Flaw vs. Leiden Refinement
# =============================================================================
def make_fuzzy_knn_and_leiden():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.22})

    # Left: Curved C-manifold (dense on upper branch, sparse on lower branch)
    rng = np.random.default_rng(605)
    # Non-uniform density along arc theta in [0.18*pi, 1.82*pi]
    t_dense = np.linspace(0.22 * math.pi, 0.95 * math.pi, 36)
    t_sparse = np.linspace(1.05 * math.pi, 1.78 * math.pi, 14)
    thetas = np.concatenate([t_dense, t_sparse])
    xs = 2.4 * np.cos(thetas) + rng.normal(0, 0.08, len(thetas))
    ys = 2.4 * np.sin(thetas) + rng.normal(0, 0.08, len(thetas))
    pts = np.column_stack([xs, ys])

    # Build k=4 nearest neighbor graph
    for i in range(len(pts)):
        ds = np.linalg.norm(pts - pts[i], axis=1)
        nbrs = np.argsort(ds)[1:5]
        for j in nbrs:
            ax1.plot([pts[i, 0], pts[j, 0]], [pts[i, 1], pts[j, 1]], color="#2563eb", lw=0.9, alpha=0.45, zorder=2)

    # Highlight geodesic path along manifold vs shortcut Euclidean chord between endpoints
    i_start, i_end = 0, len(pts) - 1
    ax1.plot(pts[:, 0], pts[:, 1], "-", color="#1d4ed8", lw=2.3, alpha=0.85, zorder=3, label="Geodesic k-NN manifold path")
    ax1.plot([pts[i_start, 0], pts[i_end, 0]], [pts[i_start, 1], pts[i_end, 1]], "--", color="#dc2626", lw=2.4, zorder=4,
             label="Euclidean chord (cuts empty space!)")

    # Show two adaptive UMAP kernel circles: dense region (tight sigma_i) vs sparse region (wide sigma_i)
    i_d, i_s = 12, 43
    d_d = np.sort(np.linalg.norm(pts - pts[i_d], axis=1))[4]
    d_s = np.sort(np.linalg.norm(pts - pts[i_s], axis=1))[4]
    ax1.add_patch(Circle(pts[i_d], d_d, facecolor="#16a34a", edgecolor="#15803d", alpha=0.18, lw=1.6, zorder=1))
    ax1.add_patch(Circle(pts[i_s], d_s, facecolor="#f59e0b", edgecolor="#b45309", alpha=0.18, lw=1.6, zorder=1))

    ax1.scatter(pts[:, 0], pts[:, 1], s=26, color="#0f172a", edgecolors="white", lw=0.5, zorder=5)
    ax1.scatter([pts[i_start, 0], pts[i_end, 0]], [pts[i_start, 1], pts[i_end, 1]], s=70, color="#dc2626", edgecolors="white", lw=1.0, zorder=6)

    ax1.annotate("Dense region:\nsmall σ_i", xy=(pts[i_d, 0], pts[i_d, 1] + d_d), xytext=(-3.25, 2.35),
                 fontsize=7.9, fontweight="bold", color="#14532d",
                 bbox=dict(boxstyle="round,pad=0.24", facecolor="#f0fdf4", edgecolor="#86efac"),
                 arrowprops=dict(arrowstyle="->", color="#16a34a", lw=1.1))
    ax1.annotate("Sparse region:\nlarge σ_i adapts!", xy=(pts[i_s, 0], pts[i_s, 1] + d_s),
                 xytext=(-1.15, -1.05), zorder=8,
                 fontsize=7.8, fontweight="bold", color="#92400e",
                 bbox=dict(boxstyle="round,pad=0.24", facecolor="#fffbeb", edgecolor="#f59e0b"),
                 arrowprops=dict(arrowstyle="->", color="#d97706", lw=1.1))

    ax1.set_xlim(-3.45, 3.45)
    ax1.set_ylim(-3.60, 3.45)
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_xlabel("State Axis 1")
    ax1.set_ylabel("State Axis 2")
    ax1.set_title("Local Adaptive k-NN vs. Global Euclidean Distance", pad=10)
    ax1.grid(True, ls=":", alpha=0.35)
    ax1.legend(loc="center left", bbox_to_anchor=(0.02, 0.50), fontsize=7.2, framealpha=0.95)

    # Right: Traag et al. (2019) Louvain Disconnection Flaw vs. Leiden Refinement Diagram
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.axis("off")
    ax2.set_title("Louvain's Disconnection Flaw vs. Leiden Refinement", pad=10)

    def draw_panel(y_off, title, sub, col_title, border_col, bg_col, a1_col, a2_col, bridge_col, is_leiden):
        card = FancyBboxPatch((0.25, y_off), 9.5, 4.35, boxstyle="round,pad=0.08,rounding_size=0.25",
                              facecolor=bg_col, edgecolor=border_col, lw=1.4)
        ax2.add_patch(card)
        ax2.text(0.55, y_off + 3.88, title, fontsize=8.6, fontweight="bold", color=col_title, va="center")
        ax2.text(0.55, y_off + 0.48, sub, fontsize=7.3, color="#334155", va="center")

        # Left clique A1, right clique A2, central bridge v* pulled into bottom-center B
        a1 = [(1.6, y_off + 2.85), (1.1, y_off + 1.85), (2.1, y_off + 1.85)]
        a2 = [(5.4, y_off + 2.85), (4.9, y_off + 1.85), (5.9, y_off + 1.85)]
        vb = (3.5, y_off + 2.35)
        b_cl = [(7.5, y_off + 2.85), (7.5, y_off + 1.70), (8.6, y_off + 2.28)]

        def draw_clique(nodes, col):
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    ax2.plot([nodes[i][0], nodes[j][0]], [nodes[i][1], nodes[j][1]], color="#475569", lw=1.4, zorder=2)
            for nx, ny in nodes:
                ax2.scatter([nx], [ny], s=140, color=col, edgecolors="white", lw=1.2, zorder=5)

        # Bridge edges
        ax2.plot([a1[2][0], vb[0]], [a1[2][1], vb[1]], color="#64748b", lw=1.3, ls="--", zorder=2)
        ax2.plot([a2[1][0], vb[0]], [a2[1][1], vb[1]], color="#64748b", lw=1.3, ls="--", zorder=2)
        ax2.plot([vb[0], b_cl[0][0]], [vb[1], b_cl[0][1]], color="#ea580c", lw=2.0, zorder=2)
        ax2.plot([vb[0], b_cl[1][0]], [vb[1], b_cl[1][1]], color="#ea580c", lw=2.0, zorder=2)

        draw_clique(a1, a1_col)
        draw_clique(a2, a2_col)
        draw_clique(b_cl, bridge_col)
        ax2.scatter([vb[0]], [vb[1]], s=175, color=bridge_col, edgecolors="#0f172a", lw=1.5, zorder=6)
        ax2.text(vb[0], vb[1] + 0.52, "bridge v*", ha="center", va="center", fontsize=7.6, fontweight="bold", color="#9a3412")

        ax2.text(1.6, y_off + 1.25, "A1" if not is_leiden else "A1 (conn.)", ha="center", fontsize=7.8, fontweight="bold", color=a1_col)
        ax2.text(5.4, y_off + 1.25, "A2" if not is_leiden else "A2 (conn.)", ha="center", fontsize=7.8, fontweight="bold", color=a2_col)
        ax2.text(8.1, y_off + 1.25, "Comm. B", ha="center", fontsize=7.8, fontweight="bold", color=bridge_col)

        if not is_leiden:
            ax2.text(3.5, y_off + 1.32, "0 internal edges\nbetween A1 & A2!", ha="center", va="center",
                     fontsize=7.3, fontweight="bold", color="#b91c1c")

    draw_panel(5.30, "Louvain: Greedy Move Leaves Community A Disconnected!",
               "Moving bridge v* into B leaves A1 & A2 with the same ID (blue) but 0 internal path!",
               "#991b1b", "#f87171", "#fef2f2", "#2563eb", "#2563eb", "#f97316", False)
    draw_panel(0.35, "Leiden: Refinement Phase Guarantees Connected Communities!",
               "Refinement splits A into connected sub-communities A1 (blue) & A2 (green) first.",
               "#14532d", "#86efac", "#f0fdf4", "#2563eb", "#16a34a", "#f97316", True)

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "fuzzy_knn_and_louvain_vs_leiden.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 4: UMAP Cross-Entropy Forces & Global Distance/Density Distortion
# =============================================================================
def make_umap_distortion():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.24})

    # Left: UMAP kernel q(d) and attractive vs repulsive cross-entropy forces
    d = np.linspace(0.08, 4.5, 350)
    a, b = 1.929, 0.7915
    q = 1.0 / (1.0 + a * (d ** (2 * b)))
    # Gradient magnitudes of CE loss w.r.t. d:
    # Attract (p=1): d/dd [-log q] = 2*a*b * d^(2b-1) / (1 + a*d^(2b))
    # Repel   (p=0): d/dd [-log(1-q)] = 2*b / (d * (1 + a*d^(2b)))
    f_att = (2.0 * a * b * (d ** (2 * b - 1))) / (1.0 + a * (d ** (2 * b)))
    f_rep = (2.0 * b) / (d * (1.0 + a * (d ** (2 * b))))

    ax1.plot(d, q, "-", color="#2563eb", lw=2.5, label=r"Low-dim affinity $q(d) = (1 + a d^{2b})^{-1}$")
    ax1.plot(d, f_att, "--", color="#16a34a", lw=2.2, label="Attractive force (neighbors, p = 1)")
    ax1.plot(d, np.minimum(f_rep, 2.5), ":", color="#dc2626", lw=2.3, label="Repulsive force (non-neighbors, p = 0)")
    ax1.axvspan(2.2, 4.5, color="#94a3b8", alpha=0.14)

    ax1.annotate(
        "Flat long-range plateau:\nrepulsion ≈ 0 beyond d > 2\n⇒ arbitrary inter-cluster gap!",
        xy=(3.1, f_rep[np.argmin(np.abs(d - 3.1))]), xytext=(1.65, 1.38),
        fontsize=7.9, fontweight="bold", color="#991b1b",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#fef2f2", edgecolor="#f87171"),
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
    )

    ax1.set_xlim(0.0, 4.5)
    ax1.set_ylim(0, 2.25)
    ax1.set_xlabel(r"2D Embedding Distance  $d = \Vert y_i - y_j \Vert_2$")
    ax1.set_ylabel("Affinity / Force Magnitude")
    ax1.set_title("UMAP Fuzzy Cross-Entropy Kernel & Forces", pad=10)
    ax1.grid(True, ls=":", alpha=0.4)
    ax1.legend(loc="upper right", fontsize=7.4, framealpha=0.95)

    # Right: True High-Dim Metric vs Distorted 2D UMAP Layout (side-by-side inside ax2)
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.axis("off")
    ax2.set_title("What UMAP Preserves vs. What It Distorts", pad=10)

    rng = np.random.default_rng(606)
    # Left sub-box: True high-dim metric (tight A near diffuse B, far C)
    ax2.add_patch(FancyBboxPatch((0.2, 0.5), 4.45, 8.9, boxstyle="round,pad=0.06,rounding_size=0.2",
                                 facecolor="#f8fafc", edgecolor="#cbd5e1", lw=1.2))
    ax2.text(2.42, 8.95, "True PCA Metric Space", ha="center", fontsize=8.8, fontweight="bold", color="#0f172a")
    A_true = rng.normal([1.3, 6.9], 0.18, (45, 2))
    B_true = rng.normal([2.8, 6.2], 0.68, (45, 2))
    C_true = rng.normal([2.4, 2.25], 0.30, (45, 2))
    ax2.scatter(A_true[:, 0], A_true[:, 1], s=14, color="#2563eb", alpha=0.82)
    ax2.scatter(B_true[:, 0], B_true[:, 1], s=14, color="#16a34a", alpha=0.78)
    ax2.scatter(C_true[:, 0], C_true[:, 1], s=14, color="#7c3aed", alpha=0.82)
    ax2.text(1.3, 7.65, "A (dense)", ha="center", fontsize=7.8, fontweight="bold", color="#1d4ed8")
    ax2.text(3.3, 7.65, "B (diffuse)", ha="center", fontsize=7.8, fontweight="bold", color="#15803d")
    ax2.text(2.4, 0.92, "C (10× farther!)", ha="center", fontsize=7.8, fontweight="bold", color="#6d28d9")
    ax2.annotate("", xy=(2.4, 2.95), xytext=(2.4, 5.2), arrowprops=dict(arrowstyle="<->", color="#dc2626", lw=1.4))
    ax2.text(2.55, 4.05, "d = 20", fontsize=7.8, fontweight="bold", color="#b91c1c")

    # Arrow between boxes
    ax2.annotate("", xy=(5.25, 5.0), xytext=(4.75, 5.0), arrowprops=dict(arrowstyle="->", color="#334155", lw=2.0))
    ax2.text(5.0, 5.55, "UMAP", ha="center", fontsize=7.6, fontweight="bold", color="#334155")

    # Right sub-box: 2D UMAP layout (equalized densities, arbitrary gaps!)
    ax2.add_patch(FancyBboxPatch((5.35, 0.5), 4.45, 8.9, boxstyle="round,pad=0.06,rounding_size=0.2",
                                 facecolor="#fef2f2", edgecolor="#fca5a5", lw=1.2))
    ax2.text(7.57, 8.95, "2D UMAP Projection", ha="center", fontsize=8.8, fontweight="bold", color="#991b1b")
    A_um = rng.normal([6.5, 6.7], 0.44, (45, 2))
    B_um = rng.normal([8.7, 6.5], 0.46, (45, 2))
    C_um = rng.normal([7.5, 3.25], 0.44, (45, 2))
    ax2.scatter(A_um[:, 0], A_um[:, 1], s=14, color="#2563eb", alpha=0.82)
    ax2.scatter(B_um[:, 0], B_um[:, 1], s=14, color="#16a34a", alpha=0.78)
    ax2.scatter(C_um[:, 0], C_um[:, 1], s=14, color="#7c3aed", alpha=0.82)
    ax2.text(6.5, 7.85, "A (expanded)", ha="center", fontsize=7.6, fontweight="bold", color="#1d4ed8")
    ax2.text(8.7, 7.85, "B (shrunk)", ha="center", fontsize=7.6, fontweight="bold", color="#15803d")
    ax2.text(7.5, 1.85, "C (pulled close!)", ha="center", fontsize=7.6, fontweight="bold", color="#6d28d9")
    ax2.text(7.57, 0.92, "Local neighbors ✓ · Global dist & area ✗", ha="center", fontsize=7.5, fontweight="bold", color="#991b1b")

    plt.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "umap_cross_entropy_and_distortion.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 5: Exact Diffusion Pseudotime (DPT) & Waddington KDE Pseudo-Potential
# =============================================================================
def make_dpt_and_waddington():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw={"wspace": 0.25})

    # Left: Exact multi-scale Diffusion Pseudotime on continuous monocyte arc + detached B-cells
    rng = np.random.default_rng(607)
    n_arc = 150
    theta = np.linspace(0.12 * math.pi, 1.15 * math.pi, n_arc)
    r_arc = 2.2 + rng.normal(0.0, 0.16, n_arc)
    arc_pts = np.column_stack([r_arc * np.cos(theta), r_arc * np.sin(theta)])

    # Build symmetric Gaussian affinity matrix on the monocyte arc
    D2 = np.sum((arc_pts[:, None, :] - arc_pts[None, :, :]) ** 2, axis=-1)
    sigma = 0.42
    A = np.exp(-D2 / (2 * sigma ** 2))
    np.fill_diagonal(A, 0.0)
    # Symmetrized k-NN sparsification (k = 12)
    mask = np.zeros_like(A, dtype=bool)
    for i in range(n_arc):
        nbrs = np.argsort(D2[i])[1:13]
        mask[i, nbrs] = True
    A = A * (mask | mask.T)

    deg = np.sum(A, axis=1)
    inv_sqrt_d = 1.0 / np.sqrt(np.maximum(deg, 1e-12))
    T_sym = (A * inv_sqrt_d[:, None]) * inv_sqrt_d[None, :]
    evals, evecs = np.linalg.eigh(T_sym)
    idx_desc = np.argsort(evals)[::-1]
    evals = evals[idx_desc]
    evecs = evecs[:, idx_desc]

    # Exact Neumann multi-scale diffusion matrix M = sum_{k=1}^{15} (lam_k / (1 - lam_k)) psi_k psi_k^T
    k_use = 15
    weights = evals[1:k_use + 1] / np.maximum(1.0 - evals[1:k_use + 1], 1e-9)
    M_coords = evecs[:, 1:k_use + 1] * weights[None, :]
    root = 0
    dpt_raw = np.linalg.norm(M_coords - M_coords[root], axis=1)
    dpt = (dpt_raw - dpt_raw.min()) / (dpt_raw.max() - dpt_raw.min())

    # Detached B-cell island (disconnected component -> dpt = inf)
    b_pts = rng.normal([0.1, -1.55], 0.34, (45, 2))

    sc = ax1.scatter(arc_pts[:, 0], arc_pts[:, 1], c=dpt, cmap="viridis", s=30, edgecolors="white", lw=0.35, zorder=5)
    ax1.scatter(b_pts[:, 0], b_pts[:, 1], color="#cbd5e1", edgecolors="#64748b", s=24, lw=0.6, zorder=4,
                label="B-cells (disconnected: d_dpt = ∞)")
    ax1.scatter([arc_pts[root, 0]], [arc_pts[root, 1]], marker="*", s=220, color="#dc2626", edgecolors="white", lw=1.0,
                zorder=8, label="Root: classical CD14++")

    cb = plt.colorbar(sc, ax=ax1, fraction=0.044, pad=0.02)
    cb.set_label("Diffusion Pseudotime (0 → 1)", fontsize=8.4)
    cb.ax.tick_params(labelsize=7.8)

    ax1.annotate("CD14++ root", xy=(arc_pts[root, 0], arc_pts[root, 1]), xytext=(0.95, 0.05),
                 fontsize=7.9, fontweight="bold", color="#991b1b",
                 arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.1))
    ax1.annotate("FCGR3A++ (CD16+)\nnon-classical end", xy=(arc_pts[-1, 0], arc_pts[-1, 1]), xytext=(-2.65, -2.05),
                 fontsize=7.9, fontweight="bold", color="#1e3a8a",
                 bbox=dict(boxstyle="round,pad=0.25", facecolor="#eff6ff", edgecolor="#93c5fd"),
                 arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.1))
    ax1.text(0.20, -2.65, "Disconnected B-cells\n(no random-walk path!)", ha="center", fontsize=7.6, color="#475569", fontweight="bold")

    ax1.set_xlim(-2.9, 2.9)
    ax1.set_ylim(-3.05, 2.85)
    ax1.set_xlabel("Diffusion / State Coord 1")
    ax1.set_ylabel("Diffusion / State Coord 2")
    ax1.set_title(r"Multi-Scale DPT: $M = (I - \tilde{T})^{-1} - I$", pad=10)
    ax1.grid(True, ls=":", alpha=0.35)
    ax1.legend(loc="upper left", fontsize=7.2, framealpha=0.95)

    # Right: Waddington KDE Pseudo-Potential U(x, y) = -log(p_hat(x, y) + eps)
    rng2 = np.random.default_rng(608)
    c_t = rng2.normal([-2.4, 1.5], 0.55, (180, 2))
    c_b = rng2.normal([2.2, 1.7], 0.45, (110, 2))
    c_m1 = rng2.normal([-1.4, -1.6], 0.50, (130, 2))
    c_m2 = rng2.normal([1.6, -1.6], 0.48, (90, 2))
    # Continuous saddle bridge of intermediate monocytes between CD14 & FCGR3A
    t_br = rng2.uniform(-1.0, 1.2, 42)
    c_br = np.column_stack([t_br, -1.6 + rng2.normal(0.0, 0.28, len(t_br))])
    all_pts = np.vstack([c_t, c_b, c_m1, c_m2, c_br])

    kde = gaussian_kde(all_pts.T, bw_method=0.34)
    xg = np.linspace(-4.2, 4.2, 110)
    yg = np.linspace(-3.5, 3.5, 100)
    XX, YY = np.meshgrid(xg, yg)
    dens = kde(np.vstack([XX.ravel(), YY.ravel()])).reshape(XX.shape)
    U = -np.log(dens + 2e-3)

    cf = ax2.contourf(XX, YY, U, levels=16, cmap="viridis_r", alpha=0.90)
    ax2.contour(XX, YY, U, levels=10, colors="#0f172a", linewidths=0.45, alpha=0.45)
    cb2 = plt.colorbar(cf, ax=ax2, fraction=0.044, pad=0.02)
    cb2.set_label(r"Pseudo-Potential  $U = -\log(\hat{p} + \epsilon)$", fontsize=8.4)
    cb2.ax.tick_params(labelsize=7.8)

    for lbl, xy in [("CD4/CD8 T\n(valley)", (-2.4, 1.5)), ("B-cells\n(valley)", (2.2, 1.7)),
                    ("CD14 Mono", (-1.6, -1.6)), ("FCGR3A\nMono", (1.8, -1.6))]:
        ax2.text(xy[0], xy[1], lbl, ha="center", va="center", fontsize=7.8, fontweight="bold", color="#ffffff",
                 bbox=dict(boxstyle="round,pad=0.18", facecolor="#0f172a", edgecolor="none", alpha=0.58))

    ax2.annotate("Low saddle pass\n(continuous maturation)", xy=(0.1, -1.6), xytext=(0.1, -3.15),
                 ha="center", fontsize=7.7, fontweight="bold", color="#fef08a",
                 bbox=dict(boxstyle="round,pad=0.22", facecolor="#0f172a", edgecolor="#fde047", alpha=0.82),
                 arrowprops=dict(arrowstyle="->", color="#fde047", lw=1.3))
    ax2.text(0.0, 0.25, "High epigenetic ridge\n(low density)", ha="center", va="center",
             fontsize=7.7, fontweight="bold", color="#0f172a",
             bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#94a3b8", alpha=0.88))

    ax2.set_xlim(-4.2, 4.2)
    ax2.set_ylim(-3.5, 3.5)
    ax2.set_xlabel("Embedding Coord 1")
    ax2.set_ylabel("Embedding Coord 2")
    ax2.set_title("Waddington KDE Pseudo-Potential Landscape", pad=10)

    plt.subplots_adjust(left=0.06, right=0.94, top=0.89, bottom=0.14)
    out = os.path.join(OUT_DIR, "diffusion_dpt_and_waddington_kde.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


if __name__ == "__main__":
    make_hvg_and_zscore()
    make_svd_and_curse()
    make_fuzzy_knn_and_leiden()
    make_umap_distortion()
    make_dpt_and_waddington()
