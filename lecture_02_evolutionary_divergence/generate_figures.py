#!/usr/bin/env python3
"""Generate 4 publication-quality figures for Lecture 2, Part 3 (Underlying Algorithms)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from Bio.Align import substitution_matrices

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
# Figure 1: Global (NW) vs. Local (SW) vs. Semi-Global DP Edit Graph
# =============================================================================
def make_dp_modes():
    seq_a = "TTACGTAC"   # rows (len 8)
    seq_b = "ACGTACGG"   # cols (len 8)
    n, m = len(seq_a), len(seq_b)
    match, mismatch, gap = 2, -2, -2

    def run_dp(mode):
        H = np.zeros((n + 1, m + 1), dtype=int)
        ptr = [[None] * (m + 1) for _ in range(n + 1)]
        if mode == "global":
            for i in range(1, n + 1):
                H[i, 0] = i * gap
                ptr[i][0] = (i - 1, 0)
            for j in range(1, m + 1):
                H[0, j] = j * gap
                ptr[0][j] = (0, j - 1)
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
                cands = [
                    (H[i - 1, j - 1] + s, (i - 1, j - 1)),
                    (H[i - 1, j] + gap, (i - 1, j)),
                    (H[i, j - 1] + gap, (i, j - 1)),
                ]
                if mode == "local":
                    cands.append((0, None))
                best_val, best_p = max(cands, key=lambda x: x[0])
                H[i, j] = best_val
                ptr[i][j] = best_p

        if mode == "global":
            cur = (n, m)
        elif mode == "local":
            idx = np.unravel_index(np.argmax(H), H.shape)
            cur = (int(idx[0]), int(idx[1]))
        else:  # semi-global: max along last row or last col
            last_row_j = int(np.argmax(H[n, :]))
            last_col_i = int(np.argmax(H[:, m]))
            if H[n, last_row_j] >= H[last_col_i, m]:
                cur = (n, last_row_j)
            else:
                cur = (last_col_i, m)

        path = [cur]
        while cur is not None:
            if mode == "local" and H[cur] == 0:
                break
            if mode == "semiglobal" and (cur[0] == 0 or cur[1] == 0):
                break
            nxt = ptr[cur[0]][cur[1]]
            if nxt is None:
                break
            path.append(nxt)
            cur = nxt
        return H, path, H[path[0]]

    modes = [
        ("global", "Global (Needleman–Wunsch)\nBorders = i·g, j·g · Corner-to-corner", "TTACGTAC--\n--ACGTACGG  (score +4)"),
        ("local", "Local (Smith–Waterman)\nBorders = 0 · Clamped H ≥ 0", "ACGTAC\nACGTAC  (score +12)"),
        ("semiglobal", "Semi-Global (Free End Gaps)\nBorders = 0 · Unclamped interior", "TTACGTAC--\n--ACGTACGG  (score +12)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 5.1), gridspec_kw={"wspace": 0.24})

    for ax, (mode, title, subtitle) in zip(axes, modes):
        H, path, best_score = run_dp(mode)
        vmin, vmax = (-16, 12) if mode != "local" else (0, 12)
        cmap = "RdYlBu" if mode != "local" else "YlGnBu"
        im = ax.imshow(H, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")

        # Grid lines and cell numbers
        path_set = set(path)
        for i in range(n + 1):
            for j in range(m + 1):
                val = H[i, j]
                in_path = (i, j) in path_set
                norm_v = (val - vmin) / max(1, (vmax - vmin))
                if mode == "local":
                    txt_col = "white" if norm_v > 0.62 else "#0f172a"
                else:
                    txt_col = "white" if (norm_v < 0.16 or norm_v > 0.84) else "#0f172a"
                if in_path:
                    edge_c = "#dc2626" if (i, j) == path[0] else ("#16a34a" if (i, j) == path[-1] else "#d97706")
                    rect = plt.Rectangle((j - 0.48, i - 0.48), 0.96, 0.96, fill=False,
                                         edgecolor=edge_c, linewidth=2.2, zorder=4)
                    ax.add_patch(rect)
                ax.text(
                    j, i, str(val), ha="center", va="center",
                    fontsize=8.6,
                    fontweight="bold" if in_path else "normal",
                    color=txt_col,
                    zorder=6,
                )

        # Draw translucent connecting traceback path behind text
        xs = [p[1] for p in path]
        ys = [p[0] for p in path]
        ax.plot(xs, ys, color="#f59e0b", linewidth=3.4, alpha=0.35, zorder=3)

        ax.set_xticks(range(m + 1))
        ax.set_yticks(range(n + 1))
        ax.set_xticklabels(["–"] + list(seq_b), fontsize=9.5, fontweight="bold")
        ax.set_yticklabels(["–"] + list(seq_a), fontsize=9.5, fontweight="bold")
        ax.xaxis.tick_top()
        ax.set_title(title, pad=14, fontsize=10.8)
        ax.set_xlabel(subtitle, fontsize=9.2, fontfamily="DejaVu Sans Mono", labelpad=8,
                      bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8fafc", edgecolor="#cbd5e1"))

        for sp in ax.spines.values():
            sp.set_color("#64748b")

    fig.text(
        0.5, 0.015,
        "Scoring: match = +2, mismatch = −2, linear gap = −2   ·   Green box = alignment start   ·   Amber boxes = traceback path   ·   Red box = optimal score",
        ha="center", fontsize=9.3, color="#475569", style="italic"
    )
    plt.subplots_adjust(left=0.04, right=0.98, top=0.81, bottom=0.14)
    out = os.path.join(OUT_DIR, "dp_edit_graph_modes.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figures 2a & 2b: BLOSUM62 Heatmap + PAM vs. BLOSUM Evolutionary Scale
# =============================================================================
def make_blosum_pam():
    B = substitution_matrices.load("BLOSUM62")
    groups = [
        ("Aliphatic", list("AVILM")),
        ("Aromatic",  list("FWY")),
        ("Polar",     list("STNQ")),
        ("Basic (+)", list("KRH")),
        ("Acidic (−)",list("DE")),
        ("Special",   list("CGP")),
    ]
    aa = [r for _, rs in groups for r in rs]
    mat = np.array([[int(B[a][b]) for b in aa] for a in aa])

    # ---- Figure 2a: Standalone BLOSUM62 Heatmap ----
    fig1, ax1 = plt.subplots(figsize=(7.4, 6.2))
    im = ax1.imshow(mat, cmap="RdBu_r", vmin=-5, vmax=11, aspect="equal")
    for i in range(20):
        for j in range(20):
            v = mat[i, j]
            txt_c = "white" if (v >= 6 or v <= -4) else "#0f172a"
            ax1.text(j, i, str(v), ha="center", va="center",
                     fontsize=7.8, fontweight="bold" if i == j else "normal", color=txt_c)

    pos = 0
    for idx, (gname, rs) in enumerate(groups):
        mid = pos + (len(rs) - 1) / 2.0
        ax1.text(-1.55, mid, gname, ha="right", va="center", fontsize=8.0,
                 fontweight="bold", color="#475569")
        pos += len(rs)
        if idx < len(groups) - 1:
            ax1.axhline(pos - 0.5, color="#1e293b", lw=1.2, ls="--", alpha=0.65)
            ax1.axvline(pos - 0.5, color="#1e293b", lw=1.2, ls="--", alpha=0.65)

    ax1.set_xticks(range(20))
    ax1.set_yticks(range(20))
    ax1.set_xticklabels(aa, fontsize=9.0, fontweight="bold")
    ax1.set_yticklabels(aa, fontsize=9.0, fontweight="bold")
    ax1.set_title("BLOSUM62 Log-Odds Matrix (Half-Bits, Grouped by Chemistry)", pad=10, fontsize=11.5)
    cbar = fig1.colorbar(im, ax=ax1, fraction=0.044, pad=0.025)
    cbar.set_label(r"Score $S(a, b) = 2\log_2(p_{ab} / q_a q_b)$", fontsize=9.5)

    w_idx = aa.index("W")
    c_idx = aa.index("C")
    for idx_hi, col_hi in [(w_idx, "#fbbf24"), (c_idx, "#22c55e")]:
        rect = plt.Rectangle((idx_hi - 0.5, idx_hi - 0.5), 1, 1, fill=False,
                             edgecolor=col_hi, linewidth=2.4, zorder=5)
        ax1.add_patch(rect)

    plt.subplots_adjust(left=0.17, right=0.92, top=0.92, bottom=0.07)
    out1 = os.path.join(OUT_DIR, "blosum62_heatmap.png")
    fig1.savefig(out1, dpi=200)
    plt.close(fig1)
    print("Saved", out1)

    # ---- Figure 2b: Standalone PAM vs. BLOSUM Evolutionary Scale ----
    fig2, ax2 = plt.subplots(figsize=(7.2, 5.4))
    pam_nums = np.array([30, 60, 100, 120, 160, 200, 250])
    pam_H    = np.array([2.57, 1.79, 1.18, 0.98, 0.70, 0.51, 0.36])

    blosum_equiv_pam = np.array([100, 120, 160, 200, 250])
    blosum_H = np.array([1.18, 0.99, 0.70, 0.52, 0.38])

    ax2.plot(pam_nums, pam_H, "-o", color="#2563eb", lw=2.6, ms=7.0, label="PAM (Markov extrapolation Mⁿ)")
    ax2.scatter(blosum_equiv_pam, blosum_H, color="#dc2626", s=95, marker="D", zorder=5,
                label="BLOSUM (Empirical BLOCKS clustering)")

    pam_ann = [
        (30, 2.57, "PAM30\n(close)", (14, -4)),
        (120, 0.98, "PAM120", (-46, -18)),
        (250, 0.36, "PAM250\n(distant)", (-50, -24)),
    ]
    for x, y, lbl, off in pam_ann:
        ax2.annotate(lbl, (x, y), xytext=off, textcoords="offset points",
                     fontsize=9.2, color="#1d4ed8", fontweight="bold")

    blos_ann = [
        (100, 1.18, "BLOSUM90", (12, 16)),
        (120, 0.99, "BLOSUM80", (14, 12)),
        (160, 0.70, "BLOSUM62\n(default)", (12, 16)),
        (250, 0.38, "BLOSUM45", (-42, 20)),
    ]
    for x, y, lbl, off in blos_ann:
        ax2.annotate(lbl, (x, y), xytext=off, textcoords="offset points",
                     fontsize=9.2, color="#b91c1c", fontweight="bold",
                     arrowprops=dict(arrowstyle="-", color="#dc2626", lw=0.9, alpha=0.65))

    ax2.annotate(
        "", xy=(255, 2.32), xytext=(65, 2.32),
        arrowprops=dict(arrowstyle="->", lw=2.2, color="#334155")
    )
    ax2.text(160, 2.41, "More Divergent Homologs  →", ha="center", fontsize=10.0, fontweight="bold", color="#0f172a")
    ax2.text(160, 2.13, "PAM number ↑   |   BLOSUM number ↓", ha="center", fontsize=9.5,
             color="#475569", bbox=dict(boxstyle="round,pad=0.30", facecolor="#f1f5f9", edgecolor="#cbd5e1"))

    ax2.set_xlim(15, 285)
    ax2.set_ylim(0.0, 2.80)
    ax2.set_xlabel("Evolutionary Divergence Scale (PAM Units)", fontsize=10.5)
    ax2.set_ylabel("Relative Entropy H (bits / aligned site)", fontsize=10.5)
    ax2.set_title("PAM vs. BLOSUM: Opposing Numbering & Information", pad=10, fontsize=11.5)
    ax2.grid(True, ls=":", alpha=0.45)
    ax2.legend(loc="center right", bbox_to_anchor=(0.99, 0.64), fontsize=8.8, framealpha=0.95)

    plt.subplots_adjust(left=0.11, right=0.97, top=0.90, bottom=0.12)
    out2 = os.path.join(OUT_DIR, "pam_vs_blosum_scale.png")
    fig2.savefig(out2, dpi=200)
    plt.close(fig2)
    print("Saved", out2)


# =============================================================================
# Figure 3: Multiple-Hits Saturation & Delta-Method Variance Explosion
# =============================================================================
def make_saturation_variance():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.0), gridspec_kw={"wspace": 0.24})

    # Left panel: p(d) saturation curves
    d = np.linspace(0, 3.0, 400)
    p_lin = np.minimum(d, 1.0)
    p_pois = 1.0 - np.exp(-d)
    p_jc = 0.75 * (1.0 - np.exp(-4.0 / 3.0 * d))
    alpha = 1.0
    p_gamma = 1.0 - (1.0 + d / alpha) ** (-alpha)

    ax1.plot(d[d <= 1.05], d[d <= 1.05], "--", color="#64748b", lw=1.6, label="No multiple hits (p = d)")
    ax1.plot(d, p_pois, color="#2563eb", lw=2.6, label=r"Poisson (protein): $p = 1 - e^{-d}$")
    ax1.plot(d, p_jc, color="#dc2626", lw=2.6, label=r"JC69 (DNA, 4 bases): $p = \frac{3}{4}(1 - e^{-4d/3})$")
    ax1.plot(d, p_gamma, color="#7c3aed", lw=2.3, ls="-.", label=r"Gamma–Poisson ($\alpha = 1$): $p = d / (1 + d)$")

    # Shade hidden hits between p=d and Poisson (up to d=1) and show JC69 ceiling
    ax1.axhline(0.75, color="#dc2626", ls=":", lw=1.4, alpha=0.8)
    ax1.text(2.96, 0.768, "DNA saturation ceiling (p = 0.75)", ha="right", va="bottom",
             fontsize=8.6, color="#b91c1c", fontweight="bold")
    ax1.axhline(0.95, color="#2563eb", ls=":", lw=1.1, alpha=0.5)
    ax1.text(2.96, 0.962, "20-aa random limit (~0.94–0.95)", ha="right", va="bottom",
             fontsize=8.2, color="#1d4ed8")

    # Vertical double arrow showing hidden hits at d = 1.2 (Poisson: p = 0.699)
    d0 = 1.204
    p0 = 1.0 - np.exp(-d0)  # 0.70
    ax1.scatter([d0], [p0], color="#2563eb", s=45, zorder=5)
    ax1.annotate(
        "Observed p = 0.70\nTrue d = 1.20 (+72% hidden hits!)",
        xy=(d0, p0), xytext=(1.42, 0.34),
        fontsize=8.7, fontweight="bold", color="#1e3a8a",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#eff6ff", edgecolor="#93c5fd"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.4, connectionstyle="arc3,rad=-0.15"),
    )

    ax1.set_xlim(0, 3.0)
    ax1.set_ylim(0, 1.04)
    ax1.set_xlabel("True Evolutionary Distance d (substitutions / site)")
    ax1.set_ylabel("Observed Fraction of Differing Sites p")
    ax1.set_title("Multiple-Hits Saturation: Observed p vs. True d", pad=10)
    ax1.grid(True, ls=":", alpha=0.45)
    ax1.legend(loc="lower right", fontsize=8.4, framealpha=0.95)

    # Right panel: Delta-method standard error SE(d_hat) vs p
    p_p = np.linspace(0.01, 0.93, 400)
    p_j = np.linspace(0.01, 0.725, 400)

    se_pois_150 = np.sqrt(p_p / ((1.0 - p_p) * 150.0))
    se_pois_500 = np.sqrt(p_p / ((1.0 - p_p) * 500.0))
    se_jc_150   = np.sqrt(p_j * (1.0 - p_j) / (((1.0 - 4.0 / 3.0 * p_j) ** 2) * 150.0))
    se_jc_500   = np.sqrt(p_j * (1.0 - p_j) / (((1.0 - 4.0 / 3.0 * p_j) ** 2) * 500.0))

    ax2.plot(p_j, se_jc_150, color="#dc2626", lw=2.5, label="JC69 DNA (L = 150)")
    ax2.plot(p_j, se_jc_500, color="#f97316", lw=2.2, ls="--", label="JC69 DNA (L = 500)")
    ax2.plot(p_p, se_pois_150, color="#2563eb", lw=2.5, label="Poisson protein (L = 150 aa, ~HBA1)")
    ax2.plot(p_p, se_pois_500, color="#0ea5e9", lw=2.2, ls="--", label="Poisson protein (L = 500 aa)")

    ax2.axvline(0.75, color="#dc2626", ls=":", lw=1.5, alpha=0.85)
    ax2.text(0.77, 0.54, "JC69 pole\np = 0.75\n(DNA limit)", ha="left", va="top",
             fontsize=8.8, color="#b91c1c", fontweight="bold")
    ax2.axvspan(0.75, 0.96, color="#fee2e2", alpha=0.35)

    ax2.annotate(
        "Delta-method SE diverges:\nshort genes lose clock precision\nlong before p = 1",
        xy=(0.84, np.sqrt(0.84 / (0.16 * 150))), xytext=(0.14, 0.37),
        fontsize=8.6, color="#1e293b",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8fafc", edgecolor="#cbd5e1"),
        arrowprops=dict(arrowstyle="->", color="#2563eb", lw=1.3),
    )

    ax2.set_xlim(0, 0.96)
    ax2.set_ylim(0, 0.58)
    ax2.set_xlabel("Observed Fraction of Differing Sites p")
    ax2.set_ylabel(r"Standard Error $\mathrm{SE}(\hat{d})$ (subs / site)")
    ax2.set_title(r"Variance Explosion Near Saturation ($\mathrm{Var}(\hat{d}) \approx \frac{p}{(1-p)L}$)", pad=10)
    ax2.grid(True, ls=":", alpha=0.45)
    ax2.legend(loc="upper left", fontsize=8.3, framealpha=0.95)

    plt.subplots_adjust(left=0.06, right=0.98, top=0.89, bottom=0.13)
    out = os.path.join(OUT_DIR, "saturation_multiple_hits.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


# =============================================================================
# Figure 4: UPGMA Failure vs. Neighbor-Joining Success on 4-Taxon Additive Tree
# =============================================================================
def make_upgma_vs_nj():
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 5.0), gridspec_kw={"wspace": 0.14})
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis("off")

    # ---- Panel 1: True Additive Tree + Distance Matrix ----
    ax = axes[0]
    ax.set_title("1. True Additive Tree (Unequal Rates)", fontsize=10.8, pad=6)

    # Unrooted tree ((A:1, B:9):3, (C:1, D:9))
    u = (3.4, 6.9)
    v = (6.6, 6.9)
    pA = (2.1, 8.6)
    pB = (0.7, 5.0)
    pC = (7.9, 8.6)
    pD = (9.3, 5.0)

    for (x1, y1), (x2, y2), lw in [(u, v, 2.8), (u, pA, 2.4), (u, pB, 2.4), (v, pC, 2.4), (v, pD, 2.4)]:
        ax.plot([x1, x2], [y1, y2], color="#334155", lw=lw, solid_capstyle="round")

    # Node markers
    ax.scatter([u[0], v[0]], [u[1], v[1]], color="#64748b", s=45, zorder=5)
    for pt, name, col, dx, dy in [
        (pA, "A", "#2563eb", -0.35, 0.35),
        (pB, "B", "#2563eb", -0.35, -0.55),
        (pC, "C", "#16a34a", 0.35, 0.35),
        (pD, "D", "#16a34a", 0.35, -0.55),
    ]:
        ax.scatter([pt[0]], [pt[1]], color=col, s=85, zorder=6, edgecolors="white", lw=1.2)
        ax.text(pt[0] + dx, pt[1] + dy, name, ha="center", va="center",
                fontsize=11.5, fontweight="bold", color=col)

    # Branch length labels
    ax.text(5.0, 7.35, "3", ha="center", fontsize=10, fontweight="bold", color="#0f172a")
    ax.text(2.95, 8.05, "1", ha="center", fontsize=10, fontweight="bold", color="#2563eb")
    ax.text(1.75, 6.25, "9", ha="center", fontsize=10, fontweight="bold", color="#dc2626")
    ax.text(7.05, 8.05, "1", ha="center", fontsize=10, fontweight="bold", color="#16a34a")
    ax.text(8.25, 6.25, "9", ha="center", fontsize=10, fontweight="bold", color="#dc2626")

    # Compact distance matrix box at bottom
    mat_txt = (
        "Pairwise Distances D:\n"
        "     A    B    C    D\n"
        "A [  0   10    5*  13 ]\n"
        "B [ 10    0   13   21 ]\n"
        "C [  5*  13    0   10 ]\n"
        "D [ 13   21   10    0 ]\n"
        "* min D is d(A,C)=5 < d(A,B)=10!"
    )
    ax.text(
        5.0, 2.05, mat_txt, ha="center", va="center",
        fontsize=8.3, fontfamily="DejaVu Sans Mono", color="#0f172a",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#f8fafc", edgecolor="#cbd5e1")
    )

    # ---- Panel 2: UPGMA (Wrong Topology!) ----
    ax = axes[1]
    ax.set_title("2. UPGMA (Greedy min D → FAILS)", fontsize=10.8, pad=6, color="#b91c1c")

    # Draw UPGMA ultrametric dendrogram:
    # Leaves at x = 8.4: A(y=8.2), C(y=6.6), D(y=5.0), B(y=3.4)
    # Scale height h in [0, 8.125] -> x = 8.4 - h * 0.85
    def hx(h):
        return 8.4 - h * 0.85

    yA, yC, yD, yB = 8.2, 6.6, 5.0, 3.4
    yAC = (yA + yC) / 2
    yACD = (yAC + yD) / 2
    yRoot = (yACD + yB) / 2

    x0 = hx(0)
    xAC = hx(2.5)
    xACD = hx(5.75)
    xRoot = hx(8.125)

    # Lines
    segs = [
        ((xAC, yA), (x0, yA)), ((xAC, yC), (x0, yC)), ((xAC, yA), (xAC, yC)),
        ((xACD, yAC), (xAC, yAC)), ((xACD, yD), (x0, yD)), ((xACD, yAC), (xACD, yD)),
        ((xRoot, yACD), (xACD, yACD)), ((xRoot, yB), (x0, yB)), ((xRoot, yACD), (xRoot, yB)),
    ]
    for (x1, y1), (x2, y2) in segs:
        ax.plot([x1, x2], [y1, y2], color="#dc2626", lw=2.3)

    for y, name, col in [(yA, "A", "#2563eb"), (yC, "C", "#16a34a"), (yD, "D", "#16a34a"), (yB, "B", "#2563eb")]:
        ax.scatter([x0], [y], color=col, s=75, zorder=5, edgecolors="white", lw=1.1)
        ax.text(x0 + 0.45, y, name, va="center", ha="left", fontsize=11.5, fontweight="bold", color=col)

    ax.text((xAC + x0) / 2, yA + 0.32, "2.5", ha="center", fontsize=8.5, color="#7f1d1d")
    ax.text((xACD + xAC) / 2, yAC + 0.32, "3.25", ha="center", fontsize=8.5, color="#7f1d1d")
    ax.text((xRoot + xACD) / 2, yACD + 0.32, "2.375", ha="center", fontsize=8.5, color="#7f1d1d")
    ax.text((xACD + x0) / 2, yD + 0.32, "5.75", ha="center", fontsize=8.5, color="#7f1d1d")
    ax.text((xRoot + x0) / 2, yB - 0.48, "8.125", ha="center", fontsize=8.5, color="#7f1d1d")

    # Wrong badge
    badge_err = (
        "WRONG TOPOLOGY: ((A, C), D), B\n"
        "UPGMA merges min D = d(A, C) = 5,\n"
        "grouping the two slow-evolving\n"
        "non-sisters instead of (A, B)!"
    )
    ax.text(
        5.0, 1.45, badge_err, ha="center", va="center",
        fontsize=8.6, color="#991b1b", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#fef2f2", edgecolor="#f87171", lw=1.3)
    )

    # ---- Panel 3: Neighbor-Joining (Exact!) ----
    ax = axes[2]
    ax.set_title("3. Neighbor-Joining (min Q → EXACT)", fontsize=10.8, pad=6, color="#15803d")

    # Draw recovered unrooted tree (same geometry, highlighted green/blue)
    for (x1, y1), (x2, y2), lw in [(u, v, 2.8), (u, pA, 2.5), (u, pB, 2.5), (v, pC, 2.5), (v, pD, 2.5)]:
        ax.plot([x1, x2], [y1, y2], color="#15803d", lw=lw, solid_capstyle="round")

    ax.scatter([u[0], v[0]], [u[1], v[1]], color="#15803d", s=50, zorder=5)
    for pt, name, col, dx, dy in [
        (pA, "A", "#2563eb", -0.35, 0.35),
        (pB, "B", "#2563eb", -0.35, -0.55),
        (pC, "C", "#16a34a", 0.35, 0.35),
        (pD, "D", "#16a34a", 0.35, -0.55),
    ]:
        ax.scatter([pt[0]], [pt[1]], color=col, s=85, zorder=6, edgecolors="white", lw=1.2)
        ax.text(pt[0] + dx, pt[1] + dy, name, ha="center", va="center",
                fontsize=11.5, fontweight="bold", color=col)

    ax.text(5.0, 7.35, "3.0", ha="center", fontsize=10, fontweight="bold", color="#15803d")
    ax.text(2.95, 8.05, "1.0", ha="center", fontsize=10, fontweight="bold", color="#15803d")
    ax.text(1.75, 6.25, "9.0", ha="center", fontsize=10, fontweight="bold", color="#15803d")
    ax.text(7.05, 8.05, "1.0", ha="center", fontsize=10, fontweight="bold", color="#15803d")
    ax.text(8.25, 6.25, "9.0", ha="center", fontsize=10, fontweight="bold", color="#15803d")

    q_txt = (
        "EXACT TOPOLOGY & BRANCHES\n"
        "R = [28, 44, 28, 44]\n"
        "Q(A, B) = 2(10) − 28 − 44 = −52  ✓\n"
        "Q(A, C) = 2(5)  − 28 − 28 = −46\n"
        "Subtracting R_i, R_j cancels fast branches!"
    )
    ax.text(
        5.0, 1.80, q_txt, ha="center", va="center",
        fontsize=8.3, fontfamily="DejaVu Sans Mono", color="#14532d",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#f0fdf4", edgecolor="#4ade80", lw=1.3)
    )

    plt.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.04)
    out = os.path.join(OUT_DIR, "upgma_vs_nj_counterexample.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


if __name__ == "__main__":
    make_dp_modes()
    make_blosum_pam()
    make_saturation_variance()
    make_upgma_vs_nj()
