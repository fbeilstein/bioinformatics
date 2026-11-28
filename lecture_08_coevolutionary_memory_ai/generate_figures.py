#!/usr/bin/env python3
"""
Generate 5 reproducible, publication-grade 200-DPI figures for
Lecture 8, Part 3: Underlying Algorithms (Co-Evolutionary Memory & AI).

Leaves all 3 existing Part 4 figures untouched:
  - attention_contact_map.png
  - dms_mutational_matrix.png
  - zero_shot_dms_correlation.png

Creates 5 new figures in lecture_08_coevolutionary_memory_ai/images/:
  1. epistasis_and_dca_couplings.png
  2. attention_scaling_and_rope.png
  3. mlm_apc_contact_emergence.png
  4. zeroshot_marginal_and_dms_metrics.png
  5. esmfold_vs_af2_and_scaling.png
"""

from pathlib import Path
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

OUT_DIR = Path(__file__).resolve().parent / "images"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11.2,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.8,
    "axes.edgecolor": "#475569",
    "axes.linewidth": 0.9,
    "xtick.labelsize": 8.8,
    "ytick.labelsize": 8.8,
    "legend.fontsize": 8.4,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

C_BLUE = "#2563EB"
C_PURPLE = "#7C3AED"
C_TEAL = "#0D9488"
C_GREEN = "#16A34A"
C_ORANGE = "#EA580C"
C_RED = "#DC2626"
C_SLATE = "#334155"
C_MUTED = "#64748B"


# ---------------------------------------------------------------------------
# Exact Statistical Physics, RoPE, Attention, APC & Correlation Helpers
# ---------------------------------------------------------------------------

def solve_potts_4site(J_val=2.2, beta=1.0, h=(0.45, 0.45, 0.45, 0.40)):
    """
    Exact 16-state Potts / Ising model on 4 binary sites (+1/-1) with
    direct 3D couplings (0,1) and (1,2), making (0,2) an indirect transitive chain,
    and site 3 uncoupled.
    """
    J = np.zeros((4, 4))
    J[0, 1] = J[1, 0] = J_val
    J[1, 2] = J[2, 1] = J_val

    states = np.array([[1 if (m >> i) & 1 else -1 for i in range(4)] for m in range(16)], dtype=float)
    h_arr = np.array(h, dtype=float)
    # E(x) = -sum h_i x_i - sum_{i<j} J_ij x_i x_j
    energies = -(states @ h_arr) - 0.5 * np.sum((states @ J) * states, axis=1)
    log_w = -beta * energies
    log_w -= np.max(log_w)
    probs = np.exp(log_w)
    probs /= np.sum(probs)

    # Single-site marginals P(x_i = +1)
    p_plus = np.array([np.sum(probs[states[:, i] == 1]) for i in range(4)])

    # Exact pairwise Mutual Information (nats)
    MI = np.zeros((4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            mi = 0.0
            for si in (-1, 1):
                pi = np.sum(probs[states[:, i] == si])
                for sj in (-1, 1):
                    pj = np.sum(probs[states[:, j] == sj])
                    pij = np.sum(probs[(states[:, i] == si) & (states[:, j] == sj)])
                    if pij > 1e-15:
                        mi += pij * np.log(pij / (pi * pj))
            MI[i, j] = MI[j, i] = mi

    # Log-odds relative to WT = (+1, +1, +1, +1) (index 15)
    def log_p(vec):
        idx = sum((1 << i) if vec[i] == 1 else 0 for i in range(4))
        return float(np.log(probs[idx]))

    def log_pssm(vec):
        return float(sum(np.log(p_plus[i] if vec[i] == 1 else (1.0 - p_plus[i])) for i in range(4)))

    wt = [1, 1, 1, 1]
    mut_a = [-1, 1, 1, 1]
    mut_b = [1, -1, 1, 1]
    mut_ab = [-1, -1, 1, 1]

    ds_a = log_p(mut_a) - log_p(wt)
    ds_b = log_p(mut_b) - log_p(wt)
    ds_ab = log_p(mut_ab) - log_p(wt)
    ds_ab_add = ds_a + ds_b
    epi = ds_ab - ds_ab_add

    pssm_a = log_pssm(mut_a) - log_pssm(wt)
    pssm_b = log_pssm(mut_b) - log_pssm(wt)
    pssm_ab = log_pssm(mut_ab) - log_pssm(wt)

    return {
        "J": J, "MI": MI, "probs": probs, "p_plus": p_plus,
        "ds_a": ds_a, "ds_b": ds_b, "ds_ab": ds_ab, "ds_ab_add": ds_ab_add, "epi": epi,
        "pssm_a": pssm_a, "pssm_b": pssm_b, "pssm_ab": pssm_ab,
    }


def rope_rotate(vec, pos, base=10000.0):
    """Apply exact 2D block Rotary Positional Embedding R_pos to vector vec (even dim)."""
    d = len(vec)
    m = np.arange(0, d, 2, dtype=float)
    inv_freq = base ** (-m / d)
    angles = pos * inv_freq
    cos_a, sin_a = np.cos(angles), np.sin(angles)
    out = np.empty_like(vec, dtype=float)
    out[0::2] = vec[0::2] * cos_a - vec[1::2] * sin_a
    out[1::2] = vec[0::2] * sin_a + vec[1::2] * cos_a
    return out


def stable_softmax(S):
    m = np.max(S, axis=-1, keepdims=True)
    e = np.exp(S - m)
    return e / np.sum(e, axis=-1, keepdims=True)


def apc_correct(A):
    """Symmetrize and apply Average Product Correction (APC) to N x N attention matrix."""
    A_sym = 0.5 * (A + A.T)
    row_m = np.mean(A_sym, axis=1, keepdims=True)
    col_m = np.mean(A_sym, axis=0, keepdims=True)
    tot_m = float(np.mean(A_sym))
    if abs(tot_m) < 1e-12:
        return A_sym.copy()
    return A_sym - (row_m @ col_m) / tot_m


# ---------------------------------------------------------------------------
# Figure 1: Epistasis (PSSM vs Potts) & Direct vs Indirect Coupling (MI vs J)
# ---------------------------------------------------------------------------

def fig1_epistasis_and_dca():
    res = solve_potts_4site(J_val=1.8, beta=1.0, h=(0.45, 0.45, 0.45, 0.40))

    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.05, 1.05], wspace=0.26)

    # Left: Double-mutant epistasis comparison
    ax0 = fig.add_subplot(gs[0])
    labels = [
        "Single Mut A\n(breaks 1–2)",
        "Single Mut B\n(breaks 1–2)",
        "Double AB\n(PSSM Additive)",
        "Double AB\n(Potts Rescue)",
    ]
    vals = [res["ds_a"], res["ds_b"], res["ds_ab_add"], res["ds_ab"]]
    cols = [C_ORANGE, C_ORANGE, C_RED, C_GREEN]
    x = np.arange(4)
    bars = ax0.bar(x, vals, width=0.56, color=cols, edgecolor=C_SLATE, lw=0.9, alpha=0.92)
    ax0.axhline(0, color=C_SLATE, lw=1.0)

    for b, v in zip(bars, vals):
        va = "top" if v < 0 else "bottom"
        dy = -0.35 if v < 0 else 0.25
        ax0.text(b.get_x() + b.get_width() / 2, v + dy, f"{v:+.2f}",
                 ha="center", va=va, fontsize=8.8, fontweight="bold", color=C_SLATE)

    # Epistasis bracket cleanly positioned between bars 2 and 3 without touching -5.40
    ax0.annotate("", xy=(2.92, res["ds_ab"] - 1.25), xytext=(2.08, res["ds_ab_add"] + 0.35),
                 arrowprops=dict(arrowstyle="<->", color=C_PURPLE, lw=1.9, ls="--"))
    ax0.text(2.50, 0.5 * (res["ds_ab"] + res["ds_ab_add"]) - 0.15,
             f"Compensatory Epistasis\nΔΔs_epi = {res['epi']:+.2f}",
             ha="center", va="center", fontsize=8.1, fontweight="bold", color=C_PURPLE,
             bbox=dict(boxstyle="round,pad=0.24", fc="#F5F3FF", ec="#DDD6FE", alpha=0.96))

    ax0.set_xticks(x)
    ax0.set_xticklabels(labels, fontsize=8.2)
    ax0.set_ylabel(r"Log-Odds Fitness Change  $\Delta s = \log\frac{P(\mathbf{x}^{\mathrm{mut}})}{P(\mathbf{x}^{\mathrm{wt}})}$")
    ax0.set_ylim(min(vals) - 2.0, 2.2)
    ax0.set_title("PSSM Additivity Failure vs. Potts Compensatory Rescue", pad=8)
    ax0.grid(True, axis="y", alpha=0.25)

    # Right: Direct Coupling |J_ij| vs Pairwise Mutual Information MI_ij
    ax1 = fig.add_subplot(gs[1])
    pairs = [(0, 1, "Pair (1, 2)\nDirect 3D"), (1, 2, "Pair (2, 3)\nDirect 3D"),
             (0, 2, "Pair (1, 3)\nIndirect Chain!"), (0, 3, "Pair (1, 4)\nUncoupled")]
    mi_vals = [max(0.0, res["MI"][i, j]) for i, j, _ in pairs]
    j_vals = [abs(res["J"][i, j]) for i, j, _ in pairs]
    mi_norm = [0.0 if v < 1e-12 else v / max(mi_vals) for v in mi_vals]
    j_norm = [0.0 if v < 1e-12 else v / max(j_vals) for v in j_vals]

    xp = np.arange(4)
    w = 0.35
    b1 = ax1.bar(xp - w / 2, mi_norm, width=w, color=C_ORANGE, edgecolor=C_SLATE, lw=0.8, label="Pairwise Mutual Info (MI / max)")
    b2 = ax1.bar(xp + w / 2, j_norm, width=w, color=C_BLUE, edgecolor=C_SLATE, lw=0.8, label="Potts Direct Coupling (|J_ij| / max)")

    for b, v in zip(b1, mi_norm):
        ax1.text(b.get_x() + b.get_width() / 2, v + 0.025, f"{v:.2f}", ha="center", va="bottom", fontsize=8.0, color="#9A3412", fontweight="bold")
    for b, v in zip(b2, j_norm):
        ax1.text(b.get_x() + b.get_width() / 2, v + 0.025, f"{v:.2f}", ha="center", va="bottom", fontsize=8.0, color="#1E3A8A", fontweight="bold")

    ax1.annotate("Transitive false positive!\n(1↔2↔3 chain fools MI,\nbut Potts J₁₃ = 0.00)",
                 xy=(2 - w / 2, mi_norm[2] + 0.06), xytext=(2.25, 0.90),
                 fontsize=8.0, fontweight="bold", color=C_RED, ha="center",
                 arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.3),
                 bbox=dict(boxstyle="round,pad=0.25", fc="#FEF2F2", ec="#FECACA", alpha=0.95))

    ax1.set_xticks(xp)
    ax1.set_xticklabels([p[2] for p in pairs], fontsize=8.2)
    ax1.set_ylabel("Normalized Coupling Score (rel. to max)")
    ax1.set_ylim(0, 1.24)
    ax1.set_title("Disentangling Direct 3D Contacts from Indirect Chains", pad=8)
    ax1.legend(loc="upper left", framealpha=0.94, fontsize=8.0)
    ax1.grid(True, axis="y", alpha=0.25)

    fig.subplots_adjust(left=0.07, right=0.985, top=0.89, bottom=0.17)
    fig.savefig(OUT_DIR / "epistasis_and_dca_couplings.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: Attention 1/sqrt(d_k) Variance Scaling & RoPE Shift Invariance
# ---------------------------------------------------------------------------

def fig2_attention_scaling_and_rope():
    rng = np.random.default_rng(123)
    dks = [4, 8, 16, 32, 64, 128, 256]
    L = 32
    ent_scaled, ent_unscaled = [], []
    for dk in dks:
        Q = rng.standard_normal((250, dk))
        K = rng.standard_normal((250, L, dk))
        dots = np.einsum("nd,nld->nl", Q, K)
        A_un = stable_softmax(dots)
        A_sc = stable_softmax(dots / np.sqrt(dk))
        h_un = -np.sum(A_un * np.log(np.clip(A_un, 1e-15, 1.0)), axis=-1) / np.log(L)
        h_sc = -np.sum(A_sc * np.log(np.clip(A_sc, 1e-15, 1.0)), axis=-1) / np.log(L)
        ent_unscaled.append(float(np.mean(h_un)))
        ent_scaled.append(float(np.mean(h_sc)))

    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 1.05], wspace=0.25)

    # Left: Softmax Entropy vs Head Dimension d_k
    ax0 = fig.add_subplot(gs[0])
    ax0.axhspan(0.0, 0.25, color="#FECACA", alpha=0.35, zorder=0)
    ax0.text(4.3, 0.05, "Gradient saturation zone (one-hot collapse!)",
             ha="left", va="bottom", fontsize=8.0, fontweight="bold", color=C_RED)

    ax0.plot(dks, ent_scaled, "-o", color=C_BLUE, lw=2.4, ms=5.5, label=r"Scaled: $\mathrm{softmax}(QK^T / \sqrt{d_k})$")
    ax0.plot(dks, ent_unscaled, "--s", color=C_RED, lw=2.4, ms=5.5, label=r"Unscaled: $\mathrm{softmax}(QK^T)$")

    ax0.set_xscale("log", base=2)
    ax0.set_xticks(dks)
    ax0.set_xticklabels([str(d) for d in dks])
    ax0.set_ylim(-0.02, 1.08)
    ax0.set_xlabel(r"Attention Head Dimension $d_k$ (log₂ scale, $L = 32$ tokens)")
    ax0.set_ylabel("Normalized Softmax Entropy  H(A) / log L")
    ax0.set_title(r"Why Divide by $\sqrt{d_k}$? Preserving Softmax Gradients", pad=8)
    ax0.legend(loc="center right", bbox_to_anchor=(0.98, 0.68), framealpha=0.95, fontsize=8.1)
    ax0.text(0.97, 0.42, r"$\mathrm{Var}(\mathbf{q}\cdot\mathbf{k}) = d_k$" + "\n" + r"$\mathrm{Var}(\mathbf{q}\cdot\mathbf{k}/\sqrt{d_k}) = 1$",
             transform=ax0.transAxes, ha="right", va="center", fontsize=8.4, color=C_SLATE,
             bbox=dict(boxstyle="round,pad=0.28", fc="#F8FAFC", ec="#CBD5E1"))
    ax0.grid(True, alpha=0.25)

    # Right: RoPE Relative Distance Kernel & Shift Invariance
    ax1 = fig.add_subplot(gs[1])
    dk_rope = 64
    q0 = np.ones(dk_rope) / np.sqrt(dk_rope)
    k0 = np.ones(dk_rope) / np.sqrt(dk_rope)
    seps = np.arange(-28, 29)
    i_base = 35
    delta = 50
    dots_base = [float(np.dot(rope_rotate(q0, i_base), rope_rotate(k0, i_base + s))) for s in seps]
    dots_shift = [float(np.dot(rope_rotate(q0, i_base + delta), rope_rotate(k0, i_base + delta + s))) for s in seps]
    max_diff = float(np.max(np.abs(np.array(dots_base) - np.array(dots_shift))))

    ax1.plot(seps, dots_base, "-", color=C_BLUE, lw=2.6, label="Original positions (i, j)")
    ax1.plot(seps[::2], dots_shift[::2], "o", color=C_ORANGE, ms=4.8, mfc="white", mew=1.6,
             label="Shifted (i + 50, j + 50) — identical!")
    ax1.axvline(0, color="#94A3B8", ls=":", lw=1.0)

    ax1.text(0.97, 0.94,
             r"$(R_{i+\Delta}\mathbf{q})^T (R_{j+\Delta}\mathbf{k}) = \mathbf{q}^T R_{j-i}\mathbf{k}$" + "\n"
             + f"max |Δ| across shift Δ = 50:\n{max_diff:.1e}  (exact invariance!)",
             transform=ax1.transAxes, ha="right", va="top", fontsize=8.3, color=C_SLATE,
             bbox=dict(boxstyle="round,pad=0.30", fc="#EFF6FF", ec="#93C5FD", alpha=0.95))

    ax1.set_xlim(-28, 28)
    ax1.set_ylim(-0.05, 1.15)
    ax1.set_xlabel("Relative Sequence Separation  (j − i)")
    ax1.set_ylabel(r"Rotary Inner Product  $\langle R_i\mathbf{q},\, R_j\mathbf{k}\rangle$")
    ax1.set_title("RoPE (ESM-2): Exact Relative-Separation Invariance", pad=8)
    ax1.legend(loc="lower center", framealpha=0.95, fontsize=8.1)
    ax1.grid(True, alpha=0.25)

    fig.subplots_adjust(left=0.06, right=0.985, top=0.89, bottom=0.15)
    fig.savefig(OUT_DIR / "attention_scaling_and_rope.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: MLM Contact Emergence — Raw Attention vs Sym+APC vs Precision
# ---------------------------------------------------------------------------

def fig3_mlm_apc_contact_emergence():
    rng = np.random.default_rng(99)
    N = 24
    # Ground-truth binary long-range contacts (|i-j| > 4): antiparallel hairpin (11..16 vs 23..18) + tertiary pack
    truth = np.zeros((N, N), dtype=int)
    signal = np.zeros((N, N))
    # Local + alpha-helix band (1..10)
    for i in range(N):
        signal[i, i] = 0.25
        if i + 1 < N:
            signal[i, i + 1] = signal[i + 1, i] = 0.20
    for i in range(0, 7):
        signal[i, i + 4] = signal[i + 4, i] = 0.34
    # Long-range antiparallel beta-hairpin contacts (|i - j| > 4)
    hairpin_pairs = [(10, 23), (11, 22), (12, 21), (13, 20), (14, 19), (2, 16), (3, 15), (4, 14), (5, 13), (6, 12)]
    for i, j in hairpin_pairs:
        truth[i, j] = truth[j, i] = 1
        signal[i, j] = signal[j, i] = 0.38

    # Add asymmetric attention + strong phylogenetic/entropy column & row streaks
    col_bias = np.zeros(N)
    col_bias[[3, 7, 17, 21]] = [0.42, 0.38, 0.44, 0.36]
    row_bias = np.zeros(N)
    row_bias[[5, 15, 20]] = [0.28, 0.32, 0.26]
    asym = rng.uniform(0.0, 0.08, size=(N, N))
    A_raw = signal + np.outer(np.ones(N), col_bias) + np.outer(row_bias, np.ones(N)) + asym
    A_raw = A_raw / np.sum(A_raw, axis=1, keepdims=True)  # row-stochastic like real softmax!

    A_sym = 0.5 * (A_raw + A_raw.T)
    A_apc = apc_correct(A_raw)

    def top_l2_precision(M):
        pairs = []
        for i in range(N):
            for j in range(i + 5, N):  # |i - j| > 4
                pairs.append((M[i, j], truth[i, j]))
        pairs.sort(key=lambda x: x[0], reverse=True)
        k = len(hairpin_pairs)  # Top-10 long-range contacts (~L/2.4)
        return float(np.mean([p[1] for p in pairs[:k]]))

    p_raw = top_l2_precision(A_raw)
    p_sym = top_l2_precision(A_sym)
    p_apc = top_l2_precision(A_apc)
    # Simulated 4-head logistic average (combining 4 noisy APC heads)
    heads_apc = [apc_correct(signal + np.outer(np.ones(N), col_bias) + rng.uniform(0, 0.11, (N, N))) for _ in range(4)]
    p_lr = top_l2_precision(np.mean(heads_apc, axis=0))

    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.0, 1.0, 1.08], wspace=0.28)

    # Panel 1: Raw Asymmetric Attention
    ax0 = fig.add_subplot(gs[0])
    ax0.imshow(A_raw, origin="lower", cmap="viridis")
    ax0.set_title("1. Raw Attention A\n(Buried in Entropy Streaks)", fontsize=10.2, pad=7)
    ax0.set_xlabel("Key Residue j")
    ax0.set_ylabel("Query Residue i")
    ax0.set_xticks([0, 11, 23]); ax0.set_xticklabels(["1", "12", "24"])
    ax0.set_yticks([0, 11, 23]); ax0.set_yticklabels(["1", "12", "24"])

    # Panel 2: Symmetrized + APC
    ax1 = fig.add_subplot(gs[1])
    ax1.imshow(np.maximum(A_apc, 0), origin="lower", cmap="viridis")
    # Overlay subtle red rings on true long-range contacts
    ti, tj = np.where(np.triu(truth, k=5) == 1)
    ax1.scatter(tj, ti, s=26, facecolors="none", edgecolors="#F87171", lw=1.2, label="True 3D contact")
    ax1.set_title(r"2. Symmetrized + APC ($A^{\mathrm{APC}}$)" + "\n(3D Contacts Pop Out!)", fontsize=10.2, pad=7)
    ax1.set_xlabel("Residue j")
    ax1.set_ylabel("Residue i")
    ax1.set_xticks([0, 11, 23]); ax1.set_xticklabels(["1", "12", "24"])
    ax1.set_yticks([0, 11, 23]); ax1.set_yticklabels(["1", "12", "24"])
    ax1.legend(loc="lower right", framealpha=0.92, fontsize=7.6)

    # Panel 3: Top-L/2 Precision Jump
    ax2 = fig.add_subplot(gs[2])
    stages = ["Raw\nAttn A", "Symm.\nĀ", "APC\nDenoised", "Multi-Head\nLogReg"]
    precs = [p_raw * 100, p_sym * 100, p_apc * 100, p_lr * 100]
    cols = [C_RED, C_ORANGE, C_BLUE, C_GREEN]
    bars = ax2.bar(np.arange(4), precs, color=cols, width=0.58, edgecolor=C_SLATE, lw=0.9)
    for b, v in zip(bars, precs):
        ax2.text(b.get_x() + b.get_width() / 2, v + 2.2, f"{v:.0f}%",
                 ha="center", va="bottom", fontsize=9.0, fontweight="bold", color=C_SLATE)
    ax2.set_xticks(np.arange(4))
    ax2.set_xticklabels(stages, fontsize=8.3)
    ax2.set_ylabel("Long-Range Contact Precision (%)")
    ax2.set_ylim(0, 118)
    ax2.set_title("3. Long-Range Precision (|i − j| > 4)", fontsize=10.2, pad=7)
    ax2.grid(True, axis="y", alpha=0.25)

    fig.subplots_adjust(left=0.055, right=0.985, top=0.84, bottom=0.16)
    fig.savefig(OUT_DIR / "mlm_apc_contact_emergence.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: Zero-Shot Masked vs WT Marginal & DMS Non-Linear Saturation
# ---------------------------------------------------------------------------

def fig4_zeroshot_and_dms():
    rng = np.random.default_rng(2026)

    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.05, 1.05], wspace=0.26)

    # Left: Masked-Marginal vs WT-Marginal probability profile at a buried hydrophobic site (WT = L)
    ax0 = fig.add_subplot(gs[0])
    aas = ["L*", "I", "V", "M", "F", "A", "Y", "W", "T", "S", "G", "P", "K", "E", "D"]
    # Context-only masked logits vs WT-self-conditioned logits
    z_masked = np.array([2.6, 2.2, 1.9, 1.4, 1.1, 0.3, 0.1, -0.4, -0.9, -1.3, -1.9, -3.2, -3.5, -3.9, -4.2])
    z_wtcond = z_masked.copy()
    z_wtcond[0] += 2.1  # Self-conditioning boost when WT token is NOT masked!
    p_masked = np.exp(z_masked - np.max(z_masked)); p_masked /= p_masked.sum()
    p_wtcond = np.exp(z_wtcond - np.max(z_wtcond)); p_wtcond /= p_wtcond.sum()

    x = np.arange(len(aas))
    w = 0.38
    ax0.bar(x - w / 2, p_masked, width=w, color=C_BLUE, edgecolor=C_SLATE, lw=0.6, label=r"Masked-Marginal $P(a \mid \mathbf{x}_{\setminus i})$")
    ax0.bar(x + w / 2, p_wtcond, width=w, color=C_ORANGE, edgecolor=C_SLATE, lw=0.6, alpha=0.88, label="WT-Marginal (unmasked, self-biased)")
    ax0.annotate("WT self-conditioning\ninflates P(WT)!", xy=(0 + w / 2, p_wtcond[0]), xytext=(3.6, 0.68),
                 fontsize=8.1, fontweight="bold", color="#9A3412",
                 arrowprops=dict(arrowstyle="->", color=C_ORANGE, lw=1.3),
                 bbox=dict(boxstyle="round,pad=0.24", fc="#FFF7ED", ec="#FDBA74"))
    ax0.text(0.97, 0.56,
             r"$\Delta s = \log P(\mathrm{mut}) - \log P(\mathrm{wt})$" + "\n"
             + r"$\equiv z_{\mathrm{mut}} - z_{\mathrm{wt}}$  ($\log Z_i$ cancels!)",
             transform=ax0.transAxes, ha="right", va="top", fontsize=8.3, color=C_SLATE,
             bbox=dict(boxstyle="round,pad=0.26", fc="#F8FAFC", ec="#CBD5E1"))

    ax0.set_xticks(x)
    ax0.set_xticklabels(aas, fontsize=8.4)
    ax0.set_xlabel("Candidate Amino Acid at Buried Core Site (L* = Wild Type)")
    ax0.set_ylabel("Predicted Probability P(a | context)")
    ax0.set_ylim(0, 0.88)
    ax0.set_title("Masked-Marginal vs. Unmasked WT-Marginal Scoring", pad=8)
    ax0.legend(loc="upper right", framealpha=0.95, fontsize=8.0)
    ax0.grid(True, axis="y", alpha=0.25)

    # Right: Synthetic DMS nonlinear Boltzmann/Hill activity curve where Spearman rho >> Pearson r
    ax1 = fig.add_subplot(gs[1])
    esm_s = np.linspace(-9.2, 1.2, 85) + rng.normal(0, 0.25, 85)
    # Cooperative monotonic fluorescence response (strongly compressed dark tail + sharp bright rise)
    true_f = np.exp(0.62 * (esm_s - 0.8)) / (1.0 + 0.15 * np.exp(0.62 * (esm_s - 0.8)))
    true_f = true_f / np.max(true_f)
    obs_f = np.clip(true_f + rng.normal(0, 0.025, size=len(esm_s)), 0.0, 1.08)

    r_val, _ = pearsonr(obs_f, esm_s)
    rho_val, _ = spearmanr(obs_f, esm_s)

    colors = np.where(obs_f < 0.20, C_RED, np.where(obs_f < 0.65, C_ORANGE, C_GREEN))
    ax1.scatter(obs_f, esm_s, c=colors, s=30, alpha=0.84, edgecolors="white", lw=0.5, zorder=3)

    s_grid = np.linspace(-9.4, 1.3, 200)
    f_raw = np.exp(0.62 * (s_grid - 0.8)) / (1.0 + 0.15 * np.exp(0.62 * (s_grid - 0.8)))
    f_grid = f_raw / np.max(f_raw)
    ax1.plot(f_grid, s_grid, "-", color=C_PURPLE, lw=2.3, label="Non-linear fluorescence curve", zorder=4)
    coef = np.polyfit(obs_f, esm_s, 1)
    x_lin = np.linspace(0.0, 1.06, 50)
    ax1.plot(x_lin, np.polyval(coef, x_lin), "--", color=C_SLATE, lw=1.6, label="Linear OLS fit", zorder=2)

    ax1.text(0.97, 0.05,
             f"Pearson linear r  = {r_val:.2f}\nSpearman rank ρ = {rho_val:.2f}\n\n"
             "Rank ρ captures monotonic\nnon-linear fitness cliffs!",
             transform=ax1.transAxes, ha="right", va="bottom", fontsize=8.3, fontweight="bold", color=C_SLATE,
             bbox=dict(boxstyle="round,pad=0.28", fc="#F5F3FF", ec="#DDD6FE", alpha=0.96))

    ax1.set_xlim(-0.06, 1.12)
    ax1.set_ylim(-9.8, 2.0)
    ax1.set_xlabel("Experimental DMS Fluorescence Fitness (0 = dark, 1 = WT)")
    ax1.set_ylabel("Zero-Shot ESM Log-Odds Score Δs")
    ax1.set_title("DMS Validation: Non-Linear Curves Favor Spearman ρ", pad=8)
    ax1.legend(loc="upper left", framealpha=0.95, fontsize=8.0)
    ax1.grid(True, alpha=0.25)

    fig.subplots_adjust(left=0.06, right=0.985, top=0.89, bottom=0.15)
    fig.savefig(OUT_DIR / "zeroshot_marginal_and_dms_metrics.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5: ESMFold vs AlphaFold 2 — Memory Chunking & MSA-Depth Trade-Off
# ---------------------------------------------------------------------------

def fig5_esmfold_and_scaling():
    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.05, 1.05], wspace=0.26)

    # Left: GPU Memory Scaling vs Sequence Length L (GB VRAM in fp16)
    ax0 = fig.add_subplot(gs[0])
    L = np.linspace(100, 950, 200)
    # Realistic ESMFold/AF2 memory footprint components (model weights ~6.5 GB fp16 + activation scaling)
    base_w = 6.5
    vram_unchunk = base_w + 1.05e-7 * (L ** 3) + 1.4e-5 * (L ** 2)
    vram_chunk64 = base_w + 64 * 1.05e-7 * (L ** 2) + 0.85e-5 * (L ** 2)
    vram_attn_only = 0.6 + 0.7e-5 * (L ** 2)  # ESM-2 LM only (no folding trunk)

    ax0.plot(L, vram_unchunk, "-", color=C_RED, lw=2.4, label=r"Unchunked Triangle Updates $O(L^3)$")
    ax0.plot(L, vram_chunk64, "-", color=C_BLUE, lw=2.4, label=r"Chunked Trunk (chunk_size=64) $O(L^2)$")
    ax0.plot(L, vram_attn_only, "--", color=C_TEAL, lw=2.1, label=r"ESM-2 Self-Attention Only $O(L^2)$")

    ax0.axhline(16.0, color=C_PURPLE, ls="-.", lw=1.6)
    ax0.text(935, 14.9, "Free Colab T4 GPU Limit (16 GB VRAM)",
             ha="right", va="top", fontsize=8.0, fontweight="bold", color=C_PURPLE,
             bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#DDD6FE", alpha=0.92))
    ax0.annotate("OOM Crash!\n(L ≈ 440 unchunked)", xy=(445, 16.0), xytext=(625, 22.8),
                 fontsize=8.1, fontweight="bold", color=C_RED,
                 arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.3),
                 bbox=dict(boxstyle="round,pad=0.24", fc="#FEF2F2", ec="#FECACA"))

    ax0.set_xlim(100, 950)
    ax0.set_ylim(0, 28)
    ax0.set_xlabel("Protein Sequence Length L (residues)")
    ax0.set_ylabel("Peak GPU Memory Footprint (GB VRAM)")
    ax0.set_title("Why Chunking & Streaming Enable 3B Folding on T4 GPUs", pad=8)
    ax0.legend(loc="upper left", bbox_to_anchor=(0.01, 0.99), framealpha=0.95, fontsize=7.8)
    ax0.grid(True, alpha=0.25)

    # Right: Accuracy (TM-Score) vs MSA Depth (N_eff) — AlphaFold 2 vs ESMFold
    ax1 = fig.add_subplot(gs[1])
    neff = np.logspace(0, 3, 120)
    log_n = np.log10(neff)
    tm_af2 = 0.44 + 0.47 / (1.0 + np.exp(-2.1 * (log_n - 0.95)))
    tm_esm = 0.56 + 0.28 / (1.0 + np.exp(-1.8 * (log_n - 1.05)))

    ax1.axvspan(1.0, 10.0, color="#DBEAFE", alpha=0.40, zorder=0)
    ax1.text(3.1, 0.34, "Orphan / Shallow MSA\n(ESMFold excels!)", ha="center", va="bottom",
             fontsize=8.0, fontweight="bold", color="#1E3A8A")

    ax1.plot(neff, tm_af2, "-", color=C_PURPLE, lw=2.5, label="AlphaFold 2 (Explicit MSA · minutes/prot)")
    ax1.plot(neff, tm_esm, "-", color=C_GREEN, lw=2.5, label="ESMFold 3B (Single Seq · ~3 sec, 60× faster)")
    ax1.axhline(0.50, color=C_SLATE, ls=":", lw=1.2)
    ax1.text(920, 0.515, "TM = 0.50 (same fold)", ha="right", va="bottom", fontsize=7.9, color=C_SLATE)

    ax1.annotate("Deep homolog family:\nexplicit MSA gives\n+0.05–0.07 TM edge",
                 xy=(520, 0.90), xytext=(15, 0.95),
                 fontsize=7.9, fontweight="bold", color=C_PURPLE, ha="center", va="top",
                 arrowprops=dict(arrowstyle="->", color=C_PURPLE, lw=1.2, connectionstyle="arc3,rad=-0.08"),
                 bbox=dict(boxstyle="round,pad=0.24", fc="#F5F3FF", ec="#DDD6FE"))

    ax1.set_xscale("log")
    ax1.set_xlim(1, 1000)
    ax1.set_ylim(0.30, 1.00)
    ax1.set_xlabel(r"Effective Evolutionary Family Depth $N_{\mathrm{eff}}$ (log scale)")
    ax1.set_ylabel("Expected Structural Accuracy (TM-Score)")
    ax1.set_title("ESMFold vs. AlphaFold 2: Speed vs. Evolutionary Depth", pad=8)
    ax1.legend(loc="lower right", framealpha=0.95, fontsize=8.0)
    ax1.grid(True, alpha=0.25)

    fig.subplots_adjust(left=0.06, right=0.985, top=0.89, bottom=0.15)
    fig.savefig(OUT_DIR / "esmfold_vs_af2_and_scaling.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    print("Generating Lecture 8 Part 3 figures...")
    fig1_epistasis_and_dca()
    fig2_attention_scaling_and_rope()
    fig3_mlm_apc_contact_emergence()
    fig4_zeroshot_and_dms()
    fig5_esmfold_and_scaling()
    print(f"All 5 figures generated cleanly in {OUT_DIR}")
