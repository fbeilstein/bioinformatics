#!/usr/bin/env python3
"""
Generate 5 reproducible, publication-grade 200-DPI figures for
Lecture 7, Part 3: Underlying Algorithms (Folding the Molecular Machine).

Leaves all 4 existing Part 4 figures untouched:
  - bfactor_mobility.png
  - distance_matrix_patterns.png
  - kabsch_superposition.png
  - ramachandran_plot.png

Creates 5 new figures in lecture_07_folding_the_molecular_machine/images/:
  1. dihedral_and_ramachandran_sterics.png
  2. kabsch_se3_and_chirality.png
  3. correspondence_and_tm_score.png
  4. distance_geometry_and_mds.png
  5. debye_waller_and_diffraction.png
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Polygon, FancyArrowPatch, Circle
from scipy.spatial.distance import cdist

OUT_DIR = Path(__file__).resolve().parent / "images"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.edgecolor": "#475569",
    "axes.linewidth": 0.9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.6,
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
# Exact Geometry Helpers (also reused for numerical cross-checks)
# ---------------------------------------------------------------------------

def place_atom(a, b, c, bond_len, bond_angle_deg, torsion_deg):
    """NeRF (Natural Extension Reference Frame) placement of atom D given A, B, C."""
    bc = c - b
    bc_u = bc / np.linalg.norm(bc)
    ab = b - a
    n = np.cross(ab, bc_u)
    n_norm = np.linalg.norm(n)
    if n_norm < 1e-12:
        # Degenerate fallback
        n = np.cross(np.array([1.0, 0.0, 0.0]), bc_u)
        if np.linalg.norm(n) < 1e-6:
            n = np.cross(np.array([0.0, 1.0, 0.0]), bc_u)
        n_norm = np.linalg.norm(n)
    n_u = n / n_norm
    m_u = np.cross(n_u, bc_u)

    th = np.radians(180.0 - bond_angle_deg)
    ph = np.radians(torsion_deg)
    d_local = np.array([
        bond_len * np.cos(th),
        bond_len * np.sin(th) * np.cos(ph),
        bond_len * np.sin(th) * np.sin(ph),
    ])
    return c + d_local[0] * bc_u + d_local[1] * m_u + d_local[2] * n_u


def dihedral_deg(r1, r2, r3, r4):
    """Signed IUPAC dihedral angle in (-180, 180] via plane normals and atan2."""
    b1 = r2 - r1
    b2 = r3 - r2
    b3 = r4 - r3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    u = b2 / np.linalg.norm(b2)
    return float(np.degrees(np.arctan2(np.dot(np.cross(n1, n2), u), np.dot(n1, n2))))


def build_dipeptide_coords(phi_deg, psi_deg, include_cb=True):
    """
    Build exact Engh-Huber Ace-Xaa-NMe capped dipeptide around (phi, psi).
    Atoms: [CA_prev, C_prev, O_prev, N, H, CA, HA, C, O, N_next, H_next, CA_next, (CB)]
    """
    # Place rigid preceding peptide plane in xy-plane with omega = 180 deg
    n_pos = np.array([0.0, 0.0, 0.0])
    ca_pos = np.array([1.458, 0.0, 0.0])
    ang_cnc = np.radians(121.7)
    c_prev = np.array([1.329 * np.cos(ang_cnc), 1.329 * np.sin(ang_cnc), 0.0])
    # H on N opposite to C_prev-N-CA exterior
    h_pos = place_atom(ca_pos, c_prev, n_pos, 1.01, 119.0, 180.0)
    # Use proper N-centered geometry for H and O_prev:
    # C_prev -> N -> CA lies in z=0; O_prev is trans to H (same side as CA across C_prev-N? Wait:
    # in trans peptide: CA_prev-C_prev-N-CA = 180, O_prev-C_prev-N-CA = 0 (O_prev cis to CA, trans to H!)
    # Check: O=C-N-H has torsion 180 deg, so O is cis to CA (torsion O-C_prev-N-CA = 0 deg) and H is cis to CA_prev.
    ca_prev = place_atom(ca_pos, n_pos, c_prev, 1.525, 116.2, 180.0)
    o_prev = place_atom(ca_pos, n_pos, c_prev, 1.231, 123.0, 0.0)
    h_pos = place_atom(o_prev, c_prev, n_pos, 1.01, 119.5, 180.0)

    # C placed by phi = dihedral(C_prev, N, CA, C)
    c_pos = place_atom(c_prev, n_pos, ca_pos, 1.525, 111.2, phi_deg)

    # L-amino acid tetrahedral substituents on CA:
    # Viewed along CA -> HA, N -> C -> CB goes clockwise (L-configuration: dihedral C_prev-N-CA-CB = phi - 122.7 deg)
    cb_pos = place_atom(c_prev, n_pos, ca_pos, 1.521, 110.5, phi_deg - 122.7)
    ha_pos = place_atom(c_prev, n_pos, ca_pos, 1.09, 109.5, phi_deg + 118.5)

    # Next peptide plane placed by psi = dihedral(N, CA, C, N_next)
    n_next = place_atom(n_pos, ca_pos, c_pos, 1.329, 116.2, psi_deg)
    o_pos = place_atom(n_pos, ca_pos, c_pos, 1.231, 120.8, psi_deg + 180.0)
    h_next = place_atom(ca_pos, c_pos, n_next, 1.01, 119.5, 180.0)
    ca_next = place_atom(ca_pos, c_pos, n_next, 1.458, 121.7, 180.0)

    atoms = [
        ("C", ca_prev), ("C", c_prev), ("O", o_prev),
        ("N", n_pos), ("H", h_pos), ("C", ca_pos), ("H", ha_pos),
        ("C", c_pos), ("O", o_pos), ("N", n_next), ("H", h_next), ("C", ca_next),
    ]
    if include_cb:
        atoms.append(("CB", cb_pos))
    return atoms


def ramachandran_steric_grid(n_grid=121):
    """
    Compute exact hard-sphere steric clash maps for L-Alanine and Glycine
    using Ramachandran (1963) outer-limit & normal contact radii.
    """
    # Non-bonded pairs that vary with (phi, psi)
    # Indices: 0:CA_prev, 1:C_prev, 2:O_prev, 3:N, 4:H, 5:CA, 6:HA, 7:C, 8:O, 9:N_next, 10:H_next, 11:CA_next, 12:CB
    # Exclude 1-2 (bonded) and 1-3 (bond-angle fixed) pairs; keep 1-4+ pairs that depend on phi/psi
    pairs_gly = [
        (1, 7), (1, 8), (1, 9), (1, 10),
        (2, 7), (2, 8), (2, 9), (2, 10), (2, 6),
        (0, 7), (0, 8), (0, 9),
        (3, 8), (3, 9), (3, 10), (3, 11),
        (4, 7), (4, 8), (4, 9), (4, 10), (4, 11),
        (5, 10), (5, 11),
        (6, 8), (6, 9), (6, 10), (6, 11),
    ]
    pairs_cb = [
        (12, 0), (12, 1), (12, 2), (12, 4),
        (12, 8), (12, 9), (12, 10), (12, 11),
    ]

    # Classic Ramachandran (1963) normal & outer-limit hard-sphere radii (Å)
    r_norm = {"C": 1.49, "CB": 1.56, "N": 1.38, "O": 1.33, "H": 0.96}
    r_outer = {"C": 1.43, "CB": 1.48, "N": 1.31, "O": 1.25, "H": 0.90}

    angles = np.linspace(-180.0, 180.0, n_grid)
    ala_state = np.zeros((n_grid, n_grid), dtype=int)  # 0=clash, 1=outer, 2=core
    gly_state = np.zeros((n_grid, n_grid), dtype=int)

    for iy, psi in enumerate(angles):
        for ix, phi in enumerate(angles):
            atoms = build_dipeptide_coords(phi, psi, include_cb=True)
            coords = np.array([a[1] for a in atoms])
            types = [a[0] for a in atoms]

            # Glycine check
            gly_core = True
            gly_outer = True
            for i, j in pairs_gly:
                d = np.linalg.norm(coords[i] - coords[j])
                hb_relax = 0.16 if {types[i], types[j]} == {"H", "O"} else 0.0
                if d < r_outer[types[i]] + r_outer[types[j]] - hb_relax:
                    gly_outer = False
                    gly_core = False
                    break
                if d < r_norm[types[i]] + r_norm[types[j]] - hb_relax:
                    gly_core = False
            gly_state[iy, ix] = 2 if gly_core else (1 if gly_outer else 0)

            # L-Alanine check adds CB pairs
            ala_core = gly_core
            ala_outer = gly_outer
            if ala_outer:
                for i, j in pairs_cb:
                    d = np.linalg.norm(coords[i] - coords[j])
                    if d < r_outer[types[i]] + r_outer[types[j]]:
                        ala_outer = False
                        ala_core = False
                        break
                    if d < r_norm[types[i]] + r_norm[types[j]]:
                        ala_core = False
            ala_state[iy, ix] = 2 if ala_core else (1 if ala_outer else 0)

    return angles, ala_state, gly_state


def kabsch_align(P, Q, enforce_so3=True):
    """Exact Kabsch SVD alignment of P onto Q (N x 3 row matrices)."""
    cP = P.mean(axis=0)
    cQ = Q.mean(axis=0)
    Pc = P - cP
    Qc = Q - cQ
    C = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(C)
    V = Vt.T
    d = 1.0
    if enforce_so3 and np.linalg.det(V @ U.T) < 0:
        d = -1.0
    D = np.diag([1.0, 1.0, d])
    R_col = V @ D @ U.T
    P_rot = Pc @ R_col.T
    diff = P_rot - Qc
    rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
    return P_rot + cQ, R_col, rmsd, float(np.linalg.det(V @ U.T))


def tm_score(P_aligned, Q):
    """Exact TM-score on paired N x 3 coordinates."""
    L = len(Q)
    d0 = 1.24 * (max(L - 15, 1) ** (1.0 / 3.0)) - 1.8
    d0 = max(d0, 0.5)
    d = np.linalg.norm(P_aligned - Q, axis=1)
    return float(np.mean(1.0 / (1.0 + (d / d0) ** 2)))


def gdt_ts(P_aligned, Q):
    """Global Distance Test - Total Score (average fraction within 1, 2, 4, 8 Å)."""
    d = np.linalg.norm(P_aligned - Q, axis=1)
    return float(0.25 * sum(np.mean(d < cut) for cut in (1.0, 2.0, 4.0, 8.0)))


def make_alpha_helix(n_res=24):
    """Exact right-handed alpha-helix CA coordinates (r=2.3 Å, rise=1.5 Å, 100°/res)."""
    idx = np.arange(n_res)
    theta = np.radians(100.0 * idx)
    return np.column_stack([2.3 * np.cos(theta), 2.3 * np.sin(theta), 1.5 * idx])


def make_beta_sheet(n_strand=12, parallel=True):
    """Two-strand beta-sheet CA coordinates (pleated z-ripple, 3.3 Å step, 4.8 Å inter-strand)."""
    i1 = np.arange(n_strand)
    s1 = np.column_stack([3.3 * i1, np.zeros(n_strand), 0.7 * ((-1) ** i1)])
    if parallel:
        # Both strands run N -> C in +x direction (connected via crossover return)
        s2 = np.column_stack([3.3 * i1, np.full(n_strand, 4.8), 0.7 * ((-1) ** i1)])
    else:
        # Anti-parallel hairpin: strand 2 runs C -> N in -x direction
        s2 = np.column_stack([3.3 * (n_strand - 1 - i1), np.full(n_strand, 4.8), 0.7 * ((-1) ** (n_strand - 1 - i1))])
    return np.vstack([s1, s2])


# ---------------------------------------------------------------------------
# Figure 1: Dihedral Geometry & Analytical Ramachandran Sterics
# ---------------------------------------------------------------------------

def fig1_dihedral_ramachandran():
    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.05, 1.05], wspace=0.22)

    # --- Left Panel: Peptide Plane & 4-Atom Dihedral Schematic ---
    ax = fig.add_subplot(gs[0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    ax.set_title("Peptide Plane Resonance (ω ≈ 180°) & Signed Torsion (φ, ψ)", pad=8)

    # Top schematic: Two rigid amide planes hinged at C_alpha
    plane1 = Polygon([(0.55, 3.45), (3.65, 3.45), (4.15, 5.45), (1.05, 5.45)],
                     closed=True, facecolor="#DBEAFE", edgecolor=C_BLUE, lw=1.4, alpha=0.85, zorder=1)
    plane2 = Polygon([(5.55, 3.45), (8.65, 3.45), (9.35, 5.35), (6.25, 5.35)],
                     closed=True, facecolor="#EDE9FE", edgecolor=C_PURPLE, lw=1.4, alpha=0.85, zorder=1)
    ax.add_patch(plane1)
    ax.add_patch(plane2)

    ax.text(2.35, 5.65, "Peptide Plane 1 (ω ≈ 180°, rigid)", ha="center", va="bottom",
            fontsize=8.6, fontweight="bold", color=C_BLUE)
    ax.text(7.55, 5.65, "Peptide Plane 2 (ω ≈ 180°, rigid)", ha="center", va="bottom",
            fontsize=8.6, fontweight="bold", color=C_PURPLE)

    # Atom coordinates in diagram space
    pts = {
        "Cm": (1.55, 4.45), "Om": (1.45, 5.10),
        "N":  (3.15, 4.30), "H":  (3.25, 3.72),
        "CA": (4.85, 4.45), "CB": (4.85, 5.45),
        "C":  (6.55, 4.35), "O":  (6.65, 5.05),
        "Np": (8.15, 4.40), "Hp": (8.05, 3.75),
    }
    bonds = [
        ("Cm", "Om", C_SLATE, 2.2), ("Cm", "N", C_BLUE, 2.6),
        ("N", "H", C_SLATE, 1.5), ("N", "CA", C_ORANGE, 2.8),
        ("CA", "CB", C_RED, 2.0), ("CA", "C", C_TEAL, 2.8),
        ("C", "O", C_SLATE, 2.2), ("C", "Np", C_PURPLE, 2.6),
        ("Np", "Hp", C_SLATE, 1.5),
    ]
    for a, b, col, lw in bonds:
        ax.plot([pts[a][0], pts[b][0]], [pts[a][1], pts[b][1]], color=col, lw=lw, zorder=2)

    atom_style = {
        "Cm": (r"$C_{i-1}$", "#475569", 0.27),
        "Om": ("O", C_RED, 0.20),
        "N":  (r"$N_i$", C_BLUE, 0.27),
        "H":  ("H", "#94A3B8", 0.17),
        "CA": (r"$C\alpha_i$", C_SLATE, 0.31),
        "CB": (r"$C\beta$", C_RED, 0.24),
        "C":  (r"$C_i$", "#475569", 0.27),
        "O":  ("O", C_RED, 0.20),
        "Np": (r"$N_{i+1}$", C_PURPLE, 0.27),
        "Hp": ("H", "#94A3B8", 0.17),
    }
    for k, (lbl, col, rad) in atom_style.items():
        x, y = pts[k]
        ax.add_patch(Circle((x, y), rad, facecolor=col, edgecolor="white", lw=1.2, zorder=4))
        ax.text(x, y, lbl, ha="center", va="center", color="white",
                fontsize=7.3 if len(lbl) > 2 else 7.8, fontweight="bold", zorder=5)

    # Curved rotation arrows for phi and psi
    ax.annotate("", xy=(4.15, 4.78), xytext=(3.85, 4.02),
                arrowprops=dict(arrowstyle="->", color=C_ORANGE, lw=2.0, connectionstyle="arc3,rad=-0.45"), zorder=6)
    ax.text(4.02, 3.32, r"$\phi$ ($N_i\!-\!C\alpha_i$)", ha="center", va="top",
            fontsize=9.0, fontweight="bold", color=C_ORANGE)

    ax.annotate("", xy=(5.85, 4.78), xytext=(5.55, 4.02),
                arrowprops=dict(arrowstyle="->", color=C_TEAL, lw=2.0, connectionstyle="arc3,rad=0.45"), zorder=6)
    ax.text(5.72, 3.32, r"$\psi$ ($C\alpha_i\!-\!C_i$)", ha="center", va="top",
            fontsize=9.0, fontweight="bold", color=C_TEAL)

    ax.text(4.85, 5.82, "L-side chain\n(absent in Gly)", ha="center", va="bottom",
            fontsize=7.6, color=C_RED, style="italic")

    # Lower schematic: 4-atom bond vectors & plane normals for atan2
    box = FancyBboxPatch((0.35, 0.20), 9.30, 2.65, boxstyle="round,pad=0.12",
                         facecolor="#F8FAFC", edgecolor="#CBD5E1", lw=1.1, zorder=1)
    ax.add_patch(box)
    ax.text(0.60, 2.66, "Signed 4-Atom Torsion via Half-Plane Normals:",
            fontsize=8.9, fontweight="bold", color=C_SLATE, va="top")

    # 4-atom chain schematic on lower left
    r1, r2, r3, r4 = (0.85, 0.95), (2.00, 1.58), (3.95, 1.58), (4.82, 0.70)
    # Shaded half-planes
    ax.add_patch(Polygon([r1, r2, r3, (2.80, 0.95)], facecolor="#BFDBFE", alpha=0.45, edgecolor=C_BLUE, lw=0.8, zorder=2))
    ax.add_patch(Polygon([r2, r3, r4, (2.88, 0.70)], facecolor="#FDE68A", alpha=0.45, edgecolor=C_ORANGE, lw=0.8, zorder=2))
    for p_a, p_b, lbl, offset in [
        (r1, r2, r"$\mathbf{b}_1$", (-0.25, 0.14)),
        (r2, r3, r"$\mathbf{b}_2\ (\hat{\mathbf{u}})$", (-0.38, -0.22)),
        (r3, r4, r"$\mathbf{b}_3$", (0.28, 0.08)),
    ]:
        ax.annotate("", xy=p_b, xytext=p_a,
                    arrowprops=dict(arrowstyle="-|>", color=C_SLATE, lw=2.0), zorder=4)
        mx, my = 0.5 * (p_a[0] + p_b[0]) + offset[0], 0.5 * (p_a[1] + p_b[1]) + offset[1]
        ax.text(mx, my, lbl, ha="center", va="center", fontsize=8.3, fontweight="bold", color=C_SLATE)

    for pt, name, dy in [(r1, r"$\mathbf{r}_1$", -0.22), (r2, r"$\mathbf{r}_2$", 0.22),
                         (r3, r"$\mathbf{r}_3$", 0.22), (r4, r"$\mathbf{r}_4$", -0.20)]:
        ax.add_patch(Circle(pt, 0.11, facecolor=C_BLUE, edgecolor="white", lw=1.0, zorder=5))
        ax.text(pt[0], pt[1] + dy, name, ha="center", va="center", fontsize=7.9, color=C_MUTED)

    # Normal vectors emanating from midpoint of b2
    ax.annotate("", xy=(2.98, 2.30), xytext=(2.98, 1.58),
                arrowprops=dict(arrowstyle="-|>", color=C_BLUE, lw=1.8), zorder=5)
    ax.text(2.55, 2.22, r"$\mathbf{n}_1$", fontsize=8.7, fontweight="bold", color=C_BLUE)
    ax.annotate("", xy=(3.62, 2.18), xytext=(2.98, 1.58),
                arrowprops=dict(arrowstyle="-|>", color=C_ORANGE, lw=1.8), zorder=5)
    ax.text(3.78, 2.20, r"$\mathbf{n}_2$", fontsize=8.7, fontweight="bold", color=C_ORANGE)

    # Formula block on lower right (100% inside box!)
    ax.text(5.28, 2.18,
            r"$\mathbf{n}_1 = \mathbf{b}_1 \times \mathbf{b}_2,\quad \mathbf{n}_2 = \mathbf{b}_2 \times \mathbf{b}_3$" + "\n"
            + r"$\theta = \mathrm{atan2}\!\left((\mathbf{n}_1 \times \mathbf{n}_2)\cdot\hat{\mathbf{u}},\ \mathbf{n}_1\cdot\mathbf{n}_2\right)$",
            fontsize=8.5, color=C_SLATE, va="top", linespacing=1.6)
    ax.text(5.28, 0.86,
            "Preserves handedness sign in (−180°, +180°]\n(unlike arccos, which collapses +θ and −θ)",
            fontsize=7.8, color=C_MUTED, va="top", linespacing=1.3)

    # --- Right Panel: Exact Hard-Sphere Ramachandran Map ---
    ax2 = fig.add_subplot(gs[1])
    angles, ala_state, gly_state = ramachandran_steric_grid(n_grid=121)
    X, Y = np.meshgrid(angles, angles)

    ax2.set_facecolor("#F8FAFC")
    ax2.contourf(X, Y, gly_state, levels=[0.5, 2.5], colors=["#FFEDD5"], alpha=0.88)
    ax2.contour(X, Y, gly_state, levels=[0.5], colors=[C_ORANGE], linewidths=1.3, linestyles="--")

    ax2.contourf(X, Y, ala_state, levels=[0.5, 1.5, 2.5], colors=["#93C5FD", "#2563EB"], alpha=0.90)
    ax2.contour(X, Y, ala_state, levels=[0.5], colors=["#1D4ED8"], linewidths=1.1)

    ax2.axvspan(-75, -50, color=C_PURPLE, alpha=0.13, zorder=2)
    ax2.axvline(-63, color=C_PURPLE, ls=":", lw=1.3, zorder=3)

    ax2.axhline(0, color="#94A3B8", lw=0.8, ls="-")
    ax2.axvline(0, color="#94A3B8", lw=0.8, ls="-")

    ax2.scatter([-135, -60, 60], [140, -45, 45], s=42, c=["white", "white", C_ORANGE],
                edgecolors=[C_SLATE, C_SLATE, "white"], lw=1.4, zorder=6)

    ax2.annotate(r"$\beta$-sheet" + "\n(−135°, +140°)", xy=(-135, 140), xytext=(-172, 55),
                 fontsize=8.4, fontweight="bold", color="#1E3A8A",
                 arrowprops=dict(arrowstyle="->", color="#1E3A8A", lw=1.1),
                 bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#93C5FD", alpha=0.94))
    ax2.annotate(r"Right-handed $\alpha$" + "\n(−60°, −45°)", xy=(-60, -45), xytext=(-172, -145),
                 fontsize=8.4, fontweight="bold", color="#1E3A8A",
                 arrowprops=dict(arrowstyle="->", color="#1E3A8A", lw=1.1),
                 bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#93C5FD", alpha=0.94))
    ax2.annotate(r"Left-handed $\alpha$" + "\n(Glycine only!)", xy=(60, 45), xytext=(55, 118),
                 fontsize=8.4, fontweight="bold", color="#9A3412",
                 arrowprops=dict(arrowstyle="->", color=C_ORANGE, lw=1.2),
                 bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#FDBA74", alpha=0.95))
    ax2.text(-63, -170, "Proline\nφ ≈ −63°", ha="center", va="bottom",
             fontsize=7.9, fontweight="bold", color=C_PURPLE,
             bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#DDD6FE", alpha=0.92))
    ax2.text(105, -112, "Forbidden for L-aa\n(Cβ···O steric clash)", ha="center", va="center",
             fontsize=8.2, color=C_RED, style="italic",
             bbox=dict(boxstyle="round,pad=0.22", fc="#FEF2F2", ec="#FECACA", alpha=0.94))

    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor="#2563EB", edgecolor="#1D4ED8", label="L-Ala core"),
        Patch(facecolor="#93C5FD", edgecolor="#1D4ED8", label="L-Ala outer"),
        Patch(facecolor="#FFEDD5", edgecolor=C_ORANGE, ls="--", label="Glycine extra"),
    ]
    # Place legend neatly in the empty central white forbidden corridor (phi in [-35, +50])
    ax2.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.48, 0.99),
               framealpha=0.95, fontsize=7.6, borderpad=0.35)

    ax2.set_xlim(-180, 180)
    ax2.set_ylim(-180, 180)
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax2.set_yticks([-180, -120, -60, 0, 60, 120, 180])
    ax2.set_xlabel(r"$\phi$ (degrees)")
    ax2.set_ylabel(r"$\psi$ (degrees)")
    ax2.set_title("Analytical Hard-Sphere Steric Map (Engh–Huber Geometry)", pad=8)

    fig.subplots_adjust(left=0.03, right=0.98, top=0.90, bottom=0.13)
    fig.savefig(OUT_DIR / "dihedral_and_ramachandran_sterics.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: Kabsch SE(3) Progression & Chirality / Transpose Traps
# ---------------------------------------------------------------------------

def fig2_kabsch_se3_and_chirality():
    rng = np.random.default_rng(42)
    # Build chiral asymmetric helix-turn peptide (18 residues)
    P = make_alpha_helix(18)
    P[12:, 0] += np.linspace(0, 5.5, 6)
    P[12:, 1] += np.linspace(0, 3.2, 6)

    # Apply true SO(3) rotation + translation + mild thermal noise
    th = np.radians(118.0)
    Rz = np.array([[np.cos(th), -np.sin(th), 0],
                   [np.sin(th),  np.cos(th), 0],
                   [0,           0,          1]])
    ph = np.radians(32.0)
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(ph), -np.sin(ph)],
                   [0, np.sin(ph),  np.cos(ph)]])
    R_true = Rz @ Rx
    t_true = np.array([15.0, -7.5, 4.0])
    Q = (P @ R_true.T) + t_true + rng.normal(0, 0.30, size=P.shape)

    raw_rmsd = float(np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1))))
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    cent_rmsd = float(np.sqrt(np.mean(np.sum((Pc - Qc) ** 2, axis=1))))
    P_al, R_col, opt_rmsd, det_same = kabsch_align(P, Q, enforce_so3=True)

    # Transpose trap: Pc @ R_col instead of Pc @ R_col.T
    P_wrong_T = Pc @ R_col
    trap_T_rmsd = float(np.sqrt(np.mean(np.sum((P_wrong_T - Qc) ** 2, axis=1))))

    # Mirror-enantiomer D-peptide trap: reflect P across z -> -z
    Q_mirror = P.copy()
    Q_mirror[:, 2] *= -1.0
    Q_mirror = (Q_mirror @ R_true.T) + t_true
    _, _, mir_so3_rmsd, det_mir = kabsch_align(P, Q_mirror, enforce_so3=True)
    _, _, mir_o3_rmsd, _ = kabsch_align(P, Q_mirror, enforce_so3=False)

    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.05, 1.0, 1.15], wspace=0.28)

    # Panel A: Raw SE(3) frames
    ax0 = fig.add_subplot(gs[0])
    ax0.plot(P[:, 0], P[:, 1], "-o", color=C_BLUE, lw=2.0, ms=4.5, label="Reference P")
    ax0.plot(Q[:, 0], Q[:, 1], "-s", color=C_ORANGE, lw=2.0, ms=4.2, label="Mobile Q (SE(3) shifted)")
    cP, cQ = P.mean(axis=0), Q.mean(axis=0)
    ax0.scatter([cP[0], cQ[0]], [cP[1], cQ[1]], s=75, marker="*", c=[C_BLUE, C_ORANGE], edgecolors="black", zorder=6)
    ax0.annotate("", xy=(cQ[0], cQ[1]), xytext=(cP[0], cP[1]),
                 arrowprops=dict(arrowstyle="-|>", color=C_RED, lw=1.6, ls="--"))
    ax0.text(0.5 * (cP[0] + cQ[0]) - 0.5, 0.5 * (cP[1] + cQ[1]) + 1.3,
             r"$\mathbf{t}^* = \mathbf{c}_Q - R\mathbf{c}_P$", color=C_RED, fontsize=8.5, fontweight="bold", ha="center")
    ax0.set_title("1. Raw Frames (Arbitrary Origin)", pad=8)
    ax0.set_xlabel("X (Å)")
    ax0.set_ylabel("Y (Å)")
    ax0.legend(loc="lower left", framealpha=0.92, fontsize=8.0)
    ax0.text(0.97, 0.96, f"Raw RMSD:\n{raw_rmsd:.1f} Å", transform=ax0.transAxes,
             ha="right", va="top", fontsize=8.8, fontweight="bold", color=C_RED,
             bbox=dict(boxstyle="round,pad=0.25", fc="#FEF2F2", ec="#FECACA"))
    ax0.grid(True, alpha=0.25)

    # Panel B: After Kabsch SVD
    ax1 = fig.add_subplot(gs[1])
    P_al_c = P_al - Q.mean(axis=0)
    ax1.plot(Qc[:, 0], Qc[:, 1], "-s", color=C_ORANGE, lw=2.2, ms=4.5, alpha=0.85, label="Centered Q_c")
    ax1.plot(P_al_c[:, 0], P_al_c[:, 1], "-o", color=C_BLUE, lw=2.0, ms=4.2, alpha=0.9, label=r"$P_c R_{\mathrm{col}}^T$ (Kabsch)")
    for i in range(0, len(Qc), 2):
        ax1.plot([Qc[i, 0], P_al_c[i, 0]], [Qc[i, 1], P_al_c[i, 1]], color=C_GREEN, lw=0.9, alpha=0.6)
    ax1.scatter([0], [0], s=75, marker="*", c="black", zorder=7, label="Shared centroid (0,0,0)")
    ax1.set_title("2. Kabsch SVD Aligned (det R = +1)", pad=8)
    ax1.set_xlabel("Centered X (Å)")
    ax1.set_ylabel("Centered Y (Å)")
    ymin1, ymax1 = ax1.get_ylim()
    ax1.set_ylim(ymin1 - 1.9, ymax1 + 1.8)
    ax1.legend(loc="lower left", framealpha=0.94, fontsize=7.6)
    ax1.text(0.97, 0.96, f"Centered: {cent_rmsd:.1f} Å\nKabsch:   {opt_rmsd:.2f} Å",
             transform=ax1.transAxes, ha="right", va="top", fontsize=8.3, fontweight="bold", color=C_GREEN,
             bbox=dict(boxstyle="round,pad=0.24", fc="#F0FDF4", ec="#BBF7D0"))
    ax1.grid(True, alpha=0.25)

    # Panel C: Failure Modes (Chirality & Transpose Trap)
    ax2 = fig.add_subplot(gs[2])
    labels = [
        "Same L-fold\nKabsch Pc@R.T",
        "Same L-fold\nBug: Pc@R",
        "L vs D-Mirror\nO(3) det=−1",
        "L vs D-Mirror\nSO(3) det=+1",
    ]
    vals = [opt_rmsd, trap_T_rmsd, mir_o3_rmsd, mir_so3_rmsd]
    cols = [C_GREEN, C_RED, C_ORANGE, C_BLUE]
    bars = ax2.bar(np.arange(4), vals, color=cols, width=0.58, edgecolor=C_SLATE, lw=0.8, alpha=0.9)
    for b, v, note in zip(bars, vals, ["True fit", "Rotated\nbackwards!", "ILLEGAL\nmirror flip!", "L-chirality\npreserved"]):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.25, f"{v:.2f} Å\n({note})",
                 ha="center", va="bottom", fontsize=7.6, fontweight="bold", color=C_SLATE)

    ax2.set_xticks(np.arange(4))
    ax2.set_xticklabels(labels, fontsize=7.8)
    ax2.set_ylabel("Computed RMSD (Å)")
    ax2.set_ylim(0, max(vals) * 1.30)
    ax2.set_title("3. Two Kabsch Traps: Transpose & Chirality", pad=8)
    ax2.grid(True, axis="y", alpha=0.25)

    fig.subplots_adjust(left=0.055, right=0.985, top=0.89, bottom=0.17)
    fig.savefig(OUT_DIR / "kabsch_se3_and_chirality.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: Frame-Shift Correspondence Pitfall & RMSD vs TM-Score / GDT-TS
# ---------------------------------------------------------------------------

def fig3_correspondence_and_tm_score():
    rng = np.random.default_rng(7)
    # Build 80-residue synthetic fold (helix + loop + beta-hairpin + helix)
    n_tot = 80
    t = np.linspace(0, 6.5 * np.pi, n_tot)
    core = np.column_stack([
        9.0 * np.cos(t * 0.35) + 2.1 * np.cos(t * 2.1),
        9.0 * np.sin(t * 0.35) + 2.1 * np.sin(t * 2.1),
        0.55 * np.arange(n_tot) + 1.8 * np.sin(t * 0.9),
    ])
    # Slightly perturbed homolog (true RMSD ~ 1.25 Å)
    homolog = core + rng.normal(0, 0.72, size=core.shape)
    # Smooth slightly along sequence to look like real backbone divergence
    for k in range(3):
        homolog[:, k] = np.convolve(homolog[:, k], [0.2, 0.6, 0.2], mode="same")

    shifts = [0, 1, 2, 3, 4]
    shift_rmsds = []
    shift_tms = []
    for sh in shifts:
        if sh == 0:
            p_sub, q_sub = core, homolog
        else:
            p_sub, q_sub = core[:-sh], homolog[sh:]
        p_al, _, r_val, _ = kabsch_align(p_sub, q_sub, enforce_so3=True)
        shift_rmsds.append(r_val)
        shift_tms.append(tm_score(p_al, q_sub))

    # Tail swing experiment: 72-res conserved core + 8-res flexible tail
    swings = np.linspace(0.0, 20.0, 25)
    tail_rmsds, tail_tms, tail_gdts = [], [], []
    base_hom = core + rng.normal(0, 0.40, size=core.shape)
    for sw in swings:
        q_mod = base_hom.copy()
        # Swing last 8 residues outward smoothly
        ramp = np.linspace(0.15, 1.0, 8)
        q_mod[-8:, 0] += sw * ramp
        q_mod[-8:, 1] -= 0.45 * sw * ramp
        p_al, _, r_val, _ = kabsch_align(core, q_mod, enforce_so3=True)
        tail_rmsds.append(r_val)
        tail_tms.append(tm_score(p_al, q_mod))
        tail_gdts.append(gdt_ts(p_al, q_mod))

    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.95, 1.15], wspace=0.28)

    # Left panel: Frame-Shift Pitfall
    ax0 = fig.add_subplot(gs[0])
    cols = [C_GREEN, C_ORANGE, C_RED, "#B91C1C", "#7F1D1D"]
    bars = ax0.bar(shifts, shift_rmsds, color=cols, width=0.58, edgecolor=C_SLATE, lw=0.9)
    for b, rv in zip(bars, shift_rmsds):
        ax0.text(b.get_x() + b.get_width() / 2, rv + 0.18, f"{rv:.2f} Å",
                 ha="center", va="bottom", fontsize=8.8, fontweight="bold", color=C_SLATE)

    ax0.annotate("Naive [:min_len] off-by-2!\n(e.g. 6LU7 1–306 vs 1UK4 3–303:\n1.35 Å → 6.21 Å illusion)",
                 xy=(2, shift_rmsds[2] + 0.14), xytext=(-0.38, max(shift_rmsds) * 0.90),
                 fontsize=8.1, fontweight="bold", color="#991B1B",
                 arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.4, connectionstyle="arc3,rad=-0.12"),
                 bbox=dict(boxstyle="round,pad=0.28", fc="#FEF2F2", ec="#FCA5A5"))

    ax0.set_xticks(shifts)
    ax0.set_xticklabels([f"Shift {s}\n({'Matched' if s == 0 else f'Off by {s}'})" for s in shifts], fontsize=8.5)
    ax0.set_ylabel("Kabsch Backbone RMSD (Å)")
    ax0.set_ylim(0, max(shift_rmsds) * 1.28)
    ax0.set_title("Correspondence Pitfall: Index Frame-Shift Explosion", pad=8)
    ax0.grid(True, axis="y", alpha=0.25)

    # Right panel: RMSD vs TM-score & GDT-TS
    ax1 = fig.add_subplot(gs[1])
    l1, = ax1.plot(swings, tail_rmsds, "-o", color=C_RED, lw=2.3, ms=4.2, label="Global RMSD (Å, left axis)")
    ax1.set_xlabel("Flexible Surface Tail Excursion (Å, 8 of 80 residues)")
    ax1.set_ylabel("Global RMSD (Å)", color=C_RED, fontweight="bold")
    ax1.tick_params(axis="y", labelcolor=C_RED)
    ax1.set_ylim(0, max(tail_rmsds) * 1.18)
    ax1.grid(True, alpha=0.25)

    ax1r = ax1.twinx()
    l2, = ax1r.plot(swings, tail_tms, "-s", color=C_BLUE, lw=2.3, ms=4.2, label="TM-score (right axis)")
    l3, = ax1r.plot(swings, tail_gdts, "--^", color=C_TEAL, lw=2.0, ms=4.0, label="GDT-TS (right axis)")
    ax1r.axhline(0.5, color=C_PURPLE, ls=":", lw=1.4)
    ax1r.text(0.35, 0.52, "TM = 0.50 (same-fold threshold)", ha="left", va="bottom",
              fontsize=8.1, fontweight="bold", color=C_PURPLE)
    ax1r.set_ylabel("Normalized Score (0 – 1)", color=C_BLUE, fontweight="bold")
    ax1r.tick_params(axis="y", labelcolor=C_BLUE)
    ax1r.set_ylim(0.0, 1.08)

    ax1.set_title("RMSD Quadratic Blow-Up vs. TM-Score & GDT-TS Stability", pad=8)
    ax1.legend(handles=[l1, l2, l3], loc="center left", bbox_to_anchor=(0.02, 0.36),
               framealpha=0.95, fontsize=8.2)
    ax1.text(0.98, 0.06, r"$\mathrm{TM} = \frac{1}{L}\sum_{i=1}^L \frac{1}{1 + (d_i/d_0)^2}$",
             transform=ax1.transAxes, ha="right", va="bottom", fontsize=9.2, color=C_SLATE,
             bbox=dict(boxstyle="round,pad=0.26", fc="#F8FAFC", ec="#CBD5E1"))

    fig.subplots_adjust(left=0.06, right=0.93, top=0.89, bottom=0.16)
    fig.savefig(OUT_DIR / "correspondence_and_tm_score.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: Contact Map Signatures & Classical MDS Eigenvalues / Chirality
# ---------------------------------------------------------------------------

def fig4_distance_geometry_and_mds():
    helix = make_alpha_helix(24)
    par_b = make_beta_sheet(12, parallel=True)
    anti_b = make_beta_sheet(12, parallel=False)

    D_h = cdist(helix, helix)
    D_p = cdist(par_b, par_b)
    D_a = cdist(anti_b, anti_b)

    # Classical MDS on 3D helix-turn fold
    fold = helix.copy()
    fold[14:, 0] += np.linspace(0, 6.0, 10)
    fold[14:, 1] += np.linspace(0, 4.0, 10)
    D_fold = cdist(fold, fold)
    N = len(fold)
    H = np.eye(N) - np.ones((N, N)) / N
    G = -0.5 * (H @ (D_fold ** 2) @ H)
    eigvals = np.sort(np.linalg.eigvalsh(G))[::-1]
    eigvals_clean = np.maximum(eigvals[:10], 0.0)

    # Verify mirror invariance
    fold_mirror = fold.copy()
    fold_mirror[:, 2] *= -1.0
    max_mir_diff = float(np.max(np.abs(D_fold - cdist(fold_mirror, fold_mirror))))

    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    outer = gridspec.GridSpec(1, 2, width_ratios=[1.48, 1.0], wspace=0.24)
    left_gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0], wspace=0.14)

    maps = [
        (D_h, r"$\alpha$-Helix" + "\n(thick i→i+4 band)"),
        (D_p, r"Parallel $\beta$-Sheet" + "\n(j − i = const)"),
        (D_a, r"Anti-Parallel $\beta$" + "\n(i + j = const)"),
    ]
    im = None
    for idx, (D, title) in enumerate(maps):
        ax = fig.add_subplot(left_gs[idx])
        im = ax.imshow(D, origin="lower", cmap="viridis", vmin=0, vmax=26)
        ax.set_title(title, fontsize=9.8, pad=6)
        ax.set_xlabel("Residue j", fontsize=8.8)
        if idx == 0:
            ax.set_ylabel("Residue i", fontsize=8.8)
        else:
            ax.set_yticklabels([])
        ax.set_xticks([0, 11, 23])
        ax.set_xticklabels(["1", "12", "24"], fontsize=8.2)
        ax.set_yticks([0, 11, 23])
        if idx == 0:
            ax.set_yticklabels(["1", "12", "24"], fontsize=8.2)

    cbar_ax = fig.add_axes([0.055, 0.09, 0.495, 0.032])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label(r"$C\alpha$–$C\alpha$ Euclidean Distance $D_{ij}$ (Å)  —  100% $\mathrm{SE}(3)$-Invariant", fontsize=8.6)
    cbar.ax.tick_params(labelsize=8.0)

    # Right panel: Classical MDS Gram Spectrum & Mirror Ambiguity
    ax_r = fig.add_subplot(outer[1])
    k_idx = np.arange(1, 11)
    bar_cols = [C_BLUE, C_BLUE, C_BLUE] + ["#CBD5E1"] * 7
    bars = ax_r.bar(k_idx, eigvals_clean, color=bar_cols, width=0.62, edgecolor=C_SLATE, lw=0.8)
    for idx in range(3):
        ax_r.text(k_idx[idx], eigvals_clean[idx] + eigvals_clean[0] * 0.025,
                  f"{eigvals_clean[idx]:.0f}", ha="center", va="bottom",
                  fontsize=8.3, fontweight="bold", color=C_BLUE)
    ax_r.text(6.8, eigvals_clean[0] * 0.06, r"$\lambda_4 = \dots = \lambda_N = 0$" + "\n(exact rank 3 for 3D coords!)",
              ha="center", va="bottom", fontsize=8.3, color=C_MUTED, fontweight="bold")

    ax_r.text(0.97, 0.94,
              r"Gram Matrix: $G = -\frac{1}{2} H D^{\circ 2} H = X_c X_c^T$" + "\n"
              + r"Reconstruct: $X_c = V_3 \Lambda_3^{1/2}$" + "\n\n"
              + f"Mirror Ambiguity (E(3) vs SE(3)):\n"
              + r"$D(X) \equiv D(-X)$" + f"  (max |ΔD| = {max_mir_diff:.1e} Å)\n"
              + "Distances alone lose L vs D chirality!",
              transform=ax_r.transAxes, ha="right", va="top", fontsize=8.2, color=C_SLATE,
              bbox=dict(boxstyle="round,pad=0.35", fc="#FFF7ED", ec="#FDBA74"))

    ax_r.set_xticks(k_idx)
    ax_r.set_xticklabels([f"λ{k}" for k in k_idx], fontsize=8.5)
    ax_r.set_ylabel("Gram Eigenvalue (Å²)")
    ax_r.set_xlabel("MDS Component Index")
    ax_r.set_ylim(0, eigvals_clean[0] * 1.18)
    ax_r.set_title("Classical MDS: Rank-3 Recovery & Mirror Limit", fontsize=10.8, pad=8)
    ax_r.grid(True, axis="y", alpha=0.25)

    fig.subplots_adjust(left=0.05, right=0.975, top=0.87, bottom=0.24)
    fig.savefig(OUT_DIR / "distance_geometry_and_mds.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5: Debye-Waller Diffraction Attenuation & Electron-Density Smearing
# ---------------------------------------------------------------------------

def fig5_debye_waller_and_diffraction():
    fig = plt.figure(figsize=(12.8, 4.6), dpi=200)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.05, 1.0], wspace=0.26)

    # Left panel: Bragg intensity decay vs resolution d (6 Å -> 1.2 Å)
    ax0 = fig.add_subplot(gs[0])
    d = np.linspace(1.2, 6.0, 300)
    s2 = 1.0 / (4.0 * d ** 2)  # (sin theta / lambda)^2

    b_cases = [
        (15,  C_BLUE,   "-",  r"Rigid core ($B=15\ \mathrm{Å}^2,\ \sigma_u=0.44\ \mathrm{Å}$)"),
        (40,  C_TEAL,   "-",  r"Domain avg ($B=40\ \mathrm{Å}^2,\ \sigma_u=0.71\ \mathrm{Å}$)"),
        (79,  C_ORANGE, "--", r"Mobile loop ($B=79\ \mathrm{Å}^2,\ \sigma_u=1.00\ \mathrm{Å}$)"),
        (140, C_RED,    ":",  r"Disordered ($B=140\ \mathrm{Å}^2,\ \sigma_u=1.33\ \mathrm{Å}$)"),
    ]
    for B, col, ls, lbl in b_cases:
        atten = np.exp(-2.0 * B * s2)
        ax0.plot(d, atten, color=col, ls=ls, lw=2.3, label=lbl)

    # Mark 6LU7 X-ray resolution = 2.16 Å
    ax0.axvline(2.16, color=C_PURPLE, ls="-.", lw=1.4)
    ax0.text(2.26, 0.45, "6LU7 X-ray\nres = 2.16 Å", ha="right", va="center",
             fontsize=8.2, fontweight="bold", color=C_PURPLE,
             bbox=dict(boxstyle="round,pad=0.22", fc="#F5F3FF", ec="#DDD6FE"))

    ax0.set_xlim(6.0, 1.2)  # Inverted: coarse -> atomic resolution
    ax0.set_ylim(-0.02, 1.06)
    ax0.set_xlabel(r"Diffraction Resolution $d = \lambda / (2\sin\theta)$ (Å)  →  Finer Detail")
    ax0.set_ylabel(r"Relative Bragg Intensity  $I/I_0 = \exp(-2Bs^2)$")
    ax0.set_title("Debye–Waller Attenuation Kills High-Resolution Reflections", pad=8)
    ax0.legend(loc="upper right", framealpha=0.95, fontsize=7.9)
    ax0.grid(True, alpha=0.25)

    # Right panel: Real-Space Gaussian Electron-Density Smearing
    ax1 = fig.add_subplot(gs[1])
    x = np.linspace(-3.6, 3.6, 400)
    sigma0 = 0.35  # Intrinsic atomic form-factor width (Å)
    peak_ref = 1.0 / np.sqrt(2.0 * np.pi * (sigma0 ** 2 + 15.0 / (8.0 * np.pi ** 2)))

    for B, col, ls, _ in b_cases:
        sig_u = np.sqrt(B / (8.0 * np.pi ** 2))
        sig_tot = np.sqrt(sigma0 ** 2 + sig_u ** 2)
        rho = (1.0 / np.sqrt(2.0 * np.pi * sig_tot ** 2)) * np.exp(-0.5 * (x / sig_tot) ** 2)
        rho_rel = rho / peak_ref
        ax1.plot(x, rho_rel, color=col, ls=ls, lw=2.3, label=f"B = {B} Å² (σ_u = {sig_u:.2f} Å)")
        if B == 15:
            ax1.fill_between(x, 0, rho_rel, color=col, alpha=0.08)

    # 1-sigma crystallographic noise floor
    noise_cut = 0.45
    ax1.axhspan(0, noise_cut, color="#94A3B8", alpha=0.18, zorder=0)
    ax1.axhline(noise_cut, color=C_SLATE, ls="--", lw=1.2)
    ax1.text(3.45, noise_cut + 0.025, "1σ crystal noise floor\n(density vanishes below!)",
             ha="right", va="bottom", fontsize=7.9, fontweight="bold", color=C_SLATE)

    ax1.annotate("High-B loop smears flat\n→ 'missing residues' in PDB",
                 xy=(0.0, 0.405), xytext=(1.05, 0.78),
                 fontsize=8.1, fontweight="bold", color=C_RED,
                 arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.3),
                 bbox=dict(boxstyle="round,pad=0.24", fc="#FEF2F2", ec="#FECACA"))

    ax1.set_xlim(-3.5, 3.5)
    ax1.set_ylim(0, 1.12)
    ax1.set_xlabel("Displacement from Equilibrium Position x (Å)")
    ax1.set_ylabel("Unit-Cell Averaged Density (rel. peak)")
    ax1.set_title(r"Real-Space Electron Density Smearing ($B = 8\pi^2\langle u^2\rangle$)", pad=8)
    ax1.legend(loc="upper left", framealpha=0.94, fontsize=8.0)
    ax1.grid(True, alpha=0.25)

    fig.subplots_adjust(left=0.06, right=0.985, top=0.89, bottom=0.15)
    fig.savefig(OUT_DIR / "debye_waller_and_diffraction.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    print("Generating Lecture 7 Part 3 figures...")
    fig1_dihedral_ramachandran()
    fig2_kabsch_se3_and_chirality()
    fig3_correspondence_and_tm_score()
    fig4_distance_geometry_and_mds()
    fig5_debye_waller_and_diffraction()
    print(f"All 5 figures generated cleanly in {OUT_DIR}")
