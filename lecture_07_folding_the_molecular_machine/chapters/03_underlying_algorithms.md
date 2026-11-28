:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Dihedral Geometry · Kabsch SVD · TM-Score · Contact Maps · MDS · Debye–Waller
:::

---

# Curriculum & Algorithmic Roadmap

Six geometric, linear-algebraic, and physical algorithms transform raw atomic coordinates (`ATOM` records $\in \mathbb{R}^{N \times 3}$) into validated folds, rigid superpositions, invariant contact maps, and crystal mobility profiles in Part 4:

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "pdb",  "label": "3D Coordinates\n(ATOM ∈ R^{N×3})",       "x": 20,   "y": 40, "width": 150, "height": 64 },
    { "id": "dih",  "label": "Internal Dihedrals\n(φ, ψ, ω via atan2)", "x": 205,  "y": 40, "width": 170, "height": 64, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "kab",  "label": "Kabsch SVD\n(SE(3), det R = +1)",         "x": 410,  "y": 40, "width": 165, "height": 64, "style": "fill: #2563EB; color: #fff;" },
    { "id": "tm",   "label": "Correspondence +\nTM-Score & GDT-TS",     "x": 610,  "y": 40, "width": 165, "height": 64, "style": "fill: #16A34A; color: #fff;" },
    { "id": "map",  "label": "Contact Maps D_ij\n& Classical MDS",      "x": 810,  "y": 40, "width": 165, "height": 64, "style": "fill: #D97706; color: #fff;" },
    { "id": "dw",   "label": "Debye–Waller\nB = 8π²⟨u²⟩ Physics",       "x": 1010, "y": 40, "width": 165, "height": 64, "style": "fill: #DC2626; color: #fff;" }
  ],
  "links": [
    { "from": "pdb", "to": "dih" },
    { "from": "dih", "to": "kab" },
    { "from": "kab", "to": "tm"  },
    { "from": "tm",  "to": "map" },
    { "from": "map", "to": "dw"  }
  ]
}
:::

:::matrix { .cards cols="1fr 1fr 1fr" gap="18px" }
[[0,0]]{ .purple }
### 1. Internal Backbone Geometry
Peptide planar resonance ($\omega \approx 180^\circ$), signed 4-atom torsion angles via `atan2`, and hard-sphere Ramachandran steric basins (L-Ala vs. Gly vs. Pro).

[[0,1]]{ .blue }
### 2. $\text{SE}(3)$ Kabsch SVD & Scoring
Centroid decoupling, $3\times 3$ cross-covariance SVD, right-handed $\det(R) = +1$ chirality guard, row/column conventions, frame-shift traps, and TM-score / GDT-TS.

[[0,2]]{ .orange }
### 3. Distance Geometry & Crystal Physics
$\text{SE}(3)$-invariant $C\alpha$ contact maps, Classical MDS rank-3 reconstruction & mirror ambiguity, and Debye–Waller diffraction attenuation of mobile loops.
:::

---

# Backbone Internal Coordinates: Peptide Planarity & $(\phi, \psi, \omega)$

A 306-residue polypeptide (`6LU7`) contains $\sim 2{,}400$ heavy atoms ($7{,}200$ Cartesian coordinates), yet bond lengths, bond angles, and amide resonance collapse its backbone conformation onto just **two torsion angles per residue**:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Rigid Planar Peptide Units ($\omega$)
- Delocalization of the amide nitrogen lone pair across $\mathrm{O=C-N}$ gives the peptide bond $\sim 40\\%$ partial double-bond character ($1.33\text{ Å}$ vs. $1.46\text{ Å}$ single bond).
- All 6 atoms $(C\alpha_i, C_i, O_i, N_{i+1}, H_{i+1}, C\alpha_{i+1})$ lie in a flat plane, locking the peptide torsion angle:
$$
\omega_i = \angle(C\alpha_i, C_i, N_{i+1}, C\alpha_{i+1}) \approx 180^\circ \text{ (trans)}
$$
- **Rare *cis* exception ($\omega \approx 0^\circ$):** `< 0.05%` of non-Proline bonds, rising to `~5%` before **Proline** where *cis* and *trans* steric clashes nearly balance.

[[0,1]]
### Two Free Backbone Rotations per Residue
Only the two single bonds flanking each tetrahedral $C\alpha_i$ hinge rotate freely:
- **$\phi_i$ (N–$C\alpha$ bond):** defined by four atoms $(C_{i-1}, N_i, C\alpha_i, C_i)$.
- **$\psi_i$ ($C\alpha$–C bond):** defined by four atoms $(N_i, C\alpha_i, C_i, N_{i+1})$.
- **Boundary condition (Termini & Chain Breaks):**
  - Residue $1$ has no $C_{i-1} \implies \phi_1$ is undefined (`None`).
  - Residue $N$ has no $N_{i+1} \implies \psi_N$ is undefined (`None`).
  - Missing crystal loop gaps reset both angles across the chain break!
:::

---

# Computing Signed Dihedral Angles via Plane Normals & `atan2`

Given four sequential atomic positions $\mathbf{r}\_1, \mathbf{r}\_2, \mathbf{r}\_3, \mathbf{r}\_4 \in \mathbb{R}^3$, how do we compute the exact signed IUPAC dihedral angle $\theta \in (-\pi, \pi]$?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Bond Vectors & Half-Plane Normals
1. Form the three consecutive directed bond vectors:
$$
\mathbf{b}\_1 = \mathbf{r}\_2 - \mathbf{r}\_1, \quad \mathbf{b}\_2 = \mathbf{r}\_3 - \mathbf{r}\_2, \quad \mathbf{b}\_3 = \mathbf{r}\_4 - \mathbf{r}\_3
$$
2. Construct the two half-plane normals perpendicular to the central rotatable bond $\mathbf{b}_2$ and unit axis $\hat{\mathbf{u}} = \mathbf{b}\_2 / \lVert \mathbf{b}\_2 \rVert_2$:
$$
\mathbf{n}\_1 = \mathbf{b}\_1 \times \mathbf{b}\_2, \qquad \mathbf{n}\_2 = \mathbf{b}\_2 \times \mathbf{b}\_3
$$
- Notice $\mathbf{n}_1 \perp \hat{\mathbf{u}}$ and $\mathbf{n}_2 \perp \hat{\mathbf{u}}$, so both normals lie strictly in the plane orthogonal to the rotation axis!

[[0,1]]
### Why `arccos` Fails & the Exact `atan2` Formula
- Using $\arccos\bigl((\mathbf{n}\_1 \cdot \mathbf{n}\_2)/(\lVert\mathbf{n}\_1\rVert \lVert\mathbf{n}\_2\rVert)\bigr)$ returns only $[0, \pi]$ — collapsing right-handed ($-\theta$) and left-handed ($+\theta$) twists and losing numerical accuracy near $0^\circ$ and $180^\circ$.
- Because $\mathbf{n}\_1 \cdot \mathbf{n}\_2 = \lVert\mathbf{n}\_1\rVert\lVert\mathbf{n}\_2\rVert \cos\theta$ and $(\mathbf{n}\_1 \times \mathbf{n}\_2) \cdot \hat{\mathbf{u}} = \lVert\mathbf{n}\_1\rVert\lVert\mathbf{n}\_2\rVert \sin\theta$, the common norm cancels identically:
$$
\theta = \operatorname{atan2}\Bigl((\mathbf{n}\_1 \times \mathbf{n}\_2) \cdot \hat{\mathbf{u}}, \ \mathbf{n}\_1 \cdot \mathbf{n}\_2\Bigr) \in (-\pi, \pi]
$$
- Sign convention: looking down $\mathbf{r}\_2 \to \mathbf{r}\_3$, clockwise rotation of the front bond $\mathbf{b}_1$ onto the rear bond $\mathbf{b}_3$ is **positive**!
:::

---

# Ramachandran Steric Physics & Amino-Acid Exceptions

![Peptide Dihedral Geometry and Analytical Hard-Sphere Ramachandran Map](lecture_07_folding_the_molecular_machine/images/dihedral_and_ramachandran_sterics.png){ width=94% center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Steric hard-sphere origin (Ramachandran et al., 1963):** Steep Pauli-repulsion walls $V_{\text{rep}}(r) \propto (r_{\min}/r)^{12}$ forbid `> 75%` of the $(\phi, \psi)$ torus where non-bonded atoms penetrate within van der Waals contact distance — carving out the canonical **$\beta$-sheet** (`−135°, +140°`) and **right-handed $\alpha$-helix** (`−60°, −45°`) valleys.

[[0,1]]
- **Chirality breaks symmetry:** Every L-amino acid places its $C\beta$ side chain asymmetrically on $C\alpha$, causing violent $C\beta \cdots O_{i-1}$ clashes across the $+\phi$ half.
- **Two structural exceptions:** Achiral **Glycine** (no $C\beta$) freely populates the **left-handed $\alpha$** basin (`+60°, +40°`), whereas **Proline's** cyclic pyrrolidine ring locks $\phi \approx -63^\circ$!
:::

---

# Interactive Lab — Signed Dihedrals & Hard-Sphere Ramachandran Explorer

Explore backbone internal coordinates live: **(1)** in **Tab 1**, drag the 3D peptide planes, switch presets ($\alpha$-helix, $\beta$-sheet, *cis*-$\omega$), and verify that `atan2` preserves the exact minus sign that `arccos` drops whenever $\theta < 0$; **(2)** in **Tab 2**, click or drag across the analytical **Hard-Sphere Ramachandran Map** while switching between **L-Alanine**, **Glycine**, and **Proline** to watch red van der Waals clashes flash in 3D!

<iframe src="lecture_07_folding_the_molecular_machine/demo/dihedral_ramachandran_lab.html" title="Interactive Signed Dihedrals and Hard-Sphere Ramachandran Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---
# Rigid-Body Motion ($\text{SE}(3)$) & the Orthogonal Procrustes Problem

Suppose we solve two crystal structures of the same protease (`6LU7` vs. `1UK4`) with paired $C\alpha$ matrices $P, Q \in \mathbb{R}^{N \times 3}$. Why does raw subtraction $\lVert P - Q \rVert_F$ fail catastrophically?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Arbitrary Crystallographic Reference Frames
- Each X-ray crystal lattice fixes its own unit-cell origin and orientation, so identical folds differ by an arbitrary **rigid-body motion** $(R, \mathbf{t}) \in \text{SE}(3)$:
$$
\text{SE}(3) = \bigl\\{(R, \mathbf{t}) : R \in \text{SO}(3), \mathbf{t} \in \mathbb{R}^3\bigr\\}
$$
- Here **$\text{SO}(3)$** is the Special Orthogonal group of **proper rotations**:
$$
R^T R = I_3 \quad \text{and} \quad \det(R) = +1
$$
- Without removing $(R, \mathbf{t})$, raw point-to-point RMSD exceeds **`100 Å`** even for identical molecules!

[[0,1]]
### Optimal Superposition Objective (RMSD)
We seek the rigid motion $(R^\*, \mathbf{t}^\*) \in \text{SE}(3)$ minimizing mean squared atomic deviation:
$$
\text{RMSD}(P, Q) = \min_{R \in \text{SO}(3), \mathbf{t} \in \mathbb{R}^3} \sqrt{\frac{1}{N} \sum_{i=1}^N \lVert R\mathbf{p}\_i + \mathbf{t} - \mathbf{q}\_i \rVert_2^2}
$$
- **Why $\det(R) = +1$ is physically sacred:** Orthogonal matrices with $\det(R) = -1$ are **mirror reflections** — physically impossible for rigid proteins because they invert L-amino acids into D-enantiomers!
- Dividing by **$N$** (atoms), **not $3N$** (scalars), keeps RMSD in true 3D Euclidean Angstroms.
:::

---

# The Kabsch Algorithm: Centroid Decoupling & Trace Maximization

Wolfgang Kabsch (1976, 1978) derived the exact closed-form global optimum over $\text{SE}(3)$ in four algebraic stages:

:::static-diagram { width=95% center }
{
  "nodes": [
    { "id": "pq",   "label": "Paired Cα Arrays\nP, Q ∈ R^{N×3}",             "x": 20,  "y": 32, "width": 160, "height": 62 },
    { "id": "cent", "label": "1. Center Centroids\nP_c, Q_c (kills t*)",      "x": 220, "y": 32, "width": 175, "height": 62, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "cov",  "label": "2. Cross-Covariance\nC = P_c^T Q_c ∈ R^{3×3}",  "x": 435, "y": 32, "width": 185, "height": 62, "style": "fill: #2563EB; color: #fff;" },
    { "id": "svd",  "label": "3. 3×3 SVD\nC = U S V^T",                       "x": 660, "y": 32, "width": 150, "height": 62, "style": "fill: #16A34A; color: #fff;" },
    { "id": "rot",  "label": "4. Chirality Guard\nd = sign(det(V U^T))",      "x": 850, "y": 32, "width": 180, "height": 62, "style": "fill: #D97706; color: #fff;" },
    { "id": "rmsd", "label": "Optimal SO(3)\nRMSD (divide by N)",             "x": 1070,"y": 32, "width": 165, "height": 62, "style": "fill: #DC2626; color: #fff;" }
  ],
  "links": [
    { "from": "pq",   "to": "cent" },
    { "from": "cent", "to": "cov"  },
    { "from": "cov",  "to": "svd"  },
    { "from": "svd",  "to": "rot"  },
    { "from": "rot",  "to": "rmsd" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Step 1: Decoupling Translation $\mathbf{t}$
Differentiate $E(R, \mathbf{t}) = \sum_{i=1}^N \lVert R\mathbf{p}\_i + \mathbf{t} - \mathbf{q}\_i \rVert_2^2$ with respect to $\mathbf{t}$ and set $\nabla_{\mathbf{t}} E = \mathbf{0}$:
$$
\mathbf{t}^* = \mathbf{c}\_Q - R\mathbf{c}\_P, \qquad \mathbf{c}\_P = \frac{1}{N}\sum_{i=1}^N \mathbf{p}\_i, \quad \mathbf{c}\_Q = \frac{1}{N}\sum_{i=1}^N \mathbf{q}\_i
$$
Substituting $\mathbf{t}^*$ centers both point clouds at the origin ($P_c = P - \mathbf{1}\mathbf{c}\_P^T$, $Q_c = Q - \mathbf{1}\mathbf{c}\_Q^T$), removing $\mathbf{t}$ completely!

[[0,1]]
### Step 2: Frobenius Expansion $\to$ Trace Maximization
For row-stored matrices $P_c, Q_c \in \mathbb{R}^{N\times 3}$, rotating column vectors $\mathbf{p}\_i \mapsto R\mathbf{p}\_i$ acts on rows as $P_c R^T$. Since $R^T R = I_3$:
$$
\lVert P_c R^T - Q_c \rVert_F^2 = \lVert P_c \rVert_F^2 + \lVert Q_c \rVert_F^2 - 2\operatorname{tr}\bigl(R \underbrace{P_c^T Q_c}\_{C \in \mathbb{R}^{3\times 3}}\bigr)
$$
Because $\lVert P_c \rVert_F^2 + \lVert Q_c \rVert_F^2$ is constant in $R$, minimizing RMSD is **exact equivalence** to maximizing $\operatorname{tr}(R C)$ over $R \in \text{SO}(3)$!
:::

---

# Kabsch SVD Proof, Chirality Guard & Matrix Conventions

Now take the Singular Value Decomposition of the tiny $3 \times 3$ cross-covariance matrix $C = P_c^T Q_c = U S V^T$ with $S = \operatorname{diag}(s_1 \ge s_2 \ge s_3 \ge 0)$:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Exact Trace Bound & Reflection Check
Using cyclic invariance of the trace ($\operatorname{tr}(AB) = \operatorname{tr}(BA)$):
$$
\operatorname{tr}(R C) = \operatorname{tr}(R U S V^T) = \operatorname{tr}\bigl(\underbrace{(V^T R U)}\_{W \in \text{O}(3)} S\bigr) = \sum_{k=1}^3 W_{kk} s_k
$$
- Since every column of orthogonal $W$ has unit norm, $\lvert W_{kk} \rvert \le 1$, so $\operatorname{tr}(W S) \le s_1 + s_2 + s_3$ — achieved at $W = I_3 \implies R = V U^T$.
- **Chirality guard:** $\det(W) = \det(V^T)\det(R)\det(U) = \det(R)\det(V U^T)$. If $\det(V U^T) = -1$, requiring $\det(R) = +1$ forces $\det(W) = -1$, whose exact optimum flips the smallest singular value $s_3$:
$$
W^* = \operatorname{diag}(1, 1, d), \qquad d = \operatorname{sign}\bigl(\det(V U^T)\bigr) \in \{+1, -1\}
$$

[[0,1]]
### Two Equivalent Conventions & Two Classic Bugs
1. **Column-vector convention** (Practice 7 Task 1):
$$
R_{\text{col}} = V \operatorname{diag}(1, 1, d) U^T \implies P_{\text{rot}} = (R_{\text{col}} P_c^T)^T = P_c R_{\text{col}}^T
$$
2. **Row-vector convention** (direct right-multiply):
$$
R_{\text{row}} = U \operatorname{diag}(1, 1, d) V^T = R_{\text{col}}^T \implies P_{\text{rot}} = P_c R_{\text{row}}
$$
- **Bug 1 (Transpose trap):** Computing `Pc @ R_col` instead of `Pc @ R_col.T` rotates by $R^{-1}$ (backwards!).
- **Bug 2 ($3N$ trap):** `np.sqrt(np.mean((P_rot - Qc)**2))` divides by $3N$ instead of $N$, underestimating RMSD by $\sqrt{3} \approx 1.732\times$!
:::

---
# Kabsch SVD vs. PCA SVD --- Same Tool, Different Geometry

![Kabsch Three-Stage Progression and Transpose vs Chirality Traps](lecture_07_folding_the_molecular_machine/images/kabsch_se3_and_chirality.png){ width=89 center }

:::matrix { cols="60/40" gap="24px" }
[[0,0]]
| Property | **PCA (Lecture 6)** | **Kabsch (Lecture 7)** |
| :--- | :--- | :--- |
| **Input Matrix** | Single matrix $X \in \mathbb{R}^{n \times p}$ | Paired frames $P_c, Q_c \in \mathbb{R}^{N \times 3}$ |
| **Matrix Decomposed** | $X^\top X \in \mathbb{R}^{p\times p}$ (auto-cov.) | $C = P_c^T Q_c \in \mathbb{R}^{3\times 3}$ (cross-cov.) |
| **Geometric Goal** | Maximize variance of 1 cloud | Align 2 frames via $R \in \text{SO}(3)$ |
| **Sign / Reflection** | Arbitrary flips $(-v_i)$ valid | $\det(R) = +1$ **mandatory** (L-aa!) |

[[0,1]]
- **Why the distinction matters:** PCA decomposes *one* dataset's self-variance where reflecting an axis leaves every pairwise distance intact. Kabsch decomposes the *cross-correlation* between two 3D frames — skipping $d = \operatorname{sign}(\det(VU^T))$ illegally mirrors an L-helix onto its D-enantiomer (faking `0.00 Å` RMSD!), while dropping `.T` on $R_{\text{col}}$ spins the chain backwards!
:::

---
# Correspondence & Length-Normalized Scores (TM, GDT-TS)

![Residue Frame-Shift Pitfall and RMSD vs TM-Score and GDT-TS Stability](lecture_07_folding_the_molecular_machine/images/correspondence_and_tm_score.png){ width=73 center }

:::matrix { cols="42/58" gap="5px" }
[[0,0]]
- **Correspondence trap (naive `[:min_len]`):** Kabsch requires exact 1-to-1 residue pairing. Because `6LU7` spans residues `1–306` while `1UK4` starts at residue `3` (`3–303`), naive slicing shifts every residue by 2 positions (`6.21 Å` artifact vs. true `1.35 Å` when paired by `res.id[1]`)!
- **RMSD's quadratic weakness:** Because $\text{RMSD} \propto \sqrt{\sum d_i^2}$, a single flexible terminal tail swinging `16 Å` inflates global RMSD dramatically even when 90% of the core superimposes cleanly (`< 0.8 Å`).

[[0,1]]
- **TM-score (Zhang & Skolnick, 2004):** Replaces quadratic error with a bounded, length-normalized Lorentzian kernel:
$$
\text{TM} = \max_{R, \mathbf{t}} \frac{1}{L_{\text{target}}} \sum_{i=1}^{L_{\text{ali}}} \frac{1}{1 + (d_i / d_0)^2}, \\; d_0 = 1.24\sqrt[3]{L_{\text{target}} - 15} - 1.8
$$
  - Random unrelated folds score `≈ 0.17`; **`TM > 0.50`** indicates the same structural fold across all protein lengths!
- **GDT-TS (CASP gold standard):** Average fraction of aligned $C\alpha$ atoms within `1, 2, 4, and 8 Å` cutoffs.
:::

---

# Interactive Lab — 3D Kabsch SVD, Chirality Trap & TM-Score

Test rigid-body superposition live: **(1)** in **Tab 1**, step from **Raw Frames** $\to$ **Centroid Centered** $\to$ **Kabsch SVD**, inspect the exact $3\times 3$ cross-covariance $C$ and singular values, then activate **Mirror D-Peptide** and uncheck **Enforce `det(R) = +1`** (or toggle **Transpose Bug `Pc @ R`**) to trigger both traps; **(2)** in **Tab 2**, introduce an **Off-by-2 Frame Shift** or swing the **Flexible Tail** to compare RMSD vs. TM-score and GDT-TS!

<iframe src="lecture_07_folding_the_molecular_machine/demo/kabsch_superposition_lab.html" title="Interactive 3D Kabsch SVD Superposition, Chirality Trap and TM-Score Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# $\text{SE}(3)$-Invariant Distance Geometry & Contact Maps

Can we compare protein folds and detect secondary structure *without* solving any alignment or choosing any global coordinate frame? Yes — move to **pairwise internal distances**!

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Exact Proof of Rigid-Body Invariance
Define the symmetric $N \times N$ $C\alpha$ Euclidean distance matrix $D_{ij} = \lVert \mathbf{x}\_i - \mathbf{x}\_j \rVert_2$. Under any rigid motion $\mathbf{x}\_i' = R\mathbf{x}\_i + \mathbf{t}$ with $R \in \text{SO}(3)$:
$$
\begin{aligned}
D_{ij}' &= \lVert (R\mathbf{x}\_i + \mathbf{t}) - (R\mathbf{x}\_j + \mathbf{t}) \rVert_2 \\
&= \lVert R(\mathbf{x}\_i - \mathbf{x}\_j) \rVert_2 = \lVert \mathbf{x}\_i - \mathbf{x}\_j \rVert_2 = D_{ij}
\end{aligned}
$$
- Translation $\mathbf{t}$ cancels identically; orthogonal $R$ preserves inner products ($(R\mathbf{v})^T(R\mathbf{v}) = \mathbf{v}^T R^T R \mathbf{v} = \lVert\mathbf{v}\rVert_2^2$).
- **Binary contact map:** $A_{ij} = \mathbb{I}(D_{ij} < r_c)$ (standard $C\alpha$ cutoff $r_c = 8.0\text{ Å}$).

[[0,1]]
### Analytical Motifs in the Distance Matrix
1. **Consecutive Backbone ($|i - j| = 1$):** Rigid *trans* peptide planes fix $D_{i, i+1} \approx 3.8\text{ Å}$ everywhere along the main diagonal.
2. **$\alpha$-Helix:** Regular $3.6$ residues/turn with $i \to i+4$ H-bonds keeps $D_{i, i+3} \approx 5.1\text{ Å}$ and $D_{i, i+4} \approx 5.5\text{ Å}$ $\implies$ **thick dark band** flanking the main diagonal.
3. **Parallel $\beta$-Sheet:** Adjacent strands run $N \to C$ together: $(i, j), (i+1, j+1), \dots \implies j - i = \text{const}$ (**stripe parallel** to diagonal, $\sim 4.8\text{ Å}$).
4. **Anti-Parallel $\beta$-Sheet:** Strands run opposite ($N \to C$ vs. $C \to N$): $(i, j), (i+1, j-1), \dots \implies i + j = \text{const}$ (**stripe perpendicular** to diagonal!).
:::

---
# From Distances Back to 3D: Classical MDS & Mirror Ambiguity

![Contact Map Secondary Motifs and Classical MDS Rank-3 Spectrum](lecture_07_folding_the_molecular_machine/images/distance_geometry_and_mds.png){ width=89 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Classical Multidimensional Scaling (MDS):** Expanding $D_{ij}^2 = \lVert\mathbf{x}\_i\rVert^2 + \lVert\mathbf{x}\_j\rVert^2 - 2\mathbf{x}\_i^T\mathbf{x}\_j$ and double-centering with $H = I_N - \frac{1}{N}\mathbf{1}\mathbf{1}^T$ recovers the exact centered Gram matrix:
$$
G = -\frac{1}{2} H D^{\circ 2} H = X_c X_c^T \quad \Longrightarrow \quad X_c = V_3 \Lambda_3^{1/2}
$$
  Exact 3D coordinates produce **precisely 3 positive eigenvalues** $(\lambda_1, \lambda_2, \lambda_3 > 0)$ and $\lambda_4 = \dots = \lambda_N = 0$!

[[0,1]]
- **The Mirror Ambiguity ($\text{E}(3)$ vs. $\text{SE}(3)$):** Because $D_{ij}$ is invariant under *all* isometries $\text{O}(3)$ — including reflections $\mathbf{x}\_i \mapsto -\mathbf{x}\_i$ — a right-handed helix and its left-handed mirror share the **exact same distance matrix** ($D(X) \equiv D(-X)$)!
- **Practical impact:** NMR NOE distance geometry and AlphaFold distograms *must* inject signed stereochemical/dihedral constraints to pick the biological L-enantiomer!
:::

---

# Crystallographic Disorder: Bragg Interference & Debye–Waller Physics

In X-ray crystallography, diffraction does not image a single molecule — it measures **coherent wave interference** summed simultaneously across $\sim 10^{12}$ periodic unit cells in a crystal:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Structure Factors & Bragg Attenuation
- For reciprocal lattice reflection $\mathbf{h}$ at resolution $d = \lambda / (2\sin\theta)$ (so $s = \frac{\sin\theta}{\lambda} = \frac{1}{2d}$), the complex structure factor sums every atom $j$ in the unit cell:
$$
F(\mathbf{h}) = \sum_{j} f_j(s) \ e^{2\pi i \ \mathbf{h} \cdot \mathbf{r}\_j} \ \underbrace{\exp\bigl(-B_j s^2\bigr)}\_{\text{Debye–Waller factor}}
$$
- Measured Bragg reflection intensity scales as $\lvert F(\mathbf{h}) \rvert^2$, attenuating by $\exp(-2 B s^2) = \exp\bigl(-B / (2d^2)\bigr)$!
- Notice the **$1/d^2$ exponent**: thermal disorder destroys fine high-resolution reflections **exponentially faster** than coarse low-resolution reflections.

[[0,1]]
### The Debye–Waller Equation ($B$-Factor)
For isotropic harmonic displacement around equilibrium with 1D mean-square displacement $\langle u^2 \rangle$ (in $\text{Å}^2$):
$$
B = 8\pi^2 \langle u^2 \rangle \qquad \Longleftrightarrow \qquad \sigma_u = \sqrt{\langle u^2 \rangle} = \sqrt{\frac{B}{8\pi^2}}
$$
- **Benchmark calibration:**
  - Rigid core ($B = 15\text{ Å}^2$): $\sigma_u \approx 0.44\text{ Å}$ (crisp atomic peak).
  - Flexible loop ($B = 79\text{ Å}^2$): $\sigma_u = \sqrt{79 / 78.96} \approx 1.00\text{ Å}$!
- **Fourier duality:** Multiplying reciprocal space by a Gaussian $\exp(-B s^2)$ convolves real-space atomic electron density with a Gaussian of width $\sigma_u$, smearing mobile atoms flat!
:::

---
# Why Mobile Surface Loops Disappear in X-Ray Crystals (vs. Cryo-EM)

![Debye-Waller Diffraction Decay and Real-Space Electron Density Smearing](lecture_07_folding_the_molecular_machine/images/debye_waller_and_diffraction.png){ width=84 center }

:::matrix { cols="45/55" gap="5px" }
[[0,0]]
**Destructive phase interference in crystals:** 
<br>
Rigid hydrophobic cores and catalytic active sites occupy identical coordinates in every unit cell, so their X-ray waves add **constructively**. Mobile surface loops adopt distinct conformations across unit cells ($\sigma_u > 1\text{ Å}$, $B > 80\text{ Å}^2$): their scattered waves cancel **destructively**, sinking density below the $1\sigma$ noise floor — leaving **missing residues** in PDB files!

[[0,1]]
**X-ray vs. Cryo-EM contrast:**
- **X-ray crystallography:** Averages over an entire periodic lattice simultaneously; flexible loops either blur out or get pinned into non-native crystal-packing contacts.
- **Cryo-EM:** Flash-freezes individual molecules in vitreous ice and applies 2D/3D single-particle classification — separating discrete conformational states computationally instead of averaging them away!
:::

---
# Interactive Lab — Contact Maps ($\text{SE}(3)$ Invariance) & Debye–Waller Mobility

Verify distance invariance and diffraction physics live: **(1)** in **Tab 1**, switch among **$\alpha$-Helix**, **Parallel $\beta$**, **Anti-Parallel $\beta$**, and **$\beta\alpha\beta$**, click any pixel $(i, j)$ to inspect the 3D chord, and hit **Random $\text{SE}(3)$ Motion** or **Mirror Reflect ($z \to -z$)** to confirm `max |ΔD| = 0.00e+00`; **(2)** in **Tab 2**, sweep **$B$-Factor** (`5 → 160 Å²`) and **Resolution $d_{\min}$** (`2.16 Å` = `6LU7`) to watch high-resolution Bragg reflections and loop density drop into noise!

<iframe src="lecture_07_folding_the_molecular_machine/demo/contact_map_bfactor_lab.html" title="Interactive SE(3)-Invariant Contact Maps and Debye-Waller Mobility Lab" style="width: 100%; height: 580px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Structural Alignment Beyond Exact Pairing (TM-Align, DALI & Foldseek)

Kabsch SVD solves alignment in $O(N)$ **only when residue correspondence is already known**. What when two remote homologs share `< 15%` sequence identity and differ by insertions/deletions?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Iterative Coordinate & Distance-Matrix Aligners
- **TM-align (Zhang & Skolnick, 2005):** Alternates two steps until convergence:
  1. Run Needleman–Wunsch **Dynamic Programming** on the current residue-pair TM similarity matrix $S_{ij} = \bigl[1 + (d_{ij}/d_0)^2\bigr]^{-1}$ to update correspondence.
  2. Re-run **Kabsch SVD** (weighted toward core pairs) to refine $(R, \mathbf{t})$.
- **DALI (Holm & Sander, 1993):** Bypasses 3D frames altogether by aligning overlapping $6 \times 6$ sub-matrices of the **$\text{SE}(3)$-invariant distance matrices** $D^{(A)}$ and $D^{(B)}$ via Monte Carlo assembly!

[[0,1]]
### Database-Scale Structural Search: Foldseek (`3Di`)
- Aligning one query against the 200M-structure AlphaFold Database via TM-align/DALI would take months on a cluster.
- **Foldseek (van Kempen et al., 2024):** Encodes each residue's **nearest 3D tertiary neighbor geometry** (relative backbone orientation, virtual dihedrals, distance) through a VQ-VAE into a **20-letter structural alphabet (`3Di`)**.
- Once 3D geometry becomes a 1D `3Di` string, classic **MMseqs2 k-mer prefiltering + banded Smith–Waterman** achieves TM-align sensitivity at **$10{,}000\times$ higher speed** — unifying Lecture 2's DP with Lecture 7's 3D geometry!
:::

---
# Algorithmic Complexity & Scaling in 3D Structural Biology

Because Kabsch SVD acts on a **fixed $3 \times 3$ matrix** regardless of protein length $N$, structural superposition is dramatically faster than pairwise distance maps or MDS reconstruction:

:::matrix { cols=100 }
[[0,0]]{ .dense }
| Algorithm | Mathematical Core | Time Complexity | Space | Key Bottleneck / Practical Note |
| :--- | :--- | :---: | :---: | :--- |
| **Backbone Dihedrals** | 4-atom cross products + `atan2` | $O(N)$ | $O(N)$ | Terminal residues & crystal loop breaks yield `None` |
| **Kabsch Superposition** | Centroids + $3\times 3$ SVD ($C = P_c^T Q_c$) | $O(N)$ | $O(N)$ | $3\times 3$ SVD is $O(1)$; requires pre-matched pairs |
| **Contact Map (`cdist`)** | All-pairs $D_{ij} = \lVert\mathbf{x}\_i - \mathbf{x}\_j\rVert_2$ | $O(N^2)$ | $O(N^2)$ | $100\\%$ $\text{SE}(3)$-invariant; loses global chirality |
| **Classical MDS** | Eigendecomposition of $G = -\frac{1}{2}HD^{\circ 2}H$ | $O(N^3)$ (exact) | $O(N^2)$ | $O(N^2)$ with Lanczos top-3; mirror ambiguity $\pm X$ |
| **TM-Align (Unknown Map)** | Iterative DP + weighted Kabsch | $O(k \cdot N_1 N_2)$ | $O(N_1 N_2)$ | Converges in $k \approx 5\text{–}15$ DP–SVD sweeps |
| **Debye–Waller Profile** | $B_i = 8\pi^2 \langle u_i^2 \rangle$ extraction | $O(N)$ | $O(N)$ | Directly stored per atom in `.pdb` / `.mmcif` records |
:::

:::matrix { .cards cols="1fr 1fr" gap="18px" }
[[0,0]]{ .green }
### Why $3 \times 3$ SVD Is Practically Free
For a 306-residue protease (`6LU7`), forming $C = P_c^T Q_c$ costs just `~2,750` multiply-adds, and the $3\times 3$ SVD takes microseconds — enabling inner-loop Kabsch calls inside TM-align, MD trajectory clustering, and AlphaFold's FAPE loss!

[[0,1]]{ .purple }
### When $O(N^2)$ Representations Shine
Modern protein structure models (AlphaFold2, ESMFold) operate on $N \times N$ **pair representations** (generalized distance/orientation maps) precisely because they are natively $\text{SE}(3)$-invariant before generating final 3D backbone frames!
:::

---

# Summary & Bridge to Part 4 (Industrial Tools) & Lecture 8

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### Backbone & Sterics
#### Internal Coords
- Planar resonance locks $\omega \approx 180^\circ$
- Signed $\phi, \psi$ via `atan2` on $\mathbf{n}\_1, \mathbf{n}\_2$
- Hard-sphere $C\beta$ clashes carve basins
- Glycine unlocks $+\phi$; Proline pins $\phi$
- **Part 4:** `Bio.PDB.PPBuilder`

[[0,1]]{ .blue }
### Kabsch & Scoring
#### Rigid $\text{SE}(3)$
- Centering decouples $\mathbf{t}^* = \mathbf{c}\_Q - R\mathbf{c}\_P$
- $3\times 3$ SVD of $C = P_c^T Q_c$ maximizes $\operatorname{tr}(RC)$
- Guard $\det(R) = +1$; watch `Pc @ R.T` & $/N$
- Match residue IDs; TM-score resists tails
- **Part 4:** `Bio.PDB.Superimposer`

[[0,2]]{ .orange }
### Maps & Mobility
#### Invariants & Physics
- $D_{ij}$ is $\text{SE}(3)$-invariant ($\alpha$ vs. $\beta$ stripes)
- MDS recovers rank-3 $X_c$ up to mirror $\pm X$
- $B = 8\pi^2 \langle u^2 \rangle$ smears high-res Bragg spots
- High-$B$ loops vanish in X-ray crystals
- **Part 4:** `scipy.cdist`, `get_bfactor()`
:::

> **Core Equations Bridge:** $\theta = \operatorname{atan2}\bigl((\mathbf{n}\_1 \times \mathbf{n}\_2)\cdot\hat{\mathbf{u}}, \ \mathbf{n}\_1\cdot\mathbf{n}\_2\bigr)$ · $R_{\text{col}} = V\operatorname{diag}(1,1,d)U^T$ · $\text{TM} = \frac{1}{L}\sum \frac{1}{1+(d_i/d_0)^2}$ · $G = -\frac{1}{2}HD^{\circ 2}H$ · $B = 8\pi^2\langle u^2\rangle$ — next in **Part 4**, apply every algorithm directly to SARS-CoV-2 $M^{\text{pro}}$ (`6LU7`) vs. SARS-CoV-1 (`1UK4`) in `Bio.PDB` & `py3Dmol`, setting up **Lecture 8** where protein language models predict contact maps and 3D folds directly from sequence!
