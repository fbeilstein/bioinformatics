:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Dihedral Geometry · Kabsch SVD · Contact Maps · Debye–Waller
:::

---

# Curriculum

- **Backbone Internal Coordinates & Dihedral Angles ($\phi, \psi, \omega$)** — planar peptide-bond resonance locking $\omega \approx 180^\circ$ (*trans*), computing signed torsion angles from four consecutive atomic positions via bond-plane normals and `atan2`, and stereochemical Ramachandran basins (including Glycine's achiral exception).
- **Rigid-Body Motion ($\text{SE}(3)$) & the Orthogonal Procrustes Problem** — why naive coordinate subtraction fails across arbitrary experimental reference frames; formulating optimal superposition as minimizing RMSD over proper rotations $R \in \text{SO}(3)$ and translations $\mathbf{t} \in \mathbb{R}^3$.
- **The Kabsch Algorithm (SVD Derivation)** — decoupling translation via centroid centering ($\mathbf{t}^* = \mathbf{c}_Q - R\mathbf{c}_P$), reducing RMSD minimization to maximizing $\operatorname{tr}(R C)$ on the $3 \times 3$ cross-covariance matrix $C = P_c^T Q_c = U S V^T$, enforcing right-handed chirality ($\det R = +1$) to forbid L→D mirror reflections, reconciling row vs. column multiplication conventions, and contrasting SVD in Kabsch vs. PCA.
- **Correspondence, Frame-Shift Pitfalls & Length-Normalized Metrics** — why naive list truncation fails when crystal structures have missing terminal or loop residues; residue-number/sequence pairing vs. iterative structural alignment; limitations of RMSD on flexible tails and length-independent scores (TM-score, GDT-TS).
- **$\text{SE}(3)$-Invariant Distance Geometry (Contact Maps)** — pairwise $C\alpha$ Euclidean distance matrix $D_{ij} = \lVert \mathbf{x}_i - \mathbf{x}_j \rVert_2$, invariance to global rotations and translations, geometric signatures of $\alpha$-helices ($i \to i+4$ diagonal band) vs. parallel and anti-parallel $\beta$-sheets, and reconstructing 3D embeddings from $D$ (classical MDS / distance geometry).
- **Crystallographic Disorder & Thermal Mobility (B-Factors)** — Debye–Waller attenuation of Bragg diffraction intensities, isotropic mean-square displacement $B = 8\pi^2\langle u^2 \rangle$, and why mobile surface loops undergo destructive phase interference across periodic crystal unit cells.

---

# Key Math Notes — Dihedrals, Contact Maps & Debye–Waller

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Signed Backbone Dihedral Angles ($\phi, \psi$)
Given 4 sequential atoms $\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3, \mathbf{r}_4$ with bond vectors $\mathbf{b}_k = \mathbf{r}_{k+1} - \mathbf{r}_k$:
- Plane normals: $\mathbf{n}_1 = \mathbf{b}_1 \times \mathbf{b}_2$, $\mathbf{n}_2 = \mathbf{b}_2 \times \mathbf{b}_3$, unit axis $\hat{\mathbf{u}} = \mathbf{b}_2 / \lVert \mathbf{b}_2 \rVert$.
- Signed torsion angle in $(-\pi, \pi]$:
  $$\theta = \operatorname{atan2}\bigl((\mathbf{n}_1 \times \mathbf{n}_2) \cdot \hat{\mathbf{u}}, \quad \mathbf{n}_1 \cdot \mathbf{n}_2\bigr)$$
- **$\phi$**: $(C_{i-1}, N_i, C\alpha_i, C_i)$ · **$\psi$**: $(N_i, C\alpha_i, C_i, N_{i+1})$ · **$\omega$**: $(C\alpha_i, C_i, N_{i+1}, C\alpha_{i+1}) \approx 180^\circ$.
- Terminal residues lack $C_{i-1}$ or $N_{i+1} \implies$ undefined (`None`).

### Ramachandran Free-Energy Basins
- **$\beta$-sheet**: $\phi \in [-140^\circ, -120^\circ]$, $\psi \in [+130^\circ, +160^\circ]$.
- **Right-handed $\alpha$-helix**: $\phi \approx -60^\circ$, $\psi \approx -45^\circ$.
- **Left-handed $\alpha$**: $\phi \approx +60^\circ, \psi \approx +40^\circ$ (sterically accessible mainly to achiral Glycine).

[[0,1]]
### $\text{SE}(3)$-Invariant Distance Matrix
$$D_{ij} = \lVert \mathbf{x}\_i - \mathbf{x}\_j \rVert_2 = \lVert (R\mathbf{x}\_i + \mathbf{t}) - (R\mathbf{x}\_j + \mathbf{t}) \rVert_2$$
- **$\alpha$-helix**: $i \to i+4$ H-bond gives $D_{i, i+4} \approx 5.5\text{ Å}$ (thick band along main diagonal).
- **Parallel $\beta$-sheet**: $(i, j), (i+1, j+1)$ → stripe **parallel** to diagonal ($j - i = \text{const}$).
- **Anti-parallel $\beta$-sheet**: $(i, j), (i+1, j-1)$ → stripe **perpendicular** to diagonal ($i + j = \text{const}$).

### Debye–Waller Crystallographic Mobility
- Atomic vibration attenuates X-ray structure factor amplitudes by $\exp(-B \sin^2\theta / \lambda^2)$:
  $$B = 8\pi^2 \langle u^2 \rangle \qquad \Longrightarrow \qquad \sqrt{\langle u^2 \rangle} = \sqrt{\frac{B}{8\pi^2}}$$
- Large $\langle u^2 \rangle$ (mobile loops) causes destructive interference across crystal unit cells, erasing electron density.
:::

---

# Key Math Notes — The Kabsch / Orthogonal Procrustes Derivation

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Objective & Centroid Decoupling
Given paired $N \times 3$ coordinate matrices $P, Q$ (rows $\mathbf{p}\_i^T, \mathbf{q}\_i^T$):
$$\text{RMSD}^2(R, \mathbf{t}) = \frac{1}{N}\sum_{i=1}^N \lVert R\mathbf{p}\_i + \mathbf{t} - \mathbf{q}\_i \rVert_2^2$$
- Setting $\nabla_{\mathbf{t}} \text{RMSD}^2 = \mathbf{0}$ yields the exact optimal translation:
  $$\mathbf{t}^* = \mathbf{c}\_Q - R\mathbf{c}\_P, \qquad \mathbf{c}\_P = \frac{1}{N}\sum_{i=1}^N \mathbf{p}\_i$$
- Centering $P_c = P - \mathbf{1}\mathbf{c}\_P^T$, $Q_c = Q - \mathbf{1}\mathbf{c}\_Q^T$ removes $\mathbf{t}$ and expands the Frobenius objective:
  $$\lVert P_c R^T - Q_c \rVert_F^2 = \lVert P_c \rVert_F^2 + \lVert Q_c \rVert_F^2 - 2\operatorname{tr}\bigl(R P_c^T Q_c\bigr)$$
- Minimizing RMSD over $R \in \text{SO}(3)$ is equivalent to **maximizing** $\operatorname{tr}(R C)$ where $C = P_c^T Q_c \in \mathbb{R}^{3 \times 3}$.

[[0,1]]
### SVD Solution & Chirality Enforcement
Decompose cross-covariance $C = P_c^T Q_c = U S V^T$ ($S = \operatorname{diag}(s_1 \ge s_2 \ge s_3 \ge 0)$):
$$\operatorname{tr}(R C) = \operatorname{tr}(R U S V^T) = \operatorname{tr}\bigl((V^T R U) S\bigr)$$
- Since $W = V^T R U \in \text{O}(3)$, each diagonal entry obeys $W_{kk} \le 1$, maximized by $W = \operatorname{diag}(1, 1, d)$ with:
  $$d = \operatorname{sign}\bigl(\det(V U^T)\bigr) \in \{+1, -1\}$$
  so $\det(R) = +1$ (forbids chirality-inverting reflections!).
- **Two equivalent conventions**:
  - Column-vector: $R_{\text{col}} = V \operatorname{diag}(1, 1, d) U^T \implies P_{\text{rot}} = P_c R_{\text{col}}^T$.
  - Row-vector: $R_{\text{row}} = U \operatorname{diag}(1, 1, d) V^T \implies P_{\text{rot}} = P_c R_{\text{row}}$.
- **Optimal RMSD**: $\sqrt{\frac{1}{N}\sum_{i=1}^N \lVert \mathbf{p}\_{\text{rot},i} - \mathbf{q}\_{c,i} \rVert_2^2}$ (sum over $x,y,z$, divide by $N$, **not** $3N$).
:::
