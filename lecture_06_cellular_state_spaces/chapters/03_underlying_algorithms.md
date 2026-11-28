:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
HVGs · Spectral PCA · k-NN Manifolds · Leiden · UMAP · PAGA · DPT
:::

---

# Curriculum & Algorithmic Roadmap

Six geometric and graph-theoretic algorithms turn the clean PBMC3k single-cell matrix (`2,608 × 13,638`) into cell-type lineages, PAGA topology, and continuous differentiation trajectories in Part 4:

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "norm", "label": "CP10k + log1p\n(freeze .raw)",         "x": 20,   "y": 40, "width": 150, "height": 64 },
    { "id": "hvg",  "label": "Binned HVGs +\nZ-Score (±10σ)",        "x": 205,  "y": 40, "width": 160, "height": 64, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "pca",  "label": "Economy SVD\n(Spectral PCA 40D)",      "x": 400,  "y": 40, "width": 165, "height": 64, "style": "fill: #2563EB; color: #fff;" },
    { "id": "knn",  "label": "Fuzzy k-NN Graph\n(obsp affinities)",  "x": 600,  "y": 40, "width": 165, "height": 64, "style": "fill: #16A34A; color: #fff;" },
    { "id": "lei",  "label": "Leiden Modularity\nQ(γ) + Refinement", "x": 800,  "y": 40, "width": 165, "height": 64, "style": "fill: #D97706; color: #fff;" },
    { "id": "umap", "label": "UMAP Cross-Entropy\n+ Wilcoxon + PAGA","x": 1000, "y": 40, "width": 175, "height": 64, "style": "fill: #DC2626; color: #fff;" },
    { "id": "dpt",  "label": "Diffusion DPT +\nKDE Potential U",     "x": 1210, "y": 40, "width": 160, "height": 64, "style": "fill: #0891B2; color: #fff;" }
  ],
  "links": [
    { "from": "norm", "to": "hvg"  },
    { "from": "hvg",  "to": "pca"  },
    { "from": "pca",  "to": "knn"  },
    { "from": "knn",  "to": "lei"  },
    { "from": "lei",  "to": "umap" },
    { "from": "umap", "to": "dpt"  }
  ]
}
:::

:::matrix { .cards cols="1fr 1fr 1fr" gap="18px" }
[[0,0]]{ .purple }
### 1. HVGs & Spectral PCA
Mean-binned dispersion normalization ($z_d$), gene z-scoring with `max_value = 10`, economy SVD ($X = USV^\top$), sign-flip invariance, and low-rank denoising.

[[0,1]]{ .green }
### 2. Manifolds, Leiden & UMAP
Curse of dimensionality, adaptive fuzzy $k$-NN symmetrization, modularity $Q(\gamma)$ & Leiden's connectedness refinement, and UMAP cross-entropy distortion.

[[0,2]]{ .blue }
### 3. Markers, PAGA & Trajectories
Wilcoxon rank-sum DE, PAGA statistical connectivity ratios, multi-scale Diffusion Pseudotime ($M = (I - \tilde{T})^{-1} - I$), and KDE pseudo-potential $U = -\log\hat{p}$.
:::

---
# Mean-Binned HVG Selection & Z-Score Clipping

![Mean-Binned HVG Selection and Z-Score Clipping Leverage](lecture_06_cellular_state_spaces/images/hvg_and_zscore_clipping.png){ width=78 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Why raw variance fails:** Count variance grows with mean expression ($\operatorname{Var} \ge \mu$), so ranking by raw variance picks abundant housekeeping genes (`RPL`/`ACTB`).
- **Mean-binned dispersion (`1,838` HVGs in PBMC3k):** Group genes into 20 expression bins, compute $d_g = \text{Var}\_g / \mu_g$, and standardize within each bin:
$$
z_d(g) = \frac{d_g - \bar{d}\_{\text{bin}}}{\sigma_{d,\text{bin}}} > 0.5
$$

[[0,1]]
- **Gene standardization & outlier clipping (`sc.pp.scale`):**
$$
z_{ij} = \operatorname{clip}\left(\frac{y_{ij} - \mu_j}{\sigma_j}, \ -10, \ 10\right)
$$
  - Prevents high-range genes from dominating PCA by $100\times$ AND stops a single $+40\sigma$ burst cell from hijacking PC1!
  - **Freeze `adata.raw` first** so discarded non-HVG marker genes remain available for dot plots and `score_genes`.
:::

---

# Spectral Decomposition: PCA via Economy SVD

Let $X \in \mathbb{R}^{n \times p}$ be the column-centered, scaled HVG matrix ($n = 2{,}608$ cells, $p = 1{,}838$ HVGs). PCA finds orthogonal axes capturing maximal variance across cells.

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Covariance Eigendecomposition vs. Economy SVD
- Sample gene–gene covariance matrix:
$$
C = \frac{1}{n - 1} X^\top X = V \Lambda V^\top \in \mathbb{R}^{p \times p}
$$
- Forming $X^\top X$ squares the condition number and wastes $O(np^2)$ work. Instead compute **economy SVD** directly:
$$
X = U S V^\top \quad (U \in \mathbb{R}^{n \times r}, \ S = \operatorname{diag}(s_1, \dots, s_r), \ V \in \mathbb{R}^{p \times r})
$$
- **Meta-gene loadings:** columns of $V$ (co-regulated gene programs).
- **Cell coordinates (scores):** $Z = X V = U S$.

[[0,1]]
### Explained Variance Ratio & Reflection Invariance
- Eigenvalues equal scaled squared singular values: $\lambda_i = \dfrac{s_i^2}{n - 1}$.
- **Explained variance ratio** must normalize by **all** singular values *before* slicing top $k$:
$$
r_i = \frac{s_i^2}{\sum_{j=1}^{\min(n, p)} s_j^2}
$$
- **Sign (reflection) invariance:** Flipping any pair of singular vectors leaves the factorization exact:
$$
(-U) S (-V)^\top = U S V^\top
$$
  Hence different SVD solvers (`arpack` vs. `np.linalg.svd`) can produce opposite signs on PC$_i$ while preserving $\lvert r \rvert = 1.000$ and all pairwise distances!
:::

---

# Interactive Demo — 2D Geometric Intuition of PCA

Before projecting 1,838 HVGs into 40 dimensions, explore PCA geometry in 2D: click to place an elongated cloud of points (**Add data point**), switch to **Move PCA coords**, drag the origin to the cloud's center, and rotate either blue axis to watch the 1D projected spread maximize along **Axis 1 (PC1)** while collapsing along **Axis 2 (PC2)**!

<iframe src="lecture_06_cellular_state_spaces/demo/demo_pca.html" title="Interactive 2D PCA Coordinate Frame and Projection Demo" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---
# PCA Denoising & the Curse of Dimensionality

![Spectral PCA Denoising and Curse of Dimensionality](lecture_06_cellular_state_spaces/images/svd_denoising_and_curse.png){ width=89 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Why PCA denoises biology:** Co-regulated gene modules add constructively across hundreds of genes, lifting biological eigenvalues cleanly above the isotropic **Marchenko–Pastur noise bulk** $\lambda_+ = \sigma^2(1 + \sqrt{p/n})^2$. Independent Poisson capture noise spreads thinly across the remaining ~1,800 axes and is discarded at $n_{\text{pcs}} = 40$.

[[0,1]]
- **Curse of dimensionality:** Under independent coordinates as $p \to \infty$, pairwise Euclidean distances concentrate on a thin shell:
$$
\frac{\text{dist}\_{\max} - \text{dist}\_{\min}}{\text{dist}\_{\min}} \to 0, \qquad \frac{\operatorname{SD}(d)}{\bar{d}} \approx \frac{1}{\sqrt{2p}}
$$
  Running $k$-NN in **40D PCA** restores crisp local contrast (`CV ≈ 11%`) vs. `1,838` raw HVGs (`CV ≈ 1.6%`)!
:::

---

# Interactive Lab — HVGs, Z-Score Leverage & Spectral PCA

Explore feature selection and spectral decomposition live: **(1)** in **Tab 1**, toggle **Raw Variance Only** vs. **Mean-Binned HVG** to see abundant housekeeping genes flood raw-variance picks, then slide `max_value` from `40σ` down to `10σ` to tame single-cell leverage; **(2)** in **Tab 2**, adjust independent noise $\sigma$ in the live SVD solver and click **Flip PC1 Sign (`−v₁`)** to verify reflection invariance (`0.00e+00` pairwise distance error)!

<iframe src="lecture_06_cellular_state_spaces/demo/hvg_pca_denoising_lab.html" title="Interactive HVG Selection, Z-Score Clipping and Spectral PCA Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Riemannian Manifolds & Fuzzy $k$-NN Graphs

Gene regulatory networks restrict cells to a curved, low-dimensional **Riemannian manifold** inside $\mathbb{R}^{20{,}000}$. A manifold is *locally* Euclidean — so we approximate it using a weighted $k$-NN graph in 40D PCA space!

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Directed $k$-NN Distances (`obsp['distances']`)
- For each cell $i$, find its $k$ closest cells in 40D PCA (e.g., $k = 10$ in Practice 6).
- Because cell $i$ itself sits at distance $0$, `obsp['distances']` excludes self and stores **$k - 1 = 9$** positive neighbor distances per row.
- **Why not a fixed radius $\epsilon$?** Cell populations vary dramatically in local density (dense CD4 T-cells vs. sparse dendritic cells): fixed $\epsilon$ either fragments sparse clusters or merges dense lineages!

[[0,1]]
### Local Adaptive Kernel & Fuzzy Symmetrization
- **UMAP adaptive bandwidth (`obsp['connectivities']`):** Subtract distance to cell $i$'s nearest neighbor $\rho_i$ (guarantees local connectivity!) and solve for local scale $\sigma_i$ by binary search such that $\sum_j p_{j \mid i} = \log_2 k$:
$$
p_{j \mid i} = \exp\left(-\frac{\max(0, d_{ij} - \rho_i)}{\sigma_i}\right)
$$
- **Fuzzy-set union symmetrization:** Combine asymmetric directed affinities $p_{i \mid j} \ne p_{j \mid i}$ probabilistically:
$$
p_{ij} = p_{i \mid j} + p_{j \mid i} - p_{i \mid j} p_{j \mid i} \in [0, 1]
$$
:::

---

# Manifold Geometry & Louvain's Disconnection Flaw

![Curved Manifold Adaptive k-NN and Louvain vs Leiden Refinement](lecture_06_cellular_state_spaces/images/fuzzy_knn_and_louvain_vs_leiden.png){ width=94% center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Geodesic vs. chord distance:** On curved differentiation manifolds, straight-line Euclidean chords jump across empty, unbiological space, whereas $k$-NN geodesic paths stay faithful to observed cell states — and adaptive $\sigma_i$ stretches cleanly across sparse branches!

[[0,1]]
- **Why Leiden replaced Louvain (Traag et al., 2019):** Greedy Louvain moves a single bridge node `v*` into Community B, leaving Community A split into **two internally disconnected islands** (`0` paths between `A1` and `A2`!). Leiden's refinement step explicitly detects and splits disconnected sub-communities before aggregation.
:::

---

# Graph Community Detection: Modularity $Q(\gamma)$ & CPM

How do we partition the symmetric affinity graph $A_{ij} = p_{ij}$ into biologically meaningful cell types?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Reichardt–Bornholdt Modularity $Q(\gamma)$
- Compares observed intra-cluster edge weight against a degree-preserving configuration null model (weighted degree $k_i = \sum_j A_{ij}$, total weight $m = \frac{1}{2}\sum_{i,j} A_{ij}$):
$$
Q(\gamma) = \frac{1}{2m} \sum_{i,j} \left(A_{ij} - \gamma \frac{k_i k_j}{2m}\right) \delta(c_i, c_j)
$$
- **Constant Potts Model (CPM)** replaces the global degree term with a local density threshold $\gamma$:
$$
\mathcal{H}_{\text{CPM}}(\gamma) = \sum_{i,j} (A_{ij} - \gamma) \delta(c_i, c_j)
$$
  eliminating modularity's *resolution limit* (where tiny rare cell types merge into larger clusters!).

[[0,1]]
### Role of the Resolution Parameter $\gamma$
- **Low $\gamma$ (`0.1–0.3`):** Weak penalty for large communities $\implies$ merges related subsets into broad parental lineages (Lymphoid vs. Myeloid).
- **High $\gamma$ (`1.0–2.0`):** Strong penalty forces communities to split into compact micro-subtypes (CD4 vs. CD8 T-cells, classical vs. non-classical monocytes).
- **No magic universal $\gamma$:** Resolution is validated biologically — stop splitting once new sub-clusters no longer express distinct, reproducible marker genes!
:::

---

# The Leiden Algorithm: Moving, Refinement & Aggregation

Leiden guarantees that every discovered community is **internally connected** and converges to an asymptotically stable partition via a three-phase loop:

:::static-diagram { width=94% center }
{
  "nodes": [
    { "id": "init", "label": "Graph G & Partition P\n(Start: singletons)",     "x": 20,  "y": 35, "width": 170, "height": 62 },
    { "id": "mov",  "label": "Phase 1: Fast Local\nMoving (queue ΔQ > 0)",     "x": 230, "y": 35, "width": 180, "height": 62, "style": "fill: #2563EB; color: #fff;" },
    { "id": "ref",  "label": "Phase 2: Partition\nRefinement (connected!)",    "x": 450, "y": 35, "width": 185, "height": 62, "style": "fill: #16A34A; color: #fff;" },
    { "id": "agg",  "label": "Phase 3: Aggregate\nRefined Super-Nodes",        "x": 675, "y": 35, "width": 180, "height": 62, "style": "fill: #D97706; color: #fff;" },
    { "id": "done", "label": "Converged Partition\n(All Comm. Connected)",     "x": 895, "y": 35, "width": 185, "height": 62, "style": "fill: #7C3AED; color: #fff;" }
  ],
  "links": [
    { "from": "init", "to": "mov"  },
    { "from": "mov",  "to": "ref"  },
    { "from": "ref",  "to": "agg"  },
    { "from": "agg",  "to": "done" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### What Each Phase Does
1. **Fast Local Moving:** Visits only nodes whose neighbors recently changed community (queue-based), moving a node to whichever neighbor community gives the largest $\Delta Q > 0$ — much faster than Louvain's full-graph sweeps!
2. **Partition Refinement:** Inside *each* community $C \in P$, restarts from singletons and merges only **well-connected** sub-nodes using randomized proportional selection — breaking up any disconnected pieces.

[[0,1]]
### Aggregation & Theoretical Guarantees
3. **Graph Aggregation:** Builds the coarser graph from the **refined** partition $P_{\text{ref}}$, while initializing super-node memberships from the **unrefined** partition $P$ (giving Leiden freedom to re-merge or split cleanly!).
- **Strict guarantees (Traag et al., 2019):**
  - Every community is **100% internally connected** after every iteration.
  - Iterating converges to a subset-optimal partition (`n_iterations = 2` in Practice 6 is standard in `igraph`).
:::

---

# Non-Linear Embedding: UMAP Fuzzy Cross-Entropy

To draw the 40D fuzzy $k$-NN graph on a 2D screen, **UMAP** (McInnes et al., 2018) optimizes low-dimensional coordinates $y_1, \dots, y_n \in \mathbb{R}^2$ so their fuzzy affinities match $p_{ij}$:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Heavy-Tailed 2D Kernel & Cross-Entropy Loss
- Define smooth low-dimensional edge probability $q_{ij}$ using a Student-$t$-like heavy-tailed kernel ($a \approx 1.93, b \approx 0.79$ for `min_dist = 0.5`):
$$
q_{ij} = \left(1 + a \lVert y_i - y_j \rVert_2^{2b}\right)^{-1}
$$
- Minimize the **fuzzy-set cross-entropy** over all pairs $(i, j)$:
$$
\mathcal{L} = \sum_{i \ne j} \left[\underbrace{p_{ij} \log\frac{p_{ij}}{q_{ij}}}\_{\text{pull neighbors together}} + \underbrace{(1 - p_{ij}) \log\frac{1 - p_{ij}}{1 - q_{ij}}}\_{\text{push non-neighbors apart}}\right]
$$

[[0,1]]
### Negative-Sampling SGD & Seed Sensitivity
- Evaluating all $O(n^2)$ non-neighbor pairs is prohibitive, so UMAP uses **Stochastic Gradient Descent (SGD)** with **negative sampling** (updating each true graph edge plus ~5 random repulsive pairs per step) initialized from a spectral layout.
- **Why `random_state` changes UMAP but NOT Leiden:**
  - `sc.tl.leiden` operates strictly on `adata.obsp['connectivities']` — **zero** dependence on 2D coordinates!
  - Changing `random_state` in `sc.tl.umap` changes the random negative samples $\implies$ clusters rotate/shift in 2D while cluster memberships stay **100% identical**.
:::

---
# Why UMAP Distorts Global Distances & Densities

![UMAP Cross-Entropy Forces and Global Distortion](lecture_06_cellular_state_spaces/images/umap_cross_entropy_and_distortion.png){ width=93 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **The flat long-range plateau:** Once non-neighbors separate beyond $d \gtrsim 2$, the repulsive gradient drops near zero — so two distant lineages experience almost no metric restoring force and settle at **arbitrary inter-cluster gaps**!
- **Density equalization:** Because each cell's local bandwidth $\sigma_i$ normalizes its $k$-NN neighborhood to identical weight sum $\log_2 k$, UMAP **expands tight clusters** and **compresses diffuse clusters** into uniform visual blobs.

[[0,1]]
- **Golden rules for interpreting UMAP:**
  - **Valid:** local neighbor proximity, cluster separation, and continuous unbroken trajectories (like monocyte maturation).
  - **Invalid:** measuring Euclidean distances between disconnected islands, comparing cluster areas/diameters, or running $k$-means/clustering directly on 2D UMAP coordinates!
:::

---

# Interactive Lab — Curse of Dim., Fuzzy $k$-NN, Leiden & UMAP

Test manifold geometry and graph clustering live: **(1)** in **Tab 1**, jump between **3D**, **40D PCA**, and **1,838 HVGs** to watch pairwise Euclidean distances collapse into a spike, and tune $k$ to inspect adaptive $\sigma_i$ on a dense-vs-sparse arc; **(2)** in **Tab 2**, sweep **Leiden resolution $\gamma$** (`0.1 → 2.2`), trigger the **Louvain Bridge Flaw**, and click **Shuffle UMAP `random_state`** to verify that Leiden communities stay 100% invariant!

<iframe src="lecture_06_cellular_state_spaces/demo/knn_leiden_umap_lab.html" title="Interactive Curse of Dimensionality, Fuzzy k-NN, Leiden and UMAP Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---
# Marker Discovery: Wilcoxon Rank-Sum & Gene-Set Scoring

Turning numerical Leiden IDs (`0, 1, 2, …`) into biological cell types (`CD4 T`, `B`, `NK`, `CD14 Mono`, `DC`) combines non-parametric testing with curated gene-set scoring:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Data-Driven Markers: Wilcoxon Rank-Sum Test
- For every gene $g$, compare expression in cluster $C$ ($n_1$ cells) vs. all other cells ($n_2$ cells).
- Single-cell `log1p` counts are **heavily zero-inflated and non-Gaussian** $\implies$ classical $t$-tests are distorted by outliers and skewness.
- **Mann–Whitney / Wilcoxon rank-sum** replaces raw values with joint ranks $R_i$ (`method='wilcoxon'`):
$$
U = \sum_{i \in C} R_i - \frac{n_1(n_1 + 1)}{2}, \qquad z = \frac{U - n_1 n_2 / 2}{\sigma_{U,\text{tie}}}
$$
  followed by tie-variance correction (crucial for tied zeros!) and Benjamini–Hochberg FDR adjustment (`padj`).

[[0,1]]
### Curated Lineage Scoring (`sc.tl.score_genes`)
- Given a canonical marker panel $\mathcal{G}$ (e.g., `MS4A1, CD79A, CD79B` for B-cells) evaluated on **`adata.raw`**:
  1. Bin all expressed genes into 25 mean-expression bins.
  2. For each marker $g \in \mathcal{G}$, sample $n_{\text{ctrl}} = 50$ random **control genes** from the *same* expression bin.
  3. Score cell $i$ as the excess mean expression above matched background:
$$
\text{Score}\_i(\mathcal{G}) = \frac{1}{\lvert\mathcal{G}\rvert} \sum_{g \in \mathcal{G}} y_{ig} - \frac{1}{\lvert\mathcal{B}\rvert} \sum_{b \in \mathcal{B}} y_{ib}
$$
  removing global sequencing-depth and expression-tier bias!
:::

---

# Coarse-Grained Manifold Topology: PAGA

Since 2D UMAP warps inter-cluster distances, how do we test which cell types are genuinely adjacent on the high-dimensional manifold?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Partition-Based Graph Abstraction (Wolf et al., 2019)
- **PAGA** collapses the full $n$-cell $k$-NN graph to a coarse graph with one super-node per cell type ($I, J$).
- Let $E_{\text{obs}}(I, J)$ be the observed number of $k$-NN edges crossing between clusters $I$ (size $n_I$) and $J$ (size $n_J$).
- Under a random graph null model with total edges $\lvert E \rvert$ over $n$ cells, the expected cross-edges are:
$$
\mathbb{E}\_{\text{rand}}[E(I, J)] = \lvert E \rvert \cdot \frac{2 \ n_I n_J}{n(n - 1)}
$$
- **PAGA connectivity score:**
$$
e_{IJ} = \frac{E_{\text{obs}}(I, J)}{\mathbb{E}\_{\text{rand}}[E(I, J)]}
$$
  (scaled and thresholded by an exact binomial/normal $p$-value to prune spurious noise edges).

[[0,1]]
### Reconciling Topology with PBMC Biology
- **Strong PAGA edges ($e_{IJ} \gg 1$):** Reflect genuine transcriptional continuity on the 40D manifold:
  - `CD14⁺ Monocytes ↔ FCGR3A⁺ Monocytes` (continuous peripheral maturation!)
  - `CD4 T ↔ CD8 T ↔ NK` (shared lymphoid/cytotoxic transcriptional programs)
- **Zero PAGA edges ($e_{IJ} = 0$):** Lineages separated by early hematopoietic divergence (`B-cells` vs. `Monocytes`) share **no** manifold bridge, regardless of how close UMAP happens to draw them!
- **Bonus:** PAGA positions serve as a topology-preserving initialization for `sc.tl.umap(init_pos='paga')`.
:::

---
# Diffusion Maps & Diffusion Pseudotime (DPT)

Discrete clusters break down during **continuous differentiation** (like classical `CD14++` $\to$ non-classical `FCGR3A++` monocytes). **Diffusion Maps** (Coifman & Lafon, 2006) and **DPT** (Haghverdi et al., 2016) order cells by random-walk geometry!

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Random-Walk Transition Matrix & Spectral Decay
- From symmetric affinity matrix $A$ and degree matrix $D_{ii} = \sum_j A_{ij}$, define the Markov random-walk transition matrix:
$$
T = D^{-1} A, \qquad \tilde{T} = D^{-1/2} A D^{-1/2}
$$
- Because $\tilde{T}$ is symmetric with the same real eigenvalues as $T$:
$$
1 = \lambda_0 > \lvert\lambda_1\rvert \ge \lvert\lambda_2\rvert \ge \dots
$$
- The stationary trivial mode ($\lambda_0 = 1$) gives equilibrium degree density; slow-decaying modes ($\lambda_1, \lambda_2, \dots \approx 1$) capture **global differentiation axes** while fast-decaying modes ($\lambda_k \approx 0$) are local noise!

[[0,1]]
### Multi-Scale Diffusion Pseudotime ($M$ Matrix)
- Rather than choosing one arbitrary walk length $t$, DPT sums transition probabilities over **all walk lengths** $t = 1, 2, \dots, \infty$ on the non-trivial operator $\tilde{T}^* = \tilde{T} - \psi_0 \psi_0^\top$:
$$
M = \sum_{t=1}^{\infty} (\tilde{T}^\*)^t = (I - \tilde{T}^*)^{-1} - I = \sum_{k=1}^{m} \frac{\lambda_k}{1 - \lambda_k} \psi_k \psi_k^\top
$$
- **Diffusion pseudotime distance** from chosen root cell $i_{\text{root}}$:
$$
d_{\text{dpt}}(i_{\text{root}}, j) = \bigl\lVert M(i_{\text{root}}, \cdot) - M(j, \cdot) \bigr\rVert_2
$$
- Cells in disconnected components (zero random-walk path from root) receive **$d_{\text{dpt}} = \infty$** (`np.inf` $\to$ `np.nan` in Practice 6!).
:::

---
# Diffusion Trajectories & Waddington KDE Landscape

![Multi-Scale Diffusion Pseudotime and Waddington KDE Pseudo-Potential](lecture_06_cellular_state_spaces/images/diffusion_dpt_and_waddington_kde.png){ width=80 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Continuous vs. discrete immune states:** Peripheral `CD14⁺ → FCGR3A⁺` monocytes actively mature in blood, forming an unbroken arc where DPT orders cells smoothly (`CD14`/`LYZ` fall while `FCGR3A`/`MS4A7` rise), whereas mature B-cells and T-cells diverged early in bone marrow/thymus and sit in separate basins.

[[0,1]]
- **KDE Waddington pseudo-potential:** Smooth 2D density with Gaussian kernel bandwidth $h$:
$$
\hat{p}(\mathbf{x}) = \frac{1}{n} \sum_{i=1}^n K_h(\mathbf{x} - \mathbf{x}\_i), \qquad U(\mathbf{x}) = -\log\bigl(\hat{p}(\mathbf{x}) + \epsilon\bigr)
$$
  High-density stable fates become **valleys**, continuous maturation forms a **saddle pass**, and empty transition zones become **ridges** *(illustrative on UMAP — coordinates are non-physical!)*.
:::

---

# Interactive Lab — PAGA, DPT & Waddington Pseudo-Potential

Explore random walks and energy landscapes live: **(1)** in **Tab 1**, move the **Root Cell** along the monocyte arc, toggle between the finite-time wave $T^t$ and **Multi-Scale DPT**, verify that zero random-walk mass leaks into disconnected B-cells, and check the **PAGA connectivity ratios**; **(2)** in **Tab 2**, adjust **KDE Bandwidth $h$** and **Transition Bridge Cells** to watch the epigenetic barrier $\Delta U$ rise or vanish!

<iframe src="lecture_06_cellular_state_spaces/demo/paga_dpt_waddington_lab.html" title="Interactive PAGA, Diffusion Pseudotime and Waddington KDE Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Algorithmic Complexity & Scaling of Manifold Learning

Why can Scanpy run PCA, $k$-NN, Leiden, UMAP, PAGA, and DPT on `2,608` cells in ~3 seconds — and scale to 1,000,000-cell atlases?

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Stage | Naive Algorithm | Modern Single-Cell Algorithm | Complexity & Scaling |
| :--- | :--- | :--- | :--- |
| **Spectral PCA** | Dense covariance $X^\top X$ + full eigendecomp. $O(np^2 + p^3)$ | Truncated Lanczos (`arpack`) / randomized SVD on HVGs | $O(n \cdot p_{\text{HVG}} \cdot k)$ ($k = 50 \ll p = 1{,}838$); **< 0.3 s** |
| **$k$-NN Graph** | Brute-force all-pairs $O(n^2 d)$ in 2,000D | Approximate **NN-Descent** (`pynndescent`) in 40D PCA | Empirical $\sim O(n^{1.14} d)$ + sparse $O(nk)$ storage (`obsp`) |
| **Community Detection** | Louvain full-graph sweeps (can disconnect!) | Queue-based **Leiden** (local moving + refinement + agg.) | Near-linear $O(\lvert E \rvert) = O(nk)$; **< 0.1 s** on PBMC3k |
| **2D Embedding (UMAP)** | Exact $O(n^2)$ repulsive force sum | Spectral init + **negative-sampling SGD** | $O(n \cdot k \cdot n_{\text{epochs}})$; linear in cell count $n$ |
| **Diffusion Maps & DPT** | Dense $n \times n$ matrix inversion $O(n^3)$ | Top-$m$ sparse symmetric ARPACK eigenvectors ($m \approx 15$) | $O(\lvert E \rvert \cdot m + n m^2)$ — avoids forming dense $M$! |
:::

:::matrix { .cards cols="1fr 1fr" gap="22px" }
[[0,0]]{ .purple }
### Why Truncating to HVGs + 40 PCs Multiplies Speed
Subsetting `13,638 → 1,838` HVGs cuts SVD memory/time by **$7.5\times$** before dense z-scoring; projecting `1,838 → 40` PCs accelerates $k$-NN search by **$46\times$** while simultaneously curing the curse of dimensionality!

[[0,1]]{ .green }
### Sparse Graph Backbone Everywhere
Once `sc.pp.neighbors` builds the sparse $n \times n$ matrix (`~99.5%` zeros), **Leiden**, **UMAP**, **PAGA**, and **Diffusion Maps / DPT** all operate strictly on sparse edges $O(nk)$ — never $O(n^2)$!
:::

---

# Summary & Bridge to Part 4 (Industrial Tools)

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### HVGs & Spectral PCA
#### Linear Denoising
- CP10k + `log1p`; freeze `adata.raw`
- 20-bin normalized dispersion $z_d$
- Z-score + clip `max_value = 10`
- Economy SVD: $X = U S V^\top$
- Scores $Z = US$; sign-flip invariant
- 40D PCA beats curse of dim.

[[0,1]]{ .green }
### $k$-NN, Leiden & UMAP
#### Manifold & Clusters
- `distances` ($k - 1$ non-self)
- Adaptive kernel $\sigma_i$ + fuzzy union
- Modularity $Q(\gamma)$ & CPM
- Leiden refinement = connected ✓
- UMAP SGD: local ✓, global ✗
- Seed alters UMAP, not Leiden

[[0,2]]{ .blue }
### Lineages & Kinetics
#### PAGA · DPT · KDE
- Wilcoxon rank-sum for zero-inflated DE
- `score_genes` control-binned on `.raw`
- PAGA ratio $e_{IJ}$ = true graph adjacency
- DPT Neumann sum $M = (I - \tilde{T}^*)^{-1} - I$
- Monocytes: `CD14⁺ → FCGR3A⁺`
- KDE $U = -\log(\hat{p} + \epsilon)$ basins
:::

> **Core Equations at a Glance:** $X = USV^\top$ · $r_i = S_i^2 / \sum_j S_j^2$ · $p_{j\mid i} = \exp(-\max(0, d_{ij} - \rho_i)/\sigma_i)$ · $p_{ij} = p_{i\mid j} + p_{j\mid i} - p_{i\mid j}p_{j\mid i}$ · $Q(\gamma) = \frac{1}{2m}\sum_{i,j}(A_{ij} - \gamma \frac{k_i k_j}{2m})\delta(c_i, c_j)$ · $M = (I - \tilde{T}^*)^{-1} - I$ · $U(\mathbf{x}) = -\log(\hat{p}(\mathbf{x}) + \epsilon)$
