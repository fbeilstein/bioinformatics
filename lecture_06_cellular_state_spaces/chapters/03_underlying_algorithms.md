:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
HVGs · Spectral PCA · k-NN Manifolds · Leiden · UMAP · PAGA · DPT
:::

---

# Curriculum

- **Normalization, Highly Variable Gene (HVG) Selection & Z-Scoring** — CP10k scaling and $\ln(1+x)$ variance stabilization; why raw variance selects abundant housekeeping genes; mean-binned normalized dispersion ($d = \text{Var}/\mu$, within-bin $z$-scores); gene standardization ($z_{ij}$) and outlier clipping (`max_value=10`) to prevent single-cell leverage in PCA.
- **Spectral Decomposition (PCA via Economy SVD)** — gene–gene covariance $C = \frac{1}{n-1}X^\top X = V\Lambda V^\top$, economy SVD $X = U S V^\top$, meta-gene loadings $V$ vs. cell scores $Z = US$, explained variance ratio over total variance, eigenvector reflection (sign-flip) invariance, and low-rank denoising (separating correlated biological modules from isotropic noise).
- **Manifold Approximation & $k$-Nearest Neighbor ($k$-NN) Graphs** — Riemannian manifold hypothesis in $\mathbb{R}^{20{,}000}$, curse of dimensionality (distance concentration as $p \to \infty$), local Euclidean neighborhoods in denoised PCA space, directed $k-1$ distance graphs vs. UMAP fuzzy-set symmetrized connectivities.
- **Graph Community Detection (Modularity & Leiden)** — modularity $Q(\gamma)$ and Constant Potts Model (CPM) with resolution parameter $\gamma$; Louvain's greedy merge flaw producing internally disconnected communities; Leiden's three-phase loop (local moving → **refinement** → aggregation) guaranteeing connected partitions.
- **Non-Linear Embedding (UMAP)** — fuzzy simplicial sets, heavy-tailed low-dimensional Student-$t$-like kernel $q_{ij}$, cross-entropy minimization via stochastic gradient descent (SGD), why local neighborhoods survive while global inter-cluster distances and densities are distorted, and seed sensitivity vs. graph invariance.
- **Marker Discovery, Graph Abstraction (PAGA) & Diffusion Trajectories (DPT)** — non-parametric Wilcoxon rank-sum testing (cluster vs. rest) and gene-set scoring, PAGA statistical inter-cluster connectivity ratio, Diffusion Maps & Diffusion Pseudotime on the random-walk transition matrix $T = D^{-1}A$, and KDE pseudo-potential $U(x,y) = -\log\hat{p}(x,y)$.

---
# Key Math Notes — HVGs, Spectral PCA & $k$-NN Manifolds

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Normalization, HVGs & Gene Scaling
- **CP10k + log1p**:
  $y_{ij} = \ln\left(1 + 10^4 \cdot \frac{x_{ij}}{\sum_k x_{ik}}\right)$
- **Binned HVG dispersion**: group genes into expression bins, compute $d_g = \text{Var}\_g / \mu_g$, and standardize within each bin:
  $$z_d(g) = {d_g - \bar{d}\_{\text{bin}}}/{\sigma_{d,\text{bin}}}$$
- **Gene z-scoring** (clipped at $\pm 10$): $z_{ij} = \text{clip}\bigl((y_{ij} - \mu_j)/\sigma_j, -10, 10\bigr)$.

### Spectral PCA via Economy SVD
- Centered matrix $X \in \mathbb{R}^{n \times p}$ ($n$ cells, $p$ HVGs):
  $$X = U S V^\top, \qquad C = \frac{1}{n-1}X^\top X = V \Lambda V^\top$$
- **Scores**: $Z = XV = US$ · **Eigenvalues**: $\lambda_i = \frac{S_i^2}{n - 1}$.
- **Explained variance ratio** (over **all** singular values):
  $$r_i = \frac{S_i^2}{\sum_{j=1}^{\min(n,p)} S_j^2}, \qquad (-U)S(-V)^\top = USV^\top$$

[[0,1]]
### Curse of Dimensionality & Local Graphs
- As dimension $p \to \infty$ under independent noise:
  $$\frac{\text{dist}\_{\max} - \text{dist}\_{\min}}{\text{dist}\_{\min}} \to 0$$
  Global Euclidean distances concentrate and cut off the manifold; restricting to the top $d \approx 40$ PCs restores local geometry.

### $k$-NN Distances & Fuzzy Connectivities
- **Directed distances** (`obsp['distances']`): each row $i$ stores Euclidean distances $d_{ij}$ to its $k - 1$ closest *other* cells (excluding self).
- **Local adaptive kernel** (distance to 1st neighbor $\rho_i$, bandwidth $\sigma_i$ solved so $\sum_j p_{j|i} = \log_2 k$):
  $$p_{j|i} = \exp\left(-\frac{\max(0, d_{ij} - \rho_i)}{\sigma_i}\right)$$
- **Fuzzy-set symmetrization** (`obsp['connectivities']`):
  $$p_{ij} = p_{i|j} + p_{j|i} - p_{i|j} p_{j|i} \in [0, 1]$$
:::

---
# Key Math Notes — Leiden, UMAP, PAGA & Diffusion Pseudotime

:::matrix { cols=50/50 }
[[0,0]]
### Leiden Modularity Optimization
- Given symmetric affinity graph $A_{ij} = p_{ij}$, weighted degree $k_i = \sum_j A_{ij}$, total weight $m = \tfrac{1}{2}\sum_{i,j} A_{ij}$:
  $$Q(\gamma) = \frac{1}{2m} \sum_{i,j} \left(A_{ij} - \gamma \frac{k_i k_j}{2m}\right)\delta(c_i, c_j)$$
- Higher resolution $\gamma$ penalizes large communities → finer clusters; Leiden's **refinement phase** splits poorly connected subsets before aggregation, guaranteeing every community is internally connected.

### UMAP Fuzzy Cross-Entropy
- Low-dimensional heavy-tailed affinity ($y_i \in \mathbb{R}^2$, default $a \approx 1.93, b \approx 0.79$):
  $q_{ij} = \left(1 + a \lVert y_i - y_j \rVert^{2b}\right)^{-1}$
- Minimized by SGD over high-dimensional affinities $p_{ij}$:
  $$\mathcal{L} = \sum_{i \ne j} \left[p_{ij}\log\frac{p_{ij}}{q_{ij}} + (1 - p_{ij})\log\frac{1 - p_{ij}}{1 - q_{ij}}\right]$$

[[0,1]]
### Markers & Coarse-Grained PAGA
- **Wilcoxon rank-sum**: non-parametric test comparing ranks of gene $g$ inside cluster $C$ vs. all other cells (with BH FDR adjustment).
- **PAGA connectivity** between clusters $I, J$:
  $e_{IJ} = \frac{E_{\text{obs}}(I, J)}{\mathbb{E}\_{\text{rand}}[E(I, J)]}$
  Quantifies true manifold adjacency free of 2D UMAP distortion.

### Diffusion Pseudotime & Pseudo-Potential
- Random-walk transition matrix $T = D^{-1}A$ ($D_{ii} = \sum_j A_{ij}$) with stationary eigenvectors $\psi_k$ and eigenvalues $1 = \lambda_0 > \lvert\lambda_1\rvert \ge \dots$:
- **Multi-scale diffusion distance** (summing over all path lengths $t \ge 1$ as $\sum_{t=1}^{\infty} \tilde{T}^t$ ):
$$
  M = (I - \tilde{T})^{-1} - I,\\;\\; d_{\text{dpt}}(x, y) = \lVert M(x,\cdot) - M(y,\cdot) \rVert\_2
$$
- **KDE Waddington pseudo-potential**:
  $$U(x, y) = -\log\bigl(\hat{p}(x, y) + \epsilon\bigr)$$
:::
