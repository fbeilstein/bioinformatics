:::titlepage
[[top]]
Part 4

[[title]]
Industrial Tools

[[subcaption]]
Scanpy · igraph · UMAP · Leiden · PAGA · DPT
:::

---

# Curriculum

- **Preprocessing & Normalization** — CP10k library size scaling, $\log(1 + x)$ variance stabilization, saving raw counts in `adata.raw`.
- **Feature Selection (HVGs)** — mean–dispersion relationships, identifying ~1,800 informative genes, gene z-scoring via `sc.pp.scale`.
- **Spectral Decomposition (PCA)** — Singular Value Decomposition ($X = USV^T$), meta-genes (loadings), cell scores, PCA as a biological denoiser, the scree/elbow plot.
- **Approximating the Manifold ($k$-NN Graph)** — Riemannian manifold hypothesis, curse of dimensionality in 2,000D vs. local Euclidean geometry in 40D PCA space, `adata.obsp['distances']` and `connectivities`.
- **Graph Partitioning (Leiden Clustering)** — modularity and Constant Potts Model (CPM) optimization, fixing Louvain's disconnected communities, the resolution parameter $\gamma$.
- **Non-Linear Embedding (UMAP)** — fuzzy simplicial sets, cross-entropy minimization, global vs. local distance distortions.
- **Marker Genes & Lineage Annotation** — data-driven Wilcoxon rank-sum differential expression, canonical immune marker panels, cell-type scoring, dot plots.
- **Coarse-Grained Topology (PAGA)** — lineage connectivity testing, reconciling graph topology with 2D projections.
- **Continuous State Transitions (Diffusion Pseudotime)** — random walks on the $k$-NN graph, transition matrices, monocyte maturation ($CD14^+ \to FCGR3A^+$).
- **Data-Driven Waddington Landscape** — pseudo-potentials $U(x, y) = -\log \hat{p}(x, y)$, attractor basins, differentiation ridges.

---

# Roadmap

The modern single-cell manifold learning, clustering, and trajectory pipeline:

:::static-diagram { width=95% center }
{
  "nodes": [
    { "id": "clean", "label": "Clean AnnData\n(from Lecture 5)",        "x": 40,   "y": 80,  "width": 150, "height": 60 },
    { "id": "norm",  "label": "normalize + log1p\nadata.raw frozen",    "x": 220,  "y": 80,  "width": 160, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "hvg",   "label": "Highly Variable Genes\n~1,800 HVGs",     "x": 410,  "y": 80,  "width": 160, "height": 60 },
    { "id": "pca",   "label": "Spectral PCA\n50 components",            "x": 600,  "y": 80,  "width": 150, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "knn",   "label": "k-NN Graph\nk=10 in 40D PCA",            "x": 780,  "y": 80,  "width": 150, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "leiden","label": "Leiden Clustering\nCommunity Detection", "x": 960,  "y": 80,  "width": 150, "height": 60, "style": "fill: #FF9800; color: #fff;" },

    { "id": "umap",  "label": "UMAP 2D Projection\nVisualization",      "x": 960,  "y": 250, "width": 160, "height": 60, "style": "fill: #009688; color: #fff;" },
    { "id": "markers","label": "Marker Genes\nWilcoxon DE + Scores",    "x": 770,  "y": 250, "width": 160, "height": 60, "style": "fill: #E91E63; color: #fff;" },
    { "id": "types", "label": "Lineage Annotation\nT, B, NK, Mono, DC", "x": 580,  "y": 250, "width": 160, "height": 60, "style": "fill: #673AB7; color: #fff;" },
    { "id": "paga",  "label": "PAGA Topology\nInter-cluster edges",     "x": 390,  "y": 250, "width": 160, "height": 60, "style": "fill: #3F51B5; color: #fff;" },
    { "id": "dpt",   "label": "Diffusion Pseudotime\nContinuous Maturation", "x": 190, "y": 250, "width": 170, "height": 60, "style": "fill: #795548; color: #fff;" },
    { "id": "wadd",  "label": "Waddington Landscape\nPotential U = -log p", "x": 10, "y": 250, "width": 150, "height": 60 }
  ],
  "links": [
    { "from": "clean", "to": "norm" },
    { "from": "norm",  "to": "hvg" },
    { "from": "hvg",   "to": "pca" },
    { "from": "pca",   "to": "knn" },
    { "from": "knn",   "to": "leiden" },
    { "from": "leiden","to": "umap" },
    { "from": "umap",  "to": "markers" },
    { "from": "markers","to": "types" },
    { "from": "types", "to": "paga" },
    { "from": "paga",  "to": "dpt" },
    { "from": "dpt",   "to": "wadd" }
  ]
}
:::

> In Lecture 5, we turned raw droplets into clean, single cells. In this lecture, we map their collective geometry: uncovering cell types, lineages, and developmental transitions.

---

# The Manifold Learning Toolbox

Analyzing cellular state spaces requires combining linear spectral decomposition with non-linear graph mechanics:

:::matrix { .cards cols="1fr 1fr 1fr 1fr" gap="20px" margin-top: 20px; }
[[0,0]]{ .purple }
### Scanpy Preprocessing
`sc.pp.normalize_total`, `sc.pp.log1p`, `sc.pp.highly_variable_genes`, and `sc.pp.scale`: standardizing variance across 20,000 gene dimensions.

[[0,1]]{ .blue }
### Spectral PCA
`sc.pp.pca`: truncated SVD finding orthogonal axes of maximal biological co-variation while discarding independent technical noise.

[[0,2]]{ .green }
### igraph & Leiden
`sc.tl.leiden`: guarantees well-connected community partitioning by optimizing modularity on the high-dimensional $k$-NN graph.

[[0,3]]{ .orange }
### UMAP & PAGA & DPT
`sc.tl.umap`, `sc.tl.paga`, `sc.tl.dpt`: non-linear visualization, coarse-grained topology, and random-walk diffusion pseudotime.
:::

```bash
# Graph libraries required for Leiden community detection
pip install -q scanpy anndata igraph leidenalg scikit-learn
```

---

# The Epigenetic Landscape in 20,000 Dimensions

In 1957, Conrad Waddington visualized cellular fate determination as a marble rolling down a bifurcating, hilly landscape:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Biological Metaphor
- **Valleys**: Stable, differentiated cellular attractor basins (mature T-cells, B-cells, monocytes).
- **Ridges**: Epigenetic and transcriptional barriers preventing transdifferentiation between lineages.
- **Slopes**: Continuous developmental and maturation trajectories.
- **The Marbles**: Single cells progressing through developmental time.

[[0,1]]
### The Mathematical Reality
- Cells do not live in 3D: each cell is a point in a **$\sim 20,000$-dimensional gene expression space**.
- However, gene regulatory networks (GRNs) enforce tight co-expression and feedback loops.
- Cells cannot occupy arbitrary states in $\mathbb{R}^{20,000}$!
- They are constrained to lie on a **low-dimensional, curved Riemannian manifold** embedded within this high-dimensional space.
:::

> **Our Computational Mission**: Reconstruct and parameterize this low-dimensional manifold directly from single-cell transcriptomic matrices.

---

# Step 1 — Rebuilding the Quality-Controlled Matrix

Starting from the raw PBMC3k 10x Chromium matrix, we apply the filtering steps established in Lecture 5:

```python
import scanpy as sc, numpy as np

# Load raw matrix and ensure unique gene symbols
adata = sc.read_10x_mtx('raw_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
adata.var_names_make_unique()

# Annotate mtDNA-encoded transcripts and calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 1. Screen empty background droplets (> 500 UMIs)
adata = adata[adata.obs.total_counts > 500, :].copy()

# 2. Dynamic MAD filter for dying cells (Median + 3 MAD)
pct_mt = adata.obs['pct_counts_mt'].values
mt_cutoff = np.median(pct_mt) + 3 * (np.median(np.abs(pct_mt - np.median(pct_mt))) * 1.4826)

# 3. Filter dying cells and potential multi-cell doublets
valid = (adata.obs['pct_counts_mt'] < mt_cutoff) & (adata.obs['n_genes_by_counts'] <= 2500)
adata = adata[valid, :].copy()

# 4. Remove unexpressed / non-informative genes
sc.pp.filter_genes(adata, min_cells=3)
```

The resulting clean dataset contains **2,608 bona fide cells** across **13,638 expressed genes**.

---

# Step 2 — Normalization and Variance Stabilization

Raw integer UMI counts confound biological expression levels with technical capture efficiency:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### 1. Library-Size Scaling (CP10k)
- Cells vary in total sequencing depth from 500 to 15,000 UMIs.
- We rescale each cell's total count to a common target ($10,000$ counts per cell):
  $$x_{ij}' = \frac{x_{ij}}{\sum_k x_{ik}} \times 10,000$$
- `sc.pp.normalize_total(adata, target_sum=1e4)`

[[0,1]]
### 2. Log-Transformation ($x \mapsto \log(1+x)$)
- Count data exhibits a heavy right tail with strong mean–variance dependence.
- Log-transformation stabilizes variance across expression tiers:
  $$y_{ij} = \ln(1 + x_{ij}')$$
- `sc.pp.log1p(adata)`
- Keeps zeros at zero ($\ln(1) = 0$).
:::

```python
# Preserve raw counts in a layer for count-based models
adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

---

# ⚠ Pitfall — Freezing `adata.raw` Before Subsetting

Before filtering genes or scaling, we must freeze the full log-normalized matrix:

```python
# Freeze the complete log-normalized matrix across all 13,638 genes!
adata.raw = adata
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Why Freezing `.raw` Is Essential
- In subsequent steps, we will discard $\sim 85\%$ of all genes to perform PCA on informative genes only.
- However, critical lineage markers (e.g. transcription factors or low-abundance cytokines) may **not** be in that subset!
- If you subset `adata` without saving `.raw`, those genes are lost forever.
- Scanpy plotting functions (`sc.pl.umap`, `sc.pl.dotplot`, `sc.pl.violin`) automatically fall back to `adata.raw` to fetch and display any gene.

[[0,1]]
### The AnnData Layer Architecture

```text
adata.layers['counts']  <-- Raw integer UMI counts
adata.raw.X            <-- Normalized log1p (ALL genes)
adata.X                <-- Scaled z-scores (HVGs ONLY)
```

- `.layers['counts']`: needed for negative binomial / Poisson differential expression.
- `.raw`: needed for marker visualization.
- `.X`: scaled, low-dimension representation used for PCA and graph math.
:::

---
# Step 3 — Highly Variable Gene (HVG) Selection

Most of the 13,638 genes are universally expressed housekeeping components (e.g. actin, basic polymerases). They add noise rather than identity signal!

![HVG Selection](lecture_06_cellular_state_spaces/images/hvg_selection.png){width=78% center}

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"HVGs selected: {adata.var['highly_variable'].sum():,} / {adata.n_vars:,}")
```

---
# The Mean–Dispersion Relationship

Why can't we simply pick the genes with the highest raw sample variance?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Poisson Baseline
- In sequencing count data, the variance naturally increases with the mean expression $\mu$:
  $$\text{Var}(X) \ge \mathbb{E}[X]$$
- Highly expressed genes (like ribosomal proteins) automatically have enormous raw variances, even when completely uniformly expressed across all cells!
- Selecting by raw variance simply selects the most abundant housekeeping genes.

[[0,1]]
### Scanpy's Dispersion Binning
1. Bins all genes by their mean $\log$-expression into 20 equal-frequency bins.
2. Inside each bin, computes the mean dispersion:
   $$d = \frac{\text{Var}}{\text{Mean}}$$
3. Standardizes dispersions within each bin to z-scores:
   $$z_d = \frac{d - \bar{d}\_{\text{bin}}}{\sigma_{d, \text{bin}}}$$
4. Genes with $z_d > 0.5$ carry excess variance **relative to their expression level** $\to$ **Highly Variable Genes**.
:::

On PBMC3k, this isolates **1,838 informative genes** carrying true biological population structure.

---

# Step 4 — Gene Scaling (Z-Scoring)

Before spectral decomposition, we center and standardize each highly variable gene:

```python
# 1. Subset AnnData matrix to Highly Variable Genes only
adata = adata[:, adata.var.highly_variable].copy()

# 2. Scale each gene to unit variance (mean = 0, std = 1)
sc.pp.scale(adata, max_value=10)
```

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Why Z-Score Genes?
$$z_{ij} = \frac{y_{ij} - \mu_j}{\sigma_j}$$
- Without scaling, a gene with an expression range of $0 \to 10$ would exert **$100\times$ more leverage** on PCA than a critical transcription factor with a range of $0 \to 1$.
- Scaling gives every selected gene an equal opportunity to contribute to the principal axes.

[[0,1]]
### Why Clip at `max_value=10`?
- In single-cell data, a single cell with an extreme burst or technical spike can produce a z-score of $+40$.
- Squaring in PCA covariance would cause that single outlier cell to hijack an entire principal component!
- Clipping at $10~\sigma$ bounds extreme outliers while preserving genuine biological variance.
:::

---

# Step 5 — Spectral Decomposition (PCA)

Let $X \in \mathbb{R}^{n \times p}$ be the scaled expression matrix ($n$ cells, $p$ HVGs).
PCA finds orthogonal linear combinations of genes that capture maximal variance:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Covariance & Eigenvectors
The sample gene–gene covariance matrix is:
$$C = \frac{1}{n-1} X^T X \in \mathbb{R}^{p \times p}$$
Its spectral eigendecomposition is:
$$C = V \Lambda V^T$$
- $V \in \mathbb{R}^{p \times p}$: **Loadings** (eigenvectors). Each column is a "meta-gene" defining a linear axis of co-regulated genes.
- $\Lambda = \text{diag}(\lambda_1, \dots, \lambda_p)$: **Eigenvalues** representing variance explained along each axis.

[[0,1]]
### Singular Value Decomposition (SVD)
In practice, computing $X^T X$ directly is inefficient. We perform truncated SVD directly on $X$:
$$X = U S V^T$$
- $U \in \mathbb{R}^{n \times k}$: Left singular vectors.
- $S = \text{diag}(s_1, \dots, s_k)$: Singular values ($s_i = \sqrt{(n-1)\lambda_i}$).
- **Cell Scores (Coordinates)**:
  $$Z = X V = U S \in \mathbb{R}^{n \times k}$$
  These are the new low-dimensional coordinates for each cell!
- **Sign (reflection) invariance**: $(-U)S(-V)^T = USV^T$ — any PC axis $v$ and $-v$ span the exact same subspace, so different SVD solvers can flip a component's sign.
:::

```python
sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
```

---

# PCA as a Biological Denoiser

Why does projecting cells into 40 or 50 principal components clean the data?

```text
High-Dimensional Space (1,838 Genes)
┌─────────────────────────────────────────────────────────────────────────────┐
│ Biological Signal: Gene modules co-vary across cells (r > 0.6)              │
│ ──► Coherent variance aligns along dominant eigenvectors (PC1 - PC40)       │
├─────────────────────────────────────────────────────────────────────────────┤
│ Technical Noise: Poisson dropouts & capture errors are independent (r ≈ 0)  │
│ ──► Uncorrelated noise spreads thinly across all 1,800 dimensions           │
└─────────────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼ Truncate at PC 40
Low-Dimensional Denoised Space (40 Dimensions)
┌─────────────────────────────────────────────────────────────────────────────┐
│ Concentrates > 80% of true biological co-variation                          │
│ Discards the unstructured technical noise tail                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

> Truncating PCA retains concerted biological programs (e.g. cell cycle, immune activation, lineage identity) while filtering out uncorrelated single-cell noise.

---
# The Scree Plot (Elbow Method)


:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
![PCA Elbow Plot](lecture_06_cellular_state_spaces/images/pca_variance_elbow.png){width=100 center}

[[0,1]]
How much variance does each principal component capture?

```python
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```

### The Elbow Curve
- PC1–PC10 exhibit steep eigenvalues: major lineage bifurcations (T-cells vs. myeloid cells vs. B-cells).
- PC10–PC30 capture subtype variation (CD4 vs. CD8 T-cells, classical vs. non-classical monocytes).

### The Asymptotic Plateau
- Beyond PC35–PC40, the variance per component flattens into an asymptotic noise floor.
- Selecting $n_{\text{pcs}} = 40$ safely captures the biological manifold while preventing noise accumulation.
:::

---

# Where PCA Information Lives in `AnnData`

Running `sc.pp.pca(adata, n_comps=50)` populates three distinct compartments:

```text
AnnData
├─ .obsm['X_pca']               <-- Array of shape (n_cells, 50)
│                                    Cell coordinates (Scores: U * S)
│
├─ .varm['PCs']                 <-- Array of shape (n_genes, 50)
│                                    Gene weights per PC (Loadings: V)
│
└─ .uns['pca']                  <-- Dictionary of run metadata:
   ├─ ['variance']              <-- Eigenvalues: s_i^2 / (n - 1)
   └─ ['variance_ratio']        <-- Fraction of total variance explained
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Cell Coordinates (`.obsm['X_pca']`)
```python
coords = adata.obsm['X_pca']  # (2608, 50)
# PC1 coordinate of first cell:
print(coords[0, 0])
```

[[0,1]]
### Top Gene Loadings (`.varm['PCs']`)
```python
# Inspect PC1–PC2 colored by markers, and top gene loadings:
sc.pl.pca(adata, color=['CST3', 'NKG7', 'MS4A1', 'PPBP'], ncols=4)
sc.pl.pca_loadings(adata, components=[1, 2, 3])
```
:::

---
# Step 6 — The $k$-Nearest Neighbor ($k$-NN) Graph

Why can't we cluster cells using global Euclidean distance in gene space?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Curse of Dimensionality
- In high dimensions ($p \ge 100$), all pairwise Euclidean distances become nearly equidistant:
  $$\lim_{p \to \infty} \frac{\text{dist}\_{\max} - \text{dist}\_{\min}}{\text{dist}_{\min}} \to 0$$
- Global straight-line Euclidean distance between two cells across a curved manifold "cuts through" empty, unbiological space.

[[0,1]]
### The Manifold Solution
- A manifold is **locally Euclidean**: nearby points on the manifold reflect true biological proximity.
- By connecting each cell only to its **$k$ nearest neighbors** in 40D PCA space, we construct a discrete graph approximation of the continuous manifold!
- Geodesic distance (shortest path along the graph) respects the biological manifold.
:::

```python
# Construct k-NN graph with k=10 neighbors using 40 PCs
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
```

---

# Inside the $k$-NN Adjacency Matrices

Scanpy stores the computed graph in two sparse $n \times n$ matrices inside `adata.obsp`:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### 1. `adata.obsp['distances']`
- **Directed sparse distance matrix** ($n \times n$).
- Row `i` stores Euclidean distances to the **$k - 1$** nearest *other* cells (here $k = 10 \implies$ **9** entries; self is excluded); $0$ elsewhere.
- Access cell `i`'s neighbors with `adata.obsp['distances'][i].indices`.

[[0,1]]
### 2. `adata.obsp['connectivities']`
- **Symmetric sparse affinity matrix** ($n \times n$).
- Converts local distances into directed probabilities $w_{j|i} = \exp\left(-\frac{\max(0, d_{ij} - \rho_i)}{\sigma_i}\right)$ ($\rho_i$ = 1st-neighbor distance; $\sum_j w_{j|i} = \log_2 k$), then symmetrizes via fuzzy union $w_{ij} = w_{i|j} + w_{j|i} - w_{i|j}w_{j|i}$.
:::

> `adata.obsp['connectivities']` serves as the direct input for both **Leiden clustering** and **UMAP projection**!

---

# Step 7 — Leiden Community Detection

How do we partition the $k$-NN graph into distinct, biologically meaningful cell types?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Modularity Maximization
Community detection seeks a partition that maximizes the density of edges within clusters compared to a null model:

$$Q = \frac{1}{2m} \sum_{i,j} \left( A_{ij} - \gamma \frac{k_i k_j}{2m} \right) \delta(c_i, c_j)$$

- $A_{ij}$: graph adjacency weight between cell $i$ and $j$.
- $k_i = \sum_j A_{ij}$: degree of cell $i$.
- $m = \frac{1}{2} \sum_{i,j} A_{ij}$: total graph edge weight.
- $\gamma$: **Resolution parameter**.
- $\delta(c_i, c_j) = 1$ if cells $i$ and $j$ share the same cluster, $0$ otherwise.

[[0,1]]
### Why Leiden Replaced Louvain
- **The Louvain Flaw**: Louvain's heuristic can greedily merge communities, often producing internally disconnected clusters (sub-graphs with zero internal paths!).
- **The Leiden Guarantee**: Traag et al. (2019) proved that Leiden's refinement step guarantees all discovered communities are **connected** and asymptotically stable.
:::

```python
sc.tl.leiden(adata, resolution=0.9, random_state=0,
             flavor='igraph', n_iterations=2, directed=False)
print(f"Discovered {adata.obs['leiden'].nunique()} distinct clusters.")
```

---

# The Resolution Parameter ($\gamma$)

The parameter `resolution` controls the granularity of the community detection:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
```text
Low Resolution (gamma = 0.2)
┌─────────────────────────────────┐
│ [   Lymphoid   ]   [  Myeloid  ]│
│   (Coarse, broad lineages)      │
└─────────────────────────────────┘
                ▲
                │ Increasing gamma
                ▼
High Resolution (gamma = 1.5)
┌─────────────────────────────────┐
│ [CD4 T] [CD8 T] [NK] [B] [Mono] │
│   (Fine-grained subtypes)       │
└─────────────────────────────────┘
```

[[0,1]]
### How Resolution Functions
- Higher $\gamma$ increases the penalty for large clusters, forcing the algorithm to split communities into smaller, denser sub-clusters.
- Lower $\gamma$ merges related sub-types into broad parental lineages.
- **There is no single "correct" resolution**:
  - Broad tissue overview: $\gamma \approx 0.2 - 0.5$
  - Immune subtype discovery: $\gamma \approx 0.8 - 1.2$
:::

> In a real analysis, resolution is guided by whether new clusters exhibit biologically validated, distinct marker genes.

---

# Step 8 — Non-Linear Embedding: UMAP

To visualize the high-dimensional graph on a 2D computer screen, we use **UMAP** (Uniform Manifold Approximation and Projection):

![UMAP Clusters](lecture_06_cellular_state_spaces/images/umap_clusters.png){width=85% center}

```python
sc.tl.umap(adata, random_state=0)
sc.pl.umap(adata, color=['leiden', 'cell_type'])
```

---

# How UMAP Operates

UMAP models the high-dimensional manifold as a fuzzy simplicial set and optimizes a low-dimensional layout to match it:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### 1. High-Dimensional Affinities ($p_{ij}$)
The probability that cell $i$ and cell $j$ share an edge on the high-dimensional Riemannian manifold:
$$p_{j|i} = \exp\left(-\frac{\max(0, d_{ij} - \rho_i)}{\sigma_i}\right)$$
$$p_{ij} = p_{i|j} + p_{j|i} - p_{i|j} p_{j|i}$$

[[0,1]]
### 2. Low-Dimensional Affinities ($q_{ij}$)
In 2D space ($y_i \in \mathbb{R}^2$), edge probabilities follow a Student-t-like heavy-tailed kernel:
$$q_{ij} = \left(1 + a \lVert y_i - y_j\rVert^{2b}\right)^{-1}$$
($a \approx 1.93, b \approx 0.79$ under default parameters).
:::

### 3. Cross-Entropy Loss Optimization
UMAP minimizes the fuzzy set cross-entropy using stochastic gradient descent (SGD):
$$\mathcal{L} = \sum_{i \ne j} \left[ p_{ij} \log\left(\frac{p_{ij}}{q_{ij}}\right) + (1 - p_{ij}) \log\left(\frac{1 - p_{ij}}{1 - q_{ij}}\right) \right]$$
- First term pulls true neighbors together (attraction); second pushes non-neighbors apart (repulsion).
- Because SGD is stochastic, `random_state` in `sc.tl.umap` alters the 2D layout, while Leiden clusters come strictly from `adata.obsp['connectivities']` and never depend on UMAP.

---

# ⚠ Pitfall — Interpreting Distances in UMAP

UMAP is a non-linear projection designed to represent local neighborhoods, **not** global Euclidean space!

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### What UMAP Faithfully Preserves
- **Local Connectivity**: If two cells are close in UMAP, they are transcriptional neighbors in 40D PCA space.
- **Cluster Compactness**: Tight clusters represent dense regions of the graph.
- **Continuous Trajectories**: Cells undergoing a continuous developmental transition form unbroken strands.

[[0,1]]
### What UMAP Distorts
- **Inter-Cluster Distances**: The distance between two separate clusters in 2D is largely arbitrary. An empty space of 5 cm on screen does not mean they are $5\times$ more distant than two clusters separated by 1 cm!
- **Cluster Densities**: UMAP expands dense clusters and compresses diffuse clusters to maintain visual readability.
:::

> **Fundamental Rule**: Never perform clustering or quantitative distance measurements on 2D UMAP coordinates! Always cluster on the $k$-NN graph in PCA space.

---

# Step 9 — Marker Genes & Lineage Determination

Clusters are merely arbitrary integers ($0, 1, 2, \dots$) until we assign biological identity:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Method A: Data-Driven Wilcoxon DE
Tests every gene in cluster $C$ against all remaining cells:
```python
sc.tl.rank_genes_groups(
    adata, 
    groupby='leiden', 
    method='wilcoxon'
)
```
- Non-parametric: ranks expression, robust against outliers.
- Extracts the top enriched transcripts per cluster.

[[0,1]]
### Method B: Curated Markers & `score_genes`
Score each cell against canonical PBMC sets (using `use_raw=True`), then label each cluster by its highest mean score:
- **CD4 / CD8 T**: `IL7R, CCR7, LDHB` · `CD8A, CD8B, GZMK`
- **NK / B**: `GNLY, NKG7, KLRD1` · `MS4A1, CD79A, CD79B`
- **CD14 / FCGR3A Mono**: `CD14, LYZ, S100A8, S100A9` · `FCGR3A, MS4A7, LST1`
- **Dendritic / Megakaryocyte**: `FCER1A, CST3, CLEC10A` · `PPBP, PF4`
```python
sc.tl.score_genes(adata, gene_list=genes, score_name=f'score_{ct}', use_raw=True)
```
:::

---
# Visualizing Expression: Dot Plots

Dot plots display two dimensions of expression simultaneously across lineages:

![Marker Dotplot](lecture_06_cellular_state_spaces/images/marker_dotplot.png){width=60 center}

```python
sc.pl.dotplot(adata, marker_genes, groupby='cell_type', standard_scale='var')
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
- **Dot Color Intensity**: Mean expression level of the gene in that cell type.

[[0,1]]
- **Dot Diameter**: Fraction of cells in that cluster expressing the gene ($> 0$ counts).
:::

---
# Step 10 — Coarse-Grained Topology: PAGA

:::matrix{cols=50/50}
[[0,0]]
![PAGA Graph](lecture_06_cellular_state_spaces/images/paga_graph.png){width=100 center}
[[0,1]]
Can we quantify true biological connections between lineages without UMAP's visual distortion?

```python
sc.tl.paga(adata, groups='cell_type')
sc.pl.paga(adata, threshold=0.05, node_size_scale=2)
```

### How PAGA Operates
- **Partition-based Graph Abstraction** (Wolf et al., 2019).
- Compares the observed number of edges between cluster $I$ and cluster $J$ against the number expected under a random graph configuration:
  $$e_{IJ} = \frac{\text{observed edges}(I, J)}{\mathbb{E}[\text{random edges}(I, J)]}$$

### Interpreting the Graph
- **Thick Edges**: Genuine biological connectivity (e.g. CD14 Monocytes $\leftrightarrow$ FCGR3A Monocytes; CD4 T $\leftrightarrow$ CD8 T).
- **No Edges**: Lineages that diverged early in hematopoiesis (e.g. B-cells and Monocytes share zero graph edges).
:::

---
# Step 11 — Continuous State Transitions (Diffusion Pseudotime)

Not all biological states are discrete clusters. In peripheral blood, monocytes undergo a continuous differentiation process:

$$\text{Classical Monocytes (CD14}^{++}\text{ CD16}^-\text{)} \quad\longrightarrow\quad \text{Non-Classical Monocytes (CD14}^+\text{ CD16}^{++}\text{)}$$

:::matrix { cols="30/70" gap="30px" }
[[0,0]]
### Diffusion Pseudotime (DPT)
- Models a random walk on the $k$-NN transition matrix $T = D^{-1} A$.
- Measures geodesic diffusion distance from a chosen root cell ($i_{\text{root}}$):
  $$d_{\text{dpt}}(x, y) = \sum_{t=1}^\infty \lVert T^t(x, \cdot) - T^t(y, \cdot)\rVert$$
[[0,1]]
![Monocyte DPT](lecture_06_cellular_state_spaces/images/monocyte_dpt_trajectory.png){width=100 center}
:::

### Scanpy DPT Workflow & Dynamics
```python
sc.pp.neighbors(mono, n_neighbors=15, n_pcs=20)
sc.tl.diffmap(mono)                          # required before dpt
mono.uns['iroot'] = int(root_cell_idx)       # most classical CD14+ cell
sc.tl.dpt(mono)                              # writes obs['dpt_pseudotime']
```
- Along pseudotime (`sc.get.obs_df(..., use_raw=True)`), `CD14` & `LYZ` decline while `FCGR3A` (CD16) & `MS4A7` rise monotonically; replace disconnected `np.inf` values with `np.nan`.

---

# Data-Driven "Waddington Landscape"

We can close the loop with Waddington's metaphor by turning cell density in the state space into a **pseudo-potential**:

$$U(x, y) = -\log \hat{p}(x, y)$$

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Kernel Density Estimation (KDE)
- We compute a 2D Gaussian kernel density estimate $\hat{p}(x, y)$ across the state space coordinates.
- Highly populated regions produce large $\hat{p}(x, y)$.
- Taking $-\log \hat{p}(x, y)$ converts high-density regions into **gravitational valleys (attractor basins)**!
- Sparsely populated transition zones become **elevated energy barriers (ridges)**.

[[0,1]]
### The Biological Interpretation
- **Valleys**: Stable cell types (e.g. CD4 T-cells, B-cells) where transcription factors hold the cell in equilibrium.
- **Valley Depth**: Relative population abundance under homeostatic conditions.
- **Passes / Saddles**: Preferred developmental routes where differentiation occurs with minimal regulatory resistance.
:::

```python
from scipy.stats import gaussian_kde
import plotly.graph_objects as go
coords = adata.obsm['X_umap']
kde = gaussian_kde(coords.T, bw_method=0.12)      # evaluate on 2D meshgrid → go.Surface
potential = -np.log(kde(coords.T) + 1e-4)         # +1e-4 floor prevents log(0)
```

---

# Homework Task 1 — Spectral PCA from Scratch

In Homework Task 1, you write a pure NumPy function `pca_svd(X, n_comps)` using SVD:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Algorithmic Requirements
1. **Column Centering**:
   $$X_{\text{centered}} = X - \mu_X$$
   Compute the column mean using `np.mean(X, axis=0)` and subtract.
2. **Economy SVD**:
   Use `np.linalg.svd(X_centered, full_matrices=False)`.
   Extract singular values $S$ and left singular vectors $U$.
3. **Cell Scores**:
   $$Z = U_{[:, :k]} \cdot S_{[:k]}$$
4. **Explained Variance Ratio** (divide by **all** singular values *before* slicing `[:k]`):
   $$\text{var}_i = \frac{S_i^2}{n - 1}, \qquad \text{ratio}_i = \frac{S_i^2}{\sum_{\text{all } j} S_j^2}$$

[[0,1]]
### Validation Checklist
Compare your top 5 PCs against `adata.obsm['X_pca']` and `adata.uns['pca']['variance_ratio']`:
```python
corr = np.abs(np.corrcoef(my_pc, scanpy_pc)[0, 1])
```
- **Dense input**: ensure `X = np.asarray(X)` (or `.toarray()` if sparse) before centering.
- **Total variance pitfall**: normalizing by `S[:k]` instead of full `S` makes ratios sum to 1.0 and breaks the match!
- **Absolute correlation**: reflect on why `np.abs(...)` is needed (see SVD geometry in Step 5).
:::

---

# Homework Task 2 — How Many PCs?

In Homework Task 2, you analyze the cumulative explained variance stored in `adata.uns['pca']`:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Cumulative Variance Math
- `adata.uns['pca']['variance_ratio']` contains the variance share for the top 50 components.
- Use `np.cumsum()` to compute the cumulative sum curve.
- Normalize relative to the total variance of the 50 PCs:
  $$C_k = \frac{\sum_{i=1}^k \text{ratio}_i}{\sum_{i=1}^{50} \text{ratio}_i}$$
- Find the first component where $C_k \ge 0.80$ (watch 0-based vs. 1-based indexing: `np.argmax(C >= 0.80) + 1`).

[[0,1]]
### Methodological Comparison
- Contrast your mathematical 80% cutoff against the empirical $n_{\text{pcs}} = 40$ used in the notebook:
  - What happens if you pick too few PCs (e.g. 5 PCs)? You lose fine subtype distinction (e.g. dendritic cells collapse into monocytes).
  - What happens if you pick too many PCs (e.g. 200 PCs)? You re-introduce uninformative Poisson technical noise into distance calculations.
:::

---

# Homework Task 3 — Building the $k$-NN Graph Yourself

In Homework Task 3, you reconstruct the neighborhood graph using `sklearn.neighbors.NearestNeighbors`:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Self-Neighbor Pitfall ⚠
When you run:
```python
from sklearn.neighbors import NearestNeighbors
nn = NearestNeighbors(n_neighbors=10).fit(X_pca)
distances, indices = nn.kneighbors(X_pca)
```
- In Scikit-Learn, query cell $i$ is evaluated against the training set, which includes cell $i$!
- **Column 0 is always cell $i$ itself** (at distance 0.0).
- Columns 1 through 9 are the 9 *other* closest cells.

[[0,1]]
### Scanpy's Storage Format
- Scanpy stores `adata.obsp['distances']` with **$k - 1 = 9$ distinct neighbors** per row (self excluded).
- Neighbor indices for cell $i$ are accessed via:
  ```python
  scanpy_neighbors = set(adata.obsp['distances'][i].indices)
  ```
- Slice off column 0 (`indices[i, 1:]`) to get the 9 other neighbors and compute overlap:
  $$\text{overlap}_i = \frac{\lvert\text{ours}_i \cap \text{scanpy}_i\rvert}{\lvert\text{scanpy}_i\rvert}$$
  *(Exact vs. approximate search can differ slightly from 1.0 on tie boundaries.)*
:::

---

# Homework Task 4 — Resolution Sweep

In Homework Task 4, you evaluate how the Leiden resolution parameter $\gamma$ governs cluster granularity:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Sweeping the Parameter
Run Leiden clustering across the designated range:
$$\gamma \in [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]$$

```python
resolutions = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
# For each r:
# 1. sc.tl.leiden(adata, resolution=r, key_added=f'leiden_r{r}',
#                 flavor='igraph', n_iterations=2, directed=False, random_state=0)
# 2. Count clusters: adata.obs[f'leiden_r{r}'].nunique()
# 3. Side-by-side UMAPs: sc.pl.umap(adata, color=['leiden_r0.1', 'leiden_r2.0'])
```

[[0,1]]
### What to Observe
- Plot the number of discovered clusters against $\gamma$.
- Observe how the curve scales monotonically:
  - At $\gamma = 0.1$, only coarse mega-clusters form.
  - At $\gamma = 2.0$, clusters split into micro-subsets.
- Which resolutions preserve major immune lineages without over-fragmenting homogeneous cell types?
:::

---

# Homework Task 5 — Identifying the Mystery Cluster

In Homework Task 5, you investigate the smallest, rarest cluster in the PBMC3k dataset:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Locating Rare Populations
1. **Find Smallest Cluster**:
   Use `adata.obs['leiden'].value_counts()` (`.idxmin()` returns its string ID).
2. **Query Marker Table**:
   Extract differential expression statistics:
   ```python
   df = sc.get.rank_genes_groups_df(adata, group=smallest_id)
   top_markers = df['names'].head(10).tolist()
   ```
3. **Cell-Type Proportions**:
   `adata.obs['cell_type'].value_counts(normalize=True) * 100` gives percentages for a bar chart.

[[0,1]]
### Biological Deductions
- Inspect the top 10 ranked markers.
- Are these genes associated with:
  - Antigen-presenting dendritic cells (`FCER1A`, `CLEC10A`)?
  - Platelet / megakaryocyte fragments (`PPBP`, `PF4`)?
  - Plasma cells (`MZB1`, `IGKC`)?
- Check how this cluster's automatic annotation aligns with cell type frequency distributions.
:::

---

# Homework Task 6 — Conceptual Review & Guidance

Task 6 presents 6 core theoretical questions examining the mathematical and biological foundations:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### 1. Why PCA Space for $k$-NN?
- Contrast the biological signal concentration in the first 40 PCs with the curse of dimensionality and unstructured technical noise across 1,800 HVGs.

### 2. Clustering on $k$-NN vs. UMAP?
- Why is the $k$-NN graph an un-distorted topological representation, while 2D UMAP intentionally warps global distances to fit a planar screen?

### 3. Skipping Gene Scaling (`sc.pp.scale`)?
- What happens to the covariance matrix if genes with raw mean expression of 100 are not standardized against genes with mean expression of 0.1?

[[0,1]]
### 4. Selecting Resolution $\gamma$?
- How does a researcher balance under-clustering (missing rare populations) vs. over-clustering (splitting homogeneous noise)? What role do marker genes play?

### 5. Waddington Landscape Realities?
- Contrast continuous peripheral monocyte maturation ($CD14 \to FCGR3A$) with early hematopoietic divergence (B-cells vs. T-cells).

### 6. UMAP Random Seeds vs. Leiden?
- Why does altering `random_state` in `sc.tl.umap` rotate and distort the 2D visual layout while leaving Leiden cluster assignments identical?
:::

---

# Putting It All Together: Complete Manifold Pipeline

```python
import numpy as np, scanpy as sc
# 1. Rebuild clean dataset (Lecture 5 QC) and preserve raw counts + .raw
adata = sc.read_10x_mtx('raw_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.total_counts > 500, :].copy()
pct_mt = adata.obs['pct_counts_mt'].values
mt_thr = np.median(pct_mt) + 3 * 1.4826 * np.median(np.abs(pct_mt - np.median(pct_mt)))
adata = adata[(adata.obs['pct_counts_mt'] < mt_thr) & (adata.obs['n_genes_by_counts'] <= 2500), :].copy()
sc.pp.filter_genes(adata, min_cells=3)
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # Crucial: freeze full transcriptome for marker lookups
# 2. Select HVGs, scale variance, and compute spectral PCA
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
# 3. Construct k-NN graph, partition via Leiden, and project onto UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=0.9, random_state=0, flavor='igraph', n_iterations=2, directed=False)
sc.tl.umap(adata, random_state=0)
# 4. Identify markers, score lineages, and compute PAGA topology
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.tl.paga(adata, groups='leiden')
```

---

# Summary

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### Preprocessing & PCA
#### Linear Reduction
- Rescale counts to CP10k
- $\log(1 + x)$ variance stabilization
- Freeze `adata.raw` before slicing
- Dispersion binning selects HVGs
- SVD: $X = U S V^T$
- PCA denoises uncorrelated noise

[[0,1]]{ .green }
### Graphs & Clustering
#### Manifold Learning
- $k$-NN in 40D PCA space ($k=10$)
- Connectivities stored in `.obsp`
- Leiden optimizes modularity
- Guarantees connected clusters
- Resolution $\gamma$ controls granularity
- UMAP embeds graph in 2D

[[0,2]]{ .blue }
### Biology & Lineages
#### Cell Identity
- Wilcoxon DE identifies markers
- Dot plots show level & fraction
- PAGA tests lineage connectivity
- DPT tracks monocyte maturation
- Potential $U = -\log \hat{p}(x, y)$
- Attractor basins = stable fates
:::

---

# Next Steps — Lecture 7

In **Lecture 7: Folding the Molecular Machine (Protein 3D Conformation)**, we transition from cell-level transcriptomes to atomic-resolution macromolecular geometry:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### What We Accomplished in Lectures 5 & 6
- Decoded 737,280 raw microfluidic droplets into clean, individual cells.
- Denoised transcriptional variance using spectral PCA decomposition.
- Reconstructed the high-dimensional cellular manifold via $k$-NN graphs.
- Partitioned immune communities with Leiden clustering and annotated lineages.
- Tracked continuous differentiation trajectories through diffusion pseudotime.

[[0,1]]
### What Comes Next (Lecture 7)
- **Peptide Geometry**: Planar resonance of peptide bonds, $\phi/\psi$ dihedral angles, and Ramachandran steric constraints.
- **Experimental Coordinates**: Parsing Protein Data Bank (`.pdb` / `.cif`) structural files with `Bio.PDB`.
- **The Kabsch Algorithm**: SVD-based optimal 3D spatial rotation and translation.
- **Structural Alignment**: Root-Mean-Square Deviation (RMSD) superposition across mutant and wild-type viral proteases.
- **In-Notebook WebGL**: 3D interactive rendering using `py3Dmol`.
:::
