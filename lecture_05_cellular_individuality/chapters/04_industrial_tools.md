:::titlepage
[[top]]
Part 4

[[title]]
Industrial Tools

[[subcaption]]
Scanpy · anndata · scipy.stats · 10x Genomics PBMC3k
:::

---

# Curriculum

- **10x Genomics PBMC3k** — raw MatrixMarket (`.mtx`), `genes.tsv`, `barcodes.tsv`, reading via `sc.read_10x_mtx`, fast `.h5ad` caching.
- **The AnnData Container** — `.X` (sparse CSC/CSR), `.obs` (droplets), `.var` (genes), memory economics (22 MB sparse vs. 90 GB dense).
- **Barcode Rank Analysis (Knee Plot)** — 737,280 droplets, the cellular plateau, the microfluidic cliff, ambient RNA soup.
- **Cell Calling Heuristics** — Cell Ranger v2 order-of-magnitude cutoff ($m/10$) vs. hard UMI threshold ($> 500$).
- **Quality Control Metrics** — `sc.pp.calculate_qc_metrics`, 13 human mitochondrial genes (`MT-`), library size, gene diversity.
- **Robust Outlier Detection** — Median Absolute Deviation (MAD), the normal scaling constant ($1.4826$), dynamic thresholding ($\text{Median} + 3~\text{MAD}$).
- **Doublet Contamination** — Poisson microfluidic encapsulation, high-gene outliers, simulation scoring with Scrublet.
- **AnnData Hygiene & Sparsity** — views vs. copies (`.copy()`), extreme sparsity, observed zeros vs. Poisson dropout $\exp(-\mu)$.
- **UMI Collision Mathematics** — occupancy statistics, barcode capacity for $L=10$ vs. $L=12$, collision inversion correction.

---

# Roadmap

The end-to-end single-cell quality control and filtering pipeline:

:::static-diagram { width=95% center }
{
  "nodes": [
    { "id": "raw",   "label": "10x Raw Files\n(.mtx, genes, barcodes)", "x": 40,   "y": 80,  "width": 150, "height": 60 },
    { "id": "read",  "label": "sc.read_10x_mtx\n737,280 droplets",      "x": 220,  "y": 80,  "width": 150, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "knee",  "label": "Knee Plot\nBarcode Rank Analysis",       "x": 400,  "y": 80,  "width": 160, "height": 60 },
    { "id": "call",  "label": "Cell Calling\ntotal_counts > 500",       "x": 590,  "y": 80,  "width": 150, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "qc",    "label": "calculate_qc_metrics\nMT- genes (13)",   "x": 770,  "y": 80,  "width": 160, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "mad",   "label": "MAD Scaling\n3 × MAD threshold",        "x": 960,  "y": 80,  "width": 150, "height": 60, "style": "fill: #FF9800; color: #fff;" },

    { "id": "mt",    "label": "Filter Dying Cells\n(Homework Task 1)",  "x": 960,  "y": 250, "width": 160, "height": 60, "style": "fill: #F44336; color: #fff;" },
    { "id": "dbt",   "label": "Filter Doublets\n(Homework Task 1)",     "x": 770,  "y": 250, "width": 150, "height": 60, "style": "fill: #F44336; color: #fff;" },
    { "id": "clean", "label": "adata_clean\nFiltered Real Cells",       "x": 580,  "y": 250, "width": 160, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "umi",   "label": "UMI Collision Math\n(Homework Task 4)",  "x": 390,  "y": 250, "width": 160, "height": 60, "style": "fill: #673AB7; color: #fff;" },
    { "id": "norm",  "label": "Normalization\nnormalize_total + log1p", "x": 190,  "y": 250, "width": 170, "height": 60, "style": "fill: #009688; color: #fff;" },
    { "id": "next",  "label": "Lecture 6 Preview\nPCA · UMAP · Leiden", "x": 10,   "y": 250, "width": 150, "height": 60 }
  ],
  "links": [
    { "from": "raw",   "to": "read" },
    { "from": "read",  "to": "knee" },
    { "from": "knee",  "to": "call" },
    { "from": "call",  "to": "qc" },
    { "from": "qc",    "to": "mad" },
    { "from": "mad",   "to": "mt" },
    { "from": "mt",    "to": "dbt" },
    { "from": "dbt",   "to": "clean" },
    { "from": "clean", "to": "umi" },
    { "from": "clean", "to": "norm" },
    { "from": "norm",  "to": "next" }
  ]
}
:::

> The practice notebook processes the raw PBMC3k count matrix, screens ambient noise, identifies apoptotic cells via robust statistics, and prepares the dataset for downstream manifold learning.

---

# The Single-Cell Toolbox

Single-cell analysis operates on massive, ultra-sparse expression matrices across thousands of individual cells:

:::matrix { .cards cols="1fr 1fr 1fr 1fr" gap="20px" margin-top: 20px; }
[[0,0]]{ .purple }
### Scanpy
The standard Python framework for single-cell transcriptomics: preprocessing, QC metrics, normalization, clustering, and trajectory inference.

[[0,1]]{ .green }
### AnnData
Annotated Data container: tracks the sparse expression matrix alongside synchronized cell metadata (`.obs`) and gene metadata (`.var`).

[[0,2]]{ .blue }
### scipy.sparse
Compressed Sparse Column (CSC) & Row (CSR) linear algebra: enables 4,000× RAM compression compared to dense matrices.

[[0,3]]{ .orange }
### scipy.stats
Robust statistical dispersion: Median Absolute Deviation (MAD), normal scale factors, and discrete occupancy distributions.
:::

```bash
# Installation in Google Colab / virtual environment
pip install -q scanpy anndata numpy pandas scipy matplotlib seaborn
```

> Numerical values and benchmarks in this lecture reflect a real run on Python 3.11 with `scanpy 1.11.5` and `anndata 0.12.19`.

---
# The Data — 10x Genomics PBMC3k

Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor sequenced on the 10x Chromium platform:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Why Start with RAW Data?
- Pre-filtered matrices include only what Cell Ranger *guessed* were cells.
- Real single-cell quality control starts at the **droplet level**:
  - Distinguishing real cells from empty droplets.
  - Measuring ambient cell-free mRNA contamination.
  - Detecting cells damaged during microfluidic loading.
- We download `pbmc3k_raw_gene_bc_matrices.tar.gz` (11.72 MB).

[[0,1]]{ .smaller-table }
### The Extracted File Tree (`hg19/`)

| File | Lines | Uncompressed Size | Role |
| :--- | :--- | :--- | :--- |
| `matrix.mtx` | 2,905,864 | **40.40 MB** | Sparse coordinate matrix |
| `barcodes.tsv` | 737,280 | **11.95 MB** | 16-nt droplet barcodes |
| `genes.tsv` | 32,738 | **798 KB** | Ensembl IDs + gene symbols |

In contrast, the pre-filtered set contains only **2,700** barcodes in `barcodes.tsv` (44.8 KB).
:::

```bash
wget "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_raw_gene_bc_matrices.tar.gz"
tar -xzf pbmc3k_raw_gene_bc_matrices.tar.gz
```
---

# Inside the Raw Files

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### `matrix.mtx` (MatrixMarket Coordinate)
```text
%%MatrixMarket matrix coordinate real general
%
32738 737280 2905861
1821 1 1
2129 1 1
```
- Line 3: `32738` genes × `737280` droplets with `2905861` non-zero entries.
- Line 4: Gene 1821 in Droplet 1 has **1** UMI.
- 1-based indexing on disk!

[[0,1]]
### `genes.tsv` & `barcodes.tsv`
```text
ENSG00000243485   MIR1302-10
ENSG00000237613   FAM138A
ENSG00000186092   OR4F5
```
Tab-separated: Ensembl ID (column 1) and HGNC symbol (column 2).

```text
AAACATACAAAACG-1
AAACATACAAAAGC-1
AAACATACAAACAG-1
```
16 nt bead barcode + `-1` (10x Chromium library index suffix).
:::

> Notice that `matrix.mtx` has genes as rows and cells as columns. Scanpy automatically transposes this to **Cells × Genes** upon ingestion!

---

# Tool 1 — Reading into AnnData

Scanpy provides an optimized reader for 10x Genomics folder structures:

```python
import scanpy as sc

# Read raw 10x directory into an AnnData object
adata = sc.read_10x_mtx(
    'raw_gene_bc_matrices/hg19/',  # path containing matrix.mtx, genes.tsv, barcodes.tsv
    var_names='gene_symbols',      # index genes by HGNC symbol instead of Ensembl ID
    cache=True                     # create an on-disk .h5ad cache for fast subsequent reads
)

# Ensure all gene symbols in var_names are unique
adata.var_names_make_unique()
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Load Times & Cache
- First read (`.mtx` text parsing): **2.79 s**
- Second read (reading cached `.h5ad`): **2.54 s**
- Generated cache file:
  `cache/...matrix.h5ad` (**33.21 MB**)

[[0,1]]
### Ingested Dimensions
```text
AnnData object with n_obs × n_vars = 737280 × 32738
    var: 'gene_ids'
```
- $n_{\text{obs}} = 737,280$ droplets (cells)
- $n_{\text{vars}} = 32,738$ annotated genes
:::

---

# Anatomy of an `AnnData` Object

An `AnnData` (Annotated Data) structure organizes high-dimensional biological measurements:

```text
               ┌─────────────────────── var ───────────────────────┐
               │    gene_ids       mt     n_cells_by_counts  ...   │
               └───────────────────────────────────────────────────┘
               ┌───────────────────────────────────────────────────┐ ▲
               │                                                   │ │
               │                                                   │ │
               │                       .X                          │ │
      obs      │             (Cells × Genes Matrix)                │ │
┌────────────┐ │            scipy.sparse.csc_matrix                │ │ n_obs
│ total_count│ │              737,280 × 32,738                     │ │ (droplets)
│ pct_mt     │ │                                                   │ │
│ n_genes    │ │                                                   │ │
└────────────┘ └───────────────────────────────────────────────────┘ ▼
               ◄──────────────────── n_vars ───────────────────────►
                                     (genes)
```

- **`.X`**: Primary count matrix (sparse `float32`).
- **`.obs`**: Observation metadata dataframe (rows = barcodes/droplets).
- **`.var`**: Variable metadata dataframe (rows = genes).
- **`.uns`**: Unstructured dictionary (colors, run parameters, global statistics).
- **`.obsm` / `.varm`**: Multidimensional embeddings (`X_pca`, `X_umap`).
- **`.obsp` / `.varp`**: Pairwise graphs (cell connectivities, gene distance matrices).

---

# Sparsity & RAM Economics

Why not just use a standard dense 2D NumPy array?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Reality of Single-Cell Data
- Total matrix entries:
  $$737,280 \times 32,738 = 24,136,873,440$$
- Non-zero entries (`nnz`): **2,905,861**
- Density:
  $$\frac{2,905,861}{24,136,873,440} \approx \mathbf{0.0120\%}$$
- **99.988%** of all entries are exact zeros!

[[0,1]]
### Memory Footprint in RAM

| Format | Representation | Memory Footprint |
| :--- | :--- | :--- |
| **Dense** | `float32` (4 bytes/val) | **89.92 GB** 💥 |
| **Dense** | `float64` (8 bytes/val) | **179.85 GB** 💥 |
| **CSC Sparse** | `data + indices + indptr` | **22.29 MB** ⚡ |

Dense allocation crashes standard laptops and Google Colab (12 GB RAM) with an instant **OOM Kill**!
:::

> `scipy.sparse` stores only non-zero coordinates and counts. AnnData seamlessly preserves sparsity during slicing, mathematical operations, and disk serialization.

---

# ⚠ Pitfall — Scipy Sparse Sums and 2D Matrices

When computing cell totals across columns (`adata.X.sum(axis=1)`):

```python
# Calculating total UMIs per droplet
s = adata.X.sum(axis=1)

print(type(s))    # <class 'numpy.matrix'>  <-- 2D MATRIX, NOT 1D NDARRAY!
print(s.shape)    # (737280, 1)
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Hidden Bug
- In SciPy sparse, `.sum()` historically returns a legacy **2D `np.matrix`**.
- Boolean filtering with a 2D matrix (`adata[s > 500]`) can cause silent index broadcasting bugs or shape mismatches.
- `np.median(s)` or `np.sort(s)` may produce unexpected 2D outputs.

[[0,1]]
### The Correct Idiom
Convert the 2D matrix into a flat 1D NumPy array:

```python
# Option A: using .A1 (NumPy matrix flattening shortcut)
umi_counts = s.A1

# Option B: explicit array flattening (recommended)
umi_counts = np.asarray(s).ravel()

# Option C: flatten()
umi_counts = np.array(s).flatten()
```
`type(umi_counts)` is now `<class 'numpy.ndarray'>` with shape `(737280,)`.
:::

---

# ⚠ Pitfall — Duplicate Gene Symbols

In human and mouse annotations, different Ensembl gene IDs often share the exact same HGNC symbol:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Problem in `genes.tsv`
- **91** gene symbols occur more than once in `genes.tsv` (95 redundant copies).
- Examples:
  - `U1`: `ENSG00000233421` and `ENSG00000228549`
  - `EFNA3`: `ENSG00000143590` and `ENSG00000251246`
  - `RPP14`: `ENSG00000163684` and `ENSG00000255154`
- If gene symbols serve as variable names (`adata.var_names`), duplicate names break pandas indexing and raise `KeyError` during slicing.

[[0,1]]
### The Scanpy Solution

```python
adata.var_names_make_unique()
```

Scanpy checks for collisions and appends `-1`, `-2`, etc. to duplicate symbols:
- `U1` $\to$ `U1-1`
- `EFNA3` $\to$ `EFNA3-1`
- `RPP14` $\to$ `RPP14-1`

> Always call `adata.var_names_make_unique()` immediately after reading 10x text outputs!
:::

---

# Tool 2 — The Barcode Rank ("Knee") Analysis

Where do the 7,071,612 sequencing reads in the raw library actually go?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### UMI Count Distribution
Out of **737,280** total bead barcodes:

| UMI Threshold | Barcode Count | Fraction of Barcodes |
| :--- | :--- | :--- |
| $\text{UMI} = 0$ | **578,390** | **78.45%** |
| $\text{UMI} \ge 1$ | **158,890** | **21.55%** |
| $\text{UMI} \ge 10$ | **14,397** | **1.95%** |
| $\text{UMI} \ge 100$ | **2,967** | **0.40%** |
| $\text{UMI} \ge 500$ | **2,708** | **0.37%** |
| $\text{UMI} \ge 1,000$ | **2,547** | **0.35%** |
| $\text{UMI} \ge 5,000$ | **69** | **0.01%** |

[[0,1]]
### Library Concentration

```text
Total Library: 7,071,612 UMIs
├─ Top 2,700 droplets (Cells):
│  └─ 6,390,631 UMIs (90.37%)
│
├─ 254 droplets (Transitional):
│  └─ 60,148 UMIs (0.85%)
│
└─ 155,928 droplets (Ambient Soup):
   └─ 616,646 UMIs (8.72%)
```

Over **90%** of all mRNA molecules are contained within fewer than **0.4%** of the droplet partitions!
:::

---
# Visualizing the Knee Plot

We rank all barcodes descending by total UMI count and plot them on a $\log_{10}$–$\log_{10}$ scale:

![Knee Plot](lecture_05_cellular_individuality/images/knee_plot.png){width=50% center}

```python
umi_counts = np.array(adata.X.sum(axis=1)).flatten()
sorted_umi = np.sort(umi_counts)[::-1]

plt.figure(figsize=(8, 4.5))
plt.loglog(sorted_umi, color='navy', lw=2)
plt.axhline(500, color='red', linestyle='--', label='UMI Cutoff = 500')
plt.xlabel('Barcode Rank (log scale)'); plt.ylabel('Total UMI Count (log scale)')
```

---

# Physics of the Three Droplet Regimes

Why does the barcode rank curve form this characteristic profile?

:::matrix { cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]
### 1. Cellular Plateau
**Ranks 1 to ~2,000**
- Intact single cells encapsulated with a hydrogel bead.
- Cell cytoplasm volume: $\sim 1\text{ pL}$.
- Delivers $1,500$ to $16,000$ distinct mRNA molecules to the droplet lysis reaction.

[[0,1]]
### 2. The Cliff (Knee)
**Ranks 2,000 to 2,700**
- Sharp transition across an order of magnitude:
  - Rank 2,000: **1,775 UMIs**
  - Rank 2,500: **1,102 UMIs**
  - Rank 2,700: **548 UMIs**
  - Rank 3,000: **89 UMIs**
- Smaller cells (platelets, low-RNA lymphocytes) and damaged cells.

[[0,2]]
### 3. Ambient Tail
**Ranks > 3,000**
- **Empty droplets**: oil bubbles containing only buffer and bead.
- Droplet volume: $\sim 1\text{ nL}$ ($1,000\times$ cell volume).
- Passively captures cell-free ambient RNA leaked from lysed cells into the loading suspension.
:::

---

# Cell Calling Algorithms

How do industrial pipelines separate real cells from ambient background?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Cell Ranger v2 Heuristic
Takes expected recovered cell count $N$ (e.g. 3,000):
1. Compute the 99th percentile of the top $N$ barcodes:
   $$m = \text{percentile}_{99}(\text{top } N) = 5,688.7$$
2. Set threshold to one-tenth of $m$:
   $$T = \frac{m}{10} = 568.9\text{ UMIs}$$
3. Call all droplets with $\text{UMI} \ge T$ as cells:
   $$\implies \mathbf{2,693}\text{ cells called}$$

[[0,1]]
### Notebook Hard Threshold (> 500)
- Simple cutoff: retain droplets with `total_counts > 500`.
- Retains **2,708** droplets out of 737,280.
- Comparing against Cell Ranger's official 2,700 filtered cells:
  - **2,700 of 2,700** (100%) Cell Ranger cells pass $> 500$.
  - Our cutoff rescues **8 extra droplets** with $506 \le \text{UMI} \le 548$ that are genuine cells.
:::

```python
# Initial cell calling: drop empty droplets (< 500 UMIs)
adata = adata[adata.obs.total_counts > 500, :].copy()
print(f"Cells remaining after basic threshold: {adata.n_obs:,}")  # 2,708
```

---

# Ambient RNA: What Is in the "Empty" Droplets?

Profiling empty droplets ($10 \le \text{UMI} \le 100$, $N = 11,435$) against real cells ($N = 2,700$):

:::matrix { cols="1.2fr 1fr" gap="30px" }
[[0,0]]{ .smaller-table }
### Top Genes in Ambient Background
| Gene | Description | % of Empty Pool | % of Cell Pool |
| :--- | :--- | :--- | :--- |
| **`MALAT1`** | Nuclear lncRNA | **3.61%** | 2.53% |
| **`B2M`** | $\beta_2$-microglobulin | **2.03%** | 1.90% |
| **`TMSB4X`** | Thymosin beta-4 | **1.84%** | 1.94% |
| **`MT-CO1`** | Cytochrome c oxidase I | **1.53%** | 0.52% |
| **`FTL`** | Ferritin light chain | **1.42%** | 1.17% |
| **`RPL13`** | Ribosomal protein L13 | **1.39%** | 1.21% |
| **`MT-CO3`** | Cytochrome c oxidase III | **1.11%** | 0.38% |

[[0,1]]
### Mitochondrial Enrichment

```text
Pooled Mitochondrial Reads (MT-):
├─ Empty Droplets:   6.62% of pool
└─ Real Intact Cells: 2.15% of pool
```

> **The Ambient Signature**: Empty droplets are enriched **$3\times$ higher** in mitochondrial RNA than intact cells! When dying cells lyse in the suspension tube before droplet encapsulation, their cytoplasm and mitochondria leak into the ambient soup.
:::

---

# Tool 3 — Quality Control (QC) Metrics

Scanpy calculates comprehensive QC metrics across cells and genes in parallel:

```python
# 1. Flag human mitochondrial genes (encoded in mtDNA, starting with 'MT-')
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# 2. Compute library sizes, detected genes, and mitochondrial percentages
sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars=['mt'], 
    percent_top=None, 
    log1p=False, 
    inplace=True
)
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Columns Added to `.obs` (Droplets)
- `total_counts`: sum of all UMIs in this droplet.
- `n_genes_by_counts`: number of genes with count $\ge 1$.
- `total_counts_mt`: sum of UMIs from `MT-` genes.
- `pct_counts_mt`: $\frac{\text{MT UMIs}}{\text{Total UMIs}} \times 100$.

[[0,1]]
### Columns Added to `.var` (Genes)
- `n_cells_by_counts`: number of cells expressing this gene.
- `total_counts`: sum of UMIs for this gene across all cells.
- `mean_counts`: average expression per cell.
- `pct_dropout_by_counts`: percentage of cells with count $= 0$.
:::

---

# The 13 Human Mitochondrial Genes

Human mitochondrial DNA (mtDNA) encodes exactly 13 respiratory chain proteins:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Complete Set in PBMC3k
- **Complex I (NADH dehydrogenase)**:
  `MT-ND1`, `MT-ND2`, `MT-ND3`, `MT-ND4`, `MT-ND4L`, `MT-ND5`, `MT-ND6`
- **Complex III (Cytochrome bc₁)**:
  `MT-CYB`
- **Complex IV (Cytochrome c oxidase)**:
  `MT-CO1`, `MT-CO2`, `MT-CO3`
- **Complex V (ATP synthase)**:
  `MT-ATP6`, `MT-ATP8`

[[0,1]]
### Why MT- and Not Nuclear Genes?
- Mitochondria also contain hundreds of nuclear-encoded genes (e.g. `MRPL30`, `NDUFA1`).
- Only transcripts encoded on the **circular mitochondrial chromosome** (`MT-`) reflect the organellar entrapment mechanism during apoptosis!
- In mouse datasets, the prefix is lowercase: `mt-` (e.g. `mt-Nd1`).
:::

```python
mt_genes = adata.var_names[adata.var['mt']].tolist()
print(f"Detected {len(mt_genes)} mitochondrial genes: {', '.join(mt_genes[:4])}...")
```

---
# Visualizing QC Distributions (Violin Plots)

Inspecting metric distributions across the called cells ($N = 2,708$):

![QC Violins](lecture_05_cellular_individuality/images/qc_violin.png){width=70% center}

```python
sc.pl.violin(
    adata, 
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4, multi_panel=True, size=1
)
```

:::matrix { cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]
**`n_genes_by_counts`**
- Median: **816 genes**
- 50% IQR: 689 – 953
- Max: 3,422 (potential doublets)

[[0,1]]
**`total_counts`**
- Median: **2,195 UMIs**
- 50% IQR: 1,753 – 2,762
- Max: 15,844

[[0,2]]
**`pct_counts_mt`**
- Median: **2.03%**
- Bulk tightly packed below 4%
- Long apoptotic tail up to **22.57%**
:::

---
# Joint QC Scatter: Dissecting Cell Pathology

Plotting `total_counts` vs. `n_genes_by_counts` colored by mitochondrial percentage:

![QC Scatter](lecture_05_cellular_individuality/images/qc_scatter.png){width=50% center}

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Diagonal Trend
- Healthy cells exhibit a tight linear correlation: more sequencing depth reveals more unique genes.

[[0,1]]
### The Two Failure Modes
- **Top Right**: Abnormally high counts ($> 6,000$) and genes ($> 2,500$) $\to$ **Doublets / Multiplets**.
- **Lower MT outliers** (orange/red points): High MT % with low library size $\to$ **Apoptotic ghosts**.
:::

---

# Biology of Mitochondrial Contamination

Why is elevated mitochondrial RNA the hallmark of dying cells?

:::matrix { cols="1.2fr 1fr" gap="30px" }
[[0,0]]
```text
Healthy Intact Cell                  Apoptotic / Ruptured Cell
┌─────────────────────────┐           ░░░░░░░░░░░░░░░░░░░░░░░░░░░░
│ Cytoplasm: [mRNA] [mRNA]│           ░ mRNA diffuses into       ░
│                         │           ░ ambient buffer           ░
│ ┌──────────────────┐    │           ░           ┌─────────┐    ░
│ │  Mitochondrion   │    │           ░           │  Mito   │    ░
│ │    [mt-mRNA]     │    │     ──►   ░           │ [mt-mRNA│    ░
│ └──────────────────┘    │           ░           └─────────┘    ░
│ Intact plasma membrane  │           ░ Cytoskeletal husk traps  ░
└─────────────────────────┘           ░ organelle inside droplet ░
                                      ░░░░░░░░░░░░░░░░░░░░░░░░░░░░
Total RNA: 98% cytoplasmic              Total RNA: 20-50% mito!
```

[[0,1]]
### The Mechanism
1. Cellular stress or mechanical shearing damages the plasma membrane.
2. Small, soluble cytoplasmic mRNAs diffuse into the buffer.
3. Mitochondria are large ($\sim 1~\mu\text{m}$), double-membraned organelles anchored to the cytoskeleton.
4. When encapsulated into a droplet, the damaged ghost retains its mitochondria but has lost its cytoplasm.
:::

---

# Tool 4 — Why NOT Standard Deviation?

Why shouldn't we filter outliers using $\text{Mean} + 3~\sigma$?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Flaw of Standard Deviation
$$\sigma = \sqrt{\frac{1}{N-1} \sum_{i=1}^N (x_i - \bar{x})^2}$$
- Squaring distances gives extreme leverage to outliers.
- The **breakdown point** of sample standard deviation is **0%**: a single extreme cell (e.g. $80\%$ MT) can inflate $\sigma$ arbitrarily high!
- On PBMC3k ($N = 2,708$):
  - $\bar{x} = 2.23\%$, $\sigma = 1.24\%$
  - $\bar{x} + 3\sigma = 2.23\% + 3(1.24\%) = \mathbf{5.94\%}$
  - Only flags **29 cells** (1.07%), missing real dying cells!

[[0,1]]
### Outlier Masking Effect

```text
Contaminated sample:
[2%, 2%, 2%, 3%, 8%, 15%, 22%]
                  ▲
Outliers inflate mean & SD!
Threshold moves far to the right,
failing to flag the 8% and 15% cells.
```

The very presence of dying cells **masks** their own detection when using non-robust statistics!
:::

---

# Robust Statistics — Median Absolute Deviation

The Median Absolute Deviation (MAD) is the robust counterpart to standard deviation:

$$\text{MAD} = \text{median}\left(\lvert x_i - \text{median}(x)\rvert\right)$$

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Key Mathematical Properties
- **50% Breakdown Point**: Up to half the observations can be corrupted without ruining the scale estimate.
- Operates on medians, completely ignoring extreme values at the tails.
- Raw MAD measures dispersion in data units.

[[0,1]]
### The Normal Consistency Factor ($k$)
To make MAD an unbiased estimator of standard deviation $\sigma$ for Gaussian data:

$$k = \frac{1}{\Phi^{-1}(0.75)} = \frac{1}{\sqrt{2}~\text{erf}^{-1}(0.5)} \approx \mathbf{1.4826022}$$

$$\text{MAD}_{\text{scaled}} = 1.4826 \times \text{MAD}$$

For $X \sim \mathcal{N}(\mu, \sigma^2)$, $\mathbb{E}[\text{MAD}_{\text{scaled}}] = \sigma$.
:::

```python
pct_mt = adata.obs['pct_counts_mt'].values
med_mt = np.median(pct_mt)
mad_mt = np.median(np.abs(pct_mt - med_mt))
mad_scaled = mad_mt * 1.4826022
```

---
# Dynamic Outlier Threshold in Action

Applying MAD scaling to PBMC3k mitochondrial percentages ($N = 2,708$):

$$\text{Threshold} = \text{Median} + 3 \times \text{MAD}_{\text{scaled}} = 2.03\% + 3(0.8015\%) = \mathbf{4.43\%}$$

![MT Histogram Thresholds](lecture_05_cellular_individuality/images/mt_hist_thresholds.png){width=50% center}

:::matrix { cols="1fr 1fr 1fr 1fr" gap="15px" }
[[0,0]]
**Median + 3 MAD**
- Cutoff: **4.43%**
- Flags: **95 cells** (3.51%)
- *Data-driven & robust*

[[0,1]]
**Fixed 5% Rule**
- Cutoff: **5.00%**
- Flags: **60 cells** (2.22%)
- *Arbitrary rule of thumb*

[[0,2]]
**Mean + 3 SD**
- Cutoff: **5.94%**
- Flags: **29 cells** (1.07%)
- *Masked by outliers*

[[0,3]]
**Fixed 10% Rule**
- Cutoff: **10.00%**
- Flags: **8 cells** (0.30%)
- *Far too permissive*
:::

---

# The Proof: Outlier Stress Test

To prove the superiority of MAD, we simulate removing the top 5% extreme outliers:

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Metric | All Cells (100%) | Drop Top 1% Outliers | Drop Top 5% Outliers | Shift (% Change) |
| :--- | :--- | :--- | :--- | :--- |
| **Mean** | 2.23% | 2.15% | 2.05% | -8.1% |
| **Std Dev ($\sigma$)** | 1.24% | 0.91% | 0.75% | **-39.5%** |
| **$\text{Mean} + 3~\sigma$** | **5.94%** | **4.89%** | **4.30%** | **-27.6% (Unstable!)** |
| **Median** | 2.03% | 2.02% | 1.98% | -2.5% |
| **$\text{MAD}_{\text{scaled}}$** | 0.80% | 0.79% | 0.74% | -7.5% |
| **$\text{Median} + 3~\text{MAD}$** | **4.43%** | **4.38%** | **4.19%** | **-5.4% (Rock solid!)** |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Standard Deviation Collapses
Removing just 5% of outliers drops the $\text{Mean} + 3~\sigma$ threshold from **5.94% to 4.30%**!
The initial threshold was inflated by the very cells it failed to exclude.

[[0,1]]
### MAD Remains Stable
The MAD threshold shifts by only **0.24%**, demonstrating that its estimate reflects the true underlying biological distribution.
:::

---
# Phenotype of the Flagged Cells

Are the 95 cells flagged by `pct_counts_mt` $\ge 4.43\%$ really dead or dying?

### Comparative Statistics

| Feature | Intact Cells ($N = 2,613$) | Flagged Cells ($N = 95$) | Change |
| :--- | :--- | :--- | :--- |
| **Median UMIs** | **2,220** | **1,208** | **-45.6%** 📉 |
| **Median Genes** | **821** | **537** | **-34.6%** 📉 |
| **Median MT %** | **1.99%** | **5.35%** | **+$2.7\times$** 📈 |
| **Max MT %** | 4.43% | **22.57%** | — |

### Biological Verdict
- Flagged cells have nearly **half the total RNA** of intact cells.
- They express **300 fewer unique genes**.
- Their remaining transcriptome is heavily skewed toward mitochondrial respiration transcripts.
- This confirms classic apoptotic cell death with extensive cytoplasmic loss.

---
# Two-Sided MAD on Library Size & Gene Diversity

Modern recommendations (`sc-best-practices.org`) apply MAD on **log-transformed** counts:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Why Log-Transform First?
- Library size (`total_counts`) and gene count (`n_genes_by_counts`) span several orders of magnitude and are strongly right-skewed.
- Transform: $y = \log(1 + x)$.
- Define bounds at **5 MADs**:
$$
  [\text{Median} - 5\\; \text{MAD}\_{\text{scaled}}, \quad \text{Median} + 5\\; \text{MAD}\_{\text{scaled}}]
$$
- Exponentiate bounds back: $\exp(y) - 1$.

[[0,1]]{ .smaller-table }
### PBMC3k 5-MAD Filter Bounds

| Metric | Median | 5-MAD Interval | Excluded |
| :--- | :--- | :--- | :--- |
| **`total_counts`** | 2,195 | [405, 11,885] | 0 low, 2 high |
| **`n_genes`** | 816 | [248, 2,679] | 5 low, 5 high |

- The upper gene boundary ($\sim 2,679$) precisely flags multi-cell doublets!
- The lower boundary ($\sim 248$) flags near-empty or damaged fragments.
:::

---

# Tool 5 — Doublet Detection: Droplet Collisions

Microfluidic droplet encapsulation is governed by a **Poisson process**:

$$P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}$$

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Microfluidic Loading Realities
- To prevent multiplet encapsulation, cell suspensions are loaded at extreme dilution ($\lambda \approx 0.05$):
  - $P(0) \approx 95.1\%$ (empty droplets)
  - $P(1) \approx 4.76\%$ (singlets)
  - $P(\ge 2) \approx 0.12\%$ (doublets)
- Among cell-containing droplets, the doublet fraction is:
  $$\frac{P(\ge 2)}{P(\ge 1)} \approx 2.5\%$$
- Doublet rate scales linearly with the number of loaded cells ($\sim 0.8\%$ per 1,000 cells recovered).

[[0,1]]
### Types of Doublets
- **Homotypic**: Two cells of the same type (e.g. 2 CD4+ T-cells). Affects library size, but not clustering.
- **Heterotypic**: Two distinct cell types (e.g. B-cell + Monocyte). Produces an artificial "hybrid" expression state that distorts trajectory inference!
:::

---

# In Silico Doublet Simulation: Scrublet

Scanpy includes **Scrublet** (`sc.pp.scrublet`) for model-based doublet detection:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### How Scrublet Works
1. Randomly pairs observed transcriptomes and sums their counts to synthesize thousands of **artificial doublets**.
2. Embeds observed cells and simulated doublets together in PCA space.
3. Constructs a $k$-Nearest Neighbors ($k$-NN) graph.
4. Computes a **Doublet Score** ($0 \to 1$) based on the local fraction of artificial doublet neighbors.
5. Sets a bimodal score threshold.

[[0,1]]
### PBMC3k Scrublet Results
- Automatic score threshold: **0.228**
- Predicted doublets: **33 cells** (1.27%)
- Median score: 0.044
- Doublet characteristics:
  - Median counts: **2,624** (vs 2,214 singlets)
  - Median genes: **900** (vs 819 singlets)
:::

```python
import scanpy as sc
# Run Scrublet on raw counts
sc.pp.scrublet(adata, random_state=0)
print(f"Doublets detected: {adata.obs['predicted_doublet'].sum()}")  # 33
```

---

# ⚠ Pitfall — AnnData Views vs. Copies

When you subset an AnnData object, Python creates a **View** to conserve memory:

```python
# Subsetting creates a VIEW
adata_subset = adata[adata.obs.total_counts > 1000]

print(adata_subset.is_view)  # True
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Dangerous Mistake
Modifying `.obs`, `.var`, or `.X` on a view triggers:
```text
ImplicitModificationWarning: Trying to modify 
attribute `.obs` of view, initializing view as actual.
```
- Causes silent deep copies in memory.
- Can invalidate pointers or cause subtle slice reference bugs.

[[0,1]]
### The Safe Solution
Always make an explicit copy when filtering:

```python
# Explicit copy detaches memory
adata_clean = adata[mask, :].copy()

print(adata_clean.is_view)  # False
```
:::

> **Golden Rule**: Always append `.copy()` when assigning a filtered subset of an AnnData object.

---

# AnnData Disk Storage: The `.h5ad` Format

AnnData serializes the entire experiment (matrix, metadata, embeddings, graphs) into a single Hierarchical Data Format (HDF5) container:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Saving and Loading

```python
# Save clean dataset with gzip compression
adata_clean.write("pbmc3k_clean.h5ad", compression="gzip")

# Reload in seconds in any Python session
adata_reloaded = sc.read_h5ad("pbmc3k_clean.h5ad")
```

- Self-contained: preserves `.obs`, `.var`, `.uns`, and PCA/UMAP layers.
- Interoperable with R via `anndata` / `zellkonverter` (converts to `SingleCellExperiment` or `Seurat`).

[[0,1]]{ .smaller-table }
### On-Disk Compression Benchmarks

| Object | Droplets | Compression | File Size |
| :--- | :--- | :--- | :--- |
| **Raw** | 737,280 | None | 58.97 MB |
| **Raw** | 737,280 | **gzip** | **33.45 MB** |
| **Clean Filtered** | ~2.6k | None | 20.92 MB |
| **Clean Filtered** | ~2.6k | **gzip** | **7.11 MB** |
:::

---

# Biological Reality — Extreme Sparsity

What does a quality-filtered single-cell dataset look like?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Matrix Metrics in Real Data
- After removing empty droplets and dying cells:
  - Total entries: $\sim 2.6\text{k cells} \times 32,738\text{ genes} \approx 85\text{ million}$
  - Non-zero entries: $\approx 2.2\text{ million}$
  - Matrix Density: **$\approx 2.6\%$**
- **Over 97% of all matrix entries remain exact zeros!**
- Median detected genes per cell: **$\sim 800$** ($< 3\%$ of the annotated human genome).
- Almost **half of all annotated genes** have zero reads across the entire cell population.

[[0,1]]
### Three Sources of Zeros
1. **Biological Silence**: PBMCs are differentiated blood cells; neuronal and liver genes are locked in heterochromatin.
2. **Transcriptional Bursting**: Genes transcribe in episodic pulses; many active genes have zero mRNA molecules at any instant.
3. **Technical Dropout**: Reverse transcription capture efficiency is only $\sim 10\% - 20\%$.
:::

---
# Zeros vs. Poisson Dropout

Is single-cell sparsity purely technical failure, or does it reflect Poisson sampling?

![Zeros vs Poisson](lecture_05_cellular_individuality/images/zeros_vs_poisson.png){width=60% center}

If molecule capture is Poisson distributed with mean expression $\mu$, the probability of observing 0 reads is:

$$P(X = 0) = \frac{\mu^0 e^{-\mu}}{0!} = \mathbf{e^{-\mu}}$$

---
# Interpreting Dropout Kinetics

Comparing observed zero fractions against the theoretical Poisson prediction $e^{-\mu}$:

| Mean Expression ($\mu$) Bin | Gene Count | Observed Zeros | Poisson $\mathbf{e^{-\mu}}$ | Median $\frac{\text{Var}}{\mu}$ |
| :--- | :--- | :--- | :--- | :--- |
| $[0.001, 0.01)$ | 4,679 | **0.9969** | 0.9966 | 1.00 |
| $[0.01, 0.1)$ | 6,612 | **0.9716** | 0.9698 | 1.05 |
| $[0.1, 1)$ | 2,129 | **0.8462** | 0.8278 | 1.97 |
| $[1, 10)$ | 222 | **0.2556** | 0.1011 | **3.78** |
| $\ge 10$ | 43 | **0.0069** | 0.0000 | **6.84** |

### Key Biological Takeaways
- **Low-expression genes ($\mu < 0.1$)**: Observed zeros match Poisson expectation almost perfectly. Dropout is driven by physical sampling limits!
- **High-expression genes ($\mu \ge 1$)**: Observed zeros significantly exceed Poisson predictions ($\text{Var}/\mu \gg 1$).
- This excess variance reflects true **biological heterogeneity**: cell-type-specific expression and transcriptional bursting!

---
# Normalization Preview (Bridge to Lecture 6)

Why can't we cluster cells directly on raw counts?

### The Library Size Artifact
- Sequencing depth varies from **506 to 15,844 UMIs** per cell due to cell size differences and technical capture variance.
- Without normalization, PCA and UMAP group cells by sequencing depth rather than biological cell type!

```python
# 1. Scale total counts per cell to 10,000 (CP10k)
sc.pp.normalize_total(adata, target_sum=1e4)

# 2. Log-transform counts: y = ln(1 + x)
sc.pp.log1p(adata)
```

### Value Transformations

| Stage | Mean Non-Zero Entry | Max Entry |
| :--- | :--- | :--- |
| **Raw Counts** | 2.79 | 419.0 |
| **CP10k Normalized** | 11.76 | 1,751.1 |
| **$\log(1+x)$ Transformed** | **2.02** | **7.47** |

Log-transformation dampens hyper-abundant structural genes and stabilizes variance across expression ranges.

---

# Combining Filters in AnnData (Homework Task 1)

In Homework Task 1, you filter `adata` to create a clean matrix `adata_clean`:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Multi-Condition Filtering Logic
- To select valid cells, combine multiple boolean masks using Python's bitwise AND operator `&`:
  - Condition 1: Keep cells whose mitochondrial percentage (`pct_counts_mt`) is strictly below the dynamic MAD threshold.
  - Condition 2: Keep cells whose gene diversity (`n_genes_by_counts`) does not exceed the doublet cutoff ($\le 2,500$).
- **Operator Precedence**: In Python, `&` has higher precedence than comparison operators (`<`, `<=`). **Always wrap each condition in parentheses!**
  ```python
  (condition_1) & (condition_2)
  ```

[[0,1]]
### Slicing & Memory Isolation
- AnnData objects use 2D indexing `[rows, columns]`:
  - `adata[cell_mask, :]`: selects cells matching the mask across all genes.
- **The `.copy()` Requirement**:
  - Slicing without `.copy()` produces an in-memory `AnnData` view.
  - You must call `.copy()` so `adata_clean` is an independent object that can be modified, normalized, and saved to disk.
:::

---

# Implementing Robust Estimators (Homework Task 2)

In Homework Task 2, you implement `compute_mad(data_array)` in pure Python/NumPy:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Step-by-Step Algorithm
$$\text{MAD} = \text{median}\left(\lvert X_i - \text{median}(X)\rvert\right)$$

1. **Center**: Compute the median of the input array using `np.median()`.
2. **Absolute Deviations**: Compute absolute residuals $\lvert X_i - M\rvert$ using `np.abs()`.
3. **Dispersion**: Compute the median of those absolute deviations.
4. **Normal Scale Factor**: Multiply by **$1.4826$** to calibrate the estimator for a standard normal distribution.

[[0,1]]
### How to Verify Your Implementation
You can verify your custom function against SciPy's built-in robust statistics module:

```python
from scipy.stats import median_abs_deviation

# Reference value using SciPy's normal scaling
expected = median_abs_deviation(data, scale='normal')

# Your implementation should match to full precision!
assert np.isclose(compute_mad(data), expected)
```

The scaling constant comes from the standard normal quantile: $k = 1 / \Phi^{-1}(0.75)$.
:::

---

# Identifying Transcriptome Drivers (Homework Task 3)

In Homework Task 3, you find which gene is the most highly expressed in `adata_clean`:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Summing Expression: Axis Geometry
- Slicing and aggregating sparse AnnData matrices:
  - `adata.X.sum(axis=1)`: sums across columns $\to$ total UMIs per **cell** (row sum).
  - `adata.X.sum(axis=0)`: sums across rows $\to$ total UMIs per **gene** (column sum).
- **Matrix Flattening**:
  - SciPy sparse operations return a 2D matrix of shape `(1, n_vars)`.
  - Use `.A1` or `np.asarray(...).ravel()` to flatten into a 1D vector of gene totals.

[[0,1]]
### Locating the Winner
- How to map array indices back to gene symbols:
  - `np.argmax(gene_totals)`: returns the index of the highest value.
  - `adata_clean.var_names[index]`: extracts the corresponding HGNC gene symbol.
- **Biological Context**:
  - Eukaryotic cells are dominated by structural lncRNAs, translation machinery (ribosomal proteins `RPS`/`RPL`), and major cytoskeletal transcripts.
  - Your code must identify which specific transcript captures the top spot!
:::

---

# The Mathematics of UMI Collisions (Homework Task 4)

Unique Molecular Identifiers (UMIs) index original cDNA molecules before PCR amplification:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Random Sampling Model
- A UMI barcode of length $L$ nucleotides has $4^L$ possible sequences:
  $$N = 4^L\text{ possible unique barcodes}$$
- For $L = 10$, $N = 4^{10} = \mathbf{1,048,576}$.
- When $k$ molecules of a specific transcript are captured in a single droplet, each molecule randomly and independently samples one UMI from $N$ with replacement.
- This is the **classical occupancy problem** (coupon collector / birthday paradox).

[[0,1]]
### The Occupancy Expectation ($E$)
The probability that a specific UMI is NOT chosen by any of the $k$ molecules is $(1 - 1/N)^k \approx e^{-k/N}$.
The expected number of distinct UMIs observed ($E$) is:

$$E = N \left(1 - \left(1 - \frac{1}{N}\right)^k\right) \approx \mathbf{N \left(1 - e^{-k/N}\right)}$$

- When $k \ll N$, $E \approx k$ (collisions are rare).
- When $k$ is large, collisions cause **under-counting**: you observe $E < k$.
- The number of lost UMIs is $k - E$.
:::

---

# UMI Capacity & Collision Correction

Why did sequencing platforms transition from 10-nt to 12-nt UMIs?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Barcode Capacity Comparison
- **10-nt UMI (10x v2)**:
  $$N = 4^{10} = \mathbf{1,048,576}$$
- **12-nt UMI (10x v3)**:
  $$N = 4^{12} = \mathbf{16,777,216}$$
- A 2-nucleotide increase expands the barcode diversity by **$16\times$**, dramatically reducing collision probabilities for highly expressed genes.

[[0,1]]
### Inverting the Collision Formula
If a sequencing experiment observes $E$ distinct UMIs for an abundant gene, we can estimate the true molecule count $k$:

$$\frac{E}{N} = 1 - e^{-k/N} \implies e^{-k/N} = 1 - \frac{E}{N}$$

$$k_{\text{corrected}} = -N \ln\left(1 - \frac{E}{N}\right)$$

> In Homework Task 4, use the formula for $E$ to compute how many *fewer* UMIs are observed when $k = 50,000$ molecules are captured with $L = 10$.
:::

---

# Biological & Instrumental Interpretation (Homework Task 5)

Task 5 asks you to explain the physical and biological mechanisms in your own words:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Key Concepts to Synthesize
1. **The Knee Phenomenon**:
   - Compare droplet volume ($\sim 1\text{ nL}$) vs. cell cytoplasm volume ($\sim 1\text{ pL}$). What happens to RNA density in empty droplets containing only cell-free ambient soup?
2. **Estimator Robustness**:
   - Why does standard deviation have a 0% breakdown point? How do squared deviations inflate $\sigma$ and mask dying cells?
3. **Mitochondrial Contamination**:
   - Trace what happens during apoptosis: cytoplasmic membrane permeabilization vs. organelle anchoring in the cytoskeleton.

[[0,1]]
### Key Concepts to Synthesize (cont.)
4. **Doublet Signatures**:
   - What happens when two distinct cells are encapsulated into a single oil droplet? How does their combined transcriptome affect `total_counts` and `n_genes_by_counts`?
5. **PCR Deduplication**:
   - Why must UMIs be attached *before* PCR amplification during cDNA reverse transcription? What would happen if barcodes were added after PCR cycles?
:::

---

# Putting It All Together: The Preprocessing Workflow

```python
import numpy as np, scanpy as sc

# 1. Ingest raw 10x matrix and cache to disk
adata = sc.read_10x_mtx('raw_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
adata.var_names_make_unique()

# 2. Initial cell calling: drop ambient empty droplets (< 500 UMIs)
adata = adata[adata.X.sum(axis=1).A1 > 500, :].copy()

# 3. Calculate mitochondrial QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 4. Compute robust outlier threshold via MAD
pct_mt = adata.obs['pct_counts_mt'].values
median_mt = np.median(pct_mt)
mad_scaled = np.median(np.abs(pct_mt - median_mt)) * 1.4826022
threshold_mt = median_mt + (3 * mad_scaled)

print(f"Dynamic MT Outlier Threshold: {threshold_mt:.2f}%")

# ---- HOMEWORK TASKS ----
# Task 1: Filter adata into adata_clean using threshold_mt and gene count limits
# Task 2: Implement your custom compute_mad() function
# Task 3: Identify the most highly expressed gene in adata_clean
# Task 4: Calculate UMI collision probabilities and lost counts for k=50,000
# Task 5: Answer the interpretation questions in markdown
```

---

# Summary

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### AnnData & Sparsity
#### Data Structure
- Transposed to Cells × Genes
- > 97% sparse matrix in RAM
- CSC memory: 22 MB vs 90 GB
- `adata.X.sum()` returns 2D matrix
- Always `.copy()` AnnData views
- Fast `.h5ad` HDF5 serialization

[[0,1]]{ .green }
### Cell Calling & QC
#### Microfluidics
- Knee plot separates cells from soup
- Top 2,700 droplets hold 90% of reads
- 13 human mitochondrial genes (`MT-`)
- Apoptosis causes cytoplasmic leakage
- Ambient soup is 6.6% MT vs 2.1% in cells
- Doublets identified by high gene count

[[0,2]]{ .blue }
### Robust Statistics
#### Inference
- Mean / SD broken by outliers
- $\text{MAD} = \text{median}(\lvert x_i - \text{median}\rvert)$
- Gaussian scaling factor: $1.4826$
- Dynamic threshold: $\text{Med} + 3~\text{MAD}$
- UMI collision occupancy correction
- Clean dataset: ready for manifolds
:::

---

# Next Steps — Lecture 6

In **Lecture 6: Cellular State Spaces (Manifolds & Differentiation)**, we take the clean AnnData object and build the geometric cellular state space:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### What We Accomplished Today (Lecture 5)
- Ingested 737,280 raw microfluidic droplets.
- Screened out empty background droplets using the barcode rank knee plot.
- Quantified mitochondrial leakage and filtered apoptotic cells via dynamic MAD scaling.
- Filtered multiplet doublets, producing a verified, quality-filtered single-cell matrix.
- Corrected for PCR amplification artifacts and understood UMI collision statistics.

[[0,1]]
### What Comes Next (Lecture 6)
- **Feature Selection**: Finding Highly Variable Genes (HVGs) using dispersion vs. mean.
- **Dimensionality Reduction**: Principal Component Analysis (PCA) spectral projection.
- **Manifold Learning**: $k$-Nearest Neighbor ($k$-NN) graph construction.
- **Non-Linear Embedding**: Uniform Manifold Approximation and Projection (UMAP).
- **Unsupervised Community Detection**: Leiden graph partitioning to identify T-cells, B-cells, Monocytes, and NK cells.
:::
