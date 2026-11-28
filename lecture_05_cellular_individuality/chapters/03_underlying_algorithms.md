:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Sparse Matrices · Cell Calling · Robust MAD · Doublets · UMI Collisions
:::

---

# Curriculum & Algorithmic Roadmap

Six quantitative algorithms turn 737,280 noisy microfluidic droplet barcodes into a clean, sparse single-cell expression matrix ready for Part 4 (`Scanpy` + `AnnData` on 10x Genomics **PBMC3k**):

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "raw",  "label": "Raw Droplets\n737k × 33k (.mtx)",  "x": 20,   "y": 40, "width": 160, "height": 64 },
    { "id": "sp",   "label": "Sparse CSC / CSR\nO(nnz) Compression","x": 215,  "y": 40, "width": 165, "height": 64, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "cc",   "label": "Knee Plot +\nEmptyDrops Test",       "x": 415,  "y": 40, "width": 160, "height": 64, "style": "fill: #2563EB; color: #fff;" },
    { "id": "mad",  "label": "Robust MAD QC\n(1.4826 Scaling)",    "x": 610,  "y": 40, "width": 160, "height": 64, "style": "fill: #16A34A; color: #fff;" },
    { "id": "dbl",  "label": "Poisson Loading +\nScrublet k-NN",   "x": 805,  "y": 40, "width": 165, "height": 64, "style": "fill: #D97706; color: #fff;" },
    { "id": "umi",  "label": "Dropout Kinetics +\nUMI Collisions", "x": 1005, "y": 40, "width": 165, "height": 64, "style": "fill: #DC2626; color: #fff;" },
    { "id": "norm", "label": "CP10k + log1p\n(Bridge to L6)",      "x": 1205, "y": 40, "width": 160, "height": 64, "style": "fill: #0891B2; color: #fff;" }
  ],
  "links": [
    { "from": "raw", "to": "sp"   },
    { "from": "sp",  "to": "cc"   },
    { "from": "cc",  "to": "mad"  },
    { "from": "mad", "to": "dbl"  },
    { "from": "dbl", "to": "umi"  },
    { "from": "umi", "to": "norm" }
  ]
}
:::

:::matrix { .cards cols="1fr 1fr 1fr" gap="18px" }
[[0,0]]{ .purple }
### 1. Sparsity & Cell Calling
COO vs. CSC/CSR indexing (`data`, `indices`, `indptr`), $4{,}000\times$ RAM savings, barcode-rank knee physics, and `EmptyDrops` ambient-profile testing.

[[0,1]]{ .green }
### 2. Robust QC & Doublets
Finite-sample breakdown ($0\\%$ SD vs. $50\\%$ MAD), deriving $k = 1.4826$, log-scale two-sided MAD, Poisson multiplet loading, and `Scrublet` $k$-NN scoring.

[[0,2]]{ .blue }
### 3. Dropout, UMIs & Scaling
Poisson ($e^{-\mu}$) vs. NB zero fractions, classical UMI occupancy & collision inversion $\hat{k}$, and `CP10k` + $\log(1+x)$ variance stabilization.
:::

---
# Sparse Matrix Compression: COO, CSC & CSR

![Sparse Matrix RAM Scaling and CSC Array Layout](lecture_05_cellular_individuality/images/sparse_storage_scaling.png){ width=90 center }

:::matrix { cols="50/50" gap="26px" }
[[0,0]]
- **Storage equations** ($n_{\text{obs}} = 737{,}280$, $n_{\text{vars}} = 32{,}738$, $\text{nnz} = 2{,}905{,}861$, **99.988%** zeros):
$$
\begin{aligned}
\text{Dense f32} = 4 \cdot n_{\text{obs}} n_{\text{vars}} \approx 89.9\text{ GB},\\\\ \text{CSC} = 8 \cdot \text{nnz} + 4(n_{\text{vars}} + 1) \approx 22.3\text{ MB}
\end{aligned}
$$
- Disk files (`matrix.mtx`) use 1-based **COO** triplets $(i, j, v)$; RAM uses **CSC** (fast per-gene column sums) or **CSR** (fast per-cell row sums) in $O(\text{nnz})$ time.

[[0,1]]
- **Access & memory semantics (`AnnData`):**
  - Gene $j$ occupies contiguous slice `data[indptr[j] : indptr[j+1]]`; empty genes have `indptr[j] == indptr[j+1]` (zero storage!).
  - Slicing an `AnnData` returns a lightweight **View** referencing parent arrays — mutating metadata requires an explicit **Copy** (`.copy()`) to detach memory safely.
:::

---

# Droplet Cell Calling: Barcode Rank ("Knee") Geometry

Why do fewer than **0.4%** of microfluidic barcodes capture over **90%** of all sequenced mRNA molecules in PBMC3k?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Physics of the Three Regimes
- **1. Cellular Plateau (Ranks 1–2,000):** An intact single cell ($\sim 1\text{ pL}$ cytoplasm) brings a dense packet of $1{,}500\text{–}16{,}000$ distinct UMIs.
- **2. The Cliff / Knee (Ranks 2,000–3,000):** Total UMIs plunge an order of magnitude across just 1,000 barcodes (Rank 2,000: `1,775` $\to$ Rank 2,700: `548` $\to$ Rank 3,000: `89`).
- **3. Ambient Soup Tail (Ranks > 3,000):** An empty oil droplet ($\sim 1\text{ nL}$, $1{,}000\times$ larger than a cell) passively traps cell-free RNA leaked from lysed cells (`78.5%` have `0` UMIs; `155,928` capture soup).

[[0,1]]
### Order-of-Magnitude vs. Fixed Cutoffs
- **Cell Ranger v2 heuristic** (given expected cell count $N = 3{,}000$):
$$
m = \text{percentile}\_{99}(\text{top } N) = 5{,}688.7, \qquad T = \frac{m}{10} = 568.9\text{ UMIs}
$$
  Retains **2,693** cells (`2,700` official filtered; Practice 5's `> 500` keeps **2,708**).
- **Fundamental limitation of rank cutoffs:** Any single horizontal line on the knee plot treats all low-UMI droplets identically — silently discarding genuinely small, low-RNA cell types (e.g., platelets, resting lymphocytes, neutrophils)!
:::

---
# Statistical Cell Calling: `EmptyDrops` & Ambient Soup

![PBMC3k Barcode Rank Curve and EmptyDrops Ambient Test](lecture_05_cellular_individuality/images/cell_calling_emptydrops.png){ width=86 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **The ambient fingerprint:** Pooling low-UMI droplets ($10 \le \text{UMI} \le 100$) estimates the background soup composition $\hat{\mathbf{p}}\_{\text{amb}}$ — in PBMC3k, empty droplets contain **$6.62\\%$ `MT-` RNA** ($3\times$ higher than intact cells' $2.15\\%$) from pre-encapsulation lysis!
- **Null hypothesis ($H_0$):** An empty droplet's count vector $\mathbf{c}$ is drawn from the ambient soup via a **Dirichlet–Multinomial** distribution $\operatorname{DM}(\alpha_0 \hat{\mathbf{p}}\_{\text{amb}})$ (accounting for overdispersion).

[[0,1]]
- **Likelihood rescue (Lun et al., 2019):**
  1. Droplets above the knee are kept immediately as unambiguous cells.
  2. Each transitional droplet below the knee is tested against $\hat{\mathbf{p}}\_{\text{amb}}$ by Monte Carlo $p$-value + Benjamini–Hochberg FDR.
  3. Droplets whose expression profile **significantly deviates** from the soup (`FDR < 0.01`) are rescued as real low-RNA cells!
:::

---

# Interactive Lab — Sparse Indexing & Droplet Cell Calling

Explore sparse matrices and cell calling live: **(1)** in **Tab 1**, click tiles in the `6 × 6` grid to watch `CSC` vs. `CSR` `indptr` update in real time, then test the **Raw PBMC3k (`737k`)** vs. **Filtered (`2.6k`)** RAM presets; **(2)** in **Tab 2**, slide **Expected Cells $N$** and **Hard Cutoff $T$** on PBMC3k's knee curve and watch **EmptyDrops** rescue low-RNA platelets and lymphocytes while rejecting ambient-rich empty droplets!

<iframe src="lecture_05_cellular_individuality/demo/sparse_and_cell_calling_lab.html" title="Interactive Sparse Matrix Indexing and Droplet Cell-Calling Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---
# Why Standard Deviation Fails: Breakdown Point & Masking

Why do modern single-cell pipelines reject classical $\bar{x} + 3s$ thresholding for quality control?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Finite-Sample Breakdown Point ($\varepsilon^*$)
- The **breakdown point** is the smallest fraction of corrupted observations that can drive an estimator to $\pm \infty$:
$$
\varepsilon^\*(\bar{x}, s) = \frac{1}{n} \to \mathbf{0\\%}, \qquad \varepsilon^*(\text{Med}, \text{MAD}) = \frac{\lfloor (n-1)/2 \rfloor}{n} \to \mathbf{50\\%}
$$
- Because $s = \sqrt{\frac{1}{n-1}\sum_i (x_i - \bar{x})^2}$ **squares** deviations, a handful of ruptured apoptotic ghosts (`20–80% MT-`) dominates the sum and drags $\bar{x} + 3s$ far rightward.

[[0,1]]
### Outlier Masking in Real PBMC3k (`N = 2,708`)
- Dying cells literally **mask** their own detection under $\bar{x} + 3s$:
  - Raw PBMC3k: $\bar{x} + 3s = 2.23\\% + 3(1.24\\%) = \mathbf{5.94\\%}$ (flags only **29** cells, missing dozens of $4.5\text{–}5.9\\%$ dying cells!).
  - Drop just the top **5%** extreme cells $\implies \bar{x} + 3s$ collapses **$-27.6\\%$** (`5.94% → 4.30%`), proving the original threshold was severely inflated!
  - Under the exact same 5% drop, $\text{Med} + 3 \ \text{MAD}_{\text{scaled}}$ moves only **$-5.4\\%$** (`4.43% → 4.19%`).
:::

---

# Median Absolute Deviation (MAD) & the $1.4826$ Constant

![Normal Consistency Factor 1.4826 and Outlier Breakdown Curve](lecture_05_cellular_individuality/images/mad_breakdown_and_normal_scale.png){ width=94% center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Definition of raw MAD:**
$$
\text{MAD}(x) = \operatorname{median}\_i \bigl\lvert x_i - \operatorname{median}(x) \bigr\rvert
$$
- For $X \sim \mathcal{N}(\mu, \sigma^2)$, the standardized variable $Z = (X - \mu)/\sigma \sim \mathcal{N}(0, 1)$ satisfies $P(\lvert Z \rvert \le c) = \frac{1}{2} \iff \Phi(c) = 0.75$:

[[0,1]]
- **Gaussian consistency factor $k$:**
$$
k = \frac{1}{\Phi^{-1}(0.75)} = \frac{1}{\sqrt{2} \ \operatorname{erf}^{-1}(0.5)} \approx \mathbf{1.4826022}
$$
- Setting $\text{MAD}\_{\text{scaled}} = k \cdot \text{MAD}(x)$ ensures $\mathbb{E}[\text{MAD}\_{\text{scaled}}] = \sigma$ on clean Gaussian data while preserving a **$50\\%$ breakdown point**!
:::

---

# One-Sided vs. Log-Scale Two-Sided MAD Filtering

Different single-cell QC metrics have distinct physical failure asymmetries — requiring matched one-sided vs. log-transformed two-sided MAD rules:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### One-Sided MAD on Mitochondrial Leakage (`pct_counts_mt`)
- **Biological mechanism:** Plasma-membrane rupture during apoptosis leaks small soluble cytoplasmic mRNAs into the buffer while large ($\sim 1 \ \mu\text{m}$), cytoskeleton-anchored mitochondria stay trapped inside the cell husk.
- Only **high** `pct_counts_mt` signals cell death (0% MT is healthy, never an outlier!) $\implies$ **one-sided upper cutoff**:
$$
x > \text{Med}(x) + n_{\text{MAD}} \cdot \text{MAD}\_{\text{scaled}}(x) \quad (n_{\text{MAD}} = 3)
$$
- On PBMC3k: $2.03\\% + 3(0.8015\\%) = \mathbf{4.43\\%}$ flags **95** cells with `-45.6%` median UMIs and `-34.6%` fewer detected genes!

[[0,1]]
### Log-Scale Two-Sided MAD (`total_counts`, `n_genes`)
- Library size and detected genes are **strictly positive and right-skewed** (spanning orders of magnitude). Applying linear MAD puts the lower bound below zero!
- Transform $y = \ln(1 + x)$ first, apply **two-sided 5-MAD**, and exponentiate back via $\exp(y) - 1$:
$$
\lvert y - \text{Med}(y) \rvert > 5 \cdot \text{MAD}\_{\text{scaled}}(y)
$$
- On PBMC3k (`sc-best-practices` rule):
  - `total_counts`: **`[405, 11,885]`** (drops tiny fragments & giant multiplets)
  - `n_genes_by_counts`: **`[248, 2,679]`** (upper bound flags doublets automatically!)
:::

---

# Microfluidic Encapsulation & Poisson Multiplet Kinetics

Why must ~95% of all droplets in a 10x Chromium run remain completely empty?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Independent Poisson Loading
- Flowing diluted cells randomly into $\sim 1\text{ nL}$ oil droplets at mean rate $\lambda$ cells/droplet gives exact **Poisson occupancy**:
$$
P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}, \qquad k = 0, 1, 2, \dots
$$
- **Conditional multiplet rate** among *non-empty* (cell-containing) droplets:
$$
P(X \ge 2 \mid X \ge 1) = \frac{1 - e^{-\lambda} - \lambda e^{-\lambda}}{1 - e^{-\lambda}} \approx \frac{\lambda}{2} \quad (\lambda \ll 1)
$$
- At 10x's target $\lambda \approx 0.05$: $P(0) = 95.1\\%$ empty, $P(1) = 4.76\\%$ singlets, and $P(X \ge 2 \mid X \ge 1) \approx \mathbf{2.5\\%}$ ($\sim 0.8\\%$ per $1{,}000$ recovered cells).

[[0,1]]
### Homotypic vs. Heterotypic Doublets
- Suppose a tissue has cell types $c = 1, \dots, C$ at proportions $\pi_c$ ($\sum_c \pi_c = 1$). When two independent cells co-encapsulate:
  - **Homotypic** (same type, probability $\sum_c \pi_c^2$): looks like a $2\times$ deeper singlet; barely distorts clustering or trajectories.
  - **Heterotypic** (two distinct types, probability $1 - \sum_c \pi_c^2$): creates a **fake hybrid state** (e.g., B-cell + Monocyte) that mimics a nonexistent transitional lineage!
- Simple `n_genes > 2,500` cutoffs miss low-UMI heterotypic doublets (e.g., two small lymphocytes totaling `1,800` genes) — motivating simulation-based detection!
:::

---

# Loading Trade-Offs & Scrublet Doublet Geometry

![Poisson Encapsulation Curves and Scrublet 2D PCA Doublet Bridges](lecture_05_cellular_individuality/images/poisson_loading_and_scrublet.png){ width=94% center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **The physical trade-off:** Pushing loading rate $\lambda$ from `0.05` to `0.25` increases singlet yield $4\times$, but inflates the conditional multiplet rate from **$2.5\\%$ to $12\\%$** — forcing super-loading protocols to rely on computational doublet removal or cell hashing.

[[0,1]]
- **Why geometry catches heterotypic hybrids:** Summing two different cell types (`T-cell + B-cell`) places the combined normalized transcriptome along a **bridge between clusters** in PCA space, completely separated from true singlets — whereas homotypic doublets land inside existing clusters.
:::

---

# In Silico Doublet Scoring Algorithm (`Scrublet`)

**Scrublet** (Wolock et al., 2019; `sc.pp.scrublet`) turns doublet detection into a supervised-by-simulation local density problem inside PCA space:

:::static-diagram { width=94% center }
{
  "nodes": [
    { "id": "obs", "label": "Observed Matrix\n(N real cells)",          "x": 20,  "y": 35, "width": 165, "height": 62 },
    { "id": "sim", "label": "Pairwise Sum A+B\n(r·N sim doublets)",      "x": 225, "y": 35, "width": 175, "height": 62, "style": "fill: #D97706; color: #fff;" },
    { "id": "pca", "label": "Joint CP10k +\nPCA Embedding",              "x": 440, "y": 35, "width": 165, "height": 62, "style": "fill: #2563EB; color: #fff;" },
    { "id": "knn", "label": "Build k-NN Graph\n(Obs + Sim joint)",       "x": 645, "y": 35, "width": 170, "height": 62, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "sc",  "label": "Bayes Score q_i\nfrom k_sim / k",           "x": 855, "y": 35, "width": 165, "height": 62, "style": "fill: #16A34A; color: #fff;" },
    { "id": "thr", "label": "Bimodal Threshold\npredicted_doublet",      "x": 1060,"y": 35, "width": 170, "height": 62, "style": "fill: #DC2626; color: #fff;" }
  ],
  "links": [
    { "from": "obs", "to": "sim" },
    { "from": "sim", "to": "pca" },
    { "from": "pca", "to": "knn" },
    { "from": "knn", "to": "sc"  },
    { "from": "sc",  "to": "thr" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Local $k$-NN Doublet Probability
- Given simulated-to-observed ratio $r = N_{\text{sim}}/N_{\text{obs}}$, expected doublet rate $\hat{\rho}$, and $k_{\text{sim}}$ simulated doublets among cell $i$'s $k$ nearest PCA neighbors:
$$
q_i = \frac{k_{\text{sim}}}{k}, \qquad \mathcal{L}_d(i) = \frac{q_i \ \hat{\rho} / r}{1 - \hat{\rho} - q_i(1 - \hat{\rho} - \hat{\rho}/r)}
$$
- Bayes-calibrates the raw neighbor fraction $q_i$ into a posterior doublet probability $\mathcal{L}_d(i) \in [0, 1]$!

[[0,1]]
### PBMC3k Real Run Diagnostics (Part 4)
- **Bimodal score histogram:** Simulated doublets split into a low-score mode (homotypic, indistinguishable from singlets) and a high-score mode (heterotypic bridges).
- Automatic valley cutoff on PBMC3k: **`0.228`** $\implies$ **33 predicted doublets** (`1.27%`), with higher median UMIs (`2,624` vs `2,214`) and genes (`900` vs `819`).
:::

---

# Interactive Lab — Robust MAD Breakdown & Scrublet Doublets

Test outlier masking and microfluidic loading live: **(1)** in **Tab 1**, drag **Apoptotic Contamination** from `5%` up to `45%` and watch `Mean + 3 SD` mask dying cells while `Median + 3×1.4826×MAD` stays rock-steady until the `50%` breakdown cliff (also toggle **Log-Scale Two-Sided**); **(2)** in **Tab 2**, sweep Poisson loading $\lambda$ from `0.05` to `0.30` and adjust **Scrublet $k$-NN** to watch heterotypic bridges light up!

<iframe src="lecture_05_cellular_individuality/demo/robust_mad_and_doublets_lab.html" title="Interactive Robust MAD Breakdown and Microfluidic Doublet Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Sparsity, Dropout Kinetics & Biological Heterogeneity

Even after removing all empty droplets and dead cells, **over 97%** of entries in PBMC3k remain exact zeros (`~800` genes detected per cell out of `32,738`). Is this a broken assay, or expected capture physics?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Pure Technical Capture (Poisson Baseline)
- Reverse transcription captures only $q \approx 10\text{–}20\\%$ of mRNA molecules in a droplet. If a homogeneous population has mean captured count $\mu = q \cdot m_{\text{true}}$, Poisson thinning predicts:
$$
P_{\text{Poisson}}(X = 0) = e^{-\mu}, \qquad \frac{\operatorname{Var}(X)}{\mu} = 1
$$
- **When $\mu < 0.1$** (PBMC3k bins `[0.001, 0.1)`): observed zeros (`0.9969`, `0.9716`) match $e^{-\mu}$ (`0.9966`, `0.9698`) and $\operatorname{Var}/\mu \approx 1.0$ — **pure Poisson shot noise**, no mystical "extra dropout" mechanism needed!

[[0,1]]
### Biological Heterogeneity & Transcriptional Bursting (NB)
- When true expression varies across cells (cell-type marker silence + episodic Pol II bursting) with Gamma dispersion $\alpha > 0$, marginalizing gives the **Negative Binomial zero probability**:
$$
P_{\text{NB}}(X = 0) = (1 + \alpha\mu)^{-1/\alpha} > e^{-\mu}, \qquad \frac{\operatorname{Var}(X)}{\mu} = 1 + \alpha\mu
$$
- **When $\mu \in [1, 10)$** in PBMC3k: observed zeros (**`0.2556`**) exceed Poisson (`0.1011`) by $2.5\times$ with $\operatorname{Var}/\mu = 3.78$ — revealing **true biological cell-type structure** that powers Lecture 6's clustering!
:::

---

# UMI Deduplication & Classical Occupancy (Birthday) Math

Why are **Unique Molecular Identifiers (UMIs)** ligated *before* PCR — and what happens when too many molecules of the same gene compete for finite barcode sequences?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Pre-PCR Tagging & Urn Occupancy
- PCR amplifies molecules exponentially with severe sequence/GC bias ($10\times$ vs. $1{,}000\times$ copies). Tagging each original cDNA with a random $L$-nt UMI *before* PCR lets us collapse all reads sharing `(cell, UMI, gene)` back to **1 molecule**.
- With barcode space $N = 4^L$ and $k$ captured molecules of one gene drawn uniformly with replacement:
  - Chance a fixed UMI stays empty: $(1 - 1/N)^k \approx e^{-k/N}$.
  - **Expected distinct UMIs observed ($E$):**
$$
E = N\left(1 - \left(1 - \frac{1}{N}\right)^k\right) \approx N\left(1 - e^{-k/N}\right)
$$
  - **Under-counting loss:** $k - E$ (fractional loss $(k - E)/k$).

[[0,1]]
### Birthday Paradox vs. Bulk Saturation
- **Two distinct collision regimes:**
  1. **First collision (Birthday Paradox):** Probability of $\ge 1$ collision anywhere among $k$ molecules:
$$
P(\ge 1\text{ collision}) = 1 - \prod_{j=0}^{k-1}\left(1 - \frac{j}{N}\right) \approx 1 - \exp\left(-\frac{k(k-1)}{2N}\right)
$$
     Hits $50\\%$ at $k \approx \sqrt{2\ln 2 \cdot N} \approx 1.18\sqrt{N}$ (`302` for 8-nt; `1,206` for 10-nt; `4,824` for 12-nt) — yet losing 1 count out of 1,200 is only `0.08%` error!
  2. **Bulk saturation ($k \sim N$):** Substantial fractional bias arises only when $k/N$ becomes non-negligible ($E \to N$ ceiling).
:::

---
# Dropout Curves, UMI Saturation & Collision Inversion

![PBMC3k Dropout Curves and UMI Occupancy Saturation](lecture_05_cellular_individuality/images/umi_occupancy_and_dropout.png){ width=78 center }

:::matrix { cols="45/55" gap="26px" }
[[0,0]]
- **Why 10x upgraded 10-nt (v2) $\to$ 12-nt (v3):** Adding just 2 nucleotides expands $N = 4^L$ by **$4^2 = 16\times$** ($1.05\text{ M} \to 16.78\text{ M}$), keeping $k/N \ll 1$ even for ultra-abundant structural transcripts.
- **Sequencing errors in UMIs:** A 1-bp error creates a fake singleton UMI! Production counters (`Cell Ranger`, `UMI-tools`) merge Hamming-distance-1 UMIs whenever a parent UMI has $\ge 2\times$ higher read support (**directional adjacency graph**).

[[0,1]]
- **Inverse collision correction (Fu et al., 2011):** Given observed distinct count $E < N$, inverting $E/N = 1 - e^{-k/N}$ recovers unbiased true counts:
$$
\hat{k} = -N \ln\left(1 - \frac{E}{N}\right)
$$
- **⚠ Saturation blow-up:** As $E \to N$, $\frac{d\hat{k}}{dE} = \frac{1}{1 - E/N} \to \infty$ — small sampling noise in $E$ explodes $\hat{k}$, making large barcode space $N$ far superior to post-hoc inversion!
:::

---

# Interactive Lab — Dropout Kinetics & UMI Collisions

Explore single-cell sparsity and UMI occupancy live: **(1)** in **Tab 1**, raise NB overdispersion $\alpha$ or lower the expressing subpopulation fraction $\pi_{\text{on}}$ to match PBMC3k's observed zero bins (`[1, 10) → 25.6%`) and inspect `CP10k + log1p` sparsity preservation; **(2)** in **Tab 2**, compare `Short 8-nt`, `Heavy 8-nt Saturation`, `Birthday Onset`, and `10x v3 12-nt` presets to see where first collisions occur vs. where bulk saturation kicks in!

<iframe src="lecture_05_cellular_individuality/demo/dropout_and_umi_collision_lab.html" title="Interactive Dropout Kinetics and UMI Collision Occupancy Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Library-Size Scaling (CP10k) & `log1p` Stabilization

In PBMC3k, sequencing depth varies **$30\times$** across intact cells (`506` to `15,844` UMIs). Why does Lecture 6 start with `normalize_total` + `log1p` before PCA?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Two-Step Normalization (`CP10k` + `log1p`)
$$
\tilde{x}\_{cg} = \ln\left(1 + 10^4 \cdot \frac{x_{cg}}{\sum_{g'} x_{cg'}}\right)
$$
1. **Library-size scaling (`CP10k`):** Equalizes total cellular capture depth to $L = 10{,}000$ so two identical cells sequenced at `1,000` vs. `10,000` UMIs share identical coordinates (otherwise PC1 merely tracks sequencing depth!).
2. **Delta-method variance stabilization:** For $X$ with $\operatorname{Var}(X) = \mu + \alpha\mu^2$, applying $g(x) = \ln(1 + x)$ gives:
$$
\operatorname{Var}\bigl(\ln(1 + X)\bigr) \approx \frac{\mu + \alpha\mu^2}{(1 + \mu)^2} \xrightarrow{\mu \gg 1} \alpha \quad (\text{constant!})
$$

[[0,1]]
### Why `log1p` Preserves Sparse Matrix RAM
- **Dynamic range compression (PBMC3k Part 4):**
  - Raw counts: mean non-zero `2.79`, max `419.0`
  - CP10k: mean non-zero `11.76`, max `1,751.1`
  - $\ln(1 + \text{CP10k})$: mean non-zero **`2.02`**, max **`7.47`** (stops a handful of hyper-abundant genes from hijacking Euclidean PCA distances!).
- **Sparsity invariance:** Because $\ln(1 + 0) = \mathbf{0}$, every exact zero stays zero! Unlike z-scoring (`sc.pp.scale`, which shifts zeros to $-\mu/\sigma$ and densifies the matrix), `CP10k` + `log1p` runs entirely in-place on CSC/CSR `.data` in $O(\text{nnz})$!
:::

---

# Algorithmic Complexity & Scaling of Single-Cell QC

How does `Scanpy` / `AnnData` process `737,280` droplets × `32,738` genes in ~3 seconds on a laptop?

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Stage | Naive Dense Algorithm | Sparse Single-Cell Algorithm | Complexity & PBMC3k Footprint |
| :--- | :--- | :--- | :--- |
| **Matrix storage & QC sums** | Dense `float32` array (`sum(axis=0/1)`) | `CSC` / `CSR` (`data`, `indices`, `indptr`) | $O(\text{nnz})$ vs. $O(n_{\text{obs}} n_{\text{vars}})$; **22.3 MB vs. 89.9 GB** |
| **Knee plot & MAD QC** | Repeated dense passes + sensitive $\sigma$ | Vectorized sort / linear-time selection (`np.median`) | $O(n_{\text{obs}} \log n_{\text{obs}})$ sort; **< 0.1 s** for 737k droplets |
| **Ambient testing (`EmptyDrops`)** | Pairwise full-matrix tests | Dirichlet–Multinomial Monte Carlo on transitional zone only | $O(\text{nnz}_{\text{trans}} \cdot B)$ over $\sim 10^4$ candidate droplets |
| **Doublet detection (`Scrublet`)** | Exact $O(N^2)$ pairwise distances | Randomized truncated SVD (PCA) + approximate $k$-NN (Annoy/PyNNDescent) | $O(\text{nnz} \cdot d + N \log N)$ — **~1 s** on 2,700 cells |
| **CP10k + `log1p` normalization** | Dense element-wise operations | In-place `.data` scaling + `np.log1p` preserving `0 → 0` | $O(\text{nnz})$ time & **zero extra RAM** |
:::

:::matrix { .cards cols="1fr 1fr" gap="22px" }
[[0,0]]{ .purple }
### Sparsity-Preserving Pipeline Design
Keeping every preprocessing step (`read_10x_mtx` $\to$ `calculate_qc_metrics` $\to$ `normalize_total` $\to$ `log1p`) strictly $O(\text{nnz})$ with `0 → 0` guarantees datasets scale cleanly from `3,000` PBMCs to **1,000,000-cell tissue atlases**.

[[0,1]]{ .green }
### When Sparsity Breaks (Preview of Lecture 6)
Mean-centering (`X - mean`) or z-scoring (`sc.pp.scale`) turns `0` into non-zero, instantly densifying the matrix! That is why Lecture 6 filters down to ~2,000 **Highly Variable Genes (HVGs)** *before* scaling and PCA.
:::

---

# Summary & Bridge to Part 4 (Industrial Tools)

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### Sparsity & Calling
#### `AnnData` · `scipy.sparse`
- Raw PBMC3k: `99.988%` zeros
- CSC (`22.3 MB`) vs. Dense (`89.9 GB`)
- Always flatten `.sum()` & `.copy()` views
- Knee plot ($m/10$) vs. `EmptyDrops`
- Ambient soup enriched $3\times$ in `MT-`

[[0,1]]{ .green }
### Robust QC & Doublets
#### `MAD` · `Scrublet`
- Breakdown: `0%` SD vs. `50%` MAD
- Consistency: $k = 1/\Phi^{-1}(0.75) \approx 1.4826$
- One-sided 3-MAD (`MT-%`: `4.43%`)
- Two-sided log-5-MAD (`counts`, `genes`)
- Poisson doublets ($\approx \lambda/2$) & $k$-NN bridges

[[0,2]]{ .blue }
### Dropout, UMIs & L6
#### `Occupancy` · `CP10k`
- Low $\mu$: Poisson $e^{-\mu}$ capture noise
- High $\mu$: NB excess zeros = biology
- UMI occupancy $E = N(1 - (1 - 1/N)^k)$
- 10-nt $\to$ 12-nt gives $16\times$ capacity
- `CP10k` + `log1p` preserves `0 → 0`
:::

> **Core Equations at a Glance:** $\text{CSC} = 8 \ \text{nnz} + 4(n_{\text{vars}} + 1)$ · $\text{MAD}\_{\text{scaled}} = 1.4826 \cdot \operatorname{median}\_i \lvert x_i - \operatorname{median}(x)\rvert$ · $P(X \ge 2 \mid X \ge 1) \approx \lambda/2$ · $P_{\text{NB}}(0) = (1+\alpha\mu)^{-1/\alpha}$ · $E = N(1 - (1 - 1/N)^k)$ · $\hat{k} = -N\ln(1 - E/N)$ · $\tilde{x}\_{cg} = \ln(1 + 10^4 x_{cg}/\sum_{g'} x_{cg'})$
