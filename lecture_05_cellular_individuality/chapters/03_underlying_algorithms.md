:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Sparse Matrices · Cell Calling · Robust MAD · Doublets · UMI Collisions
:::

---

# Curriculum

- **Sparse Matrix Representations & Compression** — why single-cell count tables exceed 97%–99.98% sparsity; MatrixMarket coordinate format (COO) vs. Compressed Sparse Column / Row (`CSC` / `CSR`) indexing (`data`, `indices`, `indptr`); $O(\text{nnz})$ memory and row/column summation complexity; view vs. copy memory semantics.
- **Droplet Cell Calling (Empty vs. Real Cells)** — barcode rank ("knee") log–log geometry separating the cellular plateau from the ambient RNA tail; order-of-magnitude cutoff (Cell Ranger v2: $m/10$) vs. fixed UMI cutoffs; statistical ambient-profile testing (`EmptyDrops` Dirichlet-multinomial test against the ambient soup to rescue low-RNA cells).
- **Robust Statistics & Dynamic Outlier Detection (MAD)** — finite-sample breakdown point ($0\%$ for mean/SD vs. $50\%$ for median/MAD) and outlier masking; derivation of the Gaussian consistency factor $k = 1 / \Phi^{-1}(0.75) \approx 1.4826$; one-sided vs. log-scale two-sided MAD filtering for mitochondrial leakage, library size, and gene diversity.
- **Microfluidic Encapsulation & Doublet Detection** — Poisson loading kinetics $P(X = k) = \lambda^k e^{-\lambda}/k!$, conditional multiplet rate among occupied droplets, homotypic vs. heterotypic doublets, and in silico doublet scoring (`Scrublet`: synthetic transcriptome summation, joint PCA, and $k$-NN doublet-neighbor fraction).
- **Sparsity, Dropout Kinetics & Normalization Preview** — Poisson sampling zero probability $P(X = 0) = e^{-\mu}$ vs. Negative Binomial overdispersion $(1 + \alpha\mu)^{-1/\alpha}$ (separating technical capture limits from biological heterogeneity and transcriptional bursting); library-size scaling (CP10k) and $\log(1+x)$ variance stabilization bridging to Lecture 6.
- **UMI Deduplication & Collision Occupancy Mathematics** — pre-PCR molecular barcoding to collapse PCR duplicates; classical occupancy (urn / birthday) model over $N = 4^L$ barcodes for $k$ molecules; expected distinct UMIs $E$, birthday collision probability, under-counting loss $k - E$, and inverse collision correction $\hat{k} = -N\ln(1 - E/N)$.

---
# Key Math Notes — Sparsity, Cell Calling & Robust MAD

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Sparse Storage Economics
- For $n_{\text{obs}}$ droplets, $n_{\text{vars}}$ genes, and $\text{nnz}$ non-zeros:
  - Dense `float32`: $4 \cdot n_{\text{obs}} \cdot n_{\text{vars}}$ bytes ($\sim 90\text{ GB}$ raw PBMC3k).
  - `CSC` (`data` + `indices` + `indptr`): $8 \cdot \text{nnz} + 4(n_{\text{vars}} + 1)$ bytes ($\sim 22\text{ MB}$).

### Droplet Cell Calling
- **Cell Ranger v2 heuristic** (expected cell count $N$):
  $$m = \text{percentile}_{99}(\text{top } N), \qquad T = \frac{m}{10}$$
- **Ambient profile test (`EmptyDrops`)**: estimate background distribution $\hat{\mathbf{p}}\_{\text{amb}}$ from low-UMI droplets; test each droplet's count vector $\mathbf{c}$ under $\text{Dirichlet-Multinomial}(\hat{\mathbf{p}}\_{\text{amb}})$ to reject the empty-droplet null hypothesis.

[[0,1]]
### Median Absolute Deviation (MAD)
$$
  \text{MAD}(x) = \operatorname{median}\_i \bigl\lvert x_i - \operatorname{median}(x) \bigr\rvert
$$
- **Breakdown point**: $\varepsilon^\*(\bar{x}, s) = 0$ vs. $\varepsilon^\*(\text{Med}, \text{MAD}) = 0.5$.
- **Normal consistency constant** ($Z \sim \mathcal{N}(0,1)$, $P(\lvert Z \rvert \le c) = 0.5$):
$$
  k = \frac{1}{\Phi^{-1}(0.75)} = \frac{1}{\sqrt{2}\\;\operatorname{erf}^{-1}(0.5)} \approx 1.4826022
$$
  $$\text{MAD}_{\text{scaled}} = k \cdot \text{MAD}(x)$$
- **Dynamic outlier rules**:
  - One-sided (`pct_counts_mt`): $x > \text{Med} + n_{\text{MAD}} \cdot \text{MAD}\_{\text{scaled}}$ ($n_{\text{MAD}} = 3$).
  - Log-scale two-sided ($y = \log(1 + x)$): $\lvert y - \text{Med}(y) \rvert > 5 \cdot \text{MAD}\_{\text{scaled}}(y)$.
:::

---
# Key Math Notes — Microfluidics, Dropout & UMI Collisions

:::matrix { cols=40/60 gap="30px" }
[[0,0]]
### Poisson Microfluidic Loading & Doublets
- Droplet cell count $X \sim \operatorname{Poisson}(\lambda)$:
  $$P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}$$
- **Conditional doublet/multiplet rate** among non-empty droplets ($\lambda \ll 1$):
$$
  P(X \ge 2 \mid X \ge 1) = \frac{1 - e^{-\lambda} - \lambda e^{-\lambda}}{1 - e^{-\lambda}} \approx \frac{\lambda}{2}
$$
- **Scrublet score**: local fraction of simulated doublets among a cell's $k$ nearest neighbors in PCA space.

### Zero Fraction & Library Normalization
- **Poisson capture baseline**: $P(X = 0) = e^{-\mu}$.
- **Negative Binomial (heterogeneity / bursting)**: $P(X = 0) = (1 + \alpha\mu)^{-1/\alpha} > e^{-\mu}$.
- **CP10k + log1p**: $\tilde{x}\_{cg} = \ln\left(1 + 10^4 \cdot \frac{x_{cg}}{\sum_{g'} x_{cg'}}\right)$.

[[0,1]]
### UMI Occupancy & Collision Correction
- Barcode space of length $L$: $N = 4^L$ distinct UMIs.
- Drawing $k$ independent molecules with replacement from $N$ bins:
  - Probability a given UMI stays empty: $(1 - 1/N)^k \approx e^{-k/N}$.
  - **Expected distinct UMIs observed**:
    $$E = N\left(1 - \left(1 - \frac{1}{N}\right)^k\right) \approx N\left(1 - e^{-k/N}\right)$$
  - **Expected lost counts**: $k - E$, fractional loss $(k - E)/k$.
  - **Birthday collision probability** ($\ge 1$ collision among $k$ molecules):
    $$P(\ge 1\text{ collision}) = 1 - \prod_{j=0}^{k-1}\left(1 - \frac{j}{N}\right) \approx 1 - \exp\left(-\frac{k(k-1)}{2N}\right)$$
  - **Inverse collision correction** (recovering true $k$ from observed $E < N$):
    $$\hat{k} = -N\ln\left(1 - \frac{E}{N}\right)$$
:::
