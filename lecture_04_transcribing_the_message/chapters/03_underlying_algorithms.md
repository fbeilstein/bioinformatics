:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Pseudoalignment · Negative Binomial GLMs · Empirical Bayes · FDR
:::

---

# Curriculum

- **Transcriptome Pseudoalignment & Equivalence-Class EM (`kallisto`)** — colored transcriptome de Bruijn graph ($k$-mer compatibility sets), equivalence-class intersection without base-level alignment, effective transcript length $\tilde{\ell}_t$, Expectation–Maximization (EM) for isoform deconvolution and fractional `est_counts`, bootstrap quantification uncertainty, and TPM normalization (and why TPM cannot replace raw counts for statistical testing).
- **Overdispersed Count Distributions (Poisson vs. Negative Binomial)** — Poisson shot noise in technical replicates ($\text{Var} = \mu$), Gamma–Poisson hierarchical mixture for biological replicates yielding $\text{NB}(\mu, \alpha)$ with quadratic variance $\text{Var} = \mu + \alpha\mu^2$, conversion to the `(n, p)` parameterization (`scipy.stats.nbinom`), and unbiased method-of-moments dispersion initialization.
- **Sequencing-Depth Normalization (Median of Ratios)** — compositional bias under extreme differential expression (why library-size / CPM scaling distorts unchanged genes), geometric-mean pseudo-reference and median-of-ratios size factors $s_j$, and keeping $s_j$ as a GLM log-offset rather than pre-dividing counts.
- **Negative Binomial Generalized Linear Models (GLMs)** — log-link regression $\log \mu_{ij} = \log s_j + \mathbf{x}_j^\top \boldsymbol{\beta}_i$, reference-level design encoding, conversion between natural-log coefficients and $\log_2\text{FC}$, and Iteratively Reweighted Least Squares (IRLS) with Fisher information standard errors.
- **Empirical Bayes Dispersion Shrinkage & Outlier Handling** — failure of gene-wise MLEs at small replicate count ($n = 3$), Cox–Reid adjusted profile likelihood, parametric mean–dispersion trend $\alpha_{\text{tr}}(\bar{\mu}) = a_0 + a_1 / \bar{\mu}$, log-normal prior and MAP shrinkage (why low-count flat likelihoods shrink strongest), dispersion-outlier protection, and Cook's distance sample-outlier flagging.
- **Hypothesis Testing, Multiplicity Control & LFC Shrinkage** — Wald test ($z = \hat{\beta}/\text{SE}$) vs. Likelihood Ratio Test (LRT), independent filtering on `baseMean` to maximize power without inflating Type I error, Benjamini–Hochberg False Discovery Rate (FDR / `padj`), and empirical Bayes LFC shrinkage (apeGLM-type priors) to stabilize noisy low-count fold-change rankings.

---

# Key Math Notes — Quantification, Normalization & Count Models

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Pseudoalignment & EM on Equivalence Classes
- Read $r$ matches k-mers $k_1, \dots, k_q$, each mapping to transcript set $C(k_m)$:
  $$E(r) = \bigcap_{m=1}^{q} C(k_m)$$
- Collapse reads into equivalence classes $c$ with counts $y_c$; maximize likelihood over **fragment probabilities** $\pi_t$ ($\sum_t \pi_t = 1$, expected counts $\hat{c}_t = N\pi_t$):
  $$\mathcal{L}(\boldsymbol{\pi}) \propto \prod_{c} \left(\sum_{t \in c} \frac{\pi_t}{\tilde{\ell}_t}\right)^{y_c}, \quad \tilde{\ell}_t \approx \ell_t - \bar{L}_{\text{frag}} + 1$$
- **EM update**: **E-step** splits class counts $y_c$ in proportion to $\pi_t / \tilde{\ell}_t$; **M-step** sums split counts $\hat{c}_t = \sum_{c \ni t} n_{c,t}$ and sets $\pi_t \leftarrow \hat{c}_t / N$ → fractional `est_counts`.
- **TPM**:
  $$\text{TPM}_t = \frac{\hat{c}_t / \tilde{\ell}_t}{\sum_{u} \hat{c}_u / \tilde{\ell}_u} \cdot 10^6$$

[[0,1]]
### Gamma–Poisson Mixture → Negative Binomial
- Biological rate $\Lambda_{ij} \sim \operatorname{Gamma}(1/\alpha_i, 1/(\alpha_i\mu_{ij}))$ and $K_{ij} \mid \Lambda_{ij} \sim \operatorname{Poisson}(\Lambda_{ij})$ gives $K_{ij} \sim \operatorname{NB}(\mu_{ij}, \alpha_i)$:
  $$\mathbb{E}[K_{ij}] = \mu_{ij}, \qquad \operatorname{Var}(K_{ij}) = \mu_{ij} + \alpha_i \mu_{ij}^2$$
- **Biological CV²**: $\text{Var}/\mu^2 = 1/\mu + \alpha_i \to \alpha_i$ as $\mu \to \infty$.
- **`scipy.stats.nbinom` parameterization**:
  $$n = \frac{1}{\alpha}, \qquad p = \frac{n}{n + \mu} \qquad (\alpha \to 0 \Rightarrow \text{Poisson})$$
- **Method-of-moments initializer** (unbiased $s^2$, `ddof=1`):
  $$\hat{\alpha}_{\text{MoM}} = \max\left(\frac{s^2 - \hat{\mu}}{\hat{\mu}^2}, \alpha_{\min}\right)$$
- **Median-of-ratios size factor** (genes with $\prod_v K_{iv} > 0$):
  $$s_j = \operatorname*{median}_i \frac{K_{ij}}{\left(\prod_{v=1}^{m} K_{iv}\right)^{1/m}}, \qquad \tilde{K}_{ij} = \frac{K_{ij}}{s_j}$$
:::

---
# Key Math Notes — NB GLMs, Empirical Bayes & Testing

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### NB GLM & Cox–Reid Dispersion Estimation
- **Log-link GLM** (offset $\log s_j$, design row $\mathbf{x}_j$):
$$
  \log \mu\_{ij} = \log s_j + \mathbf{x}\_j^\top \boldsymbol{\beta}\_i, \qquad \log_2 \text{FC} = \frac{\beta\_{i,1}}{\ln 2}
$$
- **IRLS covariance**: $\operatorname{Cov}(\hat{\boldsymbol{\beta}}\_i) = (\mathbf{X}^\top \mathbf{W}\_i \mathbf{X})^{-1}$ with weights $w_{ij} = (\mu_{ij}^{-1} + \alpha_i)^{-1}$.
- **Cox–Reid adjusted profile likelihood** (correcting bias from $p$ fitted coefficients):
$$
  \ell_{\text{CR}}(\alpha_i) = \ell(\alpha_i; \hat{\boldsymbol{\beta}}\_i) - \tfrac{1}{2}\log\det(\mathbf{X}^\top \mathbf{W}\_i \mathbf{X})
$$

### Empirical Bayes Dispersion Shrinkage
- **Parametric trend**: $\alpha_{\text{tr}}(\bar{\mu}) = a_0 + a_1 / \bar{\mu}$.
- **Log-normal prior** & **MAP estimate**:
$$
  \hat{\alpha}\_i^{\text{MAP}} = \arg\max_{\alpha} \left[\ell_{\text{CR}}(\alpha) - \frac{\bigl(\log \alpha - \log \alpha_{\text{tr}}(\bar{\mu}_i)\bigr)^2}{2\sigma_d^2}\right]
$$
- Low-count genes have flat $\ell_{\text{CR}}(\alpha)$ → prior pulls $\hat{\alpha}\_i^{\text{MAP}}$ close to $\alpha_{\text{tr}}$; high-dispersion outliers are not shrunk downward.

[[0,1]]
### Outlier Detection, Wald Test & Multiplicity
- **Cook's distance** $D_{ij}$ measures influence of sample $j$ on $\hat{\boldsymbol{\beta}}\_i$; flags extreme single-sample outliers ($D_{ij} > F_{0.99}(p, m - p)$).
- **Wald statistic** for $H_0: \beta_{i,1} = 0$:
$$
  z_i = \frac{\hat{\beta}\_{i,1}}{\operatorname{SE}(\hat{\beta}_{i,1})}, \qquad p_i = 2\bigl(1 - \Phi(\lvert z_i \rvert)\bigr)
$$
- **Independent filtering**: choose `baseMean` cutoff $\tau$ maximizing discoveries; set `padj = NaN` below $\tau$ (valid because $\bar{\mu}_i$ is independent of condition labels under $H_0$).
- **Benjamini–Hochberg FDR** (sorted $p_{(1)} \le \dots \le p_{(m)}$):
$$
  p_{(k)}^{\text{adj}} = \min_{l \ge k} \frac{m p_{(l)}}{l}, \qquad \text{FDR} = \mathbb{E}\left[\frac{V}{\max(R, 1)}\right]
$$
- **LFC shrinkage (apeGLM)**: zero-centered heavy-tailed prior on $\beta_{i,1}$ shrinks high-$\text{SE}$ (low-count) LFCs toward $0$ for ranking/MA plots while keeping Wald $p$-values unchanged.
:::
