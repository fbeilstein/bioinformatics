:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Pseudoalignment · Negative Binomial GLMs · Empirical Bayes · FDR
:::

---

# Curriculum & Algorithmic Roadmap

Six interconnected statistical algorithms turn raw bulk RNA-seq FASTQ reads into calibrated differential expression tables in Part 4 (`kallisto` + `PyDESeq2` on yeast heat-shock `GSE135568`):

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "fq",   "label": "Raw FASTQ\n(10⁷ cDNA reads)",     "x": 20,   "y": 40, "width": 150, "height": 64 },
    { "id": "kal",  "label": "Colored dBG +\nEquiv-Class EM",     "x": 210,  "y": 40, "width": 165, "height": 64, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "mor",  "label": "Median-of-Ratios\nSize Factors s_j","x": 415,  "y": 40, "width": 165, "height": 64, "style": "fill: #2563EB; color: #fff;" },
    { "id": "nb",   "label": "Gamma–Poisson\nNB(μ, α) Model",     "x": 620,  "y": 40, "width": 160, "height": 64, "style": "fill: #0891B2; color: #fff;" },
    { "id": "eb",   "label": "Cox–Reid + EB\nDispersion MAP",     "x": 820,  "y": 40, "width": 165, "height": 64, "style": "fill: #16A34A; color: #fff;" },
    { "id": "glm",  "label": "NB GLM (IRLS)\n+ Cook's Outlier",   "x": 1025, "y": 40, "width": 165, "height": 64, "style": "fill: #D97706; color: #fff;" },
    { "id": "fdr",  "label": "Wald/LRT · Filter\nBH FDR · apeGLM","x": 1230, "y": 40, "width": 170, "height": 64, "style": "fill: #DC2626; color: #fff;" }
  ],
  "links": [
    { "from": "fq",  "to": "kal" },
    { "from": "kal", "to": "mor" },
    { "from": "mor", "to": "nb"  },
    { "from": "nb",  "to": "eb"  },
    { "from": "eb",  "to": "glm" },
    { "from": "glm", "to": "fdr" }
  ]
}
:::

:::matrix { .cards cols="1fr 1fr 1fr" gap="18px" }
[[0,0]]{ .purple }
### 1. Pseudoalignment & EM
Colored transcriptome de Bruijn graph, compatibility set intersections, effective length $\tilde{\ell}_t$, fractional `est_counts`, and why TPM cannot replace raw counts.

[[0,1]]{ .blue }
### 2. Count Models & Size Factors
Gamma–Poisson mixture $\to \operatorname{NB}(\mu, \alpha)$, `scipy.stats.nbinom` $(n, p)$, compositional bias, and median-of-ratios GLM offsets.

[[0,2]]{ .green }
### 3. EB Shrinkage, GLMs & FDR
IRLS log-link GLMs, Cox–Reid likelihood, Empirical Bayes MAP dispersion shrinkage, Cook's distance, Wald/LRT, independent filtering, BH FDR, and apeGLM LFC shrinkage.
:::

---

# Pseudoalignment & Colored de Bruijn Graphs (`kallisto`)

Variant calling (Lecture 3) needed **exact base coordinates** (CIGAR, MAPQ). Expression quantification needs only **which transcripts a read is compatible with** — enabling $100\times$ faster alignment-free counting:

:::static-diagram { width=92% center }
{
  "nodes": [
    { "id": "r",   "label": "Read r\n(k = 31 bp k-mers)",             "x": 20,  "y": 35, "width": 175, "height": 62 },
    { "id": "h",   "label": "Hash Lookup into\nTranscriptome dBG",    "x": 240, "y": 35, "width": 185, "height": 62, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "c",   "label": "k-mer Color Sets\nC(k₁)={T1,T2,T3}\nC(k₂)={T1,T2}", "x": 470, "y": 25, "width": 195, "height": 82, "style": "fill: #2563EB; color: #fff;" },
    { "id": "i",   "label": "Set Intersection\nE(r) = ∩ C(k_m)\n= {T1, T2}", "x": 710, "y": 25, "width": 180, "height": 82, "style": "fill: #16A34A; color: #fff;" },
    { "id": "ec",  "label": "Equivalence Class\nIncrement y_c += 1",  "x": 935, "y": 35, "width": 185, "height": 62, "style": "fill: #D97706; color: #fff;" }
  ],
  "links": [
    { "from": "r", "to": "h" },
    { "from": "h", "to": "c" },
    { "from": "c", "to": "i" },
    { "from": "i", "to": "ec" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Colored Transcriptome de Bruijn Graph
- Index **cDNA transcripts** (not genomic DNA) as overlapping $k$-mers ($k = 31$, just like Mash in Lecture 1).
- Each graph contig stores a **color set** $C(k_m) \subseteq \{T_1, \dots, T_M\}$ of all isoforms containing $k_m$.
- Read $r$ defines a compatibility set by intersecting its hit $k$-mers:
$$
E(r) = \bigcap_{m=1}^{q} C(k_m)
$$

[[0,1]]
### Why It Is $100\times$ Faster Than STAR / BWA
- **Contig skipping:** Once a $k$-mer lands on a non-branching contig, all interior $k$-mers share the exact same color set — jump directly to the contig end!
- **Sufficient statistics:** $10^7$ raw reads collapse into a few thousand **equivalence classes** $c$ with integer counts $y_c$ — zero Smith–Waterman DP matrices and zero BAM files!
- **Error resilience:** A 1-bp sequencing error misses the hash table ($\emptyset$) and is skipped while flanking error-free $k$-mers still resolve $E(r)$.
:::

---
# Equivalence-Class EM, Effective Length & the TPM Trap

![Equivalence-Class EM Isoform Deconvolution](lecture_04_transcribing_the_message/images/em_isoform_deconvolution.png){ width=74 center }

:::matrix { cols="1.15fr 0.85fr" gap="26px" }
[[0,0]]
- **Likelihood over equivalence classes** ($N = \sum_c y_c$, fragment probabilities $\pi_t$):
$$
\mathcal{L}(\boldsymbol{\pi}) \propto \prod_{c} \left(\sum_{t \in c} \frac{\pi_t}{\tilde{\ell}\_t}\right)^{y_c}, \qquad \tilde{\ell}\_t \approx \ell_t - \bar{L}\_{\text{frag}} + 1
$$
- **EM updates:** **E-step** splits class count $y_c$ via $n_{c,t} = y_c \frac{\pi_t / \tilde{\ell}\_t}{\sum_{u \in c} \pi_u / \tilde{\ell}\_u}$; **M-step** sums $\hat{c}\_t = \sum_{c \ni t} n_{c,t}$ and sets $\pi_t \leftarrow \hat{c}\_t / N$ $\implies$ **fractional `est_counts`** (bootstrapping `-b 100` resamples $y_c$ to gauge isoform uncertainty).

[[0,1]]
- **Transcripts Per Million (TPM):**
$$
\text{TPM}\_t = \frac{\hat{c}\_t / \tilde{\ell}\_t}{\sum_u \hat{c}\_u / \tilde{\ell}\_u} \cdot 10^6
$$
- **⚠ TPM Trap:** `TPM = 100` can arise from **3 reads** on a 300 bp gene ($\text{CV} = 58\\%$) or **3,000 reads** on a 300 kb gene ($\text{CV} = 1.8\\%$). TPM erases count certainty — **never** feed TPM into DESeq2!
:::

---

# Interactive Lab — Pseudoalignment & Equivalence-Class EM

Explore `kallisto`'s core engine live: **(1)** click reads across exons and splice junctions in **Tab 1** to trace $E(r) = \bigcap C(k_m)$ (including the 1-bp error $k$-mer case); **(2)** switch to **Tab 2**, step EM from uniform $\pi_t = 1/3$ to convergence, and watch how changing unique `{T1}` reads or fragment length $\bar{L}$ shifts fractional `est_counts` and `TPM`!

<iframe src="lecture_04_transcribing_the_message/demo/pseudoalignment_em_lab.html" title="Interactive Transcriptome Pseudoalignment and Equivalence-Class EM Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# From Poisson Shot Noise to Biological Overdispersion

Why can we model **technical replicates** (re-sequencing the same library) with a Poisson distribution, yet **biological replicates** (independent yeast cultures) wildly violate $\operatorname{Var} = \mu$?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Technical Shot Noise (Poisson)
- Capturing $K$ reads from a fixed cDNA pool of tiny per-gene share $p_i \ll 1$ across $N \gg 1$ reads follows a **Poisson** law:
$$
K_{ij} \mid \Lambda_{ij} \sim \operatorname{Poisson}(\Lambda_{ij}), \qquad \operatorname{Var}(K \mid \Lambda) = \Lambda
$$
- Relative sampling noise shrinks as count grows:
$$
\operatorname{CV}^2_{\text{tech}} = \frac{\operatorname{Var}(K)}{\mu^2} = \frac{1}{\mu} \xrightarrow{\mu \to \infty} 0
$$
- At $\mu = 2{,}000$ reads, technical noise is tiny ($\text{CV} \approx 2.2\\%$).

[[0,1]]
### Gamma–Poisson Mixture $\implies$ Negative Binomial
- Between independent cultures, the true biological rate $\Lambda_{ij}$ fluctuates (promoter kinetics, TF bursts, temperature micro-gradients):
$$
\Lambda_{ij} \sim \operatorname{Gamma}\left(\text{shape} = \frac{1}{\alpha_i}, \ \text{rate} = \frac{1}{\alpha_i \mu_{ij}}\right)
$$
- Marginalizing out $\Lambda_{ij}$ yields the **Negative Binomial** $K_{ij} \sim \operatorname{NB}(\mu_{ij}, \alpha_i)$ with exact Law of Total Variance:
$$
\mathbb{E}[K_{ij}] = \mu_{ij}, \qquad \operatorname{Var}(K_{ij}) = \underbrace{\mu_{ij}}\_{\text{shot noise}} + \underbrace{\alpha_i \mu_{ij}^2}\_{\text{bio variance}}
$$
- **Total squared CV:** $\operatorname{CV}^2 = \frac{1}{\mu_{ij}} + \alpha_i \to \alpha_i$ (asymptotic biological $\text{CV} = \sqrt{\alpha_i}$; in yeast `GSE135568`, $\alpha \approx 0.014 \implies 12\\%$ culture-to-culture CV).
:::

---
# Negative Binomial & `scipy.stats.nbinom` Parameterization

![Gamma-Poisson Overdispersion and Poisson vs NB PMF](lecture_04_transcribing_the_message/images/gamma_poisson_overdispersion.png){ width=90 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **`scipy.stats.nbinom(n, p)` conversion** (Practice 4 Homework Task 3):
$$
n = \frac{1}{\alpha}, \qquad p = \frac{n}{n + \mu} \qquad (\alpha \to 0 \implies \text{Poisson})
$$
- Check `nbinom.mean(n, p) == mu` and `nbinom.var(n, p) == mu + alpha*mu**2` (swapping $p \leftrightarrow 1 - p$ silently corrupts results!).

[[0,1]]
- **Unbiased Method-of-Moments initializer** (use `ddof=1`! NumPy's `ddof=0` underestimates $s^2$ by $1/3$ at $n = 3$):
$$
\hat{\alpha}\_{\text{MoM}} = \max\left(\frac{s^2 - \hat{\mu}}{\hat{\mu}^2}, \ \alpha_{\min}\right)
$$
- At $\mu = 500, \alpha = 0.02$, a normal replicate count $K = 650$ has $P_{\text{Pois}} = 8.0\times 10^{-11}$ (fake hit!) vs. $P_{\text{NB}} = 0.029$.
:::

---
# Compositional Bias & Median-of-Ratios Normalization

![Compositional Bias and Median-of-Ratios Size Factors](lecture_04_transcribing_the_message/images/compositional_bias_size_factors.png){ width=93 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Why CPM / total-sum scaling fails:** Sequencers deliver a **fixed budget** of reads. Massive induction of chaperones (or global shut-down of ribosomal mRNAs) alters $\sum_i K_{ij}$, shifting the apparent count of *every unchanged gene*.
- **DESeq2 Median-of-Ratios** builds a geometric-mean pseudo-reference $g_i$ across all genes with $\prod_v K_{iv} > 0$:

[[0,1]]
$$
g_i = \left(\prod_{v=1}^{m} K_{iv}\right)^{1/m}, \qquad s_j = \operatorname*{median}\_i \frac{K_{ij}}{g_i}
$$
- **Robust breakdown:** Taking the **median** ignores extreme differential expression in up to $50\\%$ of genes — explaining `GSE135568`'s `HS_1` reversal (`21.1 M < WT_2 22.9 M`, yet $s_{\text{HS1}} = 1.142 > 0.937$).
:::

---

# Negative Binomial GLMs: Log Link, Offsets & IRLS

Why does DESeq2 put $s_j$ **inside the GLM** rather than pre-dividing counts $\tilde{K}_{ij} = K_{ij} / s_j$ and running regression on normalized numbers?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Log-Link GLM with Fixed Offset
- Model raw integer counts $K_{ij} \sim \operatorname{NB}(\mu_{ij}, \alpha_i)$ with:
$$
\log \mu\_{ij} = \underbrace{\log s_j}_{\text{fixed offset}} + \mathbf{x}\_j^\top \boldsymbol{\beta}\_i
$$
- For two conditions (`WT` $\to x_j = 0$, `HeatShock` $\to x_j = 1$):
  - $\beta_{i,0} = \ln q_{i,\text{WT}}$ (baseline log-expression)
  - $\beta_{i,1} = \ln(q_{i,\text{HS}} / q_{i,\text{WT}}) \implies \log_2\text{FC} = \dfrac{\beta_{i,1}}{\ln 2}$
- **Offset preserves read precision:** $K = 2$ with $s = 0.2$ ($\tilde{K} = 10$) has $5\times$ higher Poisson variance than $K = 20$ with $s = 2.0$ ($\tilde{K} = 10$). Pre-dividing destroys that information!

[[0,1]]
### IRLS Estimation & Reference-Level Pitfall
- **Iteratively Reweighted Least Squares (IRLS)** solves the score equations $\mathbf{X}^\top \mathbf{W}\_i \mathbf{z}\_i = \mathbf{X}^\top \mathbf{W}\_i \mathbf{X} \boldsymbol{\beta}\_i$ using per-sample Fisher weights:
$$
w_{ij} = \left(\frac{1}{\mu_{ij}} + \alpha_i\right)^{-1}
$$
- **Asymptotic coefficient covariance** (gives `lfcSE`):
$$
\operatorname{Cov}(\hat{\boldsymbol{\beta}}\_i) = (\mathbf{X}^\top \mathbf{W}\_i \mathbf{X})^{-1}, \qquad \operatorname{SE}(\log_2\text{FC}) = \frac{\sqrt{[(\mathbf{X}^\top \mathbf{W}\_i \mathbf{X})^{-1}]\_{11}}}{\ln 2}
$$
- **⚠ Reference encoding:** Alphabetically `"HeatShock" < "WT"` flips the sign of $\beta_{i,1}$ unless you set explicit categories or contrasts!
:::

---

# Why Gene-Wise Dispersion Fails at $n = 3$ & Cox–Reid Likelihood

With only $n = 3$ biological replicates per group ($m = 6$ total samples, $p = 2$ GLM parameters $\implies m - p = 4$ residual degrees of freedom), estimating $\alpha_i$ for each gene alone is catastrophically unstable.

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### The 3-Replicate Variance Catastrophe
- A variance estimated from 3 numbers follows a wide $\chi^2_2$ distribution: in `GSE135568`, raw gene-wise MLEs span **9 orders of magnitude** ($10^{-8}$ to $10$)!
- **False-positive factory:** When 3 replicates happen by luck to cluster tightly, $\hat{\alpha}_i \approx 0$ shrinks $\operatorname{SE}(\hat{\beta}_{i,1})$ artificially $\implies$ inflated Wald $z$ $\implies$ spurious significance!
- **Downward MLE bias:** Because $\hat{\boldsymbol{\beta}}_i$ is fitted to the *same* $m$ data points, residuals are systematically too small (the GLM analog of `ddof=0` vs. `ddof=1`).

[[0,1]]
### Cox–Reid Adjusted Profile Likelihood
- DESeq2 corrects degrees-of-freedom bias by subtracting a Fisher-information curvature penalty from the profile NB log-likelihood:
$$
\ell_{\text{CR}}(\alpha_i) = \ell(\alpha_i; \hat{\boldsymbol{\beta}}\_i) - \frac{1}{2} \log\det\bigl(\mathbf{X}^\top \mathbf{W}\_i(\alpha_i)\mathbf{X}\bigr)
$$
- For an orthogonal 2-group design, $\det(\mathbf{X}^\top \mathbf{W}\_i \mathbf{X}) = \bigl(\sum_{j \in \text{WT}} w_{ij}\bigr)\bigl(\sum_{j \in \text{HS}} w_{ij}\bigr)$.
- Yet even with Cox–Reid bias correction, **4 degrees of freedom** per gene remain far too noisy — we must **borrow strength across all 7,126 genes**.
:::

---

# Empirical Bayes Dispersion Shrinkage & Outlier Protection

![Empirical Bayes Dispersion Shrinkage and Cox-Reid Posterior Geometry](lecture_04_transcribing_the_message/images/eb_dispersion_shrinkage.png){ width=77 center }

:::matrix { cols="1.1fr 0.9fr" gap="26px" }
[[0,0]]
- **Three-stage Empirical Bayes (Love et al., 2014):**
  1. Fit parametric trend $\alpha_{\text{tr}}(\bar{\mu}) = a_0 + a_1 / \bar{\mu}$ (`GSE135568`: $0.0137 + 3.03/\bar{\mu}$).
  2. Estimate true log-dispersion prior variance $\sigma_d^2$ around the trend (`0.51` in yeast).
  3. Maximize the penalized **MAP posterior** per gene:
$$
\hat{\alpha}\_i^{\text{MAP}} = \arg\max_{\alpha} \left[\ell_{\text{CR}}(\alpha) - \frac{\bigl(\log \alpha - \log \alpha_{\text{tr}}(\bar{\mu}_i)\bigr)^2}{2\sigma_d^2}\right]
$$

[[0,1]]
- **Why low counts shrink strongest:** At low $\bar{\mu}$, Poisson noise ($\mu$) swamps $\alpha\mu^2$, making $\ell_{\text{CR}}(\alpha)$ **flat** $\implies$ the prior pulls $\hat{\alpha}\_i^{\text{MAP}}$ tightly to $\alpha_{\text{tr}}$. At high $\bar{\mu}$, $\ell_{\text{CR}}(\alpha)$ is **sharp** $\implies$ gene data dominates.
- **Upward outlier rule:** If $\log \hat{\alpha}\_i^{\text{MLE}} > \log \alpha_{\text{tr}}(\bar{\mu}_i) + 2\sigma_d$, keep $\hat{\alpha}\_i^{\text{MLE}}$ unshrunk (`94` genes in yeast) — never shrink a truly hypervariable gene into a false positive!
:::

---

# Interactive Lab — Gamma–Poisson NB & Empirical Bayes Shrinkage

Test both count-level and gene-level models live: **(1)** in **Tab 1**, inspect `Part 4 Worked` ($\mu = 500, \alpha = 0.02, k^* = 650$) to see $P_{\text{Pois}} = 8.0\times 10^{-11}$ vs. $P_{\text{NB}} = 0.029$ and check `scipy`'s $(n, p)$; **(2)** in **Tab 2**, load `YDR130C`, `HSP104`, `Low-Count Trap`, and `Dispersion Outlier`, then increase replicates $n$ from `3` to `12` to watch the likelihood sharpen!

<iframe src="lecture_04_transcribing_the_message/demo/nb_dispersion_eb_lab.html" title="Interactive Negative Binomial and Empirical Bayes Dispersion Shrinkage Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Sample Outliers & Cook's Distance

Even after empirical Bayes dispersion shrinkage, a **single corrupted sample** (e.g., PCR jackpot, RNA degradation, or sample swap) can hijack $\hat{\beta}_{i,1}$ and fabricate an extreme fold change.

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Leave-One-Out Influence in GLMs
- **Cook's distance** $D_{ij}$ quantifies how much deleting sample $j$ shifts the fitted coefficient vector $\hat{\boldsymbol{\beta}}\_i$ (scaled by its Fisher covariance):
$$
D_{ij} = \frac{R_{ij}^2}{p} \cdot \frac{h_{jj}}{(1 - h_{jj})^2}
$$
  where $R_{ij} = \dfrac{K_{ij} - \hat{\mu}\_{ij}}{\sqrt{\hat{\mu}\_{ij} + \hat{\alpha}\_i \hat{\mu}\_{ij}^2}}$ is the standardized Pearson residual and $h_{jj}$ is the **hat-matrix leverage** diagonal of $\mathbf{W}\_i^{1/2}\mathbf{X}(\mathbf{X}^\top \mathbf{W}\_i \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{W}\_i^{1/2}$.

[[0,1]]
### Flagging vs. Replacement Policy (`PyDESeq2`)
- **Flag threshold:** A gene is flagged if $\max_j D_{ij} > F_{0.99}(p, m - p)$ (the 99th percentile of the $F$-distribution).
- **Small cohorts ($3 \le n_{\text{rep}} < 7$):** Too few replicates exist to safely impute the outlier $\implies$ DESeq2 sets **`pvalue = NaN`, `padj = NaN`** (`2` genes in `GSE135568`, protecting the top-hit list).
- **Large cohorts ($n_{\text{rep}} \ge 7$):** Replace the offending count with the trimmed group mean and refit (`Replacing 0 outlier genes` in Part 4's $n = 3$ log!).
:::

---

# Hypothesis Testing: Wald Test vs. Likelihood Ratio Test (LRT)

Once $\hat{\alpha}\_i^{\text{MAP}}$ and $\hat{\boldsymbol{\beta}}\_i$ are fitted, how do we compute a $p$-value for $H_0: \beta_{i,1} = 0$ — and when should we use Wald vs. LRT?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Wald Test (Pairwise Contrasts — Default)
- Tests a single linear contrast $\mathbf{c}^\top \boldsymbol{\beta}\_i = 0$ (e.g., `HeatShock vs. WT`) using the IRLS standard error:
$$
z_i = \frac{\hat{\beta}\_{i,1}}{\operatorname{SE}(\hat{\beta}\_{i,1})}, \qquad p_i = 2\bigl(1 - \Phi(\lvert z_i \rvert)\bigr)
$$
- **Fast:** Fits only **one** GLM per gene; any contrast is extracted immediately from $\hat{\boldsymbol{\beta}}\_i$ and $(\mathbf{X}^\top \mathbf{W}\_i \mathbf{X})^{-1}$.
- **⚠ Float64 underflow:** For `SSA4` ($z = 47.1$) and `HSP12` ($z = 46.5$), exact $p \approx 10^{-470} < 2.2 \times 10^{-308}$ underflows to `0.0` $\implies$ always clip `padj.clip(lower=1e-300)` before $-\log_{10}(\text{padj})$!

[[0,1]]
### Likelihood Ratio Test (Multi-Level Designs)
- Compares a **full model** (e.g., `~ time + genotype`, $p_{\text{full}}$ params) against a **reduced model** (`~ genotype`, $p_{\text{red}}$ params):
$$
\Lambda_i = 2\bigl[\ell_{\text{full}}(\hat{\boldsymbol{\beta}}\_i^{\text{full}}) - \ell_{\text{red}}(\hat{\boldsymbol{\beta}}\_i^{\text{red}})\bigr] \xrightarrow{H_0} \chi^2_{\Delta p}
$$
  with $\Delta p = p_{\text{full}} - p_{\text{red}}$ degrees of freedom.
- **When to use LRT:** Testing whether a gene varies across **any** of 3+ timepoints (e.g., `GSE135568`'s full `0 / 20 / 120 min` series) or interaction terms in one unified test.
:::

---

# Multiple Testing: FWER vs. Benjamini–Hochberg FDR

Testing $m \approx 6{,}400$ expressed yeast genes simultaneously breaks classical $p < 0.05$: under the global null, $p_i \sim \operatorname{Uniform}(0, 1)$ still produces **$0.05 \times 6{,}400 \approx 320$ false discoveries** by pure chance!

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### FWER (Bonferroni) vs. FDR (Benjamini–Hochberg)
- **Family-Wise Error Rate (Bonferroni):** Controls $P(V \ge 1) \le \alpha$ via $p_i \le \alpha / m \approx 7.8 \times 10^{-6}$ — far too conservative for genome-wide screens (throws away hundreds of real biology hits).
- **False Discovery Rate (BH, 1995):** Controls the expected **proportion** of false positives among rejected genes:
$$
\text{FDR} = \mathbb{E}\left[\frac{V}{\max(R, 1)}\right] \le q
$$
- **Meaning of `padj < 0.05`:** Among `2,425` called yeast genes, the expected share of false positives is $\le 5\\%$ ($\lesssim 121$ genes) — **not** "each individual gene has a 95% probability of being real."

[[0,1]]
### Exact BH Step-Up Algorithm (`padj`)
1. Sort the $m$ tested $p$-values in ascending order:
$$
p_{(1)} \le p_{(2)} \le \dots \le p_{(m)}
$$
2. Find the largest rank $k^\*$ satisfying $p_{(k)} \le \dfrac{k}{m} q$, and reject $1, \dots, k^*$.
3. Equivalently, compute monotonic **BH-adjusted $p$-values** right-to-left:
$$
p^{\text{adj}}\_{(k)} = \min_{l \ge k} \left(\min\left(1, \ \frac{m \ p_{(l)}}{l}\right)\right)
$$
   and call genes with $p^{\text{adj}} < q$.
:::

---

# Independent Filtering & apeGLM LFC Shrinkage

![Benjamini-Hochberg Independent Filtering and apeGLM MA Funnel](lecture_04_transcribing_the_message/images/fdr_and_lfc_shrinkage.png){ width=94% center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Independent filtering (`baseMean` $\ge \tau$):** Ultra-low-count genes have near-zero power to reach significance yet inflate $m$ in $m p_{(l)} / l$.
- Choosing $\tau$ on `baseMean` (which is **independent** of condition labels under $H_0$, Bourgon et al. 2010) shrinks $m$, tilts the BH line steeper, and **increases true discoveries** (`126` filtered genes get `padj = NaN` in `GSE135568`).

[[0,1]]
- **apeGLM LFC shrinkage (Zhu et al., 2019):** Raw MLE fold changes flare wildly at low counts (`YIL100W`: `baseMean = 5`, `LFC = 1.16`, `SE = 1.07`).
- Placing a zero-centered **heavy-tailed Cauchy prior** $\beta \sim \operatorname{Cauchy}(0, s_0)$ shrinks noisy low-count LFCs (`1.16 → 0.21`) while leaving strong, high-count hits (`HSP12`: `9.12 → 9.11`) untouched — **for ranking/MA plots only** (Wald $p$-values stay unchanged).
:::

---

# Interactive Lab — Size Factors, BH FDR & apeGLM Shrinkage

Experiment with normalization and multiplicity live: **(1)** in **Tab 1**, drag **Chaperone induction** from `1×` to `60×` and watch CPM drag every unchanged housekeeping gene below zero while Median-of-Ratios stays locked at `0.00` (then lower **Ribosomal repression** to reproduce `GSE135568`'s `HS_1` paradox!); **(2)** in **Tab 2**, adjust `baseMean` cutoff $\tau$ to maximize BH discoveries and toggle **apeGLM Shrinkage** on the MA funnel!

<iframe src="lecture_04_transcribing_the_message/demo/glm_fdr_shrinkage_lab.html" title="Interactive Compositional Normalization, BH FDR and apeGLM Shrinkage Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Diagnostic Geometry: Reading Mean–Variance, Volcano & MA Plots

Every standard RNA-seq plot in Practice 4 is a direct visual audit of one mathematical assumption derived above:

:::matrix { .cards cols="1fr 1fr 1fr" gap="22px" }
[[0,0]]{ .purple }
### Mean–Variance Plot (Part 3)
- **Axes:** $\log_{10}(\hat{\mu})$ vs. $\log_{10}(s^2)$ across replicates of one condition.
- **What it verifies:** Points follow $y = x$ (slope 1) at low counts and bend smoothly upward to **slope 2** ($\alpha \mu^2$) at high counts.
- **Failure signal:** Points far *below* $y = x$ indicate pre-normalized/non-count input; massive scatter above the trend flags batch effects or sample swaps.

[[0,1]]{ .green }
### Volcano Plot (Part 6)
- **Axes:** `log2FoldChange` ($x$) vs. $-\log_{10}(\text{padj})$ ($y$).
- **What it verifies:** Joint **effect size** ($\hat{\beta}_1$) and **statistical evidence** ($z = \hat{\beta}_1 / \text{SE}$).
- **Failure signal:** An asymmetric right plume (`HSP12`, `SSA4` at $z \approx 47$) reflects true biology, whereas a **horizontally shifted baseline** indicates failed normalization; always clip `padj` at `1e-300` so `0.0` does not vanish to $+\infty$.

[[0,2]]{ .blue }
### MA Plot (Homework Task 2)
- **Axes:** $\log_{10}(\text{baseMean})$ ($x$) vs. `log2FoldChange` ($y$).
- **What it verifies:**
  1. Cloud centered at **$y = 0$** (confirms size factors removed compositional bias).
  2. Wide left funnel narrowing rightward ($\text{SE} \propto \sqrt{1/\mu + \alpha}$); collapses after apeGLM shrinkage.
- **Zero guard:** `593` all-zero genes give $\log_{10}(0) = -\infty$ — filter `baseMean > 0` or plot `log10(baseMean + 1)`.
:::

---

# Algorithmic Complexity & Scaling of Bulk RNA-Seq Pipelines

Why can a complete 6-sample, $130\text{ M}$-read yeast differential expression study run end-to-end on a laptop in minutes?

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Stage | Classic Approach | Modern Bulk RNA-Seq Algorithm | Complexity & Scaling |
| :--- | :--- | :--- | :--- |
| **Read quantification** | Spliced genomic alignment (`STAR`/`HISAT2`) + BAM + `featureCounts` | `kallisto` colored dBG pseudoalignment + equivalence-class EM | $O(\lvert R \rvert \cdot L/k)$ hash lookups + $O(\lvert C \rvert \cdot T_{\text{EM}})$ EM ($\lvert C \rvert \ll \lvert R \rvert$); **minutes, ~4 GB RAM** |
| **Depth normalization** | Total-sum CPM / RPKM (compositional bias) | Geometric-mean pseudo-reference + median-of-ratios | $O(G \cdot m \log G)$ vectorized matrix ops; **< 0.05 s** for 7,126 × 6 |
| **Dispersion + EB MAP** | Gene-wise MLE (unstable at $n = 3$) | MoM init $\to$ Cox–Reid profile $\ell_{\text{CR}}(\alpha)$ $\to$ trend + MAP | $O(G \cdot m \cdot T_{\text{opt}})$ embarrassingly parallel across genes $G$ |
| **GLM + Wald / LRT** | Gaussian $t$-test on $\log(\text{CPM} + 1)$ (breaks at low counts) | NB log-link IRLS + Fisher covariance + Cook's + BH + apeGLM | $O(G \cdot m \cdot p^2 \cdot T_{\text{IRLS}} + G \log G)$ — **~2 s** in `PyDESeq2` |
:::

:::matrix { .cards cols="1fr 1fr" gap="22px" }
[[0,0]]{ .purple }
### Key Computational Insight 1: Equivalence Compression
Quantification cost decouples from sequencing depth: once $10^7$ reads map into $\sim 10^4$ equivalence classes, running 100 EM iterations (and 100 bootstraps!) operates only on the tiny class count vector $y_c$, not individual reads.

[[0,1]]{ .green }
### Key Computational Insight 2: Gene-Parallel Statistical Engine
After size factors and the global trend $\alpha_{\text{tr}}(\bar{\mu})$ are computed across all $G$ genes, every gene's Cox–Reid MAP dispersion, NB GLM IRLS fit, Cook's distance, and Wald test is completely independent — scaling linearly across CPU cores.
:::

---

# Summary & Bridge to Part 4 (Industrial Tools)

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### Quantification
#### `kallisto`
- Colored transcriptome dBG ($k = 31$)
- $E(r) = \bigcap_m C(k_m)$, zero DP
- Equivalence-class EM $\to$ `est_counts`
- Effective length $\tilde{\ell}_t$ & `tpm`
- **Rule:** Feed integer counts, **not** TPM, to DESeq2

[[0,1]]{ .green }
### Count Model & EB
#### `PyDESeq2` (`dds`)
- Gamma–Poisson $\implies \text{Var} = \mu + \alpha\mu^2$
- `scipy`: $n = 1/\alpha$, $p = n/(n+\mu)$
- Median-of-ratios $s_j$ as GLM offset
- Cox–Reid + log-normal trend $\to \hat{\alpha}_i^{\text{MAP}}$
- Outlier protection ($+2\sigma_d$ & Cook's $D_{ij}$)

[[0,2]]{ .blue }
### Inference & QC
#### `DeseqStats` (`res_df`)
- IRLS log-link GLM ($\log_2\text{FC} = \beta_1/\ln 2$)
- Check reference level & positive controls
- Wald $z = \hat{\beta}_1/\text{SE}$ (or LRT for multi-level)
- Independent filtering + BH FDR (`padj`)
- apeGLM shrinkage for ranking & MA plots
:::

> **Core Equations at a Glance:** $E(r) = \bigcap_m C(k_m)$ · $\operatorname{Var}(K) = \mu + \alpha\mu^2$ · $n = 1/\alpha, \ p = n/(n+\mu)$ · $s_j = \operatorname{median}\_i(K_{ij}/g_i)$ · $\log \mu\_{ij} = \log s_j + \mathbf{x}\_j^\top \boldsymbol{\beta}\_i$ · $\hat{\alpha}\_i^{\text{MAP}} = \arg\max[\ell_{\text{CR}}(\alpha) - (\log\alpha - \log\alpha_{\text{tr}})^2/(2\sigma_d^2)]$ · $p\_{(k)}^{\text{adj}} = \min_{l \ge k} (m p_{(l)} / l)$
