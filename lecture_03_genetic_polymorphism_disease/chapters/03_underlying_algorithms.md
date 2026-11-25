:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Bayesian Genotyping · Hardy–Weinberg Priors · Haplotype Calling · Indexing
:::

---

# Curriculum

- **Read Alignment, Mapping Quality & Random-Access Indexing** — BWT / FM-index seeding (BWA-MEM overview), Phred-scaled mapping quality ($\text{MAPQ}$), BGZF block compression with hierarchical binning indexes (`.bai` / `.tbi` / `.crai`), reference-based CRAM delta encoding, and PCR duplicate marking (FLAG `0x400`) to preserve read independence.
- **Bayesian Diploid Read Classification** — Bayes' theorem for genotypes $G \in \{0/0, 0/1, 1/1\}$, independent per-read Phred error likelihoods ($\varepsilon_i = 10^{-Q_i/10}$), simplified biallelic vs. 4-base error models, and Maximum A Posteriori (MAP) calling.
- **Population Priors (Hardy–Weinberg Equilibrium)** — deriving genotype priors from population allele frequency $f$, and analyzing how strong read depth overwhelms the prior while shallow depth lets the prior decide the call.
- **Numerical Stability in Log Space** — preventing floating-point underflow via log-likelihood sums, the log-sum-exp normalization trick, and handling $Q = 0$ ($\varepsilon = 1$) boundary clamping.
- **Confidence Metrics & Haplotype Variant Calling** — site-level `QUAL` vs. sample-level `GL` and `GQ`; haplotype-based evaluation (resolving clustered SNPs/indels) and parsimonious left-alignment for canonical indel normalization.
- **Coverage & Allele-Fraction Sampling Statistics** — binomial sampling at heterozygotes ($k \sim \operatorname{Binomial}(n, \tfrac{1}{2})$), allelic dropout probability $2 \cdot (1/2)^n$ justifying depth cutoffs (`DP ≥ 15`), AAF variance, and reference mapping bias.

---

# Key Math Notes — Bayesian Genotyping & Confidence

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Bayes' Theorem & MAP Genotype
$$P(G \mid D) = \frac{P(D \mid G)\,P(G)}{\sum_{G' \in \{0/0, 0/1, 1/1\}} P(D \mid G')\,P(G')}$$
$$\hat{G}_{\text{MAP}} = \arg\max_G \bigl[\log P(D \mid G) + \log P(G)\bigr]$$

### Per-Read Likelihood ($\varepsilon_i = 10^{-Q_i/10}$)
For diploid genotype $G = a_1 / a_2$ and independent reads $b_1, \dots, b_n$:
$$P(D \mid G) = \prod_{i=1}^{n} \frac{P(b_i \mid a_1) + P(b_i \mid a_2)}{2}$$
In the biallelic model ($P(\text{match}) = 1 - \varepsilon_i$, $P(\text{mismatch}) = \varepsilon_i$):
- $P(D \mid 0/0) = \prod_{\text{ref}}(1 - \varepsilon_i)\prod_{\text{alt}}\varepsilon_i$
- $P(D \mid 0/1) = \prod_{\text{all}} \tfrac{1}{2} = \left(\tfrac{1}{2}\right)^n$
- $P(D \mid 1/1) = \prod_{\text{ref}}\varepsilon_i\prod_{\text{alt}}(1 - \varepsilon_i)$
*(Exact 4-base model uses $\varepsilon_i / 3$ for each specific wrong base.)*

[[0,1]]
### Hardy–Weinberg Priors (Allele Frequency $f$)
$$P(0/0) = (1 - f)^2, \quad P(0/1) = 2f(1 - f), \quad P(1/1) = f^2$$

### Log-Sum-Exp Normalization
With unnormalized log-posteriors $\ell_g = \log P(D \mid g) + \log P(g)$ and $m = \max_g \ell_g$:
$$P(G = g \mid D) = \frac{e^{\ell_g - m}}{\sum_{g'} e^{\ell_{g'} - m}}$$
*(Clamp $\varepsilon \in (0, 1)$ so $Q = 0$ never yields $\log 0 = -\infty$.)*

### VCF Confidence Scores & Sampling Error
- **Site `QUAL`**: $-10\log_{10} P(G = 0/0 \mid D)$
- **Genotype `GL`**: $\log_{10} P(D \mid G) - \max_{G'}\log_{10} P(D \mid G')$
- **Genotype `GQ`**: $-10\log_{10}\bigl(1 - P(\hat{G} \mid D)\bigr)$
- **Heterozygote allelic dropout** at depth $n$:
  $$P(\text{all } n \text{ reads from one copy}) = 2\left(\tfrac{1}{2}\right)^n$$
- **Alternate Allele Fraction**:
  $$\text{AAF} = \frac{AO}{RO + AO}, \qquad \operatorname{SD}(\text{AAF} \mid 0/1) = \frac{1}{2\sqrt{n}}$$
:::
