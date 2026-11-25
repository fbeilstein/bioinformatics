:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Quality Filtering · $k$-mers · MinHash · Mash Distance
:::

---

# Curriculum

- **Quality Filtering & Adapter Detection (`fastp`)** — ASCII Phred decoding, per-read & sliding-window quality trimming, paired-end overlap analysis for de novo adapter detection, low-complexity & duplication estimation.
- **$k$-mer Decomposition & Canonical $k$-mers** — extracting all length-$k$ substrings ($L - k + 1$), strand-invariant canonical $k$-mers, choosing $k$ (random collision probability vs. mutation sensitivity).
- **Streaming Error Filtering & Cardinality Estimation** — Bloom filters / multiplicity thresholding (`-m 2`) to suppress singleton sequencing errors; estimating genome size $\hat{G}$ and sequencing coverage from hash density.
- **Hash Functions & Bottom-$s$ MinHash Sketching** — uniform 64-bit hashing (MurmurHash3), maintaining the $s$ smallest hashes via a max-heap in $O(N \log s)$ time and $O(s)$ memory.
- **Jaccard Similarity, ANI & Mash Distance** — unbiased merged-sketch Jaccard estimator, standard error $O(1/\sqrt{s})$, Average Nucleotide Identity (ANI) and the $\ge 95\%$ species boundary, Poisson point-mutation model converting $J(A, B)$ to mutation rate $D \approx 1 - \text{ANI}$, and breakdown regimes (high divergence $D > 0.2$, unequal genome sizes).
- **Statistical Significance ($p$-value)** — null random-match probability in $4^k$ space and binomial/hypergeometric tail test.
- **Complexity & Streaming I/O** — $O(G^2)$ exact alignment vs. $O(N \log s)$ sketching & $O(s)$ comparison; constant-memory $O(1)$ generator streaming (`Bio.SeqIO`).

---

# Key Math Notes

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Quality Scores & Trimming
$$Q = \operatorname{ord}(c) - 33 = -10\log_{10} P_{\text{error}}$$
Sliding window mean over window size $w$:
$$\bar{Q}_i = \frac{1}{w}\sum_{j=i}^{i+w-1} Q_j < Q_{\text{threshold}}$$

### Canonical $k$-mers & Choosing $k$
$$\operatorname{canon}(x) = \min\!\bigl(x,\; \overline{x}^{\,\text{rev}}\bigr)$$
Random collision probability in genome of size $G$:
$$q = 1 - \left(1 - \frac{1}{4^k}\right)^G \approx \frac{G}{4^k} \implies k \gtrsim \log_4\!\left(\frac{G}{q}\right)$$

### Genome Size & Coverage (Bottom-$s$)
Given the $s$-th smallest hash $h_{(s)}$ in $[0, H_{\max}]$:
$$\hat{G} \approx \frac{s \cdot H_{\max}}{h_{(s)}}, \qquad \widehat{\text{Coverage}} = \frac{\text{Total bases}}{\hat{G}}$$

[[0,1]]
### Bottom-$s$ MinHash & Jaccard Index
Let $S(A) = \min_s(h(A))$ and $W = \min_s\!\bigl(S(A) \cup S(B)\bigr)$:
$$J(A, B) = \frac{|A \cap B|}{|A \cup B|} \approx \frac{|W \cap S(A) \cap S(B)|}{|W|}$$
Standard error (independent of genome length $G$):
$$\operatorname{SE}(\hat{J}) = \sqrt{\frac{J(1 - J)}{s}} = O\!\left(\frac{1}{\sqrt{s}}\right)$$

### Average Nucleotide Identity (ANI) & Mash Distance
Mean base identity across aligned homologous blocks (species boundary: $\text{ANI} \ge 0.95 \iff D \lesssim 0.05$). Under Poisson mutation rate $d \approx 1 - \text{ANI}$, $k$-mer survival is $p = (1 - d)^k \approx e^{-kd}$:
$$J \approx \frac{e^{-kd}}{2 - e^{-kd}} \implies D = -\frac{1}{k}\ln\!\left(\frac{2J}{1+J}\right) \approx 1 - \text{ANI}$$
*Breaks down when $D > 0.2$ ($\text{ANI} < 80\%$, too few surviving $k$-mers) or when $|A| \ne |B|$ (unequal genome sizes inflate $D$ vs. $1 - \text{ANI}$, since ANI scores only shared aligned regions).*

### Statistical Significance ($p$-value)
With $r_1 = \frac{n_1}{4^k}$, $r_2 = \frac{n_2}{4^k}$, and $r = \frac{r_1 r_2}{r_1 + r_2 - r_1 r_2}$:
$$p\text{-value} = P(X \ge x) = \sum_{i=x}^{s} \binom{s}{i} r^i (1 - r)^{s - i}$$
:::
