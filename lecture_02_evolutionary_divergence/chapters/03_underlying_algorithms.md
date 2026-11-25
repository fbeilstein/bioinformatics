:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Dynamic Programming · Gotoh · BLOSUM · Poisson Correction · Neighbor-Joining
:::

---

# Curriculum

- **Pairwise Alignment by Dynamic Programming** — Needleman–Wunsch (global), Smith–Waterman (local, zero-clamped recurrence), and semi-global (free end gaps); $O(NM)$ time and space, plus traceback & exponential co-optimal path multiplicity.
- **Affine Gap Penalties (Gotoh's Algorithm)** — gap cost $\text{gap}(L) = o + (L - 1)e$ (Biopython convention vs. $o + Le$); 3-state DP ($M, I_x, I_y$) reducing naive $O(N^3)$ affine DP to $O(NM)$.
- **Substitution Matrices & Log-Odds Scoring** — statistical derivation of BLOSUM ($S(a, b) = 2\log_2 \frac{p_{ab}}{q_a q_b}$ from clustered conserved blocks) vs. PAM Markov-chain extrapolation ($M^n$); DNA match/mismatch & transition/transversion scoring.
- **Evolutionary Distance & Multiple-Hits Correction** — percent identity definitions (gap-inclusive vs. gap-excluded), raw $p$-distance saturation, Poisson correction for proteins ($d = -\ln(1 - p)$), and Jukes–Cantor (JC69) / Kimura (K2P) models for DNA.
- **Distance-Based Phylogenetics** — additive vs. ultrametric tree metrics; UPGMA (assumes strict clock) vs. Neighbor-Joining ($Q$-matrix criterion, branch-length estimation, $O(n^3)$ complexity, unequal rates), and outgroup rooting.
- **Testing the Molecular Clock** — constant-rate expectation $d(t) = 2rt$, Ordinary Least Squares (OLS) regression ($\hat{\beta}_1, \hat{\beta}_0, R^2$), and statistical/biological caveats (site-rate heterogeneity, lineage rate shifts, fossil uncertainty, phylogenetic non-independence).

---

# Key Math Notes — Alignment & Scoring

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Smith–Waterman (Linear Gap $g$)
Zero-clamped local recurrence ($H_{i,0} = H_{0,j} = 0$):
$$H_{i,j} = \max\!\begin{cases} 0 \\ H_{i-1,j-1} + s(a_i, b_j) \\ H_{i-1,j} + g \\ H_{i,j-1} + g \end{cases}$$
*Global (Needleman–Wunsch) drops the $0$ branch and initializes borders with $i \cdot g$, $j \cdot g$.*

### BLOSUM Log-Odds Scoring
$$S(a, b) = 2 \log_2 \frac{p_{ab}}{q_a\,q_b}$$
$p_{ab}$: observed pair frequency in conserved blocks; $q_a q_b$: chance expectation (half-bit units).

[[0,1]]
### Gotoh's 3-Matrix Affine DP ($O(NM)$)
With gap open $o$ and extend $e$ ($\text{gap}(L) = o + (L - 1)e$):
$$I_{i,j} = \max\!\bigl(M_{i-1,j} + o,\; I_{i-1,j} + e\bigr)$$
$$D_{i,j} = \max\!\bigl(M_{i,j-1} + o,\; D_{i,j-1} + e\bigr)$$
$$M_{i,j} = \max\!\bigl([0,]\; M_{i-1,j-1},\; I_{i-1,j-1},\; D_{i-1,j-1}\bigr) + s(a_i, b_j)$$
*(Include $0$ clamp for local Gotoh; omit for global Gotoh. Both run in $O(NM)$ time and space.)*
:::

---

# Key Math Notes — Distances, Trees & Molecular Clock

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Multiple-Hits Distance Corrections
From observed fraction of differing sites $p = 1 - \text{identity}$:
- **Poisson (proteins)**: $1 - p = e^{-d}$
  $$d = -\ln(1 - p), \qquad \operatorname{Var}(\hat{d}) \approx \frac{p}{(1 - p)L}$$
- **Jukes–Cantor (JC69, DNA)** (4 equiprobable bases, saturates at $p = \frac{3}{4}$):
  $$d_{\text{JC}} = -\frac{3}{4}\ln\!\left(1 - \frac{4}{3}p\right)$$

### Molecular Clock & OLS Regression
Two lineages diverging $t$ MYA span total path $2t$:
$$d(t) = 2r \cdot t, \qquad \hat{\beta}_1 = \frac{\sum (t_i - \bar{t})(d_i - \bar{d})}{\sum (t_i - \bar{t})^2}, \qquad R^2 = r_{td}^2$$

[[0,1]]
### Neighbor-Joining (Saitou & Nei, $O(n^3)$)
1. **Select closest pair** minimizing the rate-corrected $Q$-matrix:
   $$Q(i, j) = (n - 2)\,d(i, j) - \sum_{k=1}^{n} d(i, k) - \sum_{k=1}^{n} d(j, k)$$
2. **Branch lengths** from pair $(i, j)$ to new internal node $u$:
   $$\delta(i, u) = \frac{1}{2}d(i, j) + \frac{1}{2(n - 2)}\!\left(\sum_{k} d(i, k) - \sum_{k} d(j, k)\right)$$
   $$\delta(j, u) = d(i, j) - \delta(i, u)$$
3. **Update distances** to every remaining taxon $k$:
   $$d(u, k) = \frac{1}{2}\bigl(d(i, k) + d(j, k) - d(i, j)\bigr)$$
:::
