:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Dynamic Programming · Gotoh · BLOSUM · Poisson & JC69 · Neighbor-Joining · Molecular Clock
:::

---

# Curriculum & Algorithmic Roadmap

How raw orthologous sequences become **phylogenetic trees** and **molecular clock rates** in six algorithmic stages:

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "seq",   "label": "Orthologous\nSequences",          "x": 30,  "y": 90, "width": 145, "height": 64 },
    { "id": "sub",   "label": "Log-Odds Scoring\nBLOSUM / DNA",  "x": 220, "y": 90, "width": 165, "height": 64, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "dp",    "label": "Gotoh 3-State DP\nO(NM) Affine",  "x": 430, "y": 90, "width": 170, "height": 64, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "dist",  "label": "Multiple Hits\nPoisson / JC69",   "x": 645, "y": 90, "width": 165, "height": 64, "style": "fill: #FF9800; color: #fff;" },
    { "id": "nj",    "label": "Neighbor-Joining\nQ-Matrix O(n³)","x": 855, "y": 40, "width": 175, "height": 62, "style": "fill: #2196F3; color: #fff;" },
    { "id": "clk",   "label": "Molecular Clock\nOLS d = 2rt",    "x": 855, "y": 145,"width": 175, "height": 62, "style": "fill: #009688; color: #fff;" }
  ],
  "links": [
    { "from": "seq",  "to": "sub" },
    { "from": "sub",  "to": "dp" },
    { "from": "dp",   "to": "dist" },
    { "from": "dist", "to": "nj" },
    { "from": "dist", "to": "clk" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="28px" margin-top: 10px; }
[[0,0]]
- **Pairwise DP Alignment** — Needleman–Wunsch (global), Smith–Waterman (local), and semi-global; $O(NM)$ time/space & co-optimal path explosion.
- **Affine Gap Penalties (Gotoh)** — $\text{gap}(L) = o + (L - 1)e$; 3-state $(M, I, D)$ DP collapsing naive $O(N^3)$ to $O(NM)$.
- **Substitution Matrices** — log-odds derivation of BLOSUM vs. PAM Markov extrapolation & DNA scoring.

[[0,1]]
- **Multiple-Hits Distance Corrections** — raw $p$-distance saturation, Poisson ($d = -\ln(1 - p)$), JC69 ($p \to 0.75$ ceiling), and variance explosion.
- **Distance-Based Phylogenetics** — additive vs. ultrametric metrics, why UPGMA fails, and Neighbor-Joining ($Q$-matrix + outgroup rooting).
- **Molecular Clock Testing** — $d(t) = 2rt$, OLS regression ($\hat{\beta}\_1, \hat{\beta}\_0, R^2$), and statistical/biological caveats.
:::

---

# Why Dynamic Programming? Edit Graphs & Optimal Substructure

Aligning two sequences $A = a_1 \dots a_N$ and $B = b_1 \dots b_M$ means inserting gaps `-` to maximize total column score. Why can't we just try every alignment?

:::matrix { cols="1.05fr 0.95fr" gap="30px" }
[[0,0]]
### Combinatorial Explosion (Delannoy Numbers)
Counting all valid paths from $(0, 0)$ to $(N, M)$ using diagonal, vertical, and horizontal steps gives the **Delannoy number**:
$$D(N, M) = \sum_{k=0}^{\min(N,M)} \binom{N}{k}\binom{M}{k} 2^k \sim \frac{(1 + \sqrt{2})^{2N + 1}}{\sqrt{\pi N}}$$

| Sequence Length $N = M$ | Possible Alignments $D(N, N)$ |
| :---: | :---: |
| $10$ | $8{,}097{,}453$ |
| $30$ | $1.04 \times 10^{22}$ |
| **$100$** | **$1.7 \times 10^{75}$** *(exceeds atoms in Earth)* |
| **$142$ ($\alpha$-globin)** | **$> 10^{107}$** *(exceeds atoms in universe!)* |

[[0,1]]
### Optimal Substructure on an $(N{+}1) \times (M{+}1)$ DAG
Every alignment is a path from $(0, 0)$ to $(N, M)$ on a **Directed Acyclic Graph (DAG)** whose last column must be one of three moves:

1. **Diagonal $\searrow$** $(i{-}1, j{-}1) \to (i, j)$: align $a_i$ with $b_j$, score $+s(a_i, b_j)$.
2. **Vertical $\downarrow$** $(i{-}1, j) \to (i, j)$: align $a_i$ with `-` (deletion in $B$), score $+g$.
3. **Horizontal $\rightarrow$** $(i, j{-}1) \to (i, j)$: align `-` with $b_j$ (insertion in $B$), score $+g$.

Because the optimal prefix score $H_{i,j}$ depends **only** on its three already-solved neighbors, **Dynamic Programming** finds the exact global optimum in just **$O(NM)$** operations — ~20,000 steps instead of $10^{107}$!
:::

---
# Global, Local & Semi-Global Alignment

One single recurrence family generates all three biological alignment modes by changing only **border initialization**, the **zero clamp**, and **traceback endpoints**:

![Global, Local, and Semi-Global DP Edit Graphs](lecture_02_evolutionary_divergence/images/dp_edit_graph_modes.png){width=84 center}

:::matrix { .smaller-table .dense cols="1fr" margin-top: 6px; }
[[0,0]]
| Mode | Recurrence $H_{i,j}$ (Linear Gap $g$) | Borders ($H_{i,0}, H_{0,j}$) | Traceback Start $\to$ Stop | Biological Use Case |
| :--- | :--- | :---: | :---: | :--- |
| **Global (NW)** | $\max(H_{i-1,j-1}+s,\ H_{i-1,j}+g,\ H_{i,j-1}+g)$ | $i\cdot g,\ j\cdot g$ | $(N, M) \to (0, 0)$ | Full-length orthologous proteins (`HBA1`, `CYCS`) |
| **Local (SW)** | $\max(\mathbf{0},\ H_{i-1,j-1}+s,\ H_{i-1,j}+g,\ H_{i,j-1}+g)$ | $0,\ 0$ | $\arg\max_{i,j} H_{i,j} \to H = 0$ | mRNAs with variable UTRs, conserved domains |
| **Semi-Global** | Same as NW (no interior $0$ clamp) | $0,\ 0$ | $\max(\text{last row/col}) \to \text{row/col } 0$ | Reads vs. genome, suffix–prefix overlaps |
:::

---
# Traceback, Co-Optimal Multiplicity & Complexity

:::matrix { cols="48/52" gap="10px" }
[[0,0]]
Filling $H_{i,j}$ only gives the **optimal score** --- recovering the actual gapped alignment requires tracing winning backpointers backward.

### Why Co-Optimal Alignments Explode ($2^k$)
- Whenever two transitions tie at cell $(i, j)$ (e.g. `↖` and `↑` both give the maximum score), a **diamond fork** opens in the backpointer graph.
- Across a sequence with $k$ independent tied regions (common in repeats or low-complexity DNA), path counts multiply:
  $$\\#\text{co-optimal alignments} = \prod_{r=1}^{k} c_r \ge 2^k$$
- For two random $2\text{ kb}$ sequences, $k > 63 \implies 2^k > 2^{63} - 1$ — directly causing Part 4's **`OverflowError`** on `len(aligner.align(s, t))`!
- **Rule of thumb:** Always take `aligner.align(a, b)[0]` (lazy single path) or `aligner.score(a, b)` — never materialize `list(...)`.

[[0,1]]
### Time & Space Complexity Trade-Offs

| Operation | Time | Memory | Mechanism |
| :--- | :---: | :---: | :--- |
| Full DP + traceback | $O(NM)$ | $O(NM)$ | Store full score / backpointer grid |
| Score only (`aligner.score`) | $O(NM)$ | **$O(\min(N,M))$** | Rolling 2-row buffer (previous + current) |
| **Hirschberg** (linear space) | $O(NM)$ | **$O(\min(N,M))$** | Divide-and-conquer midpoint recursion (2× time) |

### Karlin–Altschul Condition for Local Alignment
Smith–Waterman only stays **local** when random alignments have **negative expected score**:
$$\mathbb{E}[S] = \sum_{a,b} q_a q_b s(a, b) < 0$$
If $\mathbb{E}[S] > 0$, random flanking noise adds positive drift and local alignment degrades into global alignment!
:::

---
# Affine Gap Penalties & Gotoh's 3-State Automaton

Linear gaps charge $-g$ per deleted base — penalizing one 3-bp codon deletion (`ΔF508`) **3× as harshly** as a 1-bp deletion, even though replication slippage creates multi-base indels in a **single event**.

:::static-diagram { width=90 center }
{
  "nodes": [
    { "id": "I", "label": "State I (Gap in B)\nExtend: +e  |  Vertical ↑",     "x": 60,  "y": 140, "width": 230, "height": 68, "style": "fill: #FF9800; color: #fff;" },
    { "id": "M", "label": "State M (Match / Mismatch)\nScore: +s(aᵢ, bⱼ)  |  Diag ↖", "x": 405, "y": 30,  "width": 250, "height": 72, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "D", "label": "State D (Gap in A)\nExtend: +e  |  Horizontal ←",   "x": 770, "y": 140, "width": 230, "height": 68, "style": "fill: #2196F3; color: #fff;" }
  ],
  "links": [
    { "from": "M", "to": "I", "label": "open o" },
    { "from": "I", "to": "M", "label": "+s(aᵢ,bⱼ)" },
    { "from": "M", "to": "D", "label": "open o" },
    { "from": "D", "to": "M", "label": "+s(aᵢ,bⱼ)" },
    { "from": "I", "to": "D", "label": "open o" },
    { "from": "D", "to": "I", "label": "open o" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="30px" margin-top: 8px; }
[[0,0]]
### Affine Penalty Conventions ($o < e < 0$)
- **Biopython / Gotoh convention** (used in Part 4):
  $$\text{gap}(L) = o + (L - 1)e \qquad (o = -10,\ e = -1 \implies \text{gap}(3) = -12)$$
- **Alternative textbook convention:** $\text{gap}(L) = o' + L\cdot e$ (shifted by one extend: $o = o' + e$).

[[0,1]]
### Why Naive Affine DP Is $O(N^3)$ — and Gotoh Is $O(NM)$
- With one matrix $H_{i,j}$, knowing whether the next gap costs $o$ or $e$ forces scanning **all prior gap lengths** $k \in [1, i]$ and $l \in [1, j]$ $\implies O(NM(N + M)) = O(N^3)$!
- **Gotoh's insight (1982):** Maintain **3 coupled matrices** $(M, I, D)$ encoding the ending state of the prefix $\implies$ every transition is $O(1)$, restoring **$O(NM)$** time and space!
:::

---

# Gotoh's 3-Matrix Recurrence & Worked Verification

Using Biopython's convention $\text{gap}(L) = o + (L - 1)e$ (gap of length 1 costs $o$, each extra residue costs $e$):

:::matrix { cols="1.12fr 0.88fr" gap="28px" }
[[0,0]]
### Coupled Recurrence ($1 \le i \le N,\ 1 \le j \le M$)
- **Vertical gap** ($a_i$ aligned to `-`):
  $$I_{i,j} = \max\bigl(M_{i-1,j} + o,\ I_{i-1,j} + e,\ D_{i-1,j} + o\bigr)$$
- **Horizontal gap** (`-` aligned to $b_j$):
  $$D_{i,j} = \max\bigl(M_{i,j-1} + o,\ D_{i,j-1} + e,\ I_{i,j-1} + o\bigr)$$
- **Match / mismatch** ($a_i$ aligned to $b_j$):
  $$M_{i,j} = \max\bigl([0,]\ M_{i-1,j-1} + s,\ I_{i-1,j-1} + s,\ D_{i-1,j-1} + s\bigr)$$
  *(Include $0$ clamp for **local** Gotoh; omit for **global** / **semi-global**.)*

### Boundary Initialization
- **Global:** $M_{0,0} = 0$, $M_{i,0} = M_{0,j} = -\infty$; $I_{i,0} = o + (i{-}1)e$, $D_{0,j} = o + (j{-}1)e$; opposite borders $= -\infty$.
- **Local:** $M_{i,0} = M_{0,j} = 0$; $I_{i,0} = I_{0,j} = D_{i,0} = D_{0,j} = -\infty$ (local alignments never open on a gap).

[[0,1]]
### Exact Check: Part 4's Worked Example
Global Gotoh on `HEAGAWGHEE` vs. `PAWHEAE` (`BLOSUM62`, $o = -10$, $e = -1$):

```
target   H E A G A W G H E E
         - - - . | | . . . |
query    - - - P A W H E A E
```

- Leading 3-aa gap stays in state $I$:
  $$I_{1,0} = -10 \to I_{2,0} = -11 \to I_{3,0} = \mathbf{-12}$$
- Jumps $I_{3,0} \to M_{4,1}$ on `(G, P)` ($-2$), then rides diagonal $M$:
  $$-12 + (-2 + 4 + 11 - 2 + 0 - 1 + 5) = \mathbf{+3}$$
- Switch to **local Gotoh** on the exact same pair: ignores flanks and isolates `AWGHE / AW-HE` (1 internal gap $-10$, score **+18**)!
:::

---

# Interactive Lab — Pairwise DP & Gotoh Explorer

Explore Global, Local, Semi-Global, and Gotoh's 3 matrices ($M, I, D$). Try these **4 checks**: (1) verify `HEAGAWGHEE` vs. `PAWHEAE` gives **score = 3** in global and **18** in local; (2) switch to **mRNA + UTRs** and compare Global vs. Local vs. Semi-Global; (3) switch to **3-bp Slippage** and flip `Affine` $\leftrightarrow$ `Linear (-3)` to see linear gaps split the indel; (4) load **Co-Optimal Ties** to watch tied branches multiply!

<iframe src="lecture_02_evolutionary_divergence/demo/dp_alignment_explorer.html" title="Interactive Pairwise DP and Gotoh Alignment Explorer" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---
# Substitution Matrices: Statistical Log-Odds Derivation

For proteins, a single mismatch score ignores biochemistry (`Leu ↔ Ile` preserves hydrophobic packing; `Trp ↔ Gly` collapses the core). Where do substitution matrix integers actually come from?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Likelihood-Ratio Framework
For an aligned residue pair $(a, b)$, compare two statistical models:
- **$H_1$ (Homologous):** pair $(a, b)$ evolved from a shared ancestor with empirical joint probability $p_{ab}$.
- **$H_0$ (Chance):** independent draws from background frequencies $q_a q_b$.

Assuming column independence, the log-likelihood ratio of a full alignment **splits additively** into column scores (in **half-bits**):
$$S(a, b) = \operatorname{round}\left(2 \log_2 \frac{p_{ab}}{q_a q_b}\right)$$

- $S(a, b) > 0$: observed **more often** in homologs than by chance.
- $S(a, b) < 0$: observed **less often** than by chance (selected against).
- **Relative entropy** (average information per aligned homologous site):
<div style="margin-top: -44px;"></div>
$$
  H = \sum_{a, b} p_{ab} \log_2 \frac{p_{ab}}{q_a q_b} > 0\text{ bits}
$$

[[0,1]]
### Why `W–W = +11` While `A–A = +4`
- Tryptophan is rare ($q_W \approx 0.013$, so chance $q_W^2 \approx 0.00017$); Alanine is common ($q_A \approx 0.074$, $q_A^2 \approx 0.0055$).
- Matching a conserved rare Trp yields $p_{WW} / q_W^2 \approx 45\times$ ($\implies 2\log_2 45 \approx \mathbf{+11}$ half-bits!), far stronger homology evidence than matching common Ala!

### DNA Scoring & Transition/Transversion Bias
- Practice 2 uses `match = +5`, `mismatch = -4`: under uniform DNA ($q = 0.25$),
  $$\mathbb{E}[S] = \tfrac{1}{4}(+5) + \tfrac{3}{4}(-4) = \mathbf{-1.75} < 0 \quad \text{(safe for SW!)}$$
- Biologically refined DNA matrices penalize **transitions** ($\text{A}\leftrightarrow\text{G}$, $\text{C}\leftrightarrow\text{T}$) less harshly than **transversions** (purine $\leftrightarrow$ pyrimidine), reflecting the ~2:1 biochemical bias.
:::

---
# BLOSUM62: Biochemical Structure in Log-Odds Scores

Ordering the $20 \times 20$ **BLOSUM62** matrix by amino-acid chemistry reveals how evolution preserves biophysical properties:

:::matrix { cols="1.12fr 0.88fr" gap="28px" }
[[0,0]]
![BLOSUM62 Chemically Grouped Heatmap](lecture_02_evolutionary_divergence/images/blosum62_heatmap.png){width=94% center}

[[0,1]]
### Reading the Matrix
- **Positive diagonal blocks:** Substitutions *within* the same chemical family stay favorable (`I ↔ V = +3`, `F ↔ Y = +3`, `D ↔ E = +2`, `K ↔ R = +2`).
- **Rare structural anchors:** Bulky aromatic **`W–W (+11)`** and disulfide-forming **`C–C (+9)`** score highest because both are rare and structurally irreplaceable.
- **Cross-class clashes:** Swapping hydrophobic core residues for charged ones (`L ↔ D = -4`) or bulky Trp for tiny Gly (`W ↔ G = -2`) is strongly penalized.
- **Expected random score:** $\mathbb{E}[S] = \sum q_a q_b S(a,b) \approx \mathbf{-1.00}$ half-bits ($-0.50$ bits $< 0$) — valid for local Smith–Waterman alignment.
:::

---
# BLOSUM vs. PAM: Empirical Clustering vs. Markov Extrapolation


:::matrix { cols="50/50" gap="28px" }
[[0,0]]
:::matrix { .cards cols="1fr" gap="16px" }
[[0,0]]{ .blue }
### PAM (Dayhoff et al., 1978) — Extrapolation
- Estimate 1-step transition matrix $M_1$ (1 accepted mutation / 100 sites) from very close homologs ($\ge 85\\%$ identity).
- Extrapolate to distant divergence by **Markov matrix powers**: $M_n = M_1^n$ (`PAM250` $= M_1^{250}$).
- **Higher PAM number = MORE divergent** (`PAM30` close $\to$ `PAM250` distant).

[[1,0]]{ .orange }
### BLOSUM (Henikoff & Henikoff, 1992) — Direct Empirical
- Count residue pairs directly inside ungapped conserved **BLOCKS** across real divergent protein families.
- Cluster sequences sharing $\ge X\\%$ identity into one weighted profile to prevent redundant-family bias (`BLOSUM62` $= 62\\%$ cutoff).
- **Lower BLOSUM number = MORE divergent** (`BLOSUM80` close $\to$ `BLOSUM62` default $\to$ `BLOSUM45` distant).
:::

[[0,1]]
Both matrix families build half-bit log-odds tables — but from **opposite statistical philosophies** and **opposite numbering directions**!

![PAM vs BLOSUM Opposing Numbering and Relative Entropy](lecture_02_evolutionary_divergence/images/pam_vs_blosum_scale.png){width=98% center}
:::

---

# Evolutionary Distance & Multiple-Hits Correction

After aligning two sequences, we convert percent identity into an **evolutionary distance** $d$ (expected substitutions per site). Why isn't raw difference $p = 1 - \text{identity}$ enough?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### 1. Identity Conventions & Hidden Hits
- **Gap-inclusive** (`matches / aln_len`, Practice 2 local DNA) treats indels as differences; **gap-excluded** (`matches / comparable`, Homework Task 1) isolates pure residue substitution.
- **Hidden multiple hits:** As time grows, the same site mutates repeatedly:
  - **Double hit** ($\text{A} \to \text{G} \to \text{T}$): 2 events, observed as **1** difference.
  - **Back-mutation** ($\text{A} \to \text{G} \to \text{A}$): 2 events, observed as **0** differences!

### 2. Poisson Correction (Proteins — Homework)
Assume hits per site follow $K \sim \operatorname{Poisson}(d)$ and alphabet is large enough that back-mutations are rare. A site stays identical iff $K = 0$:
$$1 - p = P(K = 0) = e^{-d} \implies d = -\ln(1 - p)$$
*(Breaks down as $p \to 1$; ignores unequal amino-acid frequencies and site-rate variation.)*

[[0,1]]
### 3. Jukes–Cantor Model (JC69 — DNA 4-State CTMC)
For 4 equiprobable bases mutating symmetrically at rate $\alpha$ to each of the 3 other bases (total rate $3\alpha$ per branch), the probability $P_{\text{same}}(t)$ obeys:
$$\frac{dP_{\text{same}}}{dt} = -3\alpha P_{\text{same}} + \alpha(1 - P_{\text{same}}) = \alpha - 4\alpha P_{\text{same}}$$
Solving from $P_{\text{same}}(0) = 1$ with total two-branch distance $d = 6\alpha t$:
$$p(d) = 1 - P_{\text{same}} = \frac{3}{4}\left(1 - e^{-\frac{4}{3}d}\right) \implies d_{\text{JC}} = -\frac{3}{4}\ln\left(1 - \frac{4}{3}p\right)$$

### Richer Evolutionary Extensions
- **Kimura 2-Parameter (K2P, DNA):** Separates transition fraction $P$ and transversion fraction $Q$:
  $$d_{\text{K2P}} = -\tfrac{1}{2}\ln(1 - 2P - Q) - \tfrac{1}{4}\ln(1 - 2Q)$$
- **Gamma-Poisson (proteins, rate shape $\alpha$):** $d_\Gamma = \alpha\bigl[(1 - p)^{-1/\alpha} - 1\bigr]$.
:::

---
# Saturation Ceilings & Variance Explosion

![Multiple-Hits Saturation Curves and Delta-Method Variance Explosion](lecture_02_evolutionary_divergence/images/saturation_multiple_hits.png){width=73 center}

:::matrix { cols="1fr 1fr" gap="30px" margin-top: 8px; }
[[0,0]]
Multiple-hits correction recovers **unbiased** distances on average — but cannot defeat the **information loss** as sequences approach random saturation!

### Why Proteins Outlast DNA on Deep Clocks
- Random 4-letter DNA still matches at **$25\\%$** of sites $\implies$ **JC69 hits an impenetrable wall at $p = 0.75$** ($d_{\text{JC}} \to \infty$!).
- Random 20-letter proteins match at only **$\sim 5\text{–}6\\%$** and evolve slower under purifying selection, keeping $p$ far below saturation over $400{+}\text{ MYA}$.

[[0,1]]
### Delta-Method Variance Proof
Observed $\hat{p} = X/L$ is Binomial with $\operatorname{Var}(\hat{p}) = \frac{p(1 - p)}{L}$. Propagating through $d(p) = -\ln(1 - p)$ with derivative $\frac{dd}{dp} = \frac{1}{1 - p}$:
$$\operatorname{Var}(\hat{d}) \approx \left(\frac{dd}{dp}\right)^2 \operatorname{Var}(\hat{p}) = \frac{p}{(1 - p)L}$$
As $p \to 1$ (or $p \to 0.75$ for JC69), $(1 - p)$ vanishes and **sampling noise explodes**!
:::

---

# Interactive Lab — BLOSUM62 Log-Odds & Saturation Simulator

Two hands-on tabs: **(1)** click `W–W (+11)` vs. `A–A (+4)` in BLOSUM62 and drag the score offset $\Delta$ past `+1.0` to trigger the Karlin–Altschul $\mathbb{E}[S] \ge 0$ breakdown; **(2)** switch to **Multiple-Hits Simulator**, push DNA divergence $d \to 2.0$, and watch purple **back-mutations (`↺`)** and the $p = 0.75$ JC69 pole blow up the 95% error bar!

<iframe src="lecture_02_evolutionary_divergence/demo/substitution_distance_lab.html" title="Interactive BLOSUM62 and Multiple-Hits Saturation Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Distance-Based Phylogenetics: Ultrametric vs. Additive Metrics

Given an $n \times n$ pairwise distance matrix $D$, when does $D$ correspond **exactly** to the path lengths along an evolutionary tree?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Ultrametric Metric (Strict Molecular Clock)
A distance matrix is **ultrametric** iff for every triple of taxa $(i, j, k)$, the two largest distances are **equal**:
$$d(i, k) = d(j, k) \ge d(i, j)$$
(for some labeling of the triple).

- **Geometric meaning:** There exists a rooted tree where **every leaf has the exact same root-to-tip distance** $h$.
- **Biological meaning:** Every lineage evolved at the **exact same constant rate** since the root!
- **Algorithm:** **UPGMA** (Unweighted Pair Group Method with Arithmetic Mean) greedily merges the closest pair $\min d(i, j)$, setting node height $h(u) = \frac{1}{2} d(i, j)$ and averaging cluster distances weighted by cluster size.

[[0,1]]
### Additive Metric (Variable Branch Rates)
A distance matrix is **additive** (a general tree metric with arbitrary positive branch lengths) iff it satisfies **Buneman's 4-Point Condition**:

For every quartet $\{i, j, k, l\}$, compute the three pairings:
$$S_1 = d(i,j) + d(k,l), \quad S_2 = d(i,k) + d(j,l), \quad S_3 = d(i,l) + d(j,k)$$
Then the **two largest sums must be equal**!

- **Why?** If the true unrooted quartet splits $(i, j) \mid (k, l)$ across internal branch $w > 0$, then $S_2 = S_3 = S_1 + 2w \ge S_1$!
- Every ultrametric matrix is additive, but real genomes with fast and slow lineages are **additive (up to noise), NOT ultrametric**.
:::

---
# Why UPGMA Fails Under Unequal Rates (Clock Trap)

If sister lineages evolve at different rates, **minimum raw distance $\min d(i, j)$ does NOT guarantee true neighbors**!

![UPGMA Failure vs Neighbor-Joining Success on 4-Taxon Additive Tree](lecture_02_evolutionary_divergence/images/upgma_vs_nj_counterexample.png){width=88 center}

:::matrix { cols="1fr 1fr" gap="30px" margin-top: 8px; }
[[0,0]]
### The Greedy Failure Mode
- On the exact additive tree `((A:1, B:9):3, (C:1, D:9))`, the two slow-evolving **non-sisters** $A$ and $C$ are only $1 + 3 + 1 = \mathbf{5}$ apart, while true sisters $(A, B)$ are $1 + 9 = \mathbf{10}$ apart!
- **UPGMA** blindly grabs $\min D = d(A, C) = 5$ and produces the **wrong topology** `(((A, C), D), B)`!

[[0,1]]
### Why NJ Is Essential Before Testing a Clock
- Because UPGMA **forces** a strict clock onto the tree, using UPGMA would circularly bake clock assumptions into our phylogeny.
- **Neighbor-Joining** recovers the exact topology and exact branch lengths (`1.0, 9.0, 3.0, 1.0, 9.0`) with **zero clock assumption**!
:::

---

# Neighbor-Joining (Saitou & Nei): The $Q$-Matrix Algorithm

How does Neighbor-Joining cancel out long fast-evolving terminal branches so it always selects a true cherry?

:::matrix { cols="1.15fr 0.85fr" gap="28px" }
[[0,0]]
### Three-Step Iterative Reduction ($n \to n - 1$)
Let $R_i = \sum_{k=1}^{n} d(i, k)$ be the total distance from taxon $i$ to all other active nodes.

1. **Minimize the rate-corrected $Q$-matrix** (minimum-evolution criterion):
   $$Q(i, j) = (n - 2)\ d(i, j) - R_i - R_j$$
   *Subtracting $R_i$ and $R_j$ cancels terminal branch lengths and isolates the internal branch separating $\{i, j\}$ from all other taxa!*
2. **Estimate asymmetric branch lengths** to new parent node $u$:
   $$\delta(i, u) = \frac{1}{2}d(i, j) + \frac{R_i - R_j}{2(n - 2)}, \qquad \delta(j, u) = d(i, j) - \delta(i, u)$$
   *(The faster-evolving taxon has larger $R_i$ and automatically gets the longer branch!)*
3. **Update distances** from $u$ to every remaining taxon $k$:
   $$d(u, k) = \frac{1}{2}\bigl(d(i, k) + d(j, k) - d(i, j)\bigr)$$

[[0,1]]
### Theoretical Guarantees & Complexity
- **Exactness theorem (Studier & Keppler, 1988):** For *any* additive matrix, $\arg\min_{i \neq j} Q(i, j)$ is **mathematically guaranteed** to be a true pair of neighbors!
- **Complexity:** $n - 3$ merge iterations $\times$ $O(n^2)$ scan of $Q(i, j)$ $\implies$ **$O(n^3)$ time**, **$O(n^2)$ space** (instant for $n \le 1000$).
- **Outgroup Rooting:** Because NJ depends only on pairwise path sums, its output is an **unrooted** tree ($2n - 3$ branches).
  - Adding a known early-diverging **outgroup** (sea lamprey *Petromyzon marinus*, ~600 MYA) locates the root without assuming equal rates.
:::

---
# Interactive Lab — UPGMA vs. Neighbor-Joining Builder

Step through NJ's $Q$-matrix merges side-by-side against UPGMA: **(1)** on the default **Clock Trap**, watch orange $\min D = d(A,C) = 5$ mislead UPGMA while blue $\min Q = Q(A,B) = -52$ nails the exact tree; **(2)** switch to **Vertebrate Globin-Like (fast Mouse)** and change the **NJ Outgroup**; **(3)** slide **Noise $\sigma$** up to test Buneman 4-point robustness!

<iframe src="lecture_02_evolutionary_divergence/demo/neighbor_joining_lab.html" title="Interactive UPGMA vs Neighbor-Joining Tree Builder" style="width: 100%; height: 600px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Testing the Molecular Clock: OLS Regression & Rate Estimation

Once we have Poisson-corrected distances $d_i$ between human and species $i$ with fossil divergence time $t_i$ (MYA), how do we quantify clock speed and linearity?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Constant-Rate Expectation: Why the Factor of $2$?
- After splitting $t$ MYA ago, **both** lineages accumulate substitutions independently at per-lineage rate $r$ (substitutions / site / MYA).
- Total evolutionary path separating the two modern proteins is $2t$:
  $$d(t) = 2r \cdot t$$
- Fitting a straight line $d_i = \beta_0 + \beta_1 t_i + \varepsilon_i$ by **Ordinary Least Squares (OLS)** minimizes $\sum_{i=1}^n (d_i - \beta_0 - \beta_1 t_i)^2$:
  $$\hat{\beta}\_1 = \frac{\sum (t\_i - \bar{t})(d\_i - \bar{d})}{\sum (t\_i - \bar{t})^2}, \qquad \hat{\beta}\_0 = \bar{d} - \hat{\beta}\_1 \bar{t}$$
  $$R^2 = r_{td}^2 = \frac{\bigl[\sum (t_i - \bar{t})(d_i - \bar{d})\bigr]^2}{\sum (t_i - \bar{t})^2 \sum (d_i - \bar{d})^2}$$

[[0,1]]
### Interpreting Diagnostics (`scipy.stats.linregress`)
- **Slope $\hat{\beta}\_1 = 2\hat{r}$:** Direct measure of evolutionary rate — purifying selection makes `CYCS` shallow, `HBA1` intermediate, and `FGA` steep.
- **Intercept $\hat{\beta}\_0$:** Under an ideal clock, $\beta\_0 = 0$.
  - A **positive $\hat{\beta}\_0$** in raw $p$-distance (or negative curvature in `% identity`) is the classic signature of **multiple-hits saturation** bending deep points downward!
  - Switching from raw $p$ to Poisson $d = -\ln(1-p)$ straightens the curve, raising $R^2$ most dramatically for the **fastest-evolving protein** (`FGA`).
- **Per-lineage rate:** Remember to divide the regression slope by **2**: $\hat{r} = \hat{\beta}\_1 / 2$!
:::

---

# When the Molecular Clock Breaks: Four Key Caveats

A high $R^2$ in a simple OLS scatter plot is an informative starting point — but **real molecular-clock dating** must confront four fundamental biological and statistical pitfalls:

:::matrix { .cards cols="1fr 1fr" gap="24px" margin-top: 12px; }
[[0,0]]{ .purple }
### 1. Phylogenetic Non-Independence
Comparing every species back to Human **reuses shared internal branches** — Human–Chicken and Human–Zebrafish share the entire early vertebrate stem! Errors along one deep branch shift multiple data points together, violating OLS independence and overstating $R^2$ / $p$-values.

[[0,1]]{ .orange }
### 2. Among-Site Rate Heterogeneity
Proteins do not mutate uniformly: catalytic cores are nearly frozen while surface loops evolve rapidly. A single Poisson rate underestimates deep divergence because fast sites saturate first — requiring **Gamma-distributed site rates** ($\Gamma$ model, shape $\alpha$).

[[1,0]]{ .blue }
### 3. Lineage-Specific Rate Shifts
The molecular clock is **overdispersed** (local / relaxed clocks): short-generation organisms (e.g., rodents) accumulate more replication errors per million years than long-lived primates, and episodes of positive selection accelerate individual branches.

[[1,1]]{ .green }
### 4. Fossil Calibration Uncertainty
Divergence times $t\_i$ are **not exact numbers** — fossils provide only minimum bounds with wide geological uncertainty. Measurement error in the regressor $t\_i$ violates OLS assumptions, causing **attenuation bias** (flattening $\hat{\beta}\_1$) unless modeled with Bayesian relaxed clocks (e.g., BEAST).
:::

---

# Summary & Bridge to Part 4 (Industrial Tools)

Every algorithm from Part 3 powers a concrete tool in **Part 4** and **Practice 2**:

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .green }
### Pairwise Alignment
#### `Bio.Align.PairwiseAligner`
- Global (NW), Local (SW), Semi-Global (`end_gap_score=0`)
- Gotoh 3-state $(M, I, D)$ affine DP in $O(NM)$
- Biopython convention: $o + (L - 1)e$
- `BLOSUM62` half-bit log-odds ($\mathbb{E}[S] < 0$)
- Avoid `list(aligner.align(...))` ($2^k$ ties!)

[[0,1]]{ .orange }
### Distance & Saturation
#### `calculate_identity` + `poisson_distance`
- Gap-inclusive vs. gap-excluded identity
- Poisson correction: $d = -\ln(1 - p)$
- JC69 DNA ceiling at $p = 0.75$
- Variance $\operatorname{Var}(\hat{d}) \approx \frac{p}{(1 - p)L}$ explodes near saturation

[[0,2]]{ .blue }
### Trees & Clock
#### `Bio.Phylo` + `scipy.stats.linregress`
- Additive (4-point) vs. ultrametric (strict clock)
- NJ $Q$-matrix ($O(n^3)$) beats greedy UPGMA under unequal rates
- Outgroup rooting (*Petromyzon marinus*)
- OLS clock fit: $d(t) = 2rt$, $\hat{r} = \hat{\beta}\_1 / 2$
:::

> **Quick-Reference Core Equations:**  
> $S(a,b) = 2\log_2 \frac{p_{ab}}{q_a q_b}$ · $\text{gap}(L) = o + (L{-}1)e$ · $d_{\text{Pois}} = -\ln(1 - p)$ · $d_{\text{JC}} = -\frac{3}{4}\ln\left(1 - \frac{4}{3}p\right)$ · $Q(i,j) = (n{-}2)d(i,j) - R_i - R_j$ · $d(t) = 2rt$
