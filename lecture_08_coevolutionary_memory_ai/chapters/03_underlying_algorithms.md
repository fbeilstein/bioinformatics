:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Self-Attention · Masked Language Modeling · Contact Emergence · ESMFold
:::

---

# Curriculum

- **From Position-Independent Profiles (PSSMs) to Co-Evolution & Epistasis** — why column-independent models factorize $P(\mathbf{x}) = \prod_i P(x_i)$ and cannot represent epistasis; Potts / Direct Coupling Analysis (DCA) pairwise energy landscapes $E(\mathbf{x})$ connecting evolutionary sequence probabilities to folding free energy ($\Delta\Delta G$) via Boltzmann statistics.
- **Tokenization, Positional Geometry & Contextual Embeddings** — discrete amino-acid vocabulary (`<cls>`, `<eos>`, `<mask>`, 20 standard AAs), embedding lookup $X \in \mathbb{R}^{L \times d}$, and Rotary Positional Embeddings (RoPE in ESM-2) injecting relative distance $(i - j)$ directly into Query–Key inner products.
- **Scaled Dot-Product & Multi-Head Self-Attention** — linear projections $Q = XW_Q, K = XW_K, V = XW_V$; variance scaling by $1/\sqrt{d_k}$ preventing softmax gradient saturation; shift-invariant numerically stable softmax ($\exp(S - \max S)$); multi-head projection concatenation and $O(L^2 d)$ computational complexity.
- **Masked Language Modeling (MLM) & Unsupervised Contact Emergence** — masking 15% of amino acids across UniRef50 and minimizing cross-entropy loss; why predicting masked residues forces attention heads to attend to 3D contacting partners; extracting contact maps via symmetrization and Average Product Correction (APC).
- **Zero-Shot Variant Effect Prediction & Experimental Validation** — masked-marginal vs. wild-type-marginal log-odds scoring $\Delta s = \log P(x_i^{\text{mut}} \mid \mathbf{x}\_{\setminus i}) - \log P(x_i^{\text{wt}} \mid \mathbf{x}\_{\setminus i})$ (cancellation of the softmax partition function $\log Z$); statistical benchmarking against Deep Mutational Scanning (Pearson $r$, Spearman $\rho$).
- **MSA-Free Atomic Structure Prediction (ESMFold vs. AlphaFold 2)** — replacing explicit MSA searches and Evoformer with internal representations of a 3B-parameter ESM-2 language model; pair/single Folding Trunk, equivariant Structure Module outputting 3D coordinates and per-residue pLDDT confidence, and GPU memory chunking / streaming.

---
# Key Math Notes — Co-Evolution, Self-Attention & Stability

:::matrix { cols="52/48" gap="30px" }
[[0,0]]
### PSSMs vs. Epistatic Energy Models
- **Position-Specific Scoring Matrix (PSSM)** assumes independent sites (zero epistasis):
  $$P_{\text{PSSM}}(\mathbf{x}) = \prod_{i=1}^L P_i(x_i)$$
- **Potts / DCA co-evolutionary model** couples contacting pairs $(i, j)$:
  $$E(\mathbf{x}) = -\sum_{i=1}^L h_i(x_i) - \sum_{i < j} J_{ij}(x_i, x_j), \qquad P(\mathbf{x}) = \frac{e^{-E(\mathbf{x})}}{Z}$$
- Under evolutionary selection, surviving sequences approximate a Boltzmann distribution over stability: $\log\frac{P(\mathbf{x}^{\text{mut}})}{P(\mathbf{x}^{\text{wt}})} \propto -\Delta\Delta G$.

### Rotary & Contextual Representations
- Tokenize sequence of length $L$ → $L + 2$ tokens (`<cls>`, $x_1, \dots, x_L$, `<eos>`).
- Each Transformer layer updates contextual embeddings $X \in \mathbb{R}^{(L+2) \times d}$ via multi-head self-attention + feed-forward residual blocks; ESM-2 uses **RoPE** so $q_i^T k_j$ depends on relative separation $(i - j)$.

[[0,1]]
### Scaled Dot-Product & Multi-Head Attention
Given $Q = XW_Q$, $K = XW_K$, $V = XW_V \in \mathbb{R}^{N \times d_k}$:
$$\text{Attention}(Q, K, V) = \operatorname{softmax}\left(\frac{Q K^T}{\sqrt{d_k}}\right) V$$
- **Why divide by $\sqrt{d_k}$?** If $q_{im}, k_{jm} \sim \mathcal{N}(0, 1)$, then $\operatorname{Var}(\mathbf{q}\_i \cdot \mathbf{k}\_j) = d_k$; scaling by $1/\sqrt{d_k}$ restores unit variance and avoids vanishing softmax gradients.
- **Numerically stable row-wise softmax** ($m_i = \max_j S_{ij}$):
  $$A_{ij} = \frac{\exp(S_{ij} - m_i)}{\sum_{l=1}^N \exp(S_{il} - m_i)} = \frac{e^{-m_i} e^{S_{ij}}}{e^{-m_i}\sum_l e^{S_{il}}}$$
  Prevents `exp` floating-point overflow (`inf`/`NaN`) with zero change in mathematical value.
- **Multi-Head Attention** ($H$ heads, $d_k = d/H$):
  $$\operatorname{MultiHead}(X) = \operatorname{Concat}(\text{head}_1, \dots, \text{head}_H) W_O$$
:::

---

# Key Math Notes — MLM, Contact Emergence & Zero-Shot Fitness

:::matrix { cols="53/47"}
[[0,0]]
### Masked Language Modeling (MLM)
Given masked index set $\mathcal{M} \subset \{1, \dots, L\}$ ($\approx 15\%$ of positions):
$$\mathcal{L}\_{\text{MLM}} = -\sum_{i \in \mathcal{M}} \log P_\theta\bigl(x_i \mid \mathbf{x}\_{\setminus \mathcal{M}}\bigr)$$
- Minimizing $\mathcal{L}_{\text{MLM}}$ across millions of UniRef50 homologs forces attention heads to gather context from 3D contacting residues.

### Contact Extraction & ESMFold Confidence
- Raw attention $A^{(l,h)}$ is asymmetric; contact predictors strip `<cls>`/`<eos>`, symmetrize $\bar{A} = \tfrac{1}{2}(A + A^T)$, and apply **Average Product Correction (APC)**:
$$
  A^{\text{APC}}\_{ij} = \bar{A}\_{ij} - {\bar{A}\_{i\cdot} \bar{A}\_{\cdot j}}/{\bar{A}\_{\cdot\cdot}}
$$
before logistic regression across all $(l, h)$ heads (long-range contacts: $\lvert i - j \rvert > 4$, $C\beta$–$C\beta < 8\text{ Å}$).
- **ESMFold** feeds 3B-parameter ESM-2 representations into an equivariant Folding Trunk + Structure Module, outputting 3D coordinates and per-residue **pLDDT** $\in [0, 100]$.

[[0,1]]
### Zero-Shot Variant Effect Scoring
For a single substitution $x_i^{\text{wt}} \to x_i^{\text{mut}}$ at position $i$, **masked-marginal** log-odds is:
$$\Delta s(x_i^{\text{wt}} \to x_i^{\text{mut}}) = \log P_\theta(x_i^{\text{mut}} \mid \mathbf{x}\_{\setminus i}) - \log P_\theta(x_i^{\text{wt}} \mid \mathbf{x}\_{\setminus i})$$
- **Logit cancellation**: since $P_\theta(a \mid \cdot) = e^{z_{i,a}} / Z_i$, the partition function cancels identically:
  $$\Delta s = \bigl(z_{i,\text{mut}} - \log Z_i\bigr) - \bigl(z_{i,\text{wt}} - \log Z_i\bigr) = z_{i,\text{mut}} - z_{i,\text{wt}}$$
- **Masked vs. WT marginal**: replacing site $i$ with `<mask>` avoids self-conditioning on the WT token; unmasked WT-marginal gives a faster single-pass approximation.

### Experimental DMS Validation
Compare zero-shot scores $\mathbf{s}$ against experimental fitness $\mathbf{y}$ (e.g. GFP fluorescence) via Pearson correlation:
$$r = \frac{\sum_k (y_k - \bar{y})(s_k - \bar{s})}{\sqrt{\sum_k (y_k - \bar{y})^2 \sum_k (s_k - \bar{s})^2}}$$

:::
