:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Potts Epistasis · RoPE · Self-Attention · MLM & APC · Zero-Shot Fitness · ESMFold
:::

---

# Curriculum & Algorithmic Roadmap

Six statistical-physics and deep-learning algorithms turn 3.8 billion years of unlabelled protein sequences (`UniRef50`) into unsupervised 3D contact maps, zero-shot mutational fitness predictors, and single-sequence atomic folding in Part 4:

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "seq",  "label": "Evolutionary DB\n(UniRef50 ~65M)",         "x": 20,   "y": 40, "width": 150, "height": 64 },
    { "id": "tok",  "label": "Tokenization +\nRoPE Relative Pos.",      "x": 205,  "y": 40, "width": 165, "height": 64, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "attn", "label": "Multi-Head Attn\nsoftmax(QKᵀ/√d_k)V",      "x": 405,  "y": 40, "width": 175, "height": 64, "style": "fill: #2563EB; color: #fff;" },
    { "id": "mlm",  "label": "MLM 15% Mask +\nAPC Contact Maps",         "x": 615,  "y": 40, "width": 165, "height": 64, "style": "fill: #16A34A; color: #fff;" },
    { "id": "fit",  "label": "Zero-Shot Log-Odds\n& DMS Validation",     "x": 815,  "y": 40, "width": 165, "height": 64, "style": "fill: #D97706; color: #fff;" },
    { "id": "fold", "label": "ESMFold Trunk +\nSE(3) IPA & pLDDT",       "x": 1015, "y": 40, "width": 170, "height": 64, "style": "fill: #DC2626; color: #fff;" }
  ],
  "links": [
    { "from": "seq",  "to": "tok"  },
    { "from": "tok",  "to": "attn" },
    { "from": "attn", "to": "mlm"  },
    { "from": "mlm",  "to": "fit"  },
    { "from": "fit",  "to": "fold" }
  ]
}
:::

:::matrix { .cards cols="1fr 1fr 1fr" gap="18px" }
[[0,0]]{ .purple }
### 1. Epistasis & Positional Geometry
Why column-independent PSSMs miss epistasis, Potts / Direct Coupling energy landscapes $E(\mathbf{x})$, and Rotary Positional Embeddings (RoPE) in ESM-2.

[[0,1]]{ .blue }
### 2. Self-Attention & MLM Emergence
$1/\sqrt{d_k}$ variance scaling, max-shifted softmax overflow protection, multi-head specialization, and Average Product Correction (APC) contact denoising.

[[0,2]]{ .green }
### 3. Zero-Shot Fitness & ESMFold
Masked-marginal log-odds ($\log Z$ cancellation), DMS Pearson $r$ vs. Spearman $\rho$, Triangular pair updates, equivariant IPA, and chunked GPU folding.
:::

---

# Why Position-Independent Profiles (PSSMs) Miss Epistasis

Traditional sequence profiles score each column independently — an assumption that breaks down whenever two residues physically touch inside a folded protein:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Factorized PSSM (Zero Epistasis)
- A Position-Specific Scoring Matrix (PSSM) assumes every sequence site is statistically independent:
$$
P_{\text{PSSM}}(\mathbf{x}) = \prod_{i=1}^L P_i(x_i)
$$
- Taking logs forces double-mutant effects to be strictly **additive**:
$$
\Delta s_{\text{PSSM}}(AB) = \Delta s_{\text{PSSM}}(A) + \Delta s_{\text{PSSM}}(B)
$$
- **Physical failure:** If mutating residue $i$ (`Gly → Trp`) causes a lethal steric clash that a second mutation at contacting residue $j$ (`Leu → Ala`) rescues, a PSSM erroneously predicts an even *worse* double clash!

[[0,1]]
### Potts / DCA Pairwise Energy Landscapes
- **Direct Coupling Analysis (DCA)** introduces a coupled Hamiltonian over sequence $\mathbf{x} = (x_1, \dots, x_L)$:
$$
E(\mathbf{x}) = -\sum_{i=1}^L h_i(x_i) - \sum_{i < j} J_{ij}(x_i, x_j), \qquad P(\mathbf{x}) = \frac{e^{-E(\mathbf{x})}}{Z}
$$
  - $h_i(x_i)$: local site conservation (single-column bias).
  - $J_{ij}(x_i, x_j)$: direct pairwise epistatic coupling between sites $(i, j)$.
- **Evolutionary Boltzmann bridge:** Under purifying selection for fold stability, natural sequence frequencies approximate a Boltzmann ensemble over folding free energy $\Delta G$:
$$
\log \frac{P(\mathbf{x}^{\text{mut}})}{P(\mathbf{x}^{\text{wt}})} = -\bigl(E(\mathbf{x}^{\text{mut}}) - E(\mathbf{x}^{\text{wt}})\bigr) \propto -\Delta\Delta G_{\text{fold}}
$$
:::

---

# Direct vs. Indirect Couplings & Evolutionary Statistical Physics

![PSSM Additivity Failure vs Potts Compensatory Epistasis and Direct vs Indirect Couplings](lecture_08_coevolutionary_memory_ai/images/epistasis_and_dca_couplings.png){ width=94% center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Compensatory double-mutant rescue:** When two contacting sites carry strong coupling $J_{ij} > 0$, flipping either site alone breaks the interaction ($\Delta s < 0$), but flipping **both** restores packing — creating large positive epistasis $\Delta\Delta s_{\text{epi}} = \Delta s(AB) - \Delta s(A) - \Delta s(B) > 0$ invisible to PSSMs!

[[0,1]]
- **Why pairwise Mutual Information ($MI_{ij}$) fails:** If site 1 touches site 2 and site 2 touches site 3, correlation propagates transitively ($1 \leftrightarrow 2 \leftrightarrow 3$), giving a false $MI_{13} > 0$ peak even though 1 and 3 never touch! Global Potts couplings $J_{ij}$ (and deep Transformers trained across **all** families jointly) isolate direct 3D contacts!
:::

---

# Tokenization, Special Tokens & Rotary Positional Embeddings (RoPE)

Before self-attention can operate on a protein sequence of length $L$, discrete amino acid letters must enter a continuous geometric vector space that knows exact relative separation $(i - j)$:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Vocabulary & the `<cls>` Index Shift
- ESM-2 uses a 33-token vocabulary: 20 canonical amino acids + rare codes (`U`, `O`, `B`, `Z`, `X`) + special tokens (`<cls>`, `<eos>`, `<mask>`, `<pad>`, `<unk>`).
- Wrapping a sequence of length $L$ produces $L + 2$ token IDs mapped via embedding lookup to $X \in \mathbb{R}^{(L+2) \times d}$:
  - Token index `0` = **`<cls>`** (sequence-level summary slot)
  - Token index `i + 1` = **0-based residue `i`** (1-based residue `i + 1`!)
  - Token index `L + 1` = **`<eos>`** (end of chain)
- **Accounting rule:** Strip special tokens (`[1:-1, 1:-1]`) before mapping attention matrices onto 0-based residue coordinates!

[[0,1]]
### Rotary Positional Embeddings (RoPE in ESM-2)
- Self-attention is naturally permutation-equivariant (a bag of residues!). Instead of adding fixed absolute vectors, ESM-2 rotates each 2D pair $(2m, 2m+1)$ of Query and Key channels by angle $i \theta_m$ ($\theta_m = 10000^{-2m/d_k}$):
$$
\tilde{\mathbf{q}}\_i = R_i \mathbf{q}\_i, \qquad \tilde{\mathbf{k}}\_j = R_j \mathbf{k}\_j
$$
- Because 2D rotation matrices obey $R_i^T R_j = R_{j - i}$, their dot product depends **strictly on relative separation $(j - i)$**:
$$
(R_{i + \Delta} \mathbf{q}\_i)^T (R_{j + \Delta} \mathbf{k}\_j) = \mathbf{q}\_i^T R_{j - i} \mathbf{k}\_j
$$
  letting ESM-2 generalize cleanly to protein lengths never seen during pretraining!
:::

---
# Interactive Lab — Potts Epistasis & Rotary Positional Embeddings (RoPE)

Explore co-evolutionary physics and positional geometry live: **(1)** in **Tab 1**, drag coupling **$J$** from `0.0` (exact PSSM limit, $\Delta\Delta s_{\text{epi}} = 0.00$) up to `3.0` to watch compensatory double-mutant rescue emerge while pairwise $MI_{13}$ falsely spikes on the indirect chain `(1↔2↔3)` where direct coupling $J_{13} = 0$; **(2)** in **Tab 2**, slide the **Global Sequence Shift $+\Delta$** to verify `0.00e+00` error in RoPE's relative inner product!

<iframe src="lecture_08_coevolutionary_memory_ai/demo/epistasis_and_rope_lab.html" title="Interactive Potts Epistasis and Rotary Positional Embeddings Lab" style="width: 100%; height: 580px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Scaled Dot-Product Attention & Variance Scaling ($1/\sqrt{d_k}$)

Every Transformer layer updates contextual residue representations by allowing every residue $i$ to query all other residues $j \in \{1, \dots, N\}$ simultaneously:

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Linear Projections & Scaled Attention
Project input states $X \in \mathbb{R}^{N \times d}$ into Queries, Keys, and Values via learned matrices $W_Q, W_K, W_V \in \mathbb{R}^{d \times d_k}$:
$$
Q = X W_Q, \quad K = X W_K, \quad V = X W_V
$$
$$
\text{Attention}(Q, K, V) = \operatorname{softmax}\left(\frac{Q K^T}{\sqrt{d_k}}\right) V
$$
- **Biophysical meaning of $S_{ij} = \mathbf{q}\_i \cdot \mathbf{k}\_j / \sqrt{d_k}$:** Query $\mathbf{q}\_i$ encodes what biochemical partner residue $i$ needs (e.g., opposite charge, hydrophobic pocket, H-bond acceptor); Key $\mathbf{k}\_j$ advertises what residue $j$ offers and its RoPE separation. High dot product routes Value $\mathbf{v}\_j$ directly into residue $i$'s updated state!

[[0,1]]
### Exact Variance Proof: Why Divide by $\sqrt{d_k}$?
Suppose the coordinates $q_{im}, k_{jm}$ are independent with mean $0$ and variance $1$. Then the unscaled inner product has variance growing linearly with dimension $d_k$:
$$
\operatorname{Var}(\mathbf{q}\_i \cdot \mathbf{k}\_j) = \sum_{m=1}^{d_k} \operatorname{Var}(q_{im} k_{jm}) = d_k \implies \operatorname{SD}(\mathbf{q}\_i \cdot \mathbf{k}\_j) = \sqrt{d_k}
$$
- At $d_k = 64$, raw dot products span $\pm 16\text{–}24$, driving $\operatorname{softmax}$ into extreme saturated `0.9999 / 0.0001` one-hot vectors where the Jacobian $\frac{\partial A_{ij}}{\partial S_{il}} = A_{ij}(\delta_{jl} - A_{il}) \to 0$ (vanishing gradients!).
- Dividing by $\sqrt{d_k}$ restores $\operatorname{Var}(S_{ij}) = 1$ at **every** head width!
:::

---
# Numerically Stable Softmax & Multi-Head Attention Geometry

![Softmax Variance Scaling Collapse and RoPE Relative Separation Invariance](lecture_08_coevolutionary_memory_ai/images/attention_scaling_and_rope.png){ width=90 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Shift-invariant softmax stability (`m_i = max_j S_ij`):** In IEEE-754 float64, $e^{710}$ overflows to `Infinity`, turning `Inf / Inf` into **`NaN`**! Subtracting row max $m_i$ shifts all exponents $\le 0$ ($e^{S_{ij}-m_i} \in (0, 1]$) with zero mathematical change:
$$
A_{ij} = \frac{\exp(S_{ij} - m_i)}{\sum_{l=1}^N \exp(S_{il} - m_i)} = \frac{e^{-m_i} e^{S_{ij}}}{e^{-m_i} \sum_{l=1}^N e^{S_{il}}}
$$

[[0,1]]
- **Multi-head specialization ($H$ parallel heads, $d_k = d/H$):**
$$
\operatorname{MultiHead}(X) = \operatorname{Concat}(\text{head}\_1, \dots, \text{head}\_H) W_O
$$
  Rather than forcing one attention matrix to track every interaction, `esm2_t6_8M` splits $d = 320$ into **`20` heads ($d_k = 16$)** — letting separate heads specialize in local $i \pm 1$ backbones, $i \to i+4$ $\alpha$-helices, and distant 3D tertiary contacts!
:::

---

# Interactive Lab — Scaled Dot-Product Attention, Stability & Multi-Head Routing

Test transformer core math live: **(1)** in **Tab 1**, sweep head width **$d_k$** (`4 → 256`) with **Divide by $\sqrt{d_k}$** unchecked to watch softmax gradients collapse to `0`, then click **Trigger Overflow (`M = 750`)** to see naive `exp(S)` turn to `NaN` while **Subtract `max(S)`** stays exact; **(2)** in **Tab 2**, switch across **Local**, **$\alpha$-Helix**, and **3D Contact** heads and click tokens to inspect the `<cls>` `+1` index offset!

<iframe src="lecture_08_coevolutionary_memory_ai/demo/self_attention_lab.html" title="Interactive Scaled Dot-Product Attention, Stability and Multi-Head Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Masked Language Modeling (MLM): Objective & Information Routing

ESM-2 is never shown a single crystal structure or $\Delta\Delta G$ measurement during pretraining. How does pure text fill-in-the-blank learning discover 3D macromolecular physics?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### The MLM Objective Across UniRef50
- Sample a protein sequence $\mathbf{x} = (x_1, \dots, x_L)$ from `~65M` evolutionary clusters (UniRef50, reducing clade sampling bias).
- Randomly select a masked index set $\mathcal{M} \subset \{1, \dots, L\}$ covering **`15%`** of positions (`80%` replaced by `<mask>`, `10%` random amino acid, `10%` unchanged so representations stay calibrated at inference).
- Train transformer parameters $\theta$ by minimizing negative log-likelihood (cross-entropy) over masked positions:
$$
\mathcal{L}\_{\text{MLM}}(\theta) = -\sum_{i \in \mathcal{M}} \log P_\theta\bigl(x_i \mid \mathbf{x}\_{\setminus \mathcal{M}}\bigr)
$$

[[0,1]]
### Why 3D Contact Geometry Emerges Unsupervised
- **Information bottleneck:** To predict a hidden residue $x_i$ packed inside a hydrophobic core or catalytic cleft, the model cannot read position $i$ — it must gather evidence from whichever unmasked positions $j \ne i$ have the highest conditional mutual information with $x_i$.
- **Physics sets the signal:** In real evolution, the residues that tightly constrain $x_i$ across distant homologs are its **direct 3D physical neighbors** ($C\beta$–$C\beta < 8\text{ Å}$).
- Minimizing $\mathcal{L}_{\text{MLM}}$ across 36 layers and 3B parameters therefore forces dedicated attention heads to become explicit **3D tertiary contact detectors**!
:::

---
# Unsupervised Contact Extraction: Symmetrization & APC

![Raw Asymmetric Attention vs APC Denoised Contact Map and Precision Jump](lecture_08_coevolutionary_memory_ai/images/mlm_apc_contact_emergence.png){ width=90 center }

:::matrix { cols="1.05fr 0.95fr" gap="26px" }
[[0,0]]
- **Why raw attention maps look streaky:** Each row of $A^{(l,h)}$ sums to $1$ (making $A_{ij} \ne A_{ji}$), and high-entropy or conserved "sink" tokens attract broad background attention across whole rows and columns — masking sparse long-range contacts (`|i - j| > 4`).
- **Two-step geometric denoising (Rao et al., 2021):**
  1. Strip `<cls>`/`<eos>` and **symmetrize**: $\bar{A} = \frac{1}{2}(A + A^T)$.

[[0,1]]
  2. Apply **Average Product Correction (APC)** (borrowed directly from DCA/Potts contact prediction!): $A^{\text{APC}}\_{ij} = \bar{A}\_{ij} - {\bar{A}\_{i\cdot} \bar{A}\_{\cdot j}}/{\bar{A}\_{\cdot\cdot}}$
  subtracting the exact rank-1 outer product of row and column background means!
- A tiny **sparse logistic regression** over all $L_{\text{layers}} \times H$ APC-denoised heads reaches **`> 70%` Top-$L$ long-range precision** in ESM-2 3B — purely from 1D sequence pretraining!
:::

---
# Zero-Shot Variant Effect Scoring: Log-Odds & $\log Z$ Cancellation

Once pretrained via MLM, `EsmForMaskedLM` predicts whether an unseen point mutation $x_i^{\text{wt}} \to x_i^{\text{mut}}$ is neutral, beneficial, or lethal **without training on a single mutational assay** ("zero-shot"):

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Masked-Marginal Log-Odds & Partition Cancellation
Replace site $i$ with `<mask>` so the model scores candidate residues strictly from the surrounding sequence context $\mathbf{x}\_{\setminus i}$:
$$
\Delta s(x_i^{\text{wt}} \to x_i^{\text{mut}}) = \log P_\theta\bigl(x_i^{\text{mut}} \mid \mathbf{x}\_{\setminus i}\bigr) - \log P_\theta\bigl(x_i^{\text{wt}} \mid \mathbf{x}\_{\setminus i}\bigr)
$$
- Since the MLM head computes $P_\theta(a \mid \mathbf{x}\_{\setminus i}) = \exp(z_{i,a}) / Z_i$ where $Z_i = \sum_{b=1}^{33} e^{z_{i,b}}$, the normalizing constant $\log Z_i$ **cancels algebraically**:
$$
\Delta s = (z_{i,\text{mut}} - \log Z_i) - (z_{i,\text{wt}} - \log Z_i) = z_{i,\text{mut}} - z_{i,\text{wt}}
$$
  Log-odds probability difference equals raw logit difference identically!

[[0,1]]
### Three Practical Zero-Shot Scoring Modes
1. **Masked-Marginal (Practice 8 standard):** Mask position $i$ individually (`L` forward passes to score all $19L$ single mutants). Unbiased — avoids self-conditioning on the wild-type token!
2. **Wild-Type Marginal (Ultra-fast):** Run **1** unmasked forward pass on $\mathbf{x}^{\text{wt}}$ and read $z_{i,\text{mut}} - z_{i,\text{wt}}$ everywhere. Fast for proteome-scale screens, slightly biased toward WT.
3. **Multi-Mutant / Epistatic Scoring:** Mask *both* positions $\{i, j\}$ (or compute autoregressive pseudolikelihood $\sum_k \log P(x_k \mid \mathbf{x}\_{\setminus k})$) to capture higher-order compensatory epistasis!
:::

---
# Benchmarking Zero-Shot AI Against Deep Mutational Scanning (DMS)

![Masked vs Wild-Type Marginal Probabilities and DMS Fluorescence Cliff Correlation](lecture_08_coevolutionary_memory_ai/images/zeroshot_marginal_and_dms_metrics.png){ width=87 center }

:::matrix { cols="46/54" gap="5px" }
[[0,0]]
- **Validating evolutionary likelihood against real physics:** In **Deep Mutational Scanning** (e.g., Sarkisyan et al. *Nature* 2016: `> 50,000` GFP variants sorted by FACS-seq), fluorescence exhibits a cooperative **sigmoidal cliff** — mildly destabilizing mutations stay bright until $\Delta G$ crosses the unfolding threshold, then crash flat to `0` (dark).

[[0,1]]
- **Pearson $r$ vs. Spearman $\rho$:**
  - **Pearson $r$** measures strict *linear* proportionality between fitness $\mathbf{y}$ and log-odds $\mathbf{s}$ — penalized whenever dark/bright experimental plateaus curve away from a straight line!
  - **Spearman $\rho$** (Pearson on **ranks**) captures any monotonic ordering — making $\rho$ the gold standard on ProteinGym DMS leaderboards (`> 200` assays, 2.5M+ variants) and clinical VUS triage!
:::

---
# Interactive Lab — APC Contact Denoising, Zero-Shot Log-Odds & DMS Metrics

Experiment with unsupervised contact extraction and mutational fitness live: **(1)** in **Tab 1**, ramp **Entropy Streak Noise $\alpha$** up to `1.0` and watch **Raw Attention** lose long-range contacts while **Symmetrization + APC** recovers sharp Top-10 precision; **(2)** in **Tab 2**, drag the **Constant Logit Offset $+c$** to confirm `0.00e+00` error in $\log Z$ cancellation, click any mutant amino acid bar, and dial **DMS Assay Saturation** from `0 → 1` to watch Spearman $\rho$ outperform Pearson $r$!

<iframe src="lecture_08_coevolutionary_memory_ai/demo/mlm_apc_and_fitness_lab.html" title="Interactive APC Contact Denoising, Zero-Shot Log-Odds and DMS Lab" style="width: 100%; height: 580px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# From Contacts to 3D Coordinates: Folding Trunk & Triangular Updates

Attention contact maps give noisy 2D probabilities — how do **AlphaFold 2** and **ESMFold** refine them into physically self-consistent 3D atomic coordinates?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### Coupled Single & Pair Representations
- Recall from **Lecture 7** that pairwise distances $D_{ij}$ are $\text{SE}(3)$-invariant, yet arbitrary predicted distances can violate 3D Euclidean geometry (e.g., $d_{13} > d_{12} + d_{23}$).
- Modern folding trunks evolve two coupled tensors jointly across 48 blocks:
  - **Single representation** $\mathbf{s}\_i \in \mathbb{R}^{c_s}$ ($O(L)$): per-residue evolutionary and local chemical state.
  - **Pair representation** $\mathbf{z}\_{ij} \in \mathbb{R}^{c_z}$ ($O(L^2)$): latent relational geometry between every residue pair $(i, j)$.

[[0,1]]
### Triangular Multiplicative Updates ($O(L^3)$)
To enforce 3D Euclidean consistency across every triplet of residues $(i, j, k)$, each trunk block updates edge $(i, j)$ by summing gated products along all third-party bridges $k$:
$$
\Delta \mathbf{z}\_{ij} = \sum_{k=1}^L \mathbf{a}\_{ik} \odot \mathbf{b}\_{jk} \quad \text{(outgoing & incoming triangles)}
$$
- **Geometric intuition:** If residue $i$ strongly contacts $k$ and $k$ strongly contacts $j$, triangle updates propagate geometric constraints onto edge $(i, j)$ — resolving transitivity and enforcing the 3D triangle inequality before a single atom is placed!
:::

---

# MSA-Free Atomic Folding: ESMFold Architecture & Equivariant IPA

**ESMFold** (Lin et al., *Science* 2023) replaces AlphaFold 2's slow 100 GB database MSA search and Evoformer with the implicit co-evolutionary memory already stored inside a frozen **3B-parameter ESM-2** language model:

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "seq",   "label": "Single Sequence\n(Zero MSA Search!)",     "x": 20,   "y": 34, "width": 150, "height": 62 },
    { "id": "lm",    "label": "Frozen ESM-2 (3B)\n36 Layers · Implicit MSA","x": 205,  "y": 34, "width": 175, "height": 62, "style": "fill: #7C3AED; color: #fff;" },
    { "id": "proj",  "label": "Layer-Weighted Sum\nSingle s_i + Pair z_ij","x": 415,  "y": 34, "width": 175, "height": 62, "style": "fill: #2563EB; color: #fff;" },
    { "id": "trunk", "label": "Folding Trunk (48×)\nTriangle Pair Updates", "x": 625,  "y": 34, "width": 170, "height": 62, "style": "fill: #16A34A; color: #fff;" },
    { "id": "ipa",   "label": "Equivariant IPA\nSE(3) Residue Frames",      "x": 830,  "y": 34, "width": 165, "height": 62, "style": "fill: #D97706; color: #fff;" },
    { "id": "out",   "label": "3D Atomic PDB +\npLDDT (0–100) + PAE",       "x": 1030, "y": 34, "width": 165, "height": 62, "style": "fill: #DC2626; color: #fff;" }
  ],
  "links": [
    { "from": "seq",   "to": "lm"    },
    { "from": "lm",    "to": "proj"  },
    { "from": "proj",  "to": "trunk" },
    { "from": "trunk", "to": "ipa"   },
    { "from": "ipa",   "to": "out"   }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### $\text{SE}(3)$-Equivariant Invariant Point Attention (IPA)
- Each residue $i$ carries a rigid local backbone frame $T_i = (R_i, \mathbf{t}\_i) \in \text{SE}(3)$ centered at $C\alpha_i$ (the exact $\text{SE}(3)$ group from Lecture 7!).
- **IPA** transforms learned 3D query/key points into global space via $T_i$, computes **$\text{SE}(3)$-invariant** Euclidean distances $\lVert T_i(\mathbf{q}\_i) - T_j(\mathbf{k}\_j) \rVert_2^2$, and iteratively refines frame updates $(R_i, \mathbf{t}\_i)$ + side-chain torsions $\chi_1 \dots \chi_4$!

[[0,1]]
### FAPE Loss & Per-Residue Confidence (`pLDDT`)
- **Frame Aligned Point Error (FAPE):** Evaluates every predicted atom inside *every* residue's local frame $T_i^{-1}(\mathbf{x}\_j)$ — automatically $\text{SE}(3)$-invariant yet **chiral** (penalizes mirror L$\to$D reflections, solving MDS's mirror ambiguity!).
- **pLDDT $\in [0, 100]$:** Predicted local distance accuracy (`> 90` atomic accuracy, `70–90` confident fold, `< 50` intrinsically disordered flexible loop!).
:::

---
# ESMFold vs. AlphaFold 2: Scaling, Memory Chunking & Physical Limits

![ESMFold GPU Memory Chunking and AlphaFold 2 vs ESMFold MSA Depth Trade-Off](lecture_08_coevolutionary_memory_ai/images/esmfold_vs_af2_and_scaling.png){ width=80 center }

:::matrix { cols="47/53" gap="5px" }
[[0,0]]
- **Why chunking & `accelerate` matter on Colab T4 GPUs:**
  - Unchunked triangular updates scale as **$O(L^3)$ memory**, crashing a `16 GB` T4 GPU by `L ≈ 440` residues!
  - Setting `trunk.set_chunk_size(64)` tiles the third loop into slices of 64, reducing peak activation VRAM from $O(L^3)$ to $O(L^2)$, while Hugging Face `accelerate` streams 3B weights directly from disk to GPU VRAM (avoiding CPU RAM OOM).

[[0,1]]
- **Accuracy vs. speed trade-off & physical blind spots:**
  - **`60×` speedup** (`~3 sec` vs. `10–30 min` MSA search) enabled folding **`> 617M` metagenomic proteins** (ESM Metagenomic Atlas) and beats AF2 on low-homology orphan proteins!
  - **Physical limits:** Single-sequence models still miss explicit **post-translational modifications** (e.g., GFP `1EMA` chromophore cyclization), co-factors/ligands, and multi-state conformational switches.
:::

---

# Beyond Single-Sequence pLMs: All-Atom Diffusion & Co-Folding

How does the newest generation of structural AI overcome single-chain protein-only limits (like `6LU7 + N3` inhibitor from Lecture 7 or GFP's chromophore)?

:::matrix { cols="1fr 1fr" gap="28px" }
[[0,0]]
### All-Atom Generative Diffusion (AlphaFold 3 / Boltz / RFAA)
- Replaces rigid residue-only IPA frames with an **all-atom 3D denoising diffusion** module predicting atomic coordinates directly from Gaussian noise:
$$
\mathbf{x}\_t = \sqrt{\bar{\alpha}\_t} \mathbf{x}\_0 + \sqrt{1 - \bar{\alpha}\_t} \boldsymbol{\epsilon}, \qquad \hat{\mathbf{x}}\_0 = f_\theta(\mathbf{x}\_t, t, \mathbf{s}, \mathbf{z})
$$
- **Universal biomolecular co-folding:** Treats proteins, DNA, RNA, metal ions, covalent modifications, and drug-like small molecules on an equal footing — predicting protein–ligand complexes (such as protease + inhibitor) directly in one pass!

[[0,1]]
### Multimodal Protein Language Models (ESM3)
- Extends masked language modeling beyond 1D amino-acid tokens by jointly tokenizing **sequence**, **3D backbone structure** (VQ-VAE geometric tokens), and **functional annotations** into one unified discrete stream.
- **Programmable protein design:** Instead of only predicting *sequence $\to$ structure*, multimodal pLMs condition on desired active-site 3D geometry and function to generate novel sequences (e.g., designing distant functional GFP analogs with `< 60%` sequence identity to any natural protein!).
:::

---

# Algorithmic Complexity & Scaling of Protein AI Models

Understanding where quadratic $O(L^2)$ vs. cubic $O(L^3)$ bottlenecks live explains every engineering choice in Part 4 (`esm2_t6_8M` on CPU vs. chunked `ESMFold` on T4 GPU):

:::matrix {cols=100}
[[0,0]]{ .dense }
| Algorithm / Stage | Mathematical Core | Time Complexity | Peak Memory | Practical Bottleneck & Speedup |
| :--- | :--- | :---: | :---: | :--- |
| **PSSM Profile** | Independent column log-odds | $O(L)$ | $O(20 L)$ | Fast, zero epistasis ($J_{ij} = 0$) |
| **Potts / DCA** | Pairwise couplings $J_{ij}(a, b)$ | $O(L^2 q^2)$ | $O(L^2 q^2)$ | Needs deep MSA per family ($q = 21$) |
| **Self-Attention** | $\operatorname{softmax}(Q K^T / \sqrt{d_k}) V$ | $O(L^2 d)$ | $O(L^2)$ (`O(L)` Flash) | FlashAttention tiles $Q, K, V$ in GPU SRAM |
| **WT-Marginal DMS** | 1 unmasked MLM pass ($19L$ muts) | $O(L^2 d)$ | $O(L^2)$ | Fast screen; slight WT self-bias |
| **Masked-Marginal DMS** | $L$ masked MLM passes (`<mask>`) | $O(L^3 d)$ | $O(L^2)$ | Unbiased zero-shot standard (Practice 8) |
| **Triangular Updates** | $\sum_k \mathbf{a}\_{ik} \odot \mathbf{b}\_{jk}$ (Trunk) | $O(L^3 c_z)$ | $O(L^3)$ (`O(L^2)` chunk) | `set_chunk_size(64)` prevents T4 OOM |
:::

:::matrix { .cards cols="1fr 1fr" gap="18px" }
[[0,0]]{ .blue }
### Why `esm2_t6_8M` Runs in Milliseconds on CPU
With only `6` layers, `d = 320`, and `L = 35` (`gfp_snippet`), `L²` is just `1,225` entries per head — whereas scaling to `3B` parameters + 48 folding trunk blocks (`L = 500`) demands fp16 GPU acceleration and memory chunking!

[[0,1]]{ .purple }
### Unifying the Entire Course in One Pipeline
Every algorithmic theme across Lectures 1–8 converges here: $k$-mer / MSA homology (L1–L2), Bayesian log-odds (L3), statistical shrinkage & rank tests (L4–L6), $\text{SE}(3)$ frames & contact maps (L7), and Transformer co-evolutionary memory (L8)!
:::

---

# Summary & Bridge to Part 4 — Grand Algorithmic Synthesis

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### Epistasis & Geometry
#### Foundations
- PSSMs factorize sites (zero epistasis)
- Potts $J_{ij}$ isolates direct 3D contacts
- Boltzmann bridge: $\log\frac{P_{\text{mut}}}{P_{\text{wt}}} \propto -\Delta\Delta G$
- `<cls>` shifts 0-based indices by `+1`
- RoPE encodes exact separation $(j - i)$

[[0,1]]{ .blue }
### Attention & MLM
#### Emergence
- $1/\sqrt{d_k}$ keeps softmax gradients alive
- $\exp(S - \max S)$ prevents `inf / NaN`
- Multi-head routing specializes by range
- MLM 15% masking forces 3D attention
- Symmetrization + APC reveals contacts

[[0,2]]{ .green }
### Zero-Shot & Folding
#### Prediction
- Masked-marginal: $\Delta s = z_{\text{mut}} - z_{\text{wt}}$
- Spearman $\rho$ handles DMS cliffs
- Triangle updates enforce 3D geometry
- Equivariant IPA predicts $\text{SE}(3)$ frames
- **Part 4:** `ESM-2` & `ESMFold` in PyTorch
:::

> **Core Equations Bridge:** $E(\mathbf{x}) = -\sum h_i - \sum J_{ij}$ · $A = \operatorname{softmax}(Q K^T / \sqrt{d_k})$ · $A^{\text{APC}}\_{ij} = \bar{A}\_{ij} - \bar{A}\_{i\cdot}\bar{A}\_{\cdot j}/\bar{A}\_{\cdot\cdot}$ · $\Delta s = z_{i,\text{mut}} - z_{i,\text{wt}}$ · $\Delta\mathbf{z}\_{ij} = \sum_k \mathbf{a}\_{ik}\odot\mathbf{b}\_{jk}$ — next in **Part 4**, load Meta's `ESM-2` and `ESMFold` directly in PyTorch & Hugging Face to extract contact maps, score GFP mutations against Sarkisyan DMS fluorescence, and fold proteins directly from single sequences!
