:::titlepage
[[top]]
Part 4

[[title]]
Industrial Tools

[[subcaption]]
Hugging Face · Meta ESM-2 · ESMFold · Deep Mutational Scanning
:::

---

# Curriculum

- **Protein Language Models (pLMs)** — evolutionary scale modeling, tokenizer vocabulary, special tokens (`<cls>`, `<eos>`), contextual sequence embeddings.
- **The Self-Attention Mechanism** — Query, Key, Value projections, scaled dot-product attention $\text{softmax}(QK^T / \sqrt{d_k})V$, multi-head architecture.
- **Unsupervised Structural Emergence** — how attention heads spontaneously learn 3D contact maps and secondary structure from 1D evolutionary co-variation.
- **Masked Language Modeling (MLM)** — training objectives, masking probabilities, predicting masked amino acids from sequence context.
- **Zero-Shot Variant Effect Prediction** — log-odds scoring $\log P(\text{mut}) - \log P(\text{wt})$, evolution as a natural selection filter, predicting mutational deleteriousness.
- **Deep Mutational Scanning (DMS)** — experimental fitness landscapes (GFP mutation matrix, Sarkisyan et al.), FACS-seq validation.
- **ESMFold: MSA-Free Structure Prediction** — 3-Billion parameter transformer, bypassing Multiple Sequence Alignment searches, folding proteins directly from single sequences on GPU.

---

# Roadmap

The end-to-end protein language modeling, zero-shot fitness, and folding workflow:

:::static-diagram { width=95% center }
{
  "nodes": [
    { "id": "seq",   "label": "Protein Sequence\n(1D Amino Acids)",       "x": 40,   "y": 80,  "width": 150, "height": 60 },
    { "id": "tok",   "label": "AutoTokenizer\n<cls> + Token IDs",        "x": 220,  "y": 80,  "width": 150, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "esm",   "label": "ESM-2 Transformer\n8M to 3B Parameters",   "x": 400,  "y": 80,  "width": 160, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "attn",  "label": "Self-Attention Heads\nQ, K, V dot-products", "x": 590,  "y": 80,  "width": 160, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "cmap",  "label": "Contact Map Extraction\n3D Structural Contacts", "x": 780, "y": 80, "width": 160, "height": 60, "style": "fill: #FF9800; color: #fff;" },
    { "id": "mlm",   "label": "Masked LM Head\nLogit Probabilities",      "x": 960,  "y": 80,  "width": 150, "height": 60, "style": "fill: #E91E63; color: #fff;" },

    { "id": "dms",   "label": "DMS Fitness Matrix\nGFP Experimental Data", "x": 960,  "y": 250, "width": 160, "height": 60 },
    { "id": "score", "label": "Zero-Shot Scoring\nlog P(mut) - log P(wt)", "x": 770,  "y": 250, "width": 160, "height": 60, "style": "fill: #673AB7; color: #fff;" },
    { "id": "corr",  "label": "Benchmarking AI\nPearson Correlation",     "x": 580,  "y": 250, "width": 160, "height": 60, "style": "fill: #009688; color: #fff;" },
    { "id": "fold",  "label": "ESMFold Architecture\n3B Parameter Trunk",  "x": 390,  "y": 250, "width": 160, "height": 60, "style": "fill: #3F51B5; color: #fff;" },
    { "id": "gpu",   "label": "Accelerate Streaming\nGPU Memory Bypassing","x": 200,  "y": 250, "width": 160, "height": 60, "style": "fill: #795548; color: #fff;" },
    { "id": "pdb",   "label": "Predicted 3D PDB\nWebGL Superposition",    "x": 10,   "y": 250, "width": 160, "height": 60, "style": "fill: #4CAF50; color: #fff;" }
  ],
  "links": [
    { "from": "seq",   "to": "tok" },
    { "from": "tok",   "to": "esm" },
    { "from": "esm",   "to": "attn" },
    { "from": "attn",  "to": "cmap" },
    { "from": "cmap",  "to": "mlm" },
    { "from": "mlm",   "to": "dms" },
    { "from": "dms",   "to": "score" },
    { "from": "score", "to": "corr" },
    { "from": "corr",  "to": "fold" },
    { "from": "fold",  "to": "gpu" },
    { "from": "gpu",   "to": "pdb" }
  ]
}
:::

> In this final lecture, we explore how neural networks trained purely on the "grammar" of evolution implicitly reconstruct 3D macromolecular physics, zero-shot fitness landscapes, and atomic structures.

---

# The Protein AI Toolbox

Modern deep learning for proteins bridges transformer language models, tensor manipulation, and structural validation:

:::matrix { .cards cols="1fr 1fr 1fr 1fr" gap="20px" margin-top: 20px; }
[[0,0]]{ .purple }
### Hugging Face Transformers
Pretrained model hub and pipeline architecture: loads Meta's ESM-2 models, manages tokenizers, and exposes internal attention heads.

[[0,1]]{ .blue }
### PyTorch
Tensor algebra, GPU hardware acceleration, automatic differentiation, and memory streaming via Hugging Face `accelerate`.

[[0,2]]{ .green }
### Meta ESM-2 & ESMFold
Evolutionary Scale Modeling: 8M to 15B parameter protein language models and the 3B-parameter MSA-free atomic folding engine.

[[0,3]]{ .orange }
### SciPy & Py3Dmol
Statistical validation against experimental DMS assays (Pearson $r$) and WebGL 3D rendering of predicted atomic coordinates.
:::

```bash
# Installation in Google Colab / virtual environment
pip install -q transformers torch pandas numpy scipy matplotlib biopython py3Dmol
```

---

# The Biological Premise: Co-Evolutionary Memory

Why does statistical sequence modeling learn physical 3D structure?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Law of Compensatory Mutations
- Proteins function as densely packed 3D machines.
- If residue $i$ and residue $j$ are **spatially adjacent in 3D**, they cannot mutate independently!
- **Positive Epistasis**:
  - If a small residue at position $i$ (e.g. Glycine) mutates into a bulky Tryptophan, it causes steric clash.
  - The protein misfolds and the organism perishes **unless** position $j$ simultaneously mutates into a smaller residue to make room!
- Over billions of years, contacting pairs undergo coordinated **compensatory mutations**.

[[0,1]]
### Evolution as a Fitness Filter
- Nature has performed a **3.8-billion-year screening experiment** across trillions of organisms.
- Unviable mutations that destabilize folding or destroy catalysis are eliminated by natural selection.
- The sequences stored in databases (UniProt / BFD) represent the **survivors** of this massive experiment!
:::

> A neural network trained on millions of evolutionary sequences is not learning syntax — it is learning the **free energy landscape of protein folding**.

---
# Tool 1 — Loading Meta's ESM-2 in Transformers

Meta AI's **ESM-2** (Evolutionary Scale Modeling) is the state-of-the-art transformer family for proteins:

```python
from transformers import AutoTokenizer, EsmModel

model_name = "facebook/esm2_t6_8M_UR50D"

# 1. Load the amino acid tokenizer
tokenizer = AutoTokenizer.from_pretrained(model_name)

# 2. Load model with attention output enabled
model = EsmModel.from_pretrained(model_name, output_attentions=True)
model.eval()  # Put in evaluation mode (disables dropout)

print(f"Loaded {model_name} with {sum(p.numel() for p in model.parameters()):,} parameters.")
```

:::matrix { cols="1.3fr 1fr" gap="30px" }
[[0,0]]{ .smaller-table }
### The ESM-2 Hierarchy
| Model Identifier | Layers | Hidden Dim | Heads | Parameters |
| :--- | :--- | :--- | :--- | :--- |
| `esm2_t6_8M_UR50D` | **6** | **320** | **20** | **8 Million** |
| `esm2_t12_35M_UR50D` | 12 | 480 | 20 | 35 Million |
| `esm2_t30_150M_UR50D` | 30 | 640 | 20 | 150 Million |
| `esm2_t33_650M_UR50D` | 33 | 1280 | 20 | 650 Million |
| `esm2_t36_3B_UR50D` | 36 | 2560 | 40 | 3 Billion |

[[0,1]]
### Why Use `esm2_t6_8M` for Teaching?
- Only **~30 MB** in memory.
- Executes in milliseconds on a standard laptop CPU.
- Preserves the identical attention geometry, tokenization, and contact emergence as the 3B parameter model.
:::

---
# Tool 2 — Tokenization & Contextual Embeddings

How do raw amino acid characters enter a transformer neural network?

```python
# A 35-amino acid snippet of Green Fluorescent Protein (GFP)
gfp_snippet = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG"

# Tokenize into PyTorch tensors
inputs = tokenizer(gfp_snippet, return_tensors="pt")
print("Tokens:  ", inputs['input_ids'][0].tolist())
print("Decoded: ", tokenizer.convert_ids_to_tokens(inputs['input_ids'][0]))
```

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Special Tokens
- **`<cls>` (Index 0)**: Prepended to every sequence. Acts as a global "pooler" representing the whole protein.
- **`<eos>` (Index $N+1$)**: Appended to mark the end of the polypeptide chain.
- **Token Vocabulary**: 20 standard amino acids + rare amino acids (`U`, `O`) + special tokens (`<pad>`, `<mask>`, `<unk>`).

[[0,1]]
### Contextualized Hidden States
```python
with torch.no_grad():
    outputs = model(**inputs)
embeddings = outputs.last_hidden_state
print(embeddings.shape)  # torch.Size([1, 37, 320])
```
- Each amino acid is represented by a **320-dimensional vector**.
- Unlike static BLOSUM62 matrices or **Position-Specific Scoring Matrices (PSSMs)** (which score each sequence column independently), Transformer embeddings are **contextual** — conditioning every residue on the entire sequence so epistatic interactions can be captured!
:::

---

# Tool 3 — The Self-Attention Mechanism

The fundamental engine of the transformer is the **Scaled Dot-Product Attention**:

$$\text{Attention}(Q, K, V) = \text{softmax}\left(\frac{QK^T}{\sqrt{d_k}}\right)V$$

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Three Projections
For input embedding matrix $X \in \mathbb{R}^{N \times d}$:
1. **Query ($Q = X W_Q$)**: What each amino acid is looking for.
2. **Key ($K = X W_K$)**: What each amino acid broadcasts about its chemical properties and position.
3. **Value ($V = X W_V$)**: The information each amino acid routes across the network.
- Projection matrices: $W_Q, W_K, W_V \in \mathbb{R}^{d \times d_k}$.

[[0,1]]
### The Dot-Product Matrix ($Q K^T$)
- Pairwise dot product measures compatibility:
  $$S_{ij} = \frac{\mathbf{q}\_i \cdot \mathbf{k}\_j}{\sqrt{d_k}}$$
- If residue $i$ and residue $j$ have complementary chemical requirements, $\mathbf{q}\_i \cdot \mathbf{k}\_j$ is large and positive.
- Scaling factor $1/\sqrt{d_k}$ prevents gradient saturation in high dimensions.
:::

---
# Multi-Head Attention Architecture

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
In ESM-2 (`t6_8M`), each of the 6 layers contains **20 parallel attention heads**:

```text
Input Sequence Tokens (Length N)
                 │
  ┌──────────────┼──────────────┬──────────────┐
  ▼              ▼              ▼              ▼
Head 1         Head 2         Head 3  ...    Head 20
Q₁, K₁, V₁     Q₂, K₂, V₂     Q₃, K₃, V₃     Q₂₀, K₂₀, V₂₀
  │              │              │              │
  ▼              ▼              ▼              ▼
Attn Map 1     Attn Map 2     Attn Map 3     Attn Map 20
(N × N)        (N × N)        (N × N)        (N × N)
  │              │              │              │
  └──────────────┴──────────────┼──────────────┘
                                ▼
                   Concatenate & Project (W_O)
                                ▼
                 Contextualized Output Vectors
```
[[0,1]]
### Why Multiple Heads?
- Single attention can only look at one relationship per position.
- Multiple heads allow the network to simultaneously attend to:
  - Local backbone neighbors ($i \leftrightarrow i+1$).
  - Hydrogen bonding partners in $\alpha$-helices ($i \leftrightarrow i+4$).
  - Distant tertiary contacts across the hydrophobic core ($i \leftrightarrow j$).

### Tensor Dimensions in PyTorch
```python
attentions = outputs.attentions
# Tuple of 6 layers
layer_6 = attentions[5]  
# Shape: (batch=1, heads=20, seq_len=37, seq_len=37)
```
Each head generates an independent $37 \times 37$ attention weight matrix!
:::

---
# Tool 4 — Unsupervised Emergence of Contact Maps

:::matrix { cols="50/50" gap="30px" }
[[0,0]]
![Attention Contact Map](lecture_08_coevolutionary_memory_ai/images/attention_contact_map.png){width=95% center}

[[0,1]]
What happens when we visualize the self-attention weights of ESM-2?

```python
# Average attention across all 20 heads in Layer 6
mean_attn = outputs.attentions[5][0].mean(dim=0).numpy()
# Strip off <cls> and <eos> special tokens
clean_attn = mean_attn[1:-1, 1:-1]
# Ignore local neighbors (|i - j| <= 4) 
# to find long-range 3D contacts (< 8 Å in 1EMA)
```

> **The Deep Learning Miracle**: ESM-2 was **never shown a 3D structure** during training — only 1D sequence text. Yet its long-range attention peaks (`|i - j| > 4`) map to real physical contacts ($< 8\text{ Å}$ in crystal structure `1EMA`; note GFP's cleaved N-terminal Met1 means `1EMA` starts at residue `2`)!
:::

---

# Why Does Attention Learn 3D Structure?

How can a model learn 3D coordinates from 1D text strings?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Masked Language Modeling Game
- During pretraining, $15\%$ of amino acids are randomly hidden:
  $$\text{M S K G E E } \mathbf{\langle mask \rangle} \text{ T G V V P I L V E} \dots$$
- The transformer must predict the identity of the missing residue.
- To guess correctly, it cannot look at the hidden residue directly.
- It must query all other residues in the sequence that constrain what can fit into that position.

[[0,1]]
### Information Routing Follows Physics
- What constrains an amino acid? **The residues that touch it in 3D space!**
- If position 7 is packed against position 25 inside a $\beta$-barrel, position 25 contains the strongest predictive signal for position 7.
- To minimize cross-entropy loss, the attention mechanism is mathematically forced to assign high attention weights $A_{ij}$ to **spatially contacting pairs**!
:::

---

# Tool 5 — Predicting Mutational Fitness (Zero-Shot)

Can a protein language model predict whether a novel mutation will cause disease or loss of function?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The MLM Output Head (`EsmForMaskedLM`)
- On top of the transformer representations sits a linear projection to vocabulary logits:
  ```python
  from transformers import EsmForMaskedLM
  mlm_model = EsmForMaskedLM.from_pretrained(model_name)
  logits = mlm_model(**inputs).logits
  # Shape: (batch=1, seq_len=37, vocab_size=33)
  ```
- `logits[0, pos, :]` represents the un-normalized log-probability of every amino acid at `pos`.

[[0,1]]
### The Zero-Shot Log-Odds Metric
To score mutation $X_{\text{wt}} \to Y_{\text{mut}}$ at position $i$:
$$\Delta \text{score} = \log P(Y_{\text{mut}} \mid \text{context}) - \log P(X_{\text{wt}} \mid \text{context})$$

$$= \text{logits}[i, \text{token}\_{\text{mut}}] - \text{logits}[i, \text{token}\_{\text{wt}}]$$

- $\Delta \text{score} < 0$: Mutation is **evolutionarily unfavored** (destabilizing / deleterious).
- $\Delta \text{score} \approx 0$: Mutation is **neutral** (tolerated polymorphism).
- $\Delta \text{score} > 0$: Mutation is **beneficial** (rare evolutionary optimization).
:::

---

# Zero-Shot Scoring Function in Python

Implementing the **masked-marginal** variant effect predictor (matching notebook Part 5):

```python
import math, torch

def predict_mutation_effect(wt_sequence, position_1based, mutant_aa):
    """Masked-marginal log-odds score: log P(mut | context) - log P(wt | context)."""
    # 1. Replace target residue with <mask> so prediction uses surrounding context only
    seq_list = list(wt_sequence)
    wt_aa = seq_list[position_1based - 1]
    seq_list[position_1based - 1] = "<mask>"
    inputs = tokenizer("".join(seq_list), return_tensors="pt")

    # 2. Run MLM head; <cls> at index 0 shifts 1-based residue pos → token index pos
    with torch.no_grad():
        logits = mlm_model(**inputs).logits[0, position_1based]
    probs = torch.nn.functional.softmax(logits, dim=0)

    # 3. Log-odds ratio (identically equals logits[mut] - logits[wt] since log Z cancels!)
    wt_id  = tokenizer.convert_tokens_to_ids(wt_aa)
    mut_id = tokenizer.convert_tokens_to_ids(mutant_aa)
    return math.log(probs[mut_id].item()) - math.log(probs[wt_id].item())
```

> **Zero-shot**: requires zero training on assay labels. Masking the site (`<mask>`) prevents the model from seeing the wild-type residue; subtracting raw logits gives the exact same value because the softmax normalizing constant $\log Z$ cancels!

---

# Tool 6 — Deep Mutational Scanning (DMS)

How do experimental biologists benchmark these computational predictions?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Experimental Technology
- **Deep Mutational Scanning (DMS)**:
  - Synthesize a comprehensive DNA library containing thousands of point mutations.
  - Express mutant proteins in cell populations (e.g. yeast or bacteria).
  - Subject cells to a quantitative functional selection assay (e.g. fluorescence-activated cell sorting, FACS).
  - Sequence surviving cells before and after selection to quantify the **relative fitness** of every mutant.

[[0,1]]
### The Sarkisyan GFP Dataset
- Sarkisyan et al. (*Nature* 2016) generated over **50,000 variants** of Green Fluorescent Protein (GFP).
- Measured true physical fluorescence for each mutant:
  - $\text{Fitness} \approx 1.0$: Bright, folded green protein.
  - $\text{Fitness} \le 0.0$: Completely dark, misfolded or aggregated protein.
:::

---

# The Deep Mutational Scanning Landscape

Visualizing the functional impact of every possible amino acid substitution:

![DMS Matrix](lecture_08_coevolutionary_memory_ai/images/dms_mutational_matrix.png){width=80% center}

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Conserved Structural Cores (Red Columns)
- Buried hydrophobic residues and the catalytic fluorophore core (`Ser65–Tyr66–Gly67`).
- Almost **any** substitution is lethal ($\text{Fitness} \ll 0$), causing fluorescence extinction.

[[0,1]]
### Tolerant Surface Loops (White / Blue Columns)
- Solvent-exposed flexible residues.
- Mutating to diverse amino acids preserves wild-type fluorescence ($\text{Fitness} \approx 0$ or higher).
:::

---

# Benchmarking ESM-2 Against Experimental Data

:::matrix { cols="1.05fr 1fr" gap="30px" }
[[0,0]]
![Correlation Scatter](lecture_08_coevolutionary_memory_ai/images/zero_shot_dms_correlation.png){width=95% center}

[[0,1]]
When we correlate zero-shot ESM-2 log-odds scores against experimental DMS fluorescence:

### Strong Experimental Correlation
- Pearson correlation **$r \approx 0.65 - 0.75$** across comprehensive DMS datasets!
- Destabilizing mutations consistently receive deeply negative language model scores.

### Clinical & Industrial Significance
- Enables **in silico screening** of millions of patient variants of unknown significance (VUS).
- Guides protein engineering: introduces stabilizing mutations while avoiding lethal misfolding traps.
:::

---
# Tool 7 — ESMFold: MSA-Free 3D Structure Prediction

In 2022, Meta developed **ESMFold**: predicting full 3D atomic coordinates directly from a single sequence:

### AlphaFold 2 vs. ESMFold

| Feature | AlphaFold 2 | ESMFold |
| :--- | :--- | :--- |
| **Input Requirement** | Multiple Sequence Alignment (MSA) | **Single Sequence Only** |
| **Upstream Search** | Searches 100 GB genetic DBs (10–30 min) | **Zero database search** |
| **Prediction Speed** | Minutes per protein | **Seconds on GPU** ($60\times$ faster!) |
| **Memory Basis** | Explicit MSA alignment matrix | **Implicit weights of 3B pLM** |

### The Architecture
- Replaces AlphaFold's MSA Evoformer module with the **ESM-2 3-Billion parameter language model**.
- Language model embeddings feed into a structural trunk that outputs:
  - 3D Cartesian coordinates ($X, Y, Z$).
  - Per-residue confidence metric (**pLDDT**, $0 \to 100$).

---

# Predicting Structures with the ESMFold Web API

Meta provides a public REST API (`api.esmatlas.com`) for instant single-sequence structure prediction:

```python
import urllib.request, py3Dmol

# Fast sequence submission to the ESMFold inference server
seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

req = urllib.request.Request(
    "https://api.esmatlas.com/foldSequence/v1/pdb/",
    data=seq.encode("utf-8"),
    headers={"Content-Type": "text/plain"}
)

with urllib.request.urlopen(req) as response:
    predicted_pdb = response.read().decode("utf-8")

# Visualize true X-ray (left) vs. ESMFold prediction (right) with viewergrid
view = py3Dmol.view(width=900, height=450, viewergrid=(1, 2))
view.addModel(true_pdb, "pdb", viewer=(0, 0))
view.setStyle({'cartoon': {'color': 'white'}}, viewer=(0, 0))
view.addModel(predicted_pdb, "pdb", viewer=(0, 1))
view.setStyle({'cartoon': {'color': 'spectrum'}}, viewer=(0, 1))
view.zoomTo(); view.show()
```

> Returns a valid PDB within **2–5 seconds**! In notebook Part 7, test **1UBQ** (Ubiquitin) and **1CRN** (Crambin) — which fold accurately — vs. **1EMA** (GFP), whose barrel struggles because single-sequence ESMFold does not model post-translational chromophore cyclization.

---

# Running 3-Billion Parameter ESMFold on Colab GPUs

Can we execute the full 3B-parameter ESMFold model natively on a free Google Colab T4 GPU?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The System RAM Bottleneck
- Loading a 3-Billion parameter float32 model requires **$\sim 12\text{ GB}$ of RAM** just to read weights from disk.
- Free Colab instances provide only $12.7\text{ GB}$ of system RAM $\to$ standard loading triggers an immediate **OOM Crash**!

[[0,1]]
### The `accelerate` Streaming Solution
Hugging Face `accelerate` bypasses CPU memory by streaming model weights directly from disk into GPU VRAM:
```python
from transformers import EsmForProteinFolding

model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1",
    low_cpu_mem_usage=True,
    device_map={"": "cuda:0"},
    offload_state_dict=True
)
model.eval()
model.trunk.set_chunk_size(64)   # prevent VRAM OOM on long chains
with torch.no_grad():
    pdb_str = model.infer_pdb(my_sequence)
# Free T4 VRAM afterward: del model; gc.collect(); torch.cuda.empty_cache()
```
:::

---

# Homework Task 1 — Self-Attention from Scratch

In Homework Task 1, you write a pure NumPy function implementing the scaled dot-product attention:

$$\text{Attention}(Q, K, V) = \text{softmax}\left(\frac{QK^T}{\sqrt{d_k}}\right)V$$

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Implementation Recipe
1. **Scaled Dot Product**:
   Compute $d_k = Q.\text{shape}[-1]$ and raw compatibility scores:
   $$\text{Scores} = \frac{Q K^T}{\sqrt{d_k}}$$
   ($Q, K, V \in \mathbb{R}^{N \times d_k} \implies \text{Scores} \in \mathbb{R}^{N \times N}$).
2. **Softmax Across Rows**:
   Normalize scores into probability distributions summing to $1.0$ across columns (`axis=-1`, `keepdims=True`).
3. **Value Aggregation & Return**:
   Multiply attention weights by $V$ ($\text{Output} = \text{Weights} \cdot V \in \mathbb{R}^{N \times d_k}$) and return the tuple `(output, attention_weights)`.

[[0,1]]
### The Numerical Stability Trick ⚠
Computing $e^x$ directly can cause floating-point overflow (`inf`) if scores reach $+1000$!
To ensure numerical stability:
$$\text{Shifted} = \text{Scores} - \max(\text{Scores}, \text{axis}=-1, \text{keepdims}=\text{True})$$
$$\text{Weights} = \frac{\exp(\text{Shifted})}{\sum \exp(\text{Shifted})}$$

> Subtracting the maximum shifts all exponents $\le 0$ ($e^0 = 1$), eliminating overflow while leaving the output probabilities mathematically unchanged!
:::

---

# Homework Task 2 — Benchmarking AI Against Evolution

In Homework Task 2, you evaluate ESM-2 zero-shot predictions across experimental variants:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Input Dataset (`df_dms`)
A pandas dataframe containing experimental mutations and observed fitness values:

| mutation | fitness |
| :--- | :--- |
| `V12A` | 1.12 |
| `K3E` | 0.84 |
| `L15P` | -1.45 |
| `G20D` | -2.10 |

[[0,1]]
### Parsing Mutation Strings
- Format: `[WildType][Position][Mutant]` (e.g. `V12A`).
  - Wild-type amino acid: `mut_str[0]` (`'V'`)
  - Target residue position: `int(mut_str[1:-1])` (`12`, 1-based!)
  - Mutant amino acid: `mut_str[-1]` (`'A'`)
- Loop through the rows of `df_dms`.
- Call `predict_mutation_effect()` for each variant and store scores in `df_dms['esm_score']`.
- *Caution*: Guard bounds `1 <= pos <= len(gfp_snippet)` (snippet length is 35) before calling `predict_mutation_effect()`; assign `np.nan` for positions outside the snippet (or multi-mutants) to avoid `IndexError`.
:::

---

# Homework Task 3 — Correlation Analysis

In Homework Task 3, you quantify how accurately the AI predictions reflect experimental biology:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Computing Pearson Correlation
- Filter out `NaN` rows using `df.dropna(subset=['fitness', 'esm_score'])` before calling `pearsonr` (SciPy raises `ValueError` on `NaN`s).
- Use SciPy's stats module to compute the Pearson linear correlation coefficient $r$ and two-tailed p-value:
  ```python
  from scipy.stats import pearsonr
  corr, pval = pearsonr(df['fitness'], df['esm_score'])
  ```
- $r > 0$: Positive correlation indicating the language model successfully identifies functional vs. broken mutants.

[[0,1]]
### Plotting the Regression
- Draw a scatter plot:
  - X-axis: Experimental Fitness (Fluorescence)
  - Y-axis: Predicted ESM-2 Log-Odds Score
- Fit a linear trendline using NumPy:
  ```python
  z = np.polyfit(df['fitness'], df['esm_score'], deg=1)
  p = np.poly1d(z)
  ```
- Overlay the trendline and display the computed $r$ in the plot title.
:::

---

# Homework Task 4 — Isolating Maximum Attention Heads

In Homework Task 4, you dissect the multi-head attention tensor in Layer 6 to find a specialized contact detector:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Query
- Among the **20 attention heads** in Layer 6 (indices $0 \to 19$), which single head focuses the strongest attention weight on the connection between **residue 5 and residue 15**?

[[0,1]]
### The Index Accounting Trap ⚠
- Layer 6 attention tensor `outputs.attentions[5][0]` has shape `(20, L + 2, L + 2)` because `<cls>` sits at index `0` and `<eos>` at index `L + 1`!
  ```text
  Raw token index 0:      <cls>
  Raw token index 1:      Seq[0]
  ...
  Raw token index i + 1:  Seq[i]   (0-based residue i)
  ```
- Either index the raw tensor at `(i + 1, j + 1)` **or** slice off special tokens first (`[:, 1:-1, 1:-1]`, shape `(20, L, L)`) so 0-based sequence indices map directly.
- Convert the 20-head slice to NumPy (`.numpy()` or `.item()`) and use `np.argmax()` to find the winning Head ID and weight.
:::

---

# Homework Task 5 — Conceptual Review & Guidance

Task 5 tests your deep understanding of machine learning and structural biology principles:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### 1. Physical Meaning of $Q K^T$
- In natural language, $Q K^T$ matches verbs to direct objects.
- In proteins, what structural and physical dependencies do Query and Key vectors match?

### 2. Why MLM Predicts Thermodynamic Fitness
- The model was trained purely on evolutionary sequence survival, never on folding free energies ($\Delta\Delta G$) or fluorescence assays. Why does evolutionary likelihood correlate so strongly with physical stability?

[[0,1]]
### 3. Softmax Numerical Stability
- Why is subtracting the maximum value ($\exp(S - \max(S))$) essential for computer hardware when evaluating scaled dot products?

### 4. Epistasis: PSSMs vs. Transformers
- A Position-Specific Scoring Matrix (PSSM) evaluates columns independently.
- Why does independent column scoring fail to model epistasis? How does the transformer's full-sequence self-attention overcome this limitation?
:::

---

# Putting It All Together: Complete pLM Pipeline

```python
import torch, py3Dmol
from transformers import AutoTokenizer, EsmForMaskedLM

# 1. Load pretrained protein language model and tokenizer
model_name = "facebook/esm2_t6_8M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = EsmForMaskedLM.from_pretrained(model_name, output_attentions=True)
model.eval()

# 2. Sequence embedding & self-attention extraction
seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEG"
inputs = tokenizer(seq, return_tensors="pt")
with torch.no_grad():
    outputs = model(**inputs)
    logits = outputs.logits[0]
    attn = outputs.attentions[-1][0].mean(dim=0).numpy()[1:-1, 1:-1]

# 3. Masked-marginal zero-shot mutational scoring: V12A
wt_pos = 12
masked = seq[:wt_pos - 1] + "<mask>" + seq[wt_pos:]
with torch.no_grad():
    site_logits = model(**tokenizer(masked, return_tensors="pt")).logits[0, wt_pos]
mut_token = tokenizer.convert_tokens_to_ids("A")
wt_token  = tokenizer.convert_tokens_to_ids(seq[wt_pos - 1])
delta_score = float(site_logits[mut_token] - site_logits[wt_token])

print(f"Variant V{wt_pos}A zero-shot score: {delta_score:.2f} (negative = deleterious)")
print(f"Extracted contact map matrix shape: {attn.shape}")
```

---

# Summary

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### Language Models & Embeddings
#### Architecture
- Meta ESM-2 (8M to 15B)
- Tokenizer: `<cls>`, `<eos>`
- Contextual representations
- Multi-head self-attention
- $Q, K, V$ projections
- Denoises sequence noise

[[0,1]]{ .green }
### Attention & Structure
#### Emergence
- Co-evolution & epistasis
- Competing steric pressures
- Attention maps = contact maps
- Secondary structure ridges
- Emerges without 3D labels
- Contact detection heads

[[0,2]]{ .blue }
### Variant Fitness & ESMFold
#### Prediction
- Zero-shot log-odds scoring
- Validated by DMS assays
- Correlates with physical $\Delta\Delta G$
- ESMFold: MSA-free folding
- 3B parameters in seconds
- Fast API and GPU acceleration
:::

---

# Course Retrospective: The Computational Arc

Across these eight lectures, we traversed the entire computational scale of modern biology:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### From Genomes to Cells (Lectures 1–6)
- **Lecture 1**: $k$-mer hashing and MinHash sketching for whole-genome Mash distances.
- **Lecture 2**: Dynamic programming and affine gap penalties for comparative homology.
- **Lecture 3**: Bayesian read classification for clinical variant calling (BRCA1).
- **Lecture 4**: Negative Binomial GLMs and dispersion shrinkage for bulk RNA-seq.
- **Lecture 5**: Droplet microfluidics, sparsity filtering, and robust MAD statistics.
- **Lecture 6**: Manifold learning, spectral PCA, $k$-NN graphs, and Leiden clustering.

[[0,1]]
### To Molecules and Artificial Intelligence (Lectures 7–8)
- **Lecture 7**: 3D Cartesian coordinates, Ramachandran steric limits, the Kabsch SVD algorithm, and structural superposition.
- **Lecture 8**: Co-evolutionary memory, transformer self-attention, zero-shot mutational fitness, and MSA-free atomic folding.
:::

> **The Grand Synthesis**: From discrete characters to continuous manifolds, and from thermodynamic physics to evolutionary intelligence — computational biology unifies living systems through mathematics and algorithms.
