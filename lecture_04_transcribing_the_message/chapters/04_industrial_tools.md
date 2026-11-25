:::titlepage
[[top]]
Part 4

[[title]]
Industrial Tools

[[subcaption]]
kallisto · PyDESeq2 · pandas · scipy.stats
:::

---

# Curriculum

:::matrix{cols=50/50}
[[0,0]]
- **kallisto** — building a transcriptome index, `kallisto quant`, pseudoalignment and equivalence classes.
- **abundance.tsv** — `length`, `eff_length`, `est_counts`, `tpm`; transcripts → genes.
- **The count matrix** — genes × samples, library sizes, zeros, overdispersion in real data.
[[0,1]]
- **PyDESeq2** — `DeseqDataSet`, `deseq2()`: size factors, dispersions, GLM fit, Cook's distance.
- **DeseqStats** — Wald test, contrasts, Benjamini–Hochberg, independent filtering, LFC shrinkage.
- **scipy.stats.nbinom** — the Negative Binomial in Python, $(\mu, \alpha) \to (n, p)$.
- **Plots** — mean–variance, volcano, MA.
:::
# Roadmap

One pipeline in the practice notebook (**top**), extended in the homework (**bottom**):

:::static-diagram { center }
{
  "nodes": [
    { "id": "fq",    "label": "FASTQ\n(per sample)",          "x": 50,   "y": 80,  "width": 130, "height": 60 },
    { "id": "kal",   "label": "kallisto quant\n(+ index)",     "x": 220,  "y": 80,  "width": 160, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "geo",   "label": "GEO GSE135568\ncount matrix",   "x": 420,  "y": 80,  "width": 160, "height": 60 },
    { "id": "mv",    "label": "Mean vs variance\n(WT replicates)", "x": 620, "y": 80, "width": 170, "height": 60 },
    { "id": "dds",   "label": "PyDESeq2\ndeseq2()",            "x": 830,  "y": 80,  "width": 140, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "res",   "label": "DeseqStats\nresults_df",        "x": 1010, "y": 80,  "width": 140, "height": 60, "style": "fill: #4CAF50; color: #fff;" },

    { "id": "volc",  "label": "Volcano plot\n+ HSP check",     "x": 1010, "y": 200, "width": 150, "height": 60 },
    { "id": "top",   "label": "Top 10 up\n(padj < 0.01)",      "x": 820,  "y": 200, "width": 150, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "ma",    "label": "MA plot",                       "x": 630,  "y": 200, "width": 140, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "filt",  "label": "filter_low_counts\nsum ≥ 10",   "x": 400,  "y": 200, "width": 180, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "nb",    "label": "nb_probability\n(scipy nbinom)", "x": 150, "y": 200, "width": 190, "height": 60, "style": "fill: #FF9800; color: #fff;" }
  ],
  "links": [
    { "from": "fq",   "to": "kal" },
    { "from": "kal",  "to": "geo" },
    { "from": "geo",  "to": "mv" },
    { "from": "mv",   "to": "dds" },
    { "from": "dds",  "to": "res" },
    { "from": "res",  "to": "volc" },
    { "from": "res",  "to": "top" },
    { "from": "res",  "to": "ma" },
    { "from": "geo",  "to": "filt" },
    { "from": "mv",   "to": "nb" }
  ]
}
:::

> The practice notebook starts from a ready-made count matrix (GEO). kallisto is the step that **produces** such a matrix from raw reads — we cover it first so the matrix is not a black box.

---

# The Toolbox

Expression analysis is a **chain of tables**. Each tool reads one and writes the next:

:::matrix { .cards cols="1fr 1fr 1fr 1fr" gap="20px" margin-top: 20px; }
[[0,0]]{ .purple }
### kallisto
Pseudoaligner: FASTQ + transcriptome → estimated reads per transcript. Minutes on a laptop.

[[0,1]]{ .green }
### PyDESeq2
Python re-implementation of the R package DESeq2: Negative Binomial GLM per gene.

[[0,2]]{ .blue }
### pandas
Count matrix and results table as `DataFrame`s; filtering, sorting, joins.

[[0,3]]{ .orange }
### scipy.stats
`nbinom`, `poisson` — probabilities to check the model by hand.
:::

```python
!apt-get -y -qq install kallisto                    # command-line tool (Ubuntu package)
!pip install -q pydeseq2 pandas numpy scipy matplotlib openpyxl
#                                                     openpyxl → pd.read_excel for the GEO .xlsx
```

> Numbers in this part come from a real run: **PyDESeq2 0.5.4** on the notebook dataset. Other versions may differ in the last digits.

---

# The Data — GSE135568

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The experiment
- *Saccharomyces cerevisiae*, wild type
- **Control**: 30 °C · **Heat shock**: 39 °C for 20 min
- Bulk RNA-seq, gene-level integer counts deposited in NCBI **GEO**
- The supplementary `.xlsx` has **24** columns: WT and *mip6Δ* mutant, 30 °C / 39 °C, 0 / 20 / 120 min, 4 replicates each
- We take 3 + 3 columns

[[0,1]]
### Column → sample name

| GEO column | Our name | Condition |
| :--- | :--- | :--- |
| `wt.30.0.1` | `WT_1` | WT |
| `wt.30.0.2` | `WT_2` | WT |
| `wt.30.0.3` | `WT_3` | WT |
| `wt.39.20.1` | `HS_1` | HeatShock |
| `wt.39.20.2` | `HS_2` | HeatShock |
| `wt.39.20.3` | `HS_3` | HeatShock |

Names encode *genotype.temperature.minutes.replicate*.
:::

```python
df = pd.read_excel("GSE135568_CountMatrixRNAseq.xlsx", index_col=0)
df.columns = df.columns.str.strip()          # GEO headers carry stray spaces → KeyError without this
```

---

# The Table Chain

```
 sequencer         kallisto quant             sum per gene             PyDESeq2
   FASTQ   ──────────────────────►  abundance.tsv  ───────────►  count matrix  ───────────►  results_df
 (Lecture 1)   "which transcripts     (per sample,     (genes ×      "how many        (genes ×     "which genes
                could this read        per transcript)   samples)      reads per        statistics)    changed?"
                come from"                                              gene"
```

:::matrix { cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]
### abundance.tsv
- One file **per sample**
- One row **per transcript**
- `est_counts` (fractional!) and `tpm`

[[0,1]]
### Count matrix
- One row per **gene**, one column per **sample**
- **Raw integer** counts
- 7,126 × 6 in our notebook

[[0,2]]
### results_df
- One row per gene
- `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`
- What biologists actually read
:::

> Lecture 3 needed **where** each read sits (BAM, base by base). Expression needs only **how many** reads each gene got — this is why a much faster method than alignment is enough.

---
# Tool 1 — kallisto

## Pseudoalignment

:::matrix { cols="61/39" gap="30px" }
[[0,0]]
### Two commands

```bash
# 1. Index the TRANSCRIPTOME (cDNA), not the genome — once per organism
kallisto index -i yeast_R64.idx \
        Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# 2. Quantify each sample (paired-end)
kallisto quant -i yeast_R64.idx -o quant/WT_1 \
        -t 4 -b 100 \
        WT_1_R1.fastq.gz WT_1_R2.fastq.gz
```

[[0,1]]
### What it skips
A normal aligner (BWA, STAR, HISAT2) computes the **exact position and alignment** of every base — CIGAR, mismatches, MAPQ.

kallisto only answers: **which transcripts is this read compatible with?**

- No base-by-base alignment
- No BAM file
- $\sim 10^7$ reads in a few minutes on a laptop

:::

:::matrix { cols="55/45" gap="30px" }
[[0,0]]
| Flag | Meaning |
| :--- | :--- |
| `-i` | index file |
| `-o` | output **directory** (one per sample) |
| `-t` | threads |
| `-b` | bootstrap rounds (uncertainty of estimates) |
[[0,1]]
### Single-end reads
Fragment length cannot be measured, so you must give it:
```bash
kallisto quant -i yeast_R64.idx -o quant/WT_1 \
        --single -l 200 -s 20 WT_1.fastq.gz
#                 └ mean   └ sd of fragment length
```
:::

---

# How Pseudoalignment Works

The index is a **de Bruijn graph** of all transcript k-mers ($k = 31$ by default; k-mers as in mash, Lecture 1). Each k-mer stores the set of transcripts that contain it.

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Toy example: one gene, three isoforms

| k-mer of the read | Transcripts containing it |
| :--- | :--- |
| k₁ (shared exon) | {T1, T2, T3} |
| k₂ (shared exon) | {T1, T2, T3} |
| k₃ (exon 2) | {T1, T2} |
| k₄ (exon 2–3 junction) | {T1, T2} |

Read is compatible with
$$\{T1,T2,T3\} \cap \{T1,T2\} = \{T1, T2\}$$

This set is the read's **equivalence class**.

[[0,1]]
### From classes to counts — EM
Millions of reads collapse into a few thousand equivalence classes with counts. Then the **Expectation–Maximization** algorithm:

1. Start with equal abundances
2. **E-step**: split each class's reads between its transcripts in proportion to current abundance / length
3. **M-step**: recompute abundances from the split reads
4. Repeat until stable

Ambiguous reads go mostly to the isoform that unique reads show to be abundant. That is why `est_counts` are **not integers**.
:::

---

# kallisto Output

```
quant/WT_1/
├── abundance.tsv      ← human-readable table
├── abundance.h5       ← same + bootstraps (binary, for sleuth)
└── run_info.json      ← n_processed, n_pseudoaligned, p_pseudoaligned
```

```
target_id     length  eff_length  est_counts  tpm
YAL005C_mRNA  1929    …           …           …      ← SSA1 (642 aa + stop)
YLL026W_mRNA  2727    …           …           …      ← HSP104 (908 aa + stop)
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Columns
| Column | Meaning |
| :--- | :--- |
| `length` | transcript length, bp |
| `eff_length` | positions where a fragment can **start**: ≈ length − mean fragment + 1 |
| `est_counts` | reads assigned by EM (fractional) |
| `tpm` | transcripts per million |

[[0,1]]
### TPM

$$ \text{TPM}_i = \frac{c_i / \tilde{\ell}_i}{\sum_j c_j / \tilde{\ell}_j} \cdot 10^6 $$

$c_i$ = est_counts, $\tilde{\ell}_i$ = eff_length. Long transcripts collect more fragments at the same molecule number → divide by length.

- `p_pseudoaligned` in `run_info.json` is the first QC number: a clean poly-A library typically gives 70–90%; < 50% → wrong organism, rRNA, adapters
:::

---

# ⚠ Pitfall — TPM Is Not Input for DESeq2

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### TPM / FPKM
- Normalized for length and depth
- Good for **comparing genes within a sample**
- Lose the information *how many reads* stand behind a number

TPM = 100 can be 3 reads on a short gene or 3,000 reads on a long one — **very different certainty**.

[[0,1]]
### Raw counts
- The Negative Binomial models **counts** — the noise depends on how many reads there are (Part 3)
- DESeq2 does its own depth normalization (size factors)
- Rule: give PyDESeq2 **estimated counts**, rounded to integers
:::

### From transcripts to genes

```python
import pandas as pd
ab = pd.read_csv("quant/WT_1/abundance.tsv", sep="\t", index_col=0)
tx2gene = ab.index.str.replace("_mRNA", "")                 # yeast: 1 transcript ≈ 1 gene
gene_counts = ab["est_counts"].groupby(tx2gene).sum().round().astype(int)
```

> Yeast has few introns and almost no alternative isoforms, so the transcript → gene map is nearly 1 : 1. In human (~4 isoforms per gene) use a proper `tx2gene` table from the annotation (GTF) — or `tximport` in R.

---
# The Count Matrix


:::matrix { cols="1.3fr 1fr" gap="30px" }
[[0,0]]{ .smaller-table }
```python
counts_df = pd.read_csv("yeast_heat_shock_counts.csv", index_col=0)
counts_df.shape            # (7126, 6)
```
### First rows (real)

| Gene | WT_1 | WT_2 | WT_3 | HS_1 | HS_2 | HS_3 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| YIL100W | 7 | 0 | 2 | 14 | 4 | 8 |
| YDR130C | 293 | 365 | 253 | 343 | 317 | 316 |
| YDR473C | 293 | 162 | 179 | 270 | 208 | 238 |

Index = systematic ORF names: **Y** + chromosome (**A**–**P**) + arm (**L**/**R**) + number + strand (**W**/**C**).
`YLL026W` = chr XII, left arm, ORF 26, Watson strand = *HSP104*.

[[0,1]]
### Library sizes (column sums)

| Sample | Reads counted |
| :--- | ---: |
| WT_1 | 21.5 M |
| WT_2 | 22.9 M |
| WT_3 | 19.1 M |
| HS_1 | 21.1 M |
| HS_2 | 22.8 M |
| HS_3 | 23.1 M |

- **593** genes have zero reads in all 6 samples — and hundreds more have very low totals across all samples
- In `pandas`: `.sum(axis=0)` sums **columns** (sample library sizes); `.sum(axis=1)` sums **rows** (total reads per gene — needed for low-count filtering in Task 1)
:::

---
# Overdispersion in Real Data

Notebook Part 3: mean and variance of each gene across the 3 **WT** replicates.

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Median variance / mean ratio

| Mean count | Genes | Median $\text{Var}/\mu$ |
| :--- | ---: | ---: |
| 1 – 10 | 462 | **1.0** |
| 10 – 100 | 565 | 1.9 |
| 100 – 1,000 | 2,415 | 6.3 |
| > 1,000 | 2,764 | **32** |

88% of expressed genes have $\text{Var} > \mu$.

[[0,1]]
### Reading the table
- Poisson predicts $\text{Var}/\mu = 1$ everywhere
- **Low counts** — Poisson is fine: counting noise dominates
- **High counts** — variance grows like $\mu^2$: biological variation between cultures dominates

$$ \text{Var} = \mu + \alpha\mu^2 \\;\Rightarrow\\; \frac{\text{Var}}{\mu} = 1 + \alpha\mu $$

With $\alpha \approx 0.015$ (next slides): $\mu = 2000 \Rightarrow \text{Var}/\mu \approx 31$ ✓
:::

> On the log–log plot the points lie on the red $y = x$ line at the left and bend **above** it at the right — exactly this formula.

---

# Tool 2 — PyDESeq2

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds  import DeseqStats

dds = DeseqDataSet(
    counts=counts_df.T,        # samples × genes  ← TRANSPOSED
    metadata=metadata_df,      # samples × factors
    design="~ condition",      # formula (R-style)
)
dds.deseq2()                   # the whole fitting pipeline
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Input requirements
| Object | Shape | Must be |
| :--- | :--- | :--- |
| `counts` | 6 × 7126 | **non-negative integers** |
| `metadata` | 6 × 1 | index = same sample names |

[[0,1]]
### ⚠ Orientation
Count matrices are stored **genes × samples** (GEO, featureCounts, kallisto merge). PyDESeq2 follows the AnnData convention: **observations (samples) × variables (genes)**.

Forget `.T` → PyDESeq2 sees 7,126 "samples" and 6 "genes" and fails (the metadata index no longer matches).
:::

> `DeseqDataSet` is an **AnnData** object (the same container as Scanpy in Lectures 5–6). Results are stored inside it, not returned.

---

# What `deseq2()` Does

The log of the real run, annotated:

```
Fitting size factors...                 ← 1. sequencing-depth normalization
Fitting dispersions...                  ← 2. α per gene (maximum likelihood)
Fitting dispersion trend curve...       ← 3. α as a function of mean expression
Fitting MAP dispersions...              ← 4. empirical Bayes shrinkage toward the trend
Fitting LFCs...                         ← 5. NB GLM: log2 fold change per gene
Calculating cook's distance...          ← 6. outlier detection
Replacing 0 outlier genes.              ←    (replacement needs ≥ 7 replicates; we have 3)
```

Then `DeseqStats(...).summary()`:

```
Running Wald tests...                   ← 7. test LFC = 0, p-values, BH adjustment
```

> Steps 2–4 are the empirical Bayes procedure from Part 3. Each step stores its result inside `dds` — next slides show where.

---
# Step 1 — Size Factors (Median of Ratios)

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The algorithm
1. For each gene: geometric mean across all samples (a "pseudo-reference" sample)
2. For each sample: ratio count / pseudo-reference, for every gene with no zeros
3. Size factor $s_j$ = **median** of these ratios

$$ s_j = \operatorname*{median}\_i \frac{K_{ij}}{\left(\prod_{v} K_{iv}\right)^{1/m}} $$

The median ignores the few genes that really change — **assumes most genes do not**.

[[0,1]]
### Our run: `dds.obs["size_factors"]`

| Sample | $s_j$ | Library (M) |
| :--- | ---: | ---: |
| WT_1 | 0.901 | 21.5 |
| WT_2 | 0.937 | 22.9 |
| WT_3 | 0.821 | 19.1 |
| HS_1 | 1.142 | 21.1 |
| HS_2 | 1.103 | 22.8 |
| HS_3 | 1.167 | 23.1 |

HS_1 has a *smaller* library than WT_2 but a *larger* size factor.
:::

> Why? After heat shock a few chaperone mRNAs take a large share of all reads (SSA1 alone: ~300,000 per sample). Every other gene gets fewer reads at the same depth. Dividing by the total read count would make thousands of genes look **down-regulated**; the median of ratios does not.

---
# Normalized Counts

$$ \tilde{K}\_{ij} = K_{ij} / s_j \qquad \texttt{dds.layers["normed\_counts"]} $$

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Gene | | WT_1 | WT_2 | WT_3 | HS_1 | HS_2 | HS_3 |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| YDR130C | raw | 293 | 365 | 253 | 343 | 317 | 316 |
| | normalized | 325 | 390 | 308 | 300 | 288 | 271 |
| YLL026W (*HSP104*) | raw | 8,577 | 5,157 | 7,479 | 102,183 | 87,168 | 115,739 |
| | normalized | 9,514 | 5,503 | 9,110 | 89,500 | 79,061 | 99,164 |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Use normalized counts for
- Plots of single genes
- Clustering, heatmaps (after log)
- `baseMean` = mean of the normalized counts

[[0,1]]
### Never feed them back into DESeq2
The model takes **raw** counts and puts $s_j$ into the mean: $\mu_{ij} = s_j q_{ij}$. Dividing first destroys the information on how many reads there were.
:::

---

# Steps 2–4 — Dispersions

Stored in `dds.var` (one row per gene):

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Column | Step | Median over genes |
| :--- | :--- | ---: |
| `genewise_dispersions` | 2 — ML estimate from this gene's 6 values only | 0.0074 |
| `fitted_dispersions` | 3 — trend value at this gene's mean | 0.0172 |
| `MAP_dispersions` | 4 — posterior mode (empirical Bayes) | 0.0150 |
| `dispersions` | final value used in the test | 0.0150 |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The fitted trend (our data)

$$ \alpha_{\text{trend}}(\bar{\mu}) = a_0 + \frac{a_1}{\bar{\mu}} = 0.0137 + \frac{3.03}{\bar{\mu}} $$

`dds.uns["trend_coeffs"]` → `a0 = 0.0137, a1 = 3.03`

- Highly expressed genes: $\alpha \to 0.0137$ → biological CV ≈ $\sqrt{0.0137}$ ≈ **12%** between cultures
- Weak genes: larger $\alpha$

[[0,1]]
### Four genes

| Gene | Gene-wise | Trend | Final |
| :--- | ---: | ---: | ---: |
| YDR130C | 0.0058 | 0.0234 | 0.0152 |
| SSA1 | 0.0276 | 0.0137 | 0.0201 |
| HSP104 | 0.0490 | 0.0138 | 0.0286 |
| YIL100W (low) | 0.687 | 0.566 | 0.602 |

Every estimate is pulled **toward the trend** — up if too small, down if too large.
:::

---
# Why Shrink the Dispersion?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Three replicates are not enough
A variance from 3 numbers is very noisy. Gene-wise estimates in our run range from $10^{-8}$ (clipped lower bound) to $10$ (clipped upper bound).

A gene whose 3 replicates *happen* to agree gets $\alpha \approx 0$ → tiny standard error → huge Wald statistic → **false positive**.

YDR130C: gene-wise 0.0058 → final 0.0152 (×2.6). Its test is now less overconfident.

[[0,1]]
### How much to shrink?
DESeq2 estimates how far genes really scatter around the trend:

`dds.uns["prior_disp_var"]` = 0.51 (variance of $\log\alpha$ around the trend)

- Small scatter → strong pull to the trend
- More replicates → the gene's own data weigh more, shrinkage weakens
- **Low raw counts** → Poisson shot noise ($\mu$) dominates $\mu + \alpha\mu^2$, leaving a flat likelihood where the prior (trend) dominates; **high counts** yield a sharp likelihood that outweighs the prior

Genes whose gene-wise $\alpha$ is far **above** the trend keep their high value (dispersion-outlier rule, `dds.var["_outlier_genes"]`: **94** genes in our run) — never shrink a truly variable gene into significance.
:::

---

# Where Results Live in `dds`

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Slot | Content | Example |
| :--- | :--- | :--- |
| `dds.X` | raw counts (samples × genes) | `dds.X[0, :5]` |
| `dds.obs` | per-sample: `condition`, `size_factors` | `dds.obs["size_factors"]` |
| `dds.var` | per-gene: dispersions, `_normed_means`, outlier flags | `dds.var["dispersions"]` |
| `dds.varm["LFC"]` | GLM coefficients, **natural log** scale | columns `Intercept`, `condition[T.…]` |
| `dds.layers["normed_counts"]` | counts / size factor | for plots |
| `dds.layers["cooks"]` | Cook's distance per sample and gene | outlier check |
| `dds.uns` | dataset-wide: `trend_coeffs`, `prior_disp_var` | |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### ⚠ ln vs log₂
`dds.varm["LFC"]` is in **natural log**. `results_df["log2FoldChange"]` is converted to log₂:

$$ \log_2 FC = \beta_1 / \ln 2 $$

[[0,1]]
### The GLM behind it
$$ \log \mu_{ij} = \log s_j + \beta_{0,i} + \beta_{1,i}\, x_j $$

$x_j = 1$ for the non-reference condition, 0 otherwise. $\beta_1$ is the log fold change.
:::

---

# Step 7 — DeseqStats and the Wald Test

```python
stat_res = DeseqStats(dds, contrast=["condition", "HeatShock", "WT"])
#                                     factor      numerator    denominator
stat_res.summary()
res_df = stat_res.results_df
```

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Column | Meaning | HSP104 (YLL026W) |
| :--- | :--- | ---: |
| `baseMean` | mean of normalized counts, all 6 samples | 48,642 |
| `log2FoldChange` | $\log_2(\text{HeatShock}/\text{WT})$ from the GLM | 3.47 (×11) |
| `lfcSE` | standard error of the LFC | 0.20 |
| `stat` | Wald statistic $z = \text{LFC}/\text{SE}$ | 17.4 |
| `pvalue` | $2\,P(Z > \lvert z \rvert)$, standard normal | ~10⁻⁶⁷ → prints 0.0 |
| `padj` | Benjamini–Hochberg adjusted | 0.0 |
:::

> `pvalue = 0.0` is **underflow**, not "impossible". For a volcano plot clip before the log: `-np.log10(padj.clip(lower=1e-300))` (done in notebook Part 6).

---

# ⚠ Pitfall — Which Condition Is the Reference?

The formula needs a **reference level**. By default it is the **alphabetically first** level:

`"HeatShock" < "WT"` → reference = HeatShock → coefficient `condition[T.WT]` = **WT vs HeatShock**

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Safe: explicit contrast
```python
DeseqStats(dds, contrast=["condition", "HeatShock", "WT"])
```
HSP104: LFC = **+3.47** ✓

### Unsafe: coefficient name
```python
dds.varm["LFC"].columns
# ['Intercept', 'condition[T.WT]']
stat_res.lfc_shrink(coeff="condition[T.WT]")
```
HSP104: LFC = **−3.46** ✗ — sign flipped, every gene

[[0,1]]
### Fix: fix the level order
```python
metadata_df["condition"] = pd.Categorical(
    metadata_df["condition"],
    categories=["WT", "HeatShock"])     # first = reference
# → coefficient 'condition[T.HeatShock]'
```

Always print `dds.varm["LFC"].columns` and check **one known gene** (here: HSP104 must go **up** in heat shock).
:::

---

# Why Are Some `padj` NaN?

```python
res_df["padj"].isna().sum()      # 721 of 7126
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
| Reason | Genes | `pvalue` | `padj` |
| :--- | ---: | :--- | :--- |
| Zero counts in all samples | 593 | NaN | NaN |
| Cook's distance outlier | 2 | NaN | NaN |
| **Independent filtering** | 126 | value | NaN |
| **Total** | **721** | | |

[[0,1]]
### Independent filtering
Genes with very low `baseMean` can never reach significance — too few reads. Testing them only **adds to the multiple-testing burden**.

DESeq2 picks a `baseMean` threshold that **maximizes the number of discoveries**, and sets `padj = NaN` below it. The filter uses only `baseMean`, which does not depend on the condition labels → p-values stay valid.
:::

> NaN ≠ "not significant": `res_df[res_df.padj < 0.05]` silently drops NaN rows (correct); `res_df.padj.fillna(0)` would make them **all significant** (wrong). In the notebook: `dropna(subset=["padj", "log2FoldChange"])`.

---
# Multiple Testing — p vs padj

### 7,126 tests at once
If no gene changed at all, $p < 0.05$ would still flag ~5% by chance:

$$ 0.05 \times 6{,}400 \text{ tested} \approx 320 \text{ false positives} $$

| Threshold | Genes |
| :--- | ---: |
| `pvalue < 0.05` | 2,903 |
| `padj < 0.05` | **2,425** |

### Benjamini–Hochberg
Sort the $m$ p-values: $p_{(1)} \le \dots \le p_{(m)}$

$$ p^{\text{adj}}\_{(k)} = \min_{l \ge k} \frac{m \\, p_{(l)}}{l} $$

`padj < 0.05` → among all genes on the list, the expected share of false discoveries is ≤ **5%** (FDR) — at most ~120 of 2,425.

Not "each gene has 95% chance to be real".

---

# Results — Heat Shock vs WT

`padj < 0.05`: **2,425** genes (34% of all). With $\lvert\log_2 FC\rvert > 1$ (≥ 2-fold): **355 up**, **283 down**.

Top genes by `padj` (all printed `padj = 0.0`):

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Gene | Name | WT counts | HS counts | baseMean | log₂FC | Wald $z$ |
| :--- | :--- | :--- | :--- | ---: | ---: | ---: |
| YER103W | *SSA4* (Hsp70) | 757, 645, 697 | 85,492 · 67,834 · 86,284 | 35,452 | 6.47 | 47.1 |
| YFL014W | *HSP12* | 55, 61, 57 | 45,143 · 31,166 · 47,656 | 18,139 | **9.12** | 46.5 |
| YPL240C | *HSP82* (Hsp90) | 14,138 · 12,643 · 11,686 | 169,383 · 178,740 · 192,521 | 86,473 | 3.45 | 31.7 |
| YBL075C | *SSA3* (Hsp70) | | | 6,169 | 3.85 | 29.2 |
| YBR072W | *HSP26* | 338, 506, 289 | 92,862 · 62,873 · 96,649 | 37,073 | 7.45 | 28.8 |
| YGL009C | *LEU1* | | | 14,237 | **−2.80** | −22.4 |
:::

> HSP12: ~57 reads → ~41,000 reads, a **560-fold** induction in 20 minutes — the Hsf1 / Msn2/4 heat-shock response at work.

---

# Biology Check — the Chaperones

Notebook Part 7 checks two classic chaperones:

```python
for gene in ["YLL026W", "YAL005C"]:
    print(gene, res_df.loc[gene, "log2FoldChange"], res_df.loc[gene, "padj"])
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
| | HSP104 (YLL026W) | SSA1 (YAL005C) |
| :--- | ---: | ---: |
| WT raw | 8,577 · 5,157 · 7,479 | 35,995 · 53,803 · 31,598 |
| HS raw | 102,183 · 87,168 · 115,739 | 287,423 · 316,696 · 327,972 |
| log₂FC | **3.47** (×11) | **2.59** (×6) |
| lfcSE | 0.20 | 0.17 |
| padj | 0.0 | 0.0 |

[[0,1]]
### Why a positive control?
Before interpreting 2,425 genes, check genes whose answer **you already know**:
- wrong sign → reference level swapped
- no change → samples mislabeled
- tiny counts → wrong organism / annotation

SSA1 is already high at 30 °C (constitutive Hsp70) and still rises 6×; HSP104 is low and rises 11×.
:::

---

# ⚠ Pitfall — Ranking by Fold Change Alone

Sorting significant genes by raw `log2FoldChange` mixes well-measured biological hits with noisy low-count genes:

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Gene | baseMean | log₂FC | lfcSE | Shrunk log₂FC |
| :--- | ---: | ---: | ---: | ---: |
| YFL014W (*HSP12*, high counts) | 18,139 | 9.12 | 0.20 | 9.11 |
| YIL100W (low counts, from `head(3)`) | **5** | 1.16 | **1.07** | **0.21** |
| *Illustrative low-count gene* (2 → 20 reads) | **~11** | ~3.3 | **> 1.0** | *pulled toward 0* |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The problem
A gene going from 2 to 20 reads has the same fold change as one going from 2,000 to 20,000 — but the estimate is far noisier (SE 1.2 vs 0.2). Low-count genes dominate the extremes of any LFC ranking.

[[0,1]]
### LFC shrinkage
```python
stat_res.lfc_shrink(coeff="condition[T.HeatShock]")
```
Pulls noisy LFCs toward 0 (apeGLM-type prior); well-measured LFCs barely move. Use shrunk LFC for **ranking and plots**; p-values are not changed. Mind the coefficient name — previous pitfall!
:::

---
# scipy.stats — the Negative Binomial by Hand

scipy uses the "number of successes" form NB($n$, $p$). DESeq2 uses mean and dispersion. Conversion (homework Task 3):

$$ n = \frac{1}{\alpha}, \qquad p = \frac{n}{n + \mu} \qquad\Rightarrow\qquad \text{mean} = \mu,\\;\\; \text{Var} = \mu + \alpha\mu^2 $$

:::matrix { cols="1.2fr 1fr" gap="30px" }
[[0,0]]
```python
from scipy.stats import nbinom, poisson
mu, alpha = 500, 0.02
n = 1 / alpha
p = n / (n + mu)

nbinom.mean(n, p), nbinom.var(n, p)   # 500.0, 5500.0
nbinom.pmf(500, n, p)                 # 0.0054
poisson.pmf(500, mu)                  # 0.0178
nbinom.interval(0.95, n, p)           # (365, 655)
poisson.interval(0.95, mu)            # (457, 544)
```

Always check `nbinom.mean` and `nbinom.var` after converting — swapping $p \leftrightarrow 1-p$ gives a completely wrong mean without any error.

[[0,1]]
### What it means
A replicate with **650** reads, if the true mean is 500:

| Model | $P(K \ge 650)$ |
| :--- | ---: |
| Poisson | $8 \times 10^{-11}$ |
| NB, $\alpha = 0.02$ | **0.029** |

Poisson calls an ordinary biological fluctuation impossible. Any test built on it reports thousands of false discoveries.
:::

---

# Dispersion by Hand vs PyDESeq2

Notebook Part 8 — method of moments for one gene, 3 WT replicates (450, 620, 510):

$$ \hat{\mu} = 526.7, \quad s^2 = 7{,}433, \quad \hat{\alpha} = \frac{s^2 - \hat{\mu}}{\hat{\mu}^2} = 0.0249 $$

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### What the hand formula ignores
- **Sample variance**: use `np.var(x, ddof=1)` (unbiased $s^2$; NumPy defaults to `ddof=0`, which underestimates variance by $1/3$ at $n = 3$, whereas pandas `.var()` defaults to `ddof=1`)
- Different **size factors** per sample — must divide by $s_j$ first
- The **condition** effect — use only replicates of one group (as here) or residuals from the GLM
- Can be **negative** when $s^2 < \mu$ (by chance) → clip to a small positive value

PyDESeq2 stores exactly this rough estimate as `dds.var["_MoM_dispersions"]` and uses it only as the **starting point** for the maximum-likelihood fit.

[[0,1]]
### Compare with the real trend
At $\bar\mu \approx 527$:

$$ \alpha_{\text{trend}} = 0.0137 + 3.03/527 \approx 0.0195 $$

The hand estimate (0.0249) is a little above the trend. The MAP step would move it to between the two values.
:::

---

# Plots — Volcano and MA

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Volcano (notebook Part 6)
- x: `log2FoldChange` · y: $-\log_{10}$ `padj`
- Thresholds: $\lvert\log_2 FC\rvert = 1$, `padj = 0.05`
- **Asymmetric** here: the HSP genes form a far-right "plume" at $z \approx 30$–47
- Shows **effect size vs. evidence** together

```python
plot_df = res_df.dropna(subset=["padj", "log2FoldChange"]).copy()
plot_df["nlp"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
```

[[0,1]]
### MA plot (homework Task 2)
- x: `baseMean` on **log** scale · y: `log2FoldChange`
- Typical "funnel": wide spread at low `baseMean`, narrow at high — the SE story from the previous slides
- Cloud should be centred at **0**: if not, normalization failed
- Colour by `padj < 0.05`: at low counts only huge LFCs are significant
- **Zero-count guard**: 593 genes have `baseMean == 0` ($\log_{10} 0 = -\infty$) — drop `NaN`s and filter `baseMean > 0`, use `np.log10(baseMean + 1)`, or set `plt.xscale("log")`

Without the log x-axis 90% of the points sit in one column at the left edge.
:::

---

# Putting It All Together

```python
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds  import DeseqStats

# 0. (upstream) kallisto quant per sample → est_counts → gene sums → integer matrix
counts_df   = pd.read_csv("yeast_heat_shock_counts.csv", index_col=0)      # genes × samples
metadata_df = pd.read_csv("yeast_metadata.csv", index_col=0)
metadata_df["condition"] = pd.Categorical(metadata_df["condition"], categories=["WT", "HeatShock"])

# 1. Fit: size factors → dispersions → trend → MAP → GLM → Cook's
dds = DeseqDataSet(counts=counts_df.T, metadata=metadata_df, design="~ condition")
dds.deseq2()
print(dds.obs["size_factors"], dds.varm["LFC"].columns.tolist())

# 2. Test HeatShock vs WT and inspect top hits by padj
stat_res = DeseqStats(dds, contrast=["condition", "HeatShock", "WT"])
stat_res.summary()
res_df = stat_res.results_df
print(res_df.sort_values("padj").head(5))

# 3. Positive controls (known heat-shock chaperones)
for gene in ["YLL026W", "YAL005C"]:
    print(gene, res_df.loc[gene, ["log2FoldChange", "padj"]].to_dict())
```

> Useful `pandas` primitives for downstream filtering & ranking: row sums `df.sum(axis=1)`, boolean indexing `df[mask]`, descending sort `df.sort_values(col, ascending=False)`, and `.head(k).index.tolist()` to return gene IDs as a Python `list`.

---

# Summary

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### kallisto
#### Quantification
- Index the transcriptome
- k-mer compatibility, no alignment
- Equivalence classes + EM
- `est_counts` fractional, `tpm`
- TPM ≠ input for DESeq2
- **Output**: counts per transcript

[[0,1]]{ .green }
### PyDESeq2 `deseq2()`
#### Model
- `counts.T`: samples × genes
- Median-of-ratios size factors
- Gene-wise → trend → MAP $\alpha$
- Trend: $0.0137 + 3.03/\bar\mu$
- NB GLM, LFC in ln inside `varm`
- **Output**: fitted `dds`

[[0,2]]{ .blue }
### DeseqStats
#### Inference
- Explicit `contrast`, check the reference
- Wald $z$ = LFC / SE
- BH: `padj`, NaN ≠ not significant
- Shrink LFC before ranking
- Positive control: HSP104 ×11
- **Output**: 2,425 genes at FDR 5%
:::

---

# Next Steps

### In the Practice Session You Will:

1. **Download** the GSE135568 count matrix from GEO and select 3 WT + 3 heat-shock replicates
2. **Plot** mean vs variance and see overdispersion with your own eyes
3. **Fit** the Negative Binomial GLMs with `PyDESeq2` and extract the Heat Shock vs WT contrast
4. **Draw** a volcano plot and check the chaperones HSP104 and SSA1
5. **Compute** a dispersion by hand

### Homework — Differential Expression

1. `filter_low_counts()` — keep genes with total ≥ 10 reads
2. MA plot with a log-scaled `baseMean` axis
3. `nb_probability()` with `scipy.stats.nbinom` — convert $(\mu, \alpha) \to (n, p)$
4. Top 10 up-regulated genes (padj < 0.01, sorted by log₂FC)
5. Interpretation questions

> **Questions to keep in mind**: Why counts and not TPM? Why does the median, not the sum, normalize depth? Why do 3 replicates need help from 7,000 other genes? Which condition is your reference?
