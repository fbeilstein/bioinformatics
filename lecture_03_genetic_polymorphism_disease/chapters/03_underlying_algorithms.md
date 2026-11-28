:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
FM-Index · BGZF · Bayesian Genotyping · HWE Priors · Log-Sum-Exp · Haplotype Calling
:::

---

# Curriculum & Algorithmic Roadmap

How short Illumina reads become **statistically calibrated, normalized clinical variants** in six algorithmic stages:

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "reads", "label": "Illumina Reads\n(FASTQ)",           "x": 20,  "y": 90, "width": 145, "height": 64 },
    { "id": "fm",    "label": "BWT / FM-Index\nBWA-MEM + MAPQ",     "x": 205, "y": 90, "width": 165, "height": 64, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "idx",   "label": "BGZF / BAI / CRAM\n+ Dup Flag",      "x": 410, "y": 90, "width": 170, "height": 64, "style": "fill: #2196F3; color: #fff;" },
    { "id": "lik",   "label": "Phred Likelihood\nP(D | G) Log Space","x": 620, "y": 90, "width": 175, "height": 64, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "bayes", "label": "HWE Prior + LSE\nMAP GT · QUAL · GQ","x": 835, "y": 40, "width": 185, "height": 62, "style": "fill: #FF9800; color: #fff;" },
    { "id": "qc",    "label": "Haplotype Norm\n& Depth / AAF QC",   "x": 835, "y": 145,"width": 185, "height": 62, "style": "fill: #009688; color: #fff;" }
  ],
  "links": [
    { "from": "reads", "to": "fm" },
    { "from": "fm",    "to": "idx" },
    { "from": "idx",   "to": "lik" },
    { "from": "lik",   "to": "bayes" },
    { "from": "lik",   "to": "qc" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="28px" margin-top: 10px; }
[[0,0]]
- **Alignment, MAPQ & Indexing** — BWT / FM-index seeding (BWA-MEM), Phred $\text{MAPQ}$, BGZF + hierarchical `.bai`/`.tbi`/`.crai` indexes, CRAM delta encoding, and PCR duplicate marking.
- **Bayesian Diploid Classification** — Bayes' theorem over $G \in \{0/0, 0/1, 1/1\}$, per-read Phred error likelihoods ($\varepsilon_i = 10^{-Q_i/10}$), biallelic vs. 4-base models, and MAP calling.
- **Population Priors (HWE)** — Hardy–Weinberg priors from allele frequency $f$ and the prior-vs-depth crossover $n^*$.

[[0,1]]
- **Log-Space Numerical Stability** — preventing float64 underflow via log-likelihood sums, log-sum-exp normalization, and `Q = 0` ($\varepsilon = 1$) clamping.
- **Confidence & Haplotype Calling** — site `QUAL` vs. sample `GL` and `GQ`, haplotype evaluation for clustered variants, and Tan et al. canonical indel left-alignment.
- **Coverage & AAF Statistics** — binomial heterozygote dropout $2(1/2)^n$ (`DP ≥ 15`), AAF spread $\frac{1}{2\sqrt{n}}$, and reference/strand bias.
:::

---

# Read Alignment: BWT / FM-Index & Mapping Quality

Aligning $10^8$ reads of length $L = 76\text{ bp}$ against the $3 \times 10^9\text{ bp}$ human genome via full Smith–Waterman ($O(GL)$ per read) would require **$> 10^{19}$ operations**! Modern aligners (`BWA-MEM`) replace full DP with **seed-and-extend**:

:::static-diagram { width=94% center }
{
  "nodes": [
    { "id": "r",   "label": "76 bp Read",                     "x": 20,  "y": 40, "width": 130, "height": 58 },
    { "id": "fm",  "label": "FM-Index Backward Search\nO(|seed|) Exact SMEMs", "x": 190, "y": 40, "width": 230, "height": 58, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "ch",  "label": "Collinear Seed\nChaining",       "x": 460, "y": 40, "width": 155, "height": 58, "style": "fill: #2196F3; color: #fff;" },
    { "id": "sw",  "label": "Banded Smith–Waterman\nExtension (Local DP)", "x": 655, "y": 40, "width": 225, "height": 58, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "mq",  "label": "MAPQ from\nS₁ − S₂ Gap",         "x": 920, "y": 40, "width": 145, "height": 58, "style": "fill: #FF9800; color: #fff;" }
  ],
  "links": [
    { "from": "r",  "to": "fm" },
    { "from": "fm", "to": "ch" },
    { "from": "ch", "to": "sw" },
    { "from": "sw", "to": "mq" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="28px" margin-top: 8px; }
[[0,0]]
### Burrows–Wheeler Transform (BWT) & FM-Index
- **BWT** permutes the reference genome into the last column of all sorted cyclic rotations — clustering identical bases for high compression while remaining **invertible** via the **LF-mapping** property.
- **FM-Index** augments the BWT with two compact tables: `C[c]` (count of bases $< c$) and `Occ(c, r)` (occurrences of $c$ in prefix $0\dots r$).
- **Backward search** locates all exact matches (SMEMs) of a seed $P$ in **$O(|P|)$ time — completely independent of genome length $G$**!

[[0,1]]
### Phred-Scaled Mapping Quality (`MAPQ`)
$$\text{MAPQ} = -10\log_{10} P(\text{read placed at wrong locus})$$
- Driven by the alignment-score gap **$S_1 - S_2$** between the best hit ($S_1$) and second-best hit ($S_2$) in the genome:
  - **`MAPQ = 60`** ($P_{\text{wrong}} \le 10^{-6}$): unique placement ($S_1 \gg S_2$).
  - **`MAPQ = 0`** ($P_{\text{wrong}} \ge 0.5$): multi-mapping read inside a repeat ($S_1 = S_2$) — skipped by `freebayes -m 1`!
:::

---
# Random-Access Indexing (BGZF / BAI) & CRAM Compression

How did `samtools view` fetch **2.0 MB** for BRCA1 (`chr17:43.04–43.17 Mb`) out of a **multi-GB remote CRAM** in seconds?

![BGZF Independent Blocks and Hierarchical Binning Index](lecture_03_genetic_polymorphism_disease/images/bgzf_binning_index.png){width=79 center}

:::matrix { cols="1fr 1fr" gap="28px" margin-top: 8px; }
[[0,0]]
### BGZF + Hierarchical Binning (`.bai` / `.tbi` / `.crai`)
- **BGZF** packs coordinate-sorted records into independent $\le 64\text{ KB}$ gzip blocks addressed by 64-bit **virtual offsets** (`block_offset << 16 | intra_offset`).
- **UCSC / BAI binning** partitions each chromosome into 6 nested scales (`512 Mb` down to `16 kb`): querying `[a, b]` checks only $O(\log G)$ bins and issues an **HTTP byte-range request** for overlapping blocks!

[[0,1]]
### CRAM Reference-Based Delta Encoding
- Instead of storing every read's full DNA string, **CRAM** stores only **differences from the reference genome** (mismatches, indels, plus columnar quality streams) — cutting file size by **40–60%** relative to BAM.
- Because the reference bases are omitted, decoding **requires the exact reference FASTA** matching the header MD5 checksum (hence `samtools view -T ref.fa`).
:::

---

# PCR Duplicates & the Read Independence Assumption

In our BRCA1 BAM slice, **`samtools flagstat`** reported **4,789 PCR/optical duplicates (11%)**. Why does `freebayes` skip `FLAG 0x400` reads by default?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### The Hidden Assumption Behind $\prod_i P(b_i \mid G)$
Multiplying per-read likelihoods assumes every read $b_i$ was sampled from an **independent DNA molecule**:
$$P(D \mid G) = \prod_{i=1}^{n} P(b_i \mid G)$$
- During library prep, PCR copies the **same physical template molecule** $k$ times.
- If a polymerase error happens during an early PCR cycle (probability $\varepsilon = 10^{-3}$), **all $k$ daughter reads inherit that exact same false base**!
- Treating $k = 10$ copies of 1 error as independent shrinks the perceived error probability from **$10^{-3}$** to **$(10^{-3})^{10} = 10^{-30}$** — fabricating an astronomically confident false-positive variant!

[[0,1]]
### How Duplicate Marking Works (`FLAG 0x400`)
- **Paired-end signature:** Two read pairs that map to the **exact same 5' start coordinate, mate 5' coordinate, and strand orientation** almost certainly originated from the same sheared fragment.
- The pipeline keeps the highest base-quality representative and sets bit **`0x400` (1024)** on every remaining copy.
- **Critical distinction:**
  - Duplicates are **marked**, not deleted (preserving library complexity stats).
  - Both `freebayes` and our Homework Task 4 `pysam` loop **must skip `read.is_duplicate`** so evidence matches independent molecules!
:::

---

# Bayesian Diploid Read Classification

At any genomic locus in a diploid human genome, three biallelic genotypes compete: **`0/0`** (hom-ref), **`0/1`** (het), and **`1/1`** (hom-alt).

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Bayes' Theorem & MAP Calling
Combining the read likelihood $P(D \mid G)$ with a population prior $P(G)$ gives the **posterior genotype probability**:
$$P(G \mid D) = \frac{P(D \mid G)\ P(G)}{\sum_{G' \in \{0/0, 0/1, 1/1\}} P(D \mid G')\ P(G')}$$

The variant caller reports the **Maximum A Posteriori (MAP)** genotype:
$$\hat{G}\_{\text{MAP}} = \arg\max_{G} \bigl[\log P(D \mid G) + \log P(G)\bigr]$$

- **Likelihood $P(D \mid G)$:** *"How well does genotype $G$ explain the observed bases and their Phred quality scores?"*
- **Prior $P(G)$:** *"How common is genotype $G$ in the population before seeing these reads?"*

[[0,1]]
### Per-Read Two-Chromosome Mixture Model
Each independent read $b_i$ comes from either chromosome copy $a_1$ or $a_2$ with equal probability $\frac{1}{2}$:
$$P(b_i \mid G = a_1 / a_2) = \frac{P(b_i \mid a_1) + P(b_i \mid a_2)}{2}$$

With Phred base quality $Q_i$, the per-base sequencing error rate is:
$$\varepsilon_i = 10^{-Q_i / 10}$$

- `Q33` (`B`, 67% of BRCA1 bases): $\varepsilon = 5.0 \times 10^{-4}$ (strong evidence!)
- `Q6` (`'`, 12% of BRCA1 bases): $\varepsilon = 0.25$ (weak evidence, gently downweighted!)
- Notice: every single read's individual $Q_i$ enters the product — **no hard base-quality cutoff is required**!
:::

---

# Biallelic vs. Exact 4-Base Error Likelihoods

Plugging $P(b_i \mid a)$ into the two-chromosome mixture yields two closely related likelihood formulations:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### 1. Simplified Biallelic Model (Practice 3 Part 8)
Assumes an error on REF always flips to ALT (and vice versa): $P(\text{match}) = 1 - \varepsilon_i$, $P(\text{mismatch}) = \varepsilon_i$.

- **Hom-Ref (`0/0`):** REF reads match, ALT reads are errors:
  $$P(D \mid 0/0) = \prod_{i \in \text{ref}}(1 - \varepsilon_i) \prod_{j \in \text{alt}} \varepsilon_j$$
- **Het (`0/1`):** For *either* REF or ALT read, $\frac{(1 - \varepsilon_i) + \varepsilon_i}{2} = \frac{1}{2}$!
  $$P(D \mid 0/1) = \prod_{\text{all } i} \frac{1}{2} = \left(\frac{1}{2}\right)^n$$
- **Hom-Alt (`1/1`):** ALT reads match, REF reads are errors:
  $$P(D \mid 1/1) = \prod_{i \in \text{ref}} \varepsilon_i \prod_{j \in \text{alt}}(1 - \varepsilon_j)$$

[[0,1]]
### 2. Exact 4-Base Model (Industrial Callers)
DNA has **4 bases** — a sequencing error spreads across the **3 wrong bases**, so a *specific* mismatch has probability $\varepsilon_i / 3$:

| Genotype $G$ | $P(b_i = \text{REF} \mid G)$ | $P(b_i = \text{ALT} \mid G)$ | $P(b_i = \text{3rd base} \mid G)$ |
| :---: | :---: | :---: | :---: |
| **`0/0`** | $1 - \varepsilon_i$ | $\varepsilon_i / 3$ | $\varepsilon_i / 3$ |
| **`0/1`** | $\frac{1}{2} - \frac{\varepsilon_i}{3}$ | $\frac{1}{2} - \frac{\varepsilon_i}{3}$ | $\varepsilon_i / 3$ |
| **`1/1`** | $\varepsilon_i / 3$ | $1 - \varepsilon_i$ | $\varepsilon_i / 3$ |

- **Why the difference is small at high $Q$:** At `Q33`, $\frac{1}{2} - \frac{\varepsilon}{3} = 0.49983 \approx 0.5$, while $\varepsilon/3$ penalizes mismatching homozygotes slightly more strictly (by $\ln 3 \approx 1.1$ nats per error read).
- **3rd-base reads** cancel identically across all three genotypes!
:::

---

# Population Priors: Hardy–Weinberg Equilibrium

Without a prior (flat $P(G) = 1/3$), calling a rare homozygote `1/1` is treated as equally plausible as `0/0` — wildly inflating false positives across a 3 Gb genome where **$> 99.9\\%$ of sites are `0/0`**!

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Random-Mating HWE Derivation
If the alternate allele has population frequency $f \in (0, 1)$, drawing two independent chromosomes from the gene pool gives:

$$P(0/0) = (1 - f)^2, \quad P(0/1) = 2f(1 - f), \quad P(1/1) = f^2$$

| Frequency $f$ | $P(0/0)$ | $P(0/1)$ | $P(1/1)$ | Odds $P(0/1) / P(1/1)$ |
| :---: | :---: | :---: | :---: | :---: |
| `0.35` (common) | `0.4225` | `0.4550` | `0.1225` | **$3.7\times$** |
| `0.05` | `0.9025` | `0.0950` | `0.0025` | **$38\times$** |
| **`0.01` (Homework)** | **`0.9801`** | **`0.0198`** | **`0.0001`** | **$198\times$** |
| `0.001` (rare) | `0.9980` | `0.0020` | `10⁻⁶` | **$1{,}998\times$** |

[[0,1]]
### Exact Crossover Depth $n^*$ (`0/1` $\leftrightarrow$ `1/1`)
Suppose we observe **$0$ REF and $n$ pure ALT reads** (all high quality $\varepsilon \ll 1$). When does `1/1` overtake `0/1`?

Set $P(D \mid 1/1)P(1/1) = P(D \mid 0/1)P(0/1)$:
$$(1 - \varepsilon)^n f^2 = \left(\frac{1}{2}\right)^n \cdot 2f(1 - f)$$
$$\implies n^* = \frac{\ln\left(\frac{2(1 - f)}{f}\right)}{\ln\bigl(2(1 - \varepsilon)\bigr)} \approx \log_2 \frac{2(1 - f)}{f}$$

- At **$f = 0.01$**: $n^* \approx \log_2(198) = \mathbf{7.6}$ reads — so with only $2\text{–}7$ pure ALT reads, the HWE prior **still prefers `0/1`**!
- At **$f = 0.35$**: $n^* \approx \log_2(3.71) = \mathbf{1.9}$ reads — even $2$ ALT reads flip to **`1/1`**!
:::

---

# Prior vs. Data: Shallow vs. Deep Coverage

Strong read depth overwhelms any reasonable prior — whereas **shallow depth lets the prior decide the genotype call**:

![Prior vs Data Interplay and Exact Crossover Depth Phase Diagram](lecture_03_genetic_polymorphism_disease/images/bayesian_prior_vs_depth.png){width=96% center}

:::matrix { cols="1fr 1fr" gap="28px" margin-top: 6px; }
[[0,0]]
- **Left (`DP = 2`, Part 4 `rs12516`):** Two ALT reads give a likelihood ratio of only $2^2 = 4\times$ for `1/1` vs. `0/1` — overwhelmed by $198\times$ at $f = 0.01$ (`0/1`, 98%), yet flipped to `1/1` (52%) at $f = 0.35$!

[[0,1]]
- **Right (`DP ≥ 15`, Homework threshold):** Likelihood ratios exceed $2^{15} > 32{,}000\times$ — safely inside the data-dominated zone across all realistic frequencies $f \ge 10^{-4}$.
:::

---

# Numerical Stability: Log Space & Log-Sum-Exp

Multiplying raw probabilities across hundreds of reads causes catastrophic **IEEE-754 floating-point underflow**!

:::matrix { cols="1.05fr 0.95fr" gap="28px" }
[[0,0]]
![IEEE-754 Float64 Underflow Cliff vs Log-Likelihood](lecture_03_genetic_polymorphism_disease/images/log_underflow_lse.png){width=98% center}

[[0,1]]
### 1. Sum Log-Likelihoods Instead of Multiplying
At our BRCA1 max depth (`1,111×`), $(0.5)^{1111} \approx 3.4 \times 10^{-335}$ flushes to exact **`0.0`** in 64-bit IEEE-754 (`DBL_MIN ≈ 2.2 × 10⁻³⁰⁸`), producing **`0 / 0 = NaN`**!
$$\ell_g = \log P(D \mid g) + \log P(g)$$

### 2. The Log-Sum-Exp Normalization Trick
Subtract $m = \max_{g'} \ell_{g'}$ before exponentiating:
$$P(G = g \mid D) = \frac{e^{\ell_g - m}}{\sum_{g'} e^{\ell_{g'} - m}}$$
- The winning term is guaranteed $e^0 = \mathbf{1}$ (zero total underflow, zero overflow!).

### 3. Boundary Clamping (`Q = 0` Crash!)
When a base has `Q = 0`, $\varepsilon = 10^0 = 1 \implies \log(1 - \varepsilon) = \log(0) = -\infty$. Always clip `eps = np.clip(eps, 1e-12, 1 - 1e-12)`!
:::

---

# Interactive Lab — Bayesian Genotyper & HWE Prior

Test real BRCA1 scenarios live: **(1)** inspect `rs3092994` (`21 REF / 19 ALT`) and verify `GL` and `QUAL`; **(2)** load **Prior Flip** (`0 REF / 2 ALT`) and drag $f$ across `0.35` to watch `0/1` flip to `1/1`; **(3)** load **High QUAL / Low GQ Paradox** (`0 REF / 5 ALT`); **(4)** load **Q = 0 Trap** and toggle `Clip ε`; **(5)** load **Float64 Underflow (`1,100×`)**!

<iframe src="lecture_03_genetic_polymorphism_disease/demo/bayesian_genotyper_lab.html" title="Interactive Bayesian Diploid Genotyper and HWE Prior Lab" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Confidence Metrics: Site `QUAL` vs. Sample `GL` & `GQ`

A common variant-calling mistake is treating VCF **`QUAL`** as genotype certainty. They answer **different probability questions**:

![Site QUAL vs Genotype GQ Paradox Across Pure ALT Depth](lecture_03_genetic_polymorphism_disease/images/qual_vs_gq_paradox.png){width=96% center}

:::matrix { .smaller-table .dense cols="1fr" margin-top: 6px; }
[[0,0]]
| Metric | Exact Definition | Question Answered | Behavior at `0 REF, 5–8 ALT` ($f = 0.01$) |
| :--- | :--- | :--- | :--- |
| **Site `QUAL`** | $-10\log_{10} P(G = 0/0 \mid D)$ | *"Is there **any** variant at this site?"* | **Huge (`119 – 204`)** — virtually zero chance of `0/0`! |
| **Sample `GL`** | $\log_{10} P(D \mid G) - \max_{G'} \log_{10} P(D \mid G')$ | *"Relative data likelihood of each genotype?"* | Pure likelihoods (prior-free), best genotype scaled to `0` |
| **Sample `GQ`** | $-10\log_{10}\bigl(1 - P(\hat{G} \mid D)\bigr)$ | *"How sure is the **specific** winning genotype?"* | **Dips to `3.6`** at $n^* \approx 7.6$ where `0/1` and `1/1` tie! |
:::

---

# Haplotype-Based Calling (`freebayes`)

Why did `freebayes` emit alleles like `TAAAAAAAT > TAAAAAAT` instead of calling single bases one column at a time?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Why Column-by-Column Pileup Callers Fail
1. **Clustered SNPs (MNPs):** Two mutations 2 bp apart on the *same* chromosome (e.g., `AT → GC`) degrade alignment scores; independent single-column callers lose phase information and mis-annotate the resulting codon change!
2. **Indels near repeats:** Different reads place the same deletion at slightly different columns across a homopolymer — splitting read support across multiple columns so *no single column* passes the variant threshold!

```
Ref:     G C T A A A A A A A T C
Read 1:  G C T - A A A A A A T C   (del at col 4)
Read 2:  G C T A A A A A A - T C   (del at col 10)
             ▲ pileup sees 2 weak errors instead of 1 real deletion!
```

[[0,1]]
### How Haplotype Calling Solves It
- Instead of slicing independent 1-bp columns, **`freebayes`** scans a short **local genomic window** containing clustered mismatches/indels.
- Every spanning read votes as a **complete local haplotype observation** $h_i \in \mathcal{H}$ (e.g., full string `TAAAAAAAT` vs. `TAAAAAAT`).
- **Bayesian genotyping runs directly over haplotype pairs** $G = h_1 / h_2$:
  - Unifies inconsistently aligned indel reads into **one clean allele count** (`RO`, `AO`).
  - Preserves **cis vs. trans phasing** of neighboring variants on the same read.
  - Trade-off: complex haplotypes must be **normalized / decomposed** before matching against ClinVar!
:::

---
# Canonical Indel Normalization (Tan, Abecasis & Kang, 2015)

`freebayes` and **ClinVar** may disagree on **all 16** BRCA1 indels because the same biological mutation has many equivalent VCF strings!

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Uniqueness Theorem (Two Necessary & Sufficient Invariants)
A VCF record `(POS, REF, ALT)` is **canonical** (guaranteed unique across every caller and database!) if and only if it is:
1. **Parsimonious:**
   - The rightmost bases of `REF` and `ALT` **differ** (no redundant trailing context), AND
   - If both lengths $\ge 2$, the leftmost bases **differ** (no redundant leading context beyond a single anchor base).
2. **Left-Aligned:**
   - The variant cannot be shifted left by 1 bp while keeping the exact same mutated genome sequence.

### 3-Step Canonicalization Algorithm
1. **Right-trim & left-shift loop:** While last base of `REF` == last base of `ALT`: if either length is 1, prepend the preceding reference base (`POS -= 1`); then drop the last base from both.
2. **Left-trim loop:** While both lengths $\ge 2$ and first base of `REF` == first base of `ALT`: drop first base from both and increment `POS += 1`.

[[0,1]]
### Trace on Real BRCA1 Record (`rs8176144`)
Reference around `43097344`: `G C T A A A A A A A T C`

| Step | Operation | `POS` | `REF` | `ALT` |
| :-: | :--- | :---: | :--- | :--- |
| 0 | `freebayes` raw | `43097346` | `TAAAAAAAT` | `TAAAAAAT` |
| 1 | Right-trim `T` | `43097346` | `TAAAAAAA` | `TAAAAAA` |
| 2–7 | Right-trim `A` ($\times 6$) | `43097346` | **`TA`** | **`T`** |
| Stop | Last bases `A ≠ T`, `len(ALT)=1` | `43097346` | **`TA`** | **`T`** |

- Every intermediate record produces the **exact same mutated genome sequence** (`GCTAAAAAATC`)!
- Matches ClinVar **`rs8176144` (`Benign`, expert panel)** verbatim (`bcftools norm -f ref.fa -a -m -any`).
:::

---

# Interactive Lab — Indel Normalization & BGZF/BAI Index

Two hands-on tabs: **(1)** step through **Tan et al.'s** right-trim and repeat left-shift on real BRCA1 `rs8176144` and right-shifted tandem repeats, verifying that the mutated genome stays invariant at every step; **(2)** switch to **BGZF & BAI Lookup**, jump to **BRCA1 (`43.04–43.17 Mb`)**, and inspect the exact SAMv1 bin IDs and bytes skipped!

<iframe src="lecture_03_genetic_polymorphism_disease/demo/haplotype_indel_norm.html" title="Interactive Indel Normalization and BGZF Binning Explorer" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Coverage Statistics, Allelic Dropout & AAF Dispersion

Even with zero sequencing errors ($\varepsilon = 0$), finite read depth $n$ creates **binomial sampling noise** at every heterozygous site ($k \sim \operatorname{Binomial}(n, \tfrac{1}{2})$):

![Heterozygote Allelic Dropout and Binomial AAF Distribution](lecture_03_genetic_polymorphism_disease/images/dropout_and_aaf_spread.png){width=96% center}

:::matrix { cols="1fr 1fr" gap="28px" margin-top: 6px; }
[[0,0]]
- **Allelic dropout ($k = 0$ or $k = n$):** $P_{\text{drop}} = 2\left(\frac{1}{2}\right)^n$. At BRCA1's off-target median `DP = 4`, **$12.5\\%$** of heterozygotes drop an allele! Requiring **`DP ≥ 15`** slashes dropout to **$0.006\\%$**.

[[0,1]]
- **AAF dispersion:** $\text{AAF} = \frac{AO}{RO + AO}$ has $\operatorname{SD}(\text{AAF} \mid 0/1) = \frac{1}{2\sqrt{n}}$. At `DP = 40`, $\pm 2\sigma$ spans `[0.34, 0.66]` — so `21/19` and `25/15` are completely normal heterozygotes!
:::

---

# Beyond Ideal Binomial Sampling: Four Real-World Biases

Real exome BAMs deviate from an exact $\operatorname{Binomial}(n, 0.5)$ coin flip in four recognizable ways:

:::matrix { .cards cols="1fr 1fr" gap="24px" margin-top: 12px; }
[[0,0]]{ .purple }
### 1. Reference Mapping & Capture Bias
Reads carrying the **ALT allele** have an extra mismatch against the reference genome (lower BWA-MEM score) and hybridize slightly less tightly to reference-designed exome baits. This systematically shifts heterozygote AAF peaks **below `0.50`** (typically **`0.46 – 0.48`**, e.g. `21/40 = 0.475` at `rs3092994`).

[[0,1]]{ .orange }
### 2. Strand Bias (Systematic Artifacts)
True biological variants are sampled roughly **50/50** from forward and reverse strands. When all ALT reads come from **one strand only** (like Part 4's `chr17:43057017`: `DP = 375`, 80/80 ALT reads on the forward strand!), it flags a context-specific polymerase/sequencing artifact — Fisher's exact test rejects it (`QUAL = 0.0`).

[[1,0]]{ .blue }
### 3. Read-End & Context Error Clustering
Illumina base quality decays toward the **3' tail** (Lecture 1) and near soft-clipped indel boundaries. Callers weight by **base quality** and **read-position placement** so mismatches clustered exclusively at the last 3 bp of reads are downweighted.

[[1,1]]{ .green }
### 4. Germline vs. Somatic / Mosaic AAF Peaks
In a pure diploid germline sample (NA12878), filtered variants cluster cleanly near **`0.5`** (`0/1`) and **`1.0`** (`1/1`). Persistent well-supported variants at **`0.05 – 0.30`** signal either sample **contamination**, copy-number gains, or **somatic subclonal tumor mutations** (where diploid HWE assumptions no longer hold!).
:::

---

# Interactive Lab — Coverage, Dropout & AAF Spectrum

Explore sampling noise across depths: **(1)** click `DP = 4`, `15`, and `40` to compare allelic dropout $2(1/2)^n$ and the `freebayes -F 0.2` missed-het tail; **(2)** drag `p_alt` down to `0.45` to observe reference mapping bias; **(3)** switch to **Tab 2 (Cohort AAF Spectrum)** and toggle **Raw Off-Target (`DP ≈ 5`)** vs. **Filtered Exome (`DP ≈ 40`)** to preview Homework Task 2!

<iframe src="lecture_03_genetic_polymorphism_disease/demo/coverage_aaf_explorer.html" title="Interactive Coverage Dropout and AAF Spectrum Simulator" style="width: 100%; height: 590px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>

---

# Algorithmic Complexity & Pipeline Scaling

Every stage of the clinical variant-calling pipeline is engineered to scale across **3 billion base pairs** and **hundreds of millions of reads**:

:::matrix { .smaller-table cols="1fr" margin-top: 8px; }
[[0,0]]
| Pipeline Stage | Algorithm / Data Structure | Time Complexity | Memory / I/O Footprint | Why Naive Fails |
| :--- | :--- | :---: | :---: | :--- |
| **Read Alignment** | BWT / FM-index + banded SW (`BWA-MEM`) | $O(L)$ seed + banded ext | ~4.5 GB (human index) | Full SW is $O(GL) \approx 10^{19}$ ops |
| **Region Streaming** | BGZF 64 KB blocks + BAI/Tabix bins (`samtools`) | $O(\log G + K_{\text{slice}})$ | Fetches ~2 MB of 100 GB | Standard `.gz` must read from byte 0 |
| **Duplicate Marking** | Coordinate-sorted 5' pair hashing (`FLAG 0x400`) | $O(R)$ single pass | $O(\text{active window})$ | Correlated copies fabricate $10^{-30}$ calls |
| **Bayesian Genotyping** | Log-space HWE MAP + Log-Sum-Exp (`freebayes`) | $O(D \cdot |\mathcal{H}|^2)$ per site | $O(D)$ streaming window | Raw products underflow at $D > 1{,}024$ |
| **ClinVar Lookup** | Tabix `.tbi` index + canonical `(POS,REF,ALT)` | $O(\log B)$ per query | Constant streaming RAM | Unnormalized indels miss 16/16 matches |
:::

:::matrix { .cards cols="1fr 1fr" gap="24px" margin-top: 12px; }
[[0,0]]{ .blue }
### Why `cyvcf2` & `pysam` Wrap C (`htslib`)
Python loops over 40,000+ BAM alignments or 200 MB ClinVar VCFs would be 50–100× slower in pure Python. Both libraries delegate BGZF decompression, binary decoding, and bin-index seeking directly to compiled C (`htslib`).

[[0,1]]{ .green }
### Streaming Invariant Across All Stages
From `samtools view` $\to$ `freebayes` $\to$ `cyvcf2.VCF` $\to$ `pysam.fetch()`, every file is **coordinate-sorted**, enabling bounded-memory linear streaming across whole human chromosomes!
:::

---

# Summary & Bridge to Part 4 (Industrial Tools)

Every Part 3 algorithm directly powers a concrete tool and task in **Part 4** and **Practice 3**:

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### Alignment & Storage
#### `samtools` · `pysam`
- FM-index $O(|seed|)$ + `MAPQ` ($S_1 - S_2$)
- BGZF 64 KB blocks + BAI/Tabix bins
- CRAM reference delta (`-T ref.fa`)
- Skip `FLAG 0x400` duplicates to preserve independence

[[0,1]]{ .orange }
### Bayesian Genotyping
#### `freebayes` · `bayesian_genotyper`
- Two-chromosome mixture + $\varepsilon_i = 10^{-Q_i/10}$
- HWE priors $(1-f)^2,\ 2f(1-f),\ f^2$
- Crossover depth $n^* \approx \log_2 \frac{2(1-f)}{f}$
- Log-sum-exp + $\varepsilon$ clipping (`Q = 0`)

[[0,2]]{ .blue }
### Confidence, Indels & QC
#### `cyvcf2` · `ClinVar` · `bcftools norm`
- Site `QUAL` $\ne$ sample `GL` / `GQ`
- Haplotype calling + Tan et al. left-alignment
- Dropout $2(1/2)^n$ $\implies$ `DP ≥ 15`
- AAF $\operatorname{SD} = \frac{1}{2\sqrt{n}}$ + ref/strand bias
:::

> **Quick-Reference Core Equations:**  
> $\varepsilon_i = 10^{-Q_i/10}$ · $P(G \mid D) = \frac{e^{\ell_g - m}}{\sum_{g'} e^{\ell_{g'} - m}}$ · $n^* \approx \log_2 \frac{2(1-f)}{f}$ · $\text{QUAL} = -10\log_{10} P(0/0 \mid D)$ · $\text{GQ} = -10\log_{10}\bigl(1 - P(\hat{G} \mid D)\bigr)$ · $P_{\text{drop}} = 2(1/2)^n$
