:::titlepage
[[top]]
Part 3

[[title]]
Underlying Algorithms

[[subcaption]]
Quality Filtering · $k$-mers · MinHash · Mash Distance
:::

---

# Curriculum & Algorithmic Roadmap

How gigabytes of noisy raw reads become a 1 KB genomic fingerprint and a statistically calibrated species call:

:::static-diagram { width=96% center }
{
  "nodes": [
    { "id": "fq",    "label": "Raw FASTQ Stream\nASCII Phred+33",        "x": 20,  "y": 70,  "width": 165, "height": 62 },
    { "id": "qc",    "label": "fastp QC & Overlap\nWindow Trim + Adapter", "x": 225, "y": 70,  "width": 180, "height": 62, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "kmer",  "label": "Canonical k-mers\n2-Bit Rolling Hash",     "x": 445, "y": 70,  "width": 170, "height": 62, "style": "fill: #2196F3; color: #fff;" },
    { "id": "bloom", "label": "Bloom Filter Gate\nMultiplicity (-m 2)",   "x": 655, "y": 70,  "width": 175, "height": 62, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "heap",  "label": "MurmurHash3 + Heap\nBottom-s Sketch",      "x": 870, "y": 70,  "width": 180, "height": 62, "style": "fill: #FF9800; color: #fff;" },
    { "id": "jacc",  "label": "O(s) Sorted Merge\nJaccard + G-hat",       "x": 655, "y": 210, "width": 175, "height": 62, "style": "fill: #FF9800; color: #fff;" },
    { "id": "mash",  "label": "Poisson Mash Dist\nD ≈ 1 − ANI",           "x": 445, "y": 210, "width": 170, "height": 62, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "pval",  "label": "Binomial Null Test\np-value in 4^k Space", "x": 225, "y": 210, "width": 180, "height": 62, "style": "fill: #2196F3; color: #fff;" }
  ],
  "links": [
    { "from": "fq",    "to": "qc" },
    { "from": "qc",    "to": "kmer" },
    { "from": "kmer",  "to": "bloom" },
    { "from": "bloom", "to": "heap" },
    { "from": "heap",  "to": "jacc" },
    { "from": "jacc",  "to": "mash" },
    { "from": "mash",  "to": "pval" }
  ]
}
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
- **Preprocessing (`fastp`)**: ASCII Phred decoding, $O(L)$ sliding-window trimming, and de novo paired-end adapter overlap.
- **$k$-mers & Error Suppression**: Canonical strand invariance, the Goldilocks trade-off for $k$, Bloom filters (`-m 2`), and hash-density genome sizing.

[[0,1]]
- **MinHash Sketching (`mash`)**: 64-bit MurmurHash3, bottom-$s$ max-heap maintenance, and $O(s)$ unbiased Jaccard merging.
- **Distance, Limits & Complexity**: Poisson conversion to $D \approx 1 - \text{ANI}$, breakdown regimes, null $p$-values, and $O(1)$-RAM streaming I/O.
:::

---

# Phred Decoding & Sliding-Window Quality Trimming

Illumina chemistry accumulates **signal dephasing** across sequencing cycles, causing base-call error rates to climb sharply near the $3'$ end of each read.

:::matrix { cols="1.05fr 1fr" gap="30px" }
[[0,0]]
### ASCII Phred+33 Decoding
Every quality character $c$ encodes a logarithmic error probability $P_{\text{error}}$ in a single printable byte (offset by $33$ to skip ASCII control characters):

$$Q = \operatorname{ord}(c) - 33 = -10\log_{10} P_{\text{error}}$$

$$P_{\text{error}} = 10^{-Q / 10}$$

- `Q10` (`+`, ASCII 43): $P_{\text{error}} = 10\\%$ (1 in 10 wrong)
- `Q20` (`5`, ASCII 53): $P_{\text{error}} = 1\\%$ (99% accurate)
- `Q30` (`?`, ASCII 63): $P_{\text{error}} = 0.1\\%$ (99.9% accurate)
- `Q40` (`I`, ASCII 73): $P_{\text{error}} = 0.01\\%$ (99.99% accurate)

[[0,1]]
### Two Complementary `fastp` Filters
1. **Whole-read unqualified filter (default `-q 15 -u 40`)**:
   - Count bases with $Q_j < 15$; discard the entire read if that fraction exceeds $40\\%$ (or if `N` count $> 5$).
2. **$O(L)$ Sliding-Window Tail Trim (`--cut_right`)**:
   - Maintain a rolling sum $S_i = \sum_{j=i}^{i+w-1} Q_j$ over width $w$ (default $w = 4$).
   - Update in $O(1)$ per base: $S_{i+1} = S_i - Q_i + Q_{i+w}$.
   - Trigger cut at the **first** $i$ where:
     $$\bar{Q}\_i = \frac{S_i}{w} < Q_{\text{threshold}}$$
   - Discard positions $[i, L)$ cleanly before $k$-mer extraction.
:::

---

# De Novo Adapter Detection & Library Diagnostics

How does `fastp` strip sequencing adapters **without** being told which kit or adapter sequence was used?

:::matrix { cols="1.15fr 1fr" gap="30px" }
[[0,0]]
### Paired-End Overlap Geometry ($I < L$)
When a DNA fragment's physical insert length $I$ is shorter than the read length $L$, **both** mates read past the insert into the opposite adapter:

```text
Physical DNA:  [Adapter A]═════ Insert (I bp) ═════[Adapter B]
Read 1 (5'→3'):           ═════════════════════════[Adapter B..]
revcomp(R2):   [..Adapter A]═══════════════════════
                          ◄──── Shared Overlap ────►
```

- `fastp` reverse-complements $R_2$ and performs an **ungapped shift scan** against $R_1$ (requiring overlap $\ge 30\text{ bp}$ and mismatch rate $\le 10\\%$).
- A sharp identity peak at negative shift $I - L < 0$ reveals the exact insert boundary $I$ **de novo** — bases protruding past the shared overlap are clipped from both reads!
- **Bonus**: inside the shared overlap, mismatched bases are error-corrected toward the mate with the higher Phred score.

[[0,1]]
### Single-Pass Library Diagnostics
While streaming each read, `fastp` computes two additional quality signals in $O(L)$ time:

- **Low-Complexity Filter (`-y`)**:
  Measures the fraction of adjacent bases that differ:
  $$C = \frac{1}{L - 1}\sum_{i=1}^{L-1} \mathbf{1}[s_i \ne s_{i+1}]$$
  Homopolymer runs (e.g. `GGGGGGGG` from NovaSeq two-color dark cycles) have $C \approx 0$ and are stripped or rejected.
- **PCR Duplication Estimation**:
  Hashes fixed-length read prefixes (e.g. first 32 bp of $R_1$ and $R_2$) into a bounded hash table to estimate the fraction of clonal PCR duplicates in $O(1)$ RAM.
:::

---

# Interactive Lab — `fastp` QC & Adapter Overlap

:::matrix { cols="0.82fr 1.55fr" gap="24px" }
[[0,0]]
### Guided Exploration
1. **Sliding-Window Tab**:
   - Hover any bar to verify ASCII `ord(c) - 33 → Q → P_err`.
   - Click **▶ Sweep** to watch `--cut_right` stop at the first window with $\bar{Q}_i < 20$.
   - Switch to **Mid-Read Dip** and increase window width $W$ from `1` to `8`: see how wider windows smooth isolated single-base glitches!
2. **Adapter Overlap Tab**:
   - Drag **Insert Size $I$** from `120 bp` down to `70 bp` (below $L = 100\text{ bp}$).
   - Watch the red overhangs appear and the **shift-overlap peak** cross into negative shifts ($I - L < 0$).

<button class="demo-btn" style="margin-top: 14px; font-size: 1.1rem; padding: 10px 18px;" onclick="showDemo('lecture_01_the_molecular_blueprint/demo/fastp_qc_overlap.html')">⛶ Full-Screen Lab</button>

[[0,1]]
<iframe src="lecture_01_the_molecular_blueprint/demo/fastp_qc_overlap.html" title="fastp QC and Adapter Overlap Lab" style="width: 100%; height: 615px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>
:::

---

# $k$-mer Decomposition & Canonical $k$-mers

Instead of slow character-by-character sequence alignment, modern genomics represents a genome as the **set of all length-$k$ substrings** ($k$-mers).

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Exact Decomposition & Strand Ambiguity
- A read of length $L$ yields exactly **$L - k + 1$** overlapping $k$-mers (a genome of length $G$ yields $\approx G$ $k$-mers).
- **Double-stranded DNA problem**: a sequencer reads either the forward ($5' \to 3'$) or reverse-complement strand at random!
- Without strand normalization, a genome and its exact reverse complement would share **zero** $k$-mers.
- **Solution — Canonical $k$-mers**: map each $k$-mer $x$ and its reverse complement $\overline{x}^{\text{rev}}$ to the lexicographically smaller string:
  $$\operatorname{canon}(x) = \min\bigl(x,~ \overline{x}^{\text{rev}}\bigr)$$
  - Example ($k=5$): $x = \texttt{TGCTA}$, $\overline{x}^{\text{rev}} = \texttt{TAGCA} \implies \operatorname{canon} = \texttt{TAGCA}$.
  - Odd $k$ (e.g. $k = 21$) prevents palindromic self-complements ($x = \overline{x}^{\text{rev}}$).

[[0,1]]
### 2-Bit Rolling Encoding ($O(1)$ per Base)
For $k \le 32$, each nucleotide fits in **2 bits** (`A=00`, `C=01`, `G=10`, `T=11`), packing a full $k$-mer into a single 64-bit integer:

```text
Forward integer:  f = ((f << 2) | b) & MASK_2K
Rev-comp integer: r = (r >> 2) | ((3 ^ b) << (2k - 2))
Canonical code:   c = min(f, r)
```

- **Zero string slicing or allocations** inside the inner loop!
- Updating both forward and reverse-complement 21-mers takes **~4 CPU bit operations per base**.
- Any `N` base simply resets the rolling counter for the next $k$ bases.
:::

---

# Choosing $k$: Random Collisions vs. Mutation Sensitivity

Every $k$-mer algorithm faces a fundamental **Goldilocks trade-off** between alphabet capacity ($4^k$) and sensitivity to mutations.

:::matrix { cols="1.25fr 1fr" gap="25px" }
[[0,0]]
![Goldilocks k-mer Tradeoff](lecture_01_the_molecular_blueprint/images/kmer_tradeoff.png){width=98 center}

[[0,1]]
### 1. Lower Bound — Avoid Random Collisions
In a random genome of size $G$, the probability a given $k$-mer appears by pure chance is:
$$q = 1 - \left(1 - \frac{1}{4^k}\right)^G \approx \frac{G}{4^k} \implies k \gtrsim \log_4\left(\frac{G}{q}\right)$$

### 2. Upper Bound — Survive Mutations
One single-base substitution (or error) destroys **$k$ consecutive $k$-mers**! At divergence $d$, survival decays exponentially:
$$p = (1 - d)^k \approx e^{-kd}$$

| Organism | Genome $G$ | Min $k$ ($q=0.01$) | Typical $k$ |
| :--- | :---: | :---: | :---: |
| Virus | $10^4$ bp | 10 | **11 – 15** |
| Bacteria | $5\times 10^6$ bp | 15 | **21** (Mash) |
| Mammal | $3\times 10^9$ bp | 19 | **21 – 31** |
:::

---
# Suppressing Error $k$-mers: Bloom Filters & `-m 2`

Raw FASTQ reads contain base-calling errors — even a modest **$0.5\\%$** error rate at $k = 21$ corrupts $1 - 0.995^{21} \approx 10\\%$ of **every** extracted $k$-mer!

![Read k-mer Spectrum and Bloom Filter Gate](lecture_01_the_molecular_blueprint/images/kmer_spectrum_bloom.png){width=94% center}

---
# Suppressing Error $k$-mers: Bloom Filters & `-m 2`
:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Why Singleton Errors Break MinHash
- Each random error creates $k$ **unique** spurious $k$-mers at multiplicity $c = 1$ — often **outnumbering** true genomic $k$-mers $3:1$!
- Those spurious hashes flood the bottom-$s$ sketch, inflating the union $|A \cup B|$ and wrecking Jaccard accuracy.
- Filtering to **`-m 2`** (require $\ge 2$ sightings) cleanly slices off the red singleton spike while preserving the green genomic peak.

[[0,1]]
### The Bloom Filter Gate ($O(m)$ Fixed RAM)
- Storing exact counts for 20 million $k$-mers wastes gigabytes of RAM.
- Instead, `mash` passes each $k$-mer through a **Bloom filter** ($m$-bit array, $r$ hash functions):
  - **1st sighting**: at least one bit is `0` $\to$ flip all $r$ bits to `1`, **do not sketch**.
  - **2nd sighting**: all $r$ bits are already `1` $\to$ **admit to MinHash**!
- **Pitfall**: never pass `-m 2` on a 1× assembled FASTA genome — every unique $k$-mer appears once and gets discarded!
:::

---

# Genome Size & Coverage from Bottom-$s$ Hash Density

How can `mash sketch` print `Estimated genome size: 4.64e+06` and `Estimated coverage: 47.13` **without ever keeping a table of all $k$-mers**?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Order-Statistics Derivation
Suppose a genome contains $G$ distinct valid $k$-mers, and a uniform hash maps each independently onto $[0, H_{\max}]$ (e.g. $H_{\max} = 2^{64} - 1$).

- Normalized hashes $U_1, \dots, U_G \stackrel{\text{i.i.d.}}{\sim} \operatorname{Uniform}(0, 1)$ divide the interval $[0, 1]$ into $G + 1$ gaps of equal expected length $\frac{1}{G + 1}$!
- Therefore, the **$s$-th smallest** hash $h_{(s)}$ obeys:
  $$\mathbb{E}\left[\frac{h_{(s)}}{H_{\max}}\right] = \frac{s}{G + 1}$$
- Solving directly for the unknown number of unique $k$-mers $G$:
  $$\hat{G} \approx \frac{s \cdot H_{\max}}{h_{(s)}}$$

[[0,1]]
### Intuition on the Unit Interval
```text
0 |• • • • • • •|───────────────────────────────| 1
  ◄── s hashes ─►                                 
                ▲ h_(s) / H_max
```
- **Small genome** (virus, $G = 10^4$): hashes are sparse $\to$ $h_{(s)}$ sits far to the **right**.
- **Large genome** (bacteria, $G = 5\times 10^6$): hashes pack tightly near $0$ $\to$ $h_{(s)}$ is tiny!

### Sequencing Coverage Estimate
Once $\hat{G}$ is estimated from the $s$-th smallest hash, mean depth follows immediately from total bases streamed:
$$\widehat{\text{Coverage}} = \frac{\text{Total bases streamed}}{\hat{G}}$$
*(Relative standard error of $\hat{G}$ is $\approx 1/\sqrt{s} \approx 3\\%$ at $s = 1000$!)*
:::

---

# Uniform Hashing (MurmurHash3) & the MinHash Principle

How do we compare two sets of $5\times 10^6$ $k$-mers using only **1,000 integers** — with mathematical guarantees?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Why MurmurHash3 (64-bit)?
- Maps arbitrary canonical $k$-mers uniformly across $[0, 2^{64}-1]$ (~$1.8 \times 10^{19}$ slots $\implies$ near-zero hash collisions).
- **Strong avalanche effect**: flipping a single DNA base flips $\approx 50\\%$ of the output bits, destroying sequence correlations.
- **Non-cryptographic**: $10\times$ faster than SHA-256 (we need statistical uniformity, not cryptographic pre-image resistance).

[[0,1]]
### Broder's MinHash Theorem (1997)
For any two finite $k$-mer sets $A$ and $B$, let $h$ induce a uniform random permutation over $A \cup B$. Which element achieves the global minimum $\min h(A \cup B)$?

- By symmetry, **every element** in $A \cup B$ is equally likely to hold the minimum hash!
- The minimum of $A$ equals the minimum of $B$ **if and only if** that global minimum element lies in $A \cap B$:
  $$P\bigl(\min h(A) = \min h(B)\bigr) = \frac{|A \cap B|}{|A \cup B|} = J(A, B)$$

> One minimum hash is an unbiased **Bernoulli($J$) coin flip** for the exact Jaccard similarity!
:::

---
# Bottom-$s$ Sketching with a Max-Heap

Repeating $s = 1000$ independent hash functions would scan every genome $1000\times$ slower. **Bottom-$s$ MinHash** evaluates just **one** hash function per $k$-mer and retains its **$s$ smallest** distinct values.

:::static-diagram { width=94% center }
{
  "nodes": [
    { "id": "in",   "label": "Canonical k-mer\n(after -m 2)",    "x": 20,  "y": 60, "width": 160, "height": 58 },
    { "id": "hash", "label": "MurmurHash3\n64-bit h",            "x": 225, "y": 60, "width": 150, "height": 58, "style": "fill: #FF9800; color: #fff;" },
    { "id": "cmp",  "label": "Compare h < root?\n(99.9% skip!)", "x": 420, "y": 60, "width": 185, "height": 58, "style": "fill: #2196F3; color: #fff;" },
    { "id": "heap", "label": "Max-Heap (size s)\nRoot = h_max",  "x": 650, "y": 60, "width": 175, "height": 58, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "sort", "label": "Sort s Hashes\n.msh (~1 KB)",      "x": 870, "y": 60, "width": 155, "height": 58, "style": "fill: #4CAF50; color: #fff;" }
  ],
  "links": [
    { "from": "in",   "to": "hash" },
    { "from": "hash", "to": "cmp" },
    { "from": "cmp",  "to": "heap", "label": "h < h_max" },
    { "from": "heap", "to": "sort" }
  ]
}
:::

:::matrix { cols="1.15fr 1fr" gap="30px" }
[[0,0]]
### Why a *Max*-Heap for the *Smallest* Hashes?
```python
# Maintain the s smallest distinct hashes in O(s) RAM
for kmer in valid_canonical_kmers(stream):
    h = murmurhash3_64(kmer)
    if h in seen_in_heap: continue           # deduplicate
    if len(heap) < s:
        maxheap_push(heap, h); seen_in_heap.add(h)
    elif h < heap[0]:     # heap[0] = current s-th smallest!
        evicted = maxheap_replace_root(heap, h)
        seen_in_heap.remove(evicted); seen_in_heap.add(h)
```

[[0,1]]
### Runtime & Memory Analysis
- **Memory**: $O(s)$ — only $1000 \times 8\text{ B} = 8\text{ KB}$ (plus hash set).
- **Worst-case time**: $O(N \log s)$ if hashes arrive in decreasing order.
- **Expected time (random hash order)**: the $i$-th unique $k$-mer enters the bottom-$s$ with probability $s/i$! Total heap updates across $N$ $k$-mers:
  $$\sum_{i=s+1}^{N} \frac{s}{i} \approx s \ln\left(\frac{N}{s}\right) \approx 8{,}500 \text{ swaps for } N = 5\times 10^6$$
  Nearly **$O(N)$ linear time** — 99.8% of $k$-mers are rejected in a single $O(1)$ comparison!
:::

---
# Interactive Lab — Live MinHash & Max-Heap Sketcher

:::matrix { cols="0.82fr 1.55fr" gap="24px" }
[[0,0]]
### Guided Exploration
1. **Strand Test**:
   - Click **Reverse Complement**, then uncheck **Canonical k-mers**: Jaccard crashes to **`0.000`**!
   - Re-check **Canonical k-mers**: it snaps straight back to **`1.000`** ($D = 0$).
2. **Watch the Max-Heap**:
   - Click **↺ Reset** then **▶ Play Stream**.
   - Notice the purple root (`h_max`) shrinking steadily leftward along `[0, 1)` while most incoming $k$-mers are rejected in $O(1)$.
3. **Check Cardinality & Merge**:
   - Compare **Est. $\hat{G} = s / h_{(s)}$** vs. true unique $|A|$, and verify green shared chips in the merged bottom-$s$ set $W$.

<button class="demo-btn" style="margin-top: 12px; font-size: 1.1rem; padding: 10px 18px;" onclick="showDemo('lecture_01_the_molecular_blueprint/demo/minhash_sketcher.html')">⛶ Full-Screen Lab</button>

[[0,1]]
<iframe src="lecture_01_the_molecular_blueprint/demo/minhash_sketcher.html" title="Live MinHash and Max-Heap Sketcher" style="width: 100%; height: 650px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>
:::

---

# Unbiased Jaccard Estimation & Standard Error

Once sketches $S(A)$ and $S(B)$ are saved as sorted arrays of $s$ 64-bit integers, comparing two entire genomes takes **microseconds** — no sequences needed!

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Two-Pointer $O(s)$ Merged-Sketch Estimator
1. Merge sorted arrays $S(A)$ and $S(B)$ with two pointers (like merge-sort) and stop after collecting the **$s$ smallest distinct hashes** of their union:
   $$W = \min_s\bigl(S(A) \cup S(B)\bigr)$$
2. Count how many elements of $W$ appear in **both** $S(A)$ and $S(B)$:
   $$\hat{J}(A, B) = \frac{|W \cap S(A) \cap S(B)|}{|W|} = \frac{x}{s}$$

> **Why not just compute $|S(A) \cap S(B)| / s$?** Because $S(A)$ and $S(B)$ can have different cutoff radii $h_{(s)}$! Truncating to the common bottom-$s$ window $W$ of $A \cup B$ guarantees that every element of $W$ is a uniform sample from $A \cup B$ — making $\hat{J}$ strictly **unbiased**.

[[0,1]]
### Exact Variance & Sketch Size Scaling
Each of the $s$ draws in $W$ lands in $A \cap B$ with probability $J$ (without replacement from $|A \cup B| \gg s$, so tightly Binomial($s, J$)):

$$\operatorname{SE}(\hat{J}) = \sqrt{\frac{J(1 - J)}{s}} \le \frac{1}{2\sqrt{s}} = O\left(\frac{1}{\sqrt{s}}\right)$$

| Sketch Size $s$ | Max SE ($\hat{J}$) | `.msh` Size | Typical Use |
| :---: | :---: | :---: | :--- |
| $100$ | $\pm 5.0\\%$ | ~0.8 KB | Ultra-coarse screening |
| **$1{,}000$** | **$\pm 1.6\\%$** | **~8 KB** | **Mash default (bacteria)** |
| $10{,}000$ | $\pm 0.5\\%$ | ~80 KB | Fine strain / outbreak resolution |

**Key miracle**: accuracy depends **only on $s$**, completely **independent of genome size $G$**!
:::

---

# From Jaccard to ANI: Deriving the Mash Distance

Biologists think in **Average Nucleotide Identity (ANI)** — the mean base identity across aligned homologous blocks (where $\text{ANI} \ge 95\\%$ defines the microbial species boundary). How do we convert $k$-mer Jaccard $J$ into base-level mutation rate $d \approx 1 - \text{ANI}$?

:::matrix { cols="1.15fr 1fr" gap="30px" }
[[0,0]]
### Poisson Point-Mutation Model (Ondov et al., 2016)
Assume two equal-length genomes ($|A| \approx |B| = G$) differ by independent point mutations at rate $d$ per base:

1. **Probability a $k$-mer survives unmutated**:
   $$p = (1 - d)^k \approx e^{-kd}$$
2. **Expected intersection and union sizes**:
   $$|A \cap B| \approx p G, \qquad |A \cup B| = |A| + |B| - |A \cap B| \approx (2 - p) G$$
3. **Expected Jaccard similarity** ($G$ cancels!):
   $$J \approx \frac{p}{2 - p} \quad \Longleftrightarrow \quad p = \frac{2J}{1 + J}$$
4. **Invert $p = e^{-kd}$ to obtain the Mash distance $D$**:
   $$D = -\frac{1}{k}\ln\left(\frac{2J}{1 + J}\right) \approx d \approx 1 - \text{ANI}$$

[[0,1]]
### Worked Check (Matching Part 4)
Suppose `mash dist` reports **`91/1000`** shared hashes at $k = 21$ (*E. coli* vs. *Salmonella*):

1. $\hat{J} = \frac{91}{1000} = 0.091$
2. $\hat{p} = \frac{2 \times 0.091}{1 + 0.091} = 0.1668$
3. $D = -\frac{1}{21}\ln(0.1668) = \mathbf{0.0853}$
4. Estimated $\text{ANI} \approx 1 - 0.0853 = \mathbf{91.5\\%}$ (same genus, different species!).

### Species Boundary Mapping
- **Same species**: $\text{ANI} \ge 95\\% \iff D \lesssim \mathbf{0.05}$
  - At $k = 21$, $D = 0.05 \iff p = e^{-1.05} \approx 0.35 \implies J \approx \mathbf{0.21}$ (~212/1000 shared hashes!).
:::

---
# When Mash Distance Breaks Down

The clean identity $D \approx 1 - \text{ANI}$ relies on two hidden assumptions: **enough surviving $k$-mers** and **comparable genome sizes** ($|A| \approx |B|$).

![Mash Distance Breakdown Regimes](lecture_01_the_molecular_blueprint/images/mash_distance_regimes.png){width=60 center}

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Regime 1: High Divergence ($D > 0.20$, $\text{ANI} < 80\\%$)
- At $k = 21$ and $d = 0.25$, survival collapses to $p = e^{-5.25} \approx 0.005$ ($J \approx 0.0026 \implies \sim 2.6/1000$ hashes!).
- Relative sampling noise explodes; if $x = 0$, $D$ saturates at $1.0$ (or $\infty$).
- **Fix**: lower $k$ (e.g. $k = 13\text{–}15$) or increase $s$ (trading off random collisions).

[[0,1]]
### Regime 2: Unequal Genome Sizes ($|A| \ne |B|$)
- Suppose a small phage/plasmid $A$ is **100% identical** ($\text{ANI} = 100\\%$) to an integrated region of a large genome $B$.
- Exact Jaccard is only $J = |A| / |B| = \rho \ll 1$, so Mash reports $D = -\frac{1}{k}\ln\frac{2\rho}{1+\rho} > 0$ — falsely exceeding $0.05$ whenever $|A|/|B| < 0.21$ (at $k=21$)!
- **Fix**: asymmetric **containment** $|A \cap B| / |A|$ (`mash screen`) instead of symmetric Jaccard.
:::

---

# Interactive Lab — Mash Distance & Breakdown Explorer

:::matrix { cols="0.82fr 1.55fr" gap="24px" }
[[0,0]]
### Guided Experiments
1. **High-Divergence Cliff**:
   - Keep $k = 21, s = 1000$, drag **$d = 1 - \text{ANI}$** past `0.20`, and click **🎲 Resample Sketch** several times to feel the variance blow-up.
   - Now increase $s \to 5000$ or lower $k \to 15$ to watch the 95% CI shrink.
2. **Unequal-Size Illusion**:
   - Set $d = 0.000$ (**100% ANI!**) and drag **$|A|/|B|$** down to `0.16`: watch Mash $D$ cross `0.05` solely from size bias!
3. **Collision Phase Transition**:
   - Select **Human (3 Gb)** and slide $k$ down toward `13`: watch $-\log_{10}(p)$ collapse straight into the red noise band!

<button class="demo-btn" style="margin-top: 12px; font-size: 1.1rem; padding: 10px 18px;" onclick="showDemo('lecture_01_the_molecular_blueprint/demo/mash_regime_explorer.html')">⛶ Full-Screen Lab</button>

[[0,1]]
<iframe src="lecture_01_the_molecular_blueprint/demo/mash_regime_explorer.html" title="Mash Distance and Breakdown Explorer" style="width: 100%; height: 615px; border: 1px solid var(--border-color); border-radius: 10px; background: var(--secondary-bg);"></iframe>
:::

---
# Statistical Significance ($p$-value) in $4^k$ Space

When `mash dist` outputs `91/1000` with $p = 8.14 \times 10^{-10}$, how do we know those shared hashes didn't happen **purely by random chance**?

:::matrix { cols="50/50" gap="20px" }
[[0,0]]
![Statistical Significance Phase Transition](lecture_01_the_molecular_blueprint/images/pvalue_phase_transition.png){width=94 center}

### Null Hypothesis: Two Unrelated Random Genomes
Let two unrelated genomes contain $n_1, n_2$ unique $k$-mers drawn uniformly from the $4^k$ universe:
1. **Occupancy probabilities** in $4^k$ space:
   $$r_1 = 1 - \left(1 - \frac{1}{4^k}\right)^{n_1} \approx \frac{n_1}{4^k}, \qquad r_2 \approx \frac{n_2}{4^k}$$


[[0,1]]
2. **Expected random Jaccard** under independence:
   $$r = \frac{r_1 r_2}{r_1 + r_2 - r_1 r_2}$$
3. **Right-tail Binomial $p$-value** of seeing $\ge x$ shared hashes in a sketch of size $s$:
   $$p\text{-value} = P(X \ge x) = \sum_{i=x}^{s} \binom{s}{i} r^i (1 - r)^{s - i}$$
*(For finite sketches without replacement, the exact hypergeometric tail yields the same threshold.)*

- At **$k = 21$** ($4^{21} \approx 4.4 \times 10^{12}$), even **5/1000** hits between two 5 Mbp genomes gives $p \approx 10^{-18}$!
- At **$k = 11$**, random collisions yield ~540/1000 hits by chance ($p = 1.0$, useless).
:::

---

# Algorithmic Complexity & Constant-Memory Streaming

Why MinHash + streaming I/O turns a multi-hour cluster computation into a **3-second Colab notebook**:

| Algorithm / Approach | Preprocessing Time | Pairwise Compare | Storage per Genome | Memory Model |
| :--- | :---: | :---: | :---: | :---: |
| **Exact Pairwise Alignment (DP)** | — | $O(G^2)$ (~hours) | $O(G)$ (~5 MB) | $O(G^2)$ or $O(G)$ |
| **Exact All-$k$-mer Hash Set** | $O(N)$ | $O(G)$ (~seconds) | $O(G)$ (~40 MB) | Loads full $k$-mer table |
| **Bottom-$s$ MinHash (`mash`)** | **$O(N \log s)$** ($\approx O(N)$) | **$O(s)$ (~$\mu\text{s}$)** | **$O(s)$ (~1–8 KB)** | **Streaming $O(s + m)$** |

:::matrix { .cards cols="1fr 1fr" gap="30px" margin-top: 15px; }
[[0,0]]{ .green }
### Streaming I/O (`O(1)` RAM — Python Generators)
- Multi-gigabyte FASTQ files contain millions of reads and **cannot** be loaded into a Python `list`.
- `Bio.SeqIO.parse()` yields **one `SeqRecord` at a time** from a file buffer: peak RAM stays constant $O(L)$ whether the file has $10^3$ or $10^8$ reads!
- Same single-pass streaming discipline powers `fastp` and `mash sketch`.

[[0,1]]{ .orange }
### Scaling to All-vs-All Databases ($M$ Genomes)
- Comparing $M = 10{,}000$ bacterial genomes pairwise requires $\binom{M}{2} \approx 5\times 10^7$ comparisons:
  - Alignment: **$O(M^2 G^2)$** — computationally intractable.
  - MinHash: **$O(M N \log s + M^2 s)$** — sketch once, store 10 MB total, compare all 50 million pairs in minutes!
:::

---

# Summary & Bridge to Part 4 (Industrial Tools)

:::matrix { .summary cols="1fr 1fr 1fr" gap="22px" }
[[0,0]]{ .green }
### Quality & Cleaning
#### `fastp`
- Phred+33: $Q = \operatorname{ord}(c) - 33$
- $O(L)$ sliding window: $\bar{Q}\_i < Q_{\text{thresh}}$
- Paired-end shift overlap ($I < L$) strips adapters **de novo**
- Low-complexity & hash duplication stats

[[0,1]]{ .blue }
### $k$-mers & Sketching
#### `Bio.SeqIO` · `mash sketch`
- $\operatorname{canon}(x) = \min(x, \overline{x}^{\text{rev}})$ in 2-bit ops
- Goldilocks $k$: $q \approx G/4^k$ vs. $p = e^{-kd}$
- Bloom filter (`-m 2`) kills singleton errors
- Bottom-$s$ max-heap + $\hat{G} \approx s H_{\max} / h_{(s)}$
- $O(1)$-RAM generator streaming

[[0,2]]{ .orange }
### Distance & Statistics
#### `mash dist`
- $O(s)$ merge: $\hat{J} = |W \cap S(A) \cap S(B)| / s$
- $\operatorname{SE}(\hat{J}) = \sqrt{J(1-J)/s}$ (indep. of $G$)
- $D = -\frac{1}{k}\ln\frac{2J}{1+J} \approx 1 - \text{ANI}$ ($< 0.05$ species)
- Beware $D > 0.2$ & unequal $|A| \ne |B|$
- Binomial $p$-value in $4^k$ space
:::

> **Next up (Part 4 — Industrial Tools):** we put every algorithm above straight into practice with `fasterq-dump`, `fastp`, `Bio.SeqIO`, and `mash` on live NCBI SRA datasets.
