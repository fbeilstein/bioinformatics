:::titlepage
[[top]]
Part 4

[[title]]
Industrial Tools

[[subcaption]]
Bio.Entrez · UniProt REST · Bio.Align.PairwiseAligner · Bio.Phylo
:::

---

# Curriculum

- **Bio.Align.PairwiseAligner** (Biopython C-accelerated module) — global & local alignment with affine gaps.
- **Bio.Entrez** + **Bio.SeqIO** — downloading curated sequences from NCBI RefSeq.
- **UniProt REST API** — downloading protein sequences with `requests`.
- **Bio.Align.substitution_matrices** — BLOSUM62 for protein scoring.
- **Bio.Phylo** — distance matrices, Neighbor-Joining trees, outgroup rooting.
- **scipy.stats.linregress** — testing the molecular clock.

---

# Roadmap

Two pipelines in the practice notebook — a **tree** (practice) and a **clock** (homework):

:::static-diagram { width=95% center }
{
  "nodes": [
    { "id": "ncbi",   "label": "NCBI RefSeq\n(mRNA)",          "x": 50,   "y": 80,  "width": 140, "height": 60 },
    { "id": "entrez", "label": "Bio.Entrez\n+ SeqIO",          "x": 230,  "y": 80,  "width": 140, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "local",  "label": "PairwiseAligner\n(local, DNA)", "x": 420,  "y": 80,  "width": 170, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "dmat",   "label": "Distance\nMatrix",             "x": 640,  "y": 80,  "width": 130, "height": 60 },
    { "id": "nj",     "label": "Bio.Phylo\n(NJ + root)",       "x": 820,  "y": 80,  "width": 150, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "tree",   "label": "Rooted\nTree",                 "x": 1020, "y": 80,  "width": 120, "height": 60 },

    { "id": "uni",    "label": "UniProt\n(protein)",           "x": 50,   "y": 260, "width": 140, "height": 60 },
    { "id": "req",    "label": "requests\n(REST JSON)",        "x": 230,  "y": 260, "width": 140, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "glob",   "label": "PairwiseAligner\n(global, BLOSUM62)", "x": 420, "y": 260, "width": 170, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "pois",   "label": "Identity →\nPoisson d",        "x": 640,  "y": 260, "width": 130, "height": 60 },
    { "id": "lin",    "label": "linregress\n(d vs MYA)",       "x": 820,  "y": 260, "width": 150, "height": 60, "style": "fill: #FF9800; color: #fff;" },
    { "id": "clock",  "label": "Molecular\nClock",             "x": 1020, "y": 260, "width": 120, "height": 60 }
  ],
  "links": [
    { "from": "ncbi",   "to": "entrez" },
    { "from": "entrez", "to": "local" },
    { "from": "local",  "to": "dmat" },
    { "from": "dmat",   "to": "nj" },
    { "from": "nj",     "to": "tree" },
    { "from": "uni",    "to": "req" },
    { "from": "req",    "to": "glob" },
    { "from": "glob",   "to": "pois" },
    { "from": "pois",   "to": "lin" },
    { "from": "lin",    "to": "clock" }
  ]
}
:::

---

# Biopython — The Toolbox

**Biopython** is a free, community-maintained collection of Python modules for computational biology. We use five of its sub-modules today (plus `requests`, `pandas`, and `scipy.stats`):

:::matrix { .cards cols="1fr 1fr 1fr 1fr 1fr" gap="20px" margin-top: 20px; }
[[0,0]]{ .purple }
### Bio.Entrez
Talks to NCBI databases over the internet (E-utilities).

[[0,1]]{ .purple }
### Bio.SeqIO
Parses FASTA / GenBank / FASTQ into `SeqRecord` objects (Lecture 1).

[[0,2]]{ .green }
### Bio.Align
`PairwiseAligner` — dynamic programming written in C.

[[0,3]]{ .green }
### substitution_matrices
BLOSUM, PAM and 30+ other scoring tables.

[[0,4]]{ .blue }
### Bio.Phylo
Builds, roots and draws phylogenetic trees.
:::

```python
!pip install -q biopython requests      # Colab does not ship Biopython by default
from Bio import Entrez, SeqIO, Align, Phylo
from Bio.Align import substitution_matrices
```

---

# Tool 0a — Bio.Entrez

## Where Do the Sequences Come From?

**NCBI** hosts two nucleotide collections. We use **RefSeq**:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### GenBank
- Raw **archive** of everything anyone submitted
- Redundant (same gene many times)
- Annotation quality varies
- Example: `U03460` (lamprey globin, 1994)

[[0,1]]

### RefSeq
- **Curated**, non-redundant reference set
- One "accepted" record per transcript
- Consistent annotation
- Recognised by an underscore prefix
:::

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Prefix | Meaning | Example from notebook |
| :--- | :--- | :--- |
| `NM_` | Curated mRNA (experimentally supported) | `NM_000558.5` — human HBA1 |
| `XM_` | **Predicted** mRNA (computational gene model) | `XM_003585575.4` — cow HBA1 |
| `NP_` / `XP_` | Corresponding protein (curated / predicted) | — |
| `NC_` | Complete chromosome / genome | — |
:::

> **Versions matter**: `NM_000558.5` pins version 5. An unversioned accession (e.g., `NM_131144`) silently returns the *latest* version — your results may change next year.

---

# Bio.Entrez — Downloading a Record

```python
from Bio import Entrez, SeqIO
import time

Entrez.email = "student@example.com"      # NCBI requires you to identify yourself

handle = Entrez.efetch(db="nucleotide",    # which NCBI database
                       id="NM_000558.5",   # accession
                       rettype="fasta",    # return format
                       retmode="text")
record = SeqIO.read(handle, "fasta")       # read() = exactly ONE record
print(record.id, len(record.seq))          # NM_000558.5  577

time.sleep(0.5)                            # be polite between requests
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### `SeqIO.read` vs `SeqIO.parse`
- `read()` — expects **exactly one** record, raises an error otherwise
- `parse()` — iterator over **many** records (Lecture 1)

[[0,1]]

### Rate Limits
- **3 requests/s** without an API key
- **10 requests/s** with `Entrez.api_key = "..."`
- Exceed it → HTTP 429 / temporary IP ban
- Hence `time.sleep(0.5)` in the loop
:::

---

# Why mRNA and Not Genomic DNA?

A gene in the genome contains **introns** that can be 10–100× longer than the coding part. The mature mRNA has them spliced out — but still carries **untranslated regions (UTRs)**:

```
genomic DNA:  ──[exon1]────────intron────────[exon2]──────intron──────[exon3]──
                                     │ splicing
                                     ▼
mRNA:         5'UTR ═[  CDS (protein-coding)  ]═ 3'UTR ── AAAAA
              ▲ variable                         ▲ variable length,
                length                             poorly conserved
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Conserved
The **CDS** (coding sequence) is under strong selective pressure — mutations change the protein.

[[0,1]]

### Not conserved
**UTRs** differ in length and content between species → they *will* hurt an end-to-end (global) alignment.
:::

> This single biological fact is why the practice notebook uses **local** alignment for mRNA.

---

# Tool 0b — UniProt REST API

For the homework we skip DNA entirely and download **protein** sequences from **UniProt** — no special library, just HTTP + JSON:

:::matrix { cols="1.2fr 1fr" gap="30px" }
[[0,0]]

```python
import requests

def fetch_uniprot_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return None
    data = r.json()
    return data["sequence"]["value"]

fetch_uniprot_sequence("P69905")
# 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF...'
```

[[0,1]]

### JSON Structure (abridged)

```
{
  "primaryAccession": "P69905",
  "entryType": "UniProtKB reviewed (Swiss-Prot)",
  "organism": {"scientificName": "Homo sapiens"},
  "sequence": {
      "value":  "MVLSPADKTN...",
      "length": 142
  }
}
```
:::

:::matrix { .cards cols="1fr 1fr" gap="30px" margin-top: 10px; }
[[0,0]]{ .green }
### Swiss-Prot (reviewed)
Manually curated by experts. Typically 6-character IDs: `P69905`, `P99999`.

[[0,1]]{ .orange }
### TrEMBL (unreviewed)
Automatic translation of genomes. Often `A0A…` IDs: `A0A2J8R2J5`, `F6UZ60`. Check `entryType`, not the ID shape.
:::

---

# Why Proteins for the Molecular Clock?

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### DNA saturates quickly
- Only **4 letters** — two random sequences already match at **~25%** of positions
- **Synonymous** (silent) 3rd-codon changes accumulate fast and are invisible to selection
- After ~100 MYA many sites have mutated **more than once**

[[0,1]]

### Protein keeps the signal longer
- **20 letters** — random position-by-position identity only ~**6%**
- Only **non-synonymous** changes are recorded
- Selection slows the rate → clock "ticks" slower and stays readable across 400+ MYA
:::

### Three proteins, three speeds

| Protein | Gene | Function | Expected rate |
| :--- | :--- | :--- | :--- |
| Cytochrome c | `CYCS` | Electron transport chain — every residue matters | **slow** (highly conserved) |
| Haemoglobin α | `HBA1` | O₂ transport — surface residues can drift | **medium** |
| Fibrinogen α | `FGA` | Blood clotting — its *fibrinopeptide* is cleaved off and discarded | **fast** |

---

# Tool 1 — Bio.Align.PairwiseAligner

## Dynamic Programming in C

`PairwiseAligner` implements the DP algorithms from Part 3. The Python object only **stores the settings**; the matrix filling runs in a compiled **C extension** (`Bio.Align._pairwisealigner`).

:::matrix { .accent cols="1fr 1fr" gap="40px" margin-top: 20px; }
[[0,0]]
### One object, many alignments
- Configure **once**, reuse for thousands of pairs
- Same object does DNA, RNA and protein
- Works with `str` or `Seq` objects
- Picks the right algorithm **automatically**

[[0,1]]
### `aligner.algorithm`

| Settings | Algorithm chosen |
| :--- | :--- |
| global, linear gaps | Needleman-Wunsch |
| local, linear gaps | Smith-Waterman |
| global, affine gaps | Gotoh global |
| local, affine gaps | Gotoh local |
:::

> **Warning — the defaults are not biology!** A fresh `PairwiseAligner()` uses match = 1, mismatch = 0, gap = −1 (pure edit-distance style). Always set your scores explicitly.

---

# PairwiseAligner — Configuration

The two aligners used in the practice notebook:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Practice — mRNA (DNA)

```python
from Bio import Align

aligner = Align.PairwiseAligner()
aligner.mode = "local"         # Smith-Waterman
aligner.match_score = 5        # conserved base
aligner.mismatch_score = -4    # point mutation
aligner.open_gap_score = -10   # new indel event
aligner.extend_gap_score = -1  # longer indel

print(aligner.algorithm)
# Gotoh local alignment algorithm
```

[[0,1]]

### Homework — Protein

```python
from Bio import Align
from Bio.Align import substitution_matrices

aligner = Align.PairwiseAligner()
aligner.mode = "global"        # end-to-end
aligner.substitution_matrix = \
    substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -1

print(aligner.algorithm)
# Gotoh global alignment algorithm
```
:::

> Scores are **log-odds**: positive = "more likely in real homologs than by chance", negative = "less likely". All parameters can also be passed in the constructor: `Align.PairwiseAligner(mode="local", match_score=5, ...)`.

---

# Affine Gap Penalty — How Biopython Counts

Biopython's convention: a gap of length $L$ costs

$$ \text{gap}(L) = \text{open} + (L - 1)\cdot\text{extend} $$

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### open = −10, extend = −1

| Gap length $L$ | Affine cost | Linear cost ($-10L$) |
| :--- | :--- | :--- |
| 1 | −10 | −10 |
| 3 | −12 | −30 |
| 6 | −15 | −60 |
| 30 | −39 | −300 |

One long gap is **much** cheaper than many short ones.

[[0,1]]

### Why? Replication slippage

A single slippage event during replication can insert or delete **several bases at once**.

Example — cystic fibrosis **ΔF508** (CFTR c.1521_1523delCTT): one event deletes 3 bp and removes Phe508.

```
I    I    F    G    V
ATC  ATC  TTT  GGT  GTT
ATC  AT-  --T  GGT  GTT
       └──┘ one gap of 3, not three gaps
```

Open ≫ extend says: "*starting* an indel is rare; *making it longer* is cheap."
:::

> **Careful when comparing tools**: some textbooks and programs use $\text{open} + L\cdot\text{extend}$. Same idea, numbers shifted by one `extend`.

---

# Global vs. Local Mode

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### `mode = "global"` (Needleman-Wunsch / Gotoh)

Aligns **every residue** of both sequences, end to end.

```
5'UTR══════CDS══════3'UTR───────
───5'UTR═══CDS═══3'UTR══════════
   ▲ forced to align junk UTRs
```

✓ Orthologous **proteins** of similar length
✗ mRNAs with different UTRs

[[0,1]]

### `mode = "local"` (Smith-Waterman / Gotoh)

Finds the **best-scoring sub-region**; flanks are ignored.

```
5'UTR══[═══CDS═══]══3'UTR───────
───5'UTR[═══CDS═══]3'UTR════════
         ▲ only this is reported
```

✓ mRNAs, domains, reads vs. genome
✗ Can report a short "perfect" fragment
:::

### Middle ground: semi-global (free end gaps)

```python
aligner = Align.PairwiseAligner(mode="global", match_score=5, mismatch_score=-4,
                                open_gap_score=-10, extend_gap_score=-1,
                                end_gap_score=0)       # gaps at the ends are free
print(aligner.align("TTTTTTACGTACGTGGGGGG", "ACGTACGT")[0])
# target   0 TTTTTTACGTACGTGGGGGG 20
#          0 ------||||||||------ 20
# query    0 ------ACGTACGT------  8
```

---

# Substitution Matrices — BLOSUM62

For proteins, a single "mismatch" score is too crude: **Leu → Ile** (both hydrophobic) is harmless, **Trp → Gly** destroys the fold.

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Selected BLOSUM62 scores

| Pair | Score | Interpretation |
| :--- | :--- | :--- |
| W–W | **+11** | Trp is rare → matching it is strong evidence |
| C–C | +9 | Cys (disulfide bonds) |
| A–A | +4 | Ala is common → weak evidence |
| L–I | +2 | similar hydrophobics |
| D–E | +2 | both acidic |
| L–K | −2 | hydrophobic ↔ charged |
| W–G | −2 | bulky ↔ tiny |

[[0,1]]

### Where do the numbers come from?

Counted from **BLOCKS** of real, conserved protein alignments, clustered at **62% identity**:

$$ S(a,b) = 2\log_2 \frac{p_{ab}}{q_a\,q_b} $$

- $p_{ab}$ — observed frequency of pair $a$–$b$
- $q_a q_b$ — expected by chance

```python
from Bio.Align import substitution_matrices
B = substitution_matrices.load("BLOSUM62")
B["W"]["W"]                # 11.0
substitution_matrices.load()  # list all 30+
```

Lower number (BLOSUM45) → distant homologs; higher (BLOSUM80) → close ones.
:::

---

# Running an Alignment

```python
alignments = aligner.align(seq1, seq2)   # PairwiseAlignments object (lazy!)
best = alignments[0]                     # first optimal alignment
best.score                               # its score

aligner.score(seq1, seq2)                # score ONLY — no traceback, much faster
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Lazy evaluation

`aligner.align()` fills the DP matrix but builds alignments **on demand** during traceback. There may be many **co-optimal** alignments with the same score — ties in the matrix multiply.

Two random 2000-bp DNA sequences, local mode:

```python
len(aligner.align(s, t))
# OverflowError: number of optimal alignments
# is larger than 9223372036854775807
```

**Never** do `list(aligner.align(...))`.

[[0,1]]

### When to use what

| Need | Call | Time (2 kb × 2 kb) |
| :--- | :--- | :--- |
| Just the score | `aligner.score(a, b)` | ~0.03 s |
| One alignment | `aligner.align(a, b)[0]` | ~0.2 s |
| All alignments | iterate (rarely needed) | can be ∞ |

Both scale as $O(N \cdot M)$ — doubling sequence length → **4×** the time.
:::

---

# The Alignment Object — Anatomy

```python
best = aligner.align("HEAGAWGHEE", "PAWHEAE")[0]     # global, BLOSUM62, -10/-1
print(best)
```

```
target            0 HEAGAWGHEE 10
                  0 ---.||...| 10
query             0 ---PAWHEAE  7
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Reading the printout
- **target** = first argument, **query** = second
- Numbers = sequence coordinates (0-based, start/end)
- `|` identity  · `.` mismatch  · `-` gap

[[0,1]]

### Useful attributes

```python
best.score        # 3.0
best.shape        # (2, 10)  rows × columns
best[0]           # 'HEAGAWGHEE'  gapped target
best[1]           # '---PAWHEAE'  gapped query
best.coordinates  # aligned-block boundaries

c = best.counts()
c.identities      # 3
c.mismatches      # 4
c.gaps            # 3
```
:::

---

# Worked Example — Where Does the Score 3 Come From?

```
HEAGAWGHEE
---PAWHEAE
```

| Columns | Pairs | BLOSUM62 | Subtotal |
| :--- | :--- | :--- | :--- |
| 1–3 | `HEA` / `---` | gap, $L = 3$: $-10 + 2 \cdot (-1)$ | **−12** |
| 4 | G / P | −2 | |
| 5 | A / A | +4 | |
| 6 | W / W | +11 | |
| 7 | G / H | −2 | |
| 8 | H / E | 0 | |
| 9 | E / A | −1 | |
| 10 | E / E | +5 | **+15** |
| | | **Total** | $-12 + 15 = \mathbf{3}$ |

> The single W–W match (+11) carries most of the alignment — rare residues are strong evidence of homology.

---

# ⚠ Pitfall — Don't Parse the Pretty-Print

In modern Biopython, `format()` / `str()` **wraps long alignments every 60 columns** and adds labels and coordinates. Real output for human vs. mouse α-globin (142 aa):

```
target            0 MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
                  0 ||||..||.|.||||||.|.|..|||||||||||.|||||||||||||.|||||||||||
query             0 MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHG

target           60 KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
...
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### ✗ Fragile

```python
lines = alignment.format().split("\n")
a, b = lines[0], lines[2]
# → only the FIRST 60 columns,
#   plus "target            0 " labels!
```

[[0,1]]

### ✓ Robust

```python
a, b = alignment[0], alignment[1]
# full gapped strings, no labels

# or let Biopython count for you:
c = alignment.counts()
c.identities, c.mismatches, c.gaps
```
:::

> Counting `"|"` in `str(alignment)` (as in the practice cell) *does* work across wrapped blocks — the coordinate lines contain only digits — but `alignment[0]` / `alignment[1]` is the clean way.

---

# Defining "Percent Identity"

There is **no single definition** — the notebook uses two. For `HEAGAWGHEE` / `---PAWHEAE`:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Practice (cell "local aligner")

$$ \text{identity} = \frac{\text{matches}}{\text{alignment length}} = \frac{3}{10} = 30\% $$

Gaps count **as differences**.

```python
matches = str(aln).count("|")
length  = aln.shape[1]
```

[[0,1]]

### Homework Task 1 (`calculate_identity`)

$$ \text{identity} = \frac{\text{matches}}{\text{comparable columns}} = \frac{3}{7} \approx 42.9\% $$

Gap columns (`"-"`) are **ignored**:
- Extract full gapped strings via `aln[0]` and `aln[1]` and iterate column-by-column with `zip()`.
- Count only columns where **neither** character is `"-"`.
- Guard the edge case: if `comparable == 0`, return `np.nan`.
:::

> **Always state your definition.** And beware local mode: identity is measured only over the *reported* region — a 20-residue perfect local hit gives "100%" even if the rest of the proteins differ.

---

# From Identity to Evolutionary Distance

The raw difference $p = 1 - \text{identity}$ (the **p-distance**) *underestimates* how many substitutions really happened — a site can mutate **twice** (A→G→A looks like no change; A→G→T looks like one).

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Poisson correction (Homework Task 2)

If substitutions hit sites randomly at a constant rate, the number of hits per site is Poisson-distributed. The fraction of sites with **zero** hits is $e^{-d}$, so

$$ 1 - p = e^{-d} \quad\Longrightarrow\quad d = -\ln(1 - p) $$

- $d$ = expected substitutions **per site**.
- Use natural logarithm `np.log()` (not `log10`).
- **Edge case**: when $p \ge 1$ (zero shared residues), $\ln(0)$ is undefined — return `np.inf`.

[[0,1]]

### Correction grows with divergence

| Observed $p$ | Corrected $d$ | Hidden hits |
| :--- | :--- | :--- |
| 0.10 | 0.105 | +5% |
| 0.30 | 0.357 | +19% |
| 0.50 | 0.693 | +39% |
| 0.70 | 1.204 | +72% |
| 0.90 | 2.303 | +156% |

Close species: $d \approx p$. Distant species: $p$ **saturates**, $d$ keeps growing.

The practice cell's `get_distance` returns plain $100 - \text{identity}$ — an uncorrected p-distance in %.
:::

---

# Tool 2 — Bio.Phylo: Distance Matrix

`Bio.Phylo.TreeConstruction` wants the **lower triangle** of the (symmetric) distance matrix, diagonal included:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Full matrix → lower triangle

```
         Hs    Mm    Gg    Pm
Hs     [ 0 ]
Mm     [ d   , 0 ]
Gg     [ d   , d   , 0 ]
Pm     [ d   , d   , d   , 0 ]
```

Row $i$ has $i+1$ entries; the upper half is implied by symmetry $d(i,j) = d(j,i)$.

[[0,1]]

### Building it

```python
from Bio.Phylo.TreeConstruction import (
    DistanceMatrix, DistanceTreeConstructor)

seqs  = sequences["HBA1"]
names = list(seqs.keys())
lower = []
for i in range(len(names)):
    row = []
    for j in range(i + 1):
        row.append(0.0 if i == j else
                   get_distance(seqs[names[i]],
                                seqs[names[j]]))
    lower.append(row)

dm = DistanceMatrix(names=names, matrix=lower)
```
:::

$n$ species → $n(n-1)/2$ alignments. 8 species = 28 alignments; 100 species = 4,950.

---

# Bio.Phylo — Neighbor-Joining vs. UPGMA

```python
constructor = DistanceTreeConstructor()
nj_tree    = constructor.nj(dm)       # Neighbor-Joining
upgma_tree = constructor.upgma(dm)    # UPGMA (for comparison)
```

:::matrix { .cards cols="1fr 1fr" gap="30px" margin-top: 20px; }
[[0,0]]{ .blue }
### Neighbor-Joining (used)
- Allows **different rates** on different branches
- Branch lengths = actual change in each lineage
- Produces an **unrooted** tree
- Fast: $O(n^3)$

[[0,1]]{ .orange }
### UPGMA
- Assumes a **strict molecular clock**
- All leaves equidistant from the root
- Rooted by construction
- Wrong topology if rates differ between lineages
:::

> NJ makes **no clock assumption** — which is exactly why we can later *test* the clock independently instead of building it into the tree.

---

# Bio.Phylo — Outgroup Rooting & Drawing

NJ knows **how** species connect, but not **where the ancestor is**. We place the root on the branch of a species known to have diverged first — the **outgroup**: sea lamprey (*Petromyzon marinus*, jawless fish, ~600 MYA).

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]

```python
nj_tree.root_with_outgroup({"name": "Petromyzon marinus"})

fig, ax = plt.subplots(figsize=(10, 6))
Phylo.draw(nj_tree, axes=ax, do_show=False,
           label_func=lambda c: c.name
                      if c.is_terminal() else "")
plt.show()
```

- `root_with_outgroup` takes a clade *selector* — here a dict matching the leaf name
- `label_func` hides NJ's internal node names (`Inner1`, `Inner2`, …)
- `Phylo.draw_ascii(tree)` — quick text preview

[[0,1]]

### Toy example (`draw_ascii`)

```
                 ______ A
         _______|
  ______|       |_________ B
 |      |
_|      |___________ C
 |
 | D   ← outgroup
```

Horizontal length = **branch length** (evolutionary distance), not time.
:::

---

# Testing the Molecular Clock — `pandas` & `linregress`

If substitutions accumulate at a constant rate, distance is **linear in divergence time**:

$$ d(t) \approx r \cdot t $$

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### DataFrame & Linear Fit
```python
import pandas as pd, numpy as np
from scipy.stats import linregress

# Build table from a list of dicts (Homework Task 3)
df  = pd.DataFrame(records)
sub = df[df["Gene"] == "CYCS"]

fit = linregress(sub["MYA"], sub["Distance"])
# fit.slope        → r (subs/site/MYA)
# fit.intercept    → should be ≈ 0
# fit.rvalue ** 2  → R² (how clock-like)

# Overlay dashed regression line (Tasks 4–5)
x_line = np.linspace(sub["MYA"].min(), sub["MYA"].max(), 100)
plt.plot(x_line, fit.slope * x_line + fit.intercept, "--")
```

[[0,1]]

### Reading the plots
- **Slope** — speed of the clock. Steeper = faster-evolving protein
- **$R^2 \to 1$** — points lie on a line → clock-like
- **Curvature in identity plot** — saturation (multiple hits)
- **Scatter** — rate varies between lineages, or uncertain MYA values
- Check `len(sub) >= 3` before fitting a line; use `plt.xlim()` / `plt.ylim()` when fixed axis ranges are requested
:::

> Compare the **identity** plot (Task 4) with the **corrected distance** plot (Task 5): where does the Poisson correction change the picture most — and why?

---

# Putting It All Together

```python
from Bio import Entrez, SeqIO, Align, Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import time

# 1. Fetch curated RefSeq mRNA sequences
Entrez.email = "student@example.com"
seqs = {}
for species, acc in accession_db["HBA1"].items():
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
    seqs[species] = SeqIO.read(handle, "fasta").seq
    time.sleep(0.5)

# 2. Configure local Gotoh aligner (mRNA has variable-length UTRs)
aligner = Align.PairwiseAligner(mode="local", match_score=5, mismatch_score=-4,
                                open_gap_score=-10, extend_gap_score=-1)

# 3. Build lower-triangular distance matrix → Neighbor-Joining tree → root on outgroup
dm   = DistanceMatrix(names=names, matrix=lower)
tree = DistanceTreeConstructor().nj(dm)
tree.root_with_outgroup({"name": "Petromyzon marinus"})
Phylo.draw(tree, label_func=lambda c: c.name if c.is_terminal() else "")

# ---- HOMEWORK: THE MOLECULAR CLOCK ----
# Switch to UniProt proteins + global BLOSUM62 aligner:
#   Task 1: calculate_identity() (ignore gaps)   Task 2: poisson_distance()
#   Task 3: build DataFrame of all 3 genes       Tasks 4–5: scatter + linregress plots
```

---

# Summary

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### Entrez & UniProt
#### Data Retrieval
- RefSeq: curated, `NM_` / `XM_`
- Pin versions (`.5`)
- Respect rate limits
- UniProt REST → JSON
- Swiss-Prot vs TrEMBL
- **Output**: sequences

[[0,1]]{ .green }
### PairwiseAligner
#### Alignment
- C-accelerated DP
- `global` / `local` / free end gaps
- Affine gaps: open + (L−1)·extend
- BLOSUM62 for proteins
- `aln[0]`, `aln[1]`, `.counts()`
- **Output**: identity → distance

[[0,2]]{ .blue }
### Bio.Phylo + scipy
#### Evolution
- Lower-triangular `DistanceMatrix`
- NJ: no clock assumption
- Outgroup rooting
- Poisson: $d = -\ln(1-p)$
- `linregress`: slope & $R^2$
- **Output**: tree & clock rate
:::

---

# Next Steps

### In the Practice Session You Will:

1. **Download** α-globin (HBA1) mRNA for 8 vertebrates from NCBI RefSeq with `Bio.Entrez`
2. **Align** them with a **local** Gotoh aligner (match +5, mismatch −4, gap −10/−1)
3. **Build** a distance matrix and a **Neighbor-Joining** tree, rooted on the sea lamprey

### Homework — The Molecular Clock

1. Fetch CYCS, HBA1 and FGA **proteins** from UniProt for up to 19 species
2. Implement `calculate_identity()` (global alignment, BLOSUM62, gaps ignored)
3. Implement `poisson_distance()`
4. Build the data table, plot **identity vs. MYA** and **distance vs. MYA**, fit regression lines

> **Questions to answer**: Which protein ticks fastest? Which fits a linear clock best ($R^2$)? Where does the Poisson correction matter most? What makes a real molecular-clock analysis harder than ours?
