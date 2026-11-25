:::titlepage
[[top]]
Part 4

[[title]]
Industrial Tools

[[subcaption]]
SRA Toolkit · fastp · Bio.SeqIO · mash
:::

---
# Plan

- **SRA Toolkit** (`fasterq-dump`) for downloading raw reads from NCBI.
- **fastp** for adapter trimming and quality control.
- **Bio.SeqIO** for stream parsing of sequence files.
- **mash** for sketch-based distance calculations (`mash sketch`, `mash dist`).

# Roadmap

Four tools form a complete pipeline — from public database to species identification:

:::static-diagram { width=90% center }
{
  "nodes": [
    { "id": "sra",   "label": "NCBI SRA\nDatabase",         "x": 50,  "y": 200, "width": 140, "height": 60 },
    { "id": "fqdump","label": "fasterq-dump\n(Download)",    "x": 230, "y": 200, "width": 150, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "raw",   "label": "Raw FASTQ\nReads",            "x": 420, "y": 200, "width": 130, "height": 60 },
    { "id": "fastp", "label": "fastp\n(QC & Trim)",          "x": 590, "y": 200, "width": 130, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "clean", "label": "Clean FASTQ\nReads",           "x": 760, "y": 200, "width": 130, "height": 60 },
    { "id": "seqio", "label": "Bio.SeqIO\n(Parse & Inspect)", "x": 760, "y": 80,  "width": 160, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "mash",  "label": "mash\n(Sketch & Dist)",       "x": 940, "y": 200, "width": 140, "height": 60, "style": "fill: #FF9800; color: #fff;" },
    { "id": "dist",  "label": "Species\nIdentification",     "x": 1120, "y": 200, "width": 140, "height": 60 }
  ],
  "links": [
    { "from": "sra",    "to": "fqdump" },
    { "from": "fqdump", "to": "raw" },
    { "from": "raw",    "to": "fastp" },
    { "from": "fastp",  "to": "clean" },
    { "from": "clean",  "to": "seqio" },
    { "from": "clean",  "to": "mash" },
    { "from": "mash",   "to": "dist" }
  ]
}
:::

---
# Tool 0 — SRA Toolkit

## Getting Data: `fasterq-dump`

Before any analysis, we need raw sequencing reads. The **Sequence Read Archive (SRA)** at NCBI is the world's largest public repository of high-throughput sequencing data.

Every dataset has a unique accession number (e.g., `SRR40792821`). The **SRA Toolkit** provides `fasterq-dump` to download and convert these into FASTQ files.

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
```bash
# Download paired-end reads by SRR accession
fasterq-dump SRR40792821 --split-files --progress
```
> **Critical**: Read $N$ in the `_1` file and read $N$ in the `_2` file come from the **same physical DNA fragment**. This pairing must be preserved — never shuffle or sort one file without the other.

[[0,1]]
### What `--split-files` Produces

For paired-end data, you get **two synchronized files**:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### `SRR40792821_1.fastq`
**Forward reads** (Read 1)

The sequencer reads from the **5' → 3'** end of the fragment.

```
@SRR40792821.1
ACGTACGTACGT...
+
IIIIIIIIIIII...
```

[[0,1]]

### `SRR40792821_2.fastq`
**Reverse reads** (Read 2)

The sequencer reads the **opposite strand** of the same fragment.

```
@SRR40792821.1
TGCATGCATGCA...
+
IIIIIIIIIIII...
```
:::

:::

---

# SRA Toolkit — Spots vs. Reads

After downloading, `fasterq-dump` prints a summary:

```
spots read      : 500,000
reads read      : 1,000,000
reads written   : 1,000,000
```

:::matrix { .cards .purple cols="1fr 1fr 1fr" gap="40px" margin-top: 20px; }
[[0,0]]
### Spot
A single physical DNA fragment on the flow cell. In NCBI's terminology, one **spot** = one cluster.

[[0,1]]
### Reads
For paired-end sequencing: **reads = 2 × spots**. Each spot produces one forward read and one reverse read.

[[0,2]]
### Reads Written
Should equal reads read. If it doesn't, some reads were lost — **re-download**.
:::

---
# Tool 1 --- fastp

## What is it?

**fastp** is an ultra-fast, all-in-one FASTQ preprocessor written in C++.

It handles quality control, adapter trimming, and filtering in a **single pass** --- replacing an entire chain of older tools (FastQC → Trimmomatic → cutadapt).

:::matrix { .accent cols="1fr 1fr" gap="40px" margin-top: 30px; }
[[0,0]]
### Key Features
- Automatic adapter detection (no adapter file needed)
- Per-read quality filtering & trimming
- Poly-G / Poly-X tail trimming (NovaSeq fix)
- Overrepresented sequence analysis
- HTML report with interactive plots

[[0,1]]
### Performance
- Multi-threaded (typically 2–4 threads is optimal)
- Processes ~300 MB/s on a single core
- Written in C++ with SIMD optimizations
- **2–5× faster** than Trimmomatic
:::

---
# Why Do We Need Adapter Trimming?

During library preparation, **sequencing adapters** are ligated to both ends of the DNA fragment. If the fragment is shorter than the read length, the sequencer reads into the adapter on the opposite end:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Short Insert Problem

```
5'─[Adapter A]──ACGTACGT──[Adapter B]─3'
    ◄── Read 1 ──────────────────────►
                          ▲
                   Adapter contamination!
```

When insert size < read length, adapter sequence contaminates the 3' end of the read.

[[0,1]]

### After Trimming

```
5'──ACGTACGT──3'
    ◄── Read ►

✓ Pure genomic sequence
✓ Accurate k-mer decomposition
✓ Reliable quality scores
```

fastp detects adapters automatically using overlap analysis of paired-end reads.
:::

---
# fastp — Command Line Usage

### Basic Paired-End Trimming

```bash
fastp \
  -i reads_R1.fastq.gz \     # input forward reads
  -I reads_R2.fastq.gz \     # input reverse reads
  -o clean_R1.fastq.gz \     # output forward reads
  -O clean_R2.fastq.gz \     # output reverse reads
  --thread 4 \               # parallel threads
  --html fastp_report.html   # QC report
```

### What Happens Under the Hood

1. **Adapter detection** — overlaps paired reads to find adapter sequence
2. **Quality trimming** — slides a window across each read, trims when mean $Q < 20$
3. **Length filter** — discards reads shorter than 15 bp (default)
4. **Complexity filter** — removes low-complexity reads (e.g., poly-A)
5. **Report generation** — produces an interactive HTML summary

---
# fastp — Quality Score Filtering

Recall that Phred quality scores encode the probability of a base call error:
$$ 
  Q = -10 \log_{10} P_{\text{error}} 
$$

fastp uses a **sliding window** approach to trim low-quality tails:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Sliding Window Trim

For a window of size $w$ starting at position $i$:

$$ \bar{Q}_i = \frac{1}{w} \sum_{j=i}^{i+w-1} Q_j $$

If $\bar{Q}_i < Q_{\text{threshold}}$, trim from position $i$ onward.

**Default**: $w = 4$, $Q_{\text{threshold}} = 20$ (i.e., $P_{\text{error}} \leq 1\%$).

[[0,1]]

### Common Quality Thresholds

| Phred $Q$ | Error Rate | Accuracy |
|-----------|-----------|----------|
| 10 | 1 in 10 | 90% |
| 20 | 1 in 100 | 99% |
| 30 | 1 in 1,000 | 99.9% |
| 40 | 1 in 10,000 | 99.99% |

A threshold of $Q \geq 20$ keeps only bases with $\leq 1\%$ error probability — the standard for downstream analysis.
:::

---

# fastp — The HTML Report

fastp generates a self-contained HTML report with interactive charts:

:::matrix { .cards cols="1fr 1fr 1fr" gap="30px" margin-top: 20px; }
[[0,0]]{ .green }
### Before Filtering
- Total reads
- Total bases
- Q20 & Q30 percentages
- GC content distribution
- Per-base quality profile

[[0,1]]{ .orange }
### After Filtering
- Reads passing filter
- Adapter-trimmed reads (%)
- Quality-trimmed reads (%)
- Duplication rate estimate
- Insert size distribution

[[0,2]]{ .blue }
### Key Diagnostics
- Overrepresented sequences
- K-mer frequency profiles
- Adapter content curves
- Read length distribution
- Quality score heatmaps
:::

> **In practice**: always inspect the fastp report before proceeding. A sudden drop in Q30% or unexpected adapter content signals a problem in library preparation.

---

# Tool 2 — Bio.SeqIO

## Python's Swiss Army Knife for Sequence Files

**Bio.SeqIO** is the I/O module of the **Biopython** library. It provides a **unified interface** to read and write over 20 sequence file formats.

:::matrix { .accent cols="1fr 1fr" gap="40px" margin-top: 30px; }
[[0,0]]
### Supported Formats
- **FASTA** — reference genomes
- **FASTQ** — sequencing reads + quality
- **GenBank** — annotated sequences
- **EMBL** — European sequence format
- **ABI** — Sanger chromatograms

[[0,1]]
### Core Design Principle
**Stream parsing** — Bio.SeqIO processes records one at a time using a Python iterator, never loading the entire file into memory.

This is essential when parsing multi-gigabyte FASTQ files that contain millions of reads.
:::

---

# Bio.SeqIO — The SeqRecord Object

Every record parsed by Bio.SeqIO becomes a **SeqRecord** object:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### SeqRecord Anatomy

```python
from Bio import SeqIO

record = next(SeqIO.parse("reads.fastq", "fastq"))

record.id          # "SRR123456.1"
record.seq         # Seq('ACGTACGTACGT...')
record.description # full header line
record.letter_annotations["phred_quality"]
                   # [34, 35, 36, ...]
len(record)        # read length in bp
```

[[0,1]]

### FASTQ Record Structure

```
@SRR123456.1 (header / ID)
ACGTACGTACGTACGTACGT (sequence)
+
IIIIIIIIIIIIIIIIIII  (quality, ASCII-encoded)
```

Each character in the quality line encodes a Phred score:

$$Q = \text{ord}(\text{char}) - 33$$

For `I` → $\text{ord}(\texttt{I}) = 73$, so $Q = 40$.
:::

---

# Bio.SeqIO — Stream Parsing Pattern

The most important pattern: **iterate, don't load**.

```python
from Bio import SeqIO

# ✓ CORRECT — Stream parsing (constant memory)
for record in SeqIO.parse("large_file.fastq", "fastq"):
    process(record)  # each record is yielded one at a time

# ✗ WRONG — Loading everything into memory
all_records = list(SeqIO.parse("large_file.fastq", "fastq"))  # OOM risk!
```

### Practical Example: Quality Statistics

```python
from Bio import SeqIO
import statistics

qualities = []
for record in SeqIO.parse("clean_R1.fastq.gz", "fastq"):
    mean_q = statistics.mean(record.letter_annotations["phred_quality"])
    qualities.append(mean_q)

print(f"Reads: {len(qualities)}")
print(f"Mean quality: {statistics.mean(qualities):.1f}")
print(f"Reads with Q≥30: {sum(1 for q in qualities if q >= 30)}")
```

---

# Bio.SeqIO — Format Conversion

One of Bio.SeqIO's most powerful features: **converting between formats in one line**.

```python
from Bio import SeqIO

# FASTQ → FASTA (strips quality scores)
SeqIO.convert("reads.fastq", "fastq", "reads.fasta", "fasta")

# GenBank → FASTA (extracts sequence from annotations)
SeqIO.convert("genome.gbk", "genbank", "genome.fasta", "fasta")
```

### Filtering + Writing

```python
from Bio import SeqIO

# Write only high-quality, long reads
good_reads = (
    record for record in SeqIO.parse("raw.fastq", "fastq")
    if statistics.mean(record.letter_annotations["phred_quality"]) >= 30
    and len(record) >= 100
)

count = SeqIO.write(good_reads, "filtered.fastq", "fastq")
print(f"Kept {count} reads")
```

> **Key insight**: Because `good_reads` is a **generator**, this filters and writes millions of reads using almost no memory.

---

# Tool 3 — mash

## Sketch-Based Genomic Distance

**mash** implements the MinHash algorithm (from Part 3) as a production-grade command-line tool. It compresses entire genomes into compact **sketches** and estimates pairwise distances in seconds.

:::matrix { .accent cols="1fr 1fr" gap="40px" margin-top: 30px; }
[[0,0]]
### Why mash?
- Traditional alignment: $O(N^2)$ — hours for large genomes
- mash sketch + dist: $O(N)$ — seconds for the same genomes
- Accuracy within 1% of ANI (Average Nucleotide Identity)
- Works on raw reads *or* assembled genomes

[[0,1]]
### Key Use Cases
- Species identification from sequencing reads
- Contamination detection
- Outbreak clustering (epidemiology)
- Metagenome composition estimation
:::

---

# mash — The Two-Step Workflow

:::static-diagram { width=90% center }
{
  "nodes": [
    { "id": "genome", "label": "Genome\n(FASTA/FASTQ)", "x": 50,  "y": 150, "width": 150, "height": 60 },
    { "id": "kmer",   "label": "All k-mers\n(millions)",  "x": 270, "y": 150, "width": 140, "height": 60 },
    { "id": "hash",   "label": "Hash each\nk-mer",         "x": 470, "y": 150, "width": 130, "height": 60, "style": "fill: #FF9800; color: #fff;" },
    { "id": "sketch", "label": "Keep s smallest\n(sketch)", "x": 670, "y": 150, "width": 160, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "msh",    "label": ".msh file\n(~1 KB)",        "x": 890, "y": 150, "width": 120, "height": 60 }
  ],
  "links": [
    { "from": "genome", "to": "kmer" },
    { "from": "kmer",   "to": "hash" },
    { "from": "hash",   "to": "sketch" },
    { "from": "sketch", "to": "msh" }
  ]
}
:::

### Step 1: `mash sketch` — Compress

```bash
mash sketch -k 21 -s 1000 -m 2 ecoli_reads.fastq -o ecoli
#            │      │      │                         └─ output prefix (.msh)
#            │      │      └─ minimum k-mer copies (error filter)
#            │      └─ sketch size (keep 1000 smallest hashes)
#            └─ k-mer size
```

The `-m 2` flag requires each k-mer to appear **at least twice** before it is included. Sequencing errors produce random k-mers that typically appear only once — this filter suppresses them.

A ~5 MB genome is compressed into a **~1 KB sketch** — a 5000× reduction.

---

# mash sketch — Output Statistics

When sketching reads, mash prints diagnostic statistics before writing the `.msh` file:

```
Estimated genome size: 4.63941e+06
Estimated coverage:    47.13
Writing to ecoli_reads.fastq.msh...
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Estimated Genome Size

Mash counts distinct k-mers and estimates the total genome length in base pairs.

For *E. coli*: ~4.6 million bp (matches the known genome of 4.64 Mbp).

This is a quick sanity check — if you're sequencing *E. coli* but the estimate says 3 Gbp, something is **wrong** (contamination with human DNA?).

[[0,1]]

### Estimated Coverage

The average number of times each base in the genome was sequenced:

$$\text{Coverage} = \frac{\text{Total bases sequenced}}{\text{Genome size}}$$

A coverage of **47×** means each position was read ~47 times on average.

- **< 10×** — low confidence, many positions unsampled
- **20–50×** — good for most applications
- **> 100×** — diminishing returns, possible PCR duplicates
:::

---

# mash — Distance Computation

### Step 2: `mash dist` — Compare

```bash
mash dist ecoli.msh salmonella.msh
```

Output (5 tab-separated columns):
```
ecoli.fna    salmonella.fna    0.0853    8.14e-10    732/1000
│            │                 │         │           │
│            │                 │         │           └─ shared hashes / sketch size
│            │                 │         └─ p-value
│            │                 └─ Mash distance
│            └─ query ID
└─ reference ID
```

> **P-value**: The probability of observing this many shared hashes **by random chance** alone. A value near $0$ (e.g., $8.14 \times 10^{-10}$) means the similarity is statistically significant. A large p-value (e.g., $0.5$) means the match could be noise.

### The Math (from Part 3)

The Jaccard index is estimated from the shared hash count:

$$\hat{J} = \frac{|h_{\min}(A) \cap h_{\min}(B)|}{s} = \frac{732}{1000} = 0.732$$

Then converted to Mash distance:

$$D = -\frac{1}{k} \ln\left(\frac{2J}{1+J}\right) = -\frac{1}{21} \ln\left(\frac{2 \times 0.732}{1.732}\right) \approx 0.0853$$

---

# mash — Interpreting the Distance

### What Does the Mash Distance Mean?

The Mash distance $D$ approximates **one minus the Average Nucleotide Identity (ANI)**:

$$D \approx 1 - \text{ANI}$$

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Distance Thresholds

| Mash $D$ | ANI | Biological Meaning |
|----------|-----|--------------------|
| 0.00 | 100% | Identical genomes |
| < 0.05 | > 95% | Same species |
| 0.05–0.10 | 90–95% | Same genus |
| 0.10–0.20 | 80–90% | Same family |
| > 0.30 | < 70% | Distant relationship |

**Rule of thumb**: $D < 0.05$ → same species.

[[0,1]]

### Our Lecture Dataset

| Comparison | Expected $D$ |
|-----------|-------------|
| *E. coli* reads vs. *E. coli* K-12 | ~0.005 |
| *E. coli* reads vs. *Salmonella* | ~0.08 |

The mash distance clearly separates the query from the correct reference vs. the outgroup — exactly what we need for species identification!
:::

---

# mash — Two Ways to Filter Errors

When sketching raw reads (not assembled genomes), sequencing errors create **spurious k-mers**. Mash offers two approaches:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### `-m` flag (Manual Threshold)

```bash
mash sketch -k 21 -m 2 reads.fastq
```

Keep only k-mers seen **≥ m times**.

- Simple and explicit
- `-m 2` is sufficient for most bacterial datasets
- **Used in our practice notebook**

[[0,1]]

### `-r` flag (Automatic)

```bash
mash sketch -k 21 -r reads.fastq
```

Estimates coverage and auto-calculates the copy threshold $c$.

- Better for unknown coverage
- Also reports estimated genome size
- Both flags achieve the same goal: **suppress error k-mers**
:::

### Then compare against a reference:

```bash
mash dist reads.fastq.msh reference.msh
```

> **This is the full pipeline**: `fasterq-dump` downloads → `fastp` cleans → `mash sketch -m 2` compresses → `mash dist` identifies the species.

---

# Putting It All Together

### The Complete Pipeline (as used in the practice notebook)

```bash
# Step 0: Download mystery samples from NCBI SRA
fasterq-dump SRR40792821 --split-files --progress
mv SRR40792821_1.fastq Reference_Ecoli_1_raw.fastq
mv SRR40792821_2.fastq Reference_Ecoli_2_raw.fastq

# Step 1: Quality control with fastp
./fastp -i Reference_Ecoli_1_raw.fastq -I Reference_Ecoli_2_raw.fastq \
        -o Reference_Ecoli_1.fastq     -O Reference_Ecoli_2.fastq \
        -h Reference_Ecoli_report.html

# Step 2: Sketch the cleaned reads
mash sketch -k 21 -m 2 Reference_Ecoli_1.fastq

# Step 3: Compare against another sample
mash dist Reference_Ecoli_1.fastq.msh Mystery_Sample_01_1.fastq
```

```python
# Inspect any file with Bio.SeqIO
from Bio import SeqIO
for rec in SeqIO.parse("Reference_Ecoli_1.fastq", "fastq"):
    print(f"{rec.id}: length={len(rec)}, "
          f"mean_Q={sum(rec.letter_annotations['phred_quality'])/len(rec):.1f}")
    break  # peek at the first read
```

---

# Summary

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .green }
### fastp
#### Quality Control
- All-in-one FASTQ preprocessor
- Automatic adapter detection
- Sliding window quality trimming
- Interactive HTML report
- **Input**: Raw FASTQ
- **Output**: Clean FASTQ

[[0,1]]{ .blue }
### Bio.SeqIO
#### Stream Parsing
- Unified format interface
- Iterator-based (constant memory)
- SeqRecord objects
- Format conversion
- **Input**: Any sequence file
- **Output**: Python objects

[[0,2]]{ .orange }
### mash
#### Genomic Distance
- MinHash implementation
- 5000× genome compression
- Works on reads or assemblies
- Species-level identification
- **Input**: FASTA / FASTQ
- **Output**: Distance matrix
:::

---

# Next Steps

### In the Practice Session You Will:

1. **Download real data** from NCBI SRA with `fasterq-dump` — a reference *E. coli* and 7 mystery samples
2. **Clean the reads** with `fastp` — evaluate Q20/Q30 rates, adapter contamination, insert sizes
3. **Sketch & compare** with `mash sketch -m 2` and `mash dist` — identify which mystery samples are *E. coli*
4. **Classify unknowns** — match each mystery sample to one of 7 organisms:

| Organism | Mystery Sample # |
| :--- | :--- |
| *Bacillus subtilis* | ? |
| Ebola virus | ? |
| *Listeria monocytogenes* | ? |
| *Pseudomonas aeruginosa* | ? |
| *Salmonella enterica* | ? |
| SARS-CoV-2 | ? |
| *Staphylococcus aureus* (MRSA) | ? |

> **Homework**: Download reference genomes from the SRA database and build a classifier that identifies all 7 mystery organisms.
