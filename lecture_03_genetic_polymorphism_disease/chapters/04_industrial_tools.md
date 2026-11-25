:::titlepage
[[top]]
Part 4

[[title]]
Industrial Tools

[[subcaption]]
samtools · freebayes · cyvcf2 · ClinVar · pysam
:::

---

# Curriculum

- **samtools** — streaming a region of a remote CRAM, indexing, `flagstat`, `depth`.
- **SAM / BAM / CRAM** — anatomy of an aligned read: FLAG, CIGAR, MAPQ, base qualities.
- **freebayes** (executed in Colab) — haplotype-based Bayesian variant caller.
- **VCF** — header, records, `QUAL` vs. `GT` vs. `GL`, indel representation.
- **cyvcf2** — fast Python VCF parser; filtering and allele fractions.
- **ClinVar** + **tabix** — annotating variants with clinical significance.
- **pysam** — reading base qualities from the BAM to re-derive genotypes by hand.


# Roadmap

One pipeline in the practice notebook (**top**), extended in the homework (**bottom**):

:::static-diagram { width=95% center }
{
  "nodes": [
    { "id": "ftp",   "label": "1000 Genomes\nCRAM (NA12878)",  "x": 50,   "y": 80,  "width": 150, "height": 60 },
    { "id": "sam",   "label": "samtools view\n+ index",         "x": 240,  "y": 80,  "width": 150, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "bam",   "label": "BRCA1\nBAM slice",               "x": 430,  "y": 80,  "width": 130, "height": 60 },
    { "id": "fb",    "label": "freebayes\n(+ chr17 FASTA)",     "x": 600,  "y": 80,  "width": 160, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "vcf",   "label": "VCF\n(168 records)",             "x": 800,  "y": 80,  "width": 130, "height": 60 },
    { "id": "clin",  "label": "cyvcf2\n× ClinVar",              "x": 970,  "y": 80,  "width": 150, "height": 60, "style": "fill: #2196F3; color: #fff;" },

    { "id": "filt",  "label": "filter_variants\nQUAL ≥ 20, DP ≥ 15", "x": 800, "y": 250, "width": 170, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "aaf",   "label": "AAF\nhistogram",                 "x": 1010, "y": 250, "width": 120, "height": 60 },
    { "id": "pysam", "label": "pysam\n(base qualities)",        "x": 580,  "y": 250, "width": 170, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "bayes", "label": "bayesian_genotyper\n(HWE prior)", "x": 340, "y": 250, "width": 190, "height": 60, "style": "fill: #FF9800; color: #fff;" },
    { "id": "cmp",   "label": "Compare with\nfreebayes GT",     "x": 120,  "y": 250, "width": 170, "height": 60 }
  ],
  "links": [
    { "from": "ftp",   "to": "sam" },
    { "from": "sam",   "to": "bam" },
    { "from": "bam",   "to": "fb" },
    { "from": "fb",    "to": "vcf" },
    { "from": "vcf",   "to": "clin" },
    { "from": "vcf",   "to": "filt" },
    { "from": "filt",  "to": "aaf" },
    { "from": "filt",  "to": "pysam" },
    { "from": "bam",   "to": "pysam" },
    { "from": "pysam", "to": "bayes" },
    { "from": "bayes", "to": "cmp" }
  ]
}
:::

---

# The Toolbox

Variant calling is a **chain of file formats**. Each tool reads one format and writes the next:

:::matrix { .cards cols="1fr 1fr 1fr 1fr 1fr" gap="20px" margin-top: 20px; }
[[0,0]]{ .purple }
### samtools
Swiss-army knife for SAM/BAM/CRAM: slice, index, sort, count, depth.

[[0,1]]{ .green }
### freebayes
Bayesian variant caller: aligned reads + reference → VCF.

[[0,2]]{ .blue }
### cyvcf2
Python VCF parser built on `htslib` (C) — fast, NumPy arrays.

[[0,3]]{ .blue }
### ClinVar + tabix
NCBI archive of clinically interpreted variants; `.tbi` index for region queries.

[[0,4]]{ .purple }
### pysam
Python access to BAM/CRAM reads and pileups (also `htslib`).
:::

```python
# freebayes: static binary from GitHub (the apt package is often outdated)
!wget -q https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
!gunzip freebayes-1.3.6-linux-amd64-static.gz && chmod +x freebayes-1.3.6-linux-amd64-static
!mv freebayes-1.3.6-linux-amd64-static /usr/local/bin/freebayes
!apt-get -y -qq install samtools tabix        # command-line tools
!pip install -q cyvcf2 pysam                  # Python libraries
```

---

# The Data — NA12878 and BRCA1

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### The sample: NA12878
- Woman of Northern/Western European ancestry (CEU), **1000 Genomes Project**
- Reference sample of **Genome in a Bottle** (HG001) — the most thoroughly validated human genome
- Used to benchmark every new sequencer and caller
- Our file: **targeted exome** capture, Illumina, 76 bp paired-end reads, aligned to **GRCh38** with BWA-MEM

[[0,1]]

### The locus: BRCA1
- Tumor-suppressor gene on chromosome 17 (minus strand)
- Repairs DNA double-strand breaks by homologous recombination
- Pathogenic variants → high lifetime risk of breast & ovarian cancer
- Region used: `chr17:43044295-43170245` (≈126 kb, BRCA1 + flanks)
:::

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| File | What it is | Size in our run |
| :--- | :--- | :--- |
| `NA12878...exome.cram` (remote) | Whole exome, all chromosomes | many GB — never downloaded |
| `NA12878_BRCA1.bam` | Only reads in the BRCA1 region | **2.0 MB** |
| `chr17_ref.fa` | GRCh38 chromosome 17 | 84 MB |
| `clinvar.vcf.gz` + `.tbi` | All ClinVar variants (release 2026-10-11) | 200 MB |
:::

---

# The Format Chain

```
 sequencer        aligner (BWA-MEM)          variant caller            annotation
  FASTQ    ───────────────────►  SAM/BAM/CRAM  ─────────────────►  VCF  ─────────────►  VCF + ClinVar
 (Lecture 1)                     "which read sits where"           "where the sample      "what it means
                                                                    differs from ref"       clinically"
```

:::matrix { cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]

### SAM
- Plain **text**, one read per line
- Human-readable
- Huge (≈ size of FASTQ × 2)

[[0,1]]

### BAM
- **Binary**, BGZF block-compressed SAM
- Random access with a `.bai` index
- Standard working format

[[0,2]]

### CRAM
- Stores only **differences from the reference**
- 40–60% smaller than BAM
- **Cannot be decoded without the exact reference** FASTA
:::

> 1000 Genomes distributes CRAM. `samtools view -b` decodes the CRAM and re-encodes the slice as BAM, so the downstream tools see an ordinary BAM.

---

# Tool 0 — samtools view: Streaming a Region

```python
cram_url = "https://ftp.1000genomes.ebi.ac.uk/.../NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram"
ref_url  = "https://ftp.1000genomes.ebi.ac.uk/.../GRCh38_full_analysis_set_plus_decoy_hla.fa"

!samtools view -b -T {ref_url} {cram_url} chr17:43044295-43170245 > NA12878_BRCA1.bam
!samtools index NA12878_BRCA1.bam              # → NA12878_BRCA1.bam.bai
!samtools faidx {ref_url} chr17 > chr17_ref.fa  # extract one chromosome
!samtools faidx chr17_ref.fa                    # → chr17_ref.fa.fai
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### How can it fetch 126 kb out of many GB?
1. Download the small index (`.crai`)
2. Look up which compressed blocks overlap the region
3. Request **only those bytes** (HTTP range requests)
4. Decode the blocks using the reference

| Flag | Meaning |
| :--- | :--- |
| `-b` | write BAM (default: SAM text) |
| `-T` | reference FASTA for CRAM decoding |
| `chr:a-b` | region, **1-based, inclusive** |

[[0,1]]

### Why `-T {ref_url}`?
Without it, samtools looks up the reference by the **MD5 checksum** stored in the CRAM header and tries to download it from the EBI reference server. That download can hang without any error message.

### The `.fai` / `.bai` indexes
- **Always** index after creating a BAM or FASTA
- freebayes, pysam, IGV all refuse un-indexed files (or read them slowly)
- The index must be **newer** than the data file
:::

---

# Anatomy of a SAM Record

`samtools view NA12878_BRCA1.bam | head -n 2` — real output, columns split for readability:

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| # | Field | Read 1 | Read 2 (its mate) | Meaning |
| :--- | :--- | :--- | :--- | :--- |
| 1 | QNAME | `SRR098401.102069364` | same | read (pair) name |
| 2 | FLAG | `99` | `147` | bit field — see next slide |
| 3 | RNAME | `chr17` | `chr17` | reference sequence |
| 4 | POS | `43044281` | `43044336` | leftmost aligned base, **1-based** |
| 5 | MAPQ | `60` | `60` | mapping quality (Phred; 60 = unique) |
| 6 | CIGAR | `76M` | `76M` | 76 bases aligned, no indels/clipping |
| 7–9 | RNEXT, PNEXT, TLEN | `=`, `43044336`, `131` | `=`, `43044281`, `−131` | mate position, fragment length |
| 10 | SEQ | `ACAGTCATTCATGG…` | `ACAACATTTACAAA…` | bases (forward strand of ref) |
| 11 | QUAL | `BBBBBBBBBB…<BB` | `<<<<BBBBBB…<<B<<` | Phred + 33, one char per base |
| 12+ | tags | `NM:i:0 MD:Z:76` | `NM:i:1 MD:Z:55G20` | edit distance; mismatch map |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Reading the tags
- `NM:i:1` — one difference from the reference
- `MD:Z:55G20` — 55 matches, reference base **G** replaced, 20 matches
- `RG:Z:SRR098401` — read group (sequencing run)

[[0,1]]

### Note
Read 1 starts at 43,044,281 — **before** the region start (43,044,295). `samtools view` returns every read that **overlaps** the region, not only reads fully inside it.
:::

---

# FLAG and CIGAR

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### FLAG = sum of bits

| Bit | Value | Meaning | 99 | 147 |
| :--- | :--- | :--- | :--- | :--- |
| 0x1 | 1 | paired | ✓ | ✓ |
| 0x2 | 2 | properly paired | ✓ | ✓ |
| 0x10 | 16 | read on **reverse** strand | | ✓ |
| 0x20 | 32 | mate on reverse strand | ✓ | |
| 0x40 | 64 | first in pair | ✓ | |
| 0x80 | 128 | second in pair | | ✓ |
| 0x100 | 256 | secondary alignment | | |
| 0x400 | 1024 | **PCR/optical duplicate** | | |

99 = 1+2+32+64 · 147 = 1+2+16+128

[[0,1]]

### CIGAR operations

| Op | Meaning |
| :--- | :--- |
| `M` | aligned (match **or** mismatch!) |
| `I` | insertion in the read |
| `D` | deletion from the read |
| `S` | soft clip (bases kept, not aligned) |
| `N` | skipped region (RNA-seq introns) |

Most common CIGARs in our BAM:

```
76M      40,892 reads   (plain alignment)
71M5S        71         (adapter / low-quality tail clipped)
5S71M        62
```

`M` does not tell you whether a base matches — use `MD`/`NM` or compare with the reference.
:::

---

# samtools flagstat — Health Check

```
43196 + 0 in total (QC-passed reads + QC-failed reads)
43167 + 0 primary
29 + 0 secondary
4789 + 0 duplicates
42758 + 0 mapped (98.99% : N/A)
41871 + 0 properly paired (97.00% : N/A)
438 + 0 singletons (1.01% : N/A)
226 + 0 with mate mapped to a different chr
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### What to look at
- **Mapped %** — low → wrong reference or contamination
- **Properly paired %** — low → library or alignment problems
- **Duplicates** — 4,789 reads (**11%**) are PCR/optical copies of another fragment

[[0,1]]

### Why duplicates matter for variant calling
A PCR copy is **not independent evidence**. If one molecule with a polymerase error is amplified 10×, it looks like 10 reads supporting a variant.

→ Duplicates are **marked** (FLAG 0x400) by the aligner pipeline, and freebayes **skips them by default**.
:::

---

# Base Qualities in This BAM

Phred recap (Lecture 1): $Q = -10\log_{10} \varepsilon$, stored as `chr(Q + 33)`.

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Only 8 distinct values appear

| Char | $Q$ | Error $\varepsilon$ | Share of bases |
| :--- | :--- | :--- | :--- |
| `!` / `"` | 0 / 1 | ~1 | < 0.2% |
| `'` | 6 | 25% | 12% |
| `0` | 15 | 3.2% | 1.5% |
| `7` | 22 | 0.63% | 2% |
| `<` | 27 | 0.20% | 9% |
| `B` | 33 | 0.05% | **67%** |
| `F` | 37 | 0.02% | 8% |

(first ~3,000 reads of the slice)

[[0,1]]

### Why so few values?
To make the files smaller, the qualities were **binned**: every score is rounded to one of a few levels. Compression is much better and calling accuracy barely changes.

### Why you care
- `ε` for every read enters the **likelihood** $P(D \mid G)$ (Part 3)
- Some tools **silently drop** low-quality bases: e.g. pysam's `pileup()` defaults to `min_base_quality=13` → all **Q6** bases vanish
- freebayes default `-q 0` keeps them
:::

---

# Read Depth — samtools depth

```python
!samtools depth -a -r chr17:43044295-43170245 NA12878_BRCA1.bam > depth.txt
depth_df = pd.read_csv("depth.txt", sep="\t", header=None, names=["CHROM", "POS", "DEPTH"])
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Our BRCA1 slice

| Statistic | Value |
| :--- | :--- |
| Region length | 125,951 bp |
| Positions with ≥ 1 read | 116,757 |
| **Median** depth (covered positions) | **4×** |
| Maximum depth | 1,111× |
| Positions ≥ 10× | 18,266 (15%) |
| Positions ≥ 15× | 14,954 (12%) |

[[0,1]]

### Spikes and valleys
Hybridization capture enriches **exons** (Part 2). Depth is in the hundreds on targeted exons and falls to a few **off-target** reads in introns.

→ Only ~12% of the locus has enough depth for a confident call. Most raw "variants" sit in the low-depth valleys.

### Two `samtools depth` details
- By default it skips duplicates, secondary and QC-failed reads
- Without `-a`, it prints **only covered positions** — causing `plt.plot` to draw straight lines across zero-depth introns; `-r` restricts output to the target region
:::

---

# Tool 1 — freebayes

## Haplotype-Based Bayesian Variant Calling

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]

### What it does at each candidate locus
1. Collect all reads overlapping the locus (skip duplicates, MAPQ < 1)
2. Build candidate **alleles / short haplotypes** from what the reads show
3. For every genotype $G$ (ploidy 2 → 0/0, 0/1, 1/1): compute $P(D \mid G)$ from base qualities
4. Combine with a population **prior** $P(G)$ → posterior (Part 3)
5. Report loci where "not homozygous-reference" is plausible

"Haplotype-based": nearby SNPs and indels can be reported **together** as one *complex* allele instead of separately.

[[0,1]]

### Why does it need the reference FASTA?
- To know the **REF allele** at each position — the BAM stores read sequences, not the reference
- To construct candidate haplotypes (reference + read differences)
- To write REF / ALT columns of the VCF and **left-align** indels

```python
!freebayes -f chr17_ref.fa \
           --min-alternate-fraction 0.2 \
           NA12878_BRCA1.bam > brca1_variants.vcf
```
:::

---

# freebayes — Key Options

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Option | Default | Notebook | Meaning |
| :--- | :--- | :--- | :--- |
| `-f FILE` | — | `chr17_ref.fa` | reference FASTA (indexed) — **required** |
| `-F`, `--min-alternate-fraction` | 0.05 | **0.2** | ≥ 20% of reads at a locus must show the alt allele to evaluate it |
| `-C`, `--min-alternate-count` | 2 | 2 | ≥ 2 reads must show the alt allele |
| `-p`, `--ploidy` | 2 | 2 | diploid human autosome |
| `-m`, `--min-mapping-quality` | 1 | 1 | ignore reads with MAPQ below this |
| `-q`, `--min-base-quality` | 0 | 0 | ignore bases below this Phred value |
| `-4`, `--use-duplicate-reads` | off | off | duplicates excluded unless set |
| `-r chr:start-end` | whole file | — | restrict to a region |
| `-t FILE.bed` | — | — | restrict to capture targets (exome) |
| `--gvcf` | off | off | also report coverage at non-variant sites |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Why raise `-F` to 0.2?
A true heterozygote gives ~50% alt reads, a homozygote ~100%. At 5% you mostly pick up **sequencing errors** and **mapping artifacts**; 20% removes them before the expensive Bayesian step.

[[0,1]]

### What freebayes does *not* do
It writes every evaluated locus — the `FILTER` column stays `.` for all 168 records. **Filtering is your job** (homework Task 1).
:::

---

# The VCF File — Header

```
##fileformat=VCFv4.2
##source=freeBayes v1.3.6
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=RO,Number=1,Type=Integer,Description="Count of full observations of the reference haplotype.">
##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal ...">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods ...">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  NA12878
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Two kinds of annotations
- **INFO** — about the **site** (all samples together)
- **FORMAT** — about **each sample** (one column per sample)

Every `ID` used in the body **must** be declared here — this is the file's own documentation.

[[0,1]]

### The `Number=` field

| Value | Count of values |
| :--- | :--- |
| `1` | exactly one |
| `A` | one per **ALT** allele |
| `R` | one per allele (REF + ALTs) |
| `G` | one per **genotype** (3 for diploid bi-allelic) |
| `.` | variable |
:::

---

# The VCF File — A Record

A real record (rs3092994), INFO shortened:

```
#CHROM POS      ID REF ALT QUAL    FILTER INFO                              FORMAT              NA12878
chr17  43063808 .  C   T   312.156 .      AB=0.475;DP=40;RO=21;AO=19;...    GT:DP:AD:RO:QR:AO:QA:GL  0/1:40:21,19:21:637:19:610:-43.17,0,-45.58
```

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Field | Value | Meaning |
| :--- | :--- | :--- |
| `POS` | 43063808 | **1-based** position of the first REF base |
| `REF` / `ALT` | C / T | reference allele / alternate allele(s), comma-separated |
| `QUAL` | 312.2 | Phred-scaled confidence that the site is **variant** at all |
| `GT` | `0/1` | genotype: allele indices (0 = REF, 1 = first ALT); `/` unphased, `\|` phased |
| `DP`, `AD` | 40, `21,19` | depth; reads per allele (REF, ALT) |
| `RO`, `AO` | 21, 19 | freebayes' own ref / alt observation counts |
| `QR`, `QA` | 637, 610 | **sums** of base qualities of ref / alt observations |
| `GL` | −43.17, 0, −45.58 | $\log_{10} P(D \mid G)$ for 0/0, 0/1, 1/1 (best scaled to 0) |
:::

---

# QUAL ≠ Genotype Confidence

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Three different numbers

| Field | Question it answers | Level |
| :--- | :--- | :--- |
| `QUAL` | *Is there any variant here?* $-10\log_{10} P(\text{hom-ref} \mid D)$ | site |
| `GL` | *How well does each genotype explain the reads?* | sample |
| `GQ` | *How sure is the chosen genotype?* | sample |

A site can have **high QUAL** but an **uncertain genotype** (sure it is variant, unsure if 0/1 or 1/1).

[[0,1]]

### Reading GL by hand
For rs3092994: GL = (−43.17, 0, −45.58)

$$ \frac{P(D \mid 0/1)}{P(D \mid 0/0)} = 10^{43.17} $$

The reads are $10^{43}$ times better explained by a heterozygote. With a flat prior, $P(0/1 \mid D) \approx 1$.

### ⚠ Declared ≠ present
`GQ` is declared in the header but freebayes 1.3.6 did **not** write it in FORMAT (`GT:DP:AD:RO:QR:AO:QA:GL`):

```python
variant.format("GQ")     # → None
```
:::

---

# Tool 2 — cyvcf2

```python
from cyvcf2 import VCF

vcf = VCF("brca1_variants.vcf")
vcf.samples                          # ['NA12878']
for v in vcf:                        # streams records one by one (like SeqIO.parse)
    ...
```

Real values for the first variant that passes the filter (rs3092994):

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Site-level attributes

```python
v                    # Variant(chr17:43063808 C/T)
v.CHROM, v.POS       # 'chr17', 43063808  (1-based)
v.start, v.end       # 43063807, 43063808 (0-based)
v.REF                # 'C'
v.ALT                # ['T']        ← always a LIST
v.QUAL               # 312.156...
v.FILTER             # None         ← '.' in the file
v.INFO.get("DP")     # 40
v.INFO.get("TYPE")   # 'snp'
v.is_snp, v.is_indel # True, False
```

[[0,1]]

### Sample-level attributes

```python
v.genotypes          # [[0, 1, False]]
                     #   allele1, allele2, phased?
v.gt_types           # array([1])
                     #   0=HOM_REF 1=HET 2=UNKNOWN 3=HOM_ALT
v.gt_bases           # array(['C/T'])
v.format("RO")       # array([[21]], dtype=int32)
v.format("AO")       # array([[19]], dtype=int32)
v.format("GL")       # array([[-43.17, 0., -45.58]])
v.format("GQ")       # None  (not in this file)
```
:::

---

# cyvcf2 — Shapes and Missing Values

`v.format(key)` always returns a **2-D array**: one row per **sample**, one column per **value**.

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Indexing pattern used in the notebook

```python
ro = v.format("RO")[0][0]     # sample 0, value 0
ao = v.format("AO")[0][0]     # first ALT allele only
gl = v.format("GL")[0]        # 3 values: 0/0, 0/1, 1/1

a1, a2, phased = v.genotypes[0]
gt = f"{a1}/{a2}"             # '0/1'
```

| `Number=` | shape for 1 sample |
| :--- | :--- |
| `1` (DP, RO) | (1, 1) |
| `A` (AO) | (1, n_alt) |
| `G` (GL) | (1, 3) for bi-allelic |

[[0,1]]

### Defensive parsing
- Field absent from the record → `None` (not an error!)
- Missing allele in GT (`./.`) → `-1` in `genotypes`
- Multi-allelic site → `ALT` has several entries, `AO` several columns (none in our 168 records, but common genome-wide)

```python
ro_arr = v.format("RO")
ro = ro_arr[0][0] if ro_arr is not None else 0
```
:::

> cyvcf2 can only jump to a **region** (`vcf("chr17:1-1000")`) if the file is bgzip-compressed and indexed. Our plain-text `brca1_variants.vcf` can only be read start to end.

---

# ⚠ Pitfall — Coordinate Systems

Every tool on today's list picks its own convention. The **same** SNP, rs3092994:

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Where | Convention | How rs3092994 is addressed |
| :--- | :--- | :--- |
| VCF `POS`, `cyvcf2 v.POS` | 1-based | `43063808` |
| SAM `POS` column | 1-based | first base of read |
| `cyvcf2 v.start` / `v.end` | 0-based, half-open | `43063807` / `43063808` |
| samtools / htslib region string | 1-based, inclusive | `"chr17:43063808-43063808"` |
| cyvcf2 region query on ClinVar | 1-based, inclusive | `clinvar("17:43063808-43063808")` |
| `pysam` `fetch()` / `pileup()` | 0-based, half-open | `bam.fetch("chr17", 43063807, 43063808)` |
| `pysam` `PileupColumn.reference_pos` | 0-based | `43063807` |
| BED files (`freebayes -t`) | 0-based, half-open | `chr17  43063807  43063808` |
:::

> **Rule of thumb**: *strings* typed by humans (`"chr17:a-b"`) are 1-based; *integer arguments* of Python APIs and BED files are 0-based. Off-by-one errors here produce **plausible-looking wrong answers**: you analyse the neighbouring base.

---

# Filtering — From 168 Records to 25

Homework Task 1: keep a variant only if `QUAL ≥ 20` **and** `DP ≥ 15`.

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Where the 168 records go

| QUAL | DP | Records |
| :--- | :--- | :--- |
| ≥ 20 | ≥ 15 | **25 — kept** |
| ≥ 20 | < 15 | 68 |
| < 20 | < 15 | 72 |
| < 20 | ≥ 15 | 3 |

Median DP of all records: **5**. Most raw calls come from off-target reads in introns.

Record types: 152 SNPs, 7 deletions, 6 insertions, 3 complex.

[[0,1]]

### Why both filters?

**High QUAL, low DP** — confident in *some* variant, but genotype rests on a handful of reads:
```
43088260  T>C  QUAL=143.9  DP=5   RO=0  AO=5   1/1
```

**High DP, low QUAL** — plenty of reads, but freebayes itself rejects the variant:
```
43057017  T>G  QUAL=0.0    DP=375 RO=288 AO=80  0/0
          all 80 alt reads on the forward strand
```
21% alt reads passed the `-F 0.2` pre-filter, but the model prefers "error/artifact".
:::

---

# Why Depth Matters

For a true heterozygote, each read samples one of the two chromosomes with probability ½. The chance that **all** $n$ reads come from the **same** chromosome:

$$ P(\text{het looks homozygous}) = 2 \cdot \left(\tfrac{1}{2}\right)^n $$

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

| Depth $n$ | $P$ (one allele missed) |
| :--- | :--- |
| 2 | 50% |
| 5 | 6.3% |
| 10 | 0.2% |
| 15 | 0.006% |
| 30 | ≈ 2 × 10⁻⁹ |

[[0,1]]

### Consequences
- At 2–5× a heterozygote is often called **1/1** (or missed when the alt allele is unseen)
- Base-quality errors weigh more when there are few reads to outvote them
- `DP ≥ 15` makes the sampling error negligible — **if** the reads are independent (duplicates removed!)

Many freebayes records with `DP = 2` in our VCF are genotyped `1/1` on two reads.
:::

---

# Allele Fraction Spectrum

Homework Task 2: $\displaystyle AAF = \frac{AO}{RO + AO}$ for each of the 25 kept variants.

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Expected peaks
- **0.5** — heterozygous (one copy of two)
- **1.0** — homozygous alternate
- 0.0 — homozygous reference (freebayes rarely reports these)

### Our NA12878 BRCA1 result

| AAF bin | Variants |
| :--- | :--- |
| 0.0 – 0.2 | 0 |
| 0.2 – 0.4 | 2 |
| **0.4 – 0.6** | **22** |
| 0.6 – 0.9 | 0 |
| **≥ 0.9** | **1** |

[[0,1]]

### Reading the spectrum
- 24 of 25 calls are `0/1`: NA12878 carries one common BRCA1 haplotype on one chromosome and the reference haplotype on the other
- Peaks are **wide**: with 40 reads, 21/19 or 25/15 are both ordinary sampling outcomes
- Systematic shift **below** 0.5 → *reference bias*: reads with the alt allele align slightly worse or are captured less efficiently
- AAF near 0.1–0.3 in a germline sample → artifact, contamination — or a **somatic/mosaic** variant in tumour data
:::

---

# Tool 3 — ClinVar + tabix

**ClinVar** (NCBI) aggregates clinical interpretations of human variants submitted by laboratories and expert panels. It is distributed as a bgzipped VCF with a tabix index:

```python
!wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
!wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

clinvar_db = VCF("clinvar.vcf.gz")
for cv in clinvar_db("17:43071077-43071077"):     # region query via .tbi — no full scan
    print(cv.REF, cv.ALT, cv.INFO.get("CLNSIG"))
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### bgzip + tabix
- **bgzip**: gzip in independent 64 kB blocks → still readable by `zcat`
- **tabix**: index mapping genomic intervals → block offsets
- Query time is independent of file size

### ⚠ Contig names
Our BAM/VCF use `chr17` (UCSC style); ClinVar uses `17` (Ensembl style). Hence `chrom.replace("chr", "")`. Wrong name → **zero hits, no error**.

[[0,1]]

### Key INFO fields (real values, rs1799966)

| Field | Value |
| :--- | :--- |
| `CLNSIG` | `Benign` |
| `CLNREVSTAT` | `reviewed_by_expert_panel` |
| `CLNDN` | `Breast-ovarian_cancer,_familial,...` |
| `MC` | `SO:0001583\|missense_variant` |
| `GENEINFO` | `BRCA1:672` |
| `RS` | `1799966` (dbSNP rs ID) |
:::

---

# Match the Allele, Not Just the Position

ClinVar records at **one** position, chr17:43,071,077 — where NA12878 is heterozygous T/C:

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| REF | ALT | ClinVar CLNSIG | Consequence | Review status |
| :--- | :--- | :--- | :--- | :--- |
| T | **C** | **Benign** | missense | expert panel |
| T | A | Benign / Likely benign | missense | multiple submitters |
| T | G | Uncertain significance | missense | single submitter |
| T | TC | **Pathogenic** | frameshift | expert panel |
| TCTG | CCT | **Pathogenic** | frameshift | single submitter |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

```python
for cv_var in clinvar_db(region):
    if not (cv_var.REF == ref and alt in cv_var.ALT):
        continue          # same position, different variant
    clnsig = cv_var.INFO.get("CLNSIG")
```

[[0,1]]

Matching by position only would report this healthy reference genome as carrying **pathogenic BRCA1 frameshifts**. A position is not a variant; **(CHROM, POS, REF, ALT)** is.
:::

---

# ⚠ Pitfall — Indel Representation

The same deletion can be written in many ways. freebayes and ClinVar disagree on **all 16** of our indel/complex records:

```
freebayes:   43097346  TAAAAAAAT  > TAAAAAAT       (haplotype with extra context)
ClinVar:     43097346  TA         > T              rs8176144, Benign, expert panel
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Normalization = trim + left-align
1. Remove bases shared at the **end** of REF and ALT
2. Remove shared bases at the **start** (keep one anchor base, shift POS)
3. Shift the indel as far **left** as the repeat allows

After steps 1–2 alone, **8 of 16** indels match ClinVar exactly (all Benign). The remaining ones are complex alleles that must be split into primitives.

[[0,1]]

### In practice
```bash
bcftools norm -f chr17_ref.fa -a -m -any \
     brca1_variants.vcf -o brca1_norm.vcf
#  -f  reference, needed to left-align
#  -a  atomize complex alleles into SNPs/indels
#  -m -any  split multi-allelic records
```
Always normalize **both** files the same way before comparing call sets or annotating.
:::

---

# ClinVar Results for NA12878

The notebook loop (exact REF/ALT match, ClinVar release 2026-10-11): **58** of 168 records annotated. **56 Benign**, **2 with no germline classification** (`CLNSIG` is `None`).

Of the **25 filtered** variants, 19 are in ClinVar:

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Position | Change | GT | DP | dbSNP | Consequence | CLNSIG |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 43071077 | T>C | 0/1 | 322 | rs1799966 | missense | Benign ★★★ |
| 43082453 | A>G | 0/1 | 406 | rs1060915 | synonymous | Benign ★★★ |
| 43091983 | T>C | 0/1 | 150 | rs16942 | missense | Benign ★★★ |
| 43092418 | T>C | 0/1 | 490 | rs16941 | missense | Benign ★★★ |
| 43092919 | G>A | 0/1 | 375 | rs799917 | missense | Benign ★★★ |
| 43093220 | A>G | 0/1 | 536 | rs16940 | synonymous | Benign ★★★ |
| 43093449 | G>A | 0/1 | 276 | rs1799949 | synonymous | Benign ★★★ |
| + 10 intronic | | 0/1 | 15–134 | | intron | Benign ★★★ |
| 43138657, 43139936 | | 0/1 | 45, 35 | — | intron | `None` |
:::

> **Missense ≠ pathogenic.** Four amino-acid changes, all common in the population and classified Benign by the ENIGMA expert panel (★★★). Handle `CLNSIG = None` explicitly; `str(None)` prints as a fake classification.

---

# Tool 4 — pysam

`pysam` wraps htslib (the C library under samtools) so you can read alignments **inside Python**:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Reads overlapping a position: `fetch`

```python
import pysam
bam = pysam.AlignmentFile("NA12878_BRCA1.bam")

pos0 = 43063808 - 1                     # VCF POS → 0-based
for read in bam.fetch("chr17", pos0, pos0 + 1):
    if read.is_duplicate or read.is_secondary:
        continue
    # map reference position → position inside the read
    pairs = dict((r, q) for q, r in
                 read.get_aligned_pairs(matches_only=True))
    qpos = pairs.get(pos0)              # None → deletion here
    if qpos is None:
        continue
    base = read.query_sequence[qpos]
    q    = read.query_qualities[qpos]   # already integers
```

[[0,1]]

### Columns of bases: `pileup`

```python
for col in bam.pileup("chr17", pos0, pos0 + 1,
                      truncate=True,
                      min_base_quality=0):
    for pr in col.pileups:
        if pr.is_del or pr.is_refskip:
            continue
        aln  = pr.alignment
        base = aln.query_sequence[pr.query_position]
        q    = aln.query_qualities[pr.query_position]
```

Real result at rs3092994 (C>T):
```
ref C quals: 27 33 33 33 27 27 33 6 33 ... (21 reads)
alt T quals: 33 33 33 33 33 33 33 33 ... (19 reads)
```
:::

---

# ⚠ Pitfall — Counting Reads

Reference/alternate counts at four kept SNPs, obtained three ways:

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Position | freebayes RO / AO | naive `fetch()` (all reads) | `pileup()` defaults | `fetch()` skipping duplicates |
| :--- | :--- | :--- | :--- | :--- |
| 43063808 | 21 / 19 | 24 / 20 | 20 / 19 | **21 / 19** |
| 43071077 | 165 / 157 | 182 / 174 | 137 / 138 | **165 / 157** |
| 43082453 | 230 / 176 | 265 / 197 | 200 / 144 | **230 / 176** |
| 43093220 | 265 / 268 | 299 / 321 | 231 / 248 | **265 / 268** |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Why the differences?
- **naive fetch** counts PCR **duplicates** → too many reads
- **pileup defaults** skip duplicates, **but** also drop bases with $Q < 13$ (`min_base_quality=13`) → all the **Q6** bases disappear
- freebayes: no duplicates, **all** base qualities (`-q 0`)

[[0,1]]

### Lesson for Task 4
To compare your genotyper with freebayes **fairly**, feed it the **same evidence**: same read filters, same base-quality threshold. Otherwise disagreements reflect the input data, not the algorithm.
:::

---

# From Base Qualities to a Genotype

The notebook's likelihoods (Part 8), in **log space** to avoid underflow ($0.5^{22} \approx 2 \times 10^{-7}$, and real sites have hundreds of reads):

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]

```python
def compute_likelihoods(ref_quals, alt_quals):
    ref_eps = 10 ** (-np.array(ref_quals) / 10.0)
    alt_eps = 10 ** (-np.array(alt_quals) / 10.0)
    # 0/0: ref reads correct, alt reads are errors
    ll_00 = np.sum(np.log(1 - ref_eps)) + np.sum(np.log(alt_eps))
    # 0/1: each read comes from either copy, p = 1/2
    ll_01 = (len(ref_quals) + len(alt_quals)) * np.log(0.5)
    # 1/1: alt reads correct, ref reads are errors
    ll_11 = np.sum(np.log(ref_eps)) + np.sum(np.log(1 - alt_eps))
    return ll_00, ll_01, ll_11
```

Normalize with the **log-sum-exp** trick: subtract the maximum before `np.exp`.

[[0,1]]

### Adding the prior (Task 3)
Hardy–Weinberg with population allele frequency $f$:

$$ \log P(G \mid D) = \log P(D \mid G) + \log P(G) + \text{const} $$

| $G$ | $P(G)$ |
| :--- | :--- |
| 0/0 | $(1-f)^2$ |
| 0/1 | $2f(1-f)$ |
| 1/1 | $f^2$ |

A correct implementation, fed freebayes-equivalent read counts, reproduces freebayes' genotype at all 20 kept SNPs tested.
:::

---

# The Prior Matters When Data Are Weak

Posterior genotype probabilities, all reads Q30:

:::matrix { .smaller-table cols="1fr" }
[[0,0]]
| Evidence | Prior | $P(0/0)$ | $P(0/1)$ | $P(1/1)$ | Call |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 12 ref, 10 alt | flat | ~10⁻²⁴ | **1.000** | ~10⁻³⁰ | 0/1 |
| 0 ref, **2** alt | flat | 0.000 | 0.200 | **0.800** | 1/1 |
| 0 ref, **2** alt | HWE, $f = 0.01$ (rare) | 0.000 | **0.980** | 0.020 | **0/1** |
| 0 ref, **2** alt | HWE, $f = 0.35$ (common) | 0.000 | 0.482 | **0.518** | 1/1 (barely) |
| 0 ref, **20** alt | HWE, $f = 0.01$ | 0.000 | 0.000 | **1.000** | 1/1 |
:::

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]

### Strong data overwhelm the prior
With 20+ reads the likelihood differences are hundreds of log-units; any reasonable prior is irrelevant.

[[0,1]]

### Weak data are decided by the prior
At 2 reads the call **flips** with $f$. Our VCF has such a site: `43044391 G>A, DP=2, GT=1/1`, a common benign 3′-UTR SNP (rs12516). The DP filter exists exactly for these cases.
:::

---

# Putting It All Together

```bash
# 1. Stream BRCA1 CRAM slice → local indexed BAM + chr17 reference
samtools view -b -T "$ref_url" "$cram_url" chr17:43044295-43170245 > NA12878_BRCA1.bam
samtools index NA12878_BRCA1.bam
samtools faidx "$ref_url" chr17 > chr17_ref.fa && samtools faidx chr17_ref.fa

# 2. Inspect read statistics and per-base depth
samtools flagstat NA12878_BRCA1.bam
samtools depth -a -r chr17:43044295-43170245 NA12878_BRCA1.bam > depth.txt

# 3. Call variants with freebayes (min 20% alt fraction)
freebayes -f chr17_ref.fa --min-alternate-fraction 0.2 NA12878_BRCA1.bam > brca1_variants.vcf
```

```python
# 4. Annotate exact alleles against ClinVar (strip 'chr' prefix!)
from cyvcf2 import VCF
clinvar = VCF("clinvar.vcf.gz")
for v in VCF("brca1_variants.vcf"):
    for cv in clinvar(f"{v.CHROM.replace('chr', '')}:{v.POS}-{v.POS}"):
        if cv.REF == v.REF and v.ALT[0] in cv.ALT:
            print(v.POS, v.REF, v.ALT[0], cv.INFO.get("CLNSIG") or "no classification")

# ---- HOMEWORK TASKS ----
# Task 1: filter_variants() (QUAL >= 20, DP >= 15)    Task 2: AAF histogram
# Task 3: bayesian_genotyper() with HWE priors        Task 4: pysam re-genotyping vs freebayes
```

---

# Summary

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### samtools & BAM
#### Reads
- CRAM needs the reference (`-T`)
- Indexed region streaming
- FLAG bits, CIGAR, MD/NM
- Duplicates: 11% here
- Binned qualities (Q6…Q37)
- **Output**: BRCA1 BAM slice

[[0,1]]{ .green }
### freebayes & VCF
#### Calls
- Reference FASTA required
- `-F 0.2` pre-filter
- No filtering → you filter
- QUAL (site) vs GL/GQ (genotype)
- Normalize indels before comparing
- **Output**: 168 → 25 variants

[[0,2]]{ .blue }
### cyvcf2 · pysam · ClinVar
#### Interpretation
- `format()` → (samples × values)
- 1-based strings, 0-based integers
- Same evidence for fair comparison
- Match CHROM, POS, REF **and** ALT
- `chr17` vs `17`, `CLNSIG = None`
- **Output**: all Benign — healthy reference
:::

---

# Next Steps

### In the Practice Session You Will:

1. **Stream** the BRCA1 region of NA12878 from the 1000 Genomes CRAM with `samtools view`, and inspect it with `flagstat` and `depth`
2. **Call** variants with `freebayes` against the GRCh38 chr17 reference
3. **Parse** the VCF with `cyvcf2` and **annotate** it against ClinVar
4. **Compute** genotype posteriors by hand from Phred base qualities

### Homework — Variant Analysis

1. `filter_variants()` — QUAL ≥ 20 and DP ≥ 15
2. Allele-fraction histogram of the filtered variants
3. `bayesian_genotyper()` with Hardy–Weinberg priors
4. Re-genotype the first 20 filtered SNPs from BAM base qualities with `pysam` and compare with freebayes
5. Interpretation questions

> **Questions to keep in mind**: Which reads count as evidence? How many reads are enough? When does the prior decide the answer? Why can a missense variant be benign?
