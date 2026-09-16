Lecture-by-Lecture Breakdown

Lecture 1: The Molecular Blueprint (DNA Replication & Sequencing)
  - Biological Mechanism: DNA chemical structure, complementary base pairing, phosphodiester bond stability, and the replication mechanics of DNA polymerase.
  - Instrumental Method: Illumina Sequencing-by-Synthesis (reversible fluorescent terminators, laser scanning) and Oxford Nanopore (ionic current disruptions). Converting optical/electrical signals into Phred quality scores ($Q = -10 \log_{10} P$).
  - Underlying Algorithm: $k$-mer decomposition, hash functions, and MinHash sketching. Estimating the Jaccard similarity index:$$J(A, B) = \frac{\vert{}A \cap B\vert{}}{\vert{}A \cup B\vert{}} \approx \frac{\vert{}\min_s(h(A)) \cap \min_s(h(B))\vert{}}{s}$$ and converting to Mash distance: $$D \approx -\frac{1}{k} \ln\left(\frac{2J}{1+J}\right)$$
  - Industrial Tool: 
     - mash (installed via !apt-get install mash in Colab) for sketches and distance calculations (mash sketch, mash dist).
     - fastp (run via ! in Colab) for adapter trimming and Bio.SeqIO for stream parsing.
  - Datasets: 
     - Query reads: 10 MB Escherichia coli paired-end FASTQ subset from NCBI SRA (e.g., Run SRR...).
     - Target reference: Escherichia coli str. K-12 substr. MG1655 complete genome (.fna.gz, ~1.4 MB).
     - Outgroup reference: Salmonella enterica serovar Typhimurium complete genome (.fna.gz, ~1.5 MB).

* Lecture 2: Evolutionary Divergence (Homology & Comparative Genomics)
  - Biological Mechanism: Spontaneous point mutations, transition vs. transversion bias, replication slippage (indels), selective pressure, and the molecular clock.
  - Instrumental Method: Whole-genome reference assembly construction and high-throughput re-sequencing.
  - Underlying Algorithm: Dynamic Programming. Derivation of the Needleman-Wunsch global alignment matrix with affine gap penalties ($O(N^2)$ space/time).
  - Industrial Tool: Bio.Align.PairwiseAligner (Biopython C-accelerated module). 
  - Dataset: Orthologous $\alpha$-globin and $\beta$-globin gene sequences across 5 vertebrate species (FASTA).

* Lecture 3: Genetic Polymorphism & Disease (Variant Calling)
  - Biological Mechanism: Single Nucleotide Polymorphisms (SNPs), diploid heterozygosity vs. homozygosity, penetrance, and pathogenic missense mutations causing loss of function.
  - Instrumental Method: Hybridization capture (Targeted Exome Panels) and sequencing read depth accumulation.
  - Underlying Algorithm: Bayesian Read Classification. Calculating posterior genotype probabilities:$$P(G \mid D) \propto P(D \mid G) P(G)$$using base quality scores as error likelihoods.
  - Industrial Tool: freebayes (executed in Colab) + cyvcf2 (Python VCF parser). 
  - Dataset: 5 MB 1000 Genomes Project BAM slice covering the human BRCA1 locus.
  
* Lecture 4: Transcribing the Message (Bulk RNA-Seq Quantification)
  - Biological Mechanism: RNA Polymerase II promoter kinetics, transcription factor initiation, alternative splicing, mRNA turnover rates, and separating biological variance from technical noise.
  - Instrumental Method: Bulk RNA-seq (poly-A selection, reverse transcription to cDNA, PCR amplification, cDNA sequencing).
  - Underlying Algorithm: Negative Binomial Generalized Linear Models (GLMs). Modeling overdispersed RNA read counts and applying empirical Bayes dispersion shrinkage.
  - Industrial Tool: kallisto (pseudoalignment) + PyDESeq2 (Python DESeq2 implementation). 
  - Dataset: Saccharomyces cerevisiae wild-type vs. heat-shock response count matrix (2 MB CSV).

* Lecture 5: Cellular Individuality (Single-Cell Droplet Microfluidics)
  - Biological Mechanism: Tissue heterogeneity, cellular microenvironments, transcriptional bursting, cell lysis kinetics, and mitochondrial leakage during apoptosis.
  - Instrumental Method: Droplet Microfluidics (10x Chromium). Encapsulating single cells into oil droplets with barcoded hydrogel beads and Unique Molecular Identifiers (UMIs).
  - Underlying Algorithm: Sparsity Filtering & Robust Statistics. Median Absolute Deviation (MAD) for dynamic outlier detection and UMI collision correction.
  - Industrial Tool: Scanpy + anndata. 
  - Dataset: 10x Genomics PBMC3k raw unfiltered count matrix (5 MB).
  
* Lecture 6: Cellular State Spaces (Manifolds & Differentiation)
  - Biological Mechanism: Waddington’s epigenetic landscape, continuous cell state transitions, marker gene expression, and immune cell lineage determination.
  - Instrumental Method: Single-cell transcriptomic profiling of complex whole-tissue systems.
  - Underlying Algorithm: Graph-Based Manifold Learning. Spectral decomposition (PCA), $k$-Nearest Neighbors ($k$-NN) graph construction, and Leiden community detection.
  - Industrial Tool: Scanpy (sc.pp.pca, sc.pp.neighbors, sc.tl.umap, sc.tl.leiden). 
  - Dataset: Filtered and normalized PBMC3k AnnData object from Lecture 5.
  
* Lecture 7: Folding the Molecular Machine (Protein 3D Conformation)
  - Biological Mechanism: Ribosomal translation, peptide bond planar resonance, steric hindrance, Ramachandran angle constraints ($\phi, \psi$), hydrophobic collapse, and secondary structure nucleation ($\alpha$-helices, $\beta$-sheets).
  - Instrumental Method: Cryogenic Electron Microscopy (Cryo-EM) and X-ray crystallography (electron density phase solving).
  - Underlying Algorithm: Spatial Coordinate Transformations & Kabsch Algorithm. SVD-based rotation and translation for Root-Mean-Square Deviation (RMSD) superposition:$$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N \Vert{}\mathbf{x}_i - \mathbf{y}_i\Vert{}^2}$$
  - Industrial Tool: Bio.PDB + py3Dmol (WebGL in-notebook rendering). 
  - Dataset: SARS-CoV-2 Main Protease wild-type vs. mutant structures (6LU7.pdb).
  
* Lecture 8: Co-Evolutionary Memory & AI (Protein Language Models)
  - Biological Mechanism: Epistasis, compensatory mutations across geological timescales, sequence conservation landscapes, and predicting the functional impact of novel mutations.
  - Instrumental Method: Deep Mutational Scanning (DMS) libraries coupled with high-throughput sequencing selection assays.
  - Underlying Algorithm: Self-Attention & Masked Language Modeling (MLM). Extracting structural contact maps from transformer attention matrices:$$\text{Attention}(Q, K, V) = \text{softmax}\left(\frac{QK^T}{\sqrt{d_k}}\right)V$$
  - Industrial Tool: Hugging Face transformers + Meta's esm2_t6_8M_UR50D (using Colab T4 GPU). 
  - Dataset: Green Fluorescent Protein (GFP) deep mutational scanning fitness matrix (Sarkisyan et al., 3 MB CSV).
  
  
  
  
  
  
  
  
  
  
  
  
  
