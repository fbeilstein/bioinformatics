:::titlepage
[[top]]
Part 4

[[title]]
Industrial Tools

[[subcaption]]
Bio.PDB · py3Dmol · Kabsch Algorithm · 6LU7
:::

---

# Curriculum

- **The Protein Data Bank (`.pdb`)** — coordinate records (`ATOM` vs `HETATM`), resolution, header metadata, PDB format columns.
- **The SMCRA Hierarchy in `Bio.PDB`** — Structure $\to$ Model $\to$ Chain $\to$ Residue $\to$ Atom, filtering heteroatoms, atom coordinates.
- **The Ramachandran Plot** — peptide bond planar resonance, $\phi/\psi$ dihedral angle extraction with `PPBuilder`, sterically allowed secondary structure basins.
- **Root-Mean-Square Deviation (RMSD)** — metric coordinate difference, why simple coordinate subtraction fails due to arbitrary reference frames.
- **The Kabsch Algorithm** — centroid translation, $3 \times 3$ cross-covariance matrix, SVD decomposition, determinant reflection check, optimal rotation matrix $R$.
- **Structural Superposition in Biopython** — `Bio.PDB.Superimposer`, aligning SARS-CoV-2 Mpro (6LU7) with SARS-CoV-1 Mpro (1UK4).
- **Pairwise Distance Matrices** — $C\alpha$ contact maps, geometric signatures of $\alpha$-helices, parallel and anti-parallel $\beta$-sheets.
- **Temperature Factors (B-factors)** — crystallographic mobility, Debye-Waller equation, flexible surface loops vs rigid catalytic cores.
- **3D Molecular Visualization** — GPU-accelerated WebGL rendering with `py3Dmol` (cartoons, sticks, ligands).

---

# Roadmap

The computational structural biology workflow for atomic coordinates:

:::static-diagram { width=95% center }
{
  "nodes": [
    { "id": "pdb",   "label": "RCSB PDB Download\n6LU7.pdb (SARS-CoV-2)",  "x": 40,   "y": 80,  "width": 160, "height": 60 },
    { "id": "parse", "label": "Bio.PDB Parser\nSMCRA Hierarchy",          "x": 230,  "y": 80,  "width": 150, "height": 60, "style": "fill: #9C27B0; color: #fff;" },
    { "id": "rama",  "label": "Ramachandran Plot\nphi/psi Dihedrals",      "x": 410,  "y": 80,  "width": 160, "height": 60 },
    { "id": "dist",  "label": "C-alpha Distance Map\nSecondary Motifs",    "x": 600,  "y": 80,  "width": 160, "height": 60, "style": "fill: #2196F3; color: #fff;" },
    { "id": "bfact", "label": "B-factor Analysis\nAtomic Mobility",        "x": 790,  "y": 80,  "width": 150, "height": 60, "style": "fill: #FF9800; color: #fff;" },
    { "id": "view",  "label": "py3Dmol WebGL\nCartoon + Inhibitor",        "x": 960,  "y": 80,  "width": 150, "height": 60, "style": "fill: #009688; color: #fff;" },

    { "id": "raw2",  "label": "Target Structure\n1UK4.pdb (SARS-CoV-1)",   "x": 960,  "y": 250, "width": 160, "height": 60 },
    { "id": "pair",  "label": "Extract Paired C-alphas\nMatched Residues", "x": 770,  "y": 250, "width": 160, "height": 60 },
    { "id": "cent",  "label": "Centroid Centering\nTranslate to Origin",   "x": 580,  "y": 250, "width": 160, "height": 60, "style": "fill: #673AB7; color: #fff;" },
    { "id": "svd",   "label": "Kabsch SVD Algorithm\nCross-Covariance",    "x": 390,  "y": 250, "width": 160, "height": 60, "style": "fill: #E91E63; color: #fff;" },
    { "id": "rot",   "label": "Optimal Rotation R\ndet(R) = +1 Check",     "x": 200,  "y": 250, "width": 160, "height": 60, "style": "fill: #4CAF50; color: #fff;" },
    { "id": "rmsd",  "label": "Superposition RMSD\nMeasure Divergence",    "x": 10,   "y": 250, "width": 160, "height": 60, "style": "fill: #4CAF50; color: #fff;" }
  ],
  "links": [
    { "from": "pdb",   "to": "parse" },
    { "from": "parse", "to": "rama" },
    { "from": "rama",  "to": "dist" },
    { "from": "dist",  "to": "bfact" },
    { "from": "bfact", "to": "view" },
    { "from": "view",  "to": "raw2" },
    { "from": "raw2",  "to": "pair" },
    { "from": "pair",  "to": "cent" },
    { "from": "cent",  "to": "svd" },
    { "from": "svd",   "to": "rot" },
    { "from": "rot",   "to": "rmsd" }
  ]
}
:::

> In this lecture, we move from 1D sequence and 2D manifolds to 3D Cartesian space ($x, y, z$), solving the geometry of protein folding, conformational alignment, and ligand recognition.

---

# The Structural Biology Toolbox

Analyzing 3D macromolecules requires libraries that understand atomic coordinates, topology, and spatial transformations:

:::matrix { .cards cols="1fr 1fr 1fr 1fr" gap="20px" margin-top: 20px; }
[[0,0]]{ .purple }
### Bio.PDB
Biopython's structural module: parses PDB and mmCIF files into the hierarchical SMCRA object model, extracts backbone dihedral angles, and superimposes structures.

[[0,1]]{ .blue }
### NumPy & SciPy
Coordinate array manipulation, 3D SVD decomposition for the Kabsch algorithm, and high-performance pairwise Euclidean distance matrices (`cdist`).

[[0,2]]{ .teal }
### py3Dmol
WebGL GPU-accelerated molecular visualizer for Jupyter and web notebooks: renders cartoons, sticks, surfaces, and bound inhibitors in 3D.

[[0,3]]{ .orange }
### Matplotlib
2D projection of conformational state spaces: Ramachandran plots, contact heatmaps, and B-factor sequence profiles.
:::

```bash
# Installation in Google Colab or virtual environment
pip install -q biopython py3Dmol numpy pandas scipy matplotlib
```

---

# The Data — SARS-CoV-2 Main Protease (6LU7)

Our model system is the primary catalytic drug target of the COVID-19 pandemic:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Biological Target: $M^{\text{pro}}$
- **Main Protease ($M^{\text{pro}}$ / $3\text{CL}^{\text{pro}}$)**:
  - Essential homodimeric cysteine protease that cleaves the viral polyproteins pp1a and pp1ab at 11 conserved cleavage sites.
  - Crucial for viral replication; no human homolog with identical cleavage specificity, making it an ideal drug target.
- **PDB ID: `6LU7`** (Jin et al., *Nature* 2020):
  - Solved at **2.16 Å resolution** by X-ray crystallography.
  - Structure includes a covalently bound synthetic peptide-like inhibitor (**N3**).

[[0,1]]
### Fetching from the RCSB PDB
```python
import urllib.request

pdb_id = "6LU7"
url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
# RCSB requires a standard User-Agent header
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req) as resp, open(f"{pdb_id}.pdb", "wb") as f:
    f.write(resp.read())
```
- Downloaded size: **~239 KB** of plain-text ASCII coordinate data.
:::

---

# Inside the PDB Format: The Coordinate Line

The PDB format is a strict 80-column fixed-width legacy text standard:

```text
ATOM     42  CA  GLY A   6      -9.845  14.231   4.108  1.00 24.12           C
^^^^    ^^^  ^^  ^^^ ^ ^^^      ^^^^^^  ^^^^^^  ^^^^^^  ^^^^ ^^^^^           ^
Record   ID Name Res Ch ResNum       X       Y       Z   Occ  B-factor     Elem
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Field Definitions
- **Record**: `ATOM` (standard polypeptide residues) vs `HETATM` (ligands, ions, water).
- **Atom Name**: `CA` (Alpha Carbon), `N`, `C`, `O`, `CB`, etc.
- **Residue & Chain**: `GLY` on Chain `A`, sequence position 6.
- **Coordinates ($X, Y, Z$)**: Orthogonal Cartesian position in Angstroms (1 Å = $10^{-10}\text{ m} = 0.1\text{ nm}$).

[[0,1]]
### Physical Parameters
- **Occupancy (`1.00`)**: Fraction of crystal unit cells in which this atom occupies this specific location.
- **B-factor (`24.12`)**: Temperature factor representing isotropic thermal displacement and crystal disorder ($\text{Å}^2$).
- **Element (`C`)**: Chemical element symbol.
:::

---

# Tool 1 — The SMCRA Hierarchy in `Bio.PDB`

Biopython organizes structural macromolecular data using the **SMCRA** object architecture:

```text
Structure (6LU7)
  └── Model (0)                 <-- Model 0 (X-ray has 1 model; NMR has ensembles)
        ├── Chain ('A')         <-- Main protease enzyme polymer (306 residues)
        │     ├── Residue (1)   <-- Amino acid (Serine 1)
        │     │     ├── Atom ('N')   <-- Coordinate: [-12.2, 16.4, 3.1]
        │     │     ├── Atom ('CA')  <-- Coordinate: [-11.1, 15.2, 3.8]
        │     │     └── Atom ('C')   <-- Coordinate: [-10.4, 14.9, 2.5]
        │     └── Residue (2)   <-- Next amino acid...
        └── Chain ('C')         <-- Bound covalent inhibitor peptide (N3)
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Loading the Structure
```python
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)  # suppresses PDBConstructionWarning
structure = parser.get_structure("6LU7", "6LU7.pdb")
model = structure[0]
chain_A = model['A']
```

[[0,1]]
### Standard Residues vs Heteroatoms
Residues store a tuple ID: `(hetflag, resseq, icode)`.
- `res.id[0] == ' '`: Standard amino acid.
- `res.id[0] == 'W'`: Water molecule (`HOH`).
- `res.id[0].startswith('H_')`: Heteroatom / ligand.
:::

---

# Accessing Atoms & Coordinates in Biopython

Once inside a residue, atoms are dictionary-like items providing direct geometric queries:

```python
# Extract all standard amino acid residues in Chain A
residues = [res for res in chain_A if res.id[0] == ' ']
print(f"Chain A contains {len(residues)} standard amino acid residues.")

# Access the Alpha Carbon of the first residue
first_res = residues[0]
ca_atom = first_res['CA']

# 1. NumPy Cartesian coordinates
coord = ca_atom.get_coord()          # np.ndarray shape (3,), e.g. [-11.1, 15.2, 3.8]
# 2. Temperature factor
bfactor = ca_atom.get_bfactor()      # float, e.g. 24.12
# 3. Chemical element
element = ca_atom.element            # str, e.g. 'C'
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Vector Math
Biopython overloads the minus operator (`-`) for atoms:
```python
distance = ca_atom1 - ca_atom2
```
Computes the Euclidean distance in Angstroms directly!

[[0,1]]
### Extracting C$\alpha$ Traces
```python
ca_coords = np.array([
    res['CA'].get_coord() 
    for res in residues if 'CA' in res
])
print(ca_coords.shape)  # (306, 3)
```
:::

---
# Tool 2 — The Ramachandran Plot

Why do proteins fold into repeating secondary structures like $\alpha$-helices and $\beta$-sheets?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Rigid Peptide Bond
- The peptide bond ($C - N$) possesses $\sim 40\\%$ double-bond character due to resonance with the carbonyl oxygen:
  $$\mathrm{O=C - NH \longleftrightarrow O^- - C = N^+H}$$
- The 6 atoms of the peptide group lie in a flat, planar unit.
- The torsion angle $\omega$ is locked in the *trans* conformation ($\omega \approx 180^\circ$).

[[0,1]]
### The Two Degrees of Freedom
Each amino acid residue can rotate about only **two backbone single bonds**:
1. **$\phi$ (phi)**: angle of rotation about the $N - C\alpha$ bond.
2. **$\psi$ (psi)**: angle of rotation about the $C\alpha - C$ bond.
:::

> Most ($\phi, \psi$) combinations are sterically impossible: side-chain atoms clash violently with backbone carbonyl oxygens and amide hydrogens. The few allowed conformational zones form the **Ramachandran Plot**.

---

# Extracting Dihedral Angles with `PPBuilder`

Biopython provides `PPBuilder` (Polypeptide Builder) to navigate backbone chains and compute $(\phi, \psi)$ angles:

```python
from Bio.PDB.Polypeptide import PPBuilder
import math

builder = PPBuilder()
phi_psi_angles = []

# Build connected polypeptide segments (breaks at missing loop gaps)
for poly in builder.build_peptides(chain_A):
    # Returns a list of (phi, psi) tuples in RADIANS
    for phi, psi in poly.get_phi_psi_list():
        # First residue lacks phi; last residue lacks psi
        if phi is not None and psi is not None:
            phi_psi_angles.append((math.degrees(phi), math.degrees(psi)))
```

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Terminal Residues
- The N-terminal residue has no preceding carbon $\to \phi$ is undefined (`None`).
- The C-terminal residue has no succeeding nitrogen $\to \psi$ is undefined (`None`).
- Filter `None` before mathematical analysis or plotting!

[[0,1]]
### Coordinate Conversion
- `get_phi_psi_list()` returns angles in **radians** ($-\pi \to +\pi$).
- Standard Ramachandran plots are indexed in **degrees** ($-180^\circ \to +180^\circ$).
- Multiply by $180 / \pi$ using `math.degrees()`.
:::

---

# Visualizing the 6LU7 Ramachandran Distribution

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
![Ramachandran Plot](lecture_07_folding_the_molecular_machine/images/ramachandran_plot.png){width=95% center}

[[0,1]]
Plotting all $(\phi, \psi)$ pairs of SARS-CoV-2 Main Protease (6LU7 Chain A):

```python
phi_deg, psi_deg = zip(*phi_psi_angles)
plt.scatter(phi_deg, psi_deg, s=15, alpha=0.6)
plt.xlim(-180, 180); plt.ylim(-180, 180)
plt.xlabel(r'$\phi$ (degrees)')
plt.ylabel(r'$\psi$ (degrees)')
```

- Nearly all non-glycine residues cluster tightly in the **$\beta$-sheet** (upper-left) and **right-handed $\alpha$-helix** (lower-left) basins.
- **Glycine** (orange triangles) uniquely populates the upper-right left-handed $\alpha$ region.
:::

---

# Reading the Ramachandran Conformational Basins

The distribution of points clusters into distinct thermodynamic free-energy valleys:

:::matrix { cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]
### 1. $\beta$-Sheet Region
**Upper Left Quadrant**
$$\phi \approx -120^\circ \text{ to } -140^\circ$$
$$\psi \approx +130^\circ \text{ to } +160^\circ$$
- Extended polypeptide chain.
- Accommodates both anti-parallel and parallel $\beta$-pleated sheets without steric clash.

[[0,1]]
### 2. Right-Handed $\alpha$-Helix
**Lower Left Quadrant**
$$\phi \approx -60^\circ \text{ to } -70^\circ$$
$$\psi \approx -40^\circ \text{ to } -50^\circ$$
- Compact spiral forming regular $i \to i+4$ backbone hydrogen bonds.
- The dominant secondary structure in 6LU7.

[[0,2]]
### 3. Left-Handed $\alpha$-Helix
**Upper Right Quadrant**
$$\phi \approx +60^\circ, \quad \psi \approx +40^\circ$$
- Sterically forbidden for 19 of the 20 amino acids because side chains collide with the backbone!
- **Glycine Exception**: Having only a single hydrogen atom ($-\text{H}$), Glycine easily populates this quadrant.
:::

---

# Tool 3 — Structural Superposition & RMSD

How do we quantitatively compare two 3D protein structures (e.g. wild-type vs mutant, or SARS-CoV-2 vs SARS-CoV-1)?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Why Raw Coordinate Subtraction Fails
- Suppose protein $P$ and protein $Q$ have identical 3D shapes, but:
  - $P$ was crystallized in Paris (origin at center of unit cell).
  - $Q$ was crystallized in London (rotated by $45^\circ$ and translated by $20\text{ Å}$).
- Direct point-to-point subtraction gives huge errors:
  $$\frac{1}{N} \sum_{i=1}^N \lVert \mathbf{p}\_i - \mathbf{q}\_i\rVert^2 \gg 0$$
- Coordinate values depend on an **arbitrary global reference frame**!

[[0,1]]
### The Root-Mean-Square Deviation (RMSD)
RMSD measures true structural divergence **after** finding the optimal rigid-body translation $\mathbf{t}$ and rotation $R$:

$$\text{RMSD} = \min_{R, \mathbf{t}} \sqrt{\frac{1}{N} \sum_{i=1}^N \lVert (R \mathbf{p}\_i + \mathbf{t}) - \mathbf{q}\_i\rVert^2}$$

$$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N \lVert \mathbf{p}\_i' - \mathbf{q}\_i\rVert^2}$$
:::

---
# The Kabsch Algorithm: Step-by-Step Derivation

Wolfgang Kabsch (1976) solved the optimal rotation problem using Singular Value Decomposition:

### Step 1: Centering
Translate both coordinate sets so their centroids coincide at the origin $(0, 0, 0)$:
$$\mathbf{c}\_P = \frac{1}{N}\sum_{i=1}^N \mathbf{p}\_i, \qquad \mathbf{c}\_Q = \frac{1}{N}\sum_{i=1}^N \mathbf{q}\_i$$
$$P_c = P - \mathbf{c}\_P, \qquad Q_c = Q - \mathbf{c}\_Q$$
This eliminates the translation vector $\mathbf{t}$ entirely!

### Step 2: The Cross-Covariance Matrix
Compute the $3 \times 3$ cross-dispersion matrix:
$$C = P_c^T Q_c \in \mathbb{R}^{3 \times 3}$$
- Entry $C_{ab} = \sum_{i=1}^N P_{i, a} Q_{i, b}$.
- Encodes how the $x, y, z$ axes of structure $P$ correlate with the $x, y, z$ axes of structure $Q$.

---

# The Kabsch SVD & Determinant Check

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Step 3: SVD Decomposition
Compute the Singular Value Decomposition of the $3 \times 3$ matrix $C$:
$$C = U S V^T$$
- $U \in \mathbb{R}^{3 \times 3}$: orthogonal left singular vectors.
- $V \in \mathbb{R}^{3 \times 3}$: orthogonal right singular vectors.

[[0,1]]
### Step 4: The Optimal Rotation ($R$)
With $d = \text{sign}\left(\det(V U^T)\right)$ and $D = \begin{bmatrix} 1 & 0 & 0 \\\\ 0 & 1 & 0 \\\\ 0 & 0 & d \end{bmatrix}$:
- **Column-vector convention** (notebook Part 5 / Task 1):
  $$R_{\text{col}} = V D U^T \quad \Longrightarrow \quad P_{\text{aligned}} = (R_{\text{col}} P_c^T)^T = P_c R_{\text{col}}^T$$
- **Row-vector convention** (direct right-multiply):
  $$R_{\text{row}} = U D V^T = R_{\text{col}}^T \quad \Longrightarrow \quad P_{\text{aligned}} = P_c R_{\text{row}}$$
:::

> **⚠ Two Traps**: (1) **Chirality**: if $\det(V U^T) = -1$, failing to flip the last column of $V$ (i.e. $d = -1$) produces a mirror reflection inverting L- into D-amino acids. (2) **Transpose**: if you build $R = V D U^T$, you **must** apply `Pc @ R.T` (not `Pc @ R`, which rotates backwards!).

---

# Superposition in Biopython: `Superimposer`

Biopython provides an industrial implementation of the Kabsch algorithm in `Bio.PDB.Superimposer`:

```python
from Bio.PDB import Superimposer

# Align SARS-CoV-2 (6LU7) with SARS-CoV-1 (1UK4)
atoms_1 = [res['CA'] for res in chain_6LU7 if res.id[0] == ' ' and 'CA' in res]
atoms_2 = [res['CA'] for res in chain_1UK4 if res.id[0] == ' ' and 'CA' in res]

# Notebook Part 5 truncates by index for simplicity (real pipelines match by res.id[1] / sequence alignment!)
min_len = min(len(atoms_1), len(atoms_2))
paired_1, paired_2 = atoms_1[:min_len], atoms_2[:min_len]

super_imposer = Superimposer()
super_imposer.set_atoms(paired_1, paired_2)  # Computes centroids, C, SVD, and R
super_imposer.apply(structure_1UK4.get_atoms())
print(f"Superimposed {min_len} C-alpha atoms · RMSD: {super_imposer.rms:.3f} Å")
```

> **⚠ Pairing Pitfall**: `6LU7` spans residues `1–306` while `1UK4` starts at residue `3` (`3–303`). Naive `[:min_len]` slicing shifts every residue by 2 positions (`6.21 Å` RMSD in notebook Part 5); matching shared residue IDs (`res.id[1]`) yields the true backbone divergence (**`1.35 Å`**).

---

# Visualizing the Superposition: 6LU7 vs. 1UK4

Aligning the 2020 pandemic SARS-CoV-2 protease against the 2003 SARS-CoV-1 protease:

![Kabsch Superposition](lecture_07_folding_the_molecular_machine/images/kabsch_superposition.png){width=80% center}

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### Before Superposition (Left)
- The raw PDB files sit in completely different Cartesian coordinate spaces.
- Point-to-point distance is over **$100\text{ Å}$** purely due to crystallographic origin choices!

[[0,1]]
### After Kabsch Alignment (Right)
- Matching shared residue numbers (`3–303`, 301 $C\alpha$ atoms) brings the backbones into close coincidence.
- **$\text{RMSD} = 1.35\text{ Å}$**: confirms that the catalytic fold is highly conserved between SARS-CoV-1 and SARS-CoV-2!
:::

---

# Tool 4 — 2D Distance Matrices (Contact Maps)

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
![Distance Matrix](lecture_07_folding_the_molecular_machine/images/distance_matrix_patterns.png){width=95% center}

[[0,1]]
A 2D pairwise distance matrix maps the internal Euclidean distance between every pair of $C\alpha$ atoms:

$$D_{ij} = \lVert \mathbf{x}\_i - \mathbf{x}\_j\rVert = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}$$

```python
from scipy.spatial.distance import cdist
ca_coords = np.array([
    res['CA'].get_coord()
    for res in chain_A
    if res.id[0] == ' ' and 'CA' in res
])
dist_matrix = cdist(ca_coords, ca_coords, metric='euclidean')
```
:::

---
# Decoding Secondary Motifs on the Distance Map

The 2D distance matrix is completely **invariant to 3D rotation and translation**:

### 1. Main Diagonal Band
**$\alpha$-Helices**
- Consecutive residues ($i, i+1$) are always $3.8\text{ Å}$ apart.
- In an $\alpha$-helix, hydrogen bonding loops residues $i$ and $i+4$ within $5 - 6\text{ Å}$.
- Forms a prominent **thick dark band** along the main diagonal.

### 2. Off-Diagonal Parallel Lines
**Parallel $\beta$-Sheets**
- Two strands run in the same direction ($N \to C$).
- Residue $i$ contacts $j$, $i+1$ contacts $j+1$.
- Generates a narrow stripe running **parallel** to the main diagonal at offset $|i - j|$.

### 3. Off-Diagonal Orthogonal Lines
**Anti-Parallel $\beta$-Sheets**
- Two strands run in opposite directions ($N \to C$ alongside $C \to N$).
- Residue $i$ contacts $j$, but $i+1$ contacts $j-1$.
- Generates a narrow stripe running **perpendicular** to the main diagonal!

---

# Tool 5 — Temperature Factors (B-factors)

What does the B-factor column in the PDB coordinate record tell us about protein physics?

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Debye-Waller Equation
X-ray diffraction intensities attenuate due to thermal vibration and lattice disorder:
$$B = 8\pi^2 \langle u^2\rangle$$
- $\langle u^2\rangle$: mean-square spatial displacement of the atom around its equilibrium coordinate ($\text{Å}^2$).
- For $B = 79\text{ Å}^2$:
  $$\langle u^2\rangle = \frac{79}{8\pi^2} \approx 1.0\text{ Å}^2 \implies \sqrt{\langle u^2\rangle} \approx 1.0\text{ Å}$$
  The atom vibrates with an average amplitude of $1\text{ Å}$!

[[0,1]]
### Structural Interpretation
- **Low B-factor ($< 25\text{ Å}^2$)**: Rigid, tightly packed hydrophobic core and catalytic active site (e.g. catalytic dyad His41–Cys145).
- **High B-factor ($> 60\text{ Å}^2$)**: Highly flexible surface loops, solvent-exposed termini, or disordered regions.
- **Crystal Averaging**: X-ray diffraction sums electron density coherently across millions of periodic crystal unit cells — mobile loops adopting different conformations cancel out destructively, smearing their electron density below noise (often leaving missing residues in the PDB).
:::

---
# Visualizing 6LU7 Crystallographic Mobility

Plotting $C\alpha$ B-factors along the primary sequence of SARS-CoV-2 Main Protease:

![B-factor Profile](lecture_07_folding_the_molecular_machine/images/bfactor_mobility.png){width=71% center}

```python
std_res = [r for r in chain_A if r.id[0] == ' ' and 'CA' in r]
b_factors = [r['CA'].get_bfactor() for r in std_res]
res_nums  = [r.id[1] for r in std_res]
plt.plot(res_nums, b_factors, color='#9467bd')
plt.axhline(np.mean(b_factors), color='red', linestyle='--')
plt.xlabel('Residue Number'); plt.ylabel(r'B-factor ($\AA^2$)')
```
:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
- **Deep Valleys**: Catalytic dyad and core $\beta$-barrel domains are held rigid to preserve enzyme kinetics.

[[0,1]]
- **Sharp Peaks**: Mobile surface loops that flex during substrate binding and catalytic turnover.
:::

---

# Tool 6 — 3D Interactive Visualization with `py3Dmol`

`py3Dmol` embeds WebGL GPU-accelerated molecular graphics directly inside Jupyter notebooks:

```python
import py3Dmol

# Load PDB text data into the viewer
with open("6LU7.pdb", "r") as f:
    pdb_data = f.read()

view = py3Dmol.view(width=800, height=500)
view.addModel(pdb_data, "pdb")

# 1. Render protein backbone as cartoon colored by secondary structure
view.setStyle({'cartoon': {'color': 'spectrum'}})

# 2. Render non-standard heteroatoms (inhibitor N3) as green sticks
view.setStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon'}})

# 3. Center and display
view.zoomTo()
view.show()
```

:::matrix { cols="1fr 1fr 1fr" gap="15px" }
[[0,0]]
**Representation Styles**
- `'cartoon'`: secondary ribbons
- `'stick'`: bonds & side chains
- `'sphere'`: van der Waals spheres

[[0,1]]
**Selection Selectors**
- `{'chain': 'A'}`: specific polymer
- `{'resi': 145}`: residue sequence
- `{'hetflag': True}`: ligands

[[0,2]]
**Color Schemes**
- `'spectrum'`: rainbow $N \to C$
- `'greenCarbon'`: green carbons
- `'b': {'gradient': 'rwb'}`: by B-factor
:::

---

# Homework Task 1 — Kabsch Algorithm from Scratch

In Homework Task 1, you implement the complete Kabsch algorithm `kabsch_rmsd(P, Q)` in pure NumPy:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### The Required Mathematical Steps
1. **Centering**:
   Calculate centroids $\mathbf{c}\_P = \text{mean}(P, \text{axis}=0)$ and $\mathbf{c}\_Q = \text{mean}(Q, \text{axis}=0)$.
   Center both: $P_c = P - \mathbf{c}\_P$ and $Q_c = Q - \mathbf{c}\_Q$.
2. **Covariance Matrix**:
   $$C = P_c^T Q_c \in \mathbb{R}^{3 \times 3}$$
   Use `np.dot(Pc.T, Qc)`.
3. **SVD Decomposition**:
   `U, S, Vt = np.linalg.svd(C)`
   Notice `np.linalg.svd` returns $V^T$ already transposed!

[[0,1]]
### Reflection Check, Convention & RMSD
4. **Determinant & Rotation**:
   - Build $R = V U^T$ (with $V = (V^T)^T$). If $\det(R) < 0$, flip the sign of the last column of $V$ (equivalently `Vt[-1, :] *= -1`) and recompute $R$.
   - **⚠ Transpose trap**: with $R = V U^T$, rotate row-matrix $P_c$ via `P_rot = Pc @ R.T` (not `Pc @ R`!).
5. **RMSD (⚠ Divide by $N$, not $3N$)**:
   $$\text{RMSD} = \sqrt{\frac{1}{N}\sum_{i=1}^N \lVert \mathbf{p}\_{\text{rot}, i} - \mathbf{q}\_{c, i}\rVert^2}$$
   Sum squared diffs along `axis=1` before taking the mean across the $N$ atoms (`np.mean(diff**2)` divides by $3N$ and is off by $\sqrt{3}$!).
:::

> **Self-Check**: Your `kabsch_rmsd(P_coords, Q_coords)` on the notebook's Part 5 arrays must match `super_imposer.rms` to three decimal places!

---

# Homework Task 2 — Distance Matrix Calculation

In Homework Task 2, you compute and visualize the 2D pairwise distance matrix for 6LU7 Chain A:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Computing Pairwise Distances
- Extract the 3D coordinate array of all $C\alpha$ atoms in Chain A (shape: $N \times 3$).
- Use SciPy's spatial distance module:
  ```python
  from scipy.spatial.distance import cdist
  dist_matrix = cdist(coords, coords, metric='euclidean')
  ```
- Result is an $N \times N$ symmetric matrix with zeros along the main diagonal.

[[0,1]]
### Plotting with Matplotlib
- Use `plt.imshow(dist_matrix, origin='lower')`.
- Choose a perceptually uniform colormap (e.g. `'viridis'`, `'plasma'`, or inverted `'viridis_r'`).
- Add a colorbar with units in Angstroms (Å).
- **Interpretation Challenge**:
  - Locate the thick band along the diagonal representing $\alpha$-helices.
  - Locate off-diagonal stripes and determine whether they run parallel or anti-parallel to the diagonal.
:::

---

# Homework Task 3 — B-Factor Profile Analysis

In Homework Task 3, you analyze the crystallographic mobility of SARS-CoV-2 Main Protease:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Extraction Logic
1. Iterate through standard residues of Chain A and extract `res['CA'].get_bfactor()`.
2. Generate the sequential 1-based residue index array ($1, 2, \dots, N$).
3. Plot Residue Index vs B-factor using a clean line graph.
4. Add reference gridlines and axis labels ($\text{Å}^2$).

[[0,1]]
### Locating the Maximum Flexibility Peak
- Use NumPy's `np.argmax()` to find the index of the highest B-factor value:
  ```python
  max_idx = np.argmax(b_factors)
  ```
- Remember Python uses 0-based indexing: **add 1** to report the biological residue number!
- Check where this residue sits in the 3D structure: is it part of a rigid secondary core, or a flexible surface loop?
:::

---

# Homework Task 4 — Visualizing the Bound Ligand (N3)

In Homework Task 4, you isolate and render the covalent inhibitor N3 in 6LU7:

:::matrix { cols="1.1fr 1fr" gap="30px" }
[[0,0]]
### Why Is N3 Stored as Chain C?
- Most simple ligands (e.g. ATP, zinc ions) are stored as single `HETATM` records within Chain A.
- However, inhibitor **N3** is a synthetic **peptidomimetic compound** constructed from multiple non-natural amino acid fragments.
- The PDB annotators assigned N3 its own separate chain: **Chain C**.

[[0,1]]
### py3Dmol Visualization Recipe
1. Count total atoms in Chain C: `model['C'].get_atoms()` returns a **generator** (no `len()`), so convert first via `len(list(model['C'].get_atoms()))`.
2. Configure `py3Dmol`:
   - Hide or dim the protein (`view.setStyle({'chain': 'A'}, {'line': {'opacity': 0.3}})` or omit it).
   - Style the ligand: `view.setStyle({'chain': 'C'}, {'stick': {'colorscheme': 'greenCarbon'}})`.
   - Call `view.zoomTo({'chain': 'C'})` then `view.show()` to focus directly on N3.
:::

---

# Homework Task 5 — Conceptual Review & Guidance

Task 5 tests your understanding of structural biophysics and coordinate geometry:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### 1. Why Kabsch for RMSD?
- How does SVD isolate *internal* structural divergence from arbitrary reference frame translations and rotations?

### 2. Ramachandran Quadrants & Glycine
- What steric clashes prevent non-glycine residues from adopting positive $\phi$ and positive $\psi$ angles? Why can Glycine populate this region?

### 3. Contact Map Topology
- Why does an anti-parallel $\beta$-sheet generate a line perpendicular to the main diagonal, while a parallel $\beta$-sheet runs parallel?

[[0,1]]
### 4. Invisible Flexible Loops in X-ray
- Why does destructive X-ray interference from non-periodic, mobile loops cause regions of high B-factor to disappear from electron density maps?

### 5. SVD in PCA vs. SVD in Kabsch
- Compare the matrix being decomposed in PCA ($X^T X$, covariance of one dataset) vs. Kabsch ($P_c^T Q_c$, cross-covariance of two paired frames): what geometric problem does each SVD solve?
:::

---

# Putting It All Together: Complete Structural Pipeline

```python
import numpy as np, py3Dmol
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import PPBuilder
# 1. Parse PDB structure into SMCRA hierarchy
parser = PDBParser(QUIET=True)
s_cov2 = parser.get_structure("6LU7", "6LU7.pdb")
chain_A = s_cov2[0]['A']
# 2. Extract backbone dihedral angles (Ramachandran validation)
phi_psi = []
for poly in PPBuilder().build_peptides(chain_A):
    for phi, psi in poly.get_phi_psi_list():
        if phi is not None and psi is not None:
            phi_psi.append((np.degrees(phi), np.degrees(psi)))
# 3. Structural Superposition against SARS-CoV-1 (1UK4) — matched by residue ID
s_cov1 = parser.get_structure("1UK4", "1UK4.pdb")
c2 = {r.id[1]: r['CA'] for r in chain_A if r.id[0] == ' ' and 'CA' in r}
c1 = {r.id[1]: r['CA'] for r in s_cov1[0]['A'] if r.id[0] == ' ' and 'CA' in r}
common = sorted(set(c2) & set(c1))
sup = Superimposer()
sup.set_atoms([c2[i] for i in common], [c1[i] for i in common])
sup.apply(s_cov1.get_atoms())
print(f"Matched-residue RMSD: {sup.rms:.3f} Å across {len(common)} C-alpha atoms.")
# 4. Interactive 3D visualization of the protease and inhibitor
view = py3Dmol.view(width=800, height=450)
with open("6LU7.pdb") as f: view.addModel(f.read(), "pdb")
view.setStyle({'cartoon': {'color': 'spectrum'}})
view.setStyle({'chain': 'C'}, {'stick': {'colorscheme': 'greenCarbon'}})
view.zoomTo()
```

---

# Summary

:::matrix { .summary cols="1fr 1fr 1fr" gap="20px" }
[[0,0]]{ .purple }
### PDB & The Hierarchy
#### Structure
- SMCRA: Struct $\to$ Model $\to$ Chain $\to$ Res $\to$ Atom
- `ATOM` vs `HETATM`
- Rigid peptide planar bond
- $\phi/\psi$ dihedrals in `PPBuilder`
- Ramachandran allowed basins
- Glycine steric exception

[[0,1]]{ .green }
### Kabsch & RMSD
#### Geometry
- Coordinate subtraction fails
- Centroid origin translation
- Cross-covariance $C = P_c^T Q_c$
- SVD: $C = U S V^T$
- Reflection check: $\det(R) = +1$
- Optimal rotation: $R = U \operatorname{diag}(1,1,d) V^T$

[[0,2]]{ .blue }
### Biophysics & WebGL
#### Visualization
- Contact maps via `cdist`
- Diagonal = helices; off = sheets
- B-factors: $8\pi^2 \langle u^2\rangle$
- Core rigidity vs loop mobility
- `py3Dmol`: cartoons & sticks
- Inhibitor N3 in Chain C
:::

---

# Next Steps — Lecture 8

In **Lecture 8: Co-Evolutionary Memory & AI (Protein Language Models)**, we transition from rigid crystal coordinates to learned sequence representations and attention maps:

:::matrix { cols="1fr 1fr" gap="30px" }
[[0,0]]
### What We Accomplished in Lecture 7
- Parsed 3D atomic coordinates from PDB records using Biopython's SMCRA hierarchy.
- Evaluated backbone conformations against stereochemical Ramachandran limits.
- Mathematically derived and implemented the Kabsch algorithm for optimal 3D superposition.
- Decoded secondary structure topologies from pairwise distance contact maps.
- Profiled crystallographic flexibility using atomic temperature factors (B-factors).

[[0,1]]
### What Comes Next (Lecture 8)
- **Deep Mutational Scanning (DMS)**: High-throughput experimental fitness landscapes (GFP mutation matrix).
- **Masked Language Modeling (MLM)**: Training transformers on billions of evolutionary sequences.
- **The Self-Attention Mechanism**: $\text{Attention}(Q, K, V) = \text{softmax}\left(\frac{QK^T}{\sqrt{d_k}}\right)V$.
- **Evolutionary Contact Maps**: How transformer attention heads spontaneously learn 3D contact geometry from 1D sequence co-evolution!
- **Meta's ESM-2 & ESMFold**: Zero-shot variant effect prediction and AI structure prediction.
:::
