
# Impact of Rifampin on Gene Expression in Mycobacterium tuberculosis H37Rv: Transcriptomics and Rv0784 Molecular Modelling Analysis

This repository contains the complete, integrated computational workflow combining **RNA-Seq transcriptomics** and **all-atom Molecular Dynamics (MD) simulations** to investigate the systemic response of *Mycobacterium tuberculosis* H37Rv to subinhibitory Rifampin stress, with a specific focus on the structural modeling of the target protein **Rv0784**.
---

## Authors

* **Khaled Sabry El-Basha**
* **Adham Srour**
* **Moaaz Ahmed**
* **Sara Saeed**
* **Jana Mahmoud**
* **Habiba Eraky**
* **Alaa Mahmoud**

---

##  Repository Structure

The pipeline is modularized into distinct operational directories for maximum reproducibility:

```text
├── README.md                        # Comprehensive project documentation
├── .gitignore                      # Excludes heavy MD trajectories and binary files
│
├── 01_transcriptomics/              # RNA-Seq Differential Expression & Enrichment (R)
│   ├── 01_deseq2_analysis.R         # Data normalization, QA, and DESeq2 modeling
│   └── 02_functional_enrichment.R   # ORA and GSEA pathway analysis (clusterProfiler)
│
├── 02_molecular_modelling/          # Target & Ligand Preparation (Python)
│   └── zinc_ligand_prep.py          # Automated extraction and 3D optimization via RDKit
│
└── 03_md_simulation/                # GROMACS MD Setup & Trajectory Analysis (Bash & Python)
    ├── mdp_files/                   # Molecular Dynamics Parameter (.mdp) configurations
    ├── 01_gromacs_setup.sh          # System minimization and NVT/NPT equilibration commands
    └── 02_trajectory_analysis.py    # Structural kinetics analysis (MDAnalysis, ProLIF)

```

---

##  Scientific Modules & Methodology

### 1. Transcriptomics Pipeline (`01_transcriptomics/`)

* **Differential Expression:** Implements the `DESeq2` framework in R to model time-course gene expression variance (Control, 30 min, and 60 min) under subinhibitory Rifampin exposure ($1/4 \times \text{MIC}$). Significant DEGs are filtered via $\text{p}_{\text{adj}} < 0.05$ and $|\log_2\text{FC}| > 0.5$.
* **Quality Assurance:** Generates Relative Log Expression (RLE) and Principal Component Analysis (PCA) plots to track inter-replicate variance before and after DESeq2 normalization.
* **Functional Annotation:** Conducts Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA) via `clusterProfiler` leveraging a customized mapping file for *M. tuberculosis* H37Rv.

### 2. Molecular Modelling & Preparation (`02_molecular_modelling/`)

* Standardizes the structural data of the target protein **Rv0784**.
* Automates bulk fetching, valency validation, and 3D coordinate convergence of prospective ligand datasets from the ZINC database via `RDKit`.

### 3. All-Atom Molecular Dynamics Simulation (`03_md_simulation/`)

* **Core Engine:** Powered by **GROMACS 2024.4** leveraging native GPU acceleration.
* **Equilibration:** Employs a dual-phase equilibration protocol: 100 ps NVT (isothermal-isochoric) and 100 ps NPT (isothermal-isobaric) ensembles at 1 bar pressure.
* **Production Run:** Generates a 65 ns production trajectory tracking full atomistic coordinates.
* **Trajectory Analysis:** Quantifies structural compaction, thermodynamic stability, and local atomic flexibility via `MDAnalysis` and `mdakit-sasa`, calculating global 1D/2D RMSD, RMSF, SASA, and Radius of Gyration ($R_g$).
* **Interaction Mapping:** Extracts dynamic Protein-Ligand Interaction Fingerprints (PLIF) and binding timeline tracking via `ProLIF`.

---

##  Installation & Dependencies

To execute or reproduce the workflows, ensure the following environments are configured:

### R Environment

Execute within an active R session to install Bioconductor and CRAN prerequisites:

```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot", "EnhancedVolcano"))
install.packages(c("tidyverse", "pheatmap", "patchwork", "EDASeq"))

```

### Python Environment

Install core cheminformatics and trajectory processing libraries via `pip`:

```bash
pip install rdkit mdanalysis mdakit-sasa prolif matplotlib plotly requests

```

### Molecular Dynamics Engine

Ensure **GROMACS 2024.4** is compiled and sourced in your system path:

```bash
source /usr/local/gromacs/bin/GMXRC

```

---

##  License

This project is open-source and distributed under the terms of the **MIT License**.
