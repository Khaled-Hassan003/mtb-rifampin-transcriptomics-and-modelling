# Impact of Rifampin on Gene Expression in Mycobacterium tuberculosis H37Rv: transcriptomics and Rv0784 Molecular Modelling Analysis

This repository contains the complete, integrated computational workflow combining **RNA-Seq transcriptomics** and **all-atom Molecular Dynamics (MD) simulations** to investigate the systemic response of *Mycobacterium tuberculosis* H37Rv to subinhibitory Rifampin stress, with a specific focus on the structural modeling of the target protein **Rv0784**.

---

## 📂 Repository Architecture

The pipeline is modularized into distinct operational directories for maximum reproducibility:

```microbial-md-pipeline/
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
