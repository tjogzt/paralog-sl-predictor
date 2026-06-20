# Delta Dependency Prioritizes Paralog-Based Synthetic Lethality Candidates Across Solid Tumor Types

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**English** | **[中文](README_CN.md)**

## Overview

Delta Dependency (DD) is a simple, interpretable, discovery-stage metric for prioritizing paralog-based synthetic lethality (SL) candidates from Cancer Dependency Map (DepMap) CRISPR screening data. DD measures the shift in Chronos gene-effect scores between driver-mutant and wild-type cell lines, computed separately per cancer lineage.

**Key findings:**
- DD achieves AUROC = 0.736 using a single feature in head-to-head comparison with multi-feature ML classifiers on identical test sets
- CPTAC proteomics across 7 cohorts (n=672) reveals protein-level paralog co-variation undetectable at the RNA level
- Therapeutic window analysis nominates ARID1A→ARID1B as the leading selective candidate (TI = 4.13)
- All candidates are computationally nominated and require experimental validation

## Directory Structure

```
paralog_sl_predictor/
├── README.md                  # This file
├── LICENSE                    # MIT License
├── .gitignore                 # Excludes large data files
├── requirements.txt           # Python dependencies
├── Dockerfile                 # Reproducible environment
│
├── manuscript.tex             # Main manuscript (LaTeX)
├── supplementary.tex          # Supplementary material (LaTeX)
├── manuscript.pdf             # Compiled main manuscript
├── supplementary.pdf          # Compiled supplementary
├── cover_letter.md            # Submission cover letter
│
├── config.py                  # Driver genes, known SL pairs, parameters
├── data_loader.py             # DepMap data loading utilities
├── main.py                    # Main analysis pipeline
├── pcs.py                     # DD/PCS computation engine
├── pancancer.py               # Pan-cancer (23 solid tumor types) analysis
├── prism_analysis.py          # PRISM drug sensitivity analysis
├── msi_analysis.py            # MSI stratification analysis
├── mutation_type_analysis.py  # Truncating vs missense analysis
├── protein_features.py        # UniProt protein feature extraction
├── alphafold_analysis.py      # Structural similarity analysis
│
├── R_fig1.R ... R_fig4.R      # Main figure generation (R)
├── R_figS1.R ... R_figS9.R    # Supplementary figure generation (R)
├── R_figS8.R                  # Bootstrap + permutation (10,000 iter)
│
├── data/                      # Input data (large files in .gitignore)
│   ├── README.md              # Download instructions
│   ├── cptac_cache/           # CPTAC protein abundance (7 cohorts)
│   └── ...
│
├── output/                    # Analysis outputs
│   ├── figures/               # All PDF figures (main + supplementary)
│   ├── tables/                # TSV data tables (S1-S6)
│   └── *.csv                  # Intermediate results
│
└── R_package/                 # paralogSL R package (v1.0.0)
    ├── DESCRIPTION
    ├── R/
    ├── man/
    ├── data/
    ├── vignettes/
    └── README.md
```

## Quick Start

### Python Pipeline

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Download DepMap data (see data/README.md)

# 3. Run main analysis
python main.py

# 4. Run pan-cancer analysis
python pancancer.py

# 5. Generate figures (R)
Rscript R_fig1.R
Rscript R_fig2.R
Rscript R_fig3.R
Rscript R_fig4.R
```

### R Package

```r
# Install
devtools::install_github("tjogzt/paralogSL")

# Quick analysis
library(paralogSL)
result <- compute_dd(dep_matrix, driver_gene = "ARID1A",
                     paralog_gene = "ARID1B",
                     mut_lines = mut_ids, wt_lines = wt_ids)
```

## Reproducibility

### Docker

```bash
docker build -t paralog-sl .
docker run -v $(pwd)/data:/app/data paralog-sl python main.py
```

### Random Seeds

All analyses use `set.seed(42)` (R) or `random_state=42` (Python) for reproducibility.

## Data Availability

- **DepMap 26Q1**: https://depmap.org/portal/download/
- **CPTAC**: via cBioPortal API (cached in `data/cptac_cache/`)
- **PRISM**: https://depmap.org/portal/download/
- **Processed tables**: `output/tables/TableS1-S6`

## Citation

```
Zhu T. Delta Dependency Prioritizes Paralog-Based Synthetic Lethality
Candidates Across Solid Tumor Types. (2026).
```

## License

MIT License. See [LICENSE](LICENSE).
