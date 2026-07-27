# Delta Dependency Prioritizes Candidate Paralog Dependencies Across Solid Tumor Types

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21502030.svg)](https://doi.org/10.5281/zenodo.21502030)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**English** | **[中文](README_CN.md)**

## Overview

Delta Dependency (DD) is a simple, interpretable, discovery-stage metric for prioritizing paralog-based synthetic lethality (SL) candidates from Cancer Dependency Map (DepMap) CRISPR screening data. DD measures the shift in Chronos gene-effect scores between driver-mutant and wild-type cell lines, computed separately per cancer lineage.

**Key findings (see manuscript for full context):**
- **Gold-standard evaluation (evidence-tiered, citation-verified):** AUROC = 0.676 on the full 12-pair curated set; 0.725 excluding two DepMap-derived pairs; 0.774 on pre-DepMap evidence only; unchanged (0.676) under direction-strict relabelling. Both lineage-evaluable positives on the Tier A∪B external benchmark rank above all unlabeled controls (AUROC = 1.000; 2 of 5 benchmark pairs lineage-evaluable)
- **DD vs. published methods (contextual reference, CV3-analogous framework):** AUROC = 0.676 without training; best published CV3 result 0.790 (SLMGAE) — analogous, not identical, evaluation frameworks
- **Head-to-head (identical 72-pair test set, 6 positives):** the interpretable composite score (0.831) matches the best of four multi-feature classifiers under leave-one-pair-out CV (SVM-RBF 0.843, RF 0.722, SVM-Linear 0.114, LR 0.136; DD alone 0.566) — small-n results reported with an explicit power caveat
- DD + ≥30% sequence identity filter → AUROC = 1.000

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
└── R_package/                 # paralogSL R package (v1.1.1)
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

R figure dependencies (install once):

```r
install.packages(c("ggplot2", "dplyr", "readr", "cowplot",
                   "svglite", "ragg", "jsonlite", "tidyr"))
BiocManager::install(c("Biostrings", "pwalign"))  # FigS9 NW alignment
```

### R Package

Standalone R package repository: **[tjogzt/paralogSL](https://github.com/tjogzt/paralogSL)** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21502114.svg)](https://doi.org/10.5281/zenodo.21502114)

```r
# Install
devtools::install_github("tjogzt/paralogSL")

# Quick analysis
library(paralogSL)
result <- compute_dd(dep_matrix, driver_gene = "ARID1A",
                     paralog_gene = "ARID1B",
                     mut_lines = mut_ids, wt_lines = wt_ids)
```
```

## Reproducibility

### Docker

```bash
docker build -t paralog-sl .
docker run -v $(pwd)/data:/app/data paralog-sl python main.py
```

### Random Seeds

All analyses use `set.seed(42)` (R) or `random_state=42` (Python) for reproducibility.

### Verification checklist

Every quantitative claim in `manuscript.tex` is recomputed from raw DepMap
data by three audit scripts, each ending with an automatic claims check that
exits non-zero on mismatch. A fourth script,
`audit_manuscript_numbers.py`, closes the loop on the numbers those three do
not own (lineage-level AUROCs, TSG/ONC contrast, MSI stratification,
mutation-type analysis, direction audit, therapeutic-window module, Table S2
spot values): it recomputes 109 manuscript claims from the artifacts under
`output/` and writes `output/manuscript_number_audit.tsv` (109/109 match).
One command runs all of them plus the test suite:

```bash
./verify_all.sh              # fast: reuses cached data slices (~30 s)
VERIFY_FULL=1 ./verify_all.sh  # full: rebuilds all caches from raw CSVs
```

| Script | Recomputes | Claims check |
|---|---|---|
| `compute_headline_metrics.py` | DD values, TI, AUROC, pair counts, q-values | `output/headline_metrics.json` (16/16 match) |
| `ml_benchmark.py` | LOO-CV AUROC of LR/RF/SVM vs DD-only baseline | `output/ml_benchmark.json` |
| `regression_controls.py` | CNV/expression/TP53-adjusted regressions, CNV R² | `output/regression_controls.json` (10/10 match) |
| `audit_manuscript_numbers.py` | All remaining manuscript numbers (lineages, MSI, mutation type, direction, therapeutic window) | `output/manuscript_number_audit.tsv` (109/109 match) |

`regression_controls.py` also writes `output/cnv_independence.csv` and
`output/cnv_scatter_sample.csv`, the exact inputs of `R_figS4.R` (Fig. S4).
The test suite (`pytest tests/`, 31 tests) covers the R-package mirror and
shared pipeline utilities.

## Data Availability

- **DepMap 26Q1**: https://depmap.org/portal/download/
- **CPTAC**: via cBioPortal API (cached in `data/cptac_cache/`)
- **PRISM**: https://depmap.org/portal/download/
- **Processed tables**: `output/tables/TableS1-S10`

## Citation

```
Mo Q, Zhu T. Delta Dependency Prioritizes Candidate Paralog Dependencies
Across Solid Tumor Types. Manuscript under review at Genome Biology (2026).
Code archive DOI: 10.5281/zenodo.21502030 (concept DOI; resolves to the latest release)
```

## License

MIT License. See [LICENSE](LICENSE).
