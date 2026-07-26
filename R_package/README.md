# paralogSL

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21502113.svg)](https://doi.org/10.5281/zenodo.21502113)

**Paralog-Based Synthetic Lethality Prioritization via Delta Dependency**

An R package for prioritizing paralog-based synthetic lethality (SL) candidates using Delta Dependency (DD) analysis on DepMap Chronos gene-effect data.

## Installation

```r
# From GitHub (recommended)
devtools::install_github("tjogzt/paralogSL")

# Or from local source
install.packages("path/to/paralogSL", repos = NULL, type = "source")
```

## Quick Start (3 lines)

```r
library(paralogSL)

# Load DepMap data (25Q2+ release format: ModelID / HugoSymbol columns)
dep <- load_dependency("CRISPRGeneEffect.csv")
mut <- data.table::fread("OmicsSomaticMutations.csv",
                         select = c("ModelID", "HugoSymbol", "IsDefaultEntryForModel"))
mut <- mut[IsDefaultEntryForModel == "Yes"]
mut_df <- data.frame(DepMap_ID = mut$ModelID, Gene = mut$HugoSymbol)
mut_mat <- build_mutation_matrix(mut_df,
                                  cell_lines = rownames(dep),
                                  genes = c("ARID1A", "TP53", "PIK3CA"))

# Compute DD for a single driver-paralog pair
result <- compute_dd(dep, driver_gene = "ARID1A", paralog_gene = "ARID1B",
                     mut_lines = rownames(mut_mat)[mut_mat[,"ARID1A"] == 1],
                     wt_lines  = rownames(mut_mat)[mut_mat[,"ARID1A"] == 0])
print(result)
# $DD = 0.187, $p_value = 1.8e-10, $cohens_d = 0.69  (n_mut = 164, n_wt = 1044;
# pan-cancer "any annotated mutation" example, DepMap 26Q1 data)
```

## Workflow

```
DepMap gene-effect matrix
    → build_mutation_matrix()
    → compute_dd() per driver × paralog × lineage
    → compute_auroc() against gold standard
    → compute_therapeutic_window() for translational ranking
    → plot_pancancer_summary() / plot_therapeutic_window()
```

## Key Functions

| Function | Description |
|----------|-------------|
| `compute_dd()` | Delta Dependency: DD(D,P,c) = G_WT - G_MUT |
| `compute_pcs()` | Paralog Compensation Score |
| `compute_auroc()` | AUROC against evidence-tiered gold standard (Tier A+B by default) |
| `compute_therapeutic_window()` | DWS/TI selectivity classification |
| `run_paralog_analysis()` | Full pipeline wrapper |
| `plot_pancancer_summary()` | Cross-cancer AUROC visualization |
| `plot_therapeutic_window()` | DWS/TI bubble plot |

## Built-in Datasets

- `known_sl_pairs`: 12 evidence-tiered gold-standard pairs (Tier A/B/C + comparators; `compute_auroc()` uses Tier A+B by default)
- `gyn_drivers`: Driver gene sets for five cancer types
- `solid_tumor_summary`: Per-lineage DD AUROC, primary min>=5 frame (19 lineages)
- `cross_cancer_summary`: Per-lineage DD AUROC, sensitivity min>=3 frame (23 lineages)
- `therapeutic_window_summary`: 21-pair DWS classification (manuscript Table S6)
- `benchmark_methods`: Published CV3 contextual references + this-study values (not a head-to-head benchmark)

## Citation

Mo Q, Zhu T. Delta Dependency Prioritizes Paralog-Based Synthetic Lethality
Candidates Across Solid Tumor Types. *Genome Biology* (2026).

```bibtex
@article{Mo2026,
  title   = {Delta Dependency Prioritizes Paralog-Based Synthetic Lethality
             Candidates Across Solid Tumor Types},
  author  = {Mo, Qingqing and Zhu, Tao},
  journal = {Genome Biology},
  year    = {2026},
  doi     = {10.5281/zenodo.21502031},
}
```

## License

MIT
