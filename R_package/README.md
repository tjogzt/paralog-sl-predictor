# paralogSL

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21502114.svg)](https://doi.org/10.5281/zenodo.21502114)

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

# Load DepMap data
dep <- load_dependency("CRISPRGeneEffect.csv")
models <- load_models("Model.csv")
mut_mat <- build_mutation_matrix(mutations_file = "OmicsSomaticMutations.csv",
                                  cell_lines = rownames(dep),
                                  genes = c("ARID1A", "TP53", "PIK3CA"))

# Compute DD for a single driver-paralog pair
result <- compute_dd(dep, driver_gene = "ARID1A", paralog_gene = "ARID1B",
                     mut_lines = rownames(mut_mat)[mut_mat[,"ARID1A"] == 1],
                     wt_lines  = rownames(mut_mat)[mut_mat[,"ARID1A"] == 0])
print(result)
# $DD = 0.182, $p_value = 1.4e-26, $cohens_d = 0.59
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
| `compute_auroc()` | AUROC evaluation against gold standard |
| `compute_therapeutic_window()` | TI-based selectivity classification |
| `run_paralog_analysis()` | Full pipeline wrapper |
| `plot_pancancer_summary()` | Cross-cancer AUROC visualization |
| `plot_therapeutic_window()` | TI bubble plot |

## Built-in Datasets

- `known_sl_pairs`: 12 gold-standard paralog-SL pairs
- `gyn_drivers`: Driver gene sets for gynecological cancers
- `solid_tumor_summary`: DD AUROC across 23 solid tumor types
- `benchmark_methods`: Published method performance (CV3)

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
