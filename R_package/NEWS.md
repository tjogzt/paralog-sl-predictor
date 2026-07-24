# paralogSL News

## v1.0.0 (2026-05-25)

### Initial Release

- `compute_dd()` — Delta Dependency computation for driver-paralog pairs
- `compute_pcs()` — Paralog Compensation Score
- `compute_auroc()` — AUROC evaluation against gold-standard paralog-SL pairs
- `compute_therapeutic_window()` — Therapeutic Index (TI) selectivity classification
- `run_paralog_analysis()` — Full pipeline wrapper across multiple driver genes
- `build_mutation_matrix()` — Binary mutation matrix from DepMap annotations
- `classify_msi_status()` — MSI/dMMR classification from MMR gene mutations
- `stratify_by_mutation_type()` — Truncating vs missense stratification
- `predict_trial_response()` — Clinical trial priority scoring
- `load_dependency()`, `load_expression()`, `load_models()` — DepMap data loaders
- `filter_cancer_cell_lines()` — Cancer-type filtering by Oncotree annotations
- `get_benchmark_table()`, `plot_benchmark()` — Method comparison
- `plot_pancancer_summary()`, `plot_therapeutic_window()` — Visualization
- `plot_protein_correlation()`, `visualize_dependency_shift()` — Pair-level plots
- Built-in datasets: `gyn_drivers`, `known_sl_pairs`, `benchmark_methods`,
  `cross_cancer_summary`, `solid_tumor_summary`, `therapeutic_window_summary`
- Full vignette with reproducible workflow
