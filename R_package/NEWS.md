# paralogSL News

## v1.1.2 (2026-07-28)

### Signed-DD primary metric release (paralog-sl-predictor v1.3.2 framework)

- `compute_auroc()` now ranks by signed DD by default (the manuscript's
  primary metric: positive DD = compensation); new `score` argument with
  `score = "abs"` restoring the pre-v1.1.2 direction-agnostic behaviour,
  which the manuscript reports only as a sensitivity analysis
- All built-in datasets regenerated from the canonical artifacts:
  `solid_tumor_summary` / `cross_cancer_summary` now carry signed
  `dd_auroc` plus a `dd_auroc_abs_sensitivity` column (primary frame:
  1 of 8 evaluable lineages above AUROC 0.7; sensitivity frame: 6 of 12)
- `benchmark_methods` this-study DD row updated to the signed value
  (full lineage-level frame AUROC 0.629; |DD| sensitivity 0.676;
  Tier A + B external benchmark unchanged at 1.000)

## v1.1.0 (2026-07-26)

### Manuscript-alignment release (paralog-sl-predictor v1.1.0 framework)

- All built-in datasets regenerated from the canonical artifacts of the
  companion analysis repository; `data-raw/make_data.R` now points at
  `../paralog_sl_predictor/output/` and fails loudly instead of silently
  falling back to stale hard-coded values (root cause of earlier drift)
- `known_sl_pairs` is now evidence-tiered (Tier A/B/C + comparators with
  `direct_sl`, `inclusion`, `key_ref`), matching supplementary Table S3
- `compute_auroc()` gains a `tiers` argument and defaults to the manuscript's
  primary external benchmark (Tier A + B); `tiers = NULL` restores the
  pre-v1.1.0 all-pairs behavior
- `get_benchmark_table()` now returns the shipped `benchmark_methods`
  dataset instead of hard-coded values; this-study rows updated to the
  current frames (DD 0.676 full lineage-level frame; DD + ID >= 0.3
  1.000 on a 3-pair high-identity subset, anecdotal). Published CV3
  values are documented as contextual references, not a head-to-head
  benchmark (supersedes the 0.736/0.794 figures quoted in earlier NEWS)
- `solid_tumor_summary` / `cross_cancer_summary` now carry the primary
  (min>=5; 19 lineages) and sensitivity (min>=3; 23 lineages) pan-cancer
  frames with `dd_auroc` columns
- `therapeutic_window_summary` replaced by the 21-pair DWS classification
  (manuscript Table S5); TI terminology documented as DWS (dependency
  window score), formula unchanged and identical to the Python pipeline
- `classify_msi_status()` documented as the MMR/POLE mutation proxy
- README quick-start example re-verified against DepMap 26Q1 raw data
  (DD = 0.187, p = 1.8e-10, Cohen's d = 0.69, n_mut = 164, n_wt = 1044)
- Vignette reproducibility section rewritten to the current headline
  evaluation; citation author list corrected (Mo & Zhu)

## v1.0.1 (2026-07-25)

### Documentation & test corrections

- Vignette and testthat expectations aligned with the audited leave-one-out
  cross-validation benchmark in the companion analysis repository
  (paralog-sl-predictor v1.1.0); DD-only baseline (AUROC 0.736) outperforms
  all benchmarked classifiers
- No changes to package functions or APIs
- Added `CITATION.cff` and `inst/CITATION`

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
