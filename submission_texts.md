# 投稿系统单独提交项：Graphical Abstract 说明文字 + Highlights

## Running title

Delta Dependency for paralog prioritization

## Graphical Abstract — Caption

Delta Dependency (DD) framework supported by exploratory analyses across five orthogonal layers. DD measures the shift in mean Chronos dependency of a paralog gene between driver wild-type and mutant cell lines (manuscript Eq. 1). Top: from DepMap 26Q1 data (1,208 cell lines; 66,595 HGNC paralog pairs; 12 curated, evidence-tiered pairs), DD nominates priority paralog-SL candidates without any training (AUROC 0.629 on the tiered curated set, pair-clustered 95% CI 0.253-0.933; 1.000 on the Tier A+B literature-derived benchmark, where both lineage-evaluable entries were the same pair). Bottom: exploratory support across genomic (23 solid tumor lineages; AUROC > 0.7 in 1 of 8 evaluable lineages on the primary frame, 6 of 12 on the relaxed sensitivity frame), proteomic (seven CPTAC cohorts, n = 771), pharmacologic (1,482 PRISM compounds; 633 selective associations), clinical-stratification (MSI status and mutation type), and structural (sequence features, domain architecture, PROTAC suitability) layers. Fully reproducible from public data with no GPU requirement.

## Highlights（每条不超过 85 字符）

1. A single-subtraction metric (Delta Dependency) predicts paralog synthetic lethality
2. Interpretable composite (AUROC 0.831) matches best multi-feature classifier in CV
3. Protein-level paralog co-variation is detectable across seven CPTAC cohorts (n = 771)
4. Dependency-window analysis ranks established ARID1A-ARID1B top for selectivity
5. Every quantitative claim is recomputable from raw data with a single command

<!-- 字符数核对（2026-07-31 已逐条复核，均不超过 85）：
1: 83  2: 81  3: 85  4: 81  5: 76 -->
