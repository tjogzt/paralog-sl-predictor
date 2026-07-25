# 投稿系统单独提交项：Graphical Abstract 说明文字 + Highlights

## Graphical Abstract — Caption

Delta Dependency (DD) framework supported by exploratory analyses across five orthogonal layers. DD measures the shift in mean Chronos dependency of a paralog gene between driver wild-type and mutant cell lines (manuscript Eq. 1). Top: from DepMap 26Q1 data (1,208 cell lines; 66,595 HGNC paralog pairs; 12 gold-standard pairs), DD nominates priority paralog-SL candidates without any training and compares favorably with published methods (AUROC 0.794 vs best published CV3 result 0.790; 1.000 with a ≥30% sequence-identity filter). Bottom: exploratory support across genomic (23 solid tumor lineages; AUROC > 0.7 in 9 of 17 evaluable lineages), proteomic (seven CPTAC cohorts, n = 771), pharmacologic (1,482 PRISM compounds; 553 selective associations), clinical-stratification (MSI status and mutation type), and structural (sequence features, domain architecture, PROTAC suitability) layers. Fully reproducible from public data with no GPU requirement.

## Highlights（每条 ≤85 字符）

1. A single-subtraction metric (Delta Dependency) predicts paralog synthetic lethality
2. DD alone (AUROC 0.736) outperforms four multi-feature classifiers in cross-validation
3. Protein-level paralog compensation is confirmed across seven CPTAC cohorts (n = 771)
4. Therapeutic window analysis nominates ARID1A-ARID1B as the leading SL candidate
5. Every quantitative claim is recomputable from raw data with a single command

<!-- 字符数核对（提交前已逐条计数）：
1: 83  2: 85  3: 84  4: 79  5: 76 -->
