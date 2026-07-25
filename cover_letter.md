# Cover Letter

**Date:** July 25, 2026

**To:** The Editors
*Genome Biology*

**Re:** Submission of "Delta Dependency Prioritizes Paralog-Based Synthetic Lethality Candidates Across Solid Tumor Types"

---

Dear Editors,

We submit for consideration in *Genome Biology* our manuscript entitled **"Delta Dependency Prioritizes Paralog-Based Synthetic Lethality Candidates Across Solid Tumor Types."** This work presents Delta Dependency (DD), a simple, univariate metric for paralog-based synthetic lethality (SL) prediction, and validates it through five orthogonal layers of evidence spanning genomic, proteomic, pharmacologic, clinical-stratification, and structural analyses.

**Key contributions:**

1. **DD outperforms multi-feature classifiers in a rigorous head-to-head benchmark.** On the primary test set of 10 true sequence-paralog pairs, DD achieves AUROC = 0.837; on the full 12-pair set, AUROC = 0.794. In a rigorous head-to-head benchmark on the same 77 pairs (8 known positives), DD alone (AUROC = 0.736) outperforms four multi-feature classifiers — logistic regression (0.632), random forest (0.629), SVM-Linear (0.699), and SVM-RBF (0.563) — under leave-one-pair-out cross-validation. A DD + sequence-identity (≥30%) filter reaches AUROC = 1.000 on a high-identity subset (p = 0.0048 by exact permutation).

2. **Pan-cancer scope: 23 solid tumor types, 7 CPTAC cohorts.** This is the largest systematic evaluation of paralog-SL across human cancers (66,595 HGNC paralog pairs × 1,208 DepMap cell lines). DD exceeds AUROC = 0.7 in 9 of 17 evaluable lineages. CPTAC proteomics across seven independent cohorts (BRCA, COAD, LUAD, GBM, PDAC, UCEC, LUSC; n = 771 samples) confirms paralog co-variation at the protein level — EP300↔CREBBP significant in 5 of 7 cohorts — while RNA-level signal is indistinguishable from random (AUROC = 0.348).

3. **PRISM drug sensitivity validates the framework.** A discovery-stage scan of 1,482 compounds identifies 553 compound–paralog associations with selective killing in driver-mutant cell lines. Known biology is recapitulated (MEK inhibitors → KRAS-mutant, mTOR/AKT inhibitors → PTEN-mutant, HDAC inhibitors → EP300-mutant), and novel associations are uncovered.

4. **Clinical stratification and safety analysis.** MSI status stratifies patients for paralog-SL signal (MSS tumors show stronger signal), mutation type modulates compensation strength (truncating > missense for well-characterized TSGs), and a quantitative therapeutic window framework nominates ARID1A→ARID1B (DWS = 4.13) as the leading selective candidate, with a biomarker-enriched trial design for MSS tumors harboring truncating ARID1A mutations.

5. **Every number is one-command reproducible.** All analyses use publicly available DepMap and CPTAC data, and the method requires only subtraction: no GPU, no model training. Beyond the open-source R package (`paralogSL` v1.0.2) and Python pipeline, we provide three audit scripts that recompute every quantitative claim in the manuscript directly from raw data — headline metrics (16/16 claims verified), classifier benchmarks, and covariate-adjusted regressions (10/10 claims verified) — each ending with an automated claim-by-claim comparison that fails loudly on any mismatch. A single entry point (`verify_all.sh`) runs all audits plus the test suite (31 tests) in under one minute. The submission-version code is archived on Zenodo (https://doi.org/10.5281/zenodo.21502030; R package: https://doi.org/10.5281/zenodo.21502113).

We believe this manuscript is well-suited for *Genome Biology* because it combines methodological simplicity with exceptionally broad validation, addresses a clinically important problem (SL-based cancer target discovery), and sets a high standard for computational reproducibility that the community can immediately build on. The pan-cancer CPTAC validation and PRISM drug selectivity analysis substantially extend the evidence base beyond prior computational SL studies.

All data are from public repositories (DepMap, cBioPortal, PDC, TCGA). All code is openly available: analysis pipeline at https://github.com/tjogzt/paralog-sl-predictor and R package at https://github.com/tjogzt/paralogSL.

**Competing interests:** The authors declare no competing interests.

**Author contributions:** Q.Q.M. and T.Z. conceived the study. T.Z. developed the methodology, performed all computational analyses, and developed the R package. Q.Q.M. curated clinical datasets and provided gynecological oncology domain expertise. Q.Q.M. and T.Z. wrote the manuscript.

**Acknowledgments:** We thank the DepMap, CPTAC, TCGA, and cBioPortal teams for making their data publicly available.

We confirm that this manuscript is original, has not been published previously, and is not under consideration by any other journal. Both authors have read and approved the final manuscript and agree with its submission to *Genome Biology*.

We thank you for your consideration and look forward to your response.

Sincerely,

Qingqing Mo and Tao Zhu
$^{1}$Department of Obstetrics and Gynecology, National Clinical Research Center for Obstetrics and Gynecology, Tongji Hospital, Tongji Medical College, Huazhong University of Science and Technology, Wuhan, China
$^{2}$Key Laboratory of Cancer Invasion and Metastasis (Ministry of Education), Hubei Key Laboratory of Tumor Invasion and Metastasis, Tongji Hospital, Tongji Medical College, Huazhong University of Science and Technology, Wuhan, China

$^{\ast}$Correspondence: zhutao@tjh.tjmu.edu.cn

---

**Suggested Reviewers:**
1. Dr. Francisca Vazquez — DepMap / Broad Institute (CRISPR dependency expertise)
2. Dr. Michael P. Snyder — Stanford (proteogenomics, CPTAC expertise)
3. Dr. Jason Moffat — University of Toronto (synthetic lethality screening)
4. Dr. Rameen Beroukhim — Dana-Farber / Broad Institute (cancer genomics, paralog dependencies)
5. Dr. Bing Zhang — Baylor College of Medicine (CPTAC proteogenomics)

**Opposed Reviewers:** None.
