# Cover Letter

**Date:** June 19, 2026

**To:** The Editors
*Genome Biology*

**Re:** Submission of "Paralog Compensation Predicts Synthetic Lethality Across 27 Solid Tumor Types"

---

Dear Editors,

We submit for consideration in *Genome Biology* our manuscript entitled **"Paralog Compensation Predicts Synthetic Lethality Across 27 Solid Tumor Types."** This work presents Delta Dependency (DD)—a simple, univariate metric for paralog-based synthetic lethality (SL) prediction—and validates it through five orthogonal layers of evidence spanning genomic, proteomic, pharmacologic, clinical-stratification, and structural analyses.

**Key contributions:**

1. **DD outperforms deep learning.** Our single-subtraction metric achieves AUROC = 0.794 in gynecological cancers, exceeding eight published deep learning methods (best published CV3: DDSL, 0.720). With a ≥30% sequence identity filter, AUROC reaches 1.000. Component decomposition, bootstrap resampling, and negative controls confirm signal robustness.

2. **Pan-cancer scope: 27 solid tumor types, 7 CPTAC cohorts.** This is the largest systematic evaluation of paralog-SL across human cancers. DD is effective (AUROC > 0.7) in 8 of 21 evaluable lineages. CPTAC proteomics across seven independent cohorts (BRCA, COAD, LUAD, GBM, PDAC, UCEC, LUSC; n = 672 samples) confirms paralog compensation operates at the protein level, with EP300↔CREBBP significant in 5/7 cohorts, while RNA-level signals remain undetectable.

3. **PRISM drug sensitivity validates the framework.** A systematic scan of 1,482 compounds identifies 553 drugs with significant selective killing in driver-mutant cell lines. Known biology is recapitulated (MEK inhibitors → KRAS-mutant, mTOR inhibitors → PTEN-mutant, HDAC inhibitors → EP300-mutant), and novel associations are uncovered.

4. **Clinical stratification and safety analysis.** MSI status stratifies patients for paralog-SL signal (MSS tumors show stronger signal), mutation type modulates compensation strength (truncating > missense for well-characterized TSGs), and a quantitative therapeutic window framework identifies ARID1A→ARID1B as the safest paralog-SL target.

5. **Fully reproducible, no specialized infrastructure.** All analyses use publicly available DepMap and CPTAC data. The method requires only subtraction—no GPU, no model training. We provide an open-source R package (`paralogSL`) and a complete Python pipeline.

We believe this manuscript is well-suited for *Genome Biology* because it combines methodological simplicity with exceptionally broad validation, addresses a clinically important problem (SL-based cancer target discovery), and provides resources (code, data, R package) that the community can immediately use. The pan-cancer CPTAC validation and PRISM drug selectivity analysis substantially extend the evidence base beyond prior computational SL studies.

All data are from public repositories (DepMap, cBioPortal, PDC, TCGA). All code is available at [GitHub repository].

**Competing interests:** The authors declare no competing interests.

**Author contributions:** T.Z. conceived the study, developed the methodology, performed all analyses, wrote the manuscript, and developed the R package.

**Acknowledgments:** We thank the DepMap, CPTAC, TCGA, and cBioPortal teams for making their data publicly available.

We thank you for your consideration and look forward to your response.

Sincerely,

Tao Zhu
(on behalf of all authors)

---

**Suggested Reviewers:**
1. Dr. Francisca Vazquez — DepMap / Broad Institute (CRISPR dependency expertise)
2. Dr. Michael P. Snyder — Stanford (proteogenomics, CPTAC expertise)
3. Dr. Jason Moffat — University of Toronto (synthetic lethality screening)
4. Dr. Rameen Beroukhim — Dana-Farber / Broad Institute (cancer genomics, paralog dependencies)
5. Dr. Bing Zhang — Baylor College of Medicine (CPTAC proteogenomics)

**Opposed Reviewers:** None.
