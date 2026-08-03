#!/usr/bin/env python3
"""Section-by-section review patch (round: 主文稿分章节修改清单).

Applies text-level revisions to manuscript.tex / supplementary.tex /
cover_letter.md / submission_checklist.md / submission_texts.md /
submission_2026-07-31/README_投稿包.md, and regenerates Supplementary
Tables S5/S9/S10 LaTeX blocks from the recomputed artifacts (signed DWS).
Every replacement asserts an exact occurrence count. No numeric value is
hard-coded into table regeneration: all numbers come from output/ TSV/JSON.
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
MS = ROOT / "manuscript.tex"
SUP = ROOT / "supplementary.tex"
CL = ROOT / "cover_letter.md"
CKL = ROOT / "submission_checklist.md"
TXT = ROOT / "submission_texts.md"
RDM = ROOT.parent / "submission_2026-07-31" / "README_投稿包.md"

FAILURES = []


def sub(path, old, new, count=1, tag=""):
    text = path.read_text(encoding="utf-8")
    n = text.count(old)
    if n != count:
        FAILURES.append(f"{tag or path.name}: expected {count}x but found {n}x :: {old[:90]!r}")
        return
    text = text.replace(old, new)
    path.write_text(text, encoding="utf-8")
    print(f"OK  {tag or path.name}: {old[:60]!r} -> {new[:60]!r}")


# =========================================================================
# 1. manuscript.tex — abstract
# =========================================================================
sub(MS,
    "\\textbf{Results:} We introduce Delta Dependency (DD), an interpretable metric measuring the dependency shift on a paralog between driver-mutant and wild-type cell lines in DepMap (1,208 lines, 23 solid tumor types). On twelve curated, evidence-tiered pairs in the pre-specified primary benchmark (three gynecological lineages), signed DD achieved AUROC~$=$~0.629 (pair-clustered bootstrap 95\\% CI 0.253--0.933; 0.613 excluding two DepMap-labelled pairs; 0.525 for the eight pairs with pre-DepMap experimental evidence; the sole lineage-evaluable Tier A/B pair, ARID1A$\\rightarrow$ARID1B, ranked above all unlabeled controls in both evaluable lineages). At pair level on the same 72 pairs, the interpretable composite score reached AUROC~$=$~0.831, indistinguishable from the best classifier (SVM-RBF, 0.841; six positives, so descriptive). Exploratory analysis of seven CPTAC cohorts ($n=771$) detected protein-level paralog co-variation that is robust to tumor-purity adjustment. Dependency-window scoring recovered the established ARID1A$\\rightarrow$ARID1B pair (DWS~$=$~2.82) as the leading selective candidate and nominated KMT2D$\\rightarrow$KMT2C for experimental follow-up.",
    "\\textbf{Results:} We introduce Delta Dependency (DD), an interpretable metric measuring the dependency shift on a paralog between driver-mutant and wild-type cell lines in DepMap (1,208 lines, 23 solid tumor types). On twelve curated, evidence-tiered pairs in the primary literature-derived benchmark (three gynecological lineages), signed DD achieved AUROC~$=$~0.629 in a positive--unlabeled, internal-consistency evaluation (pair-clustered bootstrap 95\\% CI 0.253--0.933; 0.613 excluding two DepMap-labelled pairs; 0.525 for the eight pairs with pre-DepMap experimental evidence; the sole lineage-evaluable Tier A/B pair, ARID1A$\\rightarrow$ARID1B, ranked above all unlabeled controls in both evaluable lineages). Exploratory analysis of seven CPTAC cohorts ($n=771$) detected protein-level paralog co-variation. Dependency-window scoring recovered the established ARID1A$\\rightarrow$ARID1B pair (DWS~$=$~2.82) as the highest-selectivity established candidate and nominated KMT2D$\\rightarrow$KMT2C as the highest-ranked non-benchmark candidate for experimental follow-up.",
    tag="ms:abstract")

# =========================================================================
# 2. manuscript.tex — Background
# =========================================================================
sub(MS,
    "combinatorial knockout screens have begun to map paralog synthetic lethality directly \\cite{Parrish2021,Thompson2021,Dede2020,EsmaeiliAnvar2024}",
    "combinatorial knockout screens have begun to map paralog synthetic lethality directly \\cite{Parrish2021,Thompson2021,Dede2020,EsmaeiliAnvar2024,Harle2025,Flister2025}",
    tag="ms:bg-cite")

sub(MS,
    "Three practical gaps remain. Complex multi-feature models are less directly traceable to a single experimentally measurable effect: when a neural network predicts SL, it is difficult to reconstruct \\textit{why}, which hinders rational experimental design. The optimal readout for paralog compensation is unsettled, since RNA-level and protein-level signals can diverge; a systematic CPTAC-based survey of paralog protein compensation appeared only recently \\cite{Venkatesh2025}, and how such protein-level compensation relates to dependency-based prioritization remains unexplored. Finally, few studies have layered orthogonal evidence (proteomics, drug sensitivity, clinical biomarkers, sequence- and structure-derived descriptors) onto genomic prediction to move from a ranked list toward an experimentally prioritized candidate set.",
    "Two practical gaps remain. First, complex multi-feature models are less directly traceable to a single experimentally measurable effect: when a neural network predicts SL, it is difficult to reconstruct \\textit{why}, which hinders rational experimental design. Second, moving from a ranked list to an experimentally prioritized candidate set requires conservative integration of orthogonal evidence layers; the optimal readout for paralog compensation is itself unsettled, since RNA-level and protein-level signals can diverge \\cite{Venkatesh2025}, and how protein-level compensation relates to dependency-based prioritization remains unexplored.",
    tag="ms:bg-gaps")

sub(MS,
    "MSI-aware patient stratification, mutation-type analysis, dependency-window quantification, and exploratory structure-derived targetability descriptors.",
    "MSI-aware patient stratification, mutation-type analysis, dependency-window quantification, and exploratory structure-derived targetability descriptors. The primary objective was to assess whether signed DD enriches literature-supported mutation-conditioned paralog dependencies relative to unlabeled paralog pairs; secondary analyses evaluated robustness and selectivity, and the proteomic, pharmacologic, clinical, and sequence-derived analyses were exploratory.",
    tag="ms:bg-objective")

# =========================================================================
# 3. manuscript.tex — Results part 1 (benchmark, classifiers)
# =========================================================================
sub(MS,
    "reported by Feng et al.\\ (2024) \\cite{Feng2024} (Table~1).",
    "reported by Feng et al.\\ (2024) \\cite{Feng2024} (Supplementary Table~S12).",
    tag="ms:tab1-ref")

sub(MS,
    "with no detectable difference between them at this sample size (the paired-bootstrap confidence interval is too wide to establish equivalence), followed by random forest (0.722)",
    "with no detectable difference between them at this sample size (paired-bootstrap $\\Delta$AUROC~$=-0.01$, 95\\% CI $-0.21$ to $+0.18$; the interval is too wide to establish equivalence), followed by random forest (0.722)",
    tag="ms:delta-auroc")

sub(MS,
    "; in the full LR fit no individual feature reached significance (all five feature $p$-values $>0.23$)",
    "",
    tag="ms:lr-wald-del")

sub(MS,
    "the Tier A~$\\cup$~Tier B set constitutes the pre-specified primary external benchmark",
    "the Tier A~$\\cup$~Tier B set constitutes the primary literature-derived benchmark",
    tag="ms:benchmark-name")

# Fig 1 caption
sub(MS,
    "the Tier A~$\\cup$~Tier B set (five pairs) constitutes the primary external benchmark. The evaluation frame spans",
    "the Tier A~$\\cup$~Tier B set (five pairs) constitutes the primary literature-derived benchmark. The evaluation frame spans",
    tag="ms:fig1c-benchmark")
sub(MS,
    "Performance is reported under three pre-specified aggregation frameworks (see ``Evaluation frameworks'').",
    "Performance is reported under three aggregation frameworks defined before final reporting (see ``Evaluation frameworks'').",
    tag="ms:fig1c-agg")

# delete Table 1 block (moved to supplementary as Table S12)
TABLE1_BLOCK = """\\begin{table}[H]
\\centering
\\caption{\\textbf{DD performance in context of published CV3 results (not a head-to-head benchmark).} Published values are CV3 (gene-pair isolation) AUROCs from Feng et al.\\ (2024), Supplementary Data 1 (NSMRand negative sampling, 1:1 positive:negative ratio, complete dataset). The DD value from this study was evaluated on the full lineage-level frame (110 driver--paralog--lineage entries, 8 positives). \\textit{Evaluation frameworks differ: published methods were tested on general SL gene-pair universes with different datasets, tasks, and positive:negative ratios; DD was tested on paralog-SL pairs only. Direct AUROC comparison is not warranted.}}
\\label{tab:benchmark}
\\begin{tabular}{lc}
\\toprule
\\textbf{Method} & \\textbf{CV3 AUROC} \\\\
\\midrule
DD (this study) & 0.629 \\\\
SLMGAE \\cite{Hao2021} & 0.790 \\\\
NSF4SL \\cite{Wang2022} & 0.683 \\\\
GCATSL \\cite{Long2021} & 0.678 \\\\
GRSMF \\cite{Huang2019} & 0.656 \\\\
PiLSL \\cite{Liu2022} & 0.626 \\\\
KG4SL \\cite{Wang2021} & 0.563 \\\\
SLGNN \\cite{Zhu2023} & 0.530 \\\\
PTGNN \\cite{Long2022} & 0.529 \\\\
\\bottomrule
\\end{tabular}
\\end{table}
"""
sub(MS, TABLE1_BLOCK, "", tag="ms:table1-del")

# =========================================================================
# 4. manuscript.tex — CPTAC section
# =========================================================================
sub(MS,
    "\\subsection*{Orthogonal proteomic evidence: protein-level paralog co-variation}",
    "\\subsection*{Proteomic context: paralog protein co-variation}",
    tag="ms:cptac-title")

sub(MS,
    "These results serve as orthogonal mechanistic support, not causal validation.",
    "These results provide orthogonal contextual evidence of co-variation, not direct validation of conditional dependency.",
    tag="ms:cptac-framing")

sub(MS,
    "a directionally consistent but preliminary signal for compensatory protein upregulation upon driver loss.",
    "a directionally concordant but statistically inconclusive signal.",
    tag="ms:fig2d-caption")

# =========================================================================
# 5. manuscript.tex — survival sentence
# =========================================================================
sub(MS,
    "higher \\textit{ARID1B} expression was associated with shortened overall survival (multivariable Cox HR~$=$~1.31 per SD, 95\\% CI: 1.087--1.569, $p=0.004$, adjusted for age and AJCC stage; Supplementary Fig.~S10).",
    "higher \\textit{ARID1B} expression was associated with shortened overall survival (multivariable Cox HR~$=$~1.31 per SD, 95\\% CI: 1.087--1.569, $p=0.004$, adjusted for age and AJCC stage; Supplementary Fig.~S10), though the association reached FDR significance only after covariate adjustment and does not establish a relationship between ARID1B expression and ARID1A-conditioned dependency.",
    tag="ms:survival")

# =========================================================================
# 6. manuscript.tex — PRISM proliferation-rate qualifier
# =========================================================================
sub(MS,
    "the relaxed threshold is deliberate because this is a screen, not a confirmatory analysis.",
    "the relaxed threshold is deliberate because this is a screen, not a confirmatory analysis; proliferation rate was not controlled, so these associations are quality-control context only.",
    tag="ms:prism-qc")

# =========================================================================
# 7. manuscript.tex — DWS equation + surrounding text (signed formula)
# =========================================================================
sub(MS,
    "\\text{DWS} = \\frac{|\\text{DD}|}{\\max(|\\mu_{\\text{Chronos}}|,\\;f_{\\text{pan-essential}},\\;0.01)}",
    "\\text{DWS} = \\frac{\\max(\\text{DD},\\,0)}{\\max(|\\mu_{\\text{Chronos}}|,\\;f_{\\text{pan-essential}},\\;0.01)}",
    tag="ms:dws-eq")

sub(MS,
    "and selectivity is the mutant-minus-wild-type difference in that fraction between driver-mutant and wild-type groups. Each paralog was classified",
    "and selectivity is the mutant-minus-wild-type difference in that fraction between driver-mutant and wild-type groups. The signed numerator credits only compensation-direction shifts; the $|$DD$|$-numerator variant, which also credits reverse-direction shifts, is reported as a sensitivity analysis (Supplementary Table~S10). Each paralog was classified",
    tag="ms:dws-signed-note")

sub(MS,
    "within the dependency-window module ARID1A$\\rightarrow$ARID1B attained the highest mean selectivity and the second-highest mean $|$DD$|$ (0.270) of all evaluated pairs (Supplementary Table~S5).",
    "within the dependency-window module ARID1A$\\rightarrow$ARID1B attained the highest mean selectivity and the highest mean positive DD (0.270) of all evaluated pairs (Supplementary Table~S5).",
    tag="ms:arid-numerator")

sub(MS,
    "The highest raw DWS value was attained by NF1$\\rightarrow$RASA2 (5.42), but this estimate rests on a near-zero pan-essential denominator (0.1\\% of cell lines) with selectivity~$\\approx$~0, so NF1$\\rightarrow$RASA2 is not a selective candidate.",
    "NF1$\\rightarrow$RASA2 attains a moderate signed DWS (1.82) that derives entirely from a single context (Endometrial), where the 0.01 denominator floor is active; its bootstrap rank interval is wide (95\\% CI 1--18) and selectivity is~$\\approx$~0, so NF1$\\rightarrow$RASA2 is not a selective candidate.",
    tag="ms:nf1-dws")

sub(MS,
    "Nine paralog pairs maintained DWS~$>$~1.0 in $\\ge$2 cancer contexts (Supplementary Fig.~S8).",
    "Five paralog pairs maintained DWS~$>$~1.0 in $\\ge$2 cancer contexts (Supplementary Fig.~S8).",
    tag="ms:dws-n-gt1")

# =========================================================================
# 8. manuscript.tex — composite prioritization (Table S9 values)
# =========================================================================
sub(MS,
    "A pre-specified composite score (max-normalized DWS 0.40, rescaled selectivity 0.30, targetability 0.30) ranked ARID1A$\\rightarrow$ARID1B second (0.631) behind NF1$\\rightarrow$RASA2 (0.695; inflated by the near-zero pan-essential denominator noted above), leaving ARID1A$\\rightarrow$ARID1B as the leading selective candidate (Supplementary Table~S9).",
    "A composite score with heuristically fixed weights (max-normalized DWS 0.40, rescaled selectivity 0.30, targetability 0.30) ranked ARID1A$\\rightarrow$ARID1B first (0.823), followed by EP300$\\rightarrow$CREBBP (0.644) and NF1$\\rightarrow$RASA2 (0.554), making ARID1A$\\rightarrow$ARID1B the highest-priority established mutation-conditioned pair under the combined criteria (Supplementary Table~S9).",
    tag="ms:composite")

sub(MS,
    "NF1$\\rightarrow$RASA2 has high DWS (5.42) but near-zero selectivity (0.005), hence MODERATE.",
    "NF1$\\rightarrow$RASA2 has moderate signed DWS (1.82) resting on a floored denominator but near-zero selectivity (0.005), hence MODERATE.",
    tag="ms:fig4-caption")

# =========================================================================
# 9. manuscript.tex — Discussion
# =========================================================================
sub(MS,
    "We set out to test whether a single, interpretable metric — the dependency shift between driver-mutant and wild-type cells — could compete with multi-feature classifiers for paralog-SL prediction. On the 72-pair benchmark frame (three gynecological lineages), the interpretable composite score reached AUROC~$=$~0.831 in head-to-head comparison, with no detectable difference from the best multi-feature classifier (SVM-RBF, 0.841; classifier range 0.240--0.841) and a confidence interval too wide to establish equivalence,",
    "We set out to test whether a single, interpretable metric --- the dependency shift between driver-mutant and wild-type cells --- could capture the paralog-SL prioritization signal. DD provided a transparent univariate signal, whereas multi-feature performance estimates were unstable because only six positive pairs were available. On the 72-pair benchmark frame (three gynecological lineages), the interpretable composite score reached AUROC~$=$~0.831 in head-to-head comparison, with no detectable difference from the best multi-feature classifier (SVM-RBF, 0.841; classifier range 0.240--0.841; paired-bootstrap $\\Delta$AUROC~$=-0.01$, 95\\% CI $-0.21$ to $+0.18$, too wide to establish equivalence),",
    tag="ms:disc-opening")

sub(MS,
    "The mutation-conditioned analysis in UCEC remains directionally consistent with compensatory protein stabilization: ARID1B protein abundance was higher in ARID1A-mutant tumors (fold change $\\approx$1.05, $p=0.082$), though the result is not significant and needs replication in additional cohorts and pairs.",
    "The mutation-conditioned analysis in UCEC was directionally concordant but statistically inconclusive: ARID1B protein abundance was higher in ARID1A-mutant tumors (fold change $\\approx$1.05, $p=0.082$, BH $q=0.27$), and the result needs replication in additional cohorts and pairs.",
    tag="ms:disc-cptac")

sub(MS,
    "\\textbf{TSG versus oncogene contexts.}",
    "\\textbf{Directionality: TSG versus oncogene contexts.}",
    tag="ms:disc-tsg-title")

sub(MS,
    "Our DWS framework nominated ARID1A$\\rightarrow$ARID1B (DWS~$=$~2.82) as the leading selective candidate while flagging BRCA1/2 paralogs as pan-essential.",
    "Our DWS framework nominated ARID1A$\\rightarrow$ARID1B (DWS~$=$~2.82) as the highest-selectivity established candidate --- SMARCA4$\\rightarrow$SMARCA2 had the higher DWS (4.87) --- while flagging BRCA1/2 paralogs as pan-essential.",
    tag="ms:disc-dws")

sub(MS,
    "Both lineage-evaluable positive entries on the Tier A~$\\cup$~Tier B external benchmark (ARID1A$\\rightarrow$ARID1B in two lineages) ranked above all unlabeled controls;",
    "Both lineage-evaluable positive entries on the Tier A~$\\cup$~Tier B literature-derived benchmark (ARID1A$\\rightarrow$ARID1B in two lineages) ranked above all unlabeled controls;",
    tag="ms:disc-independence")

sub(MS,
    "Under the pre-specified multiple-testing framework (within-driver Benjamini--Hochberg correction across each driver's HGNC paralogs)",
    "Under the within-driver multiple-testing framework (Benjamini--Hochberg correction across each driver's HGNC paralogs)",
    tag="ms:disc-fdr")

sub(MS,
    "The strongest of these is KMT2D$\\rightarrow$KMT2C (DWS~$=$~1.77, MODERATE selectivity $+0.056$; composite score 0.477, fourth overall)",
    "The highest-ranked non-benchmark candidate is KMT2D$\\rightarrow$KMT2C (signed DWS~$=$~1.26, MODERATE selectivity $+0.056$; composite score 0.525, fourth overall)",
    tag="ms:kmt2d")

sub(MS,
    "RB1$\\rightarrow$RBL1 attains a higher raw DWS (3.44) through a near-zero essentiality denominator with selectivity~$\\approx$~0.007 and is therefore not a selective candidate.",
    "RB1$\\rightarrow$RBL1 shows a reverse-direction shift, so its signed DWS is zero and it is not a selective candidate (its former $|$DD$|$-based DWS of 3.44 rested on a near-zero essentiality denominator with selectivity~$\\approx$~0.007).",
    tag="ms:rb1")

sub(MS,
    "is the aggregation used for the per-pair classifier features (Table~1).",
    "is the aggregation used for the per-pair classifier features (Methods).",
    tag="ms:tab1-ref2")

sub(MS,
    "\\item \\textit{Benchmark scope:} The pre-specified evaluation universe comprises",
    "\\item \\textit{Benchmark scope:} The evaluation universe comprises",
    tag="ms:lim-scope")

# =========================================================================
# 10. manuscript.tex — Conclusions
# =========================================================================
sub(MS,
    "We introduce Delta Dependency, a single-subtraction metric for paralog-based SL prediction, and evaluate it across 23 solid tumor types. Signed DD achieves an internal-consistency AUROC of 0.629 for the paralog-SL task within DepMap, with a wide pair-clustered confidence interval (0.253--0.933) that reflects the small positive set (the sole lineage-evaluable Tier A/B pair, ARID1A$\\rightarrow$ARID1B, ranked above all unlabeled controls in both evaluable lineages), and the interpretable composite score reaches 0.831 as a descriptive comparison; published deep learning results in an analogous CV3 setting remain higher (SLMGAE \\cite{Hao2021}, 0.790) but are contextual reference points, not head-to-head benchmarks. We note that these evaluations use different test sets; a rigorous head-to-head benchmark would require running all methods on identical paralog-SL data. CPTAC proteomics across seven cohorts reveals protein-level paralog co-variation that is robust to tumor-purity adjustment and is paralleled by weaker, pair-dependent mRNA-level co-variation, and mutation-conditioned analysis in UCEC provides initial evidence for compensatory protein upregulation. Mutation type emerges as a candidate stratification variable for prospective study, whereas MSI status showed no significant modulation of the DD signal. Dependency-window scoring combined with PRISM drug screening and exploratory structure-derived descriptors recovers the established ARID1A$\\rightarrow$ARID1B pair as the leading selective candidate --- a framework-validating result --- and nominates KMT2D$\\rightarrow$KMT2C as the strongest unlabeled candidate for experimental follow-up. DD is an interpretable, association-based statistic for ranking candidate paralog dependencies: its current evidence supports hypothesis generation and experimental prioritization, not causal confirmation of synthetic lethality or prediction of clinical therapeutic windows.",
    "Signed DD provides a direction-aware and experimentally interpretable statistic for prioritizing mutation-conditioned paralog dependencies, evaluated here across 23 solid tumor types. Its current performance estimates remain uncertain because the literature-derived positive set is small, heterogeneous, and partly dependent on DepMap-derived evidence. The framework recovered established mutation-conditioned dependencies --- including ARID1A$\\rightarrow$ARID1B, the highest-selectivity established candidate --- and prioritized KMT2D$\\rightarrow$KMT2C as the highest-ranked non-benchmark candidate for further testing, but neither DD nor the dependency-window score establishes causal synthetic lethality or an in vivo therapeutic window. Independent combinatorial perturbation in isogenic models is the essential next step.",
    tag="ms:conclusions")

# =========================================================================
# 11. manuscript.tex — Methods
# =========================================================================
sub(MS,
    "The Tier A~$\\cup$~Tier B set (five pairs) constitutes the primary external benchmark; the two tiers differ in evidence independence",
    "The Tier A~$\\cup$~Tier B set (five pairs) constitutes the primary literature-derived benchmark; the two tiers differ in evidence independence",
    tag="ms:methods-benchmark")

sub(MS,
    "DWS defined as in Equation~\\ref{eq:ti}. The denominator uses",
    "DWS defined as in Equation~\\ref{eq:ti} with the signed numerator $\\max(\\mathrm{DD},0)$, which credits only compensation-direction shifts; the $|$DD$|$-numerator variant is reported as a sensitivity analysis (Supplementary Table~S10). The denominator uses",
    tag="ms:methods-dws")

sub(MS,
    "all weights were fixed a priori in the analysis code (\\texttt{alphafold\\_analysis.py}) and none of the descriptors was calibrated against external drug-development data.",
    "all weights were fixed heuristically before final reporting and are unchanged in the archived release (commit \\texttt{66fc633}; \\texttt{alphafold\\_analysis.py}), and none of the descriptors was calibrated against external drug-development data.",
    tag="ms:methods-weights")

sub(MS,
    "The pre-specified primary endpoint is the AUROC and AUPRC of signed DD on the Tier A~$\\cup$~Tier B external benchmark",
    "The primary endpoint is the AUROC and AUPRC of signed DD on the Tier A~$\\cup$~Tier B literature-derived benchmark",
    tag="ms:methods-endpoint")

# =========================================================================
# 12. manuscript.tex — supplementary information + additional files
# =========================================================================
sub(MS,
    "Supplementary Figures S1--S10 and Supplementary Tables S1--S11 are available as a separate PDF file (Additional file 1); machine-readable TSV mirrors of all tables are provided as Additional files 2--12 (see Additional files).",
    "Supplementary Figures S1--S10 and Supplementary Tables S1--S12 are available as a separate PDF file (Additional file 1); machine-readable TSV mirrors of the data tables are provided as Additional files 2--12 (see Additional files).",
    tag="ms:suppinfo-head")

sub(MS,
    "Supplementary Table~S5: All 21 paralog pairs ranked by dependency window score with component values and bootstrap confidence intervals.",
    "Supplementary Table~S5: All 21 paralog pairs ranked by signed dependency window score with numerator, denominator, bootstrap confidence intervals, and bootstrap rank intervals.",
    tag="ms:suppinfo-s5")

sub(MS,
    "Supplementary Table~S11: Full regression model table for confounder controls (base, CNV-, expression-, and lineage-adjusted models with robust standard errors and BH $q$-values).",
    "Supplementary Table~S11: Full regression model table for confounder controls (base, CNV-, expression-, and lineage-adjusted models with robust standard errors and BH $q$-values). Supplementary Table~S12: Contextual literature values under non-comparable evaluation settings (published CV3 AUROCs from Feng et al.\\ with the DD value from this study; evaluation frameworks differ and direct comparison is not warranted).",
    tag="ms:suppinfo-s12")

sub(MS,
    "392 claims in total",
    "405 claims in total",
    tag="ms:claims-405")

sub(MS,
    "\\noindent\\textbf{Additional file 1:} Supplementary Information (PDF). Supplementary Methods, Supplementary Figures S1--S10, and Supplementary Tables S1--S11.",
    "\\noindent\\textbf{Additional file 1:} Supplementary Information (PDF). Supplementary Methods, Supplementary Figures S1--S10, and Supplementary Tables S1--S12.",
    tag="ms:addfile1")

sub(MS,
    "\\noindent\\textbf{Additional file 6:} \\texttt{TableS5\\_DWS.tsv}. All 21 paralog pairs ranked by dependency window score with component values and bootstrap 95\\% confidence intervals.",
    "\\noindent\\textbf{Additional file 6:} \\texttt{TableS5\\_DWS.tsv}. All 21 paralog pairs ranked by signed dependency window score with numerator, denominator, floor flag, and bootstrap 95\\% confidence intervals.",
    tag="ms:addfile6")

# =========================================================================
# 13. supplementary.tex — Methods / captions terminology
# =========================================================================
sub(SUP,
    "The Tier~A~$\\cup$~Tier~B set constitutes the pre-specified primary external benchmark, with an evidence-independence asymmetry:",
    "The Tier~A~$\\cup$~Tier~B set constitutes the primary literature-derived benchmark, with an evidence-independence asymmetry:",
    tag="sup:benchmark-name")

sub(SUP,
    "Standardized LR coefficients with Wald $p$-values were obtained from a full-data maximum-likelihood fit (statsmodels 0.14). ",
    "",
    tag="sup:lr-wald-del")

sub(SUP,
    "The DWS (Equation 2 of the main text) is the ratio of $|$DD$|$ to a conservative baseline-essentiality denominator, $\\max(|\\mu_{\\text{Chronos}}|,\\;f_{\\text{pan-essential}},\\;0.01)$, where",
    "The DWS (Equation 2 of the main text) is the ratio of the signed numerator $\\max(\\mathrm{DD},0)$ --- which credits only compensation-direction shifts --- to a conservative baseline-essentiality denominator, $\\max(|\\mu_{\\text{Chronos}}|,\\;f_{\\text{pan-essential}},\\;0.01)$, where",
    tag="sup:dws-def")

sub(SUP,
    "and by bootstrap 95\\% confidence intervals for mean DWS and selectivity from 1{,}000 stratified resamples of cell lines within each context (seed 42; Tables~S5 and~S10).",
    "and by bootstrap 95\\% confidence intervals for mean DWS and selectivity, plus per-pair bootstrap rank intervals, from 1{,}000 stratified resamples of cell lines within each context (seed 42; Tables~S5 and~S10). The $|$DD$|$-numerator variant (pre-revision formula) is retained as a sensitivity analysis (Table~S10).",
    tag="sup:dws-robust")

sub(SUP,
    "with all weights fixed a priori in \\texttt{alphafold\\_analysis.py} and no calibration against drug-development data.",
    "with all weights fixed heuristically before final reporting in \\texttt{alphafold\\_analysis.py} and no calibration against drug-development data.",
    tag="sup:weights")

sub(SUP,
    "weights fixed a priori in \\texttt{pcs.py}",
    "weights fixed heuristically in \\texttt{pcs.py}",
    tag="sup:pcs-weights")

# Fig S8 caption
sub(SUP,
    "NF1$\\rightarrow$RASA2 attains a higher raw DWS through a near-zero pan-essential denominator with selectivity~$\\approx$~0.",
    "NF1$\\rightarrow$RASA2's signed DWS (1.82) derives entirely from the Endometrial context, where the 0.01 denominator floor is active, with selectivity~$\\approx$~0 and a wide bootstrap rank interval (95\\% CI 1--18).",
    tag="sup:figs8-caption")

# Fig S9 caption
sub(SUP,
    "d, Composite prioritization score (max-normalized DWS 0.40, rescaled selectivity 0.30, targetability 0.30; weights fixed a priori in \\texttt{alphafold\\_analysis.py}). NF1$\\rightarrow$RASA2 ranks first (0.695) only through a near-zero pan-essential denominator that inflates its raw DWS (selectivity~$\\approx$~0); ARID1A$\\rightarrow$ARID1B (0.631) is the leading selective candidate.",
    "d, Composite prioritization score (max-normalized signed DWS 0.40, rescaled selectivity 0.30, targetability 0.30; weights fixed heuristically before final reporting in \\texttt{alphafold\\_analysis.py}). ARID1A$\\rightarrow$ARID1B ranks first (0.823), followed by EP300$\\rightarrow$CREBBP (0.644); NF1$\\rightarrow$RASA2 ranks third (0.554), its signed DWS resting on a floored denominator with selectivity~$\\approx$~0.",
    tag="sup:figs9-caption")

# Table S3 caption
sub(SUP,
    "The Tier A~$\\cup$~Tier B set (5 pairs) constitutes the primary external benchmark.",
    "The Tier A~$\\cup$~Tier B set (5 pairs) constitutes the primary literature-derived benchmark.",
    tag="sup:s3-caption")
sub(SUP,
    "\\textbf{Inclusion}: Primary = primary external benchmark; Secondary = reported separately; Comparator = specificity reference.",
    "\\textbf{Inclusion}: Primary = primary literature-derived benchmark; Secondary = reported separately; Comparator = specificity reference.",
    tag="sup:s3-inclusion")

# TCGA survival extension
sub(SUP,
    "The BRCA association thus replicates directionally in OV but not in UCEC.",
    "The BRCA association was directionally concordant in OV (nominal age-adjusted $p=0.013$, $q=0.406$, not FDR-significant) but not in UCEC; in BRCA it reached FDR only after covariate adjustment (univariable $q=0.129$; multivariable $q=0.040$), so its significance is adjustment-dependent.",
    tag="sup:tcga-ext")

# =========================================================================
# 14. cover_letter.md
# =========================================================================
sub(CL,
    "Both lineage-evaluable positives on the Tier A + Tier B external benchmark rank above all unlabeled controls.",
    "Both lineage-evaluable positives on the Tier A + Tier B literature-derived benchmark rank above all unlabeled controls.",
    tag="cl:benchmark")

sub(CL,
    "CPTAC proteomics across seven cohorts (771 tumours) shows consistent protein-level paralog co-variation (RNA-level signal weak: AUROC = 0.547).",
    "CPTAC proteomics across seven cohorts (771 tumours) shows consistent protein-level paralog co-variation (RNA abundance shifts showed weak discrimination of curated dependencies, AUROC = 0.547, although RNA co-variation itself was detectable).",
    tag="cl:cptac")

sub(CL,
    "nominate ARID1A→ARID1B (DWS = 2.82) as the leading selective candidate",
    "nominate ARID1A→ARID1B (DWS = 2.82) as the highest-selectivity established candidate",
    tag="cl:arid")

sub(CL,
    "recomputing 392 numeric claims",
    "recomputing 405 numeric claims",
    tag="cl:405")

# =========================================================================
# 15. submission_checklist.md / README_投稿包.md / submission_texts.md
# =========================================================================
sub(CKL, "392/392 claims match", "405/405 claims match", tag="ckl:405a")
sub(CKL, "all 392 audited numeric claims", "all 405 audited numeric claims", tag="ckl:405b")

sub(RDM, "**392/392 claims match; 31 passed; ALL CHECKS PASSED**",
    "**405/405 claims match; 31 passed; ALL CHECKS PASSED**", tag="rdm:405a")
sub(RDM,
    "本轮 ARS 多视角审稿修复已写入：SVM 数字（0.841/0.240）、聚类 bootstrap 与置换检验、CPTAC RNA/纯度结论改述、TCGA UCEC/OV 生存扩展。",
    "本轮 ARS 多视角审稿修复已写入：SVM 数字（0.841/0.240）、聚类 bootstrap 与置换检验、CPTAC RNA/纯度结论改述、TCGA UCEC/OV 生存扩展。\n2. **2026-08-02 分章节评审修复**：DWS 主公式改 signed max(DD,0)（|DD| 版降为敏感性分析；NF1→RASA2 5.42→1.82，DWS>1 in ≥2 contexts 9→5）；composite 优先级分数同步切换 signed DWS（ARID1A→ARID1B 升为第一 0.823）；主文 Table 1 移入补充材料为 Table S12；摘要/结论按评审重写（删 composite/classifier 细节，标注 PU/internal-consistency）；统一术语 primary literature-derived benchmark；删除 pre-specified/a priori 表述；审计口径 392→405 claims。",
    tag="rdm:round-note")

# submission_texts.md: graphical abstract + highlights + running title
sub(TXT, "12 gold-standard pairs", "12 curated, evidence-tiered pairs", tag="txt:pairs")
sub(TXT, "AUROC 0.629 on the tiered gold-standard set", "AUROC 0.629 on the tiered curated set", tag="txt:set")
sub(TXT, "1.000 on the Tier A+B external benchmark", "1.000 on the Tier A+B literature-derived benchmark", tag="txt:benchmark")
sub(TXT,
    "4. Dependency-window analysis nominates ARID1A-ARID1B as leading selective candidate",
    "4. Dependency-window analysis ranks established ARID1A-ARID1B top for selectivity",
    tag="txt:highlight4")
sub(TXT,
    "# 投稿系统单独提交项：Graphical Abstract 说明文字 + Highlights",
    "# 投稿系统单独提交项：Graphical Abstract 说明文字 + Highlights\n\n## Running title\n\nDelta Dependency for paralog prioritization",
    tag="txt:running-title")

# =========================================================================
# 16. Regenerate supplementary Table S5 (signed DWS, from artifacts)
# =========================================================================
s5 = pd.read_csv(ROOT / "output/tables/TableS5_DWS.tsv", sep="\t")
rank_ci = {(r["driver"], r["paralog"]): r["rank_ci95"]
           for r in json.loads((ROOT / "output/dws_robustness.json").read_text())["bootstrap_rank_ci95"]}

CLASS_MAP = {"HIGH_SELECTIVITY": "HS", "MODERATE": "MOD",
             "LOW_SELECTIVITY": "LS", "PAN_ESSENTIAL": "PE"}
TYPE_MAP = {"sequence paralog": "Seq", "partial homolog": "Part",
            "functional analog": "Func"}


def num(x, nd):
    r = round(float(x), nd)
    if r == 0:
        return f"{0:.{nd}f}"
    if r < 0:
        return f"$-${abs(r):.{nd}f}"
    return f"{r:.{nd}f}"


def ci(lo, hi, nd=2):
    return f"{num(lo, nd)}--{num(hi, nd)}"


def rankfmt(v):
    return str(int(v)) if float(v) == int(v) else f"{v:.1f}"


rows5 = []
for _, r in s5.iterrows():
    key = (r["driver"], r["paralog"])
    rc = rank_ci[key]
    denom = num(r["pan_essentiality_denominator"], 3)
    if str(r["floor_0_01_active"]) == "True":
        denom += "$^{\\dagger}$"
    rows5.append(
        f"{r['driver']} & {r['paralog']} & {num(r['mean_max_dd_0'],3)} & {denom} & "
        f"{num(r['dws'],3)} & {num(r['selectivity'],3)} & {CLASS_MAP[r['classification']]} & "
        f"{TYPE_MAP[r['pair_type']]} & {ci(r['dws_ci95_lo'], r['dws_ci95_hi'])} & "
        f"{ci(r['selectivity_ci95_lo'], r['selectivity_ci95_hi'])} & "
        f"{rankfmt(rc[0])}--{rankfmt(rc[1])} \\\\")
S5_ROWS = "\n".join(rows5)

S5_OLD_START = "\\begin{tabular}{llrrrllcc}"
S5_OLD_END = "\\end{tabular}"
sup_text = SUP.read_text(encoding="utf-8")
i0 = sup_text.index(S5_OLD_START)
i1 = sup_text.index(S5_OLD_END, i0)
new_tabular5 = (
    "\\begin{tabular}{llrrrlrlll}\n"
    "\\toprule\n"
    "\\textbf{Driver} & \\textbf{Paralog} & \\textbf{Num.} & \\textbf{Denom.} & \\textbf{DWS} & "
    "\\textbf{Select.} & \\textbf{Class} & \\textbf{DWS 95\\% CI} & \\textbf{Select.\\ 95\\% CI} & "
    "\\textbf{Rank 95\\% CI} \\\\\n"
    "\\midrule\n" + S5_ROWS + "\n\\bottomrule\n"
)
sup_text = sup_text[:i0] + new_tabular5 + sup_text[i1:]

S5_CAP_OLD = ("\\caption{All 21 analyzed paralog pairs ranked by mean dependency window score (DWS) on the "
    "$\\geq$5-mutant frame. $|$DD$|$: mean absolute Delta Dependency across cancer contexts. Select.: mean selectivity "
    "(fraction mutant-essential minus fraction wild-type-essential). Class: HS=HIGH\\_SELECTIVITY, MOD=MODERATE, "
    "LS=LOW\\_SELECTIVITY, PE=PAN\\_ESSENTIAL. Type: Seq=sequence paralog, Part=partial homology, Func=functional analog. "
    "DWS 95\\% CI and Select.\\ 95\\% CI: bootstrap 95\\% percentile intervals from 1{,}000 stratified resamples of "
    "cell lines within each cancer context (seed 42; \\texttt{dws\\_robustness.py}); the DWS and Select.\\ columns "
    "report the observed point estimates, not bootstrap means. Three pairs evaluable under the previous frame "
    "(AKT1$\\rightarrow$AKT2, CCNE1$\\rightarrow$CCNE2, CDK4$\\rightarrow$CDK6) no longer meet the $\\geq$5-mutant "
    "threshold in any context and are not listed. The active denominator varies by pair (the 0.01 floor governs pairs "
    "with near-zero baseline essentiality, e.g., NF1$\\rightarrow$RASA2), and the classification thresholds are "
    "operational stratifiers, not validated clinical cutoffs (threshold sensitivity in Table~S10).}")
S5_CAP_NEW = ("\\caption{All 21 analyzed paralog pairs ranked by the signed dependency window score "
    "(DWS $=$ mean $\\max(\\mathrm{DD},0)$ $/$ denominator; Equation 2 of the main text) on the $\\geq$5-mutant frame. "
    "Num.: DWS numerator, mean $\\max(\\mathrm{DD},0)$ across cancer contexts (only compensation-direction shifts "
    "contribute). Denom.: baseline-essentiality denominator $\\max(|\\mu_{\\text{Chronos}}|, "
    "f_{\\text{pan-essential}}, 0.01)$; $^{\\dagger}$marks the only pair whose signed DWS derives from a context in "
    "which the 0.01 floor is active (NF1$\\rightarrow$RASA2, Endometrial). Select.: mean selectivity (fraction "
    "mutant-essential minus fraction wild-type-essential). Class: HS=HIGH\\_SELECTIVITY, MOD=MODERATE, "
    "LS=LOW\\_SELECTIVITY, PE=PAN\\_ESSENTIAL. DWS 95\\% CI and Select.\\ 95\\% CI: bootstrap 95\\% percentile "
    "intervals from 1{,}000 stratified resamples of cell lines within each cancer context (seed 42; "
    "\\texttt{dws\\_robustness.py}); the DWS and Select.\\ columns report the observed point estimates, not bootstrap "
    "means. Rank 95\\% CI: bootstrap 95\\% interval of each pair's DWS rank across the same resamples. Three pairs "
    "evaluable under the previous frame (AKT1$\\rightarrow$AKT2, CCNE1$\\rightarrow$CCNE2, CDK4$\\rightarrow$CDK6) "
    "do not meet the $\\geq$5-mutant threshold in any context and are not listed. The $|$DD$|$-numerator variant "
    "(pre-revision formula) is reported as a sensitivity analysis (Table~S10); classification thresholds are "
    "operational stratifiers, not validated clinical cutoffs.}")
assert sup_text.count(S5_CAP_OLD) == 1, "S5 caption anchor not found"
sup_text = sup_text.replace(S5_CAP_OLD, S5_CAP_NEW)
SUP.write_text(sup_text, encoding="utf-8")
print("OK  sup:TableS5 regenerated (21 rows, signed DWS, numerator/denominator/rank CI)")

# =========================================================================
# 17. Regenerate supplementary Table S9 (composite, signed DWS)
# =========================================================================
s9 = pd.read_csv(ROOT / "output/tables/TableS9_CompositeScore.tsv", sep="\t")
class_lookup = {(r["driver"], r["paralog"]): CLASS_MAP[r["classification"]]
                for _, r in s5.iterrows()}
rows9 = []
for _, r in s9.head(10).iterrows():
    known = "$\\star$" if str(r["is_known_sl"]) == "True" else ""
    rows9.append(
        f"{r['driver']:<10} & {r['paralog']:<8} & {r['dws']:.2f} & {class_lookup[(r['driver'], r['paralog'])]} & "
        f"{known} & {r['structural_similarity']:.3f} & {r['clinical_targetability']:.3f} \\\\")
S9_ROWS = "\n".join(rows9)

S9_OLD_BLOCK = """NF1      & RASA2    & 5.42 & MOD &         & 0.658 & 0.695 \\\\
ARID1A   & ARID1B   & 2.82 & HS  & $\\star$ & 0.952 & 0.631 \\\\
EP300    & CREBBP   & 1.82 & MOD & $\\star$ & 0.976 & 0.535 \\\\
KMT2D    & KMT2C    & 1.77 & MOD &         & 0.829 & 0.477 \\\\
PIK3CA   & PIK3CB   & 1.36 & LS  & $\\star$ & 0.847 & 0.471 \\\\
PPP2R1A  & PPP2R1B  & 0.93 & LS  & $\\star$ & 0.921 & 0.452 \\\\
TP53     & TP63     & 0.81 & MOD &         & 0.803 & 0.393 \\\\
KRAS     & HRAS     & 0.21 & LS  &         & 0.885 & 0.392 \\\\
FBXW7    & FBXW2    & 0.72 & LS  & $\\star$ & 0.666 & 0.358 \\\\
STK11    & SIK1     & 0.39 & LS  & $\\star$ & 0.710 & 0.334 \\\\"""
sup_text = SUP.read_text(encoding="utf-8")
assert sup_text.count(S9_OLD_BLOCK) == 1, "S9 row block anchor not found"
sup_text = sup_text.replace(S9_OLD_BLOCK, S9_ROWS)

S9_CAP_NOTE_OLD = ("\\textit{NF1$\\rightarrow$RASA2's top composite score is driven by a near-zero pan-essential "
    "denominator (0.1\\% of cell lines) that inflates its DWS; with selectivity~$\\approx$~0 it is not a selective "
    "candidate, and ARID1A$\\rightarrow$ARID1B remains the leading selective candidate.")
S9_CAP_NOTE_NEW = ("\\textit{ARID1A$\\rightarrow$ARID1B ranks first under the signed-DWS composite; "
    "NF1$\\rightarrow$RASA2's earlier top score was driven by a near-zero pan-essential denominator (0.1\\% of cell "
    "lines) that inflated its $|$DD$|$-based DWS, and with selectivity~$\\approx$~0 it is not a selective candidate.")
assert sup_text.count(S9_CAP_NOTE_OLD) == 1, "S9 caption note anchor not found"
sup_text = sup_text.replace(S9_CAP_NOTE_OLD, S9_CAP_NOTE_NEW)
sub_count = sup_text.count("weights fixed a priori in \\texttt{alphafold\\_analysis.py}")
if sub_count:
    sup_text = sup_text.replace("Target.: composite prioritization score (max-normalized DWS 0.40, rescaled selectivity 0.30, targetability 0.30; weights fixed a priori in \\texttt{alphafold\\_analysis.py})",
                                "Target.: composite prioritization score (max-normalized signed DWS 0.40, rescaled selectivity 0.30, targetability 0.30; weights fixed heuristically before final reporting in \\texttt{alphafold\\_analysis.py})")
SUP.write_text(sup_text, encoding="utf-8")
print("OK  sup:TableS9 regenerated (top-10 rows, signed-DWS composite)")

# =========================================================================
# 18. Regenerate supplementary Table S10 (sensitivity, incl. |DD| row)
# =========================================================================
s10 = pd.read_csv(ROOT / "output/tables/TableS10_DWS_Sensitivity.tsv", sep="\t")
S10_OLD_ROWS = """Base: $\\max(|\\mu|, f, 0.01)$ (production) & 1.000 & 5 \\\\
$|\\mu|$ only & 0.990 & 4 \\\\
$f$ (pan-essential fraction) only & 0.893 & 4 \\\\
mean$(|\\mu|, f)$ & 0.988 & 4 \\\\
Floor 0.001 & 1.000 & 5 \\\\
Floor 0.05 & 0.974 & 5 \\\\"""
v = s10.set_index("variant")
S10_NEW_ROWS = (
    "Base: signed $\\max(\\mathrm{DD},0)$ numerator, $\\max(|\\mu|, f, 0.01)$ denominator (production) & 1.000 & 5 \\\\\n"
    f"$|\\mu|$ only & {v.loc['denominator |mu| only','spearman_rho_vs_base']:.3f} & {int(v.loc['denominator |mu| only','top5_overlap_with_base'])} \\\\\n"
    f"$f$ (pan-essential fraction) only & {v.loc['denominator f only','spearman_rho_vs_base']:.3f} & {int(v.loc['denominator f only','top5_overlap_with_base'])} \\\\\n"
    f"mean$(|\\mu|, f)$ & {v.loc['denominator mean(|mu|, f)','spearman_rho_vs_base']:.3f} & {int(v.loc['denominator mean(|mu|, f)','top5_overlap_with_base'])} \\\\\n"
    f"Floor 0.001 & {v.loc['floor 0.001','spearman_rho_vs_base']:.3f} & {int(v.loc['floor 0.001','top5_overlap_with_base'])} \\\\\n"
    f"Floor 0.05 & {v.loc['floor 0.05','spearman_rho_vs_base']:.3f} & {int(v.loc['floor 0.05','top5_overlap_with_base'])} \\\\\n"
    f"$|$DD$|$ numerator (pre-revision formula, sensitivity) & {v.loc['numerator |DD| (pre-revision formula, sensitivity)','spearman_rho_vs_base']:.3f} & {int(v.loc['numerator |DD| (pre-revision formula, sensitivity)','top5_overlap_with_base'])} \\\\"
)
sup_text = SUP.read_text(encoding="utf-8")
assert sup_text.count(S10_OLD_ROWS) == 1, "S10 row block anchor not found"
sup_text = sup_text.replace(S10_OLD_ROWS, S10_NEW_ROWS)
S10_CAP_OLD = "\\caption{Sensitivity of DWS-based rankings to the denominator definition (upper block)"
S10_CAP_NEW = "\\caption{Sensitivity of DWS-based rankings to the numerator and denominator definitions (upper block)"
assert sup_text.count(S10_CAP_OLD) == 1
sup_text = sup_text.replace(S10_CAP_OLD, S10_CAP_NEW)

# =========================================================================
# 19. Insert Table S12 (moved main-text Table 1) before Supplementary References
# =========================================================================
S12_BLOCK = """
\\subsection{Table S12: Contextual literature values under non-comparable evaluation settings}
\\begin{table}[H]
\\centering
\\caption{\\textbf{Contextual literature values under non-comparable evaluation settings.} Published values are CV3 (gene-pair isolation) AUROCs from Feng et al.\\ (2024; main-text reference list), Supplementary Data 1 (NSMRand negative sampling, 1:1 positive:negative ratio, complete dataset); individual method citations appear in the main-text reference list. The DD value from this study was evaluated on the full lineage-level frame (110 driver--paralog--lineage entries, 8 positives). \\textit{Evaluation frameworks differ: published methods were tested on general SL gene-pair universes with different datasets, tasks, and positive:negative ratios; DD was tested on paralog-SL pairs only. Direct AUROC comparison is not warranted.}}
\\label{tab:benchmark}
\\begin{tabular}{lc}
\\toprule
\\textbf{Method} & \\textbf{CV3 AUROC} \\\\
\\midrule
DD (this study) & 0.629 \\\\
SLMGAE & 0.790 \\\\
NSF4SL & 0.683 \\\\
GCATSL & 0.678 \\\\
GRSMF & 0.656 \\\\
PiLSL & 0.626 \\\\
KG4SL & 0.563 \\\\
SLGNN & 0.530 \\\\
PTGNN & 0.529 \\\\
\\bottomrule
\\end{tabular}
\\end{table}

"""
anchor = "\\section{Supplementary References}"
assert sup_text.count(anchor) == 1, "Supplementary References anchor not found"
sup_text = sup_text.replace(anchor, S12_BLOCK + anchor)
SUP.write_text(sup_text, encoding="utf-8")
print("OK  sup:TableS10 regenerated + TableS12 inserted")

# =========================================================================
# 20. Verification greps
# =========================================================================
import re
checks = {
    MS: ["pre-specified", "external benchmark", "leading selective candidate",
         "a priori", "5.42", "0.695", "0.631", "Table~1", "392 claims"],
    SUP: ["pre-specified", "external benchmark", "a priori", "5.421",
          "0.695", "0.631", "leading selective candidate", "Standardized LR coefficients"],
    CL: ["external benchmark", "leading selective candidate", "392"],
}
for path, pats in checks.items():
    t = path.read_text(encoding="utf-8")
    for p in pats:
        if p in t:
            for m in re.finditer(re.escape(p), t):
                ctx = t[max(0, m.start()-60):m.end()+60].replace("\n", " ")
                print(f"LEFTOVER {path.name}: {p!r} :: ...{ctx}...")

if FAILURES:
    print("\n===== FAILURES =====")
    for f in FAILURES:
        print("FAIL:", f)
    sys.exit(1)
print("\nALL PATCHES APPLIED")
