#!/usr/bin/env python3
"""patch_gb_format.py — Genome Biology format-compliance round (2026-08-01).

Implements the format/figure checklist from the external review:
  * double line spacing (journal requirement)
  * full affiliation addresses with postal codes
  * ethics statement reworded per human-data guidance
  * Fig 1c: cross-study bar chart replaced by benchmark-selection flowchart
  * Table 1: drop Interpretability column and the anecdotal DD+ID 1.000 row
  * Fig 3c (TCGA survival forest) moved to Supplementary Fig. S10; main-text
    survival paragraph compressed to an exploratory summary
  * Fig 4a (PRISM anchors) moved to Supplementary Fig. S7b; Fig 4 becomes the
    single-panel dependency-window classification
  * figure captions: n-unit / error-bar statements; titles <=15 words
  * paralogSL version string v1.1.1 -> v1.1.2

Every replacement carries a uniqueness assertion so the script is safe to
re-run for inspection (it aborts before writing if any anchor is missing or
non-unique).
"""
from pathlib import Path

ROOT = Path(__file__).parent
MS = ROOT / "manuscript.tex"
SI = ROOT / "supplementary.tex"

ms = MS.read_text(encoding="utf-8")
si = SI.read_text(encoding="utf-8")

EDITS_MS = []
EDITS_SI = []


def add(bucket, old, new, tag):
    bucket.append((old, new, tag))


# ───────────────────────── manuscript.tex ─────────────────────────

# 1. double line spacing (journal requirement; was onehalfspacing)
add(EDITS_MS,
    r"\usepackage{setspace}\onehalfspacing",
    r"\usepackage{setspace}\doublespacing",
    "doublespacing")

# 2/3. full affiliation addresses with postal codes
add(EDITS_MS,
    "Huazhong University of Science and Technology, Wuhan, China\\\\[0.4em]",
    "Huazhong University of Science and Technology, 1095 Jiefang Avenue, Wuhan 430030, China\\\\[0.4em]",
    "affiliation-1 address")
add(EDITS_MS,
    "Tumor Invasion and Metastasis, Tongji Hospital, Tongji Medical College, Huazhong University of Science and Technology, Wuhan, China",
    "Tumor Invasion and Metastasis, Tongji Hospital, Tongji Medical College, Huazhong University of Science and Technology, 1095 Jiefang Avenue, Wuhan 430030, China",
    "affiliation-2 address")

# 4. ethics statement: explain WHY approval was not required
add(EDITS_MS,
    "Not applicable. This study used exclusively publicly available, de-identified data.",
    "No institutional ethics approval was required: this study analysed exclusively publicly available, de-identified data (DepMap, CPTAC, TCGA, PRISM).",
    "ethics statement")

# 5. paralogSL version string
add(EDITS_MS,
    r"The \texttt{paralogSL} R package (v1.1.1, 20 exported functions)",
    r"The \texttt{paralogSL} R package (v1.1.2, 20 exported functions)",
    "R package version")

# 6/7. CV3 comparison: cite Table 1 only; add exploratory-sensitivity sentence
add(EDITS_MS,
    r"reported by Feng et al.\ (2024) \cite{Feng2024} (Fig.~1c; Table~1).",
    r"reported by Feng et al.\ (2024) \cite{Feng2024} (Table~1).",
    "CV3 ref line 98")
add(EDITS_MS,
    "We therefore treat the published values as contextual reference points rather than a formal benchmark.",
    r"We therefore treat the published values as contextual reference points rather than a formal benchmark. Restricting DD to the three pairs with $\geq$30\% sequence identity (two positives) yielded AUROC~$=$~1.000; with so few pairs this estimate is anecdotal and is reported as an exploratory sensitivity analysis without a significance claim.",
    "exploratory 1.000 sentence")

# 8. evaluation-frameworks paragraph: Table 1 only
add(EDITS_MS,
    r"per-pair classifier features (Fig.~1c; Table~1).",
    r"per-pair classifier features (Table~1).",
    "CV3 ref line 232")

# 9. Fig 1 caption panel c: flowchart description replaces bar-chart description
add(EDITS_MS,
    r"""\textbf{c}, Contextual comparison: DD AUROC vs.\ published CV3 values from eight deep learning methods \cite{Feng2024}. Red bars, this study; blue bars, published values. \textit{Note: different test sets; not a head-to-head benchmark.} The interpretable composite score (AUROC~$=$~0.831) matches the best of four multi-feature classifiers tested under leave-one-pair-out cross-validation (SVM-RBF, 0.841; classifier range 0.240--0.841), and each DD-based nomination traces to a single measurable quantity (signed DD), enabling direct experimental design, whereas multi-feature classifiers produce opaque scores. The DD~+~identity-filter value (AUROC 1.000; Table~1) rests on a subset of 3 pairs with 2 positives and is anecdotal.""",
    r"""\textbf{c}, Benchmark selection and evaluation framework. Twelve curated pairs were evidence-tiered after verification of each primary citation (three Tier A, two Tier B, five Tier C, two comparators; Supplementary Table~S3); the Tier A~$\cup$~Tier B set (five pairs) constitutes the primary external benchmark. The evaluation frame spans the three gynecological lineages (Ovarian, Endometrial, Cervical) with $\geq$5 mutant and $\geq$5 wild-type cell lines per entry, yielding 110 driver$\times$paralog$\times$lineage entries (8 positive; 102 unlabeled controls). Performance is reported under three pre-specified aggregation frameworks (see ``Evaluation frameworks'').""",
    "Fig1 caption panel c -> flowchart")

# 10. Table 1 caption: drop DD+ID and Interpretability sentences
add(EDITS_MS,
    r"""\caption{\textbf{DD performance in context of published CV3 results (not a head-to-head benchmark).} Published values are CV3 (gene-pair isolation) AUROCs from Feng et al.\ (2024), Supplementary Data 1 (NSMRand negative sampling, 1:1 positive:negative ratio, complete dataset). DD and DD~+~ID values from this study: DD was evaluated on the full lineage-level frame (110 driver--paralog--lineage entries, 8 positives); DD~+~ID on the high-identity subset (3 pairs, 2 positives) and is reported as anecdotal. \textit{Evaluation frameworks differ: published methods were tested on general SL gene-pair universes; DD was tested on paralog-SL pairs only. Direct AUROC comparison is not warranted.} Interpretability: whether predictions can be directly traced to specific biological features.}""",
    r"""\caption{\textbf{DD performance in context of published CV3 results (not a head-to-head benchmark).} Published values are CV3 (gene-pair isolation) AUROCs from Feng et al.\ (2024), Supplementary Data 1 (NSMRand negative sampling, 1:1 positive:negative ratio, complete dataset). The DD value from this study was evaluated on the full lineage-level frame (110 driver--paralog--lineage entries, 8 positives). \textit{Evaluation frameworks differ: published methods were tested on general SL gene-pair universes with different datasets, tasks, and positive:negative ratios; DD was tested on paralog-SL pairs only. Direct AUROC comparison is not warranted.}}""",
    "Table 1 caption")

# 11. Table 1 header: drop Interpretability column
add(EDITS_MS,
    "\\begin{tabular}{lcc}\n\\toprule\n\\textbf{Method} & \\textbf{CV3 AUROC} & \\textbf{Interpretability} \\\\",
    "\\begin{tabular}{lc}\n\\toprule\n\\textbf{Method} & \\textbf{CV3 AUROC} \\\\",
    "Table 1 header")

# 12. delete the anecdotal DD+ID 1.000 row
add(EDITS_MS,
    "DD + ID $\\ge$ 30\\% (this study)$^{\\dagger}$ & 1.000 & High \\\\\n",
    "",
    "Table 1 DD+ID row")

# 13-21. strip " & High/Low" from the remaining nine rows
for row, val in [
    ("DD (this study)", "0.629 & High"),
    ("SLMGAE \\cite{Hao2021}", "0.790 & Low"),
    ("NSF4SL \\cite{Wang2022}", "0.683 & Low"),
    ("GCATSL \\cite{Long2021}", "0.678 & Low"),
    ("GRSMF \\cite{Huang2019}", "0.656 & Low"),
    ("PiLSL \\cite{Liu2022}", "0.626 & Low"),
    ("KG4SL \\cite{Wang2021}", "0.563 & Low"),
    ("SLGNN \\cite{Zhu2023}", "0.530 & Low"),
    ("PTGNN \\cite{Long2022}", "0.529 & Low"),
]:
    add(EDITS_MS,
        f"{row} & {val} \\\\",
        f"{row} & {val.split(' & ')[0]} \\\\",
        f"Table 1 row {row.split(' ')[0]}")

# 22. delete the dagger footnote
add(EDITS_MS,
    "\\begin{minipage}{\\textwidth}\\vspace{0.5em}\\footnotesize $^{\\dagger}$Anecdotal: high-identity subset of 3 pairs with 2 positives; reported without a significance claim (see text).\\end{minipage}\n",
    "",
    "Table 1 footnote")

# 23. subsection title: drop survival
add(EDITS_MS,
    r"\subsection*{DD signal across clinical contexts: MSI status, mutation type, and survival association}",
    r"\subsection*{DD signal across clinical contexts: MSI status and mutation type}",
    "subsection title")

# 24a. BRCA core sentence -> compressed exploratory summary pointing to S10
add(EDITS_MS,
    r"""TCGA Pan-Cancer Atlas survival analysis of breast cancer patients ($n=1{,}069$ with overall-survival and tumor mRNA expression data; 151 events) \cite{Liu2018,Cerami2012} tested whether high paralog expression tracks worse outcome. Cox proportional-hazards models were fitted on continuous (log$_2$-transformed, z-scored) expression for the 32 paralog genes of the 16 priority pairs. Higher ARID1B expression was associated with shortened overall survival (HR~$=$~1.30 per SD, 95\% CI: 1.087--1.552, $p=0.004$), and the association was robust to adjustment for age and AJCC stage (multivariable HR~$=$~1.31, 95\% CI: 1.087--1.569, $p=0.004$; proportional-hazards assumption satisfied, Schoenfeld $p=0.27$; Fig.~3c).""",
    r"""As an exploratory clinical correlation, we tested whether paralog expression tracks patient outcome in TCGA. In the BRCA cohort ($n=1{,}069$, 151 events), higher \textit{ARID1B} expression was associated with shortened overall survival (multivariable Cox HR~$=$~1.31 per SD, 95\% CI: 1.087--1.569, $p=0.004$, adjusted for age and AJCC stage; Supplementary Fig.~S10).""",
    "survival BRCA compressed")

# 24b. FDR detail + interpretation compressed
add(EDITS_MS,
    r"""Four of the 32 genes remained significant after Benjamini--Hochberg correction in the multivariable family (\textit{PIK3CA}, $q=0.035$; \textit{ARID1B}, \textit{BRCA2} and \textit{RBL1}, each $q=0.040$), whereas no association survived FDR correction in the univariable continuous family (smallest $q=0.129$, for \textit{ARID1B}; full results in Supplementary Table~S8). The signal is directionally consistent with the prioritization of ARID1A$\rightarrow$ARID1B as the leading candidate --- higher expression of the compensating paralog marks more aggressive disease --- but paralog expression alone is unlikely to serve as a strong standalone clinical biomarker, though it may contribute to multi-gene prognostic signatures. \textit{SMARCA2} --- the compensating paralog of the established SMARCA4$\rightarrow$SMARCA2 pair --- showed no association (multivariable HR~$=$~0.96), a reminder that expression-level prognosis and dependency-based synthetic lethality are distinct phenotypes.""",
    r"""Four of the 32 paralog genes passed Benjamini--Hochberg correction (full results in Supplementary Table~S8). \textit{SMARCA2} showed no association (HR~$=$~0.96), a reminder that expression-level prognosis and dependency-based synthetic lethality are distinct phenotypes; these associations are hypothesis-generating.""",
    "survival FDR compressed")

# 24c. UCEC/OV detail -> pointer (detail lives in supplementary extension)
add(EDITS_MS,
    r"""We extended the same pipeline to the two benchmark-relevant gynecological cohorts --- UCEC ($n=525$, 87 events) and OV ($n=299$, 180 events), where AJCC stage annotation is unpopulated and age-adjusted models are reported. The \textit{ARID1B} association replicated directionally in OV (univariable HR~$=$~1.26, 95\% CI: 1.092--1.461, $p=0.002$, BH $q=0.054$; age-adjusted HR~$=$~1.21, $p=0.013$) but not in UCEC (HR~$=$~1.17, 95\% CI: 0.941--1.464, $p=0.155$); within UCEC, the association was directionally consistent but non-significant in both the ARID1A-mutant (HR~$=$~1.24, 19 events) and wild-type (HR~$=$~1.15, 68 events) strata (Supplementary Information, TCGA survival extension).""",
    "Cohort-specific extensions in UCEC and OV are reported in the Supplementary TCGA survival extension.",
    "survival UCEC/OV pointer")

# 25. Fig 3 caption title
add(EDITS_MS,
    r"\caption{\textbf{Clinical stratification by MSI, mutation type, and survival.}",
    r"\caption{\textbf{Clinical stratification by MSI status and mutation type.}",
    "Fig3 caption title")

# 26. delete Fig 3 caption panel c (moved to S10)
add(EDITS_MS,
    r"""\textbf{c}, Forest plot of paralog gene expression vs.\ overall survival in TCGA PanCan Atlas BRCA ($n=1{,}069$; Cox proportional-hazards on continuous log$_2$-transformed, z-scored expression, adjusted for age and AJCC stage). Four genes withstand Benjamini--Hochberg correction across the 32-gene family (\textit{PIK3CA}, $q=0.035$; \textit{ARID1B}, \textit{BRCA2}, \textit{RBL1}, each $q=0.040$); these four are shown together with the compensating paralogs of the lead candidate pairs and \textit{ARID1A} for direct contrast. Error bars, 95\% CI. Full 32-gene results in Supplementary Table~S8.}""",
    "}",
    "Fig3 caption panel c removed")

# 27/28. PRISM anchor reference -> S7b; DWS panel reference -> Fig. 4
add(EDITS_MS,
    r"known drug-target biology was recovered (Fig.~4a):",
    r"known drug-target biology was recovered (Supplementary Fig.~S7b):",
    "Fig4a ref -> S7b")
add(EDITS_MS,
    r"selectivity tiers (Fig.~4b). Two pairs fell",
    r"selectivity tiers (Fig.~4). Two pairs fell",
    "Fig4b ref -> Fig4")

# 29. Fig 4 caption: single-panel DWS classification
add(EDITS_MS,
    r"""\caption{\textbf{Pharmacologic context and dependency-window classification.}
\textbf{a}, PRISM drug selectivity: |$\Delta$AUC| for top driver-genotype--drug sensitivity associations. Known drug-target biology is recapitulated: MEK inhibitors (AZD8330, Trametinib) selectively kill KRAS-mutant lines; mTOR/AKT inhibitors (Everolimus, Ipatasertib) kill PTEN-mutant lines; HDAC inhibitors (Panobinostat) kill EP300-mutant lines. BH $q<0.25$. The drugs shown are not necessarily direct binders of the paralog protein; the associations reflect driver-genotype-conditioned drug sensitivity, not paralog-specific targeting. See Supplementary Fig.~S7 for all associations passing the discovery-stage threshold.
\textbf{b}, Dependency-window classification. Bubble size, mean DWS; color, selectivity tier. HIGH\_SELECTIVITY requires selectivity~$>$~0.15 and DWS~$>$~1.0; NF1$\rightarrow$RASA2 has high DWS (5.42) but near-zero selectivity (0.005), hence MODERATE. See Supplementary Fig.~S8 for the full classification and Supplementary Fig.~S9 for the exploratory structure-derived descriptors.}""",
    r"""\caption{\textbf{Dependency-window classification.} Bubble size, mean DWS; color, selectivity tier ($n=21$ paralog pairs). HIGH\_SELECTIVITY requires selectivity~$>$~0.15 and DWS~$>$~1.0; NF1$\rightarrow$RASA2 has high DWS (5.42) but near-zero selectivity (0.005), hence MODERATE. See Supplementary Fig.~S8 for the full classification and Supplementary Fig.~S9 for the exploratory structure-derived descriptors; PRISM pharmacologic context is shown in Supplementary Fig.~S7.}""",
    "Fig4 caption single-panel")

# 30-32. supplementary figure list: S1-S9 -> S1-S10 + S10 description
add(EDITS_MS,
    "Supplementary Figures S1--S9 and Supplementary Tables S1--S11 are available",
    "Supplementary Figures S1--S10 and Supplementary Tables S1--S11 are available",
    "S-list count")
add(EDITS_MS,
    "domain conservation, and composite prioritization score. Supplementary Table~S1:",
    "domain conservation, and composite prioritization score. Supplementary Fig.~S10: Paralog expression and overall survival in TCGA BRCA (forest plot, multivariable Cox models). Supplementary Table~S1:",
    "S-list S10 entry")
add(EDITS_MS,
    "Supplementary Figures S1--S9, and Supplementary Tables S1--S11",
    "Supplementary Figures S1--S10, and Supplementary Tables S1--S11",
    "Additional file 1 count")

# ───────────────────────── supplementary.tex ─────────────────────────

# 33. S7 caption: two panels (heatmap + assay-validity anchors)
add(EDITS_SI,
    r"""\caption{\textbf{PRISM drug selectivity heatmap.} Heatmap of drug sensitivity differential ($\Delta$AUC = mean AUC in driver-mutant minus mean AUC in wild-type cell lines) for the top driver-genotype--drug sensitivity associations from the PRISM Repurposing screen (1,482 compounds, 727 cell lines). Negative $\Delta$AUC indicates selective sensitivity in driver-mutant cells. Associations shown met discovery-stage criteria (BH $q<0.25$, $\Delta$AUC$<0$, enrichment$>0.1$). The displayed top associations are dominated by cytotoxic chemotherapies (tubulin inhibitors, topoisomerase poisons, and mitotic-kinase inhibitors) with $\Delta$AUC$<0$, indicating genotype-conditioned cytotoxic sensitivity rather than paralog-specific targeting; the targeted-agent assay-validity anchors (MEK inhibitors--KRAS, mTOR/AKT inhibitors--PTEN, HDAC inhibitors--EP300) are shown in main-text Fig.~4a. Drugs shown are not necessarily direct binders of the paralog protein.}""",
    r"""\caption{\textbf{PRISM drug selectivity and assay-validity anchors.} \textbf{a}, Heatmap of drug sensitivity differential ($\Delta$AUC = mean AUC in driver-mutant minus mean AUC in wild-type cell lines) for the top driver-genotype--drug sensitivity associations from the PRISM Repurposing screen (1,482 compounds, 727 cell lines). Negative $\Delta$AUC indicates selective sensitivity in driver-mutant cells. Associations shown met discovery-stage criteria (BH $q<0.25$, $\Delta$AUC$<0$, enrichment$>0.1$). The displayed top associations are dominated by cytotoxic chemotherapies (tubulin inhibitors, topoisomerase poisons, and mitotic-kinase inhibitors) with $\Delta$AUC$<0$, indicating genotype-conditioned cytotoxic sensitivity rather than paralog-specific targeting. \textbf{b}, Assay-validity anchors: $|\Delta$AUC$|$ for targeted agents recapitulating known driver-genotype drug biology --- MEK inhibitors (AZD8330, Trametinib) selectively kill KRAS-mutant lines, mTOR/AKT inhibitors (Everolimus, Ipatasertib) kill PTEN-mutant lines, and the HDAC inhibitor Panobinostat kills EP300-mutant lines (BH $q<0.25$). These associations serve as an assay-validity control for the PRISM layer; they do not demonstrate paralog-specific targeting. Drugs shown are not necessarily direct binders of the paralog protein.}""",
    "S7 caption two-panel")

# 34. insert Fig S10 block at the end of the Supplementary Figures section
S10_BLOCK = r"""
% S10 (added 2026-08-01: TCGA BRCA survival forest, moved from main-text Fig. 3c)
\begin{figure}[H]
\centering
\includegraphics[width=0.85\textwidth,keepaspectratio]{output/figures/FigS10_Survival.pdf}
\caption{\textbf{Paralog expression and overall survival in TCGA BRCA.} Forest plot of paralog gene expression vs.\ overall survival in TCGA PanCan Atlas BRCA ($n=1{,}069$; Cox proportional-hazards on continuous log$_2$-transformed, z-scored expression, adjusted for age and AJCC stage). Four genes withstand Benjamini--Hochberg correction across the 32-gene family (\textit{PIK3CA}, $q=0.035$; \textit{ARID1B}, \textit{BRCA2}, \textit{RBL1}, each $q=0.040$); these four are shown together with the compensating paralogs of the lead candidate pairs and \textit{ARID1A} for direct contrast. Points, hazard ratio per SD; error bars, 95\% CI. Full 32-gene results in Supplementary Table~S8; UCEC and OV extensions in the Supplementary TCGA survival extension.}
\label{fig:survival_forest}
\end{figure}
"""
add(EDITS_SI,
    "\\end{figure}\n\n\\clearpage\n\\section{Supplementary Tables}",
    "\\end{figure}\n" + S10_BLOCK + "\n\\clearpage\n\\section{Supplementary Tables}",
    "insert Fig S10")

# 35. back-reference from the TCGA extension subsection to Fig S10
add(EDITS_SI,
    r"Machine-readable results: \texttt{output/tcga\_survival\_ucec.json}, \texttt{output/tcga\_survival\_ov.json}, and \texttt{output/tcga\_survival\_ucec\_ov\_summary.json} (\texttt{tcga\_survival\_ucec\_ov.py}).",
    r"Machine-readable results: \texttt{output/tcga\_survival\_ucec.json}, \texttt{output/tcga\_survival\_ov.json}, and \texttt{output/tcga\_survival\_ucec\_ov\_summary.json} (\texttt{tcga\_survival\_ucec\_ov.py}). The BRCA cohort forest plot is shown in Supplementary Fig.~S10.",
    "extension back-ref S10")


def apply(text, edits, name):
    for old, new, tag in edits:
        n = text.count(old)
        if n != 1:
            raise SystemExit(f"[{name}] anchor '{tag}' found {n} times (expected 1): {old[:80]!r}")
        text = text.replace(old, new, 1)
    return text


ms = apply(ms, EDITS_MS, "manuscript")
si = apply(si, EDITS_SI, "supplementary")

MS.write_text(ms, encoding="utf-8")
SI.write_text(si, encoding="utf-8")
print(f"OK: {len(EDITS_MS)} manuscript edits, {len(EDITS_SI)} supplementary edits applied.")
