#!/usr/bin/env python3
"""Humanize pass (blader/humanizer skill): purge em-dash habit, 'landscape',
'key' filler, '-ing' tails, 'notably/in particular' fillers.
Pure prose edits only — no numbers, no citations, no figure/table content.
Each replacement asserts exact-once occurrence. Table '---' cells in
supplementary (AUROC not computable) are intentionally untouched.
"""
import sys
from pathlib import Path

MS = Path("manuscript.tex")
SUP = Path("supplementary.tex")
CL = Path("cover_letter.md")

MS_FIX = [
    # ---- L69: method list appositive (literal —) -> parentheses ----
    ("the existing computational prediction methods — graph neural networks \\cite{Zhu2023,Long2021}, graph convolutional networks with dual dropout \\cite{Cai2020}, graph-regularized matrix factorization \\cite{Huang2019}, knowledge-graph neural networks \\cite{Wang2021}, and negative-sample-free contrastive learning \\cite{Wang2022} — train almost entirely",
     "the existing computational prediction methods (graph neural networks \\cite{Zhu2023,Long2021}, graph convolutional networks with dual dropout \\cite{Cai2020}, graph-regularized matrix factorization \\cite{Huang2019}, knowledge-graph neural networks \\cite{Wang2021}, and negative-sample-free contrastive learning \\cite{Wang2022}) train almost entirely"),
    # ---- L69: trailing dash -> comma ----
    ("the best published result reaches 0.790 \\cite{Hao2021} — precisely the setting where predictive models matter",
     "the best published result reaches 0.790 \\cite{Hao2021}, precisely the setting where predictive models matter"),
    # ---- L71: definition dash -> colon ----
    ("Gene duplication events leave behind paralogs — sequence-similar gene pairs",
     "Gene duplication events leave behind paralogs: sequence-similar gene pairs"),
    # ---- L75: landscape #1 ----
    ("That study established the landscape but did not benchmark",
     "That study opened the field but did not benchmark"),
    # ---- L98: classifier list appositive -> namely ----
    ("we trained four standard classifiers --- logistic regression (LR), random forest (RF), SVM-RBF, and SVM-Linear --- on five features",
     "we trained four standard classifiers, namely logistic regression (LR), random forest (RF), SVM-RBF, and SVM-Linear, on five features"),
    # ---- L98: measured quantity appositive -> commas ----
    ("maps to a single measured quantity --- a dependency shift --- which makes it straightforward",
     "maps to a single measured quantity, the dependency shift, which makes it straightforward"),
    # ---- L100: 'in particular' filler ----
    ("PCS in particular carries substantial signal",
     "Of the components, PCS carries substantial signal"),
    # ---- L102: appositive -> parentheses ----
    ("both evaluable positive entries --- the same pair, ARID1A$\\rightarrow$ARID1B, in Endometrial and Ovarian cancers --- ranked above",
     "both evaluable positive entries (the same pair, ARID1A$\\rightarrow$ARID1B, in Endometrial and Ovarian cancers) ranked above"),
    # ---- L104: power parenthesis -> parens (inner '(Methods)' -> 'see Methods') ----
    ("the evaluation is underpowered --- at the observed 0.629, power to reject a null of 0.5 is 0.34, and $\\ge$45 validated positive pairs would be required for 80\\% power at $\\alpha=0.05$ (Methods) --- so all benchmark AUROCs",
     "the evaluation is underpowered (at the observed 0.629, power to reject a null of 0.5 is 0.34, and $\\ge$45 validated positive pairs would be required for 80\\% power at $\\alpha=0.05$; see Methods), so all benchmark AUROCs"),
    # ---- L106: -ing tail -> split sentence ----
    ("reported by both studies; the single directionally consistent case",
     "reported by both studies. The single directionally consistent case"),
    ("($p=0.475$), reflecting the rarity of our driver mutations among screened lines and the low genotype penetrance of digenic interactions",
     "($p=0.475$). This reflects the rarity of our driver mutations among screened lines and the low genotype penetrance of digenic interactions"),
    # ---- L129: dash -> colon ----
    ("below DD (0.629) --- a co-variation-versus-discrimination distinction",
     "below DD (0.629): a co-variation-versus-discrimination distinction"),
    # ---- L177: 'key' filler in Fig. 3 caption ----
    ("key TSG paralog-SL pairs",
     "leading TSG paralog-SL pairs"),
    # ---- L192: metric appositive -> commas ----
    ("a single, interpretable metric --- the dependency shift between driver-mutant and wild-type cells --- could capture",
     "a single, interpretable metric, the dependency shift between driver-mutant and wild-type cells, could capture"),
    # ---- L194: 'notably' filler ----
    ("co-regulation at both layers; notably, Venkatesh et al.\\",
     "co-regulation at both layers; Venkatesh et al.\\"),
    # ---- L194: however-dash pair -> colon + comma ----
    ("Matched mRNA co-variation is also present, however --- for EP300$\\leftrightarrow$CREBBP the mRNA correlation exceeded the protein correlation in every cohort --- and across the 122",
     "Matched mRNA co-variation is also present, however: for EP300$\\leftrightarrow$CREBBP the mRNA correlation exceeded the protein correlation in every cohort, and across the 122"),
    # ---- L196: dash after parenthesis -> period ----
    ("in Breast cancer) --- cells carrying an activated \\textit{PIK3CA} allele",
     "in Breast cancer). Cells carrying an activated \\textit{PIK3CA} allele"),
    # ---- L198: DWS parenthesis -> parens, de-nest '(4.87)' ----
    ("as the highest-selectivity established candidate --- SMARCA4$\\rightarrow$SMARCA2 had the higher DWS (4.87) --- while flagging",
     "as the highest-selectivity established candidate (SMARCA4$\\rightarrow$SMARCA2 had the higher DWS, 4.87) while flagging"),
    # ---- L198: medicinal chemistry dash -> comma ----
    ("extensive medicinal chemistry --- none of which our computational targetability descriptors",
     "extensive medicinal chemistry, none of which our computational targetability descriptors"),
    # ---- L200: literal — -> comma ----
    ("in the CV3 setting — a comparison their study did not perform",
     "in the CV3 setting, a comparison their study did not perform"),
    # ---- L200: echoing dash -> comma ----
    ("hemizygous copy-number loss --- echoing the CYCLOPS paradigm",
     "hemizygous copy-number loss, echoing the CYCLOPS paradigm"),
    # ---- L200: landscape #2 ----
    ("extend this landscape at much larger scale",
     "extend this map at much larger scale"),
    # ---- L204: appositive -> commas ----
    ("the within-driver FDR test --- conditional association rather than family-wise enrichment --- and we do not use",
     "the within-driver FDR test, conditional association rather than family-wise enrichment, and we do not use"),
    # ---- L206: dash -> comma ----
    ("(NCI-H1299, MDA-MB-231) --- consistent with the low genotype penetrance",
     "(NCI-H1299, MDA-MB-231), consistent with the low genotype penetrance"),
    # ---- L206: 'Notably,' sentence opener ----
    ("Notably, this pair was included in the unbiased digenic screen",
     "This pair was also included in the unbiased digenic screen"),
    # ---- L208: benchmark universe -> colon ----
    ("all on the same benchmark universe --- the twelve curated pairs plus matched unlabeled controls",
     "all on the same benchmark universe: the twelve curated pairs plus matched unlabeled controls"),
    # ---- L208: power limitation -> comma ----
    ("the permutation null ($p=0.223$) --- a power limitation",
     "the permutation null ($p=0.223$), a power limitation"),
    # ---- L214: housekeeping appositive (literal —) -> commas ----
    ("gold-standard pairs \\cite{EsmaeiliAnvar2024} — predominantly housekeeping genes whose alterations are mostly passenger events — did not separate",
     "gold-standard pairs \\cite{EsmaeiliAnvar2024}, predominantly housekeeping genes whose alterations are mostly passenger events, did not separate"),
    # ---- L221: units dash -> colon ----
    ("computed at different units across analyses --- driver$\\times$paralog$\\times$lineage entries",
     "computed at different units across analyses: driver$\\times$paralog$\\times$lineage entries"),
    # ---- L222: landscape #3 ----
    ("are landscape descriptors rather than independent benchmarks",
     "are descriptive rather than independent benchmarks"),
    # ---- L228: candidates appositive -> commas ----
    ("unlabeled candidates --- starting with KMT2D$\\rightarrow$KMT2C in KMT2D-mutant versus wild-type isogenic lines --- with viability",
     "unlabeled candidates, starting with KMT2D$\\rightarrow$KMT2C in KMT2D-mutant versus wild-type isogenic lines, with viability"),
    # ---- L236: dependencies appositive -> commas ----
    ("mutation-conditioned dependencies --- including ARID1A$\\rightarrow$ARID1B, the highest-selectivity established candidate --- and prioritized",
     "mutation-conditioned dependencies, including ARID1A$\\rightarrow$ARID1B, the highest-selectivity established candidate, and prioritized"),
    # ---- L238: toolkit appositive (literal —) -> comma + which ----
    ("a complete open-source toolkit — the \\texttt{paralogSL} R package and a reproducible Python pipeline — that requires only",
     "a complete open-source toolkit, the \\texttt{paralogSL} R package and a reproducible Python pipeline, which together require only"),
    # ---- L258: 'key' filler in Methods ----
    ("For each key driver--paralog pair",
     "For each priority driver--paralog pair"),
    # ---- L266: independence dash -> colon ----
    ("differ in evidence independence --- Tier A derives from dual-gene perturbation",
     "differ in evidence independence: Tier A derives from dual-gene perturbation"),
    # ---- L302: thresholds list -> colon + comma ----
    ("stratification thresholds --- PAN\\_ESSENTIAL", "stratification thresholds: PAN\\_ESSENTIAL"),
    ("or LOW\\_SELECTIVITY --- which are not validated clinical cutoffs",
     "or LOW\\_SELECTIVITY, which are not validated clinical cutoffs"),
    # ---- L306: dash -> semicolon ----
    ("over the sequence features above --- it is not a three-dimensional alignment metric",
     "over the sequence features above; it is not a three-dimensional alignment metric"),
    # ---- L314: natural-unit appositive -> commas ----
    ("driver$\\times$paralog$\\times$lineage level --- the natural unit at which DD is computed --- with the pair treated",
     "driver$\\times$paralog$\\times$lineage level, the natural unit at which DD is computed, with the pair treated"),
    # ---- L318: lineage list -> parens ----
    ("the three evaluation lineages --- Ovarian, Endometrial, Cervical --- passing the minimum sample size filter",
     "the three evaluation lineages (Ovarian, Endometrial, Cervical) passing the minimum sample size filter"),
    # ---- L318 (supplementary-info heading): landscape #4 ----
    ("Evaluation landscape across 23 solid tumor types",
     "Evaluation overview across 23 solid tumor types"),
    # ---- L375: audit list -> colon; tail -> period split (literal —) ----
    ("derived artifacts — headline metrics, classifier benchmarks",
     "derived artifacts: headline metrics, classifier benchmarks"),
    ("405 claims in total — each ending with an automated claim-by-claim comparison",
     "405 claims in total. Each script ends with an automated claim-by-claim comparison"),
]

SUP_FIX = [
    # ---- L46: Tier B appositive -> commas ----
    ("natural-genotype conditional dependency --- the same data type as DD --- so Tier~B concordance",
     "natural-genotype conditional dependency, the same data type as DD, so Tier~B concordance"),
    # ---- L46: Tier C list -> colon + semicolon ----
    ("Tier~C comprises indirect evidence only --- reciprocal-direction-only evidence",
     "Tier~C comprises indirect evidence only: reciprocal-direction-only evidence"),
    ("(FBXW7$\\rightarrow$FBXW2, PPP2R1A$\\rightarrow$PPP2R1B) --- and is excluded from the primary benchmark",
     "(FBXW7$\\rightarrow$FBXW2, PPP2R1A$\\rightarrow$PPP2R1B); it is excluded from the primary benchmark"),
    # ---- L52: lineage list -> parens (distinguish from L215 by tail) ----
    ("the three evaluation lineages --- Ovarian, Endometrial, Cervical --- across the 12 curated pairs",
     "the three evaluation lineages (Ovarian, Endometrial, Cervical) across the 12 curated pairs"),
    # ---- L55: classifier list -> namely ----
    ("four standard classifiers --- logistic regression", "four standard classifiers, namely logistic regression"),
    ("SVM with linear kernel --- were trained on five features", "SVM with linear kernel, were trained on five features"),
    # ---- L67: numerator appositive -> commas ----
    ("$\\max(\\mathrm{DD},0)$ --- which credits only compensation-direction shifts --- to a conservative",
     "$\\max(\\mathrm{DD},0)$, which credits only compensation-direction shifts, to a conservative"),
    # ---- L104: Fig S1 caption landscape #5 ----
    ("\\textbf{Evaluation landscape across 23 solid tumor types.}",
     "\\textbf{Evaluation overview across 23 solid tumor types.}"),
    # ---- L215: lineage list -> parens ----
    ("the three evaluation lineages --- Ovarian, Endometrial, Cervical --- passing the minimum sample-size filter",
     "the three evaluation lineages (Ovarian, Endometrial, Cervical) passing the minimum sample-size filter"),
    # ---- L262: Tier B caption appositive -> i.e. ----
    ("natural-genotype conditional dependency --- the paralog is a demonstrated selective dependency in driver-mutant cells --- from single-gene perturbation",
     "natural-genotype conditional dependency, i.e., the paralog is a demonstrated selective dependency in driver-mutant cells, from single-gene perturbation"),
    # ---- L462: wild-type dash -> comma ----
    ("in wild-type --- directionally consistent in both strata",
     "in wild-type, directionally consistent in both strata"),
]

CL_FIX = [
    # ---- L20: screen names appositive -> colon, move to end ----
    ("Two combinatorial CRISPR screens published in 2025 — Harle et al. (*Genome Biology*; 472 nominated pairs) and Flister et al. (*Cell Reports*; the near-complete paralogome, 36,648 pairs) — were analysed as fully external tests that played no role in method development or candidate selection.",
     "Two combinatorial CRISPR screens published in 2025 were analysed as fully external tests that played no role in method development or candidate selection: Harle et al. (*Genome Biology*; 472 nominated pairs) and Flister et al. (*Cell Reports*; the near-complete paralogome, 36,648 pairs)."),
    # ---- L20: at-scale appositive -> commas ----
    ("confirming at scale — and on experimental platforms independent of DepMap — that dependency displacement",
     "confirming at scale, and on experimental platforms independent of DepMap, that dependency displacement"),
    # ---- L22: candidate dash -> colon ----
    ("as the highest-selectivity established candidate — a hypothesis-generating experimental validation strategy",
     "as the highest-selectivity established candidate: a hypothesis-generating experimental validation strategy"),
    # ---- L50-56: reviewer list name — affiliation -> comma ----
    ("Dr. Francisca Vazquez — DepMap / Broad Institute", "Dr. Francisca Vazquez, DepMap / Broad Institute"),
    ("Dr. Michael P. Snyder — Stanford University", "Dr. Michael P. Snyder, Stanford University"),
    ("Dr. Jason Moffat — University of Toronto", "Dr. Jason Moffat, University of Toronto"),
    ("Dr. Rameen Beroukhim — Dana-Farber / Broad Institute", "Dr. Rameen Beroukhim, Dana-Farber / Broad Institute"),
    ("Dr. Bing Zhang — Baylor College of Medicine", "Dr. Bing Zhang, Baylor College of Medicine"),
    ("Dr. G. Traver Hart — Department of Systems Biology", "Dr. G. Traver Hart, Department of Systems Biology"),
    ("Dr. Min Wu — Institute for Infocomm Research", "Dr. Min Wu, Institute for Infocomm Research"),
]


def apply(path, fixes):
    text = path.read_text(encoding="utf-8")
    n_ok = 0
    for old, new in fixes:
        c = text.count(old)
        if c != 1:
            print(f"FAIL [{path.name}]: occurrence count {c} (expected 1)\n  OLD: {old[:110]}")
            return False
        text = text.replace(old, new, 1)
        n_ok += 1
    path.write_text(text, encoding="utf-8")
    print(f"OK [{path.name}]: {n_ok} fixes applied")
    return True


# sup/CL already applied in the previous run; MS failed before write, rerun MS only.
ok = apply(MS, MS_FIX)
sys.exit(0 if ok else 1)
