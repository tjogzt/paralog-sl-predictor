#!/usr/bin/env python3
"""D1 batch: coordinated text edits to manuscript.tex (round-4 review, items 4.1/4.3/4.5).
Atomic: all replacements verified before a single write."""
from pathlib import Path

p = Path("manuscript.tex")
t = p.read_text()

REPS = [
    # ── 1. Intro line 73: add De Kegel / in4mer framing after D'Antonio citation ──
    ("Paralog compensation provides a simpler route into this problem \\cite{Koonin2005,DAntonio2013}.",
     "Paralog compensation provides a simpler route into this problem \\cite{Koonin2005,DAntonio2013}. "
     "Large-scale analyses of CRISPR dependency data have since established paralog pairs as the most "
     "enriched class of synthetic-lethal interactions and shown that interpretable models built on "
     "protein-interaction and evolutionary-conservation features can predict them \\cite{DeKegel2021}; "
     "combinatorial knockout screens have begun to map paralog synthetic lethality directly "
     "\\cite{Parrish2021,EsmaeiliAnvar2024}."),
    # ── 2. Intro line 73: >30% sequence identity -> proxy wording ──
    ("This is \\textit{sequence-paralog SL}: synthetic lethality between two genes that share $>$30\\% "
     "sequence identity and conserved functional domains.",
     "This is \\textit{sequence-paralog SL}: synthetic lethality between two genes that share high protein "
     "sequence identity (operationalized with a validated $k$-mer identity-enrichment proxy; Methods) and "
     "conserved functional domains."),
    # ── 3. Intro line 75: BRCA1/BRCA2 not "negative controls" ──
    ("Our analysis concentrates on sequence-paralog pairs; functional analogs including BRCA1/BRCA2 are "
     "retained as negative controls and flagged throughout (Supplementary Table~S3).",
     "Our analysis concentrates on sequence-paralog pairs; functional analogs including BRCA1/BRCA2 are "
     "retained as mechanistic comparators (specificity references, not verified negatives) and flagged "
     "throughout (Supplementary Table~S3)."),
    # ── 4. Intro line 79: PCAWG -> COSMIC CGC ──
    ("using a 40-gene driver panel (Supplementary Table~S5) drawn from TCGA \\cite{Bailey2018} and PCAWG "
     "driver catalogs \\cite{ICGC2020}",
     "using a 40-gene driver panel (Supplementary Table~S5) curated from TCGA PanCanAtlas driver analyses "
     "\\cite{Bailey2018} and the COSMIC Cancer Gene Census \\cite{Sondka2018}"),
    # ── 5. Results line 100: sequence-identity filter -> k-mer proxy ──
    ("A DD~+~sequence-identity ($\\ge$30\\%) filter achieved AUROC~$=$~1.000 on a high-identity subset of "
     "3 pairs (2 known SL)",
     "A DD~+~$k$-mer identity-enrichment filter ($\\ge$30\\% proxy) achieved AUROC~$=$~1.000 on a "
     "high-identity subset of 3 pairs (2 known SL)"),
    # ── 6. Discussion line 248: broaden to prior work (D'Antonio + De Kegel + in4mer) ──
    ("\\textbf{Relationship to D'Antonio et al.} D'Antonio et al.\\ \\cite{DAntonio2013} published the "
     "foundational analysis of paralog genetic interactions. Our work extends theirs in three ways.",
     "\\textbf{Relationship to prior work.} D'Antonio et al.\\ \\cite{DAntonio2013} published the "
     "foundational analysis of paralog genetic interactions, and De Kegel et al.\\ \\cite{DeKegel2021} "
     "systematically predicted robust paralog synthetic lethality across hundreds of cancer cell lines "
     "using interpretable features (protein interaction, evolutionary conservation). Combinatorial CRISPR "
     "platforms have since mapped paralog synthetic lethals experimentally \\cite{Parrish2021,EsmaeiliAnvar2024}, "
     "converging on 13 candidate cross-study gold standards \\cite{EsmaeiliAnvar2024} whose replication "
     "across screens remains imperfect \\cite{Chou2025}. Our work is complementary to these efforts in "
     "three ways."),
    # ── 7. Discussion line 248 continuation: "Together, the five orthogonal evidence layers..." keep ──
    # ── 8. Methods line 303: new tier system ──
    ("\\textit{Tier A} (two pairs: SMARCA4$\\rightarrow$SMARCA2, ARID1A$\\rightarrow$ARID1B) has "
     "directional experimental evidence independent of DepMap matching the direction scored here. One pair "
     "(EP300/CREBBP) is experimentally established as synthetic lethal only in the reciprocal direction "
     "(CREBBP$\\rightarrow$EP300 \\cite{Ogiwara2016,Nie2021}); it is counted as a known pair at the pair "
     "level but excluded from directional Tier A claims and relabelled as non-positive in the "
     "direction-strict sensitivity analysis. \\textit{Tier B} (five pairs) is supported by "
     "paralog-redundancy, digenic-knockout, or pharmacologic evidence only, and \\textit{Tier C} (two "
     "pairs) derives from DepMap analyses. Two mechanistic comparators (BRCA1$\\leftrightarrow$BRCA2 "
     "\\cite{Bryant2005} and STK11$\\rightarrow$SIK1, an LKB1--SIK axis reference \\cite{Hollstein2019}) "
     "are reported separately.",
     "\\textit{Tier A} (three pairs: AKT1$\\rightarrow$AKT2 \\cite{Najm2018}, CDK4$\\rightarrow$CDK6 and "
     "MAP2K1$\\rightarrow$MAP2K2 \\cite{Parrish2021}) comprises pairs with direct genetic synthetic-lethal "
     "evidence from dual-gene perturbation (combinatorial or digenic CRISPR knockout). \\textit{Tier B} "
     "(two pairs: SMARCA4$\\rightarrow$SMARCA2 \\cite{Hoffman2014}, ARID1A$\\rightarrow$ARID1B "
     "\\cite{Helming2014}) comprises pairs for which the paralog is a demonstrated selective dependency in "
     "driver-mutant cells (natural-genotype conditional dependency from single-gene perturbation with "
     "functional validation). The Tier A~$\\cup$~Tier B set (five pairs) constitutes the primary external "
     "benchmark. \\textit{Tier C} (five pairs) comprises indirect evidence only and is excluded from the "
     "primary benchmark: EP300$\\rightarrow$CREBBP, established experimentally in the reciprocal direction "
     "only (CREBBP$\\rightarrow$EP300 \\cite{Ogiwara2016,Nie2021}); PIK3CA$\\rightarrow$PIK3CB, whose "
     "primary evidence supports PTEN$\\rightarrow$PIK3CB \\cite{Wee2008}; CCNE1$\\rightarrow$CCNE2, "
     "supported by developmental redundancy in mouse double knockouts \\cite{Geng2003}; and "
     "FBXW7$\\rightarrow$FBXW2 and PPP2R1A$\\rightarrow$PPP2R1B, derived from DepMap analyses. Two "
     "mechanistic comparators (BRCA1$\\leftrightarrow$BRCA2 \\cite{Bryant2005} and STK11$\\rightarrow$SIK1, "
     "an LKB1--SIK axis reference \\cite{Hollstein2019}) are not sequence paralogs and serve as specificity "
     "references only."),
    # ── 9. Methods line 303: AUROC-inflate error -> correct prevalence statement ──
    ("This yields an imbalanced evaluation set (12 positives vs.\\ $\\sim$194 negatives per lineage), which "
     "is representative of the real-world class imbalance but may inflate AUROC; we therefore also report "
     "AUPRC where applicable.",
     "The evaluation is highly imbalanced (12 curated positives vs.\\ $\\sim$194 unlabeled controls per "
     "lineage). AUROC is insensitive to class prevalence, but precision--recall performance is not; we "
     "therefore report AUPRC (average precision) alongside AUROC, with bootstrap confidence intervals for "
     "both."),
    # ── 10. Methods line 323: ICGC2020 -> Sondka2018 for COSMIC CGC ──
    ("Forty pan-cancer driver genes were selected from TCGA pan-cancer analyses \\cite{Bailey2018} and the "
     "COSMIC Cancer Gene Census \\cite{ICGC2020}",
     "Forty pan-cancer driver genes were selected from TCGA PanCanAtlas driver analyses \\cite{Bailey2018} "
     "and the COSMIC Cancer Gene Census \\cite{Sondka2018}"),
    # ── 11. Methods line 347: tier description + primary endpoint + power ──
    ("Gold-standard pairs were stratified into evidence tiers (Supplementary Table~S3): a \\textit{Tier A} "
     "set of 2 pairs with directional external experimental evidence, one reciprocal-validated pair "
     "(EP300/CREBBP; direct evidence for CREBBP$\\rightarrow$EP300 only), \\textit{Tier B} (5 pairs; "
     "paralog-redundancy, digenic-knockout, or pharmacologic evidence only, not directional "
     "genotype-conditional), \\textit{Tier C} (2 pairs; DepMap-derived), and 2 mechanistic comparators "
     "(BRCA1$\\leftrightarrow$BRCA2, STK11$\\rightarrow$SIK1). Main AUROC results report the full curated "
     "set (12 pairs); tier-restricted and direction-strict analyses are reported as sensitivity analyses. "
     "With at most 2--3 strictly validated positive pairs, subset-level AUROC has limited power: at "
     "$\\alpha=0.05$, an estimated $\\ge$25 validated positive pairs would be needed to achieve 80\\% power "
     "for detecting AUROC~$=$~0.85 against a null of 0.5.",
     "Gold-standard pairs were stratified into evidence tiers (Supplementary Table~S3): \\textit{Tier A} "
     "(3 pairs; direct dual-perturbation genetic synthetic-lethal evidence), \\textit{Tier B} (2 pairs; "
     "natural-genotype conditional dependency with functional validation), \\textit{Tier C} (5 pairs; "
     "indirect, reciprocal-direction-only, or DepMap-derived evidence), and 2 mechanistic comparators "
     "(BRCA1$\\leftrightarrow$BRCA2, STK11$\\rightarrow$SIK1). The pre-specified primary endpoint is the "
     "AUROC and AUPRC of $|$DD$|$ on the Tier A~$\\cup$~Tier B external benchmark at the "
     "driver$\\times$paralog$\\times$lineage level --- the natural unit at which DD is computed --- with "
     "the pair treated as the grouping unit in sensitivity analyses to avoid cross-lineage information "
     "leakage. The full 12-pair curated set, tier-restricted, direction-strict, and per-pair aggregated "
     "analyses are reported as secondary and sensitivity analyses. With at most 3--5 evaluable positive "
     "pairs, subset-level AUROC has limited power: at $\\alpha=0.05$, an estimated $\\ge$25 validated "
     "positive pairs would be needed to achieve 80\\% power for detecting AUROC~$=$~0.85 against a null of "
     "0.5."),
    # ── 12. Methods line 299: min 3 -> 5 + Hedges g + sensitivity ──
    ("DD was computed as defined in Equation~\\ref{eq:dd}, using Welch's $t$-test (unequal variances) with "
     "Cohen's $d$ as the standardized effect size (pooled SD). Minimum three mutant and three wild-type "
     "cell lines per driver.",
     "DD was computed as defined in Equation~\\ref{eq:dd}, using Welch's $t$-test (unequal variances); "
     "effect sizes are reported as Cohen's $d$ (pooled SD) and the small-sample-corrected Hedges' $g$. The "
     "primary analysis required a minimum of five mutant and five wild-type cell lines per "
     "driver$\\times$lineage stratum; the more permissive $\\ge$3/$\\ge$3 threshold was re-run as a "
     "sensitivity analysis (Results)."),
    # ── 13. Methods line 307: bootstrap/AUPRC/paired bootstrap ──
    ("AUROC was computed using scikit-learn 1.3.0 \\cite{Pedregosa2011}. Bootstrap 95\\% CI: 1,000 "
     "iterations, resampling paralog pairs with replacement. For the DD~+~ID$\\ge$30\\% subset (3 pairs), "
     "bootstrap CI was computed separately to account for the reduced sample size. Negative controls: 100 "
     "permutations of shuffled known-SL labels.",
     "AUROC was computed with the rank-based Mann--Whitney statistic (identical to scikit-learn 1.3.0 "
     "\\cite{Pedregosa2011}); AUPRC is reported as average precision. Bootstrap 95\\% CIs used 10,000 "
     "resamples of evaluation entries (or pairs, for per-pair frames) with replacement; paired bootstrap "
     "resamples were used for head-to-head comparison of $|$DD$|$ against component scores on identical "
     "resamples. For the DD~+~identity-proxy subset (3 pairs), bootstrap CI was computed separately to "
     "account for the reduced sample size. Negative controls: 100 permutations of shuffled known-SL "
     "labels."),
    # ── 14. Limitations line 258: min threshold wording ──
    ("\\item \\textit{Sample size:} Several rare driver mutations appear in fewer than three DepMap cell "
     "lines and could not be analyzed;",
     "\\item \\textit{Sample size:} Several rare driver mutations appear in fewer than five DepMap cell "
     "lines and could not be analyzed in the primary framework;"),
    # ── 15. Bibliography: remove orphaned ICGC2020 ──
    ("\\bibitem{ICGC2020} ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium. Pan-cancer analysis "
     "of whole genomes. \\textit{Nature}. 2020;578(7793):82--93.\n",
     ""),
    # ── 16. Bibliography: add new refs before end ──
    ("\\bibitem{Nie2021} Nie M, Du L, Ren W, Joung J, Ye X, Shi X, et al.\\ Genome-wide CRISPR screens "
     "reveal synthetic lethal interaction between CREBBP and EP300 in diffuse large B-cell lymphoma. "
     "\\textit{Cell Death Dis}. 2021;12(5):419.",
     "\\bibitem{Nie2021} Nie M, Du L, Ren W, Joung J, Ye X, Shi X, et al.\\ Genome-wide CRISPR screens "
     "reveal synthetic lethal interaction between CREBBP and EP300 in diffuse large B-cell lymphoma. "
     "\\textit{Cell Death Dis}. 2021;12(5):419.\n\n"
     "\\bibitem{Sondka2018} Sondka Z, Bamford S, Cole CG, Ward SA, Beare DM, Gunasekaran P, et al.\\ The "
     "COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. \\textit{Nat Rev "
     "Cancer}. 2018;18(11):696--705.\n\n"
     "\\bibitem{DeKegel2021} De Kegel B, Quinn N, Thompson NA, Adams DJ, Ryan CJ. Comprehensive prediction "
     "of robust synthetic lethality between paralog pairs in cancer cell lines. \\textit{Cell Syst}. "
     "2021;12(12):1144--1159.e6.\n\n"
     "\\bibitem{EsmaeiliAnvar2024} Esmaeili Anvar N, Lin C, Ma X, Wilson LL, Steger R, Sangree AK, et al.\\ "
     "Efficient gene knockout and genetic interaction screening using the in4mer CRISPR/Cas12a multiplex "
     "knockout platform. \\textit{Nat Commun}. 2024;15:3577.\n\n"
     "\\bibitem{Chou2025} Chou J, Hart T. Z-scores outperform similar methods for analyzing CRISPR paralog "
     "synthetic lethality screens. \\textit{Genome Biol}. 2025;26(1):188."),
]

missing = []
for old, new in REPS:
    n = t.count(old)
    if n != 1:
        missing.append((n, old[:90]))
    else:
        t = t.replace(old, new)

if missing:
    for n, s in missing:
        print(f"MATCH COUNT {n}: {s}")
    raise SystemExit("ABORT: not all anchors unique; no write performed")

p.write_text(t)
print(f"OK: {len(REPS)} replacements applied to {p}")
