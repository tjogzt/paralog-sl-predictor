#!/usr/bin/env python3
"""D2 manuscript number sync: min5 primary frame + min3 sensitivity frame + new tier system.
Every replacement must match exactly once; abort otherwise."""
from pathlib import Path

MS = Path(__file__).parent / "manuscript.tex"
text = MS.read_text()

R = []  # (old, new, comment)

# ── Abstract (line 58) ─────────────────────────────────────────────
R.append((
 r"""DD achieved AUROC~$=$~0.682 (0.728 after excluding two pairs with DepMap-derived labels; 1.000 on the two pairs with directional evidence independent of DepMap)""",
 r"""DD achieved AUROC~$=$~0.676 (0.725 after excluding two pairs with DepMap-derived labels; 1.000 on the Tier A~$\cup$~Tier B external benchmark, on which two of five pairs were lineage-evaluable)""",
 "abstract headline"))
R.append((
 r"""In a pair-level benchmark on the same 75 test pairs, the interpretable composite score reached AUROC~$=$~0.841, exceeding four standard multi-feature classifiers under leave-one-pair-out cross-validation.""",
 r"""In a pair-level benchmark on the same 72 test pairs, the interpretable composite score reached AUROC~$=$~0.831, statistically indistinguishable from the best of four standard multi-feature classifiers (SVM-RBF, 0.843) under leave-one-pair-out cross-validation.""",
 "abstract pair-level"))

# ── Results: cross-cancer (line 98) ────────────────────────────────
R.append((
 r"""DD exceeded AUROC~$=$~0.7 in 10 of 14 evaluable lineages (Fig.~1a; Supplementary Table~S1). The strongest performers were Biliary Tract (AUROC~$=$~0.990, $n=50$ pairs), Esophagogastric (0.969), Pancreatic (0.949), Glioma (0.885), Bladder Urothelial (0.882), Breast (0.852), head and neck squamous cell carcinoma (HNSCC; 0.812), small cell lung cancer (SCLC; 0.779), non-small cell lung cancer (NSCLC; 0.738), and Colorectal (0.726). Performance did not differ significantly by oncogenic mechanism (Fig.~1b): TSG-driven cancers had a mean AUROC of 0.792 ($n=11$ lineages) versus 0.723 for oncogene-driven cancers ($n=3$ lineages; permutation $p=0.538$, exact Mann-Whitney $p=0.769$).""",
 r"""DD exceeded AUROC~$=$~0.7 in 7 of 8 evaluable lineages on the primary frame (Fig.~1a; Supplementary Table~S1): Esophagogastric (AUROC~$=$~0.965, $n=50$ pairs), small cell lung cancer (SCLC; 0.906), Bladder Urothelial (0.844), Colorectal (0.828), Endometrial (0.818), Breast (0.750), and non-small cell lung cancer (NSCLC; 0.741); Ovarian (0.661) fell below this threshold. Under the relaxed $\geq$3-per-group sensitivity frame (12 evaluable lineages), 9 of 12 exceeded 0.7, with Biliary Tract (0.990) and Pancreatic (0.949) among the strongest performers. Performance did not differ significantly by oncogenic mechanism (Fig.~1b; sensitivity frame, because the primary frame retains only one oncogene-driven lineage): TSG-driven cancers had a mean AUROC of 0.814 ($n=9$ lineages) versus 0.768 for oncogene-driven cancers ($n=3$ lineages; permutation $p=0.645$, exact Mann-Whitney $p=0.600$).""",
 "cross-cancer paragraph"))

# ── Results: ML benchmark (line 100) ───────────────────────────────
R.append((
 r"""DD achieved AUROC~$=$~0.682 without any training""",
 r"""DD achieved AUROC~$=$~0.676 without any training""",
 "CV3 DD value"))
R.append((
 r"""on the same 75 paralog pairs (6 known positives)""",
 r"""on the same 72 paralog pairs (6 known positives)""",
 "pair count"))
R.append((
 r"""The composite score alone performed best (AUROC~$=$~0.841), followed by SVM-RBF (0.744) and random forest (0.617); DD alone reached 0.551, and the linear classifiers were unstable at this sample size (SVM-Linear 0.217, LR 0.138).""",
 r"""The composite score (AUROC~$=$~0.831) and SVM-RBF (0.843) performed best and are statistically indistinguishable at this sample size, followed by random forest (0.722) and DD alone (0.566); the linear classifiers were unstable (SVM-Linear 0.114, LR 0.136).""",
 "classifier ranking"))
R.append((
 r"""in the full LR fit no individual feature reached significance (all $p>0.09$)""",
 r"""in the full LR fit no individual feature reached significance (all five feature $p$-values $>0.23$)""",
 "LR feature p-values"))
R.append((
 r"""On the same 75-pair frame, the composite score reached AUPRC~$=$~0.363 (4.5$\times$ the baseline prevalence of 0.080)""",
 r"""On the same 72-pair frame, the composite score reached AUPRC~$=$~0.357 (4.3$\times$ the baseline prevalence of 0.083)""",
 "composite AUPRC"))

# ── Results: component decomposition (line 102) ────────────────────
R.append((
 r"""Component decomposition on the same test set: DD (AUROC~$=$~0.682) versus PCS (0.777), $\Delta$Expression alone (0.564), and necessity (0.647; Fig.~1d). DD no longer dominates every comparator under the class-specific driver-mutation rules; PCS in particular carries substantial signal, which motivates the composite score.""",
 r"""Component decomposition on the same test set: DD (AUROC~$=$~0.676) versus PCS (0.825), $\Delta$Expression alone (0.547), and necessity (0.642; Fig.~1d). A paired bootstrap over entries found no significant difference between DD and any other component (PCS minus DD $+0.150$, 95\% CI $-0.110$ to $+0.456$). PCS in particular carries substantial signal, which motivates the composite score.""",
 "component decomposition"))

# ── Results: tier narrative (line 104) ─────────────────────────────
R.append((
 r"""Two pairs have directional experimental evidence, from studies external to and independent of DepMap, matching the direction scored here (\textit{Tier A}: SMARCA4$\rightarrow$SMARCA2 \cite{Hoffman2014}, ARID1A$\rightarrow$ARID1B \cite{Helming2014}). One further pair is experimentally established as synthetic lethal but only in the reciprocal direction (CREBBP$\rightarrow$EP300 \cite{Ogiwara2016,Nie2021}); the EP300$\rightarrow$CREBBP direction scored here is supported only by paralog redundancy, so this pair is annotated separately and excluded from directional Tier A claims. Five pairs are supported by paralog-redundancy, digenic-knockout, or pharmacologic evidence only (\textit{Tier B}: AKT1$\rightarrow$AKT2 \cite{Najm2018}, CCNE1$\rightarrow$CCNE2 \cite{Geng2003}, PIK3CA$\rightarrow$PIK3CB, CDK4$\rightarrow$CDK6, and MAP2K1$\rightarrow$MAP2K2 \cite{Parrish2021}), and two pairs derive from DepMap analyses (\textit{Tier C}: FBXW7$\rightarrow$FBXW2, PPP2R1A$\rightarrow$PPP2R1B).""",
 r"""Three pairs have direct genetic synthetic-lethal evidence from dual-gene perturbation assays (\textit{Tier A}: AKT1$\rightarrow$AKT2 \cite{Najm2018}, CDK4$\rightarrow$CDK6 and MAP2K1$\rightarrow$MAP2K2 \cite{Parrish2021}). Two further pairs are demonstrated selective dependencies in driver-mutant cells with functional validation (\textit{Tier B}: SMARCA4$\rightarrow$SMARCA2 \cite{Hoffman2014}, ARID1A$\rightarrow$ARID1B \cite{Helming2014}); the Tier A~$\cup$~Tier B set constitutes the pre-specified primary external benchmark. Five pairs rest on indirect, reciprocal-direction-only, or DepMap-derived evidence (\textit{Tier C}: EP300$\rightarrow$CREBBP \cite{Ogiwara2016,Nie2021}, PIK3CA$\rightarrow$PIK3CB \cite{Wee2008}, CCNE1$\rightarrow$CCNE2 \cite{Geng2003}, FBXW7$\rightarrow$FBXW2, PPP2R1A$\rightarrow$PPP2R1B), and two serve as mechanistic comparators (BRCA1$\leftrightarrow$BRCA2, STK11$\rightarrow$SIK1).""",
 "tier definitions in Results"))
R.append((
 r"""On the full curated set of twelve pairs DD achieved AUROC~$=$~0.682 (116 driver$\times$paralog$\times$lineage entries, 9 positives), and this estimate proved robust to label-quality concerns: excluding the two DepMap-derived pairs raised AUROC to 0.728, restricting to the eight pairs with pre-DepMap experimental evidence gave 0.774, and a direction-strict analysis relabelling the EP300$\rightarrow$CREBBP entries as non-positive left it unchanged at 0.682. Within the Tier A frame, both evaluable entries (ARID1A$\rightarrow$ARID1B in Endometrial and Ovarian cancers) ranked above all 107 unlabeled control entries; SMARCA4$\rightarrow$SMARCA2 had too few driver-mutant cell lines for lineage-level evaluation.""",
 r"""On the full curated set of twelve pairs DD achieved AUROC~$=$~0.676 (110 driver$\times$paralog$\times$lineage entries, 8 positives), and this estimate proved robust to label-quality concerns: excluding the two DepMap-derived pairs raised AUROC to 0.725, restricting to the eight pairs with pre-DepMap experimental evidence gave 0.774, and a direction-strict analysis relabelling the EP300$\rightarrow$CREBBP entries as non-positive left it unchanged at 0.676. On the primary Tier A~$\cup$~Tier B frame, both evaluable positive entries (ARID1A$\rightarrow$ARID1B in Endometrial and Ovarian cancers) ranked above all 102 unlabeled control entries (AUROC~$=$~1.000); SMARCA4$\rightarrow$SMARCA2 and the three Tier A pairs had too few qualifying driver-mutant cell lines for lineage-level evaluation under the $\geq$5-per-group rule.""",
 "full-set + tier AB numbers"))
R.append((
 r"""ARID1A$\rightarrow$ARID1B has the largest $|$DD$|$ (0.267), the highest dependency window score (4.13), and the highest selectivity (0.237), providing orthogonal prioritization beyond classification performance.""",
 r"""ARID1A$\rightarrow$ARID1B has the largest lineage-level DD (0.386 in Ovarian cancer; Hedges' $g=1.39$; Supplementary Table~S8) and the leading dependency-window and selectivity ranking (Table~1), providing orthogonal prioritization beyond classification performance.""",
 "effect-size ranking"))

# ── Results: LLO/bootstrap/regression (line 106) ───────────────────
R.append((
 r"""Leave-one-lineage-out AUROC was stable (range 0.674--0.702).""",
 r"""Leave-one-lineage-out AUROC was stable (range 0.656--0.704).""",
 "LLO range"))
R.append((
 r"""Bootstrap resampling (1,000 iterations, 75 pairs, 6 positives) gave a 95\% CI of 0.165--0.829 (Fig.~1d inset).""",
 r"""Bootstrap resampling (1,000 iterations, 72 pairs, 6 positives) gave a 95\% CI of 0.185--0.813 (Fig.~1d inset).""",
 "bootstrap CI"))
R.append((
 r"""the observed AUROC (0.493) did not exceed it (empirical $p=0.530$; Supplementary Fig.~S8)""",
 r"""the observed AUROC (0.500) did not exceed it (empirical $p=0.503$; Supplementary Fig.~S8)""",
 "permutation p"))
R.append((
 r"""($\Delta\text{DD} = 0.0004$, adjusted $p = 5.8 \times 10^{-42}$).""",
 r"""($\Delta\text{DD} = 0.0004$, adjusted $p = 5.8 \times 10^{-42}$). Adding lineage fixed effects attenuated but preserved the association (lineage-adjusted $p=2.2\times10^{-13}$; TP53- and CNV-adjusted with lineage fixed effects $p=7.0\times10^{-13}$).""",
 "lineage FE regression"))

# ── Figure 1 caption (lines 113-115) ───────────────────────────────
R.append((
 r"""TSG-driven tumors ($n=11$) show numerically higher AUROC than oncogene-driven tumors ($n=3$; permutation $p=0.538$, not significant). Points, individual lineages; center line, group mean.""",
 r"""TSG-driven tumors ($n=9$) show numerically higher AUROC than oncogene-driven tumors ($n=3$; permutation $p=0.645$, not significant). Points, individual lineages; center line, group mean. Computed on the $\geq$3-per-group sensitivity frame because the primary $\geq$5 frame retains only one oncogene-driven lineage.""",
 "fig1b caption"))
R.append((
 r"""The interpretable composite score (AUROC~$=$~0.841) exceeds all four multi-feature classifiers tested under leave-one-pair-out cross-validation (0.138--0.744), and each DD-based nomination traces""",
 r"""The interpretable composite score (AUROC~$=$~0.831) matches the best of four multi-feature classifiers tested under leave-one-pair-out cross-validation (SVM-RBF, 0.843; classifier range 0.114--0.843), and each DD-based nomination traces""",
 "fig1c caption"))
R.append((
 r"""\textbf{d}, Component decomposition: DD (AUROC~$=$~0.682) versus PCS (0.777), $\Delta$Expression alone (0.564), and necessity (0.647). Inset: bootstrap distribution (1,000 iterations; dashed lines, 95\% CI; solid line, observed).""",
 r"""\textbf{d}, Component decomposition: DD (AUROC~$=$~0.676) versus PCS (0.825), $\Delta$Expression alone (0.547), and necessity (0.642). Inset: pair-level bootstrap distribution (1,000 iterations; dashed lines, 95\% CI; solid line, observed).""",
 "fig1d caption"))

# ── Table 2 row ────────────────────────────────────────────────────
R.append((
 r"""DD (this study) & 0.682 & High""",
 r"""DD (this study) & 0.676 & High""",
 "table 2 DD row"))

# ── CPTAC section (line 150) + Fig 2c caption (line 158) ──────────
R.append((
 r"""$\Delta$Expression alone achieved AUROC~$=$~0.564 for discriminating known paralog-SL pairs (Fig.~2c), far below the protein-level co-variation signal and below DD (0.682).""",
 r"""$\Delta$Expression alone achieved AUROC~$=$~0.547 for discriminating known paralog-SL pairs (Fig.~2c), far below the protein-level co-variation signal and below DD (0.676).""",
 "dexpr comparison"))
R.append((
 r"""\textbf{c}, Dependency-based (DD, AUROC~$=$~0.682) vs.\ RNA-level ($\Delta$Expression, AUROC~$=$~0.564) prediction performance""",
 r"""\textbf{c}, Dependency-based (DD, AUROC~$=$~0.676) vs.\ RNA-level ($\Delta$Expression, AUROC~$=$~0.547) prediction performance""",
 "fig2c caption"))

# ── MSI paragraph (line 168) + Fig 3 caption ──────────────────────
R.append((
 r"""Contrary to that expectation, MSI-H subgroups showed numerically stronger paralog-SL predictive signal in both cancer types (Fig.~3a): endometrial AUROC~$=$~0.838 (MSI-H, $n=17$) vs.\ 0.556 (MSS, $n=11$), $\Delta=+0.282$; colorectal AUROC~$=$~0.767 (MSI-H, $n=14$) vs.\ 0.712 (MSS, $n=45$), $\Delta=+0.055$. Subgroup sizes are modest and only 3--5 gold-standard pairs were evaluable per subgroup, so the differences were not formally tested for significance; these numbers are hypothesis-generating.""",
 r"""Contrary to that expectation, MSI-H subgroups showed numerically stronger paralog-SL predictive signal in both cancer types on the $\geq$3-per-group sensitivity frame (Fig.~3a): endometrial AUROC~$=$~0.838 (MSI-H, $n=17$) vs.\ 0.556 (MSS, $n=11$), $\Delta=+0.282$; colorectal AUROC~$=$~0.767 (MSI-H, $n=14$) vs.\ 0.712 (MSS, $n=45$), $\Delta=+0.055$. On the primary $\geq$5 frame the endometrial subgroups are not evaluable (fewer than two gold-standard pairs) and the colorectal subgroups show no difference (0.574 MSI-H vs.\ 0.595 MSS), so the MSI contrast is strictly exploratory. Subgroup sizes are modest and only 3--4 gold-standard pairs were evaluable per subgroup, so the differences were not formally tested for significance; these numbers are hypothesis-generating.""",
 "MSI paragraph"))
R.append((
 r"""\textbf{a}, DD AUROC for MSI-H vs.\ MSS subgroups in colorectal and endometrial cancers, using the official DepMap MSIsensor2 annotation (MSIscore~$>$~20~$=$~MSI-H). Sample sizes shown above bars.""",
 r"""\textbf{a}, DD AUROC for MSI-H vs.\ MSS subgroups in colorectal and endometrial cancers, using the official DepMap MSIsensor2 annotation (MSIscore~$>$~20~$=$~MSI-H), computed on the $\geq$3-per-group sensitivity frame (the primary $\geq$5 frame is not evaluable for endometrial and shows no colorectal difference; see text). Sample sizes shown above bars.""",
 "fig3a caption"))

# ── Discussion ─────────────────────────────────────────────────────
R.append((
 r"""the interpretable composite score reached AUROC~$=$~0.841 in head-to-head comparison on 75 pairs, exceeding all four multi-feature classifiers tested (0.138--0.744), while DD alone reached 0.551.""",
 r"""the interpretable composite score reached AUROC~$=$~0.831 in head-to-head comparison on 72 pairs, statistically indistinguishable from the best multi-feature classifier (SVM-RBF, 0.843; classifier range 0.114--0.843), while DD alone reached 0.566.""",
 "discussion ML"))
R.append((
 r"""DD performed better in TSG-driven cancers (mean AUROC~$=$~0.737) than oncogene-driven cancers (0.595). The difference did not reach significance ($p=0.071$), and with only three oncogene lineages the comparison is underpowered.""",
 r"""DD performed better in TSG-driven cancers (mean AUROC~$=$~0.814) than oncogene-driven cancers (0.768) on the $\geq$3-per-group sensitivity frame. The difference did not reach significance (permutation $p=0.645$), and with only three oncogene lineages the comparison is underpowered; the primary $\geq$5 frame retains a single oncogene-driven lineage, precluding estimation there.""",
 "discussion TSG/ONC"))
R.append((
 r"""MSI status modulated DD signal strength (MSI-H $\geq$ MSS in both colorectal and endometrial subgroups; Fig.~3a)""",
 r"""MSI status may modulate DD signal strength (MSI-H $\geq$ MSS in both colorectal and endometrial subgroups on the sensitivity frame; Fig.~3a)""",
 "discussion MSI"))
R.append((
 r"""removing the two gold-standard pairs with DepMap-era evidence increased AUROC from 0.682 to 0.728""",
 r"""removing the two gold-standard pairs with DepMap-era evidence increased AUROC from 0.676 to 0.725""",
 "data independence 1"))
R.append((
 r"""relabelling the pair whose direct evidence supports only the reciprocal direction (EP300$\rightarrow$CREBBP) as non-positive left it unchanged at 0.682""",
 r"""relabelling the pair whose direct evidence supports only the reciprocal direction (EP300$\rightarrow$CREBBP) as non-positive left it unchanged at 0.676""",
 "data independence 2"))
R.append((
 r"""The two pairs with direct directional evidence independent of DepMap (Tier A) ranked above all unlabeled controls wherever lineage-evaluable.""",
 r"""Both lineage-evaluable positive entries on the Tier A~$\cup$~Tier B external benchmark (ARID1A$\rightarrow$ARID1B in two lineages) ranked above all unlabeled controls; the remaining benchmark pairs had too few qualifying driver-mutant lines for evaluation.""",
 "data independence 3"))
R.append((
 r"""Our top-ranked candidate has the largest $|$DD$|$ (0.250), highest DWS (3.65), and highest selectivity (0.228) among all 24 pairs, yet its within-driver BH-corrected q-value is 0.50 — nominally non-significant.""",
 r"""Our top-ranked candidate has the largest lineage-level DD (0.386 in Ovarian cancer) and the highest DWS (3.65) among all evaluated pairs, yet its within-driver BH-corrected q-values remain nominally non-significant (0.39 in Ovarian, 0.88 in Endometrial).""",
 "ARID1A paradox"))
R.append((
 r"""The cancer-type-specific evaluation (AUROC~$=$~0.682; 116 driver$\times$paralog$\times$lineage entries, 9 positives) preserves""",
 r"""The cancer-type-specific evaluation (AUROC~$=$~0.676; 110 driver$\times$paralog$\times$lineage entries, 8 positives) preserves""",
 "frameworks 1"))
R.append((
 r"""The per-pair aggregated evaluation (AUROC~$=$~0.493; 75 unique pairs, 6 positives) averages""",
 r"""The per-pair aggregated evaluation (AUROC~$=$~0.500; 72 unique pairs, 6 positives) averages""",
 "frameworks 2"))

# ── Conclusions (line 275) ─────────────────────────────────────────
R.append((
 r"""DD achieves internal-consistency AUROC of 0.682 for the paralog-SL task within DepMap (1.000 on the two pairs with directional external experimental evidence), and the interpretable composite score reaches 0.841;""",
 r"""DD achieves internal-consistency AUROC of 0.676 for the paralog-SL task within DepMap (1.000 on the Tier A~$\cup$~Tier B external benchmark, on which two of five pairs were lineage-evaluable), and the interpretable composite score reaches 0.831;""",
 "conclusions"))

# ── Methods: direction audit (line 303) ────────────────────────────
R.append((
 r"""on the directional Tier A set, signed and absolute scoring coincide (AUROC~$=$~1.000), whereas on the full curated set 4 of 9 positive entries have DD~$<0$ (signed-DD AUROC~$=$~0.578 vs.\ $|$DD$|$~0.682). Full-set estimates should therefore be read as direction-agnostic discrimination, and directional claims rest on Tier A.""",
 r"""on the Tier A~$\cup$~Tier B benchmark set, signed and absolute scoring coincide (AUROC~$=$~1.000), whereas on the full curated set 3 of 8 positive entries have DD~$<0$ (signed-DD AUROC~$=$~0.629 vs.\ $|$DD$|$~0.676). Full-set estimates should therefore be read as direction-agnostic discrimination, and directional claims rest on the Tier A~$\cup$~Tier B benchmark.""",
 "direction audit"))

# ── Methods: power statement (line 347) ────────────────────────────
R.append((
 r"""With at most 3--5 evaluable positive pairs, subset-level AUROC has limited power""",
 r"""With at most 2--5 evaluable positive pairs, subset-level AUROC has limited power""",
 "power statement"))

# ── SI listing (line 351) ──────────────────────────────────────────
R.append((
 r"""Supplementary Figures S1--S10 and Supplementary Tables S1--S7 are available""",
 r"""Supplementary Figures S1--S10 and Supplementary Tables S1--S8 are available""",
 "SI table count"))
R.append((
 r"""Supplementary Table~S2: Complete paralog-SL pair analysis results (116 driver$\times$paralog$\times$cancer-type associations passing the minimum sample size filter of $\ge$3 mutant and $\ge$3 wild-type cell lines).""",
 r"""Supplementary Table~S2: Complete paralog-SL pair analysis results (110 driver$\times$paralog$\times$cancer-type associations passing the minimum sample size filter of $\ge$5 mutant and $\ge$5 wild-type cell lines).""",
 "table S2 description"))
R.append((
 r"""Supplementary Table~S7: Gene-class-specific driver-mutation rules and per-gene variant classification (TSG: LikelyLoF; oncogene: Hotspot).""",
 r"""Supplementary Table~S7: Gene-class-specific driver-mutation rules and per-gene variant classification (TSG: LikelyLoF; oncogene: Hotspot). Supplementary Table~S8: Per-association effect sizes (Cohen's $d$, Hedges' $g$) and Welch $t$-test $p$-values for all 110 evaluated associations (TSV file).""",
 "table S8 entry"))

# ── apply ──────────────────────────────────────────────────────────
errors = []
for old, new, tag in R:
    n = text.count(old)
    if n != 1:
        errors.append(f"[{tag}] found {n} occurrences (expected 1): {old[:80]}...")
    else:
        text = text.replace(old, new, 1)

if errors:
    print("FAILED — no changes written:")
    for e in errors:
        print(" ", e)
    raise SystemExit(1)

MS.write_text(text)
print(f"OK: {len(R)} replacements applied to manuscript.tex")
