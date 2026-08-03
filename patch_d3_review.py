#!/usr/bin/env python3
"""
patch_d3_review.py — 2026-08-01 external-review P0 batch (cross-document items).

Every replacement is an exact string asserted to occur exactly once in the
target file; the script aborts before writing anything if any assertion
fails. No analysis numbers are fabricated here — numeric strings updated
below come from output artifacts (ml_benchmark.json fold-internal composite)
or from already-audited values; the audit suite is re-run afterwards.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
MS = ROOT / "manuscript.tex"
SU = ROOT / "supplementary.tex"
CL = ROOT / "cover_letter.md"
AU = ROOT / "audit_manuscript_numbers.py"
CK = ROOT / "submission_checklist.md"

EDITS = [
    # ── manuscript.tex ─────────────────────────────────────────────
    (MS,
     "0.525 for the eight pairs with pre-DepMap experimental evidence; 1.000 on the Tier A~$\\cup$~B benchmark, both evaluable entries being the same pair).",
     "0.525 for the eight pairs with pre-DepMap experimental evidence; the sole lineage-evaluable Tier A/B pair, ARID1A$\\rightarrow$ARID1B, ranked above all unlabeled controls in both evaluable lineages)."),
    (MS,
     "Fewer than 30,000 human SL interactions are catalogued in SynLethDB \\cite{Guo2016,SynLethDB2022}, and the existing",
     "Available SL catalogues contain tens of thousands of heterogeneous, unevenly validated interactions \\cite{Guo2016,SynLethDB2022}, and the existing"),
    (MS,
     "benchmark evaluation on the curated gold-standard set (Results, benchmark comparison)",
     "benchmark evaluation on the curated evaluation set (Results, benchmark comparison)"),
    (MS,
     "the paralog compensation score (PCS $=$ DD $\\times$ necessity, where necessity is the paralog's mean gene-effect score)",
     "the paralog compensation score (PCS $=$ $\\max(\\Delta\\mathrm{Expression},0)\\times\\max(\\mathrm{necessity},0)$, where necessity $=-\\overline{\\mathrm{Chronos}}$ so that larger values mark more essential paralogs)"),
    (MS,
     "The composite score (AUROC~$=$~0.831) and SVM-RBF (0.841) performed best and are statistically indistinguishable at this sample size, followed by",
     "The composite score (AUROC~$=$~0.831) and SVM-RBF (0.841) performed best, with no detectable difference between them at this sample size (the paired-bootstrap confidence interval is too wide to establish equivalence), followed by"),
    (MS,
     "because the same benchmark informed its development, its AUROC is descriptive context rather than independent validation.",
     "because the same benchmark informed its development, its AUROC is descriptive context rather than independent validation. Min--max scaling of the composite components is computed within each cross-validation fold on training pairs only; this fold-internal scaling leaves the composite AUROC unchanged to three decimals (0.831) and changes AUPRC from 0.356 to 0.368. The composite and the dependency-window score use magnitude-based $|$DD$|$ as an effect-size term, whereas signed DD remains the directional metric for biological interpretation."),
    (MS,
     "fewer than two gold-standard pairs",
     "fewer than two curated positive pairs"),
    (MS,
     "only 3--4 gold-standard pairs were evaluable per subgroup",
     "only 3--4 curated positive pairs were evaluable per subgroup"),
    (MS,
     "reached AUROC~$=$~0.831 in head-to-head comparison, statistically indistinguishable from the best multi-feature classifier (SVM-RBF, 0.841; classifier range 0.240--0.841)",
     "reached AUROC~$=$~0.831 in head-to-head comparison, with no detectable difference from the best multi-feature classifier (SVM-RBF, 0.841; classifier range 0.240--0.841) and a confidence interval too wide to establish equivalence"),
    (MS,
     "Both DD and the gold-standard labels come from DepMap CRISPR data",
     "Both DD and the curated labels come from DepMap CRISPR data"),
    (MS,
     "so we nominate KMT2D$\\rightarrow$KMT2C for combinatorial-CRISPR testing rather than claiming a validated dependency.",
     "so we nominate KMT2D$\\rightarrow$KMT2C for combinatorial-CRISPR testing rather than claiming a validated dependency. Notably, this pair was included in the unbiased digenic screen of Flister et al.\\ \\cite{Flister2025} and was not called lethal in either of the two models screened (NCI-H1299, MDA-MB-231) --- consistent with the low genotype penetrance of digenic interactions and with the modest dependency shift we observe, but an external caveat we report explicitly."),
    (MS,
     "\\item \\textit{Gold-standard heterogeneity:}",
     "\\item \\textit{Evaluation-set heterogeneity:}"),
    (MS,
     "DD and the gold standard both derive from DepMap data",
     "DD and the curated labels both derive from DepMap data"),
    (MS,
     "(1.000 on the Tier A~$\\cup$~Tier B external benchmark, where both lineage-evaluable entries were the same pair, ARID1A$\\rightarrow$ARID1B)",
     "(the sole lineage-evaluable Tier A/B pair, ARID1A$\\rightarrow$ARID1B, ranked above all unlabeled controls in both evaluable lineages)"),
    (MS,
     "\\textbf{Gold-standard positives and matched negatives.}",
     "\\textbf{Curated pairs and matched unlabeled controls.}"),
    (MS,
     "using the 12 gold-standard paralog-SL pairs as the positive set",
     "using the 12 curated paralog-SL pairs as the positive set"),
    (MS,
     "(two genes, CCNE1 and MAPK1, contribute zero mutant lines under the hotspot rule)",
     "(CCNE1 contributes zero mutant lines under the hotspot rule; MAPK1, audited for the external in4mer analysis rather than as a panel gene, likewise yields zero)"),
    (MS,
     "(seed 42; Supplementary Tables~S5 and~S11; \\texttt{dws\\_robustness.py})",
     "(seed 42; Supplementary Tables~S5 and~S10; \\texttt{dws\\_robustness.py})"),
    (MS,
     "Gold-standard pairs were stratified into evidence tiers",
     "Curated pairs were stratified into evidence tiers"),
    (MS,
     "archived at Zenodo (DOI: to be assigned upon acceptance) and the GitHub repository below.",
     "archived at Zenodo (concept DOI \\url{https://doi.org/10.5281/zenodo.21502030}, which resolves to the latest release; the submission snapshot is release v1.3.2 with version DOI \\url{https://doi.org/10.5281/zenodo.21634700}) and the GitHub repository below."),
    (MS,
     "\\texttt{compute\\_auroc()} against gold standard",
     "\\texttt{compute\\_auroc()} against the curated positive set"),
    (MS,
     "A frozen release with all processed data tables is archived at Zenodo (\\url{https://doi.org/10.5281/zenodo.21502030}).",
     "A frozen release with all processed data tables is archived at Zenodo (concept DOI \\url{https://doi.org/10.5281/zenodo.21502030}; submission snapshot: GitHub tag \\texttt{v1.3.2}, commit \\texttt{66fc633}, version DOI \\url{https://doi.org/10.5281/zenodo.21634700}). The R package is archived separately (concept DOI \\url{https://doi.org/10.5281/zenodo.21502113}; v1.1.2 version DOI \\url{https://doi.org/10.5281/zenodo.21634707}). Each new release mints a new version DOI under the same concept DOI, and the accepted-version snapshot will be deposited identically."),
    (MS,
     "composite-weight and in4mer seed sensitivities, the external digenic-screen hold-out",
     "composite fold-internal scaling, composite-weight and in4mer seed sensitivities, the external digenic-screen hold-out"),
    (MS,
     "Supplementary Table~S3: Gold-standard paralog-SL pair evidence.",
     "Supplementary Table~S3: Curated paralog-SL pair evidence with tier assignment."),
    (MS,
     "Machine-readable mirror of the 12 evidence-tiered gold-standard pairs.",
     "Machine-readable mirror of the 12 evidence-tiered curated pairs (file name retained for continuity)."),
    # ── supplementary.tex ──────────────────────────────────────────
    (SU,
     "\\subsection{Gold-standard curation and tier assignment}",
     "\\subsection{Curated-pair curation and tier assignment}"),
    (SU,
     "all with default scikit-learn 1.3.0 hyperparameters and no tuning",
     "all with default scikit-learn hyperparameters and no tuning (production versions pinned in the repository Dockerfile; the audit suite is additionally verified under scikit-learn 1.9.0)"),
    (SU,
     "red bars show known gold-standard positives",
     "red bars show known positives"),
    (SU,
     "for the known gold-standard SL pairs",
     "for the curated SL pairs"),
    (SU,
     "Known+: number of gold-standard positive pairs evaluable in each lineage.",
     "Known+: number of curated positive pairs evaluable in each lineage."),
    (SU,
     "\\item \\texttt{novelty}: ``Known'' if in gold-standard set (Table~S3), else ``Novel''.",
     "\\item \\texttt{novelty}: ``Known'' if in the curated set (Table~S3), else ``Novel''."),
    (SU,
     "\\item \\texttt{is\\_known\\_paralog\\_sl}: Boolean; TRUE if pair is in gold-standard set.",
     "\\item \\texttt{is\\_known\\_paralog\\_sl}: Boolean; TRUE if pair is in the curated set (Table~S3)."),
    (SU,
     "\\subsection{Table S3: Gold-standard paralog-SL pairs with evidence provenance and tier assignment}",
     "\\subsection{Table S3: Curated paralog-SL pairs with evidence provenance and tier assignment}"),
    (SU,
     "(CCNE1 and MAPK1 contribute zero mutant lines)",
     "(CCNE1 contributes zero mutant lines)"),
    (SU,
     "Only one of the 13 pairs (MAPK1--MAPK3) involves a gene in our 40-gene driver panel, and MAPK1 carries qualifying oncogenic hotspot mutations in zero DepMap 26Q1 cell lines (Table~S6), so no in4mer gold-standard pair is evaluable in our driver-mutation-conditioned framework.",
     "No in4mer pair involves a 40-panel driver gene; as the closest case we additionally audited MAPK1 (of MAPK1--MAPK3), which carries qualifying oncogenic hotspot mutations in zero DepMap 26Q1 cell lines (Table~S6), so no in4mer gold-standard pair is evaluable in our driver-mutation-conditioned framework."),
    (SU,
     "CCNE1 and MAPK1 are amplification-driven oncogenes with no qualifying hotspot mutations.}",
     "CCNE1 and MAPK1 are amplification-driven oncogenes with no qualifying hotspot mutations. The table covers the 40-gene driver panel plus genes audited outside the panel (CDK4 and MAP2K1 for the curated benchmark; MAPK1 for the external in4mer analysis).}"),
    (SU,
     "Known drug--target biology is recovered: MEK inhibitors with KRAS-mutant lines, mTOR/AKT inhibitors with PTEN-mutant lines, HDAC inhibitors with EP300-mutant lines. Drugs shown are not necessarily direct binders of the paralog protein.",
     "The displayed top associations are dominated by cytotoxic chemotherapies (tubulin inhibitors, topoisomerase poisons, and mitotic-kinase inhibitors) with $\\Delta$AUC$<0$, indicating genotype-conditioned cytotoxic sensitivity rather than paralog-specific targeting; the targeted-agent assay-validity anchors (MEK inhibitors--KRAS, mTOR/AKT inhibitors--PTEN, HDAC inhibitors--EP300) are shown in main-text Fig.~4a. Drugs shown are not necessarily direct binders of the paralog protein."),
    # ── cover_letter.md ────────────────────────────────────────────
    (CL,
     "is statistically indistinguishable from the best of four standard classifiers under leave-one-pair-out cross-validation (SVM-RBF 0.841)",
     "does not differ detectably from the best of four standard classifiers under leave-one-pair-out cross-validation (SVM-RBF 0.841; the confidence interval is too wide to establish equivalence)"),
    (CL,
     "paralogSL` v1.1.1",
     "paralogSL` v1.1.2"),
    (CL,
     "which always resolve to the latest release)",
     "which always resolve to the latest release; version-specific DOIs, release tag, and commit are listed in the manuscript)"),
    (CL,
     "a single entry point (`verify_all.sh`) runs five pipeline stages that recompute 390 numeric claims directly from the analysis artifacts plus a 31-test suite, completing in under a minute on an Apple M4 Max (128 GB RAM)",
     "a single entry point (`verify_all.sh`) runs four audit modules plus a 31-test suite, recomputing 392 numeric claims directly from the analysis artifacts in ~30 s on a standard workstation"),
    # ── audit_manuscript_numbers.py: fold-internal composite claims ─
    (AU,
     'check("ml_comp", "Composite alone AUROC", "0.831", ml["single_feature"]["composite_alone"], ML)',
     'check("ml_comp", "Composite alone AUROC", "0.831", ml["single_feature"]["composite_alone"], ML)\n'
     'check("ml_comp_lofo", "Composite fold-internal-scaling AUROC", "0.831",\n'
     '      ml["single_feature"]["composite_alone_lofo"], ML)\n'
     'check("ml_comp_lofo_auprc", "Composite fold-internal-scaling AUPRC", "0.368",\n'
     '      ml["composite_auprc_lofo"], ML)'),
    # ── manuscript claim-count self-reference ──────────────────────
    (MS,
     "390 claims in total",
     "392 claims in total"),
]

REPLACED = []
for path, old, new in EDITS:
    text = path.read_text(encoding="utf-8")
    n = text.count(old)
    if n != 1:
        print(f"ABORT: {path.name}: expected 1 occurrence, found {n}:\n  {old[:110]}")
        sys.exit(1)
    REPLACED.append((path, old, new))

for path, old, new in REPLACED:
    path.write_text(path.read_text(encoding="utf-8").replace(old, new),
                    encoding="utf-8")

# checklist: 390 -> 392 (three audited spots)
ck = CK.read_text(encoding="utf-8")
n = ck.count("390")
if n != 3:
    print(f"ABORT: submission_checklist.md: expected 3 occurrences of '390', found {n}")
    sys.exit(1)
CK.write_text(ck.replace("390", "392"), encoding="utf-8")

print(f"OK: {len(REPLACED)} exact replacements across {len(set(p for p,_,_ in REPLACED))} files + checklist 390->392 (3 spots)")
