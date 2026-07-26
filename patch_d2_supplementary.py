#!/usr/bin/env python3
"""D2 supplementary sync: Table S1/S4 rebuilt on min5 primary frame, captions updated,
Table S2 dimensions, sensitivity-frame notes. Every replacement must match exactly once."""
from pathlib import Path

SP = Path(__file__).parent / "supplementary.tex"
text = SP.read_text()

R = []

# ── Fig S2 caption (line 58) ───────────────────────────────────────
R.append((
 r"""\caption{\textbf{Cross-cancer DD AUROC.} DD AUROC for each of the 14 evaluable solid tumor types, grouped by signal strength. a, Strong signal (AUROC~$\ge$~0.7; 10 lineages including Biliary Tract 0.990, Esophagogastric 0.969, Pancreatic 0.949). b, Moderate signal (AUROC 0.5--0.7; Endometrial 0.690, Ovarian 0.651). c, Weak signal (AUROC~$<$~0.5; Melanoma 0.481, Cervical 0.476). Dashed line indicates random classifier (0.5). Bar labels show the number of tested paralog pairs per lineage. Cancer types with $<$2 known positive pairs were excluded from AUROC evaluation (see Table~S1 for the 9 non-evaluable types).}""",
 r"""\caption{\textbf{Cross-cancer DD AUROC.} DD AUROC for each of the 8 evaluable solid tumor types on the primary $\geq$5-per-group frame, grouped by signal strength. a, Strong signal (AUROC~$\ge$~0.7; 7 lineages: Esophagogastric 0.965, SCLC 0.906, Bladder Urothelial 0.844, Colorectal 0.828, Endometrial 0.818, Breast 0.750, NSCLC 0.741). b, Moderate signal (AUROC 0.5--0.7; Ovarian 0.661). Dashed line indicates random classifier (0.5). Bar labels show the number of tested paralog pairs per lineage. Cancer types with $<$2 known positive pairs were excluded from AUROC evaluation (see Table~S1 for the 15 non-evaluable types). On the $\geq$3-per-group sensitivity frame, four additional lineages become evaluable (Biliary Tract 0.990, Pancreatic 0.949, Melanoma 0.617, Cervical 0.500; 9 of 12 evaluable lineages exceed 0.7).}""",
 "figS2 caption"))

# ── Fig S6 caption ─────────────────────────────────────────────────
R.append((
 r"""overall mean near zero ($-0.005$)""",
 r"""overall mean near zero ($-0.016$)""",
 "figS6 mean diff"))
R.append((
 r"""the cross-cancer frame value (0.852; Fig.~1a)""",
 r"""the cross-cancer frame value (0.750; Fig.~1a)""",
 "figS6 cross-cancer ref"))

# ── Fig S8 caption (line 100) ──────────────────────────────────────
R.append((
 r"""\caption{\textbf{Bootstrap and negative control (per-pair evaluation, 75 pairs, 6 known positives).} a, Bootstrap distribution (1,000 iterations; 95\% CI: 0.165--0.829). b, Null distribution from 10,000 label-shuffled permutations (empirical $p = 0.530$). The observed per-pair AUROC (0.493) does not exceed the null mean (0.501): with only 6 positive pairs the aggregated per-pair evaluation has no detectable signal, and the lineage-level evaluation (0.682, reported in the benchmark comparison) serves as the primary framework because it evaluates each driver$\times$paralog$\times$lineage combination separately rather than aggregating across lineages.}""",
 r"""\caption{\textbf{Bootstrap and negative control (per-pair evaluation, 72 pairs, 6 known positives).} a, Bootstrap distribution (1,000 iterations; 95\% CI: 0.185--0.813). b, Null distribution from 10,000 label-shuffled permutations (empirical $p = 0.503$). The observed per-pair AUROC (0.500) does not exceed the null mean (0.501): with only 6 positive pairs the aggregated per-pair evaluation has no detectable signal, and the lineage-level evaluation (0.676, reported in the benchmark comparison) serves as the primary framework because it evaluates each driver$\times$paralog$\times$lineage combination separately rather than aggregating across lineages.}""",
 "figS8 caption"))

# ── Table S1: full row block + caption ─────────────────────────────
old_s1_rows = (
 r"""Biliary Tract & Biliary tract cancer & 39 & 50 & 2 & 0.990 \\
Esophagogastric & \shortstack[l]{Esophagogastric adenocarcinoma /\\esophageal SCC} & 62 & 69 & 4 & 0.969 \\
Pancreatic & Pancreatic adenocarcinoma & 44 & 41 & 2 & 0.949 \\
Glioma & Diffuse glioma & 66 & 28 & 2 & 0.885 \\
Bladder Urothelial & Bladder urothelial carcinoma & 30 & 20 & 3 & 0.882 \\
Breast & Invasive breast carcinoma & 50 & 48 & 4 & 0.852 \\
HNSCC & Head and neck squamous cell carcinoma & 60 & 18 & 2 & 0.812 \\
SCLC & Small cell lung cancer & 121 & 107 & 6 & 0.779 \\
NSCLC & Non-small cell lung cancer & 95 & 104 & 6 & 0.738 \\
Colorectal & Colorectal adenocarcinoma & 59 & 76 & 8 & 0.726 \\
\midrule
Endometrial & Endometrial carcinoma & 31 & 77 & 6 & 0.690 \\
Ovarian & Ovarian epithelial carcinoma & 55 & 65 & 8 & 0.651 \\
Melanoma & Melanoma & 78 & 31 & 4 & 0.481 \\
Cervical & Cervical carcinoma & 16 & 10 & 3 & 0.476 \\
\midrule
Mesothelioma & Pleural mesothelioma & 21 & 5 & 1 & --- \\
Hepatocellular & Hepatocellular carcinoma & 24 & 7 & 1 & --- \\
Renal Cell & Renal cell carcinoma & 26 & 9 & 0 & --- \\
Neuroblastoma & Neuroblastoma & 33 & 6 & 1 & --- \\
Osteosarcoma & Osteosarcoma & 19 & 12 & 0 & --- \\
Ewing Sarcoma & Ewing sarcoma & 18 & 2 & 0 & --- \\
Rhabdomyosarcoma & Rhabdomyosarcoma & 12 & 2 & 0 & --- \\
Other Sarcoma & Other soft-tissue sarcoma & 33 & 6 & 1 & --- \\
Thyroid & Thyroid carcinoma & 10 & 2 & 1 & --- \\""")
new_s1_rows = (
 r"""Esophagogastric & \shortstack[l]{Esophagogastric adenocarcinoma /\\esophageal SCC} & 62 & 50 & 3 & 0.965 \\
SCLC & Small cell lung cancer & 121 & 99 & 5 & 0.906 \\
Bladder Urothelial & Bladder urothelial carcinoma & 30 & 18 & 2 & 0.844 \\
Colorectal & Colorectal adenocarcinoma & 59 & 62 & 5 & 0.828 \\
Endometrial & Endometrial carcinoma & 31 & 62 & 5 & 0.818 \\
Breast & Invasive breast carcinoma & 50 & 12 & 2 & 0.750 \\
NSCLC & Non-small cell lung cancer & 95 & 96 & 5 & 0.741 \\
\midrule
Ovarian & Ovarian epithelial carcinoma & 55 & 50 & 7 & 0.661 \\
\midrule
Biliary Tract & Biliary tract cancer & 39 & 38 & 0 & --- \\
Melanoma & Melanoma & 78 & 22 & 1 & --- \\
Glioma & Diffuse glioma & 66 & 19 & 0 & --- \\
HNSCC & Head and neck squamous cell carcinoma & 60 & 13 & 1 & --- \\
Pancreatic & Pancreatic adenocarcinoma & 44 & 9 & 1 & --- \\
Hepatocellular & Hepatocellular carcinoma & 24 & 4 & 0 & --- \\
Renal Cell & Renal cell carcinoma & 26 & 4 & 0 & --- \\
Other Sarcoma & Other soft-tissue sarcoma & 33 & 4 & 0 & --- \\
Cervical & Cervical carcinoma & 16 & 3 & 1 & --- \\
Mesothelioma & Pleural mesothelioma & 21 & 2 & 0 & --- \\
Neuroblastoma & Neuroblastoma & 33 & 2 & 0 & --- \\
Osteosarcoma & Osteosarcoma & 19 & 0 & 0 & --- \\
Ewing Sarcoma & Ewing sarcoma & 18 & 0 & 0 & --- \\
Rhabdomyosarcoma & Rhabdomyosarcoma & 12 & 0 & 0 & --- \\
Thyroid & Thyroid carcinoma & 10 & 0 & 0 & --- \\""")
R.append((old_s1_rows, new_s1_rows, "table S1 rows"))

R.append((
 r"""\caption{Cell line counts and DD AUROC across all 23 analyzed solid tumor types. Top: AUROC~$>$~0.7 (10 types). Middle: AUROC~$\le$~0.7 (4 types). Bottom: insufficient known positives for AUROC (9 types). Known+: number of gold-standard positive pairs evaluable in each lineage. ``---'': AUROC not computable ($<$2 known positives). Full name: standardized cancer-type name corresponding to each short label; SCC, squamous cell carcinoma.}""",
 r"""\caption{Cell line counts and DD AUROC across all 23 analyzed solid tumor types (primary $\geq$5-per-group frame). Top: AUROC~$>$~0.7 (7 types). Middle: AUROC~$\le$~0.7 (1 type). Bottom: insufficient known positives for AUROC (15 types). Known+: number of gold-standard positive pairs evaluable in each lineage. ``---'': AUROC not computable ($<$2 known positives). Pairs: driver$\times$paralog combinations passing the $\geq$5 mutant / $\geq$5 wild-type filter. On the relaxed $\geq$3-per-group sensitivity frame, four additional lineages become evaluable (Biliary Tract 0.990, Pancreatic 0.949, Melanoma 0.617, Cervical 0.500). Full name: standardized cancer-type name corresponding to each short label; SCC, squamous cell carcinoma.}""",
 "table S1 caption"))

# ── Table S2 description + S8 entry ────────────────────────────────
R.append((
 r"""\texttt{TableS2\_FullResults.tsv} (116 rows $\times$ 13 columns, representing all driver$\times$paralog$\times$cancer-type combinations passing the minimum sample-size filter).""",
 r"""\texttt{TableS2\_FullResults.tsv} (110 rows $\times$ 17 columns, representing all driver$\times$paralog$\times$cancer-type combinations passing the minimum sample-size filter of $\geq$5 mutant and $\geq$5 wild-type cell lines). Effect sizes (Cohen's $d$, Hedges' $g$) and per-stratum Welch $t$-test $p$-values for the same 110 associations are provided as \texttt{TableS8\_EffectSizes.tsv} (Supplementary Table~S8).""",
 "table S2 description"))

# ── Table S4 rows (3-column, aligned spacing) ──────────────────────
old_s4 = (
 "Biliary Tract        & 39 & 50 & 0.990 \\\\\n"
 "Esophagogastric      & 62 & 69 & 0.969 \\\\\n"
 "Pancreatic           & 44 & 41 & 0.949 \\\\\n"
 "Glioma               & 66 & 28 & 0.885 \\\\\n"
 "Bladder Urothelial   & 30 & 20 & 0.882 \\\\\n"
 "Breast               & 50 & 48 & 0.852 \\\\\n"
 "HNSCC                & 60 & 18 & 0.812 \\\\\n"
 "SCLC                 & 121 & 107 & 0.779 \\\\\n"
 "NSCLC                & 95 & 104 & 0.738 \\\\\n"
 "Colorectal           & 59 & 76 & 0.726 \\\\\n"
 "\\midrule\n"
 "Endometrial          & 31 & 77 & 0.690 \\\\\n"
 "Ovarian              & 55 & 65 & 0.651 \\\\\n"
 "Melanoma             & 78 & 31 & 0.481 \\\\\n"
 "Cervical             & 16 & 10 & 0.476 \\\\")
new_s4 = (
 "Esophagogastric      & 62 & 50 & 0.965 \\\\\n"
 "SCLC                 & 121 & 99 & 0.906 \\\\\n"
 "Bladder Urothelial   & 30 & 18 & 0.844 \\\\\n"
 "Colorectal           & 59 & 62 & 0.828 \\\\\n"
 "Endometrial          & 31 & 62 & 0.818 \\\\\n"
 "Breast               & 50 & 12 & 0.750 \\\\\n"
 "NSCLC                & 95 & 96 & 0.741 \\\\\n"
 "\\midrule\n"
 "Ovarian              & 55 & 50 & 0.661 \\\\")
R.append((old_s4, new_s4, "table S4 rows"))

# ── min-sample mention (line ~262) ─────────────────────────────────
R.append((
 r"""with $\ge$3 mutant cell lines in at least one DepMap lineage.""",
 r"""with $\ge$5 mutant cell lines in at least one DepMap lineage.""",
 "min-sample threshold"))

# ── apply ──────────────────────────────────────────────────────────
errors = []
for old, new, tag in R:
    n = text.count(old)
    if n != 1:
        errors.append(f"[{tag}] found {n} occurrences (expected 1): {old[:80]!r}...")
    else:
        text = text.replace(old, new, 1)

if errors:
    print("FAILED — no changes written:")
    for e in errors:
        print(" ", e)
    raise SystemExit(1)

SP.write_text(text)
print(f"OK: {len(R)} replacements applied to supplementary.tex")
