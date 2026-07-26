#!/usr/bin/env python3
"""DWS unification: supplementary.tex Fig S7 caption + Table S6 rebuild (21 pairs, new frame)."""
from pathlib import Path

SP = Path(__file__).parent / "supplementary.tex"
text = SP.read_text()

R = []

# ── Fig S7 caption ─────────────────────────────────────────────────
R.append((
 r"""\caption{\textbf{Therapeutic window analysis.} a, All 24 paralog pairs ranked by mean dependency window score (DWS) across 5 cancer contexts (Ovarian, Endometrial, Breast, Colorectal, PanCancer). ARID1A$\rightarrow$ARID1B ranks first (DWS~$=$~3.65). b, In vitro selectivity tier classification: HIGH\_SELECTIVITY (selectivity~$>$~0.15 and DWS~$>$~1.0), MODERATE, LOW\_SELECTIVITY, and PAN\_ESSENTIAL ($f_{\text{pan-essential}}>0.5$). c, DWS values broken down by individual cancer context showing cross-lineage consistency for top candidates. d, Selectivity vs.\ DWS scatter plot; ARID1A$\rightarrow$ARID1B is the only pair in the upper-right quadrant (high DWS and high selectivity).}""",
 r"""\caption{\textbf{Therapeutic window analysis.} a, All 21 paralog pairs ranked by mean dependency window score (DWS) across up to 5 cancer contexts (Ovarian, Endometrial, Breast, Colorectal, PanCancer) on the $\geq$5-mutant frame. SMARCA4$\rightarrow$SMARCA2 (4.87) and ARID1A$\rightarrow$ARID1B (2.82) are the two HIGH\_SELECTIVITY pairs; NF1$\rightarrow$RASA2 attains a higher raw DWS through a near-zero pan-essential denominator with selectivity~$\approx$~0. b, In vitro selectivity tier classification: HIGH\_SELECTIVITY (selectivity~$>$~0.15 and DWS~$>$~1.0; $n=2$), MODERATE ($n=5$), LOW\_SELECTIVITY ($n=11$), and PAN\_ESSENTIAL ($f_{\text{pan-essential}}>0.5$; $n=3$). c, DWS values broken down by individual cancer context showing cross-lineage consistency for top candidates. d, Selectivity vs.\ DWS scatter plot; ARID1A$\rightarrow$ARID1B and SMARCA4$\rightarrow$SMARCA2 occupy the upper-right quadrant (high DWS and high selectivity).}""",
 "figS7 caption"))

# ── Table S6 section title ─────────────────────────────────────────
R.append((
 r"""\subsection{Table S6: All 24 paralog pairs ranked by preclinical prioritization score}""",
 r"""\subsection{Table S6: All 21 paralog pairs ranked by dependency window score}""",
 "table S6 title"))

# ── Table S6 rows ──────────────────────────────────────────────────
old_rows = r"""ARID1A & ARID1B & 0.250 & 3.647 & 0.228 & HS & Seq \\
NF1 & RASA2 & 0.062 & 3.260 & 0.002 & MOD & Seq \\
KMT2D & KMT2C & 0.112 & 2.558 & 0.131 & MOD & Seq \\
ATR & ATM & 0.049 & 2.182 & 0.031 & MOD & Seq \\
RB1 & RBL1 & 0.112 & 2.059 & 0.003 & MOD & Seq \\
AKT1 & AKT2 & 0.119 & 1.657 & 0.036 & MOD & Seq \\
PPP2R1A & PPP2R1B & 0.068 & 1.454 & 0.000 & LS & Seq \\
SMARCA4 & SMARCA2 & 0.057 & 1.268 & 0.039 & MOD & Seq \\
PIK3CA & PIK3CB & 0.137 & 1.245 & $-$0.052 & LS & Seq \\
EP300 & CREBBP & 0.187 & 1.185 & 0.120 & MOD & Seq \\
FBXW7 & FBXW2 & 0.030 & 1.032 & 0.000 & LS & Seq \\
STK11 & SIK1 & 0.037 & 1.018 & $-$0.010 & LS & Part \\
BRAF & RAF1 & 0.105 & 0.824 & $-$0.038 & LS & Seq \\
CDH1 & CDH2 & 0.057 & 0.722 & $-$0.030 & LS & Seq \\
TP53 & TP63 & 0.031 & 0.690 & 0.036 & MOD & Seq \\
MAP2K1 & MAP2K2 & 0.067 & 0.401 & $-$0.017 & LS & Seq \\
KRAS & NRAS & 0.071 & 0.358 & $-$0.055 & LS & Seq \\
BRCA1 & BRCA2 & 0.125 & 0.290 & 0.046 & PE & Func \\
BRCA2 & BRCA1 & 0.094 & 0.177 & $-$0.039 & PE & Func \\
PTEN & TNS2 & 0.041 & 0.172 & $-$0.012 & LS & Seq \\
CDK4 & CDK6 & 0.079 & 0.157 & $-$0.121 & LS & Seq \\
KRAS & HRAS & 0.028 & 0.139 & $-$0.045 & LS & Seq \\
PIK3R1 & CRKL & 0.100 & 0.111 & $-$0.112 & PE & Seq \\
CCNE1 & CCNE2 & 0.004 & 0.060 & $-$0.022 & LS & Seq \\"""
new_rows = r"""NF1 & RASA2 & 0.101 & 5.421 & 0.005 & MOD & Seq \\
SMARCA4 & SMARCA2 & 0.186 & 4.873 & 0.181 & HS & Seq \\
RB1 & RBL1 & 0.164 & 3.439 & 0.007 & MOD & Seq \\
ARID1A & ARID1B & 0.270 & 2.819 & 0.277 & HS & Seq \\
EP300 & CREBBP & 0.279 & 1.815 & 0.108 & MOD & Seq \\
KMT2D & KMT2C & 0.064 & 1.771 & 0.056 & MOD & Seq \\
ATR & ATM & 0.070 & 1.656 & $-$0.027 & LS & Seq \\
PIK3CA & PIK3CB & 0.152 & 1.356 & $-$0.065 & LS & Seq \\
BRAF & RAF1 & 0.245 & 1.174 & $-$0.178 & LS & Seq \\
PPP2R1A & PPP2R1B & 0.074 & 0.925 & 0.000 & LS & Seq \\
TP53 & TP63 & 0.035 & 0.813 & 0.037 & MOD & Seq \\
FBXW7 & FBXW2 & 0.025 & 0.720 & 0.000 & LS & Seq \\
CDH1 & CDH2 & 0.098 & 0.687 & $-$0.086 & LS & Seq \\
KRAS & NRAS & 0.083 & 0.420 & $-$0.060 & LS & Seq \\
STK11 & SIK1 & 0.017 & 0.392 & $-$0.011 & LS & Part \\
KRAS & HRAS & 0.043 & 0.206 & $-$0.040 & LS & Seq \\
BRCA2 & BRCA1 & 0.111 & 0.193 & $-$0.047 & PE & Func \\
PTEN & TNS2 & 0.042 & 0.174 & $-$0.006 & LS & Seq \\
PIK3R1 & CRKL & 0.152 & 0.167 & $-$0.161 & PE & Seq \\
BRCA1 & BRCA2 & 0.079 & 0.144 & $-$0.080 & PE & Func \\
MAP2K1 & MAP2K2 & 0.001 & 0.003 & $-$0.034 & LS & Seq \\"""
R.append((old_rows, new_rows, "table S6 rows"))

# ── Table S6 caption ───────────────────────────────────────────────
R.append((
 r"""\caption{All 24 analyzed paralog pairs ranked by mean dependency window score (DWS). $|$DD$|$: mean absolute Delta Dependency across cancer contexts. Select.: mean selectivity (fraction mutant-essential minus fraction wild-type-essential). Class: HS=HIGH\_SELECTIVITY, MOD=MODERATE, LS=LOW\_SELECTIVITY, PE=PAN\_ESSENTIAL. Type: Seq=sequence paralog, Part=partial homology, Func=functional analog.}""",
 r"""\caption{All 21 analyzed paralog pairs ranked by mean dependency window score (DWS) on the $\geq$5-mutant frame. $|$DD$|$: mean absolute Delta Dependency across cancer contexts. Select.: mean selectivity (fraction mutant-essential minus fraction wild-type-essential). Class: HS=HIGH\_SELECTIVITY, MOD=MODERATE, LS=LOW\_SELECTIVITY, PE=PAN\_ESSENTIAL. Type: Seq=sequence paralog, Part=partial homology, Func=functional analog. Three pairs evaluable under the previous frame (AKT1$\rightarrow$AKT2, CCNE1$\rightarrow$CCNE2, CDK4$\rightarrow$CDK6) no longer meet the $\geq$5-mutant threshold in any context and are not listed.}""",
 "table S6 caption"))

# ── apply ──────────────────────────────────────────────────────────
errors = []
for old, new, tag in R:
    n = text.count(old)
    if n != 1:
        errors.append(f"[{tag}] found {n}: {old[:80]!r}")
    else:
        text = text.replace(old, new, 1)

if errors:
    print("FAILED:")
    for e in errors:
        print(" ", e)
    raise SystemExit(1)

SP.write_text(text)
print(f"OK: {len(R)} replacements applied to supplementary.tex")
