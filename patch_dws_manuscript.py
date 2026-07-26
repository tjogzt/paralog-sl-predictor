#!/usr/bin/env python3
"""DWS module unification sync: manuscript.tex after TW module moved to the
driver-rule mutation matrix + min5 frame. Every replacement must match exactly once."""
from pathlib import Path

MS = Path(__file__).parent / "manuscript.tex"
text = MS.read_text()

R = []

# ── Abstract ───────────────────────────────────────────────────────
R.append((
 r"""Dependency-window scoring prioritized ARID1A$\rightarrow$ARID1B (DWS~$=$~3.65) as the leading selective candidate""",
 r"""Dependency-window scoring prioritized ARID1A$\rightarrow$ARID1B (DWS~$=$~2.82) as the leading selective candidate""",
 "abstract DWS"))

# ── Results: DWS paragraph (line 197) ──────────────────────────────
R.append((
 r"""Only ARID1A$\rightarrow$ARID1B fell into HIGH\_SELECTIVITY (mean DWS~$=$~3.65 across five cancer contexts). This is what one would expect if ARID1A loss creates a specific dependency on ARID1B-containing BAF complexes \cite{Helming2014}, and the ranking held under alternative criteria: ARID1A$\rightarrow$ARID1B also ranked first by $|$DD$|$ alone and by selectivity alone within the dependency-window module (Supplementary Table~S6). Sensitivity analysis confirmed that the DWS floor parameter (0.01 vs.\ 0.05) had no effect on this pair (denominator~$=$~0.069), whereas nine other pairs shifted rank. BRCA1$\leftrightarrow$BRCA2 and PIK3R1$\rightarrow$CRKL were classified as PAN\_ESSENTIAL (pan-essential fraction 0.55 and 0.85), suggesting that targeting them directly would cause broad cellular toxicity. Thirteen paralog pairs maintained DWS~$>$~1.0 in $\ge$2 cancer contexts (Supplementary Fig.~S6).""",
 r"""Two pairs fell into HIGH\_SELECTIVITY: ARID1A$\rightarrow$ARID1B (mean DWS~$=$~2.82 across four cancer contexts; mean selectivity~$+0.28$) and the established synthetic-lethal pair SMARCA4$\rightarrow$SMARCA2 (4.87; $+0.18$), whose recovery acts as a positive control for the framework. This is what one would expect if ARID1A loss creates a specific dependency on ARID1B-containing BAF complexes \cite{Helming2014}: within the dependency-window module ARID1A$\rightarrow$ARID1B attained the highest mean selectivity and the second-highest mean $|$DD$|$ (0.270) of all evaluated pairs (Supplementary Table~S6). The highest raw DWS value was attained by NF1$\rightarrow$RASA2 (5.42), but this estimate rests on a near-zero pan-essential denominator (0.1\% of cell lines) with selectivity~$\approx$~0, so NF1$\rightarrow$RASA2 is not a selective candidate. BRCA1$\leftrightarrow$BRCA2 and PIK3R1$\rightarrow$CRKL were classified as PAN\_ESSENTIAL (pan-essential fraction 0.53--0.55 and 0.73), suggesting that targeting them directly would cause broad cellular toxicity. Nine paralog pairs maintained DWS~$>$~1.0 in $\ge$2 cancer contexts (Supplementary Fig.~S7).""",
 "DWS paragraph"))

# ── Results: composite score sentence (line 199) ───────────────────
R.append((
 r"""A composite preclinical score integrating DWS (0.40), selectivity (0.30), structure similarity (0.15), and druggability (0.15) ranked ARID1A$\rightarrow$ARID1B first (0.815; Fig.~4d; Table~2).""",
 r"""A composite preclinical score integrating normalized DWS (0.40), normalized selectivity (0.30), and structure-based targetability (0.30) ranked NF1$\rightarrow$RASA2 first (0.695) and ARID1A$\rightarrow$ARID1B second (0.631; Fig.~4d; Table~2); the NF1 score is driven by the near-zero pan-essential denominator noted above (selectivity~$\approx$~0), so ARID1A$\rightarrow$ARID1B remains the leading selective candidate.""",
 "composite score"))

# ── Fig 4b caption (line 206) ──────────────────────────────────────
R.append((
 r"""HIGH\_SELECTIVITY requires selectivity~$>$~0.15 and DWS~$>$~1.0; NF1$\rightarrow$RASA2 has high DWS (3.26) but near-zero selectivity (0.002), hence MODERATE.""",
 r"""HIGH\_SELECTIVITY requires selectivity~$>$~0.15 and DWS~$>$~1.0; NF1$\rightarrow$RASA2 has high DWS (5.42) but near-zero selectivity (0.005), hence MODERATE.""",
 "fig4b caption"))

# ── Fig 4d caption ─────────────────────────────────────────────────
R.append((
 r"""\textbf{d}, Composite preclinical prioritization score ranking the top 10 candidates. Scores integrate dependency window score (0.40), selectivity (0.30), structural similarity (0.15), and druggability (0.15). ARID1A$\rightarrow$ARID1B ranks first (0.815).""",
 r"""\textbf{d}, Composite preclinical prioritization score ranking the top 10 candidates. Scores integrate normalized dependency window score (0.40), normalized selectivity (0.30), and structure-based targetability (0.30). NF1$\rightarrow$RASA2 ranks first (0.695) through a near-zero pan-essential denominator that inflates its raw DWS (selectivity~$\approx$~0); ARID1A$\rightarrow$ARID1B (0.631) is the leading selective candidate.""",
 "fig4d caption"))

# ── Table 2 (candidates) caption + rows ────────────────────────────
R.append((
 r"""\caption{\textbf{Top paralog-SL candidates ranked by preclinical prioritization score.} DWS: mean dependency window score across 5 cancer contexts. Class: HS~$=$~HIGH\_SELECTIVITY, MOD~$=$~MODERATE, LS~$=$~LOW\_SELECTIVITY, PE~$=$~PAN\_ESSENTIAL. Struct: structural similarity score. Target.: composite preclinical score. Known SL pairs marked with $\star$. \textit{All candidates are computationally nominated and require experimental validation before clinical translation.}""",
 r"""\caption{\textbf{Top paralog-SL candidates ranked by preclinical prioritization score.} DWS: mean dependency window score across up to 5 cancer contexts on the $\geq$5-mutant frame. Class: HS~$=$~HIGH\_SELECTIVITY, MOD~$=$~MODERATE, LS~$=$~LOW\_SELECTIVITY, PE~$=$~PAN\_ESSENTIAL. Struct: structural similarity score. Target.: composite preclinical score. Known SL pairs marked with $\star$. \textit{NF1$\rightarrow$RASA2's top composite score is driven by a near-zero pan-essential denominator (0.1\% of cell lines) that inflates its DWS; with selectivity~$\approx$~0 it is not a selective candidate, and ARID1A$\rightarrow$ARID1B remains the leading selective candidate. All candidates are computationally nominated and require experimental validation before clinical translation.}""",
 "table 2 caption"))
R.append((
 r"""ARID1A   & ARID1B   & 3.65 & HS  & $\star$ & 0.952 & 0.815 \\
NF1      & RASA2    & 3.26 & MOD &         & 0.658 & 0.652 \\
KMT2D    & KMT2C    & 2.56 & MOD &         & 0.829 & 0.638 \\
PPP2R1A  & PPP2R1B  & 1.45 & LS  & $\star$ & 0.921 & 0.543 \\
EP300    & CREBBP   & 1.19 & MOD & $\star$ & 0.976 & 0.533 \\
PIK3CA   & PIK3CB   & 1.25 & LS  & $\star$ & 0.847 & 0.509 \\
FBXW7    & FBXW2    & 1.03 & LS  & $\star$ & 0.666 & 0.418 \\
TP53     & TP63     & 0.69 & MOD &         & 0.803 & 0.409 \\
STK11    & SIK1     & 1.02 & LS  & $\star$ & 0.710 & 0.417 \\
KRAS     & HRAS     & 0.14 & LS  &         & 0.885 & 0.392 \\""",
 r"""NF1      & RASA2    & 5.42 & MOD &         & 0.658 & 0.695 \\
ARID1A   & ARID1B   & 2.82 & HS  & $\star$ & 0.952 & 0.631 \\
EP300    & CREBBP   & 1.82 & MOD & $\star$ & 0.976 & 0.535 \\
KMT2D    & KMT2C    & 1.77 & MOD &         & 0.829 & 0.477 \\
PIK3CA   & PIK3CB   & 1.36 & LS  & $\star$ & 0.847 & 0.471 \\
PPP2R1A  & PPP2R1B  & 0.93 & LS  & $\star$ & 0.921 & 0.452 \\
TP53     & TP63     & 0.81 & MOD &         & 0.803 & 0.393 \\
KRAS     & HRAS     & 0.21 & LS  &         & 0.885 & 0.392 \\
FBXW7    & FBXW2    & 0.72 & LS  & $\star$ & 0.666 & 0.358 \\
STK11    & SIK1     & 0.39 & LS  & $\star$ & 0.710 & 0.334 \\""",
 "table 2 rows"))

# ── Discussion ─────────────────────────────────────────────────────
R.append((
 r"""Our DWS framework nominated ARID1A$\rightarrow$ARID1B (DWS~$=$~3.65) as the leading candidate while flagging BRCA1/2 paralogs as pan-essential.""",
 r"""Our DWS framework nominated ARID1A$\rightarrow$ARID1B (DWS~$=$~2.82) as the leading selective candidate while flagging BRCA1/2 paralogs as pan-essential.""",
 "discussion DWS"))
R.append((
 r"""Our top-ranked candidate has the largest lineage-level DD (0.386 in Ovarian cancer) and the highest DWS (3.65) among all evaluated pairs, yet its""",
 r"""Our top-ranked candidate has the largest lineage-level DD (0.386 in Ovarian cancer) and the leading dependency-window profile among novel selective candidates (DWS~$=$~2.82; mean selectivity~$+0.28$), yet its""",
 "ARID1A paradox DWS"))

# ── fix table cross-ref introduced in D2 patch ─────────────────────
R.append((
 r"""and the leading dependency-window and selectivity ranking (Table~1), providing orthogonal prioritization""",
 r"""and a top-tier dependency-window and selectivity profile (Table~2), providing orthogonal prioritization""",
 "table xref fix"))

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

MS.write_text(text)
print(f"OK: {len(R)} replacements applied to manuscript.tex")
