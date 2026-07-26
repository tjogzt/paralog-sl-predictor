#!/usr/bin/env python3
"""D1 batch 2: supplementary.tex Table S3 rebuild + in4mer mapping subsection;
manuscript.tex Chou2025 author fix. Atomic writes."""
import re
from pathlib import Path

# ───────────────────────── supplementary.tex ─────────────────────────
sp = Path("supplementary.tex")
s = sp.read_text()

NEW_S3 = r"""\centering
\scriptsize
\resizebox{\textwidth}{!}{%
\begin{tabular}{lllp{3.6cm}p{2.6cm}p{2.2cm}cccp{1.5cm}l}
\toprule
\textbf{Tier} & \textbf{Driver} & \textbf{Paralog} & \textbf{Assay (perturbation)} & \textbf{Model system} & \textbf{Validated direction} & \textbf{Dual?} & \textbf{Indep.} & \textbf{Direct SL} & \textbf{Inclusion} & \textbf{Key ref} \\
\midrule
A & AKT1 & AKT2 & Combinatorial CRISPR digenic KO & Cancer cell lines & AKT1$\rightarrow$AKT2 & Yes & Yes & Yes & Primary & Najm 2018 \\
A & CDK4 & CDK6 & Digenic KO (pgPEN library) & Cancer cell lines & CDK4$\rightarrow$CDK6 & Yes & Yes & Yes & Primary & Parrish 2021 \\
A & MAP2K1 & MAP2K2 & Digenic KO (pgPEN library) & Cancer cell lines & MAP2K1$\rightarrow$MAP2K2 & Yes & Yes & Yes & Primary & Parrish 2021 \\
\midrule
B & SMARCA4 & SMARCA2 & CRISPR KO conditioned on natural SMARCA4 mutation & SMARCA4-mutant cancer lines & SMARCA4$\rightarrow$SMARCA2 & No & Yes & Conditional & Primary & Hoffman 2014 \\
B & ARID1A & ARID1B & shRNA knockdown conditioned on natural ARID1A mutation & ARID1A-mutant cancer lines & ARID1A$\rightarrow$ARID1B & No & Yes & Conditional & Primary & Helming 2014 \\
\midrule
C & EP300 & CREBBP & p300 degradation / CRISPR in CREBBP-deficient lines & CREBBP-mutant cancer lines & Reciprocal only (CREBBP$\rightarrow$EP300) & No & Yes & Reciprocal & Secondary & Ogiwara 2016; Nie 2021 \\
C & PIK3CA & PIK3CB & shRNA / PI3K inhibitor in PTEN-deficient lines & PTEN-deficient cancer lines & PTEN$\rightarrow$PIK3CB only & No & Yes & No (other driver) & Secondary & Wee 2008 \\
C & CCNE1 & CCNE2 & Mouse developmental double knockout & Mouse embryo & CCNE1$\leftrightarrow$CCNE2 redundancy & Yes (mouse) & Yes & No (redundancy) & Secondary & Geng 2003 \\
C & FBXW7 & FBXW2 & DepMap computational analysis & 700+ cell lines & FBXW7$\rightarrow$FBXW2 & No & No & No (computational) & Secondary & DepMap \\
C & PPP2R1A & PPP2R1B & DepMap computational analysis & 700+ cell lines & PPP2R1A$\rightarrow$PPP2R1B & No & No & No (computational) & Secondary & DepMap \\
\midrule
comp. & BRCA1 & BRCA2 & PARP inhibition / genetic screen (functional analogs) & BRCA-mutant cancer lines & BRCA1/2$\rightarrow$PARP axis & No & Yes & Functional analog & Comparator & Bryant 2005 \\
comp. & STK11 & SIK1 & Mouse genetics, LKB1--SIK axis (partial homolog) & Mouse NSCLC models & STK11$\rightarrow$SIK1/3 axis & No & Yes & Pathway axis & Comparator & Hollstein 2019 \\
\bottomrule
\end{tabular}%
}
\caption{Twelve curated pairs with evidence provenance and tier assignment after verification of the primary citations. \textbf{Tier A} (3 pairs): direct genetic synthetic-lethal evidence from dual-gene perturbation (combinatorial/digenic CRISPR knockout). \textbf{Tier B} (2 pairs): natural-genotype conditional dependency --- the paralog is a demonstrated selective dependency in driver-mutant cells --- from single-gene perturbation with functional validation. The Tier A~$\cup$~Tier B set (5 pairs) constitutes the primary external benchmark. \textbf{Tier C} (5 pairs): indirect evidence only (reciprocal-direction-only, other-driver, developmental-redundancy, or DepMap-derived); excluded from the primary benchmark and reported separately. \textbf{Comparators} (2 pairs): mechanistic reference pairs that are not sequence paralogs (functional analog or pathway axis); used as specificity references only. \textbf{Dual?}: whether the evidence comes from simultaneous dual-gene perturbation. \textbf{Indep.}: evidence independent of DepMap CRISPR data. \textbf{Direct SL}: whether the evidence directly demonstrates synthetic lethality in the direction scored in this study (``Conditional'' = genotype-conditional dependency; ``Reciprocal'' = direct evidence for the reverse direction only). \textbf{Inclusion}: Primary = primary external benchmark; Secondary = reported separately; Comparator = specificity reference. All twelve pairs are HGNC gene-group sequence paralogs except the two comparators; validated directions refer to the driver$\rightarrow$paralog direction scored here.}
\end{table}"""

# Replace the whole S3 table block: from "\centering\footnotesize\resizebox" up to first "\end{table}"
pat = re.compile(r"\\centering\n\\footnotesize\n\\resizebox\{\\textwidth\}\{!\}\{%\n\\begin\{tabular\}\{llllllc\}.*?\\end\{table\}", re.S)
n = len(pat.findall(s))
assert n == 1, f"S3 block matches: {n}"
s = pat.sub(lambda m: NEW_S3, s, count=1)

# Insert in4mer subsection after Table S3, before Table S4 subsection
IN4MER = r"""\subsection{External combinatorial-CRISPR gold standards (in4mer)}

Esmaeili Anvar et al.\ (2024) meta-analyzed five combinatorial CRISPR paralog screens (Dede 2020, Gonatopoulos-Pournatzis 2020, Parrish 2021, Thompson 2021, Ito 2021) and defined 13 candidate cross-study gold-standard paralog synthetic lethals (paralog score $>0.25$, called a hit in more than one study): CNOT7--CNOT8, PITPNA--PITPNB, TIA1--TIAL1, SAR1A--SAR1B, PTP4A1--PTP4A2, GSK3A--GSK3B, CSNK2A1--CSNK2A2, CSNK1D--CSNK1E, MAPK1--MAPK3, ARFGEF1--ARFGEF2, HDAC1--HDAC2, ASF1A--ASF1B, and SLC25A28--SLC25A37. We assessed whether these pairs could serve as an external benchmark for DD. Only one of the 13 pairs (MAPK1--MAPK3) involves a gene in our 40-gene driver panel, and MAPK1 carries qualifying oncogenic hotspot mutations in zero DepMap 26Q1 cell lines (Table~S7), so no in4mer gold-standard pair is evaluable in our driver-mutation-conditioned framework. This reflects a genuine division of scope rather than missing data: the in4mer gold standards measure background-independent paralog buffering under dual knockout, whereas DD quantifies dependency shifts conditioned on naturally occurring driver mutations. We therefore report the in4mer set as a conceptual cross-reference rather than a benchmark, and we note that only 8 of the 13 pairs replicated as hits in at least three of the four in4mer screens (Chou et al.\ 2025).

"""
anchor = "\\end{table}\n\n\\subsection{Table S4: Pan-cancer summary statistics}"
assert s.count(anchor) == 1, f"S4 anchor count: {s.count(anchor)}"
s = s.replace(anchor, "\\end{table}\n\n" + IN4MER + "\\subsection{Table S4: Pan-cancer summary statistics}")

sp.write_text(s)
print("supplementary.tex: S3 rebuilt + in4mer subsection added")

# ───────────────────────── manuscript.tex ─────────────────────────
mp = Path("manuscript.tex")
t = mp.read_text()
old = ("\\bibitem{Chou2025} Chou J, Hart T. Z-scores outperform similar methods for analyzing CRISPR paralog "
       "synthetic lethality screens. \\textit{Genome Biol}. 2025;26(1):188.")
new = ("\\bibitem{Chou2025} Chou J, Esmaeili Anvar N, Elghaish R, Chen J, Hart T. Z-scores outperform "
       "similar methods for analyzing CRISPR paralog synthetic lethality screens. \\textit{Genome Biol}. "
       "2025;26(1):188.")
assert t.count(old) == 1
t = t.replace(old, new)
mp.write_text(t)
print("manuscript.tex: Chou2025 authors corrected")
