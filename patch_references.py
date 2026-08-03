#!/usr/bin/env python3
"""patch_references.py — reference-round fixes (2026-08-01).

1. Background example list: separate citation-verified known SL pairs from
   catalogued pairs with indirect / direction-restricted evidence (fixes the
   EP300->CREBBP and PIK3CA->PIK3CB direction-citation issues).
2. Add Gonçalves, Ryan & Adams (Nat Rev Drug Discov 2026;25(1):22-38) —
   verified via PubMed PMID 40935881 — cited in Background and Discussion.
3. Methods: point to the data manifest with SHA256 checksums and download
   dates (data/README.md).
"""
from pathlib import Path

MS = Path(__file__).parent / "manuscript.tex"
ms = MS.read_text(encoding="utf-8")

EDITS = [
    # 1. Background example list: verified-known vs direction-restricted
    (
        r"Known examples include \textit{SMARCA4}$\rightarrow$\textit{SMARCA2} (SWI/SNF ATPase subunits) \cite{Hoffman2014}, \textit{ARID1A}$\rightarrow$\textit{ARID1B} (BAF complex subunits) \cite{Helming2014}, \textit{EP300}$\rightarrow$\textit{CREBBP} (histone acetyltransferases) \cite{Ogiwara2016}, \textit{PIK3CA}$\rightarrow$\textit{PIK3CB} (PI3K catalytic isoforms) \cite{SynLethDB2022}, \textit{AKT1}$\rightarrow$\textit{AKT2}, \textit{CCNE1}$\rightarrow$\textit{CCNE2}, \textit{CDK4}$\rightarrow$\textit{CDK6}, and \textit{MAP2K1}$\rightarrow$\textit{MAP2K2} \cite{SynLethDB2022}. The strength of evidence behind these catalogued pairs varies — from direct genotype-conditional experiments to redundancy and pharmacologic data — so we verified each primary citation and assigned every pair an evidence tier",
        r"Known examples include \textit{SMARCA4}$\rightarrow$\textit{SMARCA2} (SWI/SNF ATPase subunits) \cite{Hoffman2014}, \textit{ARID1A}$\rightarrow$\textit{ARID1B} (BAF complex subunits) \cite{Helming2014}, \textit{AKT1}$\rightarrow$\textit{AKT2} \cite{Najm2018}, and \textit{CDK4}$\rightarrow$\textit{CDK6} and \textit{MAP2K1}$\rightarrow$\textit{MAP2K2} \cite{Parrish2021}. Other catalogued pairs carry indirect or direction-restricted evidence: for \textit{EP300} and \textit{CREBBP} (histone acetyltransferases), direct experiments established the reciprocal CREBBP$\rightarrow$EP300 direction \cite{Ogiwara2016,Nie2021}; for \textit{PIK3CA}$\rightarrow$\textit{PIK3CB} (PI3K catalytic isoforms), the primary evidence supports PTEN$\rightarrow$PIK3CB \cite{Wee2008}; and \textit{CCNE1}$\rightarrow$\textit{CCNE2} is supported by developmental redundancy in mouse double knockouts \cite{Geng2003}. Because the evidence behind catalogued pairs varies in both type and direction, we verified each primary citation and assigned every pair an evidence tier",
        "background example list",
    ),
    # 2a. NRDD review cited in Background
    (
        "The logic is elegant; repeating it for other mutations has proven hard.",
        r"The logic is elegant; repeating it for other mutations has proven hard \cite{Goncalves2026}.",
        "NRDD cite background",
    ),
    # 2b. NRDD review cited in Discussion (translational barriers)
    (
        "That path required a dedicated bromodomain ligand, ternary-complex optimization, and extensive medicinal chemistry --- none of which our computational targetability descriptors attempt to model.",
        r"That path required a dedicated bromodomain ligand, ternary-complex optimization, and extensive medicinal chemistry --- none of which our computational targetability descriptors attempt to model \cite{Goncalves2026}.",
        "NRDD cite discussion",
    ),
    # 3. Methods: data-integrity manifest pointer
    (
        "Cell lines without expression or dependency data were excluded, yielding 1,208 cell lines.\n\n\\textbf{HGNC gene families.}",
        "Cell lines without expression or dependency data were excluded, yielding 1,208 cell lines. The exact files analysed (release filenames, sizes, SHA256 checksums, and per-file download dates) are tabulated in the repository data manifest (\\texttt{data/README.md}), so the identical inputs can be retrieved and verified.\n\n\\textbf{HGNC gene families.}",
        "data manifest sentence",
    ),
    # 4. bibliography entry (citation order: after Farmer2005)
    (
        "\\bibitem{Farmer2005} Farmer H, McCabe N, Lord CJ, Tutt ANJ, Johnson DA, Richardson TB, et al.\\ Targeting the DNA repair defect in BRCA mutant cells as a therapeutic strategy. \\textit{Nature}. 2005;434(7035):917--921.\n",
        "\\bibitem{Farmer2005} Farmer H, McCabe N, Lord CJ, Tutt ANJ, Johnson DA, Richardson TB, et al.\\ Targeting the DNA repair defect in BRCA mutant cells as a therapeutic strategy. \\textit{Nature}. 2005;434(7035):917--921.\n\n\\bibitem{Goncalves2026} Gonçalves E, Ryan CJ, Adams DJ. Synthetic lethality in cancer drug discovery: challenges and opportunities. \\textit{Nat Rev Drug Discov}. 2026;25(1):22--38.\n",
        "bibitem Goncalves2026",
    ),
]

for old, new, tag in EDITS:
    n = ms.count(old)
    if n != 1:
        raise SystemExit(f"anchor '{tag}' found {n} times (expected 1)")
    ms = ms.replace(old, new, 1)

MS.write_text(ms, encoding="utf-8")
print(f"OK: {len(EDITS)} reference edits applied.")
