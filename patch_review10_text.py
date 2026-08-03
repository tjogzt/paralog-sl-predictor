#!/usr/bin/env python3
"""Review round '10. Methods & statistics' — text-only patch (Patch A).
Items: unlabeled-controls terminology, mutation-rule explicitness, power PU
caveat, Chronos CNV note, AI-use statement, 5-part Data availability,
cover-letter Opposed Reviewers removal. No numbers changed.
"""
import sys
from pathlib import Path

MS = Path("manuscript.tex")
CL = Path("cover_letter.md")

MS_FIX = [
    # 1. controls terminology (x2) — no matching was performed (Methods L266)
    ("all on the same benchmark universe: the twelve curated pairs plus matched unlabeled controls, evaluated",
     "all on the same benchmark universe: the twelve curated pairs plus unlabeled paralog controls, evaluated"),
    ("\\textbf{Curated pairs and matched unlabeled controls.}",
     "\\textbf{Curated pairs and unlabeled paralog controls.}"),
    # 2. mutation rules: iff + at least one (matches code: any qualifying variant)
    ("for tumor suppressors (TSG), a line is mutant only if it carries a DepMap-curated likely loss-of-function variant (\\texttt{LikelyLoF}, which includes dominant-negative TP53 missense); for oncogenes, only if it carries a recurrent oncogenic hotspot variant (\\texttt{Hotspot})",
     "for tumor suppressors (TSG), a line is mutant if and only if it carries at least one DepMap-curated likely loss-of-function variant (\\texttt{LikelyLoF}, which includes dominant-negative TP53 missense); for oncogenes, if and only if it carries at least one recurrent oncogenic hotspot variant (\\texttt{Hotspot})"),
    # 2b. GATA3 override basis + severity-collapse rule
    ("(one documented override: GATA3 as TSG)",
     "(one documented override: GATA3 as TSG, resolving a 20/20 majority-vote tie given its frameshift-dominated lesion spectrum; GATA3 drives none of the curated benchmark pairs)"),
    ("The complete per-gene variant classification table (old non-silent definition vs.\\ class-specific rule, cell line counts, and variant-type distribution) is provided as Supplementary Table~S6.",
     "The complete per-gene variant classification table (old non-silent definition vs.\\ class-specific rule, cell line counts, and variant-type distribution) is provided as Supplementary Table~S6. In the mutation-type analysis, multiple variants in the same gene per line are collapsed to the most severe consequence."),
    # 3. power analysis PU caveat
    ("(binormal model, Hanley--McNeil variance; \\texttt{power\\_analysis.py}). To compensate for this limitation",
     "(binormal model, Hanley--McNeil variance; \\texttt{power\\_analysis.py}); because these calculations assume independent entries and verified true negatives, they serve as rough planning values under the positive--unlabeled design rather than exact requirements. To compensate for this limitation"),
    # 4. Chronos screen-level CNV correction note
    ("and multivariate regression controlling for CNV left DD estimates virtually unchanged (mean $|\\Delta\\text{DD}|=0.005$).",
     "and multivariate regression controlling for CNV left DD estimates virtually unchanged (mean $|\\Delta\\text{DD}|=0.005$); because Chronos scores already model out screen-level copy-number bias, this adjustment targets residual paralog-specific confounding."),
    # 5. AI-use statement (journal policy: document LLM use; LLM not an author)
    ("We thank the DepMap, CPTAC, TCGA, PRISM, and cBioPortal teams for making their data publicly available. We also acknowledge the SynLethDB, HGNC, and UniProt teams for curated resources.",
     "We thank the DepMap, CPTAC, TCGA, PRISM, and cBioPortal teams for making their data publicly available. We also acknowledge the SynLethDB, HGNC, and UniProt teams for curated resources. During the preparation of this work, the authors used Kimi (Moonshot AI) for language editing, citation cross-checking, and code-refactoring assistance. All AI-generated suggestions were reviewed, tested, and approved by the authors, who take full responsibility for the integrity and content of this publication."),
    # 6. Data availability -> five labelled blocks (all original facts/DOIs kept)
    ("\\subsection*{Data availability}\nAll data are publicly available: DepMap 26Q1 from the DepMap Portal \\cite{DepMap2026}; CPTAC proteomics from cBioPortal \\cite{Cerami2012,Gao2013} and the Proteomic Data Commons \\cite{PDC}; TCGA PanCan from UCSC Xena \\cite{Goldman2020}; PRISM from the DepMap PRISM portal \\cite{Corsello2020}. HGNC gene families from the HGNC repository \\cite{HGNC}. Processed datasets with cell line-level DD values and all supplementary tables are archived at Zenodo (concept DOI \\url{https://doi.org/10.5281/zenodo.21502030}, which resolves to the latest release; the submission snapshot is release v1.3.2 with version DOI \\url{https://doi.org/10.5281/zenodo.21634700}) and the GitHub repository below.",
     "\\subsection*{Data availability}\n\\noindent\\textbf{Raw third-party datasets.} DepMap 26Q1 from the DepMap Portal \\cite{DepMap2026} (exact release filenames, sizes, SHA256 checksums, and per-file download dates in the repository data manifest, \\texttt{data/README.md}); CPTAC proteomics from cBioPortal \\cite{Cerami2012,Gao2013} and the Proteomic Data Commons \\cite{PDC}; TCGA PanCan from UCSC Xena \\cite{Goldman2020}; PRISM from the DepMap PRISM portal \\cite{Corsello2020}; HGNC gene families from the HGNC repository \\cite{HGNC}.\n\n\\noindent\\textbf{Processed analysis-ready data.} Cell line-level DD values and all derived analysis tables are archived at Zenodo (concept DOI \\url{https://doi.org/10.5281/zenodo.21502030}, which resolves to the latest release; the submission snapshot is release v1.3.2 with version DOI \\url{https://doi.org/10.5281/zenodo.21634700}) and in the GitHub repository below.\n\n\\noindent\\textbf{Machine-readable supplementary tables.} TSV mirrors of Supplementary Tables S1--S12 are provided as Additional files 2--12, each documented in Additional file 1.\n\n\\noindent\\textbf{Code.} See Code availability below.\n\n\\noindent\\textbf{Archived manuscript-specific release.} The submission snapshot is GitHub tag \\texttt{v1.3.2} (commit \\texttt{66fc633}) with version DOI \\url{https://doi.org/10.5281/zenodo.21634700} under the concept DOI above; the accepted-version snapshot will be deposited identically."),
]

CL_FIX = [
    ("**Opposed Reviewers:** None.\n", ""),
]


def apply(path, fixes):
    text = path.read_text(encoding="utf-8")
    n = 0
    for old, new in fixes:
        c = text.count(old)
        if c != 1:
            print(f"FAIL [{path.name}] count={c}: {old[:100]}")
            return False
        text = text.replace(old, new, 1)
        n += 1
    path.write_text(text, encoding="utf-8")
    print(f"OK [{path.name}]: {n} fixes")
    return True


ok = apply(MS, MS_FIX) & apply(CL, CL_FIX)
sys.exit(0 if ok else 1)
