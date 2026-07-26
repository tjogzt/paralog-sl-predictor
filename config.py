"""
Paralog SL Predictor — Configuration
=====================================
All file paths, thresholds, and cancer-type definitions.
Edit this file to customize the analysis.
"""

from pathlib import Path

# ── Project root ──────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "output"
FIG_DIR = OUTPUT_DIR / "figures"

for d in [DATA_DIR, OUTPUT_DIR, FIG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── DepMap data (download from https://depmap.org/portal/download/) ──
DEPMAP_VERSION = "26Q1"
DEPMAP_FILES = {
    "dependency": DATA_DIR / "CRISPRGeneEffect.csv",
    "expression": DATA_DIR / "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "mutations":  DATA_DIR / "OmicsSomaticMutations.csv",
    "models":     DATA_DIR / "Model.csv",
}

# ── Ensembl paralog data ──────────────────────────────────────
ENSEMBL_PARALOG_FILE = DATA_DIR / "ensembl_paralogs.csv"

# ── SynLethDB known SL pairs ──────────────────────────────────
SYNLETHDB_FILE = DATA_DIR / "synlethdb_sl_pairs.csv"

# ── Gynecological cancer definitions (Oncotree names) ─────────
GYN_CANCER_TYPES = {
    "Ovarian": [
        "Ovarian Cancer", "Ovarian Epithelial Tumor",
        "Ovarian Serous Cystadenocarcinoma",
        "Ovarian Clear Cell Adenocarcinoma",
        "Ovarian Endometrioid Adenocarcinoma",
    ],
    "Endometrial": [
        "Endometrial Carcinoma", "Endometrial Cancer",
        "Uterine Serous Carcinoma", "Uterine Carcinosarcoma",
        "Endometrial Endometrioid Adenocarcinoma",
    ],
    "Cervical": [
        "Cervical Cancer", "Cervical Squamous Cell Carcinoma",
        "Cervical Adenocarcinoma",
    ],
    "Lung": [
        "Lung Adenocarcinoma", "Lung Squamous Cell Carcinoma",
        "Non-Small Cell Lung Cancer", "Small Cell Lung Cancer",
        "Lung Cancer",
    ],
    "Breast": [
        "Invasive Breast Carcinoma", "Breast Ductal Carcinoma In Situ",
        "Breast Neoplasm, NOS", "Breast Neoplasm, NOS",
    ],
}

# ── Driver genes of interest (from Doc1 Stage 1) ─────────────
DRIVER_GENES = {
    "Ovarian":     ["TP53", "BRCA1", "BRCA2", "ARID1A", "CCNE1",
                     "NF1", "RB1", "PTEN", "PIK3CA", "KRAS"],
    "Endometrial": ["PTEN", "ARID1A", "PIK3CA", "CTNNB1", "TP53",
                     "PPP2R1A", "KRAS", "PIK3R1", "FBXW7", "KMT2D"],
    "Cervical":    ["PIK3CA", "EP300", "FBXW7", "STK11", "ERBB2",
                     "MAPK1", "PTEN", "KRAS"],
    "Lung":        ["TP53", "KRAS", "EGFR", "STK11", "KEAP1",
                     "NF1", "BRAF", "PIK3CA", "ALK", "MET"],
    "Breast":      ["TP53", "PIK3CA", "PTEN", "BRCA1", "ERBB2",
                     "GATA3", "CDH1", "RB1", "NF1", "MAP3K1"],
}

# ── Driver mutation rules by gene class ───────────────────────
# A driver gene is "mutant" in a cell line only when the variant
# matches the gene's oncogenic mechanism:
#   TSG (tumor suppressor): LikelyLoF == True in the DepMap 26Q1
#       mutation annotation (DepMap-curated likely loss-of-function;
#       includes dominant-negative TP53 missense, excludes passengers)
#   ONC (oncogene): Hotspot == True in the same file (COSMIC/TCGA
#       recurrent oncogenic hotspots; excludes LoF events and
#       non-hotspot passengers)
# Gene classes follow the DepMap 26Q1 OncogeneHighImpact /
# TumorSuppressorHighImpact majority vote, with one documented
# override: GATA3 -> TSG (20/20 flag tie; frameshift-dominated
# loss-of-function spectrum in breast cancer, PMC PMID: 26928228).
GENE_DRIVER_CLASS = {
    # oncogenes (Hotspot rule)
    "KRAS": "ONC", "BRAF": "ONC", "PIK3CA": "ONC", "CTNNB1": "ONC",
    "ERBB2": "ONC", "MET": "ONC", "ALK": "ONC", "EGFR": "ONC",
    "MAPK1": "ONC", "CCNE1": "ONC", "CDK4": "ONC", "MAP2K1": "ONC",
    "AKT1": "ONC",
    # tumor suppressors (LikelyLoF rule)
    "TP53": "TSG", "BRCA1": "TSG", "BRCA2": "TSG", "ARID1A": "TSG",
    "NF1": "TSG", "RB1": "TSG", "PTEN": "TSG", "PIK3R1": "TSG",
    "FBXW7": "TSG", "KMT2D": "TSG", "EP300": "TSG", "STK11": "TSG",
    "KEAP1": "TSG", "GATA3": "TSG", "CDH1": "TSG", "MAP3K1": "TSG",
    "PPP2R1A": "TSG", "APC": "TSG", "SMAD4": "TSG", "ATM": "TSG",
    "ATR": "TSG", "SMARCA4": "TSG",
}


def driver_mutation_rule(gene: str) -> str:
    """Return 'TSG' or 'ONC' for a driver gene, else 'ANY'
    (gene absent from the class map -> no rule-based restriction)."""
    return GENE_DRIVER_CLASS.get(gene, "ANY")

# ── Known paralog-SL pairs (gold-standard positive set) ───────
# Sourced from: literature + SynLethDB + Doc2 case studies
KNOWN_PARALOG_SL = {
    ("SMARCA4", "SMARCA2"),   # ATPase subunits of SWI/SNF
    ("ARID1A", "ARID1B"),     # DNA-binding subunits of SWI/SNF
    ("BRCA1", "BRCA2"),       # Functional paralogs in HR repair
    ("EP300", "CREBBP"),      # Histone acetyltransferases
    ("PIK3CA", "PIK3CB"),     # p110 catalytic subunits
    ("AKT1", "AKT2"),         # AKT kinase isoforms
    ("STK11", "SIK1"),        # AMPK-family kinases (partial)
    ("FBXW7", "FBXW2"),       # E3 ubiquitin ligase F-box family
    ("PPP2R1A", "PPP2R1B"),   # PP2A scaffold subunits
    ("CCNE1", "CCNE2"),       # Cyclin E isoforms
    ("CDK4", "CDK6"),         # Cyclin-dependent kinases
    ("MEK1", "MEK2"),          # MAP2K1, MAP2K2
}

# ── Analysis thresholds ───────────────────────────────────────
# Primary analysis: >=5 mutant and >=5 WT cell lines per driver x lineage
# stratum (round-4 methods review). The former >=3/>=3 threshold is kept as
# a documented sensitivity analysis (artifacts:
# output/tables/TableS2_FullResults_min3.tsv, output/headline_metrics_min3.json).
MIN_MUT_SAMPLES = 5            # Minimum mutant cell lines per driver
MIN_WT_SAMPLES = 5             # Minimum WT cell lines per driver
MIN_DELTA_EXPR = 0.3           # Minimum |log2FC| for paralog expression change
PCS_THRESHOLD = 0.15           # PCS threshold for candidate nomination
SIGNIFICANCE_ALPHA = 0.05      # p-value threshold (after BH correction)

# ── Output ────────────────────────────────────────────────────
RESULTS_FILE = OUTPUT_DIR / "paralog_sl_candidates.csv"
SUMMARY_FILE = OUTPUT_DIR / "analysis_summary.txt"
