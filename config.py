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
MIN_MUT_SAMPLES = 3            # Minimum mutant cell lines per driver
MIN_WT_SAMPLES = 3             # Minimum WT cell lines per driver
MIN_DELTA_EXPR = 0.3           # Minimum |log2FC| for paralog expression change
PCS_THRESHOLD = 0.15           # PCS threshold for candidate nomination
SIGNIFICANCE_ALPHA = 0.05      # p-value threshold (after BH correction)

# ── Output ────────────────────────────────────────────────────
RESULTS_FILE = OUTPUT_DIR / "paralog_sl_candidates.csv"
SUMMARY_FILE = OUTPUT_DIR / "analysis_summary.txt"
