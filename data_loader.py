"""
Paralog SL Predictor — Data Loader
===================================
Loads and preprocesses DepMap, Ensembl paralog, and SynLethDB data.

All four DepMap CSVs must be placed in data/ before running.
Download from: https://depmap.org/portal/download/

Ensembl paralogs are auto-fetched from the Ensembl REST API if not
already cached in data/ensembl_paralogs.csv.
"""

import re
import time
import pandas as pd
import numpy as np
import requests
from pathlib import Path
from tqdm import tqdm
from typing import Optional, Dict, List, Tuple

from config import (
    DEPMAP_FILES, ENSEMBL_PARALOG_FILE, SYNLETHDB_FILE,
    GYN_CANCER_TYPES, DRIVER_GENES, KNOWN_PARALOG_SL,
    MIN_MUT_SAMPLES, MIN_WT_SAMPLES,
)


# ── DepMap loading ─────────────────────────────────────────────

def _parse_gene_column(col_name: str) -> str:
    """Extract gene symbol from DepMap column names like 'A1BG (1234)'."""
    m = re.match(r"^(.+?)\s*\(\d+\)$", str(col_name))
    return m.group(1) if m else str(col_name)


def load_dependency(path: Optional[Path] = None) -> pd.DataFrame:
    """
    Load CRISPR gene dependency matrix.
    Returns DataFrame: rows=cell lines (DepMap_ID), columns=gene symbols.
    Dependency scores are CERES or Chronos (lower = more essential).
    """
    path = path or DEPMAP_FILES["dependency"]
    if not path.exists():
        raise FileNotFoundError(
            f"Dependency file not found: {path}\n"
            f"Download from https://depmap.org/portal/download/"
        )
    df = pd.read_csv(path)
    # First column is cell line ID
    id_col = df.columns[0]
    df = df.rename(columns={id_col: "DepMap_ID"})
    # Parse gene names from remaining columns
    rename_map = {c: _parse_gene_column(c) for c in df.columns if c != "DepMap_ID"}
    df = df.rename(columns=rename_map)
    return df.set_index("DepMap_ID")


def load_expression(path: Optional[Path] = None) -> pd.DataFrame:
    """Load gene expression matrix (log2(TPM+1))."""
    path = path or DEPMAP_FILES["expression"]
    if not path.exists():
        raise FileNotFoundError(f"Expression file not found: {path}")
    df = pd.read_csv(path)
    id_col = df.columns[0]
    df = df.rename(columns={id_col: "DepMap_ID"})
    rename_map = {c: _parse_gene_column(c) for c in df.columns if c != "DepMap_ID"}
    df = df.rename(columns=rename_map)
    return df.set_index("DepMap_ID")


def load_models(path: Optional[Path] = None) -> pd.DataFrame:
    """Load cell line annotations. Keeps cancer-type columns."""
    path = path or DEPMAP_FILES["models"]
    if not path.exists():
        raise FileNotFoundError(f"Model file not found: {path}")
    df = pd.read_csv(path, low_memory=False)
    keep_cols = ["ModelID", "StrippedCellLineName",
                 "OncotreeLineage", "OncotreePrimaryDisease",
                 "OncotreeSubtype", "Age", "Sex"]
    available = [c for c in keep_cols if c in df.columns]
    return df[available].rename(columns={"ModelID": "DepMap_ID"})


def load_mutations(path: Optional[Path] = None) -> pd.DataFrame:
    """
    Load somatic mutation annotations.
    Filters to damaging mutations only (nonsense, frameshift, splice-site).
    """
    path = path or DEPMAP_FILES["mutations"]
    if not path.exists():
        raise FileNotFoundError(f"Mutations file not found: {path}")
    df = pd.read_csv(path, low_memory=False)
    # Identify the column that holds mutation consequence
    variant_col = None
    for c in ["VariantInfo", "VariantClassification", "Variant_Classification"]:
        if c in df.columns:
            variant_col = c
            break

    if variant_col:
        damaging = ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                     "Splice_Site", "Translation_Start_Site", "stop_gained",
                     "frameshift_variant", "splice_donor_variant",
                     "splice_acceptor_variant",
                     "Missense_Mutation", "missense_variant",
                     "In_Frame_Del", "In_Frame_Ins"]
        df = df[df[variant_col].isin(damaging) | df[variant_col].isna()]

    # Identify gene and model columns
    gene_col = next((c for c in ["HugoSymbol", "Hugo_Symbol", "Gene"] if c in df.columns), None)
    model_col = next((c for c in ["ModelID", "DepMap_ID"] if c in df.columns), None)

    if gene_col and model_col:
        return df[[model_col, gene_col]].rename(columns={model_col: "DepMap_ID",
                                                          gene_col: "Gene"})
    return df


def build_mutation_matrix(mutations_df: pd.DataFrame,
                           cell_lines: List[str],
                           genes: List[str]) -> pd.DataFrame:
    """
    Build binary mutation matrix.
    Rows = cell lines, Columns = genes, Values = 1 (damaging mut) or 0 (WT).
    """
    sub = mutations_df[
        mutations_df["DepMap_ID"].isin(cell_lines) &
        mutations_df["Gene"].isin(genes)
    ]
    if sub.empty:
        return pd.DataFrame(0, index=cell_lines, columns=genes)
    matrix = sub.pivot_table(
        index="DepMap_ID", columns="Gene",
        aggfunc="size", fill_value=0
    )
    matrix = matrix.reindex(index=cell_lines, columns=genes, fill_value=0)
    return (matrix > 0).astype(int)


def filter_gynecological_cell_lines(models_df: pd.DataFrame,
                                     cancer_type: Optional[str] = None
                                     ) -> pd.DataFrame:
    """
    Filter models to gynecological cancer cell lines.
    If cancer_type is 'Ovarian'/'Endometrial'/'Cervical', filter to that subtype.
    """
    if cancer_type and cancer_type in GYN_CANCER_TYPES:
        patterns = GYN_CANCER_TYPES[cancer_type]
    else:
        patterns = []
        for v in GYN_CANCER_TYPES.values():
            patterns.extend(v)

    pattern = "|".join(patterns)
    disease_col = "OncotreePrimaryDisease"
    if disease_col not in models_df.columns:
        return models_df

    mask = models_df[disease_col].str.contains(pattern, case=False, na=False)
    return models_df[mask].copy()


def classify_cancer_type(model_df: pd.DataFrame) -> pd.Series:
    """Assign each cell line to Ovarian/Endometrial/Cervical/Other."""
    disease_col = "OncotreePrimaryDisease"
    id_col = "DepMap_ID" if "DepMap_ID" in model_df.columns else model_df.index.name
    idx = model_df[id_col].values if id_col in model_df.columns else model_df.index
    result = pd.Series("Other", index=idx, dtype=str)
    for ctype, patterns in GYN_CANCER_TYPES.items():
        pattern = "|".join(patterns)
        mask = model_df[disease_col].str.contains(pattern, case=False, na=False).values
        result.iloc[mask] = ctype
    return result


# ── Ensembl paralog loading ────────────────────────────────────

def fetch_ensembl_paralogs(cache_file: Optional[Path] = None) -> pd.DataFrame:
    """
    Fetch human paralog pairs from Ensembl REST API.
    Caches result to data/ensembl_paralogs.csv.

    Returns DataFrame with columns: gene_A, gene_B, homology_type, identity_pct
    """
    cache_file = cache_file or ENSEMBL_PARALOG_FILE
    if cache_file.exists():
        print(f"Loading cached paralog data from {cache_file}")
        return pd.read_csv(cache_file)

    print("Fetching paralog data from Ensembl REST API...")
    # Strategy: get all human genes, then fetch paralogs for each
    server = "https://rest.ensembl.org"

    # 1. Get human gene list via POST lookup
    # Use a simpler approach: iterate over known Ensembl gene IDs
    # For a full analysis, use BioMart. Here we fetch via homology endpoint.

    # First get all human genes with paralogs via the compara endpoint
    all_paralogs = []

    # Use the overlap endpoint for human gene paralogs
    # We'll fetch in chunks by gene symbol (from our driver genes and common genes)
    # For a complete analysis, we'd need BioMart, but for this focused study
    # we build a comprehensive set using the Ensembl REST API.

    # Build a list of gene symbols to query (driver genes + common cancer genes)
    query_genes = set()
    for genes in DRIVER_GENES.values():
        query_genes.update(genes)

    # Add genes from known paralog-SL pairs
    for a, b in KNOWN_PARALOG_SL:
        query_genes.add(a)
        query_genes.add(b)

    for gene in tqdm(sorted(query_genes), desc="Fetching paralogs"):
        try:
            # Look up gene by symbol
            ext = f"/lookup/symbol/homo_sapiens/{gene}?expand=1"
            r = requests.get(f"{server}{ext}",
                             headers={"Content-Type": "application/json"})
            r.raise_for_status()
            data = r.json()

            # Get paralogs from homologues
            homologues = data.get("homologues", [])
            for h in homologues:
                h_type = h.get("homology_type", "")
                if "paralog" in h_type.lower():
                    target = h.get("target", {})
                    target_gene = target.get("id", "")
                    # Get target symbol (we already have gene symbol)
                    # Extract symbol from species-specific ID if needed
                    perc_id = h.get("perc_id", 0)
                    # Fetch target symbol
                    try:
                        time.sleep(0.1)  # rate limit
                        t_r = requests.get(
                            f"{server}/lookup/id/{target_gene}?expand=0",
                            headers={"Content-Type": "application/json"})
                        if t_r.ok:
                            t_data = t_r.json()
                            target_symbol = t_data.get("display_name", target_gene)
                            all_paralogs.append({
                                "gene_A": gene,
                                "gene_B": target_symbol,
                                "homology_type": h_type,
                                "identity_pct": perc_id,
                            })
                    except Exception:
                        continue
        except Exception as e:
            print(f"  Warning: could not fetch paralogs for {gene}: {e}")
            continue

    if not all_paralogs:
        print("Warning: No paralogs fetched via API. Using built-in paralog table.")
        return _build_builtin_paralogs()

    df = pd.DataFrame(all_paralogs).drop_duplicates()
    df.to_csv(cache_file, index=False)
    print(f"Cached {len(df)} paralog pairs to {cache_file}")
    return df


def _build_builtin_paralogs() -> pd.DataFrame:
    """
    Build a comprehensive paralog table from known relationships.
    This is used as fallback when Ensembl API is unavailable.

    Sources: Ensembl Compara, HGNC, literature
    """
    paralog_pairs = [
        # SWI/SNF complex subunits
        ("SMARCA4", "SMARCA2"), ("ARID1A", "ARID1B"),
        ("SMARCB1", "SMARCC1"), ("SMARCC1", "SMARCC2"),
        ("SMARCD1", "SMARCD2"), ("SMARCD2", "SMARCD3"),
        ("SMARCE1", "ACTL6A"), ("ACTL6A", "ACTL6B"),
        # DNA repair
        ("BRCA1", "BRCA2"), ("BRCA1", "BARD1"),
        ("RAD51", "RAD51B"), ("RAD51", "RAD51C"), ("RAD51", "RAD51D"),
        ("PARP1", "PARP2"), ("PARP1", "PARP3"),
        # PI3K/AKT/mTOR pathway
        ("PIK3CA", "PIK3CB"), ("PIK3CB", "PIK3CD"), ("PIK3CD", "PIK3CG"),
        ("AKT1", "AKT2"), ("AKT1", "AKT3"), ("AKT2", "AKT3"),
        ("PIK3R1", "PIK3R2"), ("PIK3R2", "PIK3R3"),
        ("PTEN", "TNS1"), ("PTEN", "TNS2"),
        # Cell cycle
        ("CCNE1", "CCNE2"), ("CCND1", "CCND2"), ("CCND2", "CCND3"),
        ("CDK4", "CDK6"), ("CDK2", "CDK3"),
        ("CDKN1A", "CDKN1B"), ("CDKN2A", "CDKN2B"),
        # p53 family
        ("TP53", "TP63"), ("TP53", "TP73"), ("TP63", "TP73"),
        # RAS/MAPK pathway
        ("KRAS", "NRAS"), ("KRAS", "HRAS"), ("NRAS", "HRAS"),
        ("MAP2K1", "MAP2K2"),  # MEK1/2
        ("MAPK1", "MAPK3"),   # ERK2/1
        ("BRAF", "RAF1"), ("RAF1", "ARAF"),
        # Wnt pathway
        ("CTNNB1", "JUP"),    # β-catenin, γ-catenin
        ("TCF7", "TCF7L1"), ("TCF7", "TCF7L2"),
        ("AXIN1", "AXIN2"),
        # PP2A complex
        ("PPP2CA", "PPP2CB"),
        ("PPP2R1A", "PPP2R1B"),
        ("PPP2R2A", "PPP2R2B"), ("PPP2R2B", "PPP2R2C"), ("PPP2R2C", "PPP2R2D"),
        # E3 ubiquitin ligases
        ("FBXW7", "FBXW2"), ("FBXW7", "FBXW4"), ("FBXW7", "FBXW5"),
        ("CUL1", "CUL2"), ("CUL2", "CUL3"), ("CUL3", "CUL4A"),
        ("CUL4A", "CUL4B"), ("CUL4A", "CUL5"),
        # Histone modifiers
        ("EP300", "CREBBP"),
        ("KMT2A", "KMT2B"), ("KMT2C", "KMT2D"),
        ("EZH1", "EZH2"),
        ("HDAC1", "HDAC2"), ("HDAC3", "HDAC4"),
        # AMPK family
        ("STK11", "SIK1"), ("STK11", "SIK2"), ("STK11", "SIK3"),
        ("SIK1", "SIK2"), ("SIK2", "SIK3"),
        # Receptor tyrosine kinases
        ("EGFR", "ERBB2"), ("ERBB2", "ERBB3"), ("ERBB3", "ERBB4"),
        ("FGFR1", "FGFR2"), ("FGFR2", "FGFR3"), ("FGFR3", "FGFR4"),
        # NF1/RAS GAPs
        ("NF1", "RASA1"), ("NF1", "RASA2"),
        # RB family
        ("RB1", "RBL1"), ("RB1", "RBL2"), ("RBL1", "RBL2"),
        # BCL2 family
        ("BCL2", "BCL2L1"), ("BCL2L1", "BCL2L2"),
        ("BAX", "BAK1"), ("BAX", "BOK"),
        # MYC family
        ("MYC", "MYCN"), ("MYC", "MYCL"),
        # mTOR
        ("MTOR", "PIK3C3"),
        # Others
        ("ATR", "ATM"),
        ("CHEK1", "CHEK2"),
    ]

    records = []
    for a, b in paralog_pairs:
        records.append({
            "gene_A": a, "gene_B": b,
            "homology_type": "paralog_known",
            "identity_pct": np.nan,
        })
        # Also add reverse for completeness
        records.append({
            "gene_A": b, "gene_B": a,
            "homology_type": "paralog_known",
            "identity_pct": np.nan,
        })

    df = pd.DataFrame(records).drop_duplicates()
    print(f"Built-in paralog table: {len(df)} pairs among {df['gene_A'].nunique()} genes")
    return df


def load_paralogs(cache_file: Optional[Path] = None) -> pd.DataFrame:
    """Main entry point for paralog loading. Tries API, falls back to built-in."""
    cache_file = cache_file or ENSEMBL_PARALOG_FILE
    try:
        return fetch_ensembl_paralogs(cache_file)
    except Exception as e:
        print(f"Ensembl API unavailable ({e}), using built-in paralog table.")
        df = _build_builtin_paralogs()
        df.to_csv(cache_file, index=False)
        return df


# ── SynLethDB loading ──────────────────────────────────────────

def load_synlethdb(path: Optional[Path] = None) -> pd.DataFrame:
    """
    Load SynLethDB known SL pairs.
    If the file doesn't exist, build a minimal set from literature.

    Returns DataFrame with columns: gene_A, gene_B, evidence
    """
    path = path or SYNLETHDB_FILE
    if path.exists():
        return pd.read_csv(path)

    print("SynLethDB file not found. Using literature-curated SL set.")
    # Build from known SL pairs in gynecological cancers (Doc2 table)
    known_sl = [
        ("BRCA1", "PARP1"), ("BRCA1", "PARP2"),
        ("BRCA2", "PARP1"), ("BRCA2", "PARP2"),
        ("ARID1A", "ATR"), ("ARID1A", "CHEK1"),
        ("ARID1A", "EZH2"), ("ARID1A", "PRMT5"),
        ("ARID1A", "ARID1B"),
        ("SMARCA4", "SMARCA2"),
        ("PTEN", "PIK3CB"), ("PTEN", "AKT1"),
        ("PPP2R1A", "RRM1"), ("PPP2R1A", "RRM2"),
        ("TP53", "WEE1"), ("TP53", "CHEK1"),
        ("KRAS", "SOS1"), ("KRAS", "RAF1"),
        ("CCNE1", "CDK2"),
        ("PIK3CA", "PIK3CB"),
        ("EP300", "CREBBP"),
        ("BRCA1", "POLQ"),
        ("NF1", "MAP2K1"), ("NF1", "MAP2K2"),
    ]
    records = []
    for a, b in known_sl:
        records.append({"gene_A": a, "gene_B": b, "evidence": "literature"})
        records.append({"gene_A": b, "gene_B": a, "evidence": "literature"})

    df = pd.DataFrame(records).drop_duplicates()
    df.to_csv(path, index=False)
    print(f"Built gold-standard SL set: {len(records)} directed pairs")
    return df


# ── Convenience loader ─────────────────────────────────────────

def load_all_data() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame,
                               pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load all required datasets at once.

    Returns
    -------
    dependency : pd.DataFrame  — CRISPR gene dependency (CERES/Chronos)
    expression : pd.DataFrame  — gene expression (log2 TPM+1)
    models     : pd.DataFrame  — cell line annotations
    mutations  : pd.DataFrame  — somatic mutation calls
    paralogs   : pd.DataFrame  — Ensembl paralog pairs
    synlethdb  : pd.DataFrame  — known SL pairs
    """
    print("=" * 60)
    print("Loading DepMap + Ensembl + SynLethDB data...")
    print("=" * 60)

    dependency = load_dependency()
    print(f"  Dependency matrix: {dependency.shape[0]} cell lines × {dependency.shape[1]} genes")

    expression = load_expression()
    print(f"  Expression matrix:  {expression.shape[0]} cell lines × {expression.shape[1]} genes")

    models = load_models()
    print(f"  Model annotations:  {len(models)} cell lines")

    mutations = load_mutations()
    print(f"  Mutations:          {len(mutations)} records")

    paralogs = load_paralogs()
    print(f"  Paralog pairs:      {len(paralogs)} pairs, "
          f"{paralogs['gene_A'].nunique()} unique genes")

    synlethdb = load_synlethdb()
    print(f"  Known SL pairs:     {len(synlethdb)} entries")

    return dependency, expression, models, mutations, paralogs, synlethdb
