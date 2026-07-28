"""
Mutation Type Stratification Analysis for Paralog-SL
=====================================================
Distinguishes truncating (high-impact LoF) vs missense mutations
per driver gene and evaluates whether paralog compensation signal
differs by mutation consequence type.

Hypothesis: Truncating mutations (frameshift, nonsense, splice-site)
cause complete loss of function and should elicit stronger paralog
compensation than missense mutations, which may retain partial function.
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score
from pathlib import Path
import re

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL
from data_loader import (
    load_dependency, load_expression, load_models,
    load_paralogs,
)
from pcs import ParalogCompensationScore

# ── Cancer definitions ─────────────────────────────────────────
ANALYSIS_CANCERS = {
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
    "Colorectal": [
        "Colorectal Adenocarcinoma",
    ],
    "Breast": [
        "Invasive Breast Carcinoma", "Breast Ductal Carcinoma In Situ",
        "Breast Neoplasm, NOS",
    ],
}

DRIVERS_ALL = ["TP53", "ARID1A", "BRCA1", "BRCA2", "PTEN", "PIK3CA",
               "KRAS", "PPP2R1A", "FBXW7", "CTNNB1", "RB1", "NF1",
               "ATR", "ATM", "SMARCA4", "EP300", "APC", "BRAF", "SMAD4",
               "ERBB2", "CDH1", "GATA3", "STK11", "KEAP1", "PIK3R1",
               "KMT2D"]


def load_full_mutations(path=None):
    """Load mutation data WITHOUT filtering to damaging only."""
    from config import DEPMAP_FILES
    path = path or DEPMAP_FILES["mutations"]
    df = pd.read_csv(path, low_memory=False)
    gene_col = next((c for c in ["HugoSymbol", "Hugo_Symbol", "Gene"] if c in df.columns), None)
    model_col = next((c for c in ["ModelID", "DepMap_ID"] if c in df.columns), None)
    variant_col = "VariantInfo" if "VariantInfo" in df.columns else None
    
    keep = [model_col, gene_col]
    if variant_col:
        keep.append(variant_col)
    
    df = df[keep].rename(columns={model_col: "DepMap_ID", gene_col: "Gene"})
    return df, variant_col


def classify_mutation_type(variant_str):
    """
    Classify a VEP VariantInfo string into mutation consequence category.
    
    Returns: 'truncating', 'missense', 'inframe', 'other', or 'wt'
    """
    if pd.isna(variant_str):
        return "wt"
    
    v = str(variant_str).lower()
    
    # Truncating = high-impact LoF
    truncating_terms = [
        "frameshift_variant", "stop_gained", "nonsense",
        "splice_acceptor_variant", "splice_donor_variant",
        "start_lost", "stop_lost"
    ]
    for term in truncating_terms:
        if term in v:
            return "truncating"
    
    # Missense
    if "missense_variant" in v or "protein_altering_variant" in v:
        return "missense"
    
    # Inframe indels
    if "inframe_deletion" in v or "inframe_insertion" in v:
        return "inframe"
    
    return "other"


def build_mutation_type_matrix(mutations_df, cell_lines, driver_genes):
    """
    Build mutation matrix with consequence classification.
    For each (cell_line, gene) pair, returns the most severe
    mutation type: truncating > inframe > missense > other > wt.
    """
    sub = mutations_df[
        mutations_df["DepMap_ID"].isin(cell_lines) &
        mutations_df["Gene"].isin(driver_genes)
    ]
    
    if sub.empty or "VariantInfo" not in sub.columns:
        return pd.DataFrame("wt", index=cell_lines, columns=driver_genes)
    
    # Classify each mutation
    sub = sub.copy()
    sub["mut_type"] = sub["VariantInfo"].apply(classify_mutation_type)
    
    # For each (cell_line, gene), take the most severe type
    severity_order = {"truncating": 4, "inframe": 3, "missense": 2, "other": 1, "wt": 0}
    sub["severity"] = sub["mut_type"].map(severity_order)
    
    # Group by cell line + gene, take max severity
    grouped = sub.groupby(["DepMap_ID", "Gene"])["severity"].max().reset_index()
    
    # Pivot to matrix
    matrix = grouped.pivot(index="DepMap_ID", columns="Gene", values="severity").fillna(0)
    matrix = matrix.reindex(index=cell_lines, columns=driver_genes, fill_value=0)
    
    # Convert back to labels
    severity_to_label = {0: "wt", 1: "other", 2: "missense", 3: "inframe", 4: "truncating"}
    for col in matrix.columns:
        matrix[col] = matrix[col].map(severity_to_label)
    
    return matrix


def compute_dd_per_mutation_type(dep, mut_type_matrix, driver, paralog, cell_lines):
    """
    Compute DD separately for truncating-MUT vs missense-MUT.
    Returns dict with dd_trunc, dd_miss, dd_all, n_trunc, n_miss, n_wt.
    """
    result = {"dd_trunc": 0.0, "dd_miss": 0.0, "dd_all": 0.0,
              "n_trunc": 0, "n_miss": 0, "n_wt": 0}
    
    if driver not in mut_type_matrix.columns or paralog not in dep.columns:
        return result
    
    mtypes = mut_type_matrix[driver]
    
    trunc_ids = mtypes[mtypes == "truncating"].index.tolist()
    miss_ids = mtypes[mtypes == "missense"].index.tolist()
    wt_ids = mtypes[mtypes == "wt"].index.tolist()
    
    trunc_ids = [c for c in trunc_ids if c in dep.index]
    miss_ids = [c for c in miss_ids if c in dep.index]
    wt_ids = [c for c in wt_ids if c in dep.index]
    
    result["n_trunc"] = len(trunc_ids)
    result["n_miss"] = len(miss_ids)
    result["n_wt"] = len(wt_ids)
    
    # Truncating vs WT
    if len(trunc_ids) >= 3 and len(wt_ids) >= 3:
        trunc_dep = dep.loc[trunc_ids, paralog].dropna()
        wt_dep_t = dep.loc[wt_ids, paralog].dropna()
        if len(trunc_dep) >= 3 and len(wt_dep_t) >= 3:
            # DD = mean(Chronos | WT) − mean(Chronos | MUT) (manuscript Eq. 1);
            # positive = stronger dependency in the mutant subgroup.
            result["dd_trunc"] = wt_dep_t.mean() - trunc_dep.mean()
            t_stat, p_val = stats.ttest_ind(trunc_dep, wt_dep_t, equal_var=False)
            result["dd_trunc_p"] = p_val
    
    # Missense vs WT
    if len(miss_ids) >= 3 and len(wt_ids) >= 3:
        miss_dep = dep.loc[miss_ids, paralog].dropna()
        wt_dep_m = dep.loc[wt_ids, paralog].dropna()
        if len(miss_dep) >= 3 and len(wt_dep_m) >= 3:
            result["dd_miss"] = wt_dep_m.mean() - miss_dep.mean()
            t_stat, p_val = stats.ttest_ind(miss_dep, wt_dep_m, equal_var=False)
            result["dd_miss_p"] = p_val
    
    # All MUT vs WT (standard DD)
    all_mut_ids = list(set(trunc_ids + miss_ids))
    mut_ids = [c for c in all_mut_ids if c in dep.index]
    if len(mut_ids) >= 3 and len(wt_ids) >= 3:
        mut_dep = dep.loc[mut_ids, paralog].dropna()
        wt_dep_all = dep.loc[wt_ids, paralog].dropna()
        if len(mut_dep) >= 3 and len(wt_dep_all) >= 3:
            result["dd_all"] = wt_dep_all.mean() - mut_dep.mean()
    
    return result


def run_mutation_type_analysis():
    """Main entry point."""
    print("=" * 70)
    print("  Mutation Type Stratification: Truncating vs Missense")
    print("=" * 70)
    
    # ── Load data ──
    dep = load_dependency()
    expr = load_expression()
    models = load_models()
    paralogs = load_paralogs()
    full_mutations, variant_col = load_full_mutations()
    
    print(f"\nMutation classification by {variant_col}")
    print(f"Total mutations: {len(full_mutations):,}")
    
    # ── Known SL pairs set ──
    known_set = set()
    for a, b in KNOWN_PARALOG_SL:
        known_set.add((a.upper(), b.upper()))
        known_set.add((b.upper(), a.upper()))
    
    # ── Overall mutation type distribution ──
    if "VariantInfo" in full_mutations.columns:
        full_mutations["mut_class"] = full_mutations["VariantInfo"].apply(classify_mutation_type)
        dist = full_mutations["mut_class"].value_counts()
        print("Mutation type distribution (all):")
        for tp, cnt in dist.items():
            print(f"  {tp:15s}: {cnt:>8,} ({cnt/len(full_mutations)*100:.1f}%)")
    
    all_results = {}
    
    for cancer_name, disease_patterns in ANALYSIS_CANCERS.items():
        print(f"\n{'─' * 70}")
        print(f"  {cancer_name}")
        print(f"{'─' * 70}")
        
        pat = "|".join(disease_patterns)
        mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
        cancer_models = models[mask].copy()
        
        cell_ids = cancer_models["DepMap_ID"].tolist()
        cell_ids = [c for c in cell_ids if c in dep.index and c in expr.index]
        
        if len(cell_ids) < 10:
            print(f"  Insufficient cell lines: {len(cell_ids)}")
            continue
        
        print(f"  Cell lines: {len(cell_ids)}")
        
        # Build mutation type matrix
        drivers = [d for d in DRIVERS_ALL if d in dep.columns]
        mut_type_matrix = build_mutation_type_matrix(
            full_mutations, cell_ids, drivers
        )
        
        # Count mutation types per driver
        trunc_counts = (mut_type_matrix == "truncating").sum()
        miss_counts = (mut_type_matrix == "missense").sum()
        drivers_with_trunc = trunc_counts[trunc_counts >= 3].index.tolist()
        drivers_with_miss = miss_counts[miss_counts >= 3].index.tolist()
        
        print(f"  Drivers with ≥3 truncating: {len(drivers_with_trunc)}")
        print(f"  Drivers with ≥3 missense:   {len(drivers_with_miss)}")
        
        # For each driver with sufficient samples, compute DD per mutation type
        comparison_rows = []
        detail_rows = []
        
        for driver in set(drivers_with_trunc + drivers_with_miss):
            driver_paralogs = paralogs[paralogs["gene_A"] == driver]["gene_B"].unique()
            valid_paralogs = [p for p in driver_paralogs 
                            if p in dep.columns and p in expr.columns]
            
            for paralog in valid_paralogs:
                dd_info = compute_dd_per_mutation_type(
                    dep, mut_type_matrix, driver, paralog, cell_ids
                )
                
                is_known = ((driver.upper(), paralog.upper()) in known_set)
                
                detail_rows.append({
                    "cancer": cancer_name,
                    "driver": driver,
                    "paralog": paralog,
                    "is_known_sl": is_known,
                    "dd_trunc": dd_info["dd_trunc"],
                    "dd_miss": dd_info["dd_miss"],
                    "dd_all": dd_info["dd_all"],
                    "n_trunc": dd_info["n_trunc"],
                    "n_miss": dd_info["n_miss"],
                    "n_wt": dd_info["n_wt"],
                })
        
        if not detail_rows:
            print(f"  No results")
            continue
        
        detail_df = pd.DataFrame(detail_rows)
        all_results[cancer_name] = detail_df
        
        # ── Summary: compare DD magnitude for truncating vs missense ──
        # Only consider pairs where both trunc and missense have valid DD
        valid = detail_df[
            (detail_df["dd_trunc"].abs() > 0.001) | 
            (detail_df["dd_miss"].abs() > 0.001)
        ]
        
        n_pairs = len(detail_df)
        n_valid = len(valid)
        
        # AUROC per mutation type
        yt = detail_df["is_known_sl"].astype(int).values
        
        auc_all = roc_auc_score(yt, detail_df["dd_all"].fillna(0)) if yt.sum() >= 2 else np.nan
        auc_trunc = roc_auc_score(yt, detail_df["dd_trunc"].fillna(0)) if yt.sum() >= 2 else np.nan
        auc_miss = roc_auc_score(yt, detail_df["dd_miss"].fillna(0)) if yt.sum() >= 2 else np.nan
        
        # Mean DD magnitude comparison
        mean_dd_trunc = detail_df["dd_trunc"].abs().mean()
        mean_dd_miss = detail_df["dd_miss"].abs().mean()
        
        print(f"  Total pairs: {n_pairs}")
        print(f"  Known SL:    {int(yt.sum())}")
        print(f"  DD AUROC — All:       {auc_all:.3f}" if not np.isnan(auc_all) else "  DD AUROC — All:       N/A")
        print(f"  DD AUROC — Truncating: {auc_trunc:.3f}" if not np.isnan(auc_trunc) else "  DD AUROC — Truncating: N/A")
        print(f"  DD AUROC — Missense:   {auc_miss:.3f}" if not np.isnan(auc_miss) else "  DD AUROC — Missense:   N/A")
        print(f"  Mean |DD| — Truncating: {mean_dd_trunc:.4f}")
        print(f"  Mean |DD| — Missense:   {mean_dd_miss:.4f}")
        
        # ── Top pairs where truncating >> missense ──
        if n_valid > 0:
            valid_e = valid.copy()
            valid_e["dd_diff"] = valid_e["dd_trunc"].abs() - valid_e["dd_miss"].abs()
            top_trunc = valid_e.nlargest(10, "dd_diff")
            
            print(f"\n  Top 10 pairs where truncating DD > missense DD:")
            for _, r in top_trunc.iterrows():
                flag = "★" if r["is_known_sl"] else "·"
                print(f"    {flag} {r['driver']:10s}→{r['paralog']:10s}  "
                      f"trunc={r['dd_trunc']:+.3f}  miss={r['dd_miss']:+.3f}  "
                      f"Δ={r['dd_diff']:+.3f}")
        
        # Save
        out_path = OUTPUT_DIR / f"muttype_{cancer_name.lower()}_results.csv"
        detail_df.to_csv(out_path, index=False)
    
    # ── Pan-cancer summary ──
    print(f"\n{'=' * 70}")
    print(f"  Pan-Cancer Mutation Type Summary")
    print(f"{'=' * 70}")
    
    # Combine all
    if all_results:
        combined = pd.concat(all_results.values(), ignore_index=True)
        combined.to_csv(OUTPUT_DIR / "muttype_all_results.csv", index=False)
        
        # Per cancer AUC comparison
        print(f"\n{'Cancer':20s} {'Pairs':>6s} {'Known':>6s} {'All':>7s} {'Trunc':>7s} {'Miss':>7s}")
        print(f"{'-' * 60}")
        for cancer_name, detail_df in all_results.items():
            yt = detail_df["is_known_sl"].astype(int).values
            ys_all = detail_df["dd_all"].fillna(0)
            ys_trunc = detail_df["dd_trunc"].fillna(0)
            ys_miss = detail_df["dd_miss"].fillna(0)
            
            auc_all = roc_auc_score(yt, ys_all) if yt.sum() >= 2 else np.nan
            auc_trunc = roc_auc_score(yt, ys_trunc) if yt.sum() >= 2 else np.nan
            auc_miss = roc_auc_score(yt, ys_miss) if yt.sum() >= 2 else np.nan
            
            def fmt(a): return f"{a:.3f}" if not np.isnan(a) else "N/A"
            print(f"{cancer_name:20s} {len(detail_df):>6d} {int(yt.sum()):>6d} "
                  f"{fmt(auc_all):>7s} {fmt(auc_trunc):>7s} {fmt(auc_miss):>7s}")
        
        # Overall: truncating vs missense DD magnitude
        print(f"\n  Overall DD magnitude comparison:")
        print(f"  Mean |DD truncating|: {combined['dd_trunc'].abs().mean():.4f}")
        print(f"  Mean |DD missense|:   {combined['dd_miss'].abs().mean():.4f}")
        t_stat, p_val = stats.ttest_rel(
            combined["dd_trunc"].abs(), combined["dd_miss"].abs(), 
            nan_policy="omit"
        )
        print(f"  Paired t-test: t={t_stat:.3f}, p={p_val:.4f}")
    
    print(f"\nResults saved to {OUTPUT_DIR}/muttype_*.csv")
    return all_results


if __name__ == "__main__":
    run_mutation_type_analysis()
