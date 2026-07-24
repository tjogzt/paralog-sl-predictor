"""
PRISM Drug Sensitivity Systematic Analysis for Paralog-SL
==========================================================
For each high-confidence paralog-SL candidate pair, scans all 1,482
PRISM drugs to identify compounds that selectively kill driver-mutant
cell lines (lower log2AUC in MUT vs WT).

Key metric: ΔAUC = mean(log2AUC | MUT) - mean(log2AUC | WT)
Negative ΔAUC = drug is more effective in MUT lines (selective killing)
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import false_discovery_control
from pathlib import Path

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL, GYN_CANCER_TYPES
from data_loader import (
    load_dependency, load_expression, load_models,
    load_mutations, load_paralogs, build_mutation_matrix,
)

# ── Top paralog-SL candidate pairs to analyze (from existing results) ──
TOP_CANDIDATES = [
    # Known SL pairs (gold standard)
    ("ARID1A", "ARID1B"),
    ("BRCA1", "BRCA2"),
    ("BRCA2", "BRCA1"),
    ("EP300", "CREBBP"),
    ("PIK3CA", "PIK3CB"),
    ("PPP2R1A", "PPP2R1B"),
    ("FBXW7", "FBXW2"),
    ("STK11", "SIK1"),
    
    # De novo candidates
    ("KRAS", "HRAS"),
    ("PIK3R1", "CRKL"),
    ("PTEN", "TNS2"),
    ("TP53", "TP63"),
    ("KRAS", "NRAS"),
    ("RB1", "RBL1"),
    ("NF1", "RASA2"),
    ("ATR", "ATM"),
    ("SMARCA4", "SMARCA2"),
    ("KMT2D", "KMT2C"),
    ("CDH1", "CDH2"),
]

# ── Cancer contexts for analysis ──
CANCER_CONTEXTS = {
    "PanCancer": None,  # All cell lines
    "Ovarian": GYN_CANCER_TYPES["Ovarian"],
    "Endometrial": GYN_CANCER_TYPES["Endometrial"],
    "Breast": GYN_CANCER_TYPES["Breast"],
    "Colorectal": ["Colorectal Adenocarcinoma"],
}

# MMR genes (to filter by MSI status if needed)
MMR_GENES = ["MLH1", "MSH2", "MSH6", "PMS2"]


def load_prism(path=None):
    """Load PRISM drug sensitivity data (log2AUC)."""
    path = path or DATA_DIR / "PRISM_log2AUC.csv"
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "DepMap_ID"})
    df = df.set_index("DepMap_ID")
    return df


def get_driver_mutation_status(mutations_df, driver_gene, cell_ids):
    """
    Return the set of cell lines with damaging mutations in driver_gene.
    Uses the same damaging mutation filter as the main pipeline.
    """
    # Use existing load_mutations to get filtered damaging mutations
    filtered = mutations_df[mutations_df["Gene"] == driver_gene]
    mut_ids = filtered["DepMap_ID"].unique().tolist()
    return [c for c in mut_ids if c in cell_ids]


def compute_drug_selectivity(prism_df, cell_ids, mut_ids, wt_ids, min_samples=3):
    """
    For each drug, compute differential sensitivity between MUT and WT cell lines.
    
    Returns DataFrame with columns: drug, mean_MUT, mean_WT, delta_AUC, p_value, ...
    """
    if len(mut_ids) < min_samples or len(wt_ids) < min_samples:
        return pd.DataFrame()
    
    # Filter to drugs with data in both groups
    prism_sub = prism_df.loc[prism_df.index.isin(cell_ids)]
    drugs = prism_sub.columns.tolist()
    
    results = []
    for drug in drugs:
        mut_vals = prism_sub.loc[prism_sub.index.isin(mut_ids), drug].dropna()
        wt_vals = prism_sub.loc[prism_sub.index.isin(wt_ids), drug].dropna()
        
        if len(mut_vals) < min_samples or len(wt_vals) < min_samples:
            continue
        
        mut_mean = mut_vals.mean()
        wt_mean = wt_vals.mean()
        delta = mut_mean - wt_mean  # negative = more sensitive in MUT
        
        # Welch's t-test
        t_stat, p_val = stats.ttest_ind(mut_vals, wt_vals, equal_var=False)
        
        # Effect size (Cohen's d)
        pooled_std = np.sqrt((mut_vals.var() + wt_vals.var()) / 2)
        cohens_d = delta / pooled_std if pooled_std > 0 else 0.0
        
        # Also compute rank-based selectivity
        # Fraction of MUT lines in bottom 25% sensitivity vs WT
        all_vals = pd.concat([mut_vals, wt_vals])
        q25 = all_vals.quantile(0.25)
        mut_sensitive_frac = (mut_vals <= q25).mean()
        wt_sensitive_frac = (wt_vals <= q25).mean()
        enrichment = mut_sensitive_frac - wt_sensitive_frac
        
        results.append({
            "drug": drug,
            "mean_mut": mut_mean,
            "mean_wt": wt_mean,
            "delta_auc": delta,
            "cohens_d": cohens_d,
            "p_value": p_val,
            "t_stat": t_stat,
            "n_mut": len(mut_vals),
            "n_wt": len(wt_vals),
            "mut_sensitive_frac": mut_sensitive_frac,
            "wt_sensitive_frac": wt_sensitive_frac,
            "enrichment": enrichment,
        })
    
    if not results:
        return pd.DataFrame()
    
    result_df = pd.DataFrame(results)
    result_df["abs_delta"] = result_df["delta_auc"].abs()
    result_df["abs_cohens_d"] = result_df["cohens_d"].abs()
    
    # BH correction
    if len(result_df) >= 2:
        result_df["bh_q"] = false_discovery_control(result_df["p_value"].fillna(1.0).values)
    
    return result_df.sort_values("delta_auc")  # Most negative = most selective


def run_prism_analysis():
    """Main entry point."""
    print("=" * 70)
    print("  PRISM Drug Sensitivity: Paralog-SL Selective Killing Scan")
    print(f"  {len(TOP_CANDIDATES)} candidate pairs × {len(CANCER_CONTEXTS)} contexts")
    print("=" * 70)
    
    # ── Load data ──
    prism = load_prism()
    dep = load_dependency()
    expr = load_expression()
    models = load_models()
    mutations = load_mutations()
    paralogs = load_paralogs()
    
    print(f"\nPRISM data: {len(prism)} cell lines × {len(prism.columns)} drugs")
    print(f"Dependency data: {len(dep)} cell lines")
    
    # ── Known SL pairs set ──
    known_set = set()
    for a, b in KNOWN_PARALOG_SL:
        known_set.add((a.upper(), b.upper()))
        known_set.add((b.upper(), a.upper()))
    
    all_hits = []
    
    for context_name, disease_patterns in CANCER_CONTEXTS.items():
        print(f"\n{'─' * 70}")
        print(f"  Context: {context_name}")
        print(f"{'─' * 70}")
        
        # Filter cell lines for this context
        if disease_patterns is None:
            # PanCancer: all cell lines with data in both dep+prism
            model_subset = models.copy()
        else:
            pat = "|".join(disease_patterns)
            mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
            model_subset = models[mask].copy()
        
        all_cell_ids = model_subset["DepMap_ID"].tolist()
        # Intersect with dep, expr, and prism
        valid_ids = [c for c in all_cell_ids 
                    if c in dep.index and c in expr.index and c in prism.index]
        
        if len(valid_ids) < 20:
            print(f"  Insufficient cell lines with PRISM data: {len(valid_ids)}")
            continue
        
        print(f"  Cell lines with PRISM data: {len(valid_ids)}")
        
        context_hits = []
        
        for driver_gene, paralog_gene in TOP_CANDIDATES:
            if driver_gene not in dep.columns or paralog_gene not in dep.columns:
                continue
            
            # Get mutation status
            mut_ids = get_driver_mutation_status(mutations, driver_gene, valid_ids)
            wt_ids = [c for c in valid_ids if c not in mut_ids]
            
            if len(mut_ids) < 3 or len(wt_ids) < 3:
                continue
            
            is_known = ((driver_gene.upper(), paralog_gene.upper()) in known_set)
            
            # Compute DD for reference
            dep_mut = dep.loc[dep.index.isin(mut_ids), paralog_gene].dropna()
            dep_wt = dep.loc[dep.index.isin(wt_ids), paralog_gene].dropna()
            dd = dep_mut.mean() - dep_wt.mean() if len(dep_mut) >= 3 and len(dep_wt) >= 3 else 0.0
            
            # Compute drug selectivity
            drug_results = compute_drug_selectivity(prism, valid_ids, mut_ids, wt_ids)
            
            if drug_results.empty:
                continue
            
            # Add pair info
            drug_results["driver"] = driver_gene
            drug_results["paralog"] = paralog_gene
            drug_results["context"] = context_name
            drug_results["is_known_sl"] = is_known
            drug_results["dd"] = dd
            drug_results["n_mut_total"] = len(mut_ids)
            drug_results["n_wt_total"] = len(wt_ids)
            
            # ── Significant hits (BH q < 0.25 and negative delta) ──
            sig_hits = drug_results[
                (drug_results["bh_q"] < 0.25) & 
                (drug_results["delta_auc"] < 0) & 
                (drug_results["enrichment"] > 0.1)
            ]
            
            if len(sig_hits) > 0:
                for _, hit in sig_hits.sort_values("delta_auc").head(5).iterrows():
                    context_hits.append({
                        "driver": driver_gene,
                        "paralog": paralog_gene,
                        "context": context_name,
                        "is_known": is_known,
                        "dd": dd,
                        "drug": hit["drug"],
                        "delta_auc": hit["delta_auc"],
                        "cohens_d": hit["cohens_d"],
                        "bh_q": hit["bh_q"],
                        "enrichment": hit["enrichment"],
                        "n_mut": int(hit["n_mut"]),
                        "n_wt": int(hit["n_wt"]),
                        "mut_sens_frac": hit["mut_sensitive_frac"],
                        "wt_sens_frac": hit["wt_sensitive_frac"],
                    })
            
            # Save full results per pair
            all_hits.append(drug_results)
            
            n_sig = len(sig_hits)
            top_drug = drug_results.iloc[0]["drug"] if not drug_results.empty else "N/A"
            top_delta = drug_results.iloc[0]["delta_auc"] if not drug_results.empty else 0
            flag = "★" if is_known else "·"
            print(f"  {flag} {driver_gene:10s}→{paralog_gene:10s}  "
                  f"DD={dd:+.3f}  MUT={len(mut_ids)} WT={len(wt_ids)}  "
                  f"sig_drugs={n_sig}  top={top_drug[:25]:25s} Δ={top_delta:+.3f}")
        
        # ── Top hits summary for this context ──
        if context_hits:
            hits_df = pd.DataFrame(context_hits)
            hits_df = hits_df.sort_values("delta_auc")
            print(f"\n  Top selective drugs in {context_name}:")
            for _, h in hits_df.head(15).iterrows():
                flag = "★" if h["is_known"] else "·"
                print(f"    {flag} {h['driver']:10s}→{h['paralog']:10s}  "
                      f"{h['drug']:30s}  ΔAUC={h['delta_auc']:+.3f}  "
                      f"d={h['cohens_d']:+.2f}  q={h['bh_q']:.3f}")
    
    # ── Pan-cancer summary ──
    print(f"\n{'=' * 70}")
    print(f"  PRISM Analysis: Global Hits Summary")
    print(f"{'=' * 70}")
    
    if all_hits:
        combined = pd.concat(all_hits, ignore_index=True)
        
        # Save full results
        combined.to_csv(OUTPUT_DIR / "prism_full_results.csv", index=False)
        
        # Top hits across all pairs and contexts
        sig_all = combined[
            (combined["bh_q"] < 0.25) & 
            (combined["delta_auc"] < 0) &
            (combined["enrichment"] > 0.1)
        ]
        
        if len(sig_all) > 0:
            sig_summary = sig_all.sort_values("delta_auc").head(30)
            
            print(f"\nTop 30 selective drugs (pan-cancer):")
            print(f"{'Pair':25s} {'Context':15s} {'Drug':30s} {'ΔAUC':>8s} {'d':>6s} {'q':>6s}")
            print("-" * 100)
            for _, h in sig_summary.iterrows():
                pair = f"{h['driver']}→{h['paralog']}"
                flag = "★" if h["is_known_sl"] else "·"
                print(f"{flag} {pair:23s} {h['context']:15s} "
                      f"{h['drug'][:30]:30s} {h['delta_auc']:>+8.3f} "
                      f"{h['cohens_d']:>+6.2f} {h['bh_q']:>6.3f}")
            
            sig_summary.to_csv(OUTPUT_DIR / "prism_top_hits.csv", index=False)
        
        print(f"\nTotal pairs analyzed: {combined[['driver','paralog','context']].drop_duplicates().shape[0]}")
        print(f"Significant selective drugs: {len(sig_all)}")
    else:
        print("  No results generated.")
    
    print(f"\nResults saved to {OUTPUT_DIR}/prism_*.csv")
    return all_hits


if __name__ == "__main__":
    run_prism_analysis()
