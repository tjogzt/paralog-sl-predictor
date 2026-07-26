"""
Therapeutic Window Analysis for Paralog-SL Candidates
======================================================
Quantifies the therapeutic index for each paralog-SL candidate pair
by assessing paralog essentiality in driver-MUT vs driver-WT contexts.

Key metrics:
  1. Therapeutic Index (TI) = |DD| / (pan-essentiality + ε)
     - Higher TI = more selective killing in MUT context
  2. Selectivity Score = fraction of cell lines where paralog is
     essential (CERES < -0.5) in MUT vs WT
  3. Window Width = quantile difference in CERES distribution

This filters out paralogs that are pan-essential (poor drug targets)
and highlights those with genuinely selective essentiality.
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL, GYN_CANCER_TYPES, MIN_MUT_SAMPLES, MIN_WT_SAMPLES
from data_loader import (
    load_dependency, load_expression, load_models,
    load_mutations, load_paralogs, build_mutation_matrix,
)

# ── Essentiality thresholds ────────────────────────────────────
CERES_ESSENTIAL_THRESHOLD = -0.5  # CERES < -0.5 = essential
CERES_DEMETER_THRESHOLD = -1.0    # Strong essentiality

# ── Cancer contexts ────────────────────────────────────────────
THERAPEUTIC_CANCERS = {
    "Ovarian": GYN_CANCER_TYPES["Ovarian"],
    "Endometrial": GYN_CANCER_TYPES["Endometrial"],
    "Breast": GYN_CANCER_TYPES["Breast"],
    "Colorectal": ["Colorectal Adenocarcinoma"],
    "PanCancer": None,
}

# ── Key driver genes and their paralogs ────────────────────────
DRIVER_PARALOG_PAIRS = [
    ("ARID1A", "ARID1B"), ("PIK3CA", "PIK3CB"),
    ("BRCA1", "BRCA2"), ("BRCA2", "BRCA1"),
    ("EP300", "CREBBP"), ("PPP2R1A", "PPP2R1B"),
    ("FBXW7", "FBXW2"), ("STK11", "SIK1"),
    ("SMARCA4", "SMARCA2"), ("CCNE1", "CCNE2"),
    ("CDK4", "CDK6"), ("AKT1", "AKT2"),
    ("MAP2K1", "MAP2K2"),
    # De novo
    ("KRAS", "HRAS"), ("KRAS", "NRAS"),
    ("PIK3R1", "CRKL"), ("PTEN", "TNS2"),
    ("TP53", "TP63"), ("RB1", "RBL1"),
    ("NF1", "RASA2"), ("ATR", "ATM"),
    ("KMT2D", "KMT2C"), ("CDH1", "CDH2"),
    ("BRAF", "RAF1"),
]


def compute_pan_essentiality(dep, gene, all_cell_ids):
    """Compute pan-essentiality metrics for a gene."""
    vals = dep.loc[dep.index.isin(all_cell_ids), gene].dropna()
    if len(vals) == 0:
        return {"mean_ceres": 0, "frac_essential": 0, "frac_strong_essential": 0,
                "median_ceres": 0, "n": 0}
    
    return {
        "mean_ceres": vals.mean(),
        "median_ceres": vals.median(),
        "frac_essential": (vals < CERES_ESSENTIAL_THRESHOLD).mean(),
        "frac_strong_essential": (vals < CERES_DEMETER_THRESHOLD).mean(),
        "n": len(vals),
    }


def compute_therapeutic_window(dep, driver, paralog, mut_ids, wt_ids, all_ids):
    """
    Compute comprehensive therapeutic window metrics.
    """
    mut_ids_valid = [c for c in mut_ids if c in dep.index]
    wt_ids_valid = [c for c in wt_ids if c in dep.index]
    
    if len(mut_ids_valid) < MIN_MUT_SAMPLES or len(wt_ids_valid) < MIN_WT_SAMPLES:
        return None
    
    mut_vals = dep.loc[mut_ids_valid, paralog].dropna()
    wt_vals = dep.loc[wt_ids_valid, paralog].dropna()
    all_vals = dep.loc[dep.index.isin(all_ids), paralog].dropna()
    
    if len(mut_vals) < MIN_MUT_SAMPLES or len(wt_vals) < MIN_WT_SAMPLES:
        return None
    
    # Basic DD (manuscript Eq. 1: WT − MUT; positive = compensation in
    # driver-mutant lines). All downstream ranking uses |dd| (dd_abs), so
    # this sign convention is presentation-only but kept consistent with
    # pcs.py and the paralogSL R package (compute_dd).
    dd = wt_vals.mean() - mut_vals.mean()
    
    # Pan-essentiality
    pan_mean = all_vals.mean()
    pan_essential_frac = (all_vals < CERES_ESSENTIAL_THRESHOLD).mean()
    
    # Context-specific essentiality
    mut_essential_frac = (mut_vals < CERES_ESSENTIAL_THRESHOLD).mean()
    wt_essential_frac = (wt_vals < CERES_ESSENTIAL_THRESHOLD).mean()
    selectivity = mut_essential_frac - wt_essential_frac
    
    # Therapeutic Index: how much MORE essential in MUT, normalized by pan-essentiality
    # Add epsilon to avoid division by zero
    epsilon = 0.01
    pan_essentiality = max(abs(pan_mean), pan_essential_frac, epsilon)
    ti = abs(dd) / pan_essentiality
    
    # Alternative: Window Width
    # Difference in the 25th percentile between MUT and WT
    mut_q25 = mut_vals.quantile(0.25)
    wt_q25 = wt_vals.quantile(0.25)
    window_width = mut_q25 - wt_q25
    
    # Variance ratio (MUT vs WT)
    mut_var = mut_vals.var()
    wt_var = wt_vals.var()
    var_ratio = mut_var / wt_var if wt_var > 0 else 999
    
    # Welch's t-test for DD
    t_stat, p_val = stats.ttest_ind(mut_vals, wt_vals, equal_var=False)
    
    # Cohen's d
    pooled_std = np.sqrt((mut_var + wt_var) / 2)
    cohens_d = dd / pooled_std if pooled_std > 0 else 0
    
    return {
        "driver": driver,
        "paralog": paralog,
        "dd": dd,
        "dd_abs": abs(dd),
        "cohens_d": cohens_d,
        "p_value": p_val,
        
        # Pan-essentiality
        "paralog_mean_ceres": pan_mean,
        "paralog_pan_essential_frac": pan_essential_frac,
        
        # Context essentiality
        "mut_essential_frac": mut_essential_frac,
        "wt_essential_frac": wt_essential_frac,
        "selectivity": selectivity,  # MUT - WT essential fraction
        
        # Therapeutic indices
        "therapeutic_index": ti,
        "window_width": window_width,
        
        # Sample sizes
        "n_mut": len(mut_ids_valid),
        "n_wt": len(wt_ids_valid),
        "n_all": len(all_ids),
        
        # Variance
        "mut_variance": mut_var,
        "wt_variance": wt_var,
        "var_ratio": var_ratio,
        
        # Distribution quantiles
        "mut_q25": mut_q25,
        "wt_q25": wt_q25,
        "mut_median": mut_vals.median(),
        "wt_median": wt_vals.median(),
    }


def run_therapeutic_window_analysis():
    """Main entry point."""
    print("=" * 70)
    print("  Therapeutic Window Analysis for Paralog-SL")
    print(f"  {len(DRIVER_PARALOG_PAIRS)} pairs × {len(THERAPEUTIC_CANCERS)} contexts")
    print("=" * 70)
    
    # ── Load data ──
    dep = load_dependency()
    expr = load_expression()
    models = load_models()
    mutations = load_mutations()
    paralogs_df = load_paralogs()
    
    # Driver-mutant sets under the same gene-class-specific rules as the main
    # pipeline (TSG: LikelyLoF; oncogene: Hotspot) — unifies the TW/DD frame
    # with Table S2 (round-4 review).
    drivers = sorted({d for d, _ in DRIVER_PARALOG_PAIRS})
    mut_matrix = build_mutation_matrix(mutations, dep.index.tolist(), drivers)
    mut_sets = {g: set(mut_matrix.index[mut_matrix[g] == 1]) for g in mut_matrix.columns}
    
    # ── Known SL pairs ──
    known_set = set()
    for a, b in KNOWN_PARALOG_SL:
        known_set.add((a.upper(), b.upper()))
        known_set.add((b.upper(), a.upper()))
    
    all_results = []
    context_summaries = []
    
    for context_name, disease_patterns in THERAPEUTIC_CANCERS.items():
        print(f"\n{'─' * 70}")
        print(f"  Context: {context_name}")
        print(f"{'─' * 70}")
        
        if disease_patterns is None:
            model_subset = models.copy()
        else:
            pat = "|".join(disease_patterns)
            mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
            model_subset = models[mask].copy()
        
        all_cell_ids = model_subset["DepMap_ID"].tolist()
        valid_ids = [c for c in all_cell_ids if c in dep.index and c in expr.index]
        
        if len(valid_ids) < 20:
            print(f"  Insufficient cell lines: {len(valid_ids)}")
            continue
        
        print(f"  Cell lines: {len(valid_ids)}")
        
        context_results = []
        
        for driver, paralog in DRIVER_PARALOG_PAIRS:
            if driver not in dep.columns or paralog not in dep.columns:
                continue
            
            # Mutation status from the driver-rule matrix (same as main pipeline)
            driver_mut_set = mut_sets.get(driver, set())
            mut_ids = [c for c in valid_ids if c in driver_mut_set]
            wt_ids = [c for c in valid_ids if c not in driver_mut_set]
            
            if len(mut_ids) < MIN_MUT_SAMPLES:
                continue
            
            tw = compute_therapeutic_window(dep, driver, paralog, mut_ids, wt_ids, valid_ids)
            if tw is None:
                continue
            
            tw["context"] = context_name
            tw["is_known_sl"] = ((driver.upper(), paralog.upper()) in known_set)
            context_results.append(tw)
        
        if not context_results:
            print("  No results")
            continue
        
        context_df = pd.DataFrame(context_results)
        all_results.append(context_df)
        
        # Summary stats
        n_total = len(context_df)
        n_known = context_df["is_known_sl"].sum()
        ti_median = context_df["therapeutic_index"].median()
        n_selective = (context_df["selectivity"] > 0).sum()
        n_good_ti = (context_df["therapeutic_index"] > 0.5).sum()
        
        print(f"  Pairs analyzed: {n_total} ({int(n_known)} known)")
        print(f"  Median TI: {ti_median:.3f}")
        print(f"  Selectivity > 0: {n_selective}/{n_total} ({n_selective/n_total*100:.0f}%)")
        print(f"  TI > 0.5 (good window): {n_good_ti}/{n_total}")
        
        # Top pairs by therapeutic index
        top_ti = context_df.nlargest(15, "therapeutic_index")
        print(f"\n  Top therapeutic windows:")
        for _, r in top_ti.iterrows():
            flag = "★" if r["is_known_sl"] else "·"
            ti_stars = "★" if r["therapeutic_index"] > 2.0 else "☆" if r["therapeutic_index"] > 1.0 else "·"
            print(f"    {flag} {r['driver']:10s}→{r['paralog']:10s}  "
                  f"TI={r['therapeutic_index']:.3f}  "
                  f"DD={r['dd']:+.3f}  d={r['cohens_d']:+.2f}  "
                  f"sel={r['selectivity']:+.2f}  "
                  f"pan_ess={r['paralog_pan_essential_frac']:.2f}  "
                  f"{ti_stars}")
        
        # Save per-context results
        out_path = OUTPUT_DIR / f"therapeutic_window_{context_name.lower()}_results.csv"
        context_df.to_csv(out_path, index=False)
        
        context_summaries.append({
            "context": context_name,
            "n_lines": len(valid_ids),
            "n_pairs": n_total,
            "n_known": int(n_known),
            "median_ti": ti_median,
            "n_selective": n_selective,
            "n_good_window": n_good_ti,
        })
    
    # ── Pan-context synthesis ──
    print(f"\n{'=' * 70}")
    print(f"  Therapeutic Window: Pan-Context Summary")
    print(f"{'=' * 70}")
    
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(OUTPUT_DIR / "therapeutic_window_all_results.csv", index=False)
        
        # Summary table
        summary_df = pd.DataFrame(context_summaries)
        print(f"\n{'Context':20s} {'Lines':>6s} {'Pairs':>6s} {'Known':>6s} {'TI':>7s} {'Sel%':>6s} {'Good':>6s}")
        print("-" * 65)
        for _, r in summary_df.iterrows():
            sel_pct = r["n_selective"] / r["n_pairs"] * 100 if r["n_pairs"] > 0 else 0
            print(f"{r['context']:20s} {int(r['n_lines']):>6d} {int(r['n_pairs']):>6d} "
                  f"{int(r['n_known']):>6d} {r['median_ti']:>7.3f} {sel_pct:>5.0f}% "
                  f"{int(r['n_good_window']):>6d}")
        summary_df.to_csv(OUTPUT_DIR / "therapeutic_window_summary.csv", index=False)
        
        # ── Consistent high-TI pairs across contexts ──
        # Find pairs with TI > 1.0 in ≥2 contexts
        high_ti = combined[combined["therapeutic_index"] > 1.0]
        pair_context_counts = high_ti.groupby(["driver", "paralog"])["context"].nunique().reset_index()
        pair_context_counts.columns = ["driver", "paralog", "n_contexts"]
        pair_mean_ti = high_ti.groupby(["driver", "paralog"])["therapeutic_index"].mean().reset_index()
        
        consistent = pair_context_counts.merge(pair_mean_ti, on=["driver", "paralog"])
        consistent = consistent[consistent["n_contexts"] >= 2].sort_values(
            ["n_contexts", "therapeutic_index"], ascending=[False, False]
        )
        
        if len(consistent) > 0:
            print(f"\n  Consistent high-TI pairs (TI > 1 in ≥2 contexts):")
            for _, r in consistent.head(15).iterrows():
                is_known = ((r["driver"].upper(), r["paralog"].upper()) in known_set)
                flag = "★" if is_known else "·"
                print(f"    {flag} {r['driver']:10s}→{r['paralog']:10s}  "
                      f"contexts={int(r['n_contexts'])}  mean_TI={r['therapeutic_index']:.3f}")
        
        # ── Paralog toxicity classification ──
        # Pan-essential paralogs: essential in >50% of lines
        # Context-selective: essential mostly in MUT lines
        # Safe paralogs: never essential (good drug target)
        paralog_summary = combined.groupby(["driver", "paralog"]).agg(
            mean_ti=("therapeutic_index", "mean"),
            mean_dd=("dd_abs", "mean"),
            mean_selectivity=("selectivity", "mean"),
            mean_pan_essential=("paralog_pan_essential_frac", "mean"),
            n_contexts=("context", "nunique"),
        ).reset_index()
        
        # Classify
        def classify_paralog(row):
            if row["mean_pan_essential"] > 0.5:
                return "PAN_ESSENTIAL"  # Too toxic
            elif row["mean_selectivity"] > 0.15 and row["mean_ti"] > 1.0:
                return "HIGH_SELECTIVITY"  # Best candidates
            elif row["mean_selectivity"] > 0:
                return "MODERATE"
            else:
                return "LOW_SELECTIVITY"
        
        paralog_summary["classification"] = paralog_summary.apply(classify_paralog, axis=1)
        
        print(f"\n  Paralog safety classification:")
        for cls in ["HIGH_SELECTIVITY", "MODERATE", "LOW_SELECTIVITY", "PAN_ESSENTIAL"]:
            subset = paralog_summary[paralog_summary["classification"] == cls]
            print(f"    {cls:20s}: {len(subset)} paralogs")
            
            if cls == "HIGH_SELECTIVITY" and len(subset) > 0:
                for _, r in subset.iterrows():
                    is_known = ((r["driver"].upper(), r["paralog"].upper()) in known_set)
                    flag = "★" if is_known else "·"
                    print(f"      {flag} {r['driver']:10s}→{r['paralog']:10s}  "
                          f"TI={r['mean_ti']:.2f}  sel={r['mean_selectivity']:+.2f}  "
                          f"pan_ess={r['mean_pan_essential']:.2f}")
        
        paralog_summary.to_csv(OUTPUT_DIR / "therapeutic_window_paralog_classification.csv", index=False)
    
    print(f"\nResults saved to {OUTPUT_DIR}/therapeutic_window_*.csv")
    return all_results, context_summaries


if __name__ == "__main__":
    run_therapeutic_window_analysis()
