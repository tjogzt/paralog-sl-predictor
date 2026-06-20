"""
CPTAC Multi-Cancer Proteomics Validation
=========================================
Expands protein-level paralog co-variation analysis across all
available CPTAC cohorts on cBioPortal.

For each cohort, fetches protein abundance data for all genes
in the paralog-SL analysis and computes Pearson correlations
between paralog pairs.

Cohorts covered:
  - BRCA (Breast Invasive Carcinoma)
  - COAD/READ (Colorectal)
  - LUAD (Lung Adenocarcinoma)
  - HNSCC (Head and Neck)
  - ccRCC (Clear Cell Renal)
  - GBM (Glioblastoma)
  - PDAC (Pancreatic)
  - PRAD (Prostate)
  - LSCC (Lung Squamous)
  - UCEC (Endometrial) — if not already cached

Also includes existing OV (PDC) and UCEC results for comparison.
"""

import requests
import json
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import false_discovery_control
from pathlib import Path
import time
import sys

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL

# ── cBioPortal CPTAC cohorts ───────────────────────────────────
# Verified 2025-05-26 against live API at https://www.cbioportal.org/api
# POST method required; sampleListId matches molecular profile
CPTAC_STUDIES = {
    "BRCA": {
        "study": "brca_cptac_2020",
        "profile": "brca_cptac_2020_protein_quantification",
        "sample_list": "brca_cptac_2020_protein_quantification",
        "n_expected": 122,
    },
    "COAD": {
        "study": "coad_cptac_2019",
        "profile": "coad_cptac_2019_protein_quantification",
        "sample_list": "coad_cptac_2019_protein_quantification",
        "n_expected": 95,
    },
    "LUAD": {
        "study": "luad_cptac_2020",
        "profile": "luad_cptac_2020_protein_quantification",
        "sample_list": "luad_cptac_2020_protein_quantification",
        "n_expected": 84,
    },
    "GBM": {
        "study": "gbm_cptac_2021",
        "profile": "gbm_cptac_2021_protein_quantification",
        "sample_list": "gbm_cptac_2021_protein_quantification",
        "n_expected": 99,
    },
    "PDAC": {
        "study": "paad_cptac_2021",
        "profile": "paad_cptac_2021_protein_quantification",
        "sample_list": "paad_cptac_2021_protein_quantification",
        "n_expected": 81,
    },
    "UCEC": {
        "study": "ucec_cptac_2020",
        "profile": "ucec_cptac_2020_protein_quantification",
        "sample_list": "ucec_cptac_2020_custom_protein",
        "n_expected": 83,
    },
    "LUSC": {
        "study": "lusc_cptac_2021",
        "profile": "lusc_cptac_2021_protein_quantification",
        "sample_list": "lusc_cptac_2021_protein",
        "n_expected": 108,
    },
}

# ── Gene mapping (symbol → Entrez ID) ──
# All key genes from the analysis + paralog pairs
GENE_ENTREZ = {
    # Known SL paralogs
    "ARID1A": 8289, "ARID1B": 57492,
    "PIK3CA": 5290, "PIK3CB": 5291,
    "BRCA1": 672, "BRCA2": 675,
    "EP300": 2033, "CREBBP": 1387,
    "PPP2R1A": 5518, "PPP2R1B": 5519,
    "FBXW7": 55294, "FBXW2": 26190,
    "STK11": 6794, "SIK1": 150094,
    "SMARCA4": 6597, "SMARCA2": 6595,
    "CCNE1": 898, "CCNE2": 9134,
    "CDK4": 1019, "CDK6": 1021,
    "AKT1": 207, "AKT2": 208,
    "MAP2K1": 5604, "MAP2K2": 5605,  # MEK1/MEK2
    
    # De novo candidates
    "KRAS": 3845, "HRAS": 3265, "NRAS": 4893,
    "PIK3R1": 5295, "CRKL": 1399, "CRK": 1398,
    "PTEN": 5728, "TNS2": 23371, "TNS1": 7145,
    "TP53": 7157, "TP63": 8626, "TP73": 7161,
    "RB1": 5925, "RBL1": 5933, "RBL2": 5934,
    "NF1": 4763, "RASA1": 5921, "RASA2": 5922,
    "ATR": 545, "ATM": 472,
    "KMT2D": 8085, "KMT2C": 58508,
    "CDH1": 999, "CDH2": 1000,
    "CTNNB1": 1499,
    "BRAF": 673, "RAF1": 5894,
    "APC": 324, "SMAD4": 4089,
    "ERBB2": 2064, "EGFR": 1956,
    "GATA3": 2625, "KEAP1": 9817,
    "NFE2L2": 4780,  # NRF2
}

# ── Paralog pairs to test ──
PARALOG_PAIRS = [
    # Known SL (gold standard)
    ("ARID1A", "ARID1B"), ("PIK3CA", "PIK3CB"),
    ("BRCA1", "BRCA2"), ("EP300", "CREBBP"),
    ("PPP2R1A", "PPP2R1B"), ("FBXW7", "FBXW2"),
    ("STK11", "SIK1"), ("SMARCA4", "SMARCA2"),
    ("CCNE1", "CCNE2"), ("CDK4", "CDK6"),
    ("AKT1", "AKT2"), ("MAP2K1", "MAP2K2"),
    
    # De novo & exploratory
    ("KRAS", "HRAS"), ("KRAS", "NRAS"), ("HRAS", "NRAS"),
    ("PIK3R1", "CRKL"), ("PIK3R1", "CRK"),
    ("PTEN", "TNS2"), ("PTEN", "TNS1"),
    ("TP53", "TP63"), ("TP53", "TP73"),
    ("RB1", "RBL1"), ("RB1", "RBL2"),
    ("NF1", "RASA2"), ("NF1", "RASA1"),
    ("ATR", "ATM"),
    ("KMT2D", "KMT2C"),
    ("CDH1", "CDH2"),
    ("BRAF", "RAF1"),
    ("KEAP1", "NFE2L2"),
]

BASE_URL = "https://www.cbioportal.org/api"


def fetch_protein_data(profile_id, sample_list_id, gene_symbol, entrez_id, timeout=60):
    """Fetch protein abundance for one gene from cBioPortal API (POST)."""
    url = f"{BASE_URL}/molecular-profiles/{profile_id}/molecular-data/fetch"
    body = {"sampleListId": sample_list_id, "entrezGeneIds": [entrez_id]}
    
    try:
        r = requests.post(url, json=body, timeout=timeout,
                          headers={"Content-Type": "application/json"})
        if not r.ok:
            return None, f"HTTP {r.status_code}"
        
        if not r.text.strip():
            return None, "empty response"
        
        data = r.json()
        values = {}
        for d in data:
            try:
                sample = d.get("sampleId", "")
                val = d.get("value", "")
                if sample and val != "" and val != "NA" and val != "NaN":
                    values[sample] = float(val)
            except (ValueError, KeyError, TypeError):
                pass
        
        if values:
            return pd.Series(values, name=gene_symbol), None
        else:
            return None, "no valid values"
    
    except requests.exceptions.Timeout:
        return None, "timeout"
    except Exception as e:
        return None, str(e)[:100]


def compute_pair_correlations(prot_data, pairs, min_samples=10):
    """
    Compute Pearson correlations for all paralog pairs.
    Returns list of dicts with correlation results.
    """
    results = []
    for a, b in pairs:
        if a not in prot_data or b not in prot_data:
            results.append({
                "gene_a": a, "gene_b": b,
                "n": 0, "r": np.nan, "p": np.nan,
                "status": "missing_data"
            })
            continue
        
        common = prot_data[a].index.intersection(prot_data[b].index)
        n = len(common)
        
        if n < min_samples:
            results.append({
                "gene_a": a, "gene_b": b,
                "n": n, "r": np.nan, "p": np.nan,
                "status": "insufficient_samples"
            })
            continue
        
        r_val, p_val = stats.pearsonr(prot_data[a][common], prot_data[b][common])
        results.append({
            "gene_a": a, "gene_b": b,
            "n": n, "r": r_val, "p": p_val,
            "status": "ok"
        })
    
    return results


def run_cptac_expansion(use_cache=True):
    """Main entry point for CPTAC multi-cancer analysis."""
    cache_dir = DATA_DIR / "cptac_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("  CPTAC Multi-Cancer Proteomics Paralogs Validation")
    print(f"  {len(CPTAC_STUDIES)} CPTAC cohorts × {len(PARALOG_PAIRS)} paralog pairs")
    print("=" * 70)
    
    all_cohort_results = {}
    cohort_summaries = []
    
    for cohort_name, cfg in CPTAC_STUDIES.items():
        print(f"\n{'─' * 70}")
        print(f"  {cohort_name} CPTAC ({cfg['study']})")
        print(f"{'─' * 70}")
        
        cache_file = cache_dir / f"{cohort_name}_protein_data.json"
        
        prot_data = {}
        
        # Try loading from cache
        if use_cache and cache_file.exists():
            try:
                with open(cache_file, "r") as f:
                    raw = json.load(f)
                for gene, vals in raw.items():
                    if vals:
                        prot_data[gene] = pd.Series(vals, name=gene)
                print(f"  Loaded {len(prot_data)} genes from cache")
            except Exception as e:
                print(f"  Cache load failed: {e}")
                cache_file.unlink(missing_ok=True)
        
        # Fetch missing genes from API
        genes_to_fetch = [g for g in GENE_ENTREZ if g not in prot_data]
        
        if genes_to_fetch:
            print(f"  Fetching {len(genes_to_fetch)} genes from API...")
            n_success = 0
            
            for gene in genes_to_fetch:
                eid = GENE_ENTREZ[gene]
                series, error = fetch_protein_data(
                    cfg["profile"], cfg["sample_list"], gene, eid
                )
                
                if series is not None:
                    prot_data[gene] = series
                    n_success += 1
                    # Progress indicator
                    if n_success % 10 == 0:
                        print(f"    {n_success}/{len(genes_to_fetch)}...")
                else:
                    # Don't print every failure to avoid noise
                    pass
                
                # Rate limiting: 300ms between API calls
                time.sleep(0.3)
            
            print(f"  Fetched: {n_success}/{len(genes_to_fetch)} genes")
            
            # Save to cache
            try:
                cache_dict = {}
                for gene, series in prot_data.items():
                    cache_dict[gene] = series.to_dict()
                with open(cache_file, "w") as f:
                    json.dump(cache_dict, f)
            except Exception:
                pass
        
        # Report samples per gene
        sample_counts = {g: len(s) for g, s in prot_data.items()}
        if sample_counts:
            print(f"  Sample range: {min(sample_counts.values())}–{max(sample_counts.values())}")
        
        # Compute correlations
        corr_results = compute_pair_correlations(prot_data, PARALOG_PAIRS)
        all_cohort_results[cohort_name] = corr_results
        
        # Summary stats
        n_sig = sum(1 for r in corr_results if r["status"] == "ok" and r["p"] < 0.05)
        n_tested = sum(1 for r in corr_results if r["status"] == "ok")
        median_r = np.median([r["r"] for r in corr_results if r["status"] == "ok" and not np.isnan(r["r"])])
        
        print(f"  Pairs tested: {n_tested}/{len(PARALOG_PAIRS)}")
        print(f"  Significant (p<0.05): {n_sig}")
        print(f"  Median r: {median_r:.3f}")
        
        # Show top correlations
        sig_results = [r for r in corr_results if r["status"] == "ok" and r["p"] < 0.05]
        sig_results.sort(key=lambda x: abs(x["r"]), reverse=True)
        
        if sig_results:
            print(f"\n  Top correlations:")
            for r in sig_results[:10]:
                star = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*"
                is_known = ((r["gene_a"], r["gene_b"]) in 
                           {(a, b) for a, b in KNOWN_PARALOG_SL} | 
                           {(b, a) for a, b in KNOWN_PARALOG_SL})
                flag = "★" if is_known else "·"
                print(f"    {flag} {r['gene_a']:8s}↔{r['gene_b']:8s}  "
                      f"n={r['n']:>4d}  r={r['r']:+.4f}  p={r['p']:.4f} {star}")
        
        cohort_summaries.append({
            "cohort": cohort_name,
            "study": cfg["study"],
            "n_genes": len(prot_data),
            "n_pairs_tested": n_tested,
            "n_sig": n_sig,
            "median_r": median_r,
        })
    
    # ── Cross-cohort synthesis ──
    print(f"\n{'=' * 70}")
    print(f"  Cross-Cohort Synthesis")
    print(f"{'=' * 70}")
    
    cohort_df = pd.DataFrame(cohort_summaries)
    
    print(f"\n{'Cohort':10s} {'Genes':>6s} {'Pairs':>6s} {'Sig':>6s} {'Median r':>9s}")
    print("-" * 45)
    for _, r in cohort_df.iterrows():
        print(f"{r['cohort']:10s} {int(r['n_genes']):>6d} "
              f"{int(r['n_pairs_tested']):>6d} {int(r['n_sig']):>6d} "
              f"{r['median_r']:>+8.3f}")
    
    # ── Build cross-cohort matrix: for each paralog pair, show r across cohorts ──
    print(f"\n{'─' * 70}")
    print(f"  Paralogs with Consistent Cross-Cohort Co-Variation")
    print(f"{'─' * 70}")
    
    # Collect per-pair data
    pair_across_cohorts = {}
    for cohort_name, results in all_cohort_results.items():
        for r in results:
            if r["status"] != "ok":
                continue
            pair_key = f"{r['gene_a']}↔{r['gene_b']}"
            if pair_key not in pair_across_cohorts:
                pair_across_cohorts[pair_key] = {}
            pair_across_cohorts[pair_key][cohort_name] = r
    
    # Find pairs with significant correlation in ≥2 cohorts
    consistent_pairs = []
    for pair_key, cohort_data in pair_across_cohorts.items():
        sig_cohorts = [c for c in cohort_data if cohort_data[c]["p"] < 0.05]
        if len(sig_cohorts) >= 2:
            mean_r = np.mean([cohort_data[c]["r"] for c in cohort_data])
            consistent_pairs.append({
                "pair": pair_key,
                "n_cohorts": len(sig_cohorts),
                "total_cohorts": len(cohort_data),
                "mean_r": mean_r,
                "sig_cohorts": ", ".join(sig_cohorts),
            })
    
    consistent_pairs.sort(key=lambda x: (x["n_cohorts"], abs(x["mean_r"])), reverse=True)
    
    if consistent_pairs:
        print(f"\n  Pairs with significant correlation in ≥2 cohorts:")
        for p in consistent_pairs[:15]:
            gene_a, gene_b = p["pair"].split("↔")
            is_known = ((gene_a, gene_b) in KNOWN_PARALOG_SL or 
                       (gene_b, gene_a) in KNOWN_PARALOG_SL)
            flag = "★" if is_known else "·"
            print(f"    {flag} {p['pair']:25s}  sig in {p['n_cohorts']}/{p['total_cohorts']} "
                  f"cohorts  mean r={p['mean_r']:+.3f}  [{p['sig_cohorts']}]")
    else:
        print("  No pairs found with consistent cross-cohort significance.")
    
    # Save results
    # Per-cohort details
    for cohort_name, results in all_cohort_results.items():
        df = pd.DataFrame(results)
        df.to_csv(OUTPUT_DIR / f"cptac_{cohort_name.lower()}_correlations.csv", index=False)
    
    # Summary
    cohort_df.to_csv(OUTPUT_DIR / "cptac_cohort_summary.csv", index=False)
    
    # Cross-cohort matrix
    if pair_across_cohorts:
        matrix_rows = []
        cohorts_list = list(CPTAC_STUDIES.keys())
        for pair_key, cohort_data in pair_across_cohorts.items():
            row = {"pair": pair_key}
            for c in cohorts_list:
                if c in cohort_data:
                    row[f"{c}_r"] = cohort_data[c]["r"]
                    row[f"{c}_p"] = cohort_data[c]["p"]
                    row[f"{c}_n"] = cohort_data[c]["n"]
            matrix_rows.append(row)
        matrix_df = pd.DataFrame(matrix_rows)
        matrix_df.to_csv(OUTPUT_DIR / "cptac_pair_matrix.csv", index=False)
    
    print(f"\nResults saved to {OUTPUT_DIR}/cptac_*.csv")
    return all_cohort_results, cohort_df


if __name__ == "__main__":
    run_cptac_expansion()
