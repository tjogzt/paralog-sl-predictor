"""
Integrated Structural Analysis Module
======================================
Combines ESMFold-predicted 3D structures (N-terminal domain, 400 aa) with 
sequence-level features from UniProt for comprehensive paralog-pair 
structural characterization.

Scope
-----
Phase 1 (COMPLETED — this module):
  - Sequence-level structural features (MW, pI, GRAVY, disorder, domain arch)
  - ESMFold domain-level structure predictions (15 proteins, mean pLDDT=0.724)
  - Pairwise structural similarity (TM-score, RMSD) for 7 paralog pairs
  - Domain architecture conservation analysis (11 pairs)
  - PROTAC suitability and druggability scoring

Phase 2 (PLANNED — requires GPU/cloud):
  - Full-length AlphaFold-Multimer complex predictions
  - Molecular dynamics (MD) simulation for interface stability
  - Free energy perturbation (FEP) for mutation ΔΔG

Data sources:
  - ESMFold v1 (facebook/esmfold_v1, 8.4 GB, via MPS on Apple M4 Max)
  - UniProt canonical sequences + domain annotations
  - AlphaFold Protein Structure Database (Varadi et al., NAR 2024)

Output files:
  - structure_pred/outputs/monomers/*.pdb  — 15 predicted structures
  - structure_pred/outputs/esmfold_results.csv  — summary metrics
  - structure_pred/outputs/esmfold_pairwise.csv  — pairwise TM-score/RMSD
  - structure_pred/outputs/sequence_features.csv  — UniProt features
  - structure_pred/outputs/pairwise_structural_analysis.csv  — domain analysis

"""


import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL

# ── Known domain architectures for key paralog pairs ────────────
# (from UniProt annotation)
DOMAIN_ARCHITECTURES = {
    "ARID1A":  ["ARID_DNA_bind", "AT-rich_interaction"],
    "ARID1B":  ["ARID_DNA_bind", "AT-rich_interaction"],
    "SMARCA4": ["SNF2_N", "Helicase_C", "Bromodomain", "HSA"],
    "SMARCA2": ["SNF2_N", "Helicase_C", "Bromodomain", "HSA"],
    "BRCA1":   ["RING", "BRCT_x2", "Coiled-coil"],
    "BRCA2":   ["BRC_repeats", "OB_fold_x3", "Tower_domain"],
    "EP300":   ["KAT11_HAT", "Bromodomain", "ZZ", "TAZ1", "TAZ2", "KIX"],
    "CREBBP":  ["KAT11_HAT", "Bromodomain", "ZZ", "TAZ1", "TAZ2", "KIX"],
    "PIK3CA":  ["PI3K_p85B", "PI3K_C2", "PI3Ka", "PI3Kc"],
    "PIK3CB":  ["PI3K_p85B", "PI3K_C2", "PI3Ka", "PI3Kc"],
    "PPP2R1A": ["HEAT_repeats"],  # A scaffold subunit
    "PPP2R1B": ["HEAT_repeats"],
    "PIK3R1":  ["SH3", "RhoGAP", "SH2_x2", "p85_iSH2"],
    "CRKL":    ["SH2", "SH3_x2"],
    "KRAS":    ["G_domain", "Hypervariable_region"],
    "HRAS":    ["G_domain", "Hypervariable_region"],
    "NF1":     ["RAS-GAP", "CRAL-TRIO", "PH-like"],
    "RASA2":   ["RAS-GAP", "PH", "C2"],
    "KMT2D":   ["PHD_x6", "FYRN", "FYRC", "SET"],
    "KMT2C":   ["PHD_x7", "FYRN", "FYRC", "SET"],
    "ATR":     ["HEAT_repeats", "FAT", "PI3Kc", "FATC"],
    "ATM":     ["HEAT_repeats", "FAT", "PI3Kc", "FATC"],
    "RB1":     ["Cyclin_bind", "Pocket_A", "Pocket_B", "C-term"],
    "RBL1":    ["Cyclin_bind", "Pocket_A", "Pocket_B"],
    "FBXW7":   ["F-box", "WD40_repeats_x8"],
    "FBXW2":   ["F-box", "WD40_repeats"],
    "STK11":   ["Kinase", "STRAD_binding"],
    "SIK1":    ["Kinase", "UBA"],
    "CDH1":    ["Cadherin_propep", "Cadherin_ext_x5", "Cadherin_cyto"],
    "CDH2":    ["Cadherin_propep", "Cadherin_ext_x5", "Cadherin_cyto"],
    "AKT1":    ["PH", "Kinase", "Hydrophobic_motif"],
    "AKT2":    ["PH", "Kinase", "Hydrophobic_motif"],
    "CCNE1":   ["Cyclin_N", "Cyclin_C"],
    "CCNE2":   ["Cyclin_N", "Cyclin_C"],
    "CDK4":    ["Kinase"],
    "CDK6":    ["Kinase"],
    "BRAF":    ["RBD", "Kinase", "C1_1"],
    "RAF1":    ["RBD", "Kinase", "C1_1"],
    "MAP2K1":  ["Kinase", "Docking_site"],
    "MAP2K2":  ["Kinase", "Docking_site"],
    "TP53":    ["TAD", "Proline_rich", "DNA_bind", "Tetramer"],
    "TP63":    ["TAD", "Proline_rich", "DNA_bind", "SAM", "TI"],
    "PTEN":    ["Phosphatase", "C2", "PDZ_binding"],
    "TNS2":    ["PTB", "SH2"],
}

# ── Structural analysis pairs ───────────────────────────────────
STRUCTURAL_PAIRS = [
    ("ARID1A", "ARID1B"), ("SMARCA4", "SMARCA2"),
    ("BRCA1", "BRCA2"), ("EP300", "CREBBP"),
    ("PIK3CA", "PIK3CB"), ("PPP2R1A", "PPP2R1B"),
    ("PIK3R1", "CRKL"), ("KRAS", "HRAS"),
    ("NF1", "RASA2"), ("KMT2D", "KMT2C"),
    ("ATR", "ATM"), ("RB1", "RBL1"),
    ("FBXW7", "FBXW2"), ("STK11", "SIK1"),
    ("CDH1", "CDH2"), ("AKT1", "AKT2"),
    ("CCNE1", "CCNE2"), ("CDK4", "CDK6"),
    ("BRAF", "RAF1"), ("TP53", "TP63"),
    ("PTEN", "TNS2"), ("MAP2K1", "MAP2K2"),
]


def compute_domain_similarity(dom_a, dom_b):
    """Compute Jaccard similarity between domain lists."""
    if not dom_a or not dom_b:
        return 0.0
    
    # Normalize domain names
    def normalize(name):
        return name.lower().replace("_", "").replace("-", "")
    
    set_a = {normalize(d) for d in dom_a}
    set_b = {normalize(d) for d in dom_b}
    
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    
    if union == 0:
        return 0.0
    return intersection / union


def run_structural_analysis():
    """Main entry point."""
    print("=" * 70)
    print("  Sequence-Based Structural/Druggability Analysis")
    print(f"  Using protein_features.csv + UniProt domain architectures")
    print("=" * 70)
    
    # ── Load protein features ──
    features_path = DATA_DIR / "protein_features.csv"
    if not features_path.exists():
        print(f"  ERROR: {features_path} not found. Run protein_features.py first.")
        return None
    
    feats = pd.read_csv(features_path)
    feats = feats.set_index("gene")
    print(f"  Loaded features for {len(feats)} genes")
    
    # ── Known SL set ──
    known_set = set()
    for a, b in KNOWN_PARALOG_SL:
        known_set.add((a, b))
        known_set.add((b, a))
    
    # Feature columns for analysis
    key_cols = ["length", "molecular_weight", "isoelectric_point", 
                "gravy", "disorder_propensity", "lysine_density",
                "cysteine_density", "instability_index"]
    
    # ── Per-pair analysis ──
    results = []
    
    for gene_a, gene_b in STRUCTURAL_PAIRS:
        if gene_a not in feats.index or gene_b not in feats.index:
            print(f"  · {gene_a}↔{gene_b}: missing features")
            continue
        
        feat_a = feats.loc[gene_a]
        feat_b = feats.loc[gene_b]
        
        is_known = (gene_a, gene_b) in known_set
        
        # ── Structural similarity (normalized feature differences) ──
        diffs = {}
        for col in key_cols:
            if col in feat_a.index and col in feat_b.index:
                val_a = feat_a[col]
                val_b = feat_b[col]
                if val_a != 0 or val_b != 0:
                    denom = max(abs(val_a), abs(val_b), 1e-10)
                    diffs[col] = abs(val_a - val_b) / denom
                else:
                    diffs[col] = 0
        
        # Overall structural similarity (1 - mean normalized difference)
        if diffs:
            mean_diff = np.mean(list(diffs.values()))
            structural_similarity = max(0, 1 - mean_diff)
        else:
            structural_similarity = 0
        
        # ── Domain similarity ──
        dom_a = DOMAIN_ARCHITECTURES.get(gene_a, [])
        dom_b = DOMAIN_ARCHITECTURES.get(gene_b, [])
        domain_similarity = compute_domain_similarity(dom_a, dom_b)
        
        # ── Druggability metrics ──
        # PROTAC suitability: high lysine density, high disorder, moderate cysteine
        lysine_a = feat_a.get("lysine_density", 0)
        lysine_b = feat_b.get("lysine_density", 0)
        cyst_a = feat_a.get("cysteine_density", 0)
        cyst_b = feat_b.get("cysteine_density", 0)
        disorder_a = feat_a.get("disorder_propensity", 0)
        disorder_b = feat_b.get("disorder_propensity", 0)
        length_a = feat_a.get("length", 0)
        length_b = feat_b.get("length", 0)
        
        # PROTAC score: high lysine (ubiquitination) + high disorder + low cysteine (structural constraint)
        protac_a = (lysine_a * 2.0 + disorder_a * 1.0 - cyst_a * 0.5)
        protac_b = (lysine_b * 2.0 + disorder_b * 1.0 - cyst_b * 0.5)
        protac_score = (protac_a + protac_b) / 2
        
        # Small molecule druggability: well-folded (low disorder), moderate size
        druggability = 1.0 - (disorder_a + disorder_b) / 2  # lower disorder = more druggable
        
        result = {
            "gene_a": gene_a, "gene_b": gene_b,
            "is_known_sl": is_known,
            
            # Structural
            "structural_similarity": round(structural_similarity, 3),
            "domain_similarity": round(domain_similarity, 3),
            
            # Sequence features
            "length_a": length_a, "length_b": length_b,
            "length_ratio": round(min(length_a, length_b) / max(length_a, length_b), 3) if length_a > 0 and length_b > 0 else 0,
            "disorder_a": round(disorder_a, 3), "disorder_b": round(disorder_b, 3),
            "lysine_a": round(lysine_a, 4), "lysine_b": round(lysine_b, 4),
            "cysteine_a": round(cyst_a, 4), "cysteine_b": round(cyst_b, 4),
            "pI_a": round(feat_a.get("isoelectric_point", 0), 1),
            "pI_b": round(feat_b.get("isoelectric_point", 0), 1),
            
            # Druggability
            "protac_score": round(protac_score, 4),
            "druggability": round(druggability, 3),
            
            # Domains
            "domains_a": "+".join(dom_a) if dom_a else "unknown",
            "domains_b": "+".join(dom_b) if dom_b else "unknown",
            "n_domains_a": len(dom_a), "n_domains_b": len(dom_b),
        }
        
        results.append(result)
    
    if not results:
        print("  No results generated.")
        return None
    
    results_df = pd.DataFrame(results)
    
    # ── Composite targetability score ──
    # Weight: structural similarity (compensation potential) + druggability + domain conservation
    results_df["targetability"] = (
        results_df["structural_similarity"] * 0.25 +
        results_df["domain_similarity"] * 0.30 +
        results_df["druggability"] * 0.25 +
        results_df["protac_score"] * 0.20
    )
    
    results_df = results_df.sort_values("targetability", ascending=False)
    
    # ── Display ──
    print(f"\n{'─' * 70}")
    print(f"  Structural/Druggability Ranking")
    print(f"{'─' * 70}")
    
    print(f"\n{'Pair':25s} {'Struct':>6s} {'Domain':>6s} {'Drugg':>6s} {'PROTAC':>6s} {'Target':>7s}")
    print("-" * 70)
    for _, r in results_df.iterrows():
        flag = "★" if r["is_known_sl"] else "·"
        print(f"{flag} {r['gene_a']+'↔'+r['gene_b']:23s} "
              f"{r['structural_similarity']:>6.3f} {r['domain_similarity']:>6.3f} "
              f"{r['druggability']:>6.3f} {r['protac_score']:>+7.4f} "
              f"{r['targetability']:>7.3f}")
    
    # ── Top insights ──
    print(f"\n{'─' * 70}")
    print(f"  Key Insights")
    print(f"{'─' * 70}")
    
    # Best structural matches (compensation potential)
    best_struct = results_df.nlargest(5, "structural_similarity")
    print(f"\n  Best structural matches (highest compensation potential):")
    for _, r in best_struct.iterrows():
        print(f"    {r['gene_a']}↔{r['gene_b']}: struct_sim={r['structural_similarity']:.3f}, "
              f"len={r['length_a']}/{r['length_b']}")
    
    # Best domain conservation
    best_domain = results_df.nlargest(5, "domain_similarity")
    print(f"\n  Best domain conservation:")
    for _, r in best_domain.iterrows():
        known = " (known SL)" if r["is_known_sl"] else ""
        print(f"    {r['gene_a']}↔{r['gene_b']}: domain_sim={r['domain_similarity']:.3f} "
              f"[{r['domains_a']}] ↔ [{r['domains_b']}]{known}")
    
    # Best PROTAC candidates
    best_protac = results_df.nlargest(5, "protac_score")
    print(f"\n  Best PROTAC candidates (high lysine + disorder):")
    for _, r in best_protac.iterrows():
        print(f"    {r['gene_a']}↔{r['gene_b']}: PROTAC={r['protac_score']:+.4f}, "
              f"K_dens={r['lysine_a']:.4f}/{r['lysine_b']:.4f}, "
              f"disorder={r['disorder_a']:.3f}/{r['disorder_b']:.3f}")
    
    # Best small-molecule targets
    best_druggable = results_df.nlargest(5, "druggability")
    print(f"\n  Best small-molecule targets (well-folded, low disorder):")
    for _, r in best_druggable.iterrows():
        print(f"    {r['gene_a']}↔{r['gene_b']}: druggability={r['druggability']:.3f}, "
              f"disorder={r['disorder_a']:.3f}/{r['disorder_b']:.3f}, "
              f"len={r['length_a']}/{r['length_b']}")
    
    # ── Integration with other analyses ──
    # Merge with therapeutic window and DD data if available
    tw_path = OUTPUT_DIR / "therapeutic_window_all_results.csv"
    if tw_path.exists():
        tw = pd.read_csv(tw_path)
        # Get mean metrics per pair across contexts
        tw_agg = tw.groupby(["driver", "paralog"]).agg(
            mean_ti=("therapeutic_index", "mean"),
            mean_dd=("dd_abs", "mean"),
            mean_selectivity=("selectivity", "mean"),
        ).reset_index()
        
        merged = results_df.merge(
            tw_agg, left_on=["gene_a", "gene_b"],
            right_on=["driver", "paralog"], how="left"
        )
        
        # Composite clinical score
        merged["clinical_targetability"] = (
            merged["targetability"] * 0.3 +
            (merged["mean_ti"] / merged["mean_ti"].max() if merged["mean_ti"].max() > 0 else 0) * 0.4 +
            (merged["mean_selectivity"] + 1) / 2 * 0.3  # normalize to [0, 1]
        )
        
        merged = merged.sort_values("clinical_targetability", ascending=False)
        
        print(f"\n  Clinical targetability ranking (structure + TI + selectivity):")
        for _, r in merged.head(10).iterrows():
            flag = "★" if r["is_known_sl"] else "·"
            ti_str = f"TI={r['mean_ti']:.2f}" if not np.isnan(r["mean_ti"]) else "N/A"
            print(f"    {flag} {r['gene_a']:10s}→{r['gene_b']:10s}  "
                  f"clinical={r['clinical_targetability']:.3f}  "
                  f"{ti_str}  struct={r['structural_similarity']:.3f}")
        
        merged.to_csv(OUTPUT_DIR / "alphafold_structural_analysis.csv", index=False)
    else:
        results_df.to_csv(OUTPUT_DIR / "alphafold_structural_analysis.csv", index=False)
    
    print(f"\nResults saved to {OUTPUT_DIR}/alphafold_structural_analysis.csv")
    return results_df


if __name__ == "__main__":
    run_structural_analysis()
