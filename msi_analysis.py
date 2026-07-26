"""
MSI/dMMR Subgroup Analysis for Paralog-SL
==========================================
Classifies UCEC and CRC cell lines by microsatellite instability status
using the official DepMap 26Q1 annotation (OmicsGlobalSignatures.csv,
MSIScore computed with MSIsensor2 by the DepMap consortium; lines with
MSIScore > 20 are labelled MSI-H, otherwise MSS), and evaluates
paralog-SL signal separately for MSI-H vs MSS subgroups.

The earlier mutation-based proxy (damaging MLH1/MSH2/MSH6/PMS2/POLE
mutations) is retained ONLY as a sensitivity cross-check in the audit
outputs (output/msi_classification_audit.csv). POLE-mutant lines are
reported separately there: POLE-ultramutated tumours are typically
microsatellite-stable, so the official MSIsensor2 annotation — unlike
the proxy — does not pool them with dMMR/MSI-H lines.

Rationale: Endometrial (UCEC ~30%) and colorectal (CRC ~15%) cancers
have substantial MSI-H populations with distinct mutational landscapes.
MMR status profoundly affects dependency relationships and may
modulate paralog compensation signals.
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score
from pathlib import Path

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL
from data_loader import (
    load_dependency, load_expression, load_models, load_mutations, load_paralogs,
    build_mutation_matrix,
)
from pcs import ParalogCompensationScore

# ── MMR genes and related ──────────────────────────────────────
MMR_GENES = ["MLH1", "MSH2", "MSH6", "PMS2"]
POLE_GENE = "POLE"   # POLE exonuclease domain mutations also cause hypermutation

# ── Cancer types with MSI relevance ────────────────────────────
MSI_CANCERS = {
    "Endometrial": [
        "Endometrial Carcinoma", "Endometrial Cancer",
        "Uterine Serous Carcinoma", "Uterine Carcinosarcoma",
        "Endometrial Endometrioid Adenocarcinoma",
    ],
    "Colorectal": [
        "Colorectal Adenocarcinoma",
    ],
}

# ── Drivers for each cancer type ───────────────────────────────
MSI_DRIVERS = {
    "Endometrial": ["PTEN", "ARID1A", "PIK3CA", "CTNNB1", "TP53",
                     "PPP2R1A", "KRAS", "PIK3R1", "FBXW7", "KMT2D",
                     "BRCA1", "BRCA2", "RB1", "NF1", "ATR", "ATM"],
    "Colorectal": ["TP53", "APC", "KRAS", "PIK3CA", "BRAF", "SMAD4",
                    "FBXW7", "PTEN", "ARID1A", "CTNNB1", "BRCA1", "BRCA2",
                    "RB1", "NF1", "ATM", "ATR"],
}


# ── Official MSI annotation (DepMap 26Q1) ──────────────────────
SIGNATURES_FILE = DATA_DIR / "OmicsGlobalSignatures.csv"
MSI_SCORE_THRESHOLD = 20.0   # MSIsensor2 MSIscore > 20 → MSI (DepMap convention)


def load_official_msi(path: Path = SIGNATURES_FILE) -> pd.DataFrame:
    """
    Load the official DepMap 26Q1 MSI annotation.

    OmicsGlobalSignatures.csv is profile-level; we keep the default
    profile per model and binarise the MSIsensor2 MSIscore at > 20
    (DepMap consortium convention: MSIscore > 20 = MSI, else MSS).

    Returns DataFrame with columns: DepMap_ID, msi_score, msi_official
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Official MSI annotation not found: {path}. "
            "Download OmicsGlobalSignatures.csv (DepMap Public 26Q1) into data/."
        )
    sig = pd.read_csv(path)
    sig = sig[sig["IsDefaultEntryForModel"] == "Yes"].copy()
    sig["msi_official"] = np.where(
        sig["MSIScore"] > MSI_SCORE_THRESHOLD, "MSI_H", "MSS"
    )
    return sig[["ModelID", "MSIScore", "msi_official"]].rename(
        columns={"ModelID": "DepMap_ID", "MSIScore": "msi_score"}
    )


def classify_official_msi(models_df, sig_df, mutations_df):
    """
    Classify cell lines by the official MSIsensor2 annotation.
    The mutation proxy is used only for lines missing from the
    signatures file (flagged source='proxy_fallback').

    Returns DataFrame: DepMap_ID, msi_status, msi_score, msi_source
    """
    proxy = classify_mmr_status(models_df, mutations_df)
    merged = models_df[["DepMap_ID"]].merge(sig_df, on="DepMap_ID", how="left")
    merged = merged.merge(proxy, on="DepMap_ID", how="left")

    has_official = merged["msi_official"].notna()
    merged["msi_status"] = np.where(
        has_official, merged["msi_official"], merged["mmr_status"]
    )
    merged["msi_source"] = np.where(
        has_official, "MSIsensor2_26Q1", "proxy_fallback"
    )
    return merged[["DepMap_ID", "msi_status", "msi_score", "msi_source"]]


def build_classification_audit(models_df, mutations_df, sig_df):
    """
    Per-line audit table: official MSI call, MMR-mutation proxy call,
    POLE-mutation flag, and their concordance. POLE-mutant lines are
    expected to be microsatellite-stable under the official annotation.
    """
    proxy = classify_mmr_status(models_df, mutations_df)
    pole_ids = set(
        mutations_df.loc[mutations_df["Gene"] == POLE_GENE, "DepMap_ID"]
    )
    audit = models_df[["DepMap_ID"]].merge(sig_df, on="DepMap_ID", how="left")
    audit = audit.merge(proxy, on="DepMap_ID", how="left")
    audit["pole_mutant"] = audit["DepMap_ID"].isin(pole_ids)
    audit["concordant"] = audit["msi_official"] == audit["mmr_status"]
    return audit


def classify_mmr_status(models_df, mutations_df):
    """
    Mutation-based proxy for MMR deficiency (damaging mutations in
    MLH1/MSH2/MSH6/PMS2 or POLE). Retained for sensitivity audit only —
    not used for the primary subgroup classification.
    
    Returns DataFrame with columns: DepMap_ID, mmr_status
    """
    all_cell_lines = models_df["DepMap_ID"].tolist()
    
    # Get damaging mutations in MMR genes
    mut_sub = mutations_df[
        mutations_df["Gene"].isin(MMR_GENES + [POLE_GENE])
    ]
    
    if mut_sub.empty:
        print("  WARNING: No MMR gene mutations found")
        result = pd.DataFrame({
            "DepMap_ID": all_cell_lines,
            "mmr_status": "MSS"
        })
        return result
    
    # Cell lines with damaging MMR/POLE mutations
    mmr_mutant_ids = mut_sub["DepMap_ID"].unique().tolist()
    
    classifications = []
    for cl in all_cell_lines:
        status = "MSI_H" if cl in mmr_mutant_ids else "MSS"
        classifications.append({"DepMap_ID": cl, "mmr_status": status})
    
    return pd.DataFrame(classifications)


def compute_mutation_type_distribution(mutations_df, cell_ids, driver_genes):
    """
    For diagnostics: report what types of mutations are present
    in driver genes of MSI-H vs MSS cell lines.
    """
    mut_sub = mutations_df[
        mutations_df["DepMap_ID"].isin(cell_ids) &
        mutations_df["Gene"].isin(driver_genes)
    ]
    
    if "VariantType" in mut_sub.columns:
        dist = mut_sub["VariantType"].value_counts().to_dict()
    elif "VariantInfo" in mut_sub.columns:
        dist = mut_sub["VariantInfo"].value_counts().to_dict()
    else:
        dist = {}
    
    n_cells_with_mut = mut_sub["DepMap_ID"].nunique()
    return n_cells_with_mut, dist


def run_msi_analysis():
    """Main entry point."""
    print("=" * 65)
    print("  MSI/dMMR Subgroup Paralog-SL Analysis")
    print("  Classification: official DepMap 26Q1 MSIsensor2 MSIscore")
    print(f"  (OmicsGlobalSignatures.csv; MSI-H = score > {MSI_SCORE_THRESHOLD:g})")
    print("=" * 65)
    
    # ── Load data ──
    dep = load_dependency()
    expr = load_expression()
    models = load_models()
    mutations = load_mutations()
    paralogs = load_paralogs()
    sig_df = load_official_msi()
    print(f"  Official MSI annotation loaded for {len(sig_df)} default-entry models")
    
    # ── Known SL pairs set ──
    known_set = set()
    for a, b in KNOWN_PARALOG_SL:
        known_set.add((a.upper(), b.upper()))
        known_set.add((b.upper(), a.upper()))
    
    all_results = {}
    summary_rows = []
    audit_frames = []
    
    for cancer_name, disease_patterns in MSI_CANCERS.items():
        print(f"\n{'─' * 65}")
        print(f"  {cancer_name}")
        print(f"{'─' * 65}")
        
        # Filter cell lines for this cancer type
        pat = "|".join(disease_patterns)
        mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
        cancer_models = models[mask].copy()
        
        # Get valid cell lines (with dependency + expression data)
        cell_ids = cancer_models["DepMap_ID"].tolist()
        cell_ids = [c for c in cell_ids if c in dep.index and c in expr.index]
        
        if len(cell_ids) < 10:
            print(f"  Insufficient cell lines: {len(cell_ids)}")
            continue
        
        # ── Classify MSI status (official MSIsensor2 annotation) ──
        msi_class = classify_official_msi(cancer_models, sig_df, mutations)
        n_fallback = int((msi_class["msi_source"] == "proxy_fallback").sum())
        msi_ids = msi_class.loc[msi_class["msi_status"] == "MSI_H", "DepMap_ID"].tolist()
        mss_ids = msi_class.loc[msi_class["msi_status"] == "MSS", "DepMap_ID"].tolist()
        
        msi_valid = [c for c in msi_ids if c in dep.index and c in expr.index]
        mss_valid = [c for c in mss_ids if c in dep.index and c in expr.index]
        
        print(f"  Total lines: {len(cell_ids)}")
        print(f"  MSI-H (MSIsensor2 score > {MSI_SCORE_THRESHOLD:g}): {len(msi_valid)} cells")
        print(f"  MSS (score <= {MSI_SCORE_THRESHOLD:g}):              {len(mss_valid)} cells")
        if n_fallback:
            print(f"  NOTE: {n_fallback} lines lacked official annotation and used the mutation proxy fallback")
        
        # ── Classification audit (official vs proxy, POLE separate) ──
        audit = build_classification_audit(cancer_models, mutations, sig_df)
        audit.insert(0, "cancer", cancer_name)
        audit_frames.append(audit)
        both = audit.dropna(subset=["msi_official"])
        if len(both):
            conc = both["concordant"].mean()
            n_pole = int(both["pole_mutant"].sum())
            n_pole_msi = int((both["pole_mutant"] & (both["msi_official"] == "MSI_H")).sum())
            proxy_msi = both["mmr_status"] == "MSI_H"
            flipped = int((proxy_msi & (both["msi_official"] == "MSS")).sum())
            print(f"  Audit: proxy↔official concordance = {conc:.1%} "
                  f"({flipped} proxy-MSI lines are officially MSS); "
                  f"POLE-mutant lines: {n_pole}, of which officially MSI-H: {n_pole_msi}")
        
        # ── Mutation profile diagnostics ──
        drivers = MSI_DRIVERS.get(cancer_name, [])
        msi_n_muts, msi_vtypes = compute_mutation_type_distribution(
            mutations, msi_valid, drivers
        )
        mss_n_muts, mss_vtypes = compute_mutation_type_distribution(
            mutations, mss_valid, drivers
        )
        print(f"  Driver mutations — MSI-H: {msi_n_muts} cells, MSS: {mss_n_muts} cells")
        
        # ── Run analysis for each subgroup ──
        for subgroup_name, subgroup_ids in [("MSI_H", msi_valid), ("MSS", mss_valid)]:
            if len(subgroup_ids) < 6:
                print(f"  {subgroup_name}: insufficient cell lines ({len(subgroup_ids)}) — skipping")
                continue
            
            pcs = ParalogCompensationScore(dep, expr, models, mutations, paralogs)
            results_list = []
            
            for driver in drivers:
                if driver not in dep.columns:
                    continue
                try:
                    result = pcs.compute_pcs_for_driver(
                        driver, subgroup_ids,
                        cancer_label=f"{cancer_name}_{subgroup_name}"
                    )
                    if result is not None and not result.empty:
                        results_list.append(result)
                except Exception as e:
                    pass
            
            if not results_list:
                print(f"  {subgroup_name}: no drivers with ≥3 MUT lines — skipping")
                continue
            
            results_df = pd.concat(results_list, ignore_index=True)
            results_df["is_known_paralog_sl"] = results_df.apply(
                lambda r: (r["driver_gene"].upper(), r["paralog_gene"].upper()) in known_set,
                axis=1
            )
            
            key = f"{cancer_name}_{subgroup_name}"
            all_results[key] = results_df
            
            yt = results_df["is_known_paralog_sl"].astype(int).values
            ys = results_df["dependency_dd"].abs().fillna(0).values
            nk = int(yt.sum())
            n_pairs = len(results_df)
            
            auc = roc_auc_score(yt, ys) if nk >= 2 else float("nan")
            
            summary_rows.append({
                "cancer": cancer_name, "subgroup": subgroup_name,
                "n_lines": len(subgroup_ids), "n_pairs": n_pairs,
                "n_known": nk, "dd_auroc": auc
            })
            
            astr = f"{auc:.3f}" if not np.isnan(auc) else "N/A"
            print(f"  {subgroup_name:6s}: {n_pairs:>4d} pairs, {nk} known, DD AUROC={astr}")
            
            # Save subgroup results
            out_path = OUTPUT_DIR / f"msi_{key.lower()}_results.csv"
            results_df.to_csv(out_path, index=False)
    
    # ── Save classification audit (official vs proxy, POLE separate) ──
    if audit_frames:
        audit_all = pd.concat(audit_frames, ignore_index=True)
        audit_all.to_csv(OUTPUT_DIR / "msi_classification_audit.csv", index=False)
        print(f"  Classification audit saved ({len(audit_all)} lines)")

    # ── Summary ──
    print(f"\n{'=' * 65}")
    print(f"  MSI Subgroup Analysis Summary")
    print(f"{'=' * 65}")
    print(f"{'Cancer':20s} {'Subgroup':8s} {'Lines':>6s} {'Pairs':>6s} {'Known':>6s} {'AUROC':>8s}")
    print(f"{'-' * 60}")
    
    if summary_rows:
        summary = pd.DataFrame(summary_rows)
        summary.to_csv(OUTPUT_DIR / "msi_subgroup_summary.csv", index=False)
        
        for _, r in summary.iterrows():
            astr = f"{r['dd_auroc']:.3f}" if not np.isnan(r['dd_auroc']) else "N/A"
            print(f"{r['cancer']:20s} {r['subgroup']:8s} "
                  f"{int(r['n_lines']):>6d} {int(r['n_pairs']):>6d} "
                  f"{int(r['n_known']):>6d} {astr:>8s}")
        
        # Machine-readable key numbers for manuscript traceability
        import json
        key_numbers = {
            "classification": {
                "source": "DepMap Public 26Q1 OmicsGlobalSignatures.csv",
                "method": "MSIsensor2 MSIscore",
                "threshold": f"MSI-H = MSIscore > {MSI_SCORE_THRESHOLD:g}",
            },
            "subgroups": {
                f"{r['cancer']}_{r['subgroup']}": {
                    "n_lines": int(r["n_lines"]),
                    "n_pairs": int(r["n_pairs"]),
                    "n_known": int(r["n_known"]),
                    "dd_auroc": (None if np.isnan(r["dd_auroc"]) else round(float(r["dd_auroc"]), 4)),
                }
                for _, r in summary.iterrows()
            },
        }
        with open(OUTPUT_DIR / "msi_key_numbers.json", "w") as fh:
            json.dump(key_numbers, fh, indent=2)
        
        # ── Within-cancer comparison ──
        print(f"\n{'─' * 65}")
        print(f"  Within-Cancer MSI-H vs MSS Comparison")
        print(f"{'─' * 65}")
        for cancer_name in MSI_CANCERS:
            msi_row = summary[(summary["cancer"] == cancer_name) & (summary["subgroup"] == "MSI_H")]
            mss_row = summary[(summary["cancer"] == cancer_name) & (summary["subgroup"] == "MSS")]
            if not msi_row.empty and not mss_row.empty:
                msi_auc = msi_row["dd_auroc"].values[0]
                mss_auc = mss_row["dd_auroc"].values[0]
                if not np.isnan(msi_auc) and not np.isnan(mss_auc):
                    delta = msi_auc - mss_auc
                    direction = "MSI-H > MSS" if delta > 0 else "MSS > MSI-H"
                    print(f"  {cancer_name:20s}: MSI_H={msi_auc:.3f}  MSS={mss_auc:.3f}  "
                          f"Δ={delta:+.3f}  ({direction})")
    else:
        print("  No results generated.")
    
    print(f"\nResults saved to {OUTPUT_DIR}/msi_*.csv")
    return summary, all_results


if __name__ == "__main__":
    run_msi_analysis()
