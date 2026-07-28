"""
Fine-grained pan-cancer analysis: specific solid tumor types (≥10 cell lines).
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL
from data_loader import load_dependency, load_expression, load_models, load_mutations, load_paralogs
from pcs import ParalogCompensationScore

# ── Fine-grained solid tumor types (≥10 cell lines with data) ──
SOLID_TUMORS = {
    # Thoracic
    "NSCLC":               ["Non-Small Cell Lung Cancer"],
    "SCLC":                ["Lung Neuroendocrine Tumor", "Small Cell Lung Cancer"],
    "Mesothelioma":        ["Pleural Mesothelioma"],

    # GI tract
    "Colorectal":          ["Colorectal Adenocarcinoma"],
    "Esophagogastric":     ["Esophagogastric Adenocarcinoma", "Esophageal Squamous Cell Carcinoma"],
    "Pancreatic":          ["Pancreatic Adenocarcinoma"],
    "Hepatocellular":      ["Hepatocellular Carcinoma", "Hepatoblastoma"],
    "Biliary Tract":       ["Intraductal Papillary Neoplasm of the Bile Duct",
                             "Intracholecystic Papillary Neoplasm", "Ampullary Carcinoma"],

    # Breast
    "Breast":              ["Invasive Breast Carcinoma", "Breast Ductal Carcinoma In Situ",
                             "Breast Neoplasm, NOS"],

    # GYN
    "Ovarian":             ["Ovarian Epithelial Tumor"],
    "Endometrial":         ["Endometrial Carcinoma", "Uterine Sarcoma/Mesenchymal"],
    "Cervical":            ["Cervical Squamous Cell Carcinoma", "Cervical Adenocarcinoma"],

    # GU
    "Renal Cell":          ["Renal Cell Carcinoma", "Wilms' Tumor"],
    "Bladder Urothelial":  ["Bladder Urothelial Carcinoma"],
    "Prostate":            ["Prostate Adenocarcinoma"],

    # Skin
    "Melanoma":            ["Melanoma", "Ocular Melanoma", "Merkel Cell Carcinoma"],

    # CNS
    "Glioma":              ["Adult-Type Diffuse Glioma", "Diffuse Glioma", "Meningothelial Tumor"],
    "Neuroblastoma":       ["Neuroblastoma"],

    # Sarcoma
    "Osteosarcoma":        ["Osteosarcoma"],
    "Ewing Sarcoma":       ["Ewing Sarcoma"],
    "Rhabdomyosarcoma":    ["Rhabdomyosarcoma"],
    "Liposarcoma":         ["Liposarcoma"],
    "Leiomyosarcoma":      ["Leiomyosarcoma"],
    "Other Sarcoma":       ["Undifferentiated Pleomorphic Sarcoma/Malignant Fibrous Histiocytoma/"
                            "High-Grade Spindle Cell Sarcoma", "Synovial Sarcoma",
                            "Rhabdoid Cancer", "SMARCA4-deficient undifferentiated tumor",
                            "Nerve Sheath Tumor", "Chondrosarcoma", "Chordoma",
                            "Fibrosarcoma", "Epithelioid Sarcoma"],

    # Head and Neck
    "HNSCC":               ["Head and Neck Squamous Cell Carcinoma", "Head and Neck Carcinoma, Other"],
    "Thyroid":             ["Anaplastic Thyroid Cancer", "Well-Differentiated Thyroid Cancer",
                             "Medullary Thyroid Cancer"],

    # Other
    "Testicular":          ["Non-Seminomatous Germ Cell Tumor"],
}

SOLID_DRIVERS = [
    "TP53", "PIK3CA", "PTEN", "ARID1A", "KRAS", "BRAF", "CTNNB1",
    "BRCA1", "BRCA2", "RB1", "NF1", "CDKN2A", "APC", "KMT2D",
    "ATM", "ATR", "SMARCA4", "FBXW7", "EP300", "PPP2R1A",
    "EGFR", "ERBB2", "MET", "ALK", "STK11", "KEAP1", "NFE2L2",
    "IDH1", "IDH2", "VHL", "SETD2", "BAP1", "PBRM1",
    "GATA3", "CDH1", "MYC", "CCNE1", "NOTCH1", "NOTCH2",
    "SMAD4", "ARID2", "KMT2C", "SPOP", "CIC",
]


def main():
    print("█"*60)
    print("  Fine-Grained Solid Tumor Paralog-SL Analysis")
    print(f"  {len(SOLID_TUMORS)} specific cancer types")
    print("█"*60 + "\n")

    dep = load_dependency(); expr = load_expression()
    mod = load_models(); mut = load_mutations(); para = load_paralogs()

    known_set = set()
    for a, b in KNOWN_PARALOG_SL:
        known_set.add((a.upper(), b.upper()))
        known_set.add((b.upper(), a.upper()))

    # ── Performance: one shared PCS object (necessity is global) and one
    # mutation matrix per driver over the union of valid cell lines, sliced
    # per lineage (build_mutation_matrix is row-independent per cell line).
    pcs = ParalogCompensationScore(dep, expr, mod, mut, para)

    lineage_cells = {}
    for cancer_name, patterns in SOLID_TUMORS.items():
        pat = "|".join(patterns)
        mask = mod["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
        lin_mod = mod[mask]
        if len(lin_mod) < 6:
            continue
        cell_ids = [c for c in lin_mod["DepMap_ID"].tolist()
                    if c in dep.index and c in expr.index]
        if len(cell_ids) >= 6:
            lineage_cells[cancer_name] = cell_ids

    from data_loader import build_mutation_matrix
    all_cells = sorted({c for ids in lineage_cells.values() for c in ids})
    drivers = [d for d in SOLID_DRIVERS if d in dep.columns]
    print(f"  Precomputing mutation matrices for {len(drivers)} drivers "
          f"over {len(all_cells)} cell lines ...")
    mut_mats = {d: build_mutation_matrix(mut, all_cells, [d]) for d in drivers}

    all_results = {}
    summary_rows = []

    for cancer_name, cell_ids in lineage_cells.items():
        print(f"  {cancer_name:25s}: {len(cell_ids):>3d} valid lines")

        results_list = []
        for driver in drivers:
            try:
                result = pcs.compute_pcs_for_driver(
                    driver, cell_ids, cancer_label=cancer_name,
                    mut_matrix=mut_mats[driver])
                if result is not None and not result.empty:
                    results_list.append(result)
            except Exception:
                pass

        if not results_list:
            print(f"    No driver with ≥3 MUT lines")
            continue

        results_df = pd.concat(results_list, ignore_index=True)
        results_df["is_known_paralog_sl"] = results_df.apply(
            lambda r: (r["driver_gene"].upper(), r["paralog_gene"].upper()) in known_set, axis=1
        )

        all_results[cancer_name] = results_df

        yt = results_df["is_known_paralog_sl"].astype(int).values
        ys = results_df["dependency_dd"].fillna(0).values
        nk = int(yt.sum())
        auc = roc_auc_score(yt, ys) if nk >= 2 else float("nan")

        summary_rows.append({
            "cancer": cancer_name, "n_lines": len(cell_ids),
            "n_pairs": len(results_df), "n_known": nk, "dd_auroc": auc
        })

        if np.isnan(auc):
            print(f"    {len(results_df):>4d} pairs, {nk} known (insufficient for AUC)")
        else:
            print(f"    {len(results_df):>4d} pairs, {nk} known, DD AUROC={auc:.3f}")

    summary = pd.DataFrame(summary_rows).sort_values("dd_auroc", ascending=False, na_position="last")
    summary.to_csv(OUTPUT_DIR / "solid_tumor_summary.csv", index=False)

    print(f"\n{'='*70}")
    print(f"Solid Tumor DD AUROC Summary")
    print(f"{'='*70}")
    print(f"{'Cancer Type':25s}  {'Lines':>5s}  {'Pairs':>5s}  {'Known':>5s}  {'AUROC':>7s}")
    print("-"*60)
    for _, r in summary.iterrows():
        astr = f"{r['dd_auroc']:.3f}" if not np.isnan(r['dd_auroc']) else "N/A"
        flag = "★" if (not np.isnan(r['dd_auroc']) and r['dd_auroc'] > 0.7) else "·"
        print(f"  {flag} {r['cancer']:23s}  {int(r['n_lines']):>5d}  {int(r['n_pairs']):>5d}  "
              f"{int(r['n_known']):>5d}  {astr:>7s}")

    for name, df in all_results.items():
        safe = name.replace("/", "_").replace(" ", "_")
        df.to_csv(OUTPUT_DIR / f"solid_{safe}_results.csv", index=False)

    return summary, all_results


if __name__ == "__main__":
    main()
