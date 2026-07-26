"""
Task 3 — PRISM delta-AUC for the five manuscript drug associations
===================================================================
Manuscript (Results, PRISM paragraph; Fig. 4a):
  - AZD8330     d=-0.62, q=0.0007 (KRAS-mutant)
  - Trametinib  d=-0.57, q=0.0006 (KRAS-mutant)
  - Everolimus  d=-0.61, q=0.0009 (PTEN-mutant)
  - Ipatasertib d=-0.99, q=0.0001 (PTEN-mutant)
  - Panobinostat d=-1.45, q=0.010 (EP300-mutant, Ovarian)

All numbers are extracted from the existing analysis product
output/prism_full_results.csv (produced by prism_analysis.py from
data/PRISM_log2AUC.csv; delta_auc = mean(log2AUC|MUT) - mean(log2AUC|WT),
negative = selective killing of mutant lines). No recomputation from raw
data is needed because the product stores delta_auc alongside Cohen's d.

Output: output/prism_delta_auc.json
"""

import json

import pandas as pd

from config import OUTPUT_DIR

# (drug_substring, driver, context) — paralog recorded from the matched row
TARGETS = [
    ("AZD8330", "KRAS", "PanCancer"),
    ("TRAMETINIB", "KRAS", "PanCancer"),
    ("EVEROLIMUS", "PTEN", "PanCancer"),
    ("IPATASERTIB", "PTEN", "PanCancer"),
    ("PANOBINOSTAT", "EP300", "Ovarian"),
]

MANUSCRIPT = {
    "AZD8330": {"d": -0.62, "q": 0.0007},
    "TRAMETINIB": {"d": -0.57, "q": 0.0006},
    "EVEROLIMUS": {"d": -0.61, "q": 0.0009},
    "IPATASERTIB": {"d": -0.99, "q": 0.0001},
    "PANOBINOSTAT": {"d": -1.45, "q": 0.010},
}


def main():
    df = pd.read_csv(OUTPUT_DIR / "prism_full_results.csv")
    df["drug_upper"] = df["drug"].str.upper()

    records = []
    for drug, driver, context in TARGETS:
        sub = df[(df["drug_upper"] == drug)
                 & (df["driver"] == driver)
                 & (df["context"] == context)]
        if sub.empty:
            records.append({"drug": drug, "driver": driver, "context": context,
                            "error": "not found in prism_full_results.csv"})
            continue
        # Multiple paralog rows can share the same driver-genotype groups
        # (e.g. KRAS->HRAS and KRAS->NRAS are identical); keep the first and
        # note duplicates.
        sub = sub.sort_values("paralog")
        row = sub.iloc[0]
        dup_pairs = sorted(f"{r.driver}->{r.paralog}" for r in sub.itertuples())
        m = MANUSCRIPT[drug]
        rec = {
            "drug": row["drug"],
            "driver": row["driver"],
            "paralog": row["paralog"],
            "context": row["context"],
            "identical_rows_same_driver_genotype": dup_pairs,
            "cohens_d": round(float(row["cohens_d"]), 6),
            "delta_auc": round(float(row["delta_auc"]), 6),
            "q": float(row["bh_q"]),
            "p_value": float(row["p_value"]),
            "n_mut": int(row["n_mut"]),
            "n_wt": int(row["n_wt"]),
            "mean_log2auc_mut": float(row["mean_mut"]),
            "mean_log2auc_wt": float(row["mean_wt"]),
            "manuscript_d": m["d"],
            "manuscript_q": m["q"],
            "d_matches_manuscript": abs(round(float(row["cohens_d"]), 2) - m["d"]) < 1e-9,
            "q_matches_manuscript_rounded": abs(round(float(row["bh_q"]), 4) - m["q"]) < 5e-5,
        }
        records.append(rec)

    result = {
        "source": "output/prism_full_results.csv (prism_analysis.py; "
                  "delta_auc = mean log2AUC MUT - mean log2AUC WT)",
        "note": "Drug selectivity is conditioned on driver genotype only; the "
                "paralog column identifies the pair row in the product. "
                "n_mut/n_wt differ across drugs because PRISM compound "
                "coverage differs per cell line.",
        "associations": records,
    }

    out = OUTPUT_DIR / "prism_delta_auc.json"
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)

    for r in records:
        if "error" not in r:
            print(f"{r['drug']:14s} {r['driver']:6s} ({r['context']:9s}) "
                  f"deltaAUC={r['delta_auc']:+.4f}  d={r['cohens_d']:+.3f}  "
                  f"q={r['q']:.2e}  n={r['n_mut']}/{r['n_wt']}  "
                  f"match d/q: {r['d_matches_manuscript']}/{r['q_matches_manuscript_rounded']}")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
