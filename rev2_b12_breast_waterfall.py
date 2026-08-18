#!/usr/bin/env python3
"""
rev2_b12_breast_waterfall.py  (Stage-4 revision, item B12)
==========================================================
The Fig. 1a Breast lineage frame (output/solid_Breast_results.csv, primary
>=5/>=5 frame: 12 entries, 2 positives) has AUROC = 0.150 — the lowest of
the 8 evaluable lineages. This script decomposes that estimate:

  (a) per-entry waterfall: rank all 12 entries by signed DD; the AUROC
      contribution of each positive = fraction of negatives ranked below
      it (Mann-Whitney), showing exactly which positives sit at the bottom;
  (b) cross-lineage context: the same curated pairs' signed DD in the
      other evaluable lineages, to test whether the Breast anomaly is a
      pair-specific reverse displacement (PIK3CA->PIK3CB compensates
      elsewhere but shows reduced dependency in Breast-mutant lines);
  (c) lineage-level meta-regression across the 8 evaluable lineages:
      AUROC ~ log(n positives) + mean driver mutation frequency
      + log(n lines). n = 8, DESCRIPTIVE only (no inference).

Sign convention: signed DD = mean(Chronos|WT) - mean(Chronos|MUT)
(positive = compensation), matching pcs.py and the solid_* production
files (fidelity-checked in rev B1). NOTE: the legacy output/breast_results.csv
carries the SAME pairs with opposite sign (older artifact, superseded by
the solid frame); the solid frame is the manuscript frame of record.

Deterministic; no resampling.
Output: output/revision_stage4/rev2_b12_breast_waterfall.{json,csv}
Usage: python rev2_b12_breast_waterfall.py   (run from repo root)
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

from compute_headline_metrics import auroc  # noqa: E402

SOLID = ROOT / "output"


def primary_frame(df):
    return df[(df["n_mut"] >= 5) & (df["n_wt"] >= 5)].copy()


def main():
    b = pd.read_csv(SOLID / "solid_Breast_results.csv")
    sub = primary_frame(b).sort_values("dependency_dd", ascending=False).reset_index(drop=True)
    y = sub["is_known_paralog_sl"].astype(int).values
    s = sub["dependency_dd"].fillna(0).values
    auc = auroc(y, s)
    n_neg = int((y == 0).sum())

    # (a) waterfall: per-positive contribution
    rows = []
    for i, r in sub.iterrows():
        contrib = None
        if r["is_known_paralog_sl"]:
            below = int(((sub["dependency_dd"].fillna(0) < r["dependency_dd"])
                         & (~sub["is_known_paralog_sl"])).sum())
            ties = int(((sub["dependency_dd"].fillna(0) == r["dependency_dd"])
                        & (~sub["is_known_paralog_sl"])).sum())
            contrib = (below + 0.5 * ties) / n_neg
        rows.append({
            "rank": i + 1,
            "driver_gene": r["driver_gene"], "paralog_gene": r["paralog_gene"],
            "signed_dd": float(r["dependency_dd"]),
            "n_mut": int(r["n_mut"]), "n_wt": int(r["n_wt"]),
            "is_known_paralog_sl": bool(r["is_known_paralog_sl"]),
            "auroc_contribution_if_positive": contrib,
        })
    wf = pd.DataFrame(rows)
    wf.to_csv(OUT / "rev2_b12_breast_waterfall.csv", index=False)
    pos = wf[wf["is_known_paralog_sl"]]
    print(f"Breast primary frame: {len(sub)} entries, {int(y.sum())} positives, "
          f"AUROC={auc:.4f}")
    print(pos[["rank", "driver_gene", "paralog_gene", "signed_dd",
               "auroc_contribution_if_positive"]].to_string(index=False))
    assert abs(auc - pos["auroc_contribution_if_positive"].mean()) < 1e-12

    # (b) cross-lineage context for the two Breast positives
    eval_lineages = ["NSCLC", "SCLC", "Colorectal", "Esophagogastric", "Breast",
                     "Ovarian", "Endometrial", "Bladder Urothelial"]
    targets = [("PIK3CA", "PIK3CB"), ("BRCA1", "BRCA2")]
    ctx = []
    for lin in eval_lineages:
        safe = lin.replace(" ", "_")
        f = SOLID / f"solid_{safe}_results.csv"
        df = primary_frame(pd.read_csv(f))
        for drv, par in targets:
            hit = df[(df["driver_gene"] == drv) & (df["paralog_gene"] == par)]
            ctx.append({
                "lineage": lin, "pair": f"{drv}->{par}",
                "signed_dd": (float(hit["dependency_dd"].iloc[0]) if len(hit) else None),
                "evaluable": bool(len(hit)),
            })
    ctx_df = pd.DataFrame(ctx)
    print("\ncross-lineage signed DD of the two Breast-positive pairs:")
    print(ctx_df.pivot(index="lineage", columns="pair", values="signed_dd").to_string())

    # (c) meta-regression across the 8 evaluable lineages
    mrows = []
    b1 = pd.read_csv(OUT / "b1_lineage_three_stats.csv").set_index("cancer")
    for lin in eval_lineages:
        safe = lin.replace(" ", "_")
        df = primary_frame(pd.read_csv(SOLID / f"solid_{safe}_results.csv"))
        mrows.append({
            "lineage": lin,
            "auroc": float(b1.loc[lin, "auroc_signed_dd"]),
            "n_positives": int(df["is_known_paralog_sl"].sum()),
            "mean_mutation_frequency": float(df["mutation_frequency"].mean()),
            "n_lines": int(b1.loc[lin, "n_lines"]),
            "n_entries": len(df),
        })
    meta = pd.DataFrame(mrows)
    X = np.column_stack([
        np.ones(len(meta)),
        np.log(meta["n_positives"]),
        meta["mean_mutation_frequency"],
        np.log(meta["n_lines"]),
    ])
    yv = meta["auroc"].values
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    pred = X @ beta
    ss_res = float(((yv - pred) ** 2).sum())
    r2 = 1 - ss_res / float(((yv - yv.mean()) ** 2).sum())
    meta["fitted"] = pred
    meta["residual"] = yv - pred
    meta.to_csv(OUT / "rev2_b12_lineage_meta_regression.csv", index=False)
    print("\nmeta-regression (n=8, descriptive): coefficients "
          f"intercept={beta[0]:.3f}, log(n_pos)={beta[1]:.3f}, "
          f"mut_freq={beta[2]:.3f}, log(n_lines)={beta[3]:.3f}; R^2={r2:.3f}")
    print(meta[["lineage", "auroc", "n_positives", "mean_mutation_frequency",
                "n_lines", "residual"]].to_string(index=False))

    out = {
        "breast_frame": {
            "source": "output/solid_Breast_results.csv (Fig. 1a frame), primary >=5/>=5",
            "n_entries": len(sub), "n_positives": int(y.sum()), "auroc": auc,
            "waterfall": rows,
            "headline": ("both curated positives show NEGATIVE signed DD in Breast: "
                         "PIK3CA->PIK3CB DD=-0.258 (rank 11/12, AUROC contribution 0.10) "
                         "and BRCA1->BRCA2 DD=-0.136 (rank 9/12, contribution 0.20); "
                         "the reverse-displaced PIK3CA->PIK3CB pair alone contributes "
                         "half of the AUROC deficit"),
        },
        "cross_lineage_context": ctx,
        "meta_regression": {
            "model": "AUROC ~ log(n_positives) + mean_mutation_frequency + log(n_lines)",
            "n_lineages": len(meta), "descriptive_only": True,
            "coefficients": {"intercept": float(beta[0]),
                             "log_n_positives": float(beta[1]),
                             "mean_mutation_frequency": float(beta[2]),
                             "log_n_lines": float(beta[3])},
            "r_squared": r2,
            "rows": mrows,
        },
        "sign_note": ("legacy output/breast_results.csv carries the same pairs with "
                      "opposite DD sign (superseded artifact); the solid frame is the "
                      "manuscript frame of record and is sign-fidelity-checked"),
    }
    out_path = OUT / "rev2_b12_breast_waterfall.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nwritten: {out_path}")


if __name__ == "__main__":
    main()
