"""Build per-lineage DD AUROC summaries from the 23 per-cancer result files.

Outputs
-------
output/solid_tumor_summary.csv        primary frame (>=5 mutant & >=5 WT per stratum)
output/solid_tumor_summary_min3.csv   sensitivity frame (all tested pairs, >=3 rule)

Scoring convention (2026-07-28): AUROC of SIGNED DD is the primary metric
(positive DD = paralog more essential in driver-mutant lines). The |DD|
magnitude convention is retained as a sensitivity analysis column
(dd_auroc_abs_sensitivity) because unsigned ranking also promotes
reverse-direction effects (e.g. PIK3CA->PIK3CB, where mutant lines are
*less* dependent on the paralog).

Inputs are the per-cancer files output/solid_<Cancer>_results.csv written by
pancancer.py (full DepMap reanalysis). Cell-line counts (n_lines) are carried
over from the previous summary files; pair/positive counts are recomputed.
"""
import glob
import sys

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"


def per_cancer_auroc(df, score_col="dependency_dd"):
    y = df["is_known_paralog_sl"].astype(int).values
    if y.sum() < 2:
        return float("nan")
    return float(roc_auc_score(y, df[score_col].fillna(0).values))


def main():
    old_primary = OUT / "solid_tumor_summary.csv"
    old_min3 = OUT / "solid_tumor_summary_min3.csv"
    n_lines_map = {}
    for old in (old_primary, old_min3):
        if old.exists():
            for r in pd.read_csv(old).itertuples():
                n_lines_map.setdefault(r.cancer, r.n_lines)

    files = sorted(glob.glob(str(OUT / "solid_*_results.csv")))
    if not files:
        sys.exit("ERROR: no output/solid_*_results.csv files — run pancancer.py first")

    rows_primary, rows_min3 = [], []
    for f in files:
        d = pd.read_csv(f)
        ct = d["cancer_type"].iloc[0]
        n_lines = n_lines_map.get(ct, float("nan"))

        # Sensitivity (min3) frame: all tested pairs
        rows_min3.append({
            "cancer": ct, "n_lines": n_lines, "n_pairs": len(d),
            "n_known": int(d["is_known_paralog_sl"].sum()),
            "dd_auroc": per_cancer_auroc(d),
            "dd_auroc_abs_sensitivity": per_cancer_auroc(d.assign(dependency_dd=d["dependency_dd"].abs())),
        })

        # Primary frame: >=5 mutant and >=5 WT cell lines per stratum
        p = d[(d["n_mut"] >= 5) & (d["n_wt"] >= 5)]
        rows_primary.append({
            "cancer": ct, "n_lines": n_lines, "n_pairs": len(p),
            "n_known": int(p["is_known_paralog_sl"].sum()),
            "dd_auroc": per_cancer_auroc(p),
            "dd_auroc_abs_sensitivity": per_cancer_auroc(p.assign(dependency_dd=p["dependency_dd"].abs())),
        })

    prim = pd.DataFrame(rows_primary).sort_values("dd_auroc", ascending=False, na_position="last")
    min3 = pd.DataFrame(rows_min3).sort_values("dd_auroc", ascending=False, na_position="last")
    prim.to_csv(old_primary, index=False)
    min3.to_csv(old_min3, index=False)

    n_eval_p = int(prim["dd_auroc"].notna().sum())
    n_eval_m = int(min3["dd_auroc"].notna().sum())
    print(f"Wrote {old_primary}  ({n_eval_p} evaluable lineages, primary >=5 frame)")
    print(f"Wrote {old_min3}  ({n_eval_m} evaluable lineages, min3 sensitivity frame)")
    print(f"  primary: {(prim['dd_auroc'] > 0.7).sum()}/{n_eval_p} lineages with signed DD AUROC > 0.7")
    print(f"  min3   : {(min3['dd_auroc'] > 0.7).sum()}/{n_eval_m} lineages with signed DD AUROC > 0.7")


if __name__ == "__main__":
    main()
