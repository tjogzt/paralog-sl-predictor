"""
Task 4 — Stratified bootstrap test of MSI-H vs MSS DD-AUROC differences
========================================================================
Manuscript (MSI section, >=3-per-group sensitivity frame, signed DD):
  endometrial AUROC = 0.303 (MSI-H, n=17 lines) vs 0.500 (MSS, n=11),
  colorectal  AUROC = 0.558 (MSI-H, n=14) vs 0.545 (MSS, n=45);
  the differences were not formally tested. This script performs that test.

Data: pair-level subgroup results written by the min3 sensitivity run
(output/msi_<cancer>_<subgroup>_results_min3.csv). Score = signed
dependency_dd (primary metric; positive = compensation), as produced by
msi_analysis.py.

Procedure (seed 42, 10,000 iterations): within each subgroup, resample the
positive pairs and the negative pairs separately with replacement
(stratified, preserving the within-subgroup pos/neg structure), recompute
both AUROCs, and record delta* = AUROC(MSI-H) - AUROC(MSS). Report the
percentile 95% CI and a two-sided empirical p-value against delta = 0:
p = 2 * min(P(delta* <= 0), P(delta* >= 0)) with the (+1)/(B+1) correction.

Output: output/msi_bootstrap_test.json
"""

import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import OUTPUT_DIR

B = 10_000
SEED = 42

CANCERS = {
    "Endometrial": {
        "msi_h": OUTPUT_DIR / "msi_endometrial_msi_h_results_min3.csv",
        "mss": OUTPUT_DIR / "msi_endometrial_mss_results_min3.csv",
    },
    "Colorectal": {
        "msi_h": OUTPUT_DIR / "msi_colorectal_msi_h_results_min3.csv",
        "mss": OUTPUT_DIR / "msi_colorectal_mss_results_min3.csv",
    },
}


def load_scores(path):
    df = pd.read_csv(path)
    y = df["is_known_paralog_sl"].astype(int).values
    s = df["dependency_dd"].fillna(0).values
    return y, s


def auroc_from(y, s):
    return roc_auc_score(y, s)


def stratified_auroc(y, s, rng):
    """Resample positives and negatives separately, then AUROC."""
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    pi = rng.choice(pos, size=len(pos), replace=True)
    ni = rng.choice(neg, size=len(neg), replace=True)
    idx = np.concatenate([pi, ni])
    return roc_auc_score(y[idx], s[idx])


def main():
    rng = np.random.default_rng(SEED)
    results = []

    for cancer, files in CANCERS.items():
        y_h, s_h = load_scores(files["msi_h"])
        y_s, s_s = load_scores(files["mss"])

        auc_h = auroc_from(y_h, s_h)
        auc_s = auroc_from(y_s, s_s)
        delta = auc_h - auc_s

        deltas = np.empty(B)
        for b in range(B):
            deltas[b] = (stratified_auroc(y_h, s_h, rng)
                         - stratified_auroc(y_s, s_s, rng))

        ci_low, ci_high = np.percentile(deltas, [2.5, 97.5])
        p_le = (np.sum(deltas <= 0) + 1) / (B + 1)
        p_ge = (np.sum(deltas >= 0) + 1) / (B + 1)
        p_emp = min(1.0, 2.0 * min(p_le, p_ge))

        rec = {
            "cancer": cancer,
            "auroc_msi_h": round(float(auc_h), 4),
            "auroc_mss": round(float(auc_s), 4),
            "delta": round(float(delta), 4),
            "boot_ci_low": round(float(ci_low), 4),
            "boot_ci_high": round(float(ci_high), 4),
            "p_empirical": round(float(p_emp), 4),
            "n_pairs_evaluable": {
                "msi_h": int(len(y_h)), "mss": int(len(y_s)),
                "msi_h_positives": int(y_h.sum()), "mss_positives": int(y_s.sum()),
            },
        }
        results.append(rec)
        print(f"{cancer:12s} MSI-H={auc_h:.3f} MSS={auc_s:.3f} "
              f"delta={delta:+.3f}  95%CI [{ci_low:+.3f}, {ci_high:+.3f}]  "
              f"p_emp={p_emp:.4f}")

    out_obj = {
        "method": {
            "frame": ">=3-per-group sensitivity frame (min3 products)",
            "score": "signed dependency_dd (positive = compensation)", "label": "is_known_paralog_sl",
            "bootstrap": f"stratified (pos/neg preserved within subgroup), "
                         f"B={B}, seed={SEED}",
            "p_value": "two-sided empirical vs delta=0: "
                       "2*min(P(delta*<=0), P(delta*>=0)), (+1)/(B+1) corrected",
        },
        "results": results,
    }
    out = OUTPUT_DIR / "msi_bootstrap_test.json"
    with open(out, "w") as fh:
        json.dump(out_obj, fh, indent=2)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
