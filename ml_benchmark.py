#!/usr/bin/env python3
"""
ml_benchmark.py
===============
Reproducible head-to-head classifier benchmark for the manuscript's
"true head-to-head comparison on the same 77 paralog pairs (8 known
positives)" paragraph.

Universe and aggregation (same as compute_headline_metrics.py):
  * input: output/tables/TableS2_FullResults.tsv (118 gyn3 lineage entries)
  * aggregate to 77 unique driver->paralog pairs by MEAN across lineages
    (this is the frame on which "DD alone (0.736, using only |DD|)"
    reproduces as 0.7355)
  * 8 known positives (gold-standard paralog-SL pairs)

Features (per manuscript): |DD|, PCS, |deltaExpression|, necessity,
composite score, mutation frequency. Features are z-score standardized.

Classifiers: logistic regression, random forest (seeded), SVM-RBF,
SVM-Linear. Evaluation: leave-one-pair-out cross-validation (77 folds),
AUROC over out-of-fold decision scores.

The LR coefficient test ("standardized beta = 2.64, p = 0.009" in the
manuscript) is a statsmodels Logit fit on all 77 pairs with standardized
features (Wald z-test).

Outputs:
  output/ml_benchmark.json
  output/tables/ml_benchmark.tsv

Usage: python ml_benchmark.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.svm import SVC

ROOT = Path(__file__).resolve().parent
TABLES2 = ROOT / "output" / "tables" / "TableS2_FullResults.tsv"
JSON_OUT = ROOT / "output" / "ml_benchmark.json"
TSV_OUT = ROOT / "output" / "tables" / "ml_benchmark.tsv"

# composite is a deterministic function of (dd_abs, pcs, dexpr_abs,
# necessity); including it makes the linear design matrix collinear and
# destabilizes LR / linear-SVM fits. It is therefore evaluated only as a
# single-feature baseline (composite_alone), never as a classifier input.
FEATURES = ["dd_abs", "pcs", "dexpr_abs", "necessity", "mut_freq"]
SEED = 42

# Values stated in the manuscript, for the automated claims check.
# Updated 2026-07-26 to the leave-one-pair-out results after the C7
# class-specific driver-mutation rules shrank the evaluable positive set
# to 6 aggregated pairs; composite removed from classifier features
# (deterministic function of the other features -> collinearity) and kept
# as single-feature baseline only. With n_pos = 6, linear-classifier LOO
# AUROCs are unstable and are reported with an explicit small-n caveat.
CLAIMS = {
    "LR": 0.138,
    "RF": 0.617,
    "SVM_RBF_low": 0.744,
    "SVM_Linear_high": 0.217,
    "dd_alone": 0.551,
    "composite_alone": 0.841,
    "lr_beta_dd": 0.40,
    "lr_p_dd": 0.253,
}
TOL = 0.005
TOL_COEF = 0.02


def auroc(labels, scores):
    """Rank-based AUROC (average ranks for ties); NaN if a class is absent."""
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    ok = ~np.isnan(scores)
    labels, scores = labels[ok], scores[ok]
    n_pos, n_neg = int(labels.sum()), int((1 - labels).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = pd.Series(scores).rank(method="average").values
    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def build_pair_frame(df):
    """Aggregate lineage entries to unique pairs (mean across lineages)."""
    g = df.groupby(["driver_gene", "paralog_gene"], as_index=False).agg(
        known=("is_known_paralog_sl", "max"),
        dd_abs=("dependency_dd", lambda s: s.abs().mean()),
        pcs=("pcs", "mean"),
        dexpr_abs=("delta_expression", lambda s: s.abs().mean()),
        necessity=("necessity", "mean"),
        composite=("composite_score", "mean"),
        mut_freq=("mutation_frequency", "mean"),
    )
    return g


def loo_scores(clf, X, y):
    """Leave-one-pair-out CV; returns out-of-fold decision scores."""
    loo = LeaveOneOut()
    scores = np.full(len(y), np.nan)
    for train_idx, test_idx in loo.split(X):
        clf.fit(X[train_idx], y[train_idx])
        if hasattr(clf, "decision_function"):
            s = clf.decision_function(X[test_idx])
        else:
            s = clf.predict_proba(X[test_idx])[:, 1]
        scores[test_idx] = s[0]
    return scores


def main():
    if not TABLES2.exists():
        sys.exit(f"ERROR: {TABLES2} not found — run main.py + tables.py first")

    df = pd.read_csv(TABLES2, sep="\t")
    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)
    pairs = build_pair_frame(df)

    y = pairs["known"].astype(int).values
    X = pairs[FEATURES].fillna(0.0).values.astype(float)
    n_pairs, n_pos = len(pairs), int(y.sum())
    print(f"Universe: {n_pairs} unique pairs, {n_pos} positives "
          f"(expected 77/8)")

    # Standardize features
    mu, sd = X.mean(axis=0), X.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd

    classifiers = {
        "LR": LogisticRegression(max_iter=10000),
        "RF": RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1),
        "SVM_RBF": SVC(kernel="rbf"),
        "SVM_Linear": SVC(kernel="linear"),
    }

    results = {}
    for name, clf in classifiers.items():
        scores = loo_scores(clf, Xs, y)
        results[name] = {"auroc": auroc(y, scores)}
        print(f"  {name:12s} LOO AUROC = {results[name]['auroc']:.4f}")

    # Single-feature references (no training)
    dd_alone = auroc(y, pairs["dd_abs"].fillna(0).values)
    comp_alone = auroc(y, pairs["composite"].fillna(0).values)
    print(f"  {'DD alone':12s} AUROC = {dd_alone:.4f}")
    print(f"  {'Composite alone':12s} AUROC = {comp_alone:.4f}")

    # LR standardized coefficients + Wald p-values (statsmodels MLE fit)
    import statsmodels.api as sm
    coef_info = {}
    try:
        model = sm.Logit(y, sm.add_constant(Xs)).fit(disp=0)
        params = model.params
        pvals = model.pvalues
        for i, feat in enumerate(["const"] + FEATURES):
            coef_info[feat] = {"beta": float(params[i]), "p_value": float(pvals[i])}
        print(f"  LR beta(|DD|) = {coef_info['dd_abs']['beta']:.3f}, "
              f"p = {coef_info['dd_abs']['p_value']:.4f}")
    except Exception as e:
        coef_info = {"status": f"fit failed: {e}"}
        print(f"  [WARN] statsmodels Logit failed: {e}")

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "ml_benchmark.py",
        "source": str(TABLES2.relative_to(ROOT)),
        "universe": {"n_pairs": n_pairs, "n_positives": n_pos,
                     "aggregation": "mean across gyn3 lineages per pair"},
        "features": FEATURES,
        "cv": "leave-one-pair-out (77 folds), AUROC on out-of-fold scores",
        "seed": SEED,
        "classifiers": results,
        "single_feature": {"dd_alone": dd_alone, "composite_alone": comp_alone},
        "lr_coefficients": coef_info,
    }

    # Claims check
    computed = {
        "LR": results["LR"]["auroc"],
        "RF": results["RF"]["auroc"],
        "SVM_RBF_low": results["SVM_RBF"]["auroc"],
        "SVM_Linear_high": results["SVM_Linear"]["auroc"],
        "dd_alone": dd_alone,
        "composite_alone": comp_alone,
        "lr_beta_dd": (coef_info.get("dd_abs") or {}).get("beta"),
        "lr_p_dd": (coef_info.get("dd_abs") or {}).get("p_value"),
    }
    checks = []
    print("\n=== Manuscript claims check ===")
    for name, claimed in CLAIMS.items():
        got = computed.get(name)
        if got is None or (isinstance(got, float) and np.isnan(got)):
            status = "not_reproducible"
        else:
            tol = TOL_COEF if name.startswith("lr_") else TOL
            status = "match" if abs(got - claimed) <= tol else "MISMATCH"
        checks.append({"metric": name, "claimed": claimed,
                       "computed": None if got is None else round(float(got), 4),
                       "status": status})
        gs = "—" if got is None else f"{got:.4f}"
        print(f"  [{status:>17s}] {name}: claimed {claimed} vs computed {gs}")
    out["manuscript_claims_check"] = checks

    JSON_OUT.write_text(json.dumps(out, indent=2, allow_nan=False, default=str))
    print(f"\nWrote {JSON_OUT}")

    rows = [{"metric": k, "value": (None if v is None else f"{float(v):.4f}"),
             "claimed": CLAIMS[k]}
            for k, v in computed.items()]
    TSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TSV_OUT, sep="\t", index=False)
    print(f"Wrote {TSV_OUT}")

    n_bad = sum(1 for c in checks if c["status"] == "MISMATCH")
    if n_bad:
        sys.exit(f"FAILED: {n_bad} manuscript claim(s) MISMATCH")


if __name__ == "__main__":
    main()
