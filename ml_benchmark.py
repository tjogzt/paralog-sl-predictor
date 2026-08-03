#!/usr/bin/env python3
"""
ml_benchmark.py
===============
Reproducible head-to-head classifier benchmark for the manuscript's
"true head-to-head comparison on the same 72 paralog pairs (6 known
positives)" paragraph.

Universe and aggregation (same as compute_headline_metrics.py):
  * input: output/tables/TableS2_FullResults.tsv (gyn3 lineage entries)
  * aggregate to 72 unique driver->paralog pairs by MEAN across lineages
    (this is the frame on which "DD alone (0.736, using only |DD|)"
    reproduces as 0.7355)
  * 6 known positives (gold-standard paralog-SL pairs)

Features (per manuscript): |DD|, PCS, |deltaExpression|, necessity,
composite score, mutation frequency. Features are z-score standardized
INSIDE each CV fold (scaler fit on training pairs only; fixed 2026-07-31 —
previously standardized once on all pairs before CV, leaking full-data
moments).

The composite score shipped in TableS2 (pcs.py) min-max normalizes its
four components on the FULL dataset (mild leakage flagged in review).
composite_fold_scores() therefore recomputes the composite leave-one-pair-out:
for each held-out pair, component min/max come from the 71 TRAINING pairs
only. Reported as single_feature.composite_alone_lofo (AUROC) and
composite_auprc_lofo (AUPRC); the full-data composite_alone is kept
unchanged for comparison.

Classifiers: logistic regression, random forest (seeded), SVM-RBF,
SVM-Linear. Evaluation: leave-one-pair-out cross-validation (72 folds),
AUROC over out-of-fold decision scores.

The LR coefficient test ("standardized beta = 2.64, p = 0.009" in the
manuscript) is a statsmodels Logit fit on all 72 pairs with standardized
features (Wald z-test).

Method-difference uncertainty (round-7 review): a paired bootstrap over
the 72 pairs (10,000 resamples with replacement, same resample for both
methods per iteration, seed 42) yields the distribution of
ΔAUROC = composite LOPO − SVM-RBF LOPO; mean Δ, the 95% percentile CI,
and P(Δ<0) are written to the `paired_bootstrap_delta_auroc` key.

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
from sklearn.metrics import average_precision_score
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
# Updated 2026-07-26 to the min>=5/mutation-rule framework: 72 unique
# pairs, 6 positives; SVM-RBF and the composite are statistically
# indistinguishable at this sample size (paired resamples overlap).
# Updated 2026-07-31 after the fold-internal-standardization fix in
# loo_scores() (scaler now fitted on the 71 training pairs of each fold
# instead of once on all 72 pairs before CV, which leaked full-data
# moments). RF and LR out-of-fold scores are unchanged (RF is invariant
# to per-fold monotone feature rescaling); SVM_RBF moved 0.8434 -> 0.8409
# (within tolerance); SVM_Linear moved 0.1136 -> 0.2399, so its CLAIMS
# entry is updated 0.114 -> 0.240. Single-feature baselines and the
# full-data LR coefficient fit are unaffected by the fix.
CLAIMS = {
    "LR": 0.136,
    "RF": 0.722,
    "SVM_RBF_low": 0.843,
    "SVM_Linear_high": 0.240,
    "dd_alone": 0.672,
    "composite_alone": 0.831,
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
        dd_signed=("dependency_dd", "mean"),
        pcs=("pcs", "mean"),
        dexpr_abs=("delta_expression", lambda s: s.abs().mean()),
        necessity=("necessity", "mean"),
        composite=("composite_score", "mean"),
        mut_freq=("mutation_frequency", "mean"),
        q_value=("q_value", "mean"),
    )
    return g


def composite_fold_scores(pairs):
    """Leakage-free composite: min-max refit leave-one-pair-out.

    The shipped composite_score column (pcs.py) is
      0.50*minmax(pcs) + 0.20*minmax(|DD|) + 0.15*minmax(-log10(q+1e-10))
      + 0.15*minmax(mutation_frequency),
    with each min-max taken over the FULL dataset, so every pair's score
    embeds global extrema (mild leakage flagged in review). Here, for each
    held-out pair i, component min/max are computed on the 71 TRAINING
    pairs only and pair i is scored as the weighted sum of
    (v_i - train_min) / (train_max - train_min) per component (zero-range
    components guarded with denominator 1.0). Returns the 72 held-out
    composite scores in the row order of `pairs`."""
    comp = pd.DataFrame({
        "pcs": pairs["pcs"].astype(float),
        "dd_abs": pairs["dd_abs"].astype(float),
        "neglog10_q": -np.log10(pairs["q_value"].astype(float) + 1e-10),
        "mut_freq": pairs["mut_freq"].astype(float),
    }).fillna(0.0)
    weights = np.array([0.50, 0.20, 0.15, 0.15])  # pcs, dd_abs, neglog10_q, mut_freq
    V = comp.values
    n = len(comp)
    scores = np.full(n, np.nan)
    for i in range(n):
        train = np.delete(V, i, axis=0)
        cmin = train.min(axis=0)
        cmax = train.max(axis=0)
        rng = cmax - cmin
        rng[rng == 0] = 1.0  # zero-range guard
        scores[i] = float((((V[i] - cmin) / rng) * weights).sum())
    return scores


def loo_scores(clf, X, y):
    """Leave-one-pair-out CV; returns out-of-fold decision scores.

    Fold-internal z-score standardization (fixed 2026-07-31): the scaler is
    fitted on the 71 TRAINING pairs of each fold only and applied to both
    train and the held-out pair. Previously features were standardized once
    on all 72 pairs before CV, leaking full-data moments (mean/variance of
    the held-out pair) into every fold."""
    loo = LeaveOneOut()
    scores = np.full(len(y), np.nan)
    for train_idx, test_idx in loo.split(X):
        mu = X[train_idx].mean(axis=0)
        sd = X[train_idx].std(axis=0, ddof=1)
        sd[sd == 0] = 1.0
        Xtr = (X[train_idx] - mu) / sd
        Xte = (X[test_idx] - mu) / sd
        clf.fit(Xtr, y[train_idx])
        if hasattr(clf, "decision_function"):
            s = clf.decision_function(Xte)
        else:
            s = clf.predict_proba(Xte)[:, 1]
        scores[test_idx] = s[0]
    return scores


BOOT_DELTA_ITERS = 10000
BOOT_DELTA_SEED = 42


def paired_bootstrap_delta_auroc(y, scores_a, scores_b,
                                 iters=BOOT_DELTA_ITERS, seed=BOOT_DELTA_SEED):
    """Paired bootstrap CI for an AUROC difference between two methods.

    Pairs are resampled with replacement; the SAME resample is used for
    both score vectors in every iteration, so the distribution of
    delta = AUROC(a) - AUROC(b) is properly paired. Resamples containing
    only one class are skipped (AUROC undefined). Returns mean delta,
    the 95% percentile CI, and P(delta < 0)."""
    rng = np.random.default_rng(seed)
    n = len(y)
    deltas = []
    n_skipped = 0
    for _ in range(iters):
        idx = rng.integers(0, n, n)
        yb = y[idx]
        if yb.sum() == 0 or yb.sum() == n:
            n_skipped += 1
            continue
        deltas.append(auroc(yb, scores_a[idx]) - auroc(yb, scores_b[idx]))
    d = np.asarray(deltas)
    return {
        "n_iterations": iters,
        "n_valid": int(len(d)),
        "n_skipped_single_class": int(n_skipped),
        "seed": seed,
        "mean_delta": float(d.mean()),
        "ci95": [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
        "p_delta_lt_0": float((d < 0).mean()),
    }


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
          f"(expected 72/6 under min>=5)")

    # Full-data standardization, used ONLY for the statsmodels LR coefficient
    # fit below (a full-sample fit, so no leakage concern). Cross-validation
    # uses fold-internal standardization inside loo_scores().
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
    loo_score_vectors = {}
    for name, clf in classifiers.items():
        scores = loo_scores(clf, X, y)
        results[name] = {"auroc": auroc(y, scores)}
        loo_score_vectors[name] = scores
        print(f"  {name:12s} LOO AUROC = {results[name]['auroc']:.4f}")

    # Single-feature references (no training)
    dd_alone = auroc(y, pairs["dd_signed"].fillna(0).values)
    dd_abs_alone = auroc(y, pairs["dd_abs"].fillna(0).values)
    comp_alone = auroc(y, pairs["composite"].fillna(0).values)
    print(f"  {'DD alone':12s} AUROC = {dd_alone:.4f}")
    print(f"  {'|DD| (sensitivity)':12s} AUROC = {dd_abs_alone:.4f}")
    print(f"  {'Composite alone':12s} AUROC = {comp_alone:.4f}")

    # Fold-internal composite (leakage-free): min-max refit per LOFO fold
    comp_lofo_scores = composite_fold_scores(pairs)
    comp_lofo = auroc(y, comp_lofo_scores)
    comp_lofo_auprc = float(average_precision_score(y, comp_lofo_scores))
    print(f"  {'Composite LOFO':12s} AUROC = {comp_lofo:.4f}, "
          f"AUPRC = {comp_lofo_auprc:.4f}")

    # Paired bootstrap: composite (LOPO) vs SVM-RBF (LOPO) AUROC difference.
    # The same pair resample is applied to both out-of-fold score vectors in
    # each of 10,000 iterations, so the Δ distribution is properly paired.
    delta_boot = paired_bootstrap_delta_auroc(
        y, comp_lofo_scores, loo_score_vectors["SVM_RBF"])
    delta_boot["comparison"] = "composite LOPO minus SVM-RBF LOPO (out-of-fold scores)"
    print(f"  ΔAUROC (composite − SVM-RBF): mean {delta_boot['mean_delta']:+.4f}, "
          f"95% CI [{delta_boot['ci95'][0]:+.4f}, {delta_boot['ci95'][1]:+.4f}], "
          f"P(Δ<0) = {delta_boot['p_delta_lt_0']:.4f} "
          f"({delta_boot['n_valid']} valid resamples)")

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
        "cv": "leave-one-pair-out (72 folds), fold-internal z-score standardization "
              "(scaler fit on training pairs only), AUROC on out-of-fold scores",
        "seed": SEED,
        "classifiers": results,
        "single_feature": {"dd_alone": dd_alone, "dd_abs_alone_sensitivity": dd_abs_alone,
                           "composite_alone": comp_alone,
                           "composite_alone_lofo": comp_lofo},
        "composite_auprc_lofo": comp_lofo_auprc,
        "paired_bootstrap_delta_auroc": delta_boot,
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
