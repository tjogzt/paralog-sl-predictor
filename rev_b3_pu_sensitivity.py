#!/usr/bin/env python3
"""
rev_b3_pu_sensitivity.py  (Stage-4 revision, item B3)
======================================================
Positive-unlabeled (PU) quantitative sensitivity of the primary-frame
signed-DD AUROC (0.629, 110 entries: 8 labeled positives / 102 unlabeled
controls treated as negatives).

Design (pre-specified here):
  * Assumed hidden-positive rate pi in {0.01, 0.05, 0.10} among the 102
    unlabeled entries -> k = round(pi * 102) = {1, 5, 10} hidden positives.
  * Random-flip resampling: 2,000 replicates per pi (seed 42); in each
    replicate k distinct unlabeled entries are chosen uniformly at random
    and relabeled positive; AUROC is recomputed. Report mean, 2.5-97.5%
    percentile interval, and min.
  * Adversarial worst case (deterministic): the k unlabeled entries with
    the HIGHEST signed DD are flipped (hidden positives concentrate at the
    top of the ranking -> maximal AUROC damage).
  * Favourable case (deterministic): the k LOWEST-scoring unlabeled entries
    are flipped.
  * Elkan-Noto correction (optional, SCAR assumption): with labeled
    positives representative of all positives, e = P(labeled | positive)
    ~ 8 / (8 + k); corrected AUROC = 0.5 + (AUROC_obs - 0.5) / e.
    (Elkan & Noto 2008, eq. on class-prior recovery; labelled optional.)

No fabricated numbers: the only randomness is the seeded relabeling of
real TableS2 entries.

Outputs (output/revision_stage4/):
  b3_pu_sensitivity.csv     one row per pi: resampling distribution summaries
  b3_pu_sensitivity.json    full detail incl. per-replicate settings

Usage: python rev_b3_pu_sensitivity.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

from compute_headline_metrics import auroc  # noqa: E402

SEED = 42
N_REP = 2000
PIS = (0.01, 0.05, 0.10)


def main():
    print("=" * 72)
    print("  rev B3: PU quantitative sensitivity of the primary AUROC")
    print("=" * 72)

    df = pd.read_csv(ROOT / "output" / "tables" / "TableS2_FullResults.tsv", sep="\t")
    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)
    y = df["is_known_paralog_sl"].astype(int).to_numpy()
    s = df["dependency_dd"].fillna(0).to_numpy()
    obs = auroc(y, s)
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    neg_idx = np.where(y == 0)[0]
    print(f"  frame: {len(y)} entries, {n_pos} positives, {n_neg} unlabeled; "
          f"observed AUROC = {obs:.4f}")

    rng = np.random.default_rng(SEED)
    rows = []
    for pi in PIS:
        k = int(round(pi * n_neg))
        vals = np.empty(N_REP)
        for i in range(N_REP):
            flip = rng.choice(neg_idx, size=k, replace=False)
            y2 = y.copy()
            y2[flip] = 1
            vals[i] = auroc(y2, s)
        # deterministic worst / best case: AUROC is hurt most when the
        # hidden positives sit at the BOTTOM of the ranking (lowest-scoring
        # unlabeled entries flipped positive), and helped when they sit at
        # the top.
        order_neg = neg_idx[np.argsort(-s[neg_idx], kind="mergesort")]
        y_worst = y.copy(); y_worst[order_neg[-k:]] = 1
        auc_worst = auroc(y_worst, s)
        y_best = y.copy(); y_best[order_neg[:k]] = 1
        auc_best = auroc(y_best, s)
        e = n_pos / (n_pos + k)
        en = 0.5 + (obs - 0.5) / e
        rows.append({
            "pi": pi, "k_hidden_positives": k,
            "n_replicates": N_REP, "seed": SEED,
            "observed_auroc": obs,
            "random_flip_mean_auroc": float(vals.mean()),
            "random_flip_sd": float(vals.std()),
            "random_flip_ci95_lo": float(np.percentile(vals, 2.5)),
            "random_flip_ci95_hi": float(np.percentile(vals, 97.5)),
            "random_flip_min": float(vals.min()),
            "worst_case_bottomk_auroc": float(auc_worst),
            "best_case_topk_auroc": float(auc_best),
            "elkan_noto_e": float(e),
            "elkan_noto_corrected_auroc": float(en),
            "frac_replicates_above_0_5": float((vals > 0.5).mean()),
        })
        print(f"  pi={pi:.2f} (k={k}): random flip mean {vals.mean():.4f} "
              f"[{np.percentile(vals, 2.5):.4f}, {np.percentile(vals, 97.5):.4f}], "
              f"worst-case {auc_worst:.4f}, Elkan-Noto {en:.4f}")

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "b3_pu_sensitivity.csv", index=False)
    meta = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "frame": "PRIMARY lineage-level (TableS2): 110 entries, 8 labeled positives, "
                 "102 unlabeled controls; score = signed DD (NaN->0)",
        "design": {
            "random_flip": f"{N_REP} replicates per pi, seed {SEED}; k=round(pi*102) "
                           "uniformly chosen unlabeled entries relabeled positive",
            "adversarial_worst_case": "flip the k LOWEST-scoring unlabeled entries "
                                      "(hidden positives at the bottom of the ranking "
                                      "maximally degrade AUROC)",
            "best_case": "flip the k HIGHEST-scoring unlabeled entries (they were "
                         "already ranked like positives; AUROC increases)",
            "elkan_noto": "e = 8/(8+k) under SCAR; AUROC_corrected = 0.5 + (AUROC_obs-0.5)/e "
                          "(Elkan & Noto 2008; optional)",
        },
        "observed_auroc": float(obs),
        "table": "b3_pu_sensitivity.csv",
        "results": rows,
    }
    (OUT / "b3_pu_sensitivity.json").write_text(json.dumps(meta, indent=2))
    print(f"\n  wrote {OUT}/b3_pu_sensitivity.csv/.json")


if __name__ == "__main__":
    main()
