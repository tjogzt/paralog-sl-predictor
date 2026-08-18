#!/usr/bin/env python3
"""
rev_b4_necessity_increment.py  (Stage-4 revision, item B4)
===========================================================
Conditional-increment analysis of `necessity` (pan-lineage mean essentiality
of the paralog, -mean Chronos) vs signed DD on the PRIMARY lineage-level
frame (TableS2, 110 entries, 8 positives):

  (a) necessity-decile stratification: per-decile n entries / n positives /
      mean and median signed DD by class. Per-decile AUROC is reported only
      where both classes are present (with 8 positives most deciles are not
      estimable -- recorded explicitly).
  (b) Head-to-head AUROC: necessity alone vs signed DD alone vs the
      {necessity + signed DD} combination (two pre-specified combos:
      equal-weight z-score sum; logistic regression apparent fit, flagged as
      in-sample). Uncertainty: paired bootstrap (10,000 entry resamples,
      same draws for all scores, seed 42) with delta-AUROC 95% CI vs DD.
  (c) Composite-score ablation WITHOUT necessity: PCS' = max(delta_expr, 0)
      (necessity factor removed); composite' = 0.50*minmax(PCS') +
      0.20*minmax(|DD|) + 0.15*minmax(-log10 q) + 0.15*minmax(mut_freq),
      entry-level global min-max exactly as pcs.py run_full_analysis.
      Evaluated on the per-pair mean frame (72 pairs; shipped composite =
      0.831) and the lineage-level frame. Fidelity: rebuilt shipped
      composite must reproduce TableS2 composite_score (1e-9) and the 0.831
      per-pair AUROC.

Outputs (output/revision_stage4/):
  b4_necessity_decile_stratification.csv
  b4_necessity_vs_dd.json
  b4_composite_ablation.json

Usage: python rev_b4_necessity_increment.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

from compute_headline_metrics import auroc  # noqa: E402

SEED = 42
N_BOOT = 10_000
W = {"pcs": 0.50, "dd": 0.20, "q": 0.15, "freq": 0.15}  # pcs.py composite weights
GYN3 = ["Ovarian", "Endometrial", "Cervical"]


def minmax(s: pd.Series) -> pd.Series:
    if s.max() == s.min():
        return pd.Series(0.0, index=s.index)
    return (s - s.min()) / (s.max() - s.min())


def main():
    print("=" * 72)
    print("  rev B4: necessity conditional increment")
    print("=" * 72)

    df = pd.read_csv(ROOT / "output" / "tables" / "TableS2_FullResults.tsv", sep="\t")
    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)
    y = df["is_known_paralog_sl"].astype(int).to_numpy()
    dd = df["dependency_dd"].fillna(0).to_numpy()
    nec = df["necessity"].fillna(0).to_numpy()

    # ── (a) necessity-decile stratification ───────────────────────
    df["_nec_decile"] = pd.qcut(df["necessity"], 10, labels=False, duplicates="drop")
    rows = []
    for d in sorted(df["_nec_decile"].unique()):
        sub = df[df["_nec_decile"] == d]
        yp = sub[sub["is_known_paralog_sl"]]
        yn = sub[~sub["is_known_paralog_sl"]]
        row = {
            "nec_decile": int(d) + 1,
            "necessity_range": f"[{sub['necessity'].min():.4f}, {sub['necessity'].max():.4f}]",
            "n_entries": int(len(sub)), "n_positives": int(len(yp)),
            "mean_dd_pos": float(yp["dependency_dd"].mean()) if len(yp) else None,
            "mean_dd_neg": float(yn["dependency_dd"].mean()) if len(yn) else None,
            "median_dd_pos": float(yp["dependency_dd"].median()) if len(yp) else None,
            "median_dd_neg": float(yn["dependency_dd"].median()) if len(yn) else None,
        }
        row["dd_separation_pos_minus_neg"] = (
            row["mean_dd_pos"] - row["mean_dd_neg"]
            if row["mean_dd_pos"] is not None and row["mean_dd_neg"] is not None else None)
        if len(yp) >= 1 and len(yn) >= 1:
            row["auroc_signed_dd_within_decile"] = float(
                auroc(sub["is_known_paralog_sl"].astype(int).to_numpy(),
                      sub["dependency_dd"].fillna(0).to_numpy()))
        else:
            row["auroc_signed_dd_within_decile"] = None
        rows.append(row)
    dec = pd.DataFrame(rows)
    dec.to_csv(OUT / "b4_necessity_decile_stratification.csv", index=False)
    n_estimable = dec["auroc_signed_dd_within_decile"].notna().sum()
    print(f"[a] deciles written; within-decile AUROC estimable in {n_estimable}/"
          f"{len(dec)} deciles (8 positives total)")

    # ── (b) necessity vs DD vs combination ────────────────────────
    z = lambda v: (v - v.mean()) / (v.std(ddof=1) if v.std(ddof=1) > 0 else 1.0)
    combo_z = z(nec) + z(dd)
    X = np.column_stack([z(nec), z(dd)])
    lr = LogisticRegression(max_iter=10000).fit(X, y)
    combo_lr = lr.predict_proba(X)[:, 1]
    scores = {"necessity": nec, "signed_dd": dd,
              "combo_equal_weight_z": combo_z, "combo_lr_apparent": combo_lr}
    aucs = {k: auroc(y, v) for k, v in scores.items()}
    print(f"[b] AUROC: necessity={aucs['necessity']:.4f} dd={aucs['signed_dd']:.4f} "
          f"z-combo={aucs['combo_equal_weight_z']:.4f} lr-apparent={aucs['combo_lr_apparent']:.4f}")

    rng = np.random.default_rng(SEED)
    n = len(y)
    deltas = {k: [] for k in scores if k != "signed_dd"}
    boots = {k: [] for k in scores}
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        yb = y[idx]
        if yb.sum() == 0 or yb.sum() == n:
            continue
        a = {k: auroc(yb, v[idx]) for k, v in scores.items()}
        for k in boots:
            boots[k].append(a[k])
        for k in deltas:
            deltas[k].append(a[k] - a["signed_dd"])
    pb = {}
    for k, v in deltas.items():
        arr = np.asarray(v)
        pb[f"{k}_minus_signed_dd"] = {
            "mean_delta": float(arr.mean()),
            "ci95": [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))],
            "frac_above_0": float((arr > 0).mean())}
    out_b = {
        "frame": "PRIMARY lineage-level (TableS2): 110 entries, 8 positives",
        "auroc": {k: float(v) for k, v in aucs.items()},
        "auroc_ci95_paired_bootstrap": {
            k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
            for k, v in boots.items()},
        "paired_bootstrap_delta_vs_signed_dd": pb,
        "combo_definitions": {
            "combo_equal_weight_z": "z(necessity) + z(signed DD), pre-specified equal weights",
            "combo_lr_apparent": "logistic regression label ~ z(necessity)+z(DD), apparent "
                                 "(in-sample) score; optimistic by construction",
        },
        "n_boot": N_BOOT, "seed": SEED,
        "decile_table": "b4_necessity_decile_stratification.csv",
        "within_decile_auroc_estimable_deciles": f"{n_estimable}/{len(dec)}",
    }
    (OUT / "b4_necessity_vs_dd.json").write_text(json.dumps(out_b, indent=2))

    # ── (c) composite ablation without necessity ──────────────────
    # The shipped composite_score was min-max normalized GLOBALLY over all 197
    # candidate entries (incl. Lung/Breast) in pcs.py run_full_analysis;
    # TableS2 is the gyn3 subset of that run. Rebuild on the full candidate
    # file, verify fidelity, then subset to gyn3.
    cand = pd.read_csv(ROOT / "output" / "paralog_sl_candidates.csv")
    cand["is_known_paralog_sl"] = cand["is_known_paralog_sl"].astype(bool)

    def build_composite(d, pcs_col):
        return (W["pcs"] * minmax(d[pcs_col]) + W["dd"] * minmax(d["dependency_dd"].abs())
                + W["q"] * minmax(-np.log10(d["q_value"].clip(lower=1e-10)))
                + W["freq"] * minmax(d["mutation_frequency"]))

    rebuilt = build_composite(cand, "pcs")
    if (rebuilt - cand["composite_score"]).abs().max() > 1e-9:
        raise RuntimeError("FIDELITY FAIL: rebuilt composite != candidates composite_score")
    cand["_pcs_no_nec"] = cand["delta_expression"].clip(lower=0.0)  # necessity factor removed
    cand["_composite_no_nec"] = build_composite(cand, "_pcs_no_nec")

    df = cand[cand["cancer_type"].isin(GYN3)].copy().reset_index(drop=True)
    # fidelity: TableS2 composite matches the candidates subset
    ts2 = pd.read_csv(ROOT / "output" / "tables" / "TableS2_FullResults.tsv", sep="\t")
    merged = df.merge(ts2, on=["driver_gene", "paralog_gene", "cancer_type"],
                      suffixes=("", "_ts2"))
    if len(merged) != len(ts2) or (merged["composite_score"] - merged["composite_score_ts2"]).abs().max() > 1e-9:
        raise RuntimeError("FIDELITY FAIL: candidates gyn3 subset != TableS2")
    y = df["is_known_paralog_sl"].astype(int).to_numpy()

    g = df.groupby(["driver_gene", "paralog_gene"], as_index=False).agg(
        known=("is_known_paralog_sl", "max"),
        comp=("composite_score", "mean"),
        comp_no_nec=("_composite_no_nec", "mean"),
        dd=("dependency_dd", "mean"),
        nec=("necessity", "mean"))
    yg = g["known"].astype(int).to_numpy()
    per_pair = {
        "composite_shipped": float(auroc(yg, g["comp"].fillna(0).to_numpy())),
        "composite_without_necessity": float(auroc(yg, g["comp_no_nec"].fillna(0).to_numpy())),
        "signed_dd_alone": float(auroc(yg, g["dd"].fillna(0).to_numpy())),
        "necessity_alone": float(auroc(yg, g["nec"].fillna(0).to_numpy())),
    }
    if abs(per_pair["composite_shipped"] - 0.8308) > 5e-3:
        raise RuntimeError(f"FIDELITY FAIL: per-pair composite {per_pair['composite_shipped']} != 0.831")
    lineage_level = {
        "composite_shipped": float(auroc(y, df["composite_score"].fillna(0).to_numpy())),
        "composite_without_necessity": float(auroc(y, df["_composite_no_nec"].fillna(0).to_numpy())),
    }
    # paired bootstrap on the per-pair frame for the ablation delta
    rng = np.random.default_rng(SEED)
    npair = len(yg)
    d_abl = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, npair, npair)
        yb = yg[idx]
        if yb.sum() == 0 or yb.sum() == npair:
            continue
        d_abl.append(auroc(yb, g["comp_no_nec"].fillna(0).to_numpy()[idx])
                     - auroc(yb, g["comp"].fillna(0).to_numpy()[idx]))
    d_abl = np.asarray(d_abl)
    out_c = {
        "ablation": "PCS' = max(delta_expression, 0) (necessity factor removed); composite "
                    "rebuilt with the same 0.50/0.20/0.15/0.15 weights and entry-level "
                    "global min-max normalization as pcs.py",
        "fidelity": "rebuilt shipped composite == TableS2 composite_score (max abs diff < 1e-9); "
                    "per-pair shipped composite AUROC reproduces 0.831",
        "per_pair_mean_frame_72_pairs": {**per_pair,
            "n_pairs": int(npair), "n_positives": int(yg.sum()),
            "ablation_delta": per_pair["composite_without_necessity"] - per_pair["composite_shipped"],
            "ablation_delta_paired_bootstrap": {
                "mean": float(d_abl.mean()),
                "ci95": [float(np.percentile(d_abl, 2.5)), float(np.percentile(d_abl, 97.5))],
                "frac_above_0": float((d_abl > 0).mean()), "n_boot": N_BOOT, "seed": SEED}},
        "lineage_level_frame_110_entries": {**lineage_level,
            "ablation_delta": lineage_level["composite_without_necessity"]
                                - lineage_level["composite_shipped"]},
    }
    (OUT / "b4_composite_ablation.json").write_text(json.dumps(out_c, indent=2))
    print(f"[c] composite ablation: per-pair shipped {per_pair['composite_shipped']:.4f} "
          f"-> no-necessity {per_pair['composite_without_necessity']:.4f} "
          f"(delta {out_c['per_pair_mean_frame_72_pairs']['ablation_delta']:+.4f})")
    print("\nDone: b4_necessity_decile_stratification.csv, b4_necessity_vs_dd.json, "
          "b4_composite_ablation.json")


if __name__ == "__main__":
    main()
