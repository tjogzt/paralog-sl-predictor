#!/usr/bin/env python3
"""
direction_audit.py (C6)
========================
Audit the use of |DD| as the ranking score across all evaluation frames.

DD is defined (manuscript Eq. 1) as mean(Chronos|WT) - mean(Chronos|MUT);
positive DD means increased paralog dependency after driver mutation.
Scoring by |DD| additionally promotes REVERSE effects (DD < 0: dependency
decreases after mutation, e.g. oncogene-addiction displacement). This audit
quantifies, for every evaluation frame used in the manuscript:

  * AUROC with |DD| (primary metric), signed DD, and -DD;
  * how many gold-standard positives have DD < 0 ("reverse-signal
    positives") and how many of those rank in the top |DD| quartile
    (reverse-signal contamination of the primary metric).

Frames audited:
  1. pair-level primary set (TableS2, gyn3) — overall and Tier A only;
  2. each pan-cancer lineage (solid_*_results.csv);
  3. MSI subgroups (msi_*_results.csv);
  4. mutation-type frames (muttype_*_results.csv).

Output: output/direction_audit.json (+ console summary).
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
TABLES2 = OUT / "tables" / "TableS2_FullResults.tsv"
JSON_OUT = OUT / "direction_audit.json"

# Directional Tier A gold pairs (compute_headline_metrics.py)
TIER_A = {("SMARCA4", "SMARCA2"), ("ARID1A", "ARID1B")}


def auroc(labels, scores):
    """Rank-based AUROC with average ranks for ties; NaN if a class absent."""
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    ok = ~np.isnan(scores)
    labels, scores = labels[ok], scores[ok]
    if labels.sum() == 0 or labels.sum() == len(labels):
        return np.nan
    order = scores.argsort()
    ranks = np.empty(len(scores))
    ranks[order] = np.arange(1, len(scores) + 1)
    # average ranks for ties
    _, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    csum = np.cumsum(counts)
    avg = (csum - counts + 1 + csum) / 2.0
    ranks = avg[inv]
    n_pos = labels.sum()
    n_neg = len(labels) - n_pos
    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def audit_frame(df, name, restrict_tier_a=False):
    df = df.dropna(subset=["dependency_dd"])
    if restrict_tier_a:
        mask = df.apply(lambda r: (r["driver_gene"], r["paralog_gene"]) in TIER_A, axis=1)
        df = df[mask | ~df["is_known_paralog_sl"]]
    labels = df["is_known_paralog_sl"].astype(int).values
    dd = df["dependency_dd"].values
    n_pos = int(labels.sum())
    res = {
        "frame": name,
        "n_pairs": int(len(df)),
        "n_pos": n_pos,
        "auroc_abs": round(auroc(labels, np.abs(dd)), 4),
        "auroc_signed": round(auroc(labels, dd), 4),
        "auroc_negated": round(auroc(labels, -dd), 4),
    }
    if n_pos:
        pos = df[df["is_known_paralog_sl"]]
        rev = pos[pos["dependency_dd"] < 0]
        res["n_pos_dd_negative"] = int(len(rev))
        res["pos_dd_negative_pairs"] = [
            f"{r.driver_gene}->{r.paralog_gene} (DD={r.dependency_dd:+.3f})"
            for r in rev.itertuples()
        ]
        if len(df) >= 4:
            thr = df["dependency_dd"].abs().quantile(0.75)
            res["n_reverse_in_top_abs_quartile"] = int((rev["dependency_dd"].abs() >= thr).sum())
    return res


def main():
    report = {"frames": []}

    ts2 = pd.read_csv(TABLES2, sep="\t")
    report["frames"].append(audit_frame(ts2, "pair_level_primary_gyn3"))
    report["frames"].append(audit_frame(ts2, "pair_level_primary_gyn3_TIER_A",
                                        restrict_tier_a=True))

    for f in sorted(OUT.glob("solid_*_results.csv")):
        df = pd.read_csv(f)
        if df["is_known_paralog_sl"].sum() >= 1:
            report["frames"].append(audit_frame(
                df, f"lineage:{f.stem.replace('solid_', '').replace('_results', '')}"))

    for f in sorted(OUT.glob("msi_*_results.csv")):
        df = pd.read_csv(f)
        if "is_known_paralog_sl" in df.columns and df["is_known_paralog_sl"].sum() >= 1:
            report["frames"].append(audit_frame(df, f"msi:{f.stem}"))

    for f in sorted(OUT.glob("muttype_*_results.csv")):
        df = pd.read_csv(f)
        if "is_known_paralog_sl" in df.columns and df["is_known_paralog_sl"].sum() >= 1:
            report["frames"].append(audit_frame(df, f"muttype:{f.stem}"))

    # ── Console summary ──
    print(f"{'frame':40s} {'n':>5s} {'pos':>4s} {'|DD|':>7s} {'signed':>7s} "
          f"{'-DD':>7s} {'rev+':>5s} {'revTopQ':>8s}")
    for r in report["frames"]:
        print(f"{r['frame']:40s} {r['n_pairs']:5d} {r['n_pos']:4d} "
              f"{r['auroc_abs']:7.3f} {r['auroc_signed']:7.3f} "
              f"{r['auroc_negated']:7.3f} {r.get('n_pos_dd_negative', 0):5d} "
              f"{r.get('n_reverse_in_top_abs_quartile', 0):8d}")

    JSON_OUT.write_text(json.dumps(report, indent=2, allow_nan=False))
    print(f"\nWrote {JSON_OUT}")


if __name__ == "__main__":
    main()
