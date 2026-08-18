#!/usr/bin/env python3
"""
rev2_b15_dws_threshold_stability.py  (Stage-4 revision, item B15)
=================================================================
The production DWS classification (dws_robustness.py, Table S5) labels a pair
HIGH_SELECTIVITY when mean selectivity > 0.15 AND mean signed DWS > 1.0
(and PAN_ESSENTIAL when pan-essential fraction > 0.5, which overrides).
dws_robustness.py already reports flips at selectivity thresholds
0.10/0.15/0.20; this item completes the stability analysis by perturbing
BOTH classification thresholds jointly:

  selectivity threshold s in {0.10, 0.125, 0.15, 0.175, 0.20}   (base 0.15)
  DWS threshold         d in {0.80, 0.90, 1.00, 1.10, 1.20}     (base 1.00)

For every (s, d) cell the HIGH_SELECTIVITY set is recomputed from the
production pair-level means (output/therapeutic_window_all_results.csv,
same aggregation as dws_robustness.py) and compared with the base set by
Jaccard index and retention (fraction of base members retained). Per-pair
distances to both thresholds (margins) identify the borderline members.

Deterministic; no resampling. Fidelity: the recomputed base classification
must reproduce the shipped threshold-0.15 HIGH_SELECTIVITY membership of
dws_robustness.json exactly.
Output: output/revision_stage4/rev2_b15_dws_threshold_stability.{json,csv}
Usage: python rev2_b15_dws_threshold_stability.py   (run from repo root)
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

from config import OUTPUT_DIR  # noqa: E402

SEL_GRID = [0.10, 0.125, 0.15, 0.175, 0.20]
DWS_GRID = [0.80, 0.90, 1.00, 1.10, 1.20]
BASE_SEL, BASE_DWS = 0.15, 1.00
PAN_THR = 0.5


def pair_table(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby(["driver", "paralog"], as_index=False).agg(
        mean_dws=("dws_signed", "mean"),
        mean_sel=("selectivity", "mean"),
        mean_pan=("paralog_pan_essential_frac", "mean"))
    g["pair"] = g["driver"] + "->" + g["paralog"]
    return g


def classify(tab, sel_thr, dws_thr):
    cls = []
    for r in tab.itertuples():
        if r.mean_pan > PAN_THR:
            cls.append("PAN_ESSENTIAL")
        elif r.mean_sel > sel_thr and r.mean_dws > dws_thr:
            cls.append("HIGH_SELECTIVITY")
        elif r.mean_sel > 0:
            cls.append("MODERATE")
        else:
            cls.append("LOW_SELECTIVITY")
    return np.array(cls)


def main():
    df = pd.read_csv(OUTPUT_DIR / "therapeutic_window_all_results.csv")
    with open(OUTPUT_DIR / "dws_robustness.json") as fh:
        shipped = json.load(fh)
    ship_base = None
    for c in shipped["sensitivity"]:
        if c.get("variant") == "HIGH_SELECTIVITY selectivity threshold 0.15":
            ship_base = set(c["high_selectivity_pairs"])
    assert ship_base is not None, "shipped 0.15 threshold entry missing"

    tab = pair_table(df)
    base_cls = classify(tab, BASE_SEL, BASE_DWS)
    base_set = set(tab.loc[base_cls == "HIGH_SELECTIVITY", "pair"])
    assert base_set == ship_base, (
        f"FIDELITY FAIL: recomputed base set {sorted(base_set)} != "
        f"shipped {sorted(ship_base)}")
    print(f"base HIGH_SELECTIVITY (s>{BASE_SEL}, DWS>{BASE_DWS}): "
          f"{len(base_set)} pairs — fidelity OK vs dws_robustness.json")

    # ── (s, d) grid ─────────────────────────────────────────────────
    rows = []
    for s in SEL_GRID:
        for d in DWS_GRID:
            cls = classify(tab, s, d)
            cur = set(tab.loc[cls == "HIGH_SELECTIVITY", "pair"])
            inter = len(cur & base_set)
            union = len(cur | base_set)
            rows.append({
                "selectivity_threshold": s, "dws_threshold": d,
                "n_high_selectivity": len(cur),
                "jaccard_vs_base": inter / union if union else 1.0,
                "retention_of_base": inter / len(base_set) if base_set else 1.0,
                "gained_vs_base": sorted(cur - base_set),
                "lost_vs_base": sorted(base_set - cur),
                "members": sorted(cur),
            })
            if (s, d) != (BASE_SEL, BASE_DWS):
                print(f"  s>{s:<5} DWS>{d:<4} n={len(cur):2d} "
                      f"J={rows[-1]['jaccard_vs_base']:.3f} "
                      f"ret={rows[-1]['retention_of_base']:.3f} "
                      f"lost={rows[-1]['lost_vs_base']} gained={rows[-1]['gained_vs_base']}")

    grid_df = pd.DataFrame([{k: v for k, v in r.items()
                             if k not in ("members", "gained_vs_base", "lost_vs_base")}
                            for r in rows])
    grid_df.to_csv(OUT / "rev2_b15_dws_threshold_stability.csv", index=False)

    # ── per-member margins to the thresholds ───────────────────────
    margins = []
    for pair in sorted(base_set):
        r = tab.loc[tab["pair"] == pair].iloc[0]
        margins.append({
            "pair": pair,
            "mean_selectivity": float(r["mean_sel"]),
            "mean_dws": float(r["mean_dws"]),
            "margin_selectivity": float(r["mean_sel"] - BASE_SEL),
            "margin_dws": float(r["mean_dws"] - BASE_DWS),
        })
    # nearest non-members below the thresholds (would enter on relaxation)
    non = tab.loc[~tab["pair"].isin(base_set) & (tab["mean_pan"] <= PAN_THR)].copy()
    non["dist_below"] = np.maximum(BASE_SEL - non["mean_sel"], 0) \
        + np.maximum(BASE_DWS - non["mean_dws"], 0)
    nearest = non.sort_values("dist_below").head(5)[
        ["pair", "mean_sel", "mean_dws", "dist_below"]]
    nearest = [{"pair": r.pair, "mean_selectivity": float(r.mean_sel),
                "mean_dws": float(r.mean_dws),
                "distance_below_thresholds": float(r.dist_below)}
               for r in nearest.itertuples()]

    out = {
        "method": ("deterministic reclassification of production pair-level means "
                   "(output/therapeutic_window_all_results.csv) on a joint grid of "
                   "the selectivity and DWS thresholds; PAN_ESSENTIAL override "
                   "(pan-essential fraction > 0.5) kept fixed; Jaccard and "
                   "retention vs the base HIGH_SELECTIVITY set"),
        "base": {"selectivity_threshold": BASE_SEL, "dws_threshold": BASE_DWS,
                 "n_high_selectivity": len(base_set), "members": sorted(base_set)},
        "grid": rows,
        "member_margins": margins,
        "nearest_non_members": nearest,
        "fidelity": "recomputed base set == dws_robustness.json 0.15-threshold set",
    }
    out_path = OUT / "rev2_b15_dws_threshold_stability.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nwritten: {out_path}")


if __name__ == "__main__":
    main()
