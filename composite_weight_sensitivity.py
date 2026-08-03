#!/usr/bin/env python3
"""
composite_weight_sensitivity.py
===============================
Dirichlet weight-perturbation sensitivity analysis for the manuscript's
composite paralog-compensation score:

    composite = 0.50*minmax(PCS) + 0.20*minmax(|DD|)
              + 0.15*minmax(-log10 q) + 0.15*minmax(mut_freq)

on the per-pair MEAN frame (72 unique pairs, 6 positives): each component
is first aggregated across a pair's lineage entries by the mean
(build_pair_frame-style; q_value aggregated by the mean), then min-max
normalized across the 72 pairs, then combined with the weights above —
the pair-frame analogue of pcs.py lines 266–273. The base-weight AUROC on
this frame reproduces 0.828 (~= the 0.831 manuscript/ml_benchmark value,
which instead averages the entry-level composite_score column; both are
reported for provenance).

Perturbation: 200 weight vectors sampled from
Dirichlet(concentration = base_weights x 20 = [10, 4, 3, 3]), seed 42;
the composite AUROC is recomputed on the same fixed 72-pair frame for
each draw (labels and components never change — only the weights).

Output: output/composite_weight_sensitivity.json
Usage:  python3 composite_weight_sensitivity.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compute_headline_metrics import auroc  # noqa: E402

ROOT = Path(__file__).resolve().parent
TABLES2 = ROOT / "output" / "tables" / "TableS2_FullResults.tsv"
JSON_OUT = ROOT / "output" / "composite_weight_sensitivity.json"

BASE_WEIGHTS = np.array([0.50, 0.20, 0.15, 0.15])  # PCS, |DD|, -log10 q, mut_freq
CONCENTRATION_MULT = 20
N_SAMPLES = 200
SEED = 42


def minmax(s: pd.Series) -> pd.Series:
    """Min-max normalize to [0, 1]; 0 if constant (same rule as pcs.py)."""
    if s.max() == s.min():
        return pd.Series(0.0, index=s.index)
    return (s - s.min()) / (s.max() - s.min())


def main():
    if not TABLES2.exists():
        sys.exit(f"ERROR: {TABLES2} not found — run main.py + tables.py first")

    df = pd.read_csv(TABLES2, sep="\t")
    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)

    # Per-pair MEAN frame with the four composite ingredients
    g = df.groupby(["driver_gene", "paralog_gene"], as_index=False).agg(
        known=("is_known_paralog_sl", "max"),
        pcs=("pcs", "mean"),
        dd_abs=("dependency_dd", lambda s: s.abs().mean()),
        q_value=("q_value", "mean"),
        mut_freq=("mutation_frequency", "mean"),
        composite_col=("composite_score", "mean"),
    )
    y = g["known"].astype(int).values

    comps = np.column_stack([
        minmax(g["pcs"]).values,
        minmax(g["dd_abs"]).values,
        minmax(-np.log10(g["q_value"] + 1e-10)).values,
        minmax(g["mut_freq"]).values,
    ])

    base_auroc = auroc(y, comps @ BASE_WEIGHTS)
    col_auroc = auroc(y, g["composite_col"].fillna(0).values)
    print(f"Universe: {len(g)} pairs, {int(y.sum())} positives")
    print(f"Base-weight composite AUROC (pair-frame minmax): {base_auroc:.4f}")
    print(f"Reference: mean of entry-level composite column: {col_auroc:.4f}")

    # Dirichlet perturbation around the base weights
    rng = np.random.default_rng(SEED)
    conc = BASE_WEIGHTS * CONCENTRATION_MULT
    draws = rng.dirichlet(conc, size=N_SAMPLES)
    aucs = np.array([auroc(y, comps @ w) for w in draws])

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "composite_weight_sensitivity.py",
        "source": str(TABLES2.relative_to(ROOT)),
        "frame": {
            "description": "per-pair MEAN frame (72 unique pairs, 6 positives); components "
                           "aggregated across lineage entries by mean, min-max normalized "
                           "across the 72 pairs, then weighted",
            "n_pairs": int(len(g)),
            "n_positives": int(y.sum()),
        },
        "composite_formula": "w0*minmax(PCS) + w1*minmax(|DD|) + w2*minmax(-log10(q+1e-10)) "
                             "+ w3*minmax(mut_freq); component order [PCS, |DD|, -log10 q, mut_freq]",
        "base_weights": BASE_WEIGHTS.tolist(),
        "perturbation": {
            "distribution": "Dirichlet",
            "concentration": (BASE_WEIGHTS * CONCENTRATION_MULT).tolist(),
            "concentration_rule": "base_weights x 20",
            "n_samples": N_SAMPLES,
            "seed": SEED,
        },
        "base_weight_auroc": float(base_auroc),
        "reference_auroc_composite_column_mean": float(col_auroc),
        "reference_note": "0.8308 = ml_benchmark.py composite_alone (mean of the entry-level "
                          "composite_score column; manuscript claim 0.831)",
        "auroc_distribution": {
            "min": float(aucs.min()),
            "p25": float(np.percentile(aucs, 25)),
            "p50": float(np.percentile(aucs, 50)),
            "p75": float(np.percentile(aucs, 75)),
            "max": float(aucs.max()),
            "mean": float(aucs.mean()),
            "std": float(aucs.std()),
        },
        "auroc_samples": [round(float(a), 6) for a in aucs],
    }

    JSON_OUT.write_text(json.dumps(out, indent=2, allow_nan=False, default=str))
    print(f"\nDirichlet AUROC: min={aucs.min():.4f} p25={np.percentile(aucs, 25):.4f} "
          f"p50={np.percentile(aucs, 50):.4f} p75={np.percentile(aucs, 75):.4f} max={aucs.max():.4f}")
    print(f"Wrote {JSON_OUT}")


if __name__ == "__main__":
    main()
