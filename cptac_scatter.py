#!/usr/bin/env python3
"""
cptac_scatter.py
================
Fig S4 inventory + full per-cohort scatter rendering.

INVENTORY FINDING (R_figS4_cptac.R): the figure script reads per-sample
protein abundances from data/cptac_cache/UCEC_protein_data.json (a dict of
gene -> {sample_id: log2 abundance}) and draws three UCEC scatter panels
(EP300/CREBBP, PIK3CA/PIK3CB, PIK3R1/CRKL). Inspection of
data/cptac_cache/*.json shows all 7 cohort caches (BRCA, COAD, GBM, LUAD,
LUSC, PDAC, UCEC) contain per-sample values (34-45 genes x 88-140
samples), so NO network fetching is needed. The ccrcc/hnscc/lscc cohorts
have correlation CSVs but NO per-sample cache and are skipped.

This script renders, for every driver x paralog pair of the CPTAC
correlation universe (30 pairs, output/cptac_*_correlations.csv) and
every cached cohort, a scatter PNG (driver vs paralog log2 protein
abundance, OLS regression line, Pearson r + p in the title) whenever
>=10 samples have both genes measured.

Outputs:
  output/cptac_scatter/{cohort}_{driver}_{paralog}.png
  output/cptac_scatter/manifest.csv   (per-plot n, Pearson r, p — numeric artifact)
Usage: python3 cptac_scatter.py   (run from repo root)
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "data" / "cptac_cache"
PAIR_SOURCE = ROOT / "output" / "cptac_ucec_correlations.csv"
OUT_DIR = ROOT / "output" / "cptac_scatter"
MIN_PAIRED = 10

# 7 cohorts with per-sample caches (ccrcc/hnscc/lscc have correlation CSVs
# only, no per-sample cache -> not renderable offline)
COHORTS = ["BRCA", "COAD", "GBM", "LUAD", "LUSC", "PDAC", "UCEC"]


def main():
    pairs = pd.read_csv(PAIR_SOURCE)[["gene_a", "gene_b"]].drop_duplicates()
    pair_list = [tuple(r) for r in pairs.itertuples(index=False)]
    print(f"Pairs: {len(pair_list)} (from {PAIR_SOURCE.name})")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    n_plotted = 0
    for cohort in COHORTS:
        f = CACHE / f"{cohort}_protein_data.json"
        if not f.exists():
            rows.append({"cohort": cohort, "driver": None, "paralog": None,
                         "status": "missing_cache", "n": 0, "r": np.nan, "p": np.nan,
                         "file": None})
            continue
        d = json.loads(f.read_text())
        for a, b in pair_list:
            tag = f"{cohort}_{a}_{b}"
            if a not in d or b not in d:
                rows.append({"cohort": cohort, "driver": a, "paralog": b,
                             "status": "gene_not_in_cache", "n": 0,
                             "r": np.nan, "p": np.nan, "file": None})
                continue
            xa = pd.Series(d[a], dtype=float)
            xb = pd.Series(d[b], dtype=float)
            both = pd.concat([xa, xb], axis=1, keys=["x", "y"]).dropna()
            n = len(both)
            if n < MIN_PAIRED:
                rows.append({"cohort": cohort, "driver": a, "paralog": b,
                             "status": f"fewer_than_{MIN_PAIRED}_paired_samples",
                             "n": n, "r": np.nan, "p": np.nan, "file": None})
                continue
            r, p = stats.pearsonr(both["x"], both["y"])
            slope, intercept = np.polyfit(both["x"], both["y"], 1)

            fig, ax = plt.subplots(figsize=(3.2, 3.0), dpi=150)
            ax.scatter(both["x"], both["y"], s=8, alpha=0.45,
                       color="#2171B5", edgecolors="none")
            xs = np.linspace(both["x"].min(), both["x"].max(), 50)
            ax.plot(xs, slope * xs + intercept, color="#CB181D", lw=1.2)
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "(ns)"
            ax.set_title(f"{cohort}: {a}/{b} | r={r:.3f} {sig}", fontsize=8,
                         fontweight="bold")
            ax.set_xlabel(f"{a} log2 abundance", fontsize=7)
            ax.set_ylabel(f"{b} log2 abundance", fontsize=7)
            ax.tick_params(labelsize=6)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            fig.tight_layout()
            png = OUT_DIR / f"{tag}.png"
            fig.savefig(png)
            plt.close(fig)
            n_plotted += 1
            rows.append({"cohort": cohort, "driver": a, "paralog": b,
                         "status": "plotted", "n": n, "r": r, "p": p,
                         "file": png.name})

    man = pd.DataFrame(rows)
    man.to_csv(OUT_DIR / "manifest.csv", index=False)
    summary = {
        "inventory": "R_figS4_cptac.R reads data/cptac_cache/UCEC_protein_data.json "
                     "(per-sample dict gene -> {sample: log2 abundance}); all 7 cached "
                     "cohorts have per-sample values -> rendered offline, no fetching",
        "cohorts_with_cache": COHORTS,
        "cohorts_without_cache": ["CCRCC", "HNSCC", "LSCC"],
        "pair_source": str(PAIR_SOURCE.relative_to(ROOT)),
        "n_pairs": len(pair_list),
        "min_paired_samples": MIN_PAIRED,
        "n_plots": n_plotted,
        "n_skipped": int((man["status"] != "plotted").sum()),
        "skipped_by_reason": man.loc[man["status"] != "plotted", "status"]
                              .value_counts().to_dict(),
        "manifest": str((OUT_DIR / "manifest.csv").relative_to(ROOT)),
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Plots written: {n_plotted} -> {OUT_DIR}")
    print(f"Skipped: {summary['n_skipped']} {summary['skipped_by_reason']}")
    print(f"Manifest: {OUT_DIR/'manifest.csv'}")


if __name__ == "__main__":
    main()
