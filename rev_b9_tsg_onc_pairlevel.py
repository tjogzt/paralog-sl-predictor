#!/usr/bin/env python3
"""
rev_b9_tsg_onc_pairlevel.py  (Stage-4 revision, item B9)
=========================================================
TSG-vs-oncogene mechanism contrast redone at the PAIR level, replacing the
lineage-level AUROC test of Fig. 1b (which compared 9 TSG-driven vs 3
oncogene-driven lineage AUROCs; exact permutation p = 0.127, MW p = 0.145
on the min3 sensitivity frame -- audit_manuscript_numbers.py).

Design:
  * Unit: one row per driver->paralog pair (72 pairs of the primary frame),
    score = mean signed DD across the pair's evaluable gyn3 lineages (the
    0.672 per-pair mean frame of ml_benchmark.py / compute_headline_metrics.py).
  * Group: driver class from output/driver_mutation_rules.csv (TSG = LoF
    rule, ONC = hotspot rule). Pairs whose driver has no TSG/ONC class
    (rule 'ANY', e.g. CDKN2A not in the class map) are excluded and listed.
  * Test: permutation of the pair-level class labels (10,000 shuffles,
    seed 42), statistic = mean(DD | TSG) - mean(DD | ONC); two-sided
    empirical p. Effect sizes: mean/median difference + Cohen's d.
  * Sensitivity: lineage-stratified permutation -- labels shuffled
    independently within each of the 3 lineage strata; statistic = mean of
    the 3 per-lineage TSG-minus-ONC differences (10,000, seed 42).
  * Post-hoc: this pair-level operationalization is declared post hoc, as
    required by the review (R1-M3 / R4-M4/A7).

Note on the parent-specified variant "positive pairs only": the 6 positive
pairs split TSG {ARID1A, FBXW7, BRCA1, BRCA2, STK11} vs ONC {PIK3CA} --
n_ONC = 1, so a positives-only permutation test is NOT COMPUTABLE; it is
recorded as such and the all-pairs test above is the primary pair-level
analysis.

Outputs (output/revision_stage4/):
  b9_tsg_onc_pairlevel.json
  b9_pair_level_dd_by_class.csv

Usage: python rev_b9_tsg_onc_pairlevel.py   (run from repo root)
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

SEED = 42
N_PERM = 10_000


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    sp = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / sp) if sp > 0 else float("nan")


def main():
    print("=" * 72)
    print("  rev B9: TSG vs oncogene pair-level mechanism contrast")
    print("=" * 72)

    df = pd.read_csv(ROOT / "output" / "tables" / "TableS2_FullResults.tsv", sep="\t")
    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)
    rules = pd.read_csv(ROOT / "output" / "driver_mutation_rules.csv")
    cls_of = dict(zip(rules["gene"], rules["driver_class"]))

    # lineage-level entries with class
    df["driver_class"] = df["driver_gene"].map(cls_of)
    excluded = sorted(df.loc[~df["driver_class"].isin(["TSG", "ONC"]), "driver_gene"].unique())

    # pair-level frame (mean signed DD across lineages)
    g = (df[df["driver_class"].isin(["TSG", "ONC"])]
         .groupby(["driver_gene", "paralog_gene"], as_index=False)
         .agg(dd=("dependency_dd", "mean"),
              known=("is_known_paralog_sl", "max"),
              driver_class=("driver_class", "first"),
              n_lineages=("cancer_type", "count")))
    g.to_csv(OUT / "b9_pair_level_dd_by_class.csv", index=False)

    tsg = g.loc[g["driver_class"] == "TSG", "dd"].to_numpy()
    onc = g.loc[g["driver_class"] == "ONC", "dd"].to_numpy()
    n_tsg, n_onc = len(tsg), len(onc)
    obs = tsg.mean() - onc.mean()
    labels = (g["driver_class"] == "TSG").to_numpy()
    scores = g["dd"].to_numpy()

    rng = np.random.default_rng(SEED)
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        lab = rng.permutation(labels)
        null[i] = scores[lab].mean() - scores[~lab].mean()
    p_2s = float((1 + np.sum(np.abs(null) >= abs(obs))) / (1 + N_PERM))
    p_gt = float((1 + np.sum(null >= obs)) / (1 + N_PERM))
    print(f"  pair-level: TSG n={n_tsg} (mean DD {tsg.mean():.4f}) vs "
          f"ONC n={n_onc} (mean DD {onc.mean():.4f}); diff {obs:.4f}; "
          f"perm p2s={p_2s:.4f}")

    # lineage-stratified sensitivity
    sub = df[df["driver_class"].isin(["TSG", "ONC"])].copy()
    lineages = sorted(sub["cancer_type"].unique())
    def strat_stat(d):
        diffs = []
        for lin in lineages:
            s = d[d["cancer_type"] == lin]
            t = s.loc[s["driver_class"] == "TSG", "dependency_dd"]
            o = s.loc[s["driver_class"] == "ONC", "dependency_dd"]
            if len(t) and len(o):
                diffs.append(t.mean() - o.mean())
        return np.mean(diffs) if diffs else np.nan
    obs_strat = strat_stat(sub)
    null_strat = np.empty(N_PERM)
    for i in range(N_PERM):
        parts = []
        for lin in lineages:
            s = sub[sub["cancer_type"] == lin].copy()
            s["driver_class"] = rng.permutation(s["driver_class"].to_numpy())
            parts.append(s)
        null_strat[i] = strat_stat(pd.concat(parts))
    ok = ~np.isnan(null_strat)
    p_strat_2s = float((1 + np.sum(np.abs(null_strat[ok]) >= abs(obs_strat))) / (1 + ok.sum()))
    print(f"  lineage-stratified: combined diff {obs_strat:.4f}; stratified perm p2s={p_strat_2s:.4f}")

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "post_hoc_declaration": "This pair-level operationalization of the TSG/oncogene "
                                "contrast was defined after seeing the data (review "
                                "R1-M3/R4-M4); it replaces the lineage-level AUROC test.",
        "unit_and_score": "72 driver->paralog pairs; mean signed DD across evaluable gyn3 "
                          "lineages (the per-pair mean frame)",
        "class_source": "output/driver_mutation_rules.csv (TSG=LikelyLoF rule, ONC=Hotspot rule)",
        "excluded_drivers_no_class": excluded,
        "n_pairs_tsg": int(n_tsg), "n_pairs_onc": int(n_onc),
        "mean_dd_tsg": float(tsg.mean()), "mean_dd_onc": float(onc.mean()),
        "median_dd_tsg": float(np.median(tsg)), "median_dd_onc": float(np.median(onc)),
        "mean_diff_tsg_minus_onc": float(obs),
        "median_diff_tsg_minus_onc": float(np.median(tsg) - np.median(onc)),
        "cohens_d": cohens_d(tsg, onc),
        "permutation_test": {
            "n_permutations": N_PERM, "seed": SEED,
            "scheme": "shuffle pair-level class labels; stat = mean(DD|TSG)-mean(DD|ONC)",
            "p_two_sided": p_2s, "p_one_sided_greater": p_gt,
            "null_mean": float(null.mean()), "null_std": float(null.std())},
        "lineage_stratified_sensitivity": {
            "combined_stat_mean_of_lineage_diffs": float(obs_strat),
            "n_permutations": N_PERM, "seed": SEED,
            "scheme": "class labels shuffled independently within each lineage stratum",
            "p_two_sided": p_strat_2s,
            "per_lineage_diffs": {lin: float(
                sub.loc[(sub["cancer_type"] == lin) & (sub["driver_class"] == "TSG"),
                        "dependency_dd"].mean()
                - sub.loc[(sub["cancer_type"] == lin) & (sub["driver_class"] == "ONC"),
                          "dependency_dd"].mean()) for lin in lineages}},
        "positives_only_variant": {
            "status": "NOT COMPUTABLE",
            "reason": "the 6 positive pairs split TSG {ARID1A->ARID1B, FBXW7->FBXW2, "
                      "BRCA1->BRCA2, BRCA2->BRCA1, STK11->SIK1} vs ONC {PIK3CA->PIK3CB}: "
                      "n_ONC = 1 positive pair; a two-group permutation test on positive "
                      "pairs only is undefined.",
        },
        "reference_lineage_level_values": {
            "source": "audit_manuscript_numbers.py / solid_tumor_summary_min3.csv",
            "tsg_mean_auroc_n9": 0.649, "onc_mean_auroc_n3": 0.834,
            "exact_permutation_p": 0.127, "exact_mannwhitney_p": 0.145,
            "note": "existing Fig. 1b lineage-level contrast retained as scatter without "
                    "test annotation per the review",
        },
    }
    (OUT / "b9_tsg_onc_pairlevel.json").write_text(json.dumps(out, indent=2))
    print(f"\n  wrote {OUT}/b9_tsg_onc_pairlevel.json + b9_pair_level_dd_by_class.csv")


if __name__ == "__main__":
    main()
