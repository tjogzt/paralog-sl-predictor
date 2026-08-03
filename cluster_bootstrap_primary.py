#!/usr/bin/env python3
"""
cluster_bootstrap_primary.py
============================
Cluster-aware uncertainty quantification for the PRIMARY lineage-level
DD benchmark (manuscript AUROC = 0.629), plus the never-permutation-tested
per-pair MEAN aggregation (AUROC = 0.672).

Frame (same evaluation universe as compute_headline_metrics.py /
validation_report.json):
  output/tables/TableS2_FullResults.tsv — 110 driver×paralog×lineage
  entries in the 3 gynecological lineages (Ovarian/Endometrial/Cervical),
  72 unique pairs, 8 positive entries, 6 positive pairs; labels in
  `is_known_paralog_sl`, score = signed `dependency_dd`.

Entries are clustered within pairs (1–2 lineage entries per pair), so the
naive entry-level bootstrap used elsewhere in the pipeline understates
uncertainty. This script:

  (a) pair-clustered bootstrap 95% CI of the primary AUROC: resample the
      72 pairs with replacement (10,000 iterations, seed 42), carry ALL
      lineage entries of each resampled pair, recompute the rank-based
      AUROC over the resulting entry multiset; iterations whose resample
      contains <2 positives are skipped and their fraction reported.
  (b) pair-level permutation test of the same primary statistic: permute
      the pair-level labels (6 positive / 66 control), carry labels down
      to entries, 10,000 permutations (seed 42); empirical
      p = (1 + #{null >= observed}) / (1 + n_perm).
  (c) effective sample size: average cluster size m, ANOVA-style ICC of
      dependency_dd by pair (one-way random effects, m0 correction for
      unequal cluster sizes), NEFF = 110 / (1 + (m-1)*ICC).
  (d) permutation test (10,000, seed 42) for the MEAN-aggregated per-pair
      DD frame (72 pairs, mean signed DD across lineages; AUROC = 0.672
      per ml_benchmark.py / compute_headline_metrics.py).
  (e) AUPRC (average precision) of signed DD on the primary lineage-level
      frame + pair-clustered bootstrap 95% CI (same resample scheme as
      (a), same skipped iterations); baseline prevalence 8/110.
  (f) precision@k on the primary frame: fraction of known positives among
      the top-10 and top-20 entries ranked by signed DD.
  (g) driver-block bootstrap 95% CI of the primary AUROC: resample the
      driver genes with replacement (10,000 iterations, seed 42), carrying
      ALL pairs and lineage entries under each resampled driver. Controls
      and curated pairs share driver genes, so multiple pairs under one
      driver are non-independent; this block scheme accounts for that
      coarser clustering level. Same skip rule as (a).

AUROC/AUPRC definitions are imported from compute_headline_metrics.py so
they are byte-identical to the manuscript's headline-metrics code
(rank-based Mann-Whitney AUROC with average ties; average precision =
mean of precision@k over the ranks of the positives).

Output: output/cluster_bootstrap_primary.json
Usage:  python3 cluster_bootstrap_primary.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compute_headline_metrics import auroc, auprc  # noqa: E402

ROOT = Path(__file__).resolve().parent
TABLES2 = ROOT / "output" / "tables" / "TableS2_FullResults.tsv"
JSON_OUT = ROOT / "output" / "cluster_bootstrap_primary.json"

SEED = 42
N_BOOT = 10_000
N_PERM = 10_000


def main():
    if not TABLES2.exists():
        sys.exit(f"ERROR: {TABLES2} not found — run main.py + tables.py first")

    df = pd.read_csv(TABLES2, sep="\t")
    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)

    scores = df["dependency_dd"].fillna(0).values.astype(float)
    labels = df["is_known_paralog_sl"].astype(int).values
    n_entries = len(df)
    n_pos_entries = int(labels.sum())

    # Pair structure
    pair_keys = list(zip(df["driver_gene"], df["paralog_gene"]))
    pair_ids = sorted(set(pair_keys))
    pair_to_int = {p: i for i, p in enumerate(pair_ids)}
    entry_pair = np.array([pair_to_int[k] for k in pair_keys])
    n_pairs = len(pair_ids)
    pair_label = np.array(
        [int(df.loc[entry_pair == i, "is_known_paralog_sl"].max()) for i in range(n_pairs)]
    )
    n_pos_pairs = int(pair_label.sum())

    # Per-pair entry index lists (cluster membership)
    pair_members = [np.where(entry_pair == i)[0] for i in range(n_pairs)]

    observed_auroc = auroc(labels, scores)
    observed_auprc = auprc(labels, scores)
    prevalence = n_pos_entries / n_entries
    print(f"Universe: {n_entries} entries / {n_pairs} pairs, "
          f"{n_pos_entries} positive entries / {n_pos_pairs} positive pairs")
    print(f"Observed primary AUROC = {observed_auroc:.4f}, AUPRC = {observed_auprc:.4f}, "
          f"baseline prevalence = {prevalence:.4f}")

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "cluster_bootstrap_primary.py",
        "source": str(TABLES2.relative_to(ROOT)),
        "frame": {
            "description": "PRIMARY lineage-level frame (TableS2): signed DD over 110 "
                           "driver×paralog×lineage entries (Ovarian/Endometrial/Cervical)",
            "n_entries": n_entries,
            "n_unique_pairs": n_pairs,
            "n_positive_entries": n_pos_entries,
            "n_positive_pairs": n_pos_pairs,
            "score": "dependency_dd (signed)",
            "label": "is_known_paralog_sl",
        },
        "seed": SEED,
        "definitions": "auroc/auprc imported from compute_headline_metrics.py "
                       "(rank-based Mann-Whitney AUROC, average ties; AUPRC = average precision)",
    }

    # ── (a) + (e) Pair-clustered bootstrap: AUROC and AUPRC ─────────────
    rng = np.random.default_rng(SEED)
    boot_auroc = np.full(N_BOOT, np.nan)
    boot_auprc = np.full(N_BOOT, np.nan)
    n_skipped = 0
    for b in range(N_BOOT):
        draw = rng.integers(0, n_pairs, n_pairs)  # resample pairs with replacement
        idx = np.concatenate([pair_members[j] for j in draw])
        yb = labels[idx]
        if yb.sum() < 2:
            n_skipped += 1
            continue
        sb = scores[idx]
        boot_auroc[b] = auroc(yb, sb)
        boot_auprc[b] = auprc(yb, sb)
    n_used = N_BOOT - n_skipped
    auroc_vals = boot_auroc[~np.isnan(boot_auroc)]
    auprc_vals = boot_auprc[~np.isnan(boot_auprc)]
    out["pair_clustered_bootstrap_auroc"] = {
        "observed": float(observed_auroc),
        "boot_mean": float(auroc_vals.mean()),
        "ci95_percentile": [float(np.percentile(auroc_vals, 2.5)),
                            float(np.percentile(auroc_vals, 97.5))],
        "n_boot": N_BOOT,
        "n_used": int(n_used),
        "n_skipped_fewer_than_2_positives": int(n_skipped),
        "frac_skipped": float(n_skipped / N_BOOT),
        "scheme": "resample 72 pairs with replacement; carry all lineage entries "
                  "of each resampled pair; AUROC over resulting entry multiset",
    }
    out["pair_clustered_bootstrap_auprc"] = {
        "observed": float(observed_auprc),
        "boot_mean": float(auprc_vals.mean()),
        "ci95_percentile": [float(np.percentile(auprc_vals, 2.5)),
                            float(np.percentile(auprc_vals, 97.5))],
        "baseline_prevalence": float(prevalence),
        "n_boot": N_BOOT,
        "n_used": int(n_used),
        "n_skipped_fewer_than_2_positives": int(n_skipped),
        "frac_skipped": float(n_skipped / N_BOOT),
        "scheme": "same resampled draws and skip rule as the AUROC bootstrap",
    }
    print(f"[a] clustered bootstrap AUROC: mean={auroc_vals.mean():.4f}, "
          f"95% CI [{np.percentile(auroc_vals, 2.5):.4f}, "
          f"{np.percentile(auroc_vals, 97.5):.4f}], skipped {n_skipped}/{N_BOOT}")
    print(f"[e] clustered bootstrap AUPRC: mean={auprc_vals.mean():.4f}, "
          f"95% CI [{np.percentile(auprc_vals, 2.5):.4f}, "
          f"{np.percentile(auprc_vals, 97.5):.4f}] (baseline {prevalence:.4f})")

    # ── (g) Driver-block bootstrap: resample driver genes ──────────────
    # Controls and curated pairs share driver genes, so pairs under one
    # driver are non-independent. Resample driver blocks (each block = all
    # lineage entries of every pair under that driver) with replacement.
    driver_keys = df["driver_gene"].values
    driver_ids = sorted(set(driver_keys))
    driver_to_int = {d: i for i, d in enumerate(driver_ids)}
    entry_driver = np.array([driver_to_int[d] for d in driver_keys])
    n_drivers = len(driver_ids)
    driver_members = [np.where(entry_driver == i)[0] for i in range(n_drivers)]

    rng = np.random.default_rng(SEED)
    boot_db = np.full(N_BOOT, np.nan)
    n_skipped_db = 0
    for b in range(N_BOOT):
        draw = rng.integers(0, n_drivers, n_drivers)  # resample driver blocks
        idx = np.concatenate([driver_members[j] for j in draw])
        yb = labels[idx]
        if yb.sum() < 2:
            n_skipped_db += 1
            continue
        boot_db[b] = auroc(yb, scores[idx])
    db_vals = boot_db[~np.isnan(boot_db)]
    out["driver_block_bootstrap_auroc"] = {
        "observed": float(observed_auroc),
        "boot_mean": float(db_vals.mean()),
        "ci95_percentile": [float(np.percentile(db_vals, 2.5)),
                            float(np.percentile(db_vals, 97.5))],
        "n_boot": N_BOOT,
        "n_used": int(N_BOOT - n_skipped_db),
        "n_skipped_fewer_than_2_positives": int(n_skipped_db),
        "frac_skipped": float(n_skipped_db / N_BOOT),
        "n_driver_blocks": int(n_drivers),
        "scheme": "resample driver genes with replacement; carry all lineage "
                  "entries of every pair under each resampled driver; AUROC "
                  "over the resulting entry multiset; same skip rule as the "
                  "pair-clustered bootstrap",
    }
    print(f"[g] driver-block bootstrap AUROC ({n_drivers} driver blocks): "
          f"mean={db_vals.mean():.4f}, 95% CI [{np.percentile(db_vals, 2.5):.4f}, "
          f"{np.percentile(db_vals, 97.5):.4f}], skipped {n_skipped_db}/{N_BOOT}")

    # ── (b) Pair-level permutation test of the primary statistic ────────
    rng = np.random.default_rng(SEED)
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        pl = rng.permutation(pair_label)
        null[i] = auroc(pl[entry_pair], scores)
    emp_p = (1 + float(np.sum(null >= observed_auroc))) / (1 + N_PERM)
    out["pair_level_permutation_primary"] = {
        "observed": float(observed_auroc),
        "null_mean": float(null.mean()),
        "null_std": float(null.std()),
        "empirical_p": float(emp_p),
        "n_permutations": N_PERM,
        "scheme": "permute pair-level labels (6 positive / 66 control); carry labels "
                  "down to that pair's lineage entries; AUROC over the 110 entries",
        "p_formula": "(1 + #{null >= observed}) / (1 + n_perm)",
    }
    print(f"[b] pair-label permutation: null mean={null.mean():.4f}, p={emp_p:.5f}")

    # ── (c) Effective sample size (ICC of DD by pair) ───────────────────
    k = n_pairs
    cluster_sizes = np.array([len(m) for m in pair_members], dtype=float)
    m_bar = n_entries / k  # average cluster size (instruction's m)
    grand = scores.mean()
    ssb = float(sum(len(m) * (scores[m].mean() - grand) ** 2 for m in pair_members))
    ssw = float(sum(((scores[m] - scores[m].mean()) ** 2).sum() for m in pair_members))
    msb = ssb / (k - 1)
    msw = ssw / (n_entries - k)
    m0 = (n_entries - float((cluster_sizes ** 2).sum()) / n_entries) / (k - 1)
    icc = (msb - msw) / (msb + (m0 - 1) * msw)
    neff = n_entries / (1 + (m_bar - 1) * icc)
    out["effective_sample_size"] = {
        "n_entries": n_entries,
        "n_clusters": int(k),
        "avg_cluster_size_m": float(m_bar),
        "m0_unequal_size_correction": float(m0),
        "cluster_size_distribution": {str(int(s)): int((cluster_sizes == s).sum())
                                      for s in np.unique(cluster_sizes)},
        "anova_icc_dependency_dd_by_pair": float(icc),
        "msb": float(msb), "msw": float(msw),
        "neff_formula": "NEFF = 110 / (1 + (m-1)*ICC)",
        "neff": float(neff),
        "note": "one-way random-effects ANOVA ICC = (MSB-MSW)/(MSB+(m0-1)*MSW); "
                "negative ICC implies within-pair variance exceeds between-pair "
                "variance (no clustering penalty)",
    }
    print(f"[c] m={m_bar:.4f}, ICC={icc:.4f}, NEFF={neff:.1f}")

    # ── (d) Permutation test for the MEAN-aggregated per-pair frame ──────
    g = df.groupby(["driver_gene", "paralog_gene"], as_index=False).agg(
        known=("is_known_paralog_sl", "max"),
        dd=("dependency_dd", "mean"),
    )
    yg = g["known"].astype(int).values
    sg = g["dd"].values.astype(float)
    obs_mean_agg = auroc(yg, sg)
    rng = np.random.default_rng(SEED)
    null_mean_agg = np.empty(N_PERM)
    for i in range(N_PERM):
        null_mean_agg[i] = auroc(rng.permutation(yg), sg)
    emp_p_mean = (1 + float(np.sum(null_mean_agg >= obs_mean_agg))) / (1 + N_PERM)
    out["per_pair_mean_aggregation_permutation"] = {
        "observed": float(obs_mean_agg),
        "null_mean": float(null_mean_agg.mean()),
        "null_std": float(null_mean_agg.std()),
        "empirical_p": float(emp_p_mean),
        "n_permutations": N_PERM,
        "n_pairs": int(len(g)),
        "n_positives": int(yg.sum()),
        "aggregation": "mean signed DD across lineages per pair (the 0.672 frame of "
                       "ml_benchmark.py / compute_headline_metrics.py; never previously "
                       "permutation-tested)",
        "p_formula": "(1 + #{null >= observed}) / (1 + n_perm)",
    }
    print(f"[d] mean-agg permutation: observed={obs_mean_agg:.4f}, "
          f"null mean={null_mean_agg.mean():.4f}, p={emp_p_mean:.5f}")

    # ── (f) Precision@k on the primary frame ─────────────────────────────
    order = np.argsort(-scores, kind="mergesort")
    y_sorted = labels[order]
    prec = {}
    for kk in (10, 20):
        top = y_sorted[:kk]
        prec[f"precision_at_{kk}"] = {
            "n_positives_in_top_k": int(top.sum()),
            "k": kk,
            "precision": float(top.mean()),
        }
    out["precision_at_k_primary"] = {
        **prec,
        "ranking": "entries ranked by signed dependency_dd, descending",
        "n_positives_total": n_pos_entries,
        "n_entries": n_entries,
    }
    print(f"[f] precision@10={prec['precision_at_10']['precision']:.3f}, "
          f"precision@20={prec['precision_at_20']['precision']:.3f}")

    JSON_OUT.write_text(json.dumps(out, indent=2, allow_nan=False, default=str))
    print(f"\nWrote {JSON_OUT}")


if __name__ == "__main__":
    main()
