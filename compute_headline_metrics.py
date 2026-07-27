#!/usr/bin/env python3
"""
compute_headline_metrics.py
===========================
Single source of truth for the manuscript's headline AUROC metrics.

Recomputes every headline number that can be derived from pipeline output
artifacts, and records provenance for those that cannot. Figure and table
scripts (tables.py, validation_viz.py, R_fig1.R) MUST read the generated
files instead of using hard-coded literals:

    output/headline_metrics.json        full detail + claims check
    output/tables/headline_metrics.tsv  flat metric/value/provenance table (for R)
    output/tables/TableS7_EffectSizes.tsv  per-entry DD/Cohen's d/Hedges' g

Inputs (all small artifacts; no DepMap re-download, no sklearn required —
AUROC is computed with the rank-based Mann-Whitney formula, identical to
sklearn.metrics.roc_auc_score for binary labels; AUPRC is average precision,
identical to sklearn.metrics.average_precision_score):

    output/tables/TableS2_FullResults.tsv   driver x paralog x lineage entries
    output/validation_report.json           per-pair negative control + bootstrap
    output/paralog_identity.csv             k-mer Jaccard identity (optional;
                                            run compute_sequence_identity.R first)

Manuscript definitions (see manuscript.tex) — evidence tiers after the
round-4 methods review (Supplementary Table S3):
  * Tier A:      3 pairs with DIRECT genetic synthetic-lethal evidence from
                 dual-gene perturbation: AKT1->AKT2 (Najm 2018 combinatorial
                 CRISPR), CDK4->CDK6 and MAP2K1->MAP2K2 (Parrish 2021 pgPEN).
  * Tier B:      2 pairs with natural-genotype conditional dependency +
                 functional validation: SMARCA4->SMARCA2 (Hoffman 2014),
                 ARID1A->ARID1B (Helming 2014).
  * PRIMARY external benchmark: Tier A ∪ Tier B (5 pairs).
  * Tier C:      5 pairs with indirect evidence only — EP300->CREBBP
                 (reciprocal direction only), PIK3CA->PIK3CB (Wee 2008
                 supports PTEN->PIK3CB only), CCNE1->CCNE2 (mouse
                 developmental double-KO redundancy), and the two
                 DepMap-derived pairs (FBXW7->FBXW2, PPP2R1A->PPP2R1B).
                 Excluded from the primary benchmark.
  * Comparators: BRCA1<->BRCA2, STK11->SIK1 — not sequence paralogs;
                 specificity references only.
  * Full set:    all 12 curated pairs (secondary analysis).
  * Leave-out:   excludes the two DepMap-era pairs.
  * Pre-DepMap:  keeps only the 8 pairs with pre-DepMap experimental evidence
                 (12 - 2 comparators - 2 DepMap-era).
  * Scoring:     AUROC of |DD| (manuscript: "DD alone, using only |DD|");
                 AUPRC (average precision) reported alongside.
  * Bootstrap:   10,000 resamples of evaluation entries (or pairs for
                 per-pair frames), percentile 95% CI; paired resamples for
                 component head-to-head comparisons.

Usage: python compute_headline_metrics.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
TABLES2 = ROOT / "output" / "tables" / "TableS2_FullResults.tsv"
VALIDATION_REPORT = ROOT / "output" / "validation_report.json"
IDENTITY_CSV = ROOT / "output" / "paralog_identity.csv"
JSON_OUT = ROOT / "output" / "headline_metrics.json"
TSV_OUT = ROOT / "output" / "tables" / "headline_metrics.tsv"
EFFECTS_OUT = ROOT / "output" / "tables" / "TableS7_EffectSizes.tsv"

N_BOOT = 10_000
BOOT_SEED = 42

# ── Evidence tiers (round-4 review; Supplementary Table S3) ────────────
# Comparators: mechanistic reference pairs, NOT sequence paralogs; scored
# as UNORDERED (BRCA1<->BRCA2 appears in both directions in TableS2).
FUNCTIONAL_ANALOGS = {("BRCA1", "BRCA2"), ("STK11", "SIK1")}
# Tier A: direct dual-perturbation genetic SL evidence.
TIER_A = {
    ("AKT1", "AKT2"),      # Najm 2018 combinatorial CRISPR digenic KO
    ("CDK4", "CDK6"),      # Parrish 2021 digenic KO (pgPEN)
    ("MAP2K1", "MAP2K2"),  # Parrish 2021 digenic KO (pgPEN)
}
# Tier B: natural-genotype conditional dependency + functional validation.
TIER_B = {
    ("SMARCA4", "SMARCA2"),  # Hoffman 2014 CRISPR in SMARCA4-mutant lines
    ("ARID1A", "ARID1B"),    # Helming 2014 shRNA in ARID1A-mutant lines
}
# Tier C: indirect evidence only (excluded from the primary benchmark).
TIER_C_INDIRECT = {
    ("EP300", "CREBBP"),   # Ogiwara 2016 + Nie 2021, reciprocal direction only
    ("PIK3CA", "PIK3CB"),  # Wee 2008 supports PTEN->PIK3CB only
    ("CCNE1", "CCNE2"),    # Geng 2003 mouse developmental double-KO redundancy
}
DEPMAP_ERA = {("FBXW7", "FBXW2"), ("PPP2R1A", "PPP2R1B")}  # DepMap-derived
TIER_C = TIER_C_INDIRECT | DEPMAP_ERA
# Reciprocal-direction set retained for the direction-strict sensitivity
# analysis (EP300->CREBBP relabelled non-positive).
TIER_RECIPROCAL = {("EP300", "CREBBP")}

# Published CV3 AUROC values — literature constants, NOT recomputed here.
# Source: Feng et al. (2024) Nat Commun 15:9058, Supplementary Data 1,
# CV3 (gene-pair isolation), NSMRand negative sampling, 1:1 pos:neg ratio,
# complete dataset. Values cross-checked against main-text Table 3 F1 scores.
PUBLISHED_BENCHMARKS = {
    "SLMGAE": 0.790, "NSF4SL": 0.683, "GCATSL": 0.678, "GRSMF": 0.656,
    "PiLSL": 0.626, "KG4SL": 0.563, "SLGNN": 0.530, "PTGNN": 0.529,
}

# Values stated in the manuscript, for the automated claims check.
# Updated 2026-07-26 to the round-4 primary framework: minimum 5 mutant +
# 5 WT cell lines per driver x lineage stratum, Tier A∪B primary benchmark
# (Tier A pairs have no >=5/>=5 stratum in the gyn3 frame; primary carried
# by the two Tier B pairs — see manuscript "Evaluation frameworks").
MANUSCRIPT_CLAIMS = {
    "dd_auroc_lineage_full": 0.676,
    "auprc_lineage_full": 0.386,
    "dd_auroc_lineage_tier_ab": 1.000,
    "auprc_lineage_tier_ab": 1.000,
    "dd_auroc_lineage_leave_out_depmap_era": 0.725,
    "dd_auroc_lineage_pre_depmap_only": 0.774,
    "dd_auroc_lineage_full_direction_strict": 0.676,
    "dd_auroc_id_filter_0.2": 0.583,
    "dd_auroc_id_filter_0.3": 1.000,
    "component_dd": 0.676,
    "component_pcs": 0.825,
    "component_delta_expression": 0.547,
    "component_necessity": 0.642,
    "per_pair_auroc": 0.500,
    "dd_auroc_per_pair_mean": 0.566,
    "composite_auroc_per_pair_mean": 0.831,
    "auprc_composite_per_pair": 0.357,
    "llo_auroc_min": 0.656,
    "llo_auroc_max": 0.704,
}
TOL = 0.005  # claims match if |computed - claimed| <= TOL


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


def auprc(labels, scores):
    """Average precision (identical to sklearn.metrics.average_precision_score):
    mean of precision@k over the ranks k of the positives, scores sorted desc."""
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    ok = ~np.isnan(scores)
    labels, scores = labels[ok], scores[ok]
    n_pos = int(labels.sum())
    if n_pos == 0 or n_pos == len(labels):
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    return float(np.mean([y[: i + 1].mean() for i in range(len(y)) if y[i] == 1]))


def bootstrap_ci(labels, scores, stat_fn, n_boot=N_BOOT, seed=BOOT_SEED):
    """Percentile 95% CI of stat_fn under entry resampling with replacement."""
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    ok = ~np.isnan(scores)
    labels, scores = labels[ok], scores[ok]
    n = len(labels)
    if n == 0:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        v = stat_fn(labels[idx], scores[idx])
        if not np.isnan(v):
            vals.append(v)
    if len(vals) < 100:
        return (float("nan"), float("nan"))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return (float(lo), float(hi))


def pair_key(df):
    return list(zip(df["driver_gene"], df["paralog_gene"]))


def lineage_metrics(df, label, with_ci=True):
    yt = df["is_known_paralog_sl"].astype(int).values
    ys = df["dependency_dd"].abs().fillna(0).values
    out = {
        "auroc": auroc(yt, ys),
        "auprc": auprc(yt, ys),
        "n_entries": int(len(df)),
        "n_positives": int(yt.sum()),
        "label": label,
    }
    if with_ci and len(df) >= 20:
        out["auroc_ci95"] = list(bootstrap_ci(yt, ys, auroc))
        out["auprc_ci95"] = list(bootstrap_ci(yt, ys, auprc))
    return out


def main():
    if not TABLES2.exists():
        sys.exit(f"ERROR: {TABLES2} not found. Run the pipeline (main.py + tables) first.")

    df = pd.read_csv(TABLES2, sep="\t")
    required = {"driver_gene", "paralog_gene", "dependency_dd", "is_known_paralog_sl"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: TableS2 missing columns: {missing}")

    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)
    keys = pair_key(df)

    metrics = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "compute_headline_metrics.py",
        "source_files": [str(TABLES2.relative_to(ROOT))],
        "scoring": "AUROC of |DD| (manuscript convention); AUPRC = average precision",
        "bootstrap": f"percentile 95% CI, {N_BOOT} entry resamples, seed {BOOT_SEED}",
        "sign_convention": "DD = mean(Chronos|WT) - mean(Chronos|MUT); positive = compensation",
        "min_samples": "primary analysis: >=5 mutant and >=5 WT cell lines per stratum",
    }

    def exclude(df_, pairs, unordered=None):
        unordered = unordered or set()
        ku = [frozenset(k) for k in pair_key(df_)]
        mask = ~pd.Series([(k in pairs) or (u in unordered) for k, u in zip(pair_key(df_), ku)],
                          index=df_.index)
        return df_[mask]

    COMPARATOR_UNORDERED = {frozenset(p) for p in FUNCTIONAL_ANALOGS}

    # ── 1. Lineage-level frames ────────────────────────────────────────
    full = lineage_metrics(df, "Full set (12 curated pairs; secondary)")
    metrics["lineage_full"] = full

    # PRIMARY external benchmark: Tier A ∪ Tier B positives; Tier C and
    # comparator entries removed from the frame.
    tier_ab_df = exclude(df, TIER_C, COMPARATOR_UNORDERED)
    metrics["lineage_tier_ab"] = lineage_metrics(
        tier_ab_df, "PRIMARY: Tier A ∪ Tier B (5 pairs, external direct/conditional evidence)")
    metrics["lineage_tier_ab"]["tier_a_pairs"] = sorted(f"{a}->{b}" for a, b in TIER_A)
    metrics["lineage_tier_ab"]["tier_b_pairs"] = sorted(f"{a}->{b}" for a, b in TIER_B)

    metrics["lineage_tier_a"] = lineage_metrics(
        exclude(df, TIER_B | TIER_C, COMPARATOR_UNORDERED),
        "Tier A only (3 pairs, direct dual-perturbation evidence)")

    metrics["lineage_tier_b"] = lineage_metrics(
        exclude(df, TIER_A | TIER_C, COMPARATOR_UNORDERED),
        "Tier B only (2 pairs, genotype-conditional dependency)")

    mask_primary = ~pd.Series([k in FUNCTIONAL_ANALOGS for k in keys], index=df.index)
    metrics["lineage_primary"] = lineage_metrics(
        df[mask_primary], "Curated sequence paralogs (10 pairs; comparators excluded)")
    metrics["lineage_primary"]["excluded_pairs"] = sorted(f"{a}->{b}" for a, b in FUNCTIONAL_ANALOGS)

    mask_leave = ~pd.Series([k in DEPMAP_ERA for k in keys], index=df.index)
    metrics["lineage_leave_out_depmap_era"] = lineage_metrics(
        df[mask_leave], "Leave-out (DepMap-era pairs removed)")
    metrics["lineage_leave_out_depmap_era"]["excluded_pairs"] = sorted(f"{a}->{b}" for a, b in DEPMAP_ERA)

    mask_pre = ~pd.Series([(k in DEPMAP_ERA) or (k in FUNCTIONAL_ANALOGS) for k in keys], index=df.index)
    metrics["lineage_pre_depmap_only"] = lineage_metrics(
        df[mask_pre], "Pre-DepMap evidence only (8 pairs)")

    # Direction-strict full set: EP300->CREBBP relabelled non-positive because
    # direct experimental evidence supports only the reciprocal direction.
    yt_ds = df["is_known_paralog_sl"].astype(int).copy()
    yt_ds[pd.Series([k in TIER_RECIPROCAL for k in keys], index=df.index)] = 0
    metrics["lineage_full_direction_strict"] = {
        "auroc": auroc(yt_ds.values, df["dependency_dd"].abs().fillna(0).values),
        "auprc": auprc(yt_ds.values, df["dependency_dd"].abs().fillna(0).values),
        "n_entries": int(len(df)),
        "n_positives": int(yt_ds.sum()),
        "label": "Full set, direction-strict (EP300->CREBBP relabelled non-positive)",
    }

    # ── 2. Component decomposition + paired bootstrap (same universe) ──
    yt = df["is_known_paralog_sl"].astype(int).values
    comp_scores = {
        "dd": df["dependency_dd"].abs().fillna(0).values,
        "pcs": df["pcs"].fillna(0).values,
        "delta_expression_abs": df["delta_expression"].abs().fillna(0).values,
        "delta_expression_signed": df["delta_expression"].fillna(0).values,
        "necessity": df["necessity"].fillna(0).values,
    }
    comp = {k: auroc(yt, v) for k, v in comp_scores.items()}
    comp_auprc = {k: auprc(yt, v) for k, v in comp_scores.items()}
    metrics["component_decomposition_lineage"] = comp
    metrics["component_decomposition_lineage_auprc"] = comp_auprc

    # Paired bootstrap: same resampled entries for every component per draw.
    n = len(yt)
    rng = np.random.default_rng(BOOT_SEED)
    boot_stats = {k: [] for k in ("dd", "pcs", "delta_expression_abs", "necessity")}
    boot_delta = {k: [] for k in ("pcs", "delta_expression_abs", "necessity")}
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        yb = yt[idx]
        if yb.sum() == 0 or yb.sum() == len(yb):
            continue
        vals = {k: auroc(yb, v[idx]) for k, v in comp_scores.items() if k != "delta_expression_signed"}
        for k, v in vals.items():
            if not np.isnan(v):
                boot_stats[k].append(v)
        for k in boot_delta:
            if not np.isnan(vals.get(k, float("nan"))) and not np.isnan(vals.get("dd", float("nan"))):
                boot_delta[k].append(vals[k] - vals["dd"])
    paired = {}
    for k, vals in boot_stats.items():
        arr = np.asarray(vals)
        paired[k] = {"ci95": [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))],
                     "n_boot": int(len(arr))}
    for k, vals in boot_delta.items():
        arr = np.asarray(vals)
        paired[f"{k}_minus_dd"] = {
            "mean_delta": float(arr.mean()),
            "ci95": [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))],
            "frac_above_dd": float((arr > 0).mean()),
        }
    metrics["component_paired_bootstrap"] = paired

    # Persist the paired delta distributions (long format) for the Fig. 1d
    # inset histogram drawn by R_fig1.R — same draws as the CIs above.
    deltas_out = pd.concat(
        [pd.DataFrame({"component": f"{k}_minus_dd", "delta": vals})
         for k, vals in boot_delta.items() if len(vals)],
        ignore_index=True)
    deltas_out.to_csv(ROOT / "output" / "component_paired_bootstrap_deltas.csv",
                      index=False)

    # ── 3. Sequence-identity filter (needs output/paralog_identity.csv) ──
    if IDENTITY_CSV.exists():
        ident = pd.read_csv(IDENTITY_CSV)
        id_map = {(r.driver_gene, r.paralog_gene): r.kmer_jaccard for r in ident.itertuples()}
        df["_kmer_id"] = [id_map.get(k, np.nan) for k in keys]
        metrics["identity_filter"] = {
            "source": str(IDENTITY_CSV.relative_to(ROOT)),
            "metric": "k-mer Jaccard (k=3), Methods-validated identity-enrichment proxy",
            "pairs_with_identity": int(df['_kmer_id'].notna().sum()),
        }
        for thr in (0.2, 0.3):
            sub = df[df["_kmer_id"] >= thr]
            metrics["identity_filter"][f"id_ge_{thr}"] = {
                "auroc": auroc(sub["is_known_paralog_sl"].astype(int).values,
                               sub["dependency_dd"].abs().fillna(0).values),
                "n_entries": int(len(sub)),
                "n_unique_pairs": int(sub[["driver_gene", "paralog_gene"]].drop_duplicates().shape[0]),
                "n_positives": int(sub["is_known_paralog_sl"].sum()),
            }
    else:
        metrics["identity_filter"] = {
            "status": "not_reproducible_from_artifacts",
            "reason": "output/paralog_identity.csv missing; run compute_sequence_identity.R",
        }

    # ── 4. Per-pair framework (from pipeline's own validation report) ──
    if VALIDATION_REPORT.exists():
        vr = json.loads(VALIDATION_REPORT.read_text())
        metrics["per_pair_framework"] = {
            "source": str(VALIDATION_REPORT.relative_to(ROOT)),
            "observed_auroc": vr.get("negative_control", {}).get("observed_auroc"),
            "null_auroc_mean": vr.get("negative_control", {}).get("null_auroc_mean"),
            "empirical_p_value": vr.get("negative_control", {}).get("empirical_p_value"),
            "n_known": vr.get("negative_control", {}).get("n_known"),
            "n_total": vr.get("negative_control", {}).get("n_total"),
            "n_permutations": vr.get("negative_control", {}).get("n_permutations"),
            "bootstrap": vr.get("bootstrap", {}),
            "expression_only_auroc": vr.get("component_decomposition", {}).get("expression_only"),
        }
        metrics["source_files"].append(str(VALIDATION_REPORT.relative_to(ROOT)))
    else:
        metrics["per_pair_framework"] = {"status": "validation_report.json missing"}

    # ── 4b. Per-pair mean framework recomputed from TableS2 ──
    g = df.groupby(["driver_gene", "paralog_gene"], as_index=False).agg(
        known=("is_known_paralog_sl", "max"),
        dd=("dependency_dd", lambda s: s.abs().mean()),
        pcs=("pcs", "mean"),
        dexpr=("delta_expression", lambda s: s.abs().mean()),
        necessity=("necessity", "mean"),
    )
    yt_g = g["known"].astype(int).values
    metrics["per_pair_mean_from_tables2"] = {
        "auroc_dd": auroc(yt_g, g["dd"].values),
        "auprc_dd": auprc(yt_g, g["dd"].values),
        "auroc_pcs": auroc(yt_g, g["pcs"].values),
        "auroc_delta_expression": auroc(yt_g, g["dexpr"].values),
        "auroc_necessity": auroc(yt_g, g["necessity"].values),
        "auroc_dd_ci95": list(bootstrap_ci(yt_g, g["dd"].values, auroc)),
        "n_pairs": int(len(g)),
        "n_positives": int(yt_g.sum()),
        "aggregation": "mean across lineages per pair, from TableS2",
    }

    # ── 4c. Per-pair MAX aggregation ──
    gmax = df.groupby(["driver_gene", "paralog_gene"], as_index=False).agg(
        known=("is_known_paralog_sl", "max"),
        dd=("dependency_dd", lambda s: s.abs().max()),
    )
    yt_m = gmax["known"].astype(int).values
    metrics["per_pair_max_from_tables2"] = {
        "auroc_dd": auroc(yt_m, gmax["dd"].values),
        "auprc_dd": auprc(yt_m, gmax["dd"].values),
        "n_pairs": int(len(gmax)),
        "n_positives": int(yt_m.sum()),
        "aggregation": "max |DD| across lineages per pair, from TableS2",
    }

    # ── 4d. Composite score on the per-pair mean frame ──
    gcomp = df.groupby(["driver_gene", "paralog_gene"], as_index=False).agg(
        known=("is_known_paralog_sl", "max"),
        comp=("composite_score", "mean"),
    )
    yt_c = gcomp["known"].astype(int).values
    metrics["per_pair_composite_mean"] = {
        "auroc": auroc(yt_c, gcomp["comp"].values),
        "auprc": auprc(yt_c, gcomp["comp"].values),
        "auroc_ci95": list(bootstrap_ci(yt_c, gcomp["comp"].values, auroc)),
        "baseline_prevalence": float(yt_c.mean()),
        "n_pairs": int(len(gcomp)),
        "n_positives": int(yt_c.sum()),
    }

    # ── 4e. Leave-one-lineage-out (lineage-level frame) ──
    llo = {}
    for ct in sorted(df["cancer_type"].unique()):
        sub = df[df["cancer_type"] != ct]
        yt_s = sub["is_known_paralog_sl"].astype(int).values
        llo[f"without_{ct}"] = auroc(yt_s, sub["dependency_dd"].abs().fillna(0).values)
    metrics["leave_one_lineage_out"] = {
        "values": llo,
        "range": [float(np.nanmin(list(llo.values()))), float(np.nanmax(list(llo.values())))],
    }

    # ── 5. Published benchmarks (literature constants, labelled) ──
    metrics["published_benchmarks"] = {
        "values": PUBLISHED_BENCHMARKS,
        "provenance": "Literature constants (Feng et al. 2024, Suppl. Data 1, CV3 NSMRand 1:1); not recomputed",
    }

    # ── 5b. Per-entry effect sizes (Hedges' g) — Table S8 ──
    if "hedges_g" in df.columns:
        eff_cols = ["driver_gene", "paralog_gene", "cancer_type", "dependency_dd",
                    "cohens_d", "hedges_g", "dd_p_value", "q_value",
                    "n_mut", "n_wt", "is_known_paralog_sl"]
        eff = df[[c for c in eff_cols if c in df.columns]].copy()
        eff = eff.sort_values(["is_known_paralog_sl", "dependency_dd"],
                              ascending=[False, False], key=lambda s: s.abs() if s.name == "dependency_dd" else s)
        eff.to_csv(EFFECTS_OUT, sep="\t", index=False)
        metrics["effect_sizes_table"] = str(EFFECTS_OUT.relative_to(ROOT))
    else:
        metrics["effect_sizes_table"] = "TableS2 lacks hedges_g — rerun tables.py after pcs.py update"

    # ── 6. Automated claims check against manuscript values ──
    computed_map = {
        "dd_auroc_lineage_full": metrics["lineage_full"]["auroc"],
        "auprc_lineage_full": metrics["lineage_full"]["auprc"],
        "dd_auroc_lineage_tier_ab": metrics["lineage_tier_ab"]["auroc"],
        "auprc_lineage_tier_ab": metrics["lineage_tier_ab"]["auprc"],
        "dd_auroc_lineage_tier_a": metrics["lineage_tier_a"]["auroc"],
        "dd_auroc_lineage_tier_b": metrics["lineage_tier_b"]["auroc"],
        "dd_auroc_lineage_primary": metrics["lineage_primary"]["auroc"],
        "dd_auroc_lineage_leave_out_depmap_era": metrics["lineage_leave_out_depmap_era"]["auroc"],
        "dd_auroc_lineage_pre_depmap_only": metrics["lineage_pre_depmap_only"]["auroc"],
        "dd_auroc_lineage_full_direction_strict": metrics["lineage_full_direction_strict"]["auroc"],
        "component_dd": comp["dd"],
        "component_pcs": comp["pcs"],
        "component_delta_expression": comp["delta_expression_abs"],
        "component_necessity": comp["necessity"],
        "per_pair_auroc": (metrics["per_pair_framework"] or {}).get("observed_auroc"),
        "dd_auroc_per_pair_mean": metrics["per_pair_mean_from_tables2"]["auroc_dd"],
        "composite_auroc_per_pair_mean": metrics["per_pair_composite_mean"]["auroc"],
        "auprc_composite_per_pair": metrics["per_pair_composite_mean"]["auprc"],
        "llo_auroc_min": metrics["leave_one_lineage_out"]["range"][0],
        "llo_auroc_max": metrics["leave_one_lineage_out"]["range"][1],
    }
    idf = metrics["identity_filter"]
    if "id_ge_0.2" in idf:
        computed_map["dd_auroc_id_filter_0.2"] = idf["id_ge_0.2"]["auroc"]
        computed_map["dd_auroc_id_filter_0.3"] = idf["id_ge_0.3"]["auroc"]

    checks = []
    for name, claimed in MANUSCRIPT_CLAIMS.items():
        got = computed_map.get(name)
        if got is None or (isinstance(got, float) and np.isnan(got)):
            status = "not_reproducible_from_artifacts"
        else:
            status = "match" if abs(got - claimed) <= TOL else "MISMATCH"
        checks.append({"metric": name, "claimed": claimed,
                       "computed": None if got is None or (isinstance(got, float) and np.isnan(got)) else round(got, 4),
                       "status": status})
    metrics["manuscript_claims_check"] = checks

    def _jsonable(obj):
        if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
            return None
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_jsonable(v) for v in obj]
        return obj

    JSON_OUT.write_text(json.dumps(_jsonable(metrics), indent=2, allow_nan=False, default=str))
    print(f"Wrote {JSON_OUT}")

    # Key-number console summary (always visible, even when claims mismatch)
    for key in ("lineage_full", "lineage_tier_ab", "lineage_tier_a", "lineage_tier_b",
                "lineage_leave_out_depmap_era", "lineage_pre_depmap_only"):
        fr = metrics[key]
        print(f"  {key:34s} AUROC={fr['auroc']!s:>7} AUPRC={fr['auprc']!s:>7} "
              f"n={fr['n_entries']:4d} pos={fr['n_positives']:3d}")
    print(f"  per_pair_mean auroc={metrics['per_pair_mean_from_tables2']['auroc_dd']:.4f} "
          f"composite={metrics['per_pair_composite_mean']['auroc']:.4f} "
          f"LLO={metrics['leave_one_lineage_out']['range']}")

    # Flat TSV for R consumption
    rows = []

    def add(metric, value, provenance):
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return
        rows.append({"metric": metric, "value": f"{float(value):.4f}", "provenance": provenance})

    def add_ci(prefix, frame, provenance):
        for stat in ("auroc", "auprc"):
            ci = frame.get(f"{stat}_ci95")
            if ci and not any(np.isnan(ci)):
                add(f"{prefix}_{stat}_ci_lo", ci[0], provenance)
                add(f"{prefix}_{stat}_ci_hi", ci[1], provenance)

    add("dd_auroc_lineage_full", metrics["lineage_full"]["auroc"], "recomputed:TableS2")
    add("auprc_lineage_full", metrics["lineage_full"]["auprc"], "recomputed:TableS2")
    add_ci("dd_lineage_full", metrics["lineage_full"], "recomputed:TableS2(bootstrap)")
    add("dd_auroc_lineage_tier_ab", metrics["lineage_tier_ab"]["auroc"], "recomputed:TableS2")
    add("auprc_lineage_tier_ab", metrics["lineage_tier_ab"]["auprc"], "recomputed:TableS2")
    add("n_entries_tier_ab", metrics["lineage_tier_ab"]["n_entries"], "recomputed:TableS2")
    add("n_positives_tier_ab", metrics["lineage_tier_ab"]["n_positives"], "recomputed:TableS2")
    add_ci("dd_lineage_tier_ab", metrics["lineage_tier_ab"], "recomputed:TableS2(bootstrap)")
    add("dd_auroc_lineage_tier_a", metrics["lineage_tier_a"]["auroc"], "recomputed:TableS2")
    add("dd_auroc_lineage_tier_b", metrics["lineage_tier_b"]["auroc"], "recomputed:TableS2")
    add("n_positives_tier_a", metrics["lineage_tier_a"]["n_positives"], "recomputed:TableS2")
    add("n_positives_tier_b", metrics["lineage_tier_b"]["n_positives"], "recomputed:TableS2")
    add("dd_auroc_lineage_primary", metrics["lineage_primary"]["auroc"], "recomputed:TableS2")
    add("dd_auroc_lineage_leave_out_depmap_era", metrics["lineage_leave_out_depmap_era"]["auroc"], "recomputed:TableS2")
    add("dd_auroc_lineage_pre_depmap_only", metrics["lineage_pre_depmap_only"]["auroc"], "recomputed:TableS2")
    add("dd_auroc_lineage_full_direction_strict", metrics["lineage_full_direction_strict"]["auroc"], "recomputed:TableS2")
    add("component_dd", comp["dd"], "recomputed:TableS2")
    add("component_pcs", comp["pcs"], "recomputed:TableS2")
    add("component_delta_expression", comp["delta_expression_abs"], "recomputed:TableS2")
    add("component_necessity", comp["necessity"], "recomputed:TableS2")
    add("auprc_component_dd", comp_auprc["dd"], "recomputed:TableS2")
    for k in ("dd", "pcs", "delta_expression_abs", "necessity"):
        pb = paired.get(k, {})
        if pb.get("ci95"):
            add(f"component_{k}_ci_lo", pb["ci95"][0], "recomputed:TableS2(paired bootstrap)")
            add(f"component_{k}_ci_hi", pb["ci95"][1], "recomputed:TableS2(paired bootstrap)")
    for k in ("pcs", "delta_expression_abs", "necessity"):
        pd_ = paired.get(f"{k}_minus_dd", {})
        if pd_.get("ci95"):
            add(f"paired_{k}_minus_dd", pd_["mean_delta"], "recomputed:TableS2(paired bootstrap)")
            add(f"paired_{k}_minus_dd_ci_lo", pd_["ci95"][0], "recomputed:TableS2(paired bootstrap)")
            add(f"paired_{k}_minus_dd_ci_hi", pd_["ci95"][1], "recomputed:TableS2(paired bootstrap)")
    if "id_ge_0.3" in idf:
        add("dd_auroc_id_filter_0.3", idf["id_ge_0.3"]["auroc"], "recomputed:TableS2+paralog_identity")
        add("dd_auroc_id_filter_0.2", idf["id_ge_0.2"]["auroc"], "recomputed:TableS2+paralog_identity")
    pp = metrics.get("per_pair_framework", {})
    add("per_pair_auroc", pp.get("observed_auroc"), "artifact:validation_report.json")
    add("per_pair_null_mean", pp.get("null_auroc_mean"), "artifact:validation_report.json")
    if isinstance(pp.get("empirical_p_value"), (int, float)):
        add("per_pair_empirical_p", pp["empirical_p_value"], "artifact:validation_report.json")
    ppm = metrics["per_pair_mean_from_tables2"]
    add("dd_auroc_per_pair_mean", ppm["auroc_dd"], "recomputed:TableS2(per-pair mean)")
    add("auprc_dd_per_pair_mean", ppm["auprc_dd"], "recomputed:TableS2(per-pair mean)")
    add("component_pcs_per_pair_mean", ppm["auroc_pcs"], "recomputed:TableS2(per-pair mean)")
    add("component_delta_expression_per_pair_mean", ppm["auroc_delta_expression"], "recomputed:TableS2(per-pair mean)")
    add("component_necessity_per_pair_mean", ppm["auroc_necessity"], "recomputed:TableS2(per-pair mean)")
    ppx = metrics["per_pair_max_from_tables2"]
    add("dd_auroc_per_pair_max", ppx["auroc_dd"], "recomputed:TableS2(per-pair max)")
    add("auprc_dd_per_pair_max", ppx["auprc_dd"], "recomputed:TableS2(per-pair max)")
    ppc = metrics["per_pair_composite_mean"]
    add("composite_auroc_per_pair_mean", ppc["auroc"], "recomputed:TableS2(per-pair mean)")
    add("auprc_composite_per_pair", ppc["auprc"], "recomputed:TableS2(per-pair mean)")
    add("llo_auroc_min", metrics["leave_one_lineage_out"]["range"][0], "recomputed:TableS2(LLO)")
    add("llo_auroc_max", metrics["leave_one_lineage_out"]["range"][1], "recomputed:TableS2(LLO)")
    for m, v in PUBLISHED_BENCHMARKS.items():
        add(f"published_{m}", v, "literature:Feng2024_SuppData1_CV3_NSMRand_1to1")

    tsv = pd.DataFrame(rows)
    TSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    tsv.to_csv(TSV_OUT, sep="\t", index=False)
    print(f"Wrote {TSV_OUT}")

    # Console summary
    print("\n=== Manuscript claims check ===")
    for c in checks:
        comp_str = "—" if c["computed"] is None else f"{c['computed']:.4f}"
        print(f"  [{c['status']:>33s}] {c['metric']}: claimed {c['claimed']} vs computed {comp_str}")
    n_bad = sum(1 for c in checks if c["status"] == "MISMATCH")
    print(f"\n{sum(1 for c in checks if c['status'] == 'match')}/{len(checks)} claims match, "
          f"{n_bad} mismatch, {sum(1 for c in checks if c['status'] == 'not_reproducible_from_artifacts')} not reproducible.")
    if n_bad:
        sys.exit(f"FAILED: {n_bad} manuscript claim(s) MISMATCH")


if __name__ == "__main__":
    main()
