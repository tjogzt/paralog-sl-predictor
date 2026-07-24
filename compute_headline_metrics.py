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

Inputs (all small artifacts; no DepMap re-download, no sklearn required —
AUROC is computed with the rank-based Mann-Whitney formula, identical to
sklearn.metrics.roc_auc_score for binary labels):

    output/tables/TableS2_FullResults.tsv   118 driver x paralog x lineage entries
    output/validation_report.json           per-pair negative control + bootstrap
    output/paralog_identity.csv             k-mer Jaccard identity (optional;
                                            run compute_sequence_identity.R first)

Manuscript definitions (see manuscript.tex):
  * Full set:    12 gold-standard positives (TableS3), lineage-level entries.
  * Primary set: 10 true sequence paralogs — excludes the two functional
                 analogs BRCA1<->BRCA2 and STK11->SIK1 (Methods: "Secondary
                 set of 2 functional analogs").
  * Leave-out:   excludes the two pairs whose evidence is DepMap-era
                 (FBXW7->FBXW2, PPP2R1A->PPP2R1B).
  * Pre-DepMap:  keeps only the 8 pairs with pre-DepMap experimental evidence
                 (12 - 2 functional analogs - 2 DepMap-era).
  * Scoring:     AUROC of |DD| (manuscript: "DD alone, using only |DD|").

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

# Pairs excluded from the Primary set (manuscript Methods: functional analogs)
FUNCTIONAL_ANALOGS = {("BRCA1", "BRCA2"), ("STK11", "SIK1")}
# Pairs whose gold-standard evidence is DepMap-era (manuscript Results)
DEPMAP_ERA = {("FBXW7", "FBXW2"), ("PPP2R1A", "PPP2R1B")}

# Published CV3 AUROC values — literature constants, NOT recomputed here.
# Source: Feng et al. (2024) Nat Commun 15:9058, Supplementary Data 1,
# CV3 (gene-pair isolation), NSMRand negative sampling, 1:1 pos:neg ratio,
# complete dataset. Values cross-checked against main-text Table 3 F1 scores.
PUBLISHED_BENCHMARKS = {
    "SLMGAE": 0.790, "NSF4SL": 0.683, "GCATSL": 0.678, "GRSMF": 0.656,
    "PiLSL": 0.626, "KG4SL": 0.563, "SLGNN": 0.530, "PTGNN": 0.529,
}

# Values stated in the manuscript, for the automated claims check.
# Updated 2026-07-25 to the recomputed values after the DD-sign/BH fixes and
# the manuscript text alignment (see manuscript.tex git history).
MANUSCRIPT_CLAIMS = {
    "dd_auroc_lineage_full": 0.794,
    "dd_auroc_lineage_primary": 0.837,
    "dd_auroc_lineage_leave_out_depmap_era": 0.833,
    "dd_auroc_lineage_pre_depmap_only": 0.900,
    "dd_auroc_id_filter_0.2": 0.778,
    "dd_auroc_id_filter_0.3": 1.000,
    "component_dd": 0.794,
    "component_pcs": 0.720,
    "component_delta_expression": 0.348,
    "component_necessity": 0.597,
    "per_pair_auroc": 0.668,
    "dd_auroc_per_pair_mean": 0.736,
    "composite_auroc_per_pair_mean": 0.730,
    "auprc_composite_per_pair": 0.271,
    "llo_auroc_min": 0.763,
    "llo_auroc_max": 0.821,
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


def pair_key(df):
    return list(zip(df["driver_gene"], df["paralog_gene"]))


def lineage_metrics(df, label):
    yt = df["is_known_paralog_sl"].astype(int).values
    ys = df["dependency_dd"].abs().fillna(0).values
    return {
        "auroc": auroc(yt, ys),
        "n_entries": int(len(df)),
        "n_positives": int(yt.sum()),
        "label": label,
    }


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
        "scoring": "AUROC of |DD| (manuscript convention)",
        "sign_convention": "DD = mean(Chronos|WT) - mean(Chronos|MUT); positive = compensation",
    }

    # ── 1. Lineage-level AUROC: full / primary / leave-out / pre-DepMap ──
    full = lineage_metrics(df, "Full set (12 gold-standard pairs)")
    metrics["lineage_full"] = full

    mask_primary = ~pd.Series([k in FUNCTIONAL_ANALOGS for k in keys], index=df.index)
    metrics["lineage_primary"] = lineage_metrics(
        df[mask_primary], "Primary set (10 true sequence paralogs; functional analogs excluded)")
    metrics["lineage_primary"]["excluded_pairs"] = sorted(f"{a}->{b}" for a, b in FUNCTIONAL_ANALOGS)
    # Variant for transparency: functional analogs kept but relabelled as negatives
    yt_variant = df["is_known_paralog_sl"].astype(int).copy()
    yt_variant[pd.Series([k in FUNCTIONAL_ANALOGS for k in keys], index=df.index)] = 0
    metrics["lineage_primary_variant_analogs_as_negatives"] = {
        "auroc": auroc(yt_variant.values, df["dependency_dd"].abs().fillna(0).values),
        "n_entries": int(len(df)),
        "n_positives": int(yt_variant.sum()),
    }

    mask_leave = ~pd.Series([k in DEPMAP_ERA for k in keys], index=df.index)
    metrics["lineage_leave_out_depmap_era"] = lineage_metrics(
        df[mask_leave], "Leave-out (DepMap-era pairs removed)")
    metrics["lineage_leave_out_depmap_era"]["excluded_pairs"] = sorted(f"{a}->{b}" for a, b in DEPMAP_ERA)

    mask_pre = ~pd.Series([(k in DEPMAP_ERA) or (k in FUNCTIONAL_ANALOGS) for k in keys], index=df.index)
    metrics["lineage_pre_depmap_only"] = lineage_metrics(
        df[mask_pre], "Pre-DepMap evidence only (8 pairs)")

    # ── 2. Component decomposition (same lineage-level universe) ──
    yt = df["is_known_paralog_sl"].astype(int).values
    comp = {
        "dd": auroc(yt, df["dependency_dd"].abs().fillna(0).values),
        "pcs": auroc(yt, df["pcs"].fillna(0).values),
        "delta_expression_abs": auroc(yt, df["delta_expression"].abs().fillna(0).values),
        "delta_expression_signed": auroc(yt, df["delta_expression"].fillna(0).values),
        "necessity": auroc(yt, df["necessity"].fillna(0).values),
    }
    metrics["component_decomposition_lineage"] = comp

    # ── 3. Sequence-identity filter (needs output/paralog_identity.csv) ──
    if IDENTITY_CSV.exists():
        ident = pd.read_csv(IDENTITY_CSV)
        id_map = {(r.driver_gene, r.paralog_gene): r.kmer_jaccard for r in ident.itertuples()}
        df["_kmer_id"] = [id_map.get(k, np.nan) for k in keys]
        metrics["identity_filter"] = {
            "source": str(IDENTITY_CSV.relative_to(ROOT)),
            "metric": "k-mer Jaccard (k=3), Methods-validated proxy for >=30% identity",
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
    # This is the "head-to-head" universe of the manuscript (77 unique pairs,
    # |DD| averaged across lineages per pair). NOTE: validation_report.json's
    # per-pair value (0.6685) derives from a different aggregation/source;
    # both are reported — see manuscript "Evaluation frameworks".
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
        "n_pairs": int(len(g)),
        "n_positives": int(yt_g.sum()),
        "aggregation": "mean across lineages per pair, from TableS2",
    }

    # ── 4c. Per-pair MAX aggregation — reproduces validation_report.json's
    # observed per-pair AUROC (0.6685) exactly; documented framework of the
    # bootstrap/negative-control analyses.
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

    # ── 4d. Composite score on the per-pair mean frame (head-to-head universe).
    # Manuscript: "the composite score (AUROC)" and "AUPRC reached 0.271
    # (2.6x baseline prevalence of 0.104)" — both reproduce on this frame.
    gcomp = df.groupby(["driver_gene", "paralog_gene"], as_index=False).agg(
        known=("is_known_paralog_sl", "max"),
        comp=("composite_score", "mean"),
    )
    yt_c = gcomp["known"].astype(int).values
    metrics["per_pair_composite_mean"] = {
        "auroc": auroc(yt_c, gcomp["comp"].values),
        "auprc": auprc(yt_c, gcomp["comp"].values),
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

    # ── 6. Automated claims check against manuscript values ──
    computed_map = {
        "dd_auroc_lineage_full": metrics["lineage_full"]["auroc"],
        "dd_auroc_lineage_primary": metrics["lineage_primary"]["auroc"],
        "dd_auroc_lineage_leave_out_depmap_era": metrics["lineage_leave_out_depmap_era"]["auroc"],
        "dd_auroc_lineage_pre_depmap_only": metrics["lineage_pre_depmap_only"]["auroc"],
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

    JSON_OUT.write_text(json.dumps(metrics, indent=2, allow_nan=False, default=str))
    print(f"Wrote {JSON_OUT}")

    # Flat TSV for R consumption
    rows = []

    def add(metric, value, provenance):
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return
        rows.append({"metric": metric, "value": f"{float(value):.4f}", "provenance": provenance})

    add("dd_auroc_lineage_full", metrics["lineage_full"]["auroc"], "recomputed:TableS2")
    add("dd_auroc_lineage_primary", metrics["lineage_primary"]["auroc"], "recomputed:TableS2")
    add("dd_auroc_lineage_leave_out_depmap_era", metrics["lineage_leave_out_depmap_era"]["auroc"], "recomputed:TableS2")
    add("dd_auroc_lineage_pre_depmap_only", metrics["lineage_pre_depmap_only"]["auroc"], "recomputed:TableS2")
    add("component_dd", comp["dd"], "recomputed:TableS2")
    add("component_pcs", comp["pcs"], "recomputed:TableS2")
    add("component_delta_expression", comp["delta_expression_abs"], "recomputed:TableS2")
    add("component_necessity", comp["necessity"], "recomputed:TableS2")
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
