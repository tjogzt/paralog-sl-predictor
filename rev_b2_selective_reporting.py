#!/usr/bin/env python3
"""
rev_b2_selective_reporting.py  (Stage-4 revision, item B2)
==========================================================
Selective-reporting repair: assembles (and where missing, computes) the
numbers the reviewers asked to be reported together. No re-run of the main
analysis — everything derives from frozen artifacts; only the uncertainty
quantification of the comparator-excluded frame is newly computed (the
manuscript never reported CI/permutation p for that frame).

  (a) Comparator-excluded frame ("Curated sequence paralogs (10 pairs)",
      compute_headline_metrics.lineage_primary: TableS2 minus the two
      ordered comparator tuples (BRCA1->BRCA2, STK11->SIK1)):
      AUROC + pair-clustered bootstrap 95% CI + pair-level permutation p
      (identical schemes to cluster_bootstrap_primary.py, seed 42, 10,000).
  (b) Label-subset AUROC overview table (full / sequence-paralog-only
      comparator-excluded / Tier A u B / leave-out DepMap-era / pre-DepMap,
      plus tier A, tier B, direction-strict for completeness) from
      output/headline_metrics.json.
  (c) Composition table of the 8 positive entries of the primary frame
      (driver / paralog / lineage / tier / evidence source / signed DD /
      n_mut / n_wt).
  (d) STK11->SIK1 label-identity audit: where the pair is treated as a
      positive vs as a comparator, with the exact mechanism of each
      treatment and its numerical consequence.

Outputs (output/revision_stage4/):
  b2_comparator_excluded_frame.json
  b2_subset_auroc_overview.csv
  b2_positive_entries_composition.csv
  b2_stk11_sik1_label_audit.json

Usage: python rev_b2_selective_reporting.py   (run from repo root)
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

from compute_headline_metrics import auroc, auprc, FUNCTIONAL_ANALOGS  # noqa: E402

SEED = 42
N_BOOT = 10_000
N_PERM = 10_000

TIER_OF = {
    ("AKT1", "AKT2"): "A", ("CDK4", "CDK6"): "A", ("MAP2K1", "MAP2K2"): "A",
    ("SMARCA4", "SMARCA2"): "B", ("ARID1A", "ARID1B"): "B",
    ("EP300", "CREBBP"): "C", ("PIK3CA", "PIK3CB"): "C", ("CCNE1", "CCNE2"): "C",
    ("FBXW7", "FBXW2"): "C", ("PPP2R1A", "PPP2R1B"): "C",
    ("BRCA1", "BRCA2"): "Comparator", ("STK11", "SIK1"): "Comparator",
}
KEY_REF_OF = {
    ("AKT1", "AKT2"): "Najm 2018", ("CDK4", "CDK6"): "Parrish 2021",
    ("MAP2K1", "MAP2K2"): "Parrish 2021", ("SMARCA4", "SMARCA2"): "Hoffman 2014",
    ("ARID1A", "ARID1B"): "Helming 2014", ("EP300", "CREBBP"): "Ogiwara 2016; Nie 2021",
    ("PIK3CA", "PIK3CB"): "Wee 2008", ("CCNE1", "CCNE2"): "Geng 2003",
    ("FBXW7", "FBXW2"): "DepMap 26Q1 release", ("PPP2R1A", "PPP2R1B"): "DepMap 26Q1 release",
    ("BRCA1", "BRCA2"): "Bryant 2005", ("STK11", "SIK1"): "Hollstein 2019",
}
DIRECT_SL_OF = {
    ("AKT1", "AKT2"): "Yes", ("CDK4", "CDK6"): "Yes", ("MAP2K1", "MAP2K2"): "Yes",
    ("SMARCA4", "SMARCA2"): "Conditional", ("ARID1A", "ARID1B"): "Conditional",
    ("EP300", "CREBBP"): "Reciprocal", ("PIK3CA", "PIK3CB"): "No (other driver)",
    ("CCNE1", "CCNE2"): "No (redundancy)", ("FBXW7", "FBXW2"): "No (computational)",
    ("PPP2R1A", "PPP2R1B"): "No (computational)",
    ("BRCA1", "BRCA2"): "Functional analog", ("STK11", "SIK1"): "Pathway axis",
}


def pair_clustered_bootstrap(labels, scores, entry_pair, n_pairs, members,
                             n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    vals, n_skip = [], 0
    for _ in range(n_boot):
        draw = rng.integers(0, n_pairs, n_pairs)
        idx = np.concatenate([members[j] for j in draw])
        yb = labels[idx]
        if yb.sum() < 2:
            n_skip += 1
            continue
        vals.append(auroc(yb, scores[idx]))
    v = np.asarray(vals)
    return {"boot_mean": float(v.mean()),
            "ci95_percentile": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))],
            "n_boot": n_boot, "n_used": int(len(v)), "n_skipped": int(n_skip)}


def pair_level_permutation(labels, scores, entry_pair, pair_label,
                           n_perm=N_PERM, seed=SEED):
    rng = np.random.default_rng(seed)
    obs = auroc(labels, scores)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = auroc(rng.permutation(pair_label)[entry_pair], scores)
    return {"observed": float(obs), "null_mean": float(null.mean()),
            "null_std": float(null.std()),
            "empirical_p": float((1 + np.sum(null >= obs)) / (1 + n_perm)),
            "n_permutations": n_perm}


def main():
    print("=" * 72)
    print("  rev B2: selective-reporting repair (assembly + comparator-excluded UQ)")
    print("=" * 72)

    df = pd.read_csv(ROOT / "output" / "tables" / "TableS2_FullResults.tsv", sep="\t")
    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)
    hm = json.loads((ROOT / "output" / "headline_metrics.json").read_text())

    # ── (a) comparator-excluded frame ─────────────────────────────
    keys = list(zip(df["driver_gene"], df["paralog_gene"]))
    mask = ~pd.Series([k in FUNCTIONAL_ANALOGS for k in keys], index=df.index)
    sub = df[mask].reset_index(drop=True)
    labels = sub["is_known_paralog_sl"].astype(int).to_numpy()
    scores = sub["dependency_dd"].fillna(0).to_numpy()
    obs = auroc(labels, scores)
    ref = hm["lineage_primary"]["auroc"]
    if abs(obs - ref) > 1e-9:
        raise RuntimeError(f"FIDELITY FAIL: comparator-excluded AUROC {obs} vs headline_metrics {ref}")

    pair_keys = list(zip(sub["driver_gene"], sub["paralog_gene"]))
    pair_ids = sorted(set(pair_keys))
    p2i = {p: i for i, p in enumerate(pair_ids)}
    entry_pair = np.array([p2i[k] for k in pair_keys])
    n_pairs = len(pair_ids)
    pair_label = np.array([int(sub.loc[entry_pair == i, "is_known_paralog_sl"].max())
                           for i in range(n_pairs)])
    members = [np.where(entry_pair == i)[0] for i in range(n_pairs)]

    boot = pair_clustered_bootstrap(labels, scores, entry_pair, n_pairs, members)
    perm = pair_level_permutation(labels, scores, entry_pair, pair_label)

    pos_pairs = sorted({f"{d}->{p}" for d, p in zip(sub["driver_gene"], sub["paralog_gene"])
                        if (d, p) not in FUNCTIONAL_ANALOGS and
                        bool(sub.loc[(sub["driver_gene"] == d) & (sub["paralog_gene"] == p),
                                     "is_known_paralog_sl"].max())})
    out_a = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "frame": "comparator-excluded (compute_headline_metrics.lineage_primary): TableS2 "
                 "minus ordered tuples (BRCA1->BRCA2, STK11->SIK1); note BRCA2->BRCA1 "
                 "(reverse-direction comparator entry) REMAINS as a positive",
        "n_entries": int(len(sub)),
        "n_unique_pairs": int(n_pairs),
        "n_positive_entries": int(labels.sum()),
        "n_positive_pairs": int(pair_label.sum()),
        "positive_pairs_in_frame": pos_pairs,
        "auroc": float(obs),
        "auprc": float(auprc(labels, scores)),
        "pair_clustered_bootstrap_auroc": boot,
        "pair_level_permutation": perm,
        "headline_metrics_entry_bootstrap_ci95": hm["lineage_primary"].get("auroc_ci95"),
        "seed": SEED,
    }
    (OUT / "b2_comparator_excluded_frame.json").write_text(json.dumps(out_a, indent=2))
    print(f"[a] comparator-excluded: AUROC={obs:.4f}, pair-clustered CI "
          f"[{boot['ci95_percentile'][0]:.4f}, {boot['ci95_percentile'][1]:.4f}], "
          f"perm p={perm['empirical_p']:.5f} ({len(sub)} entries, {int(labels.sum())} pos, "
          f"{int(pair_label.sum())} pos pairs)")

    # ── (b) label-subset AUROC overview ───────────────────────────
    def row(key, label, curated_pairs):
        f = hm[key]
        def r4(x):
            return round(x, 4) if isinstance(x, (int, float)) else None
        ci = f.get("auroc_ci95")
        return {"subset": label, "curated_pairs_in_label_set": curated_pairs,
                "n_entries": f["n_entries"], "n_positive_entries": f["n_positives"],
                "auroc_signed_dd": r4(f["auroc"]), "auprc": r4(f["auprc"]),
                "auroc_abs_dd_sensitivity": r4(f["auroc_abs_dd_sensitivity"]),
                "entry_bootstrap_ci95": [r4(x) for x in ci] if ci else None,
                "note": ("AUROC not estimable (0 positive entries in frame)"
                         if f["n_positives"] == 0 else None),
                "source": "output/headline_metrics.json"}
    overview = pd.DataFrame([
        row("lineage_full", "Full 12-pair curated set (secondary)", "12"),
        row("lineage_primary", "Sequence paralogs only (comparator-excluded, 10 pairs)", "10"),
        row("lineage_tier_ab", "PRIMARY Tier A u B (5 pairs)", "5"),
        row("lineage_leave_out_depmap_era", "Leave-out DepMap-era (10 pairs)", "10"),
        row("lineage_pre_depmap_only", "Pre-DepMap evidence only (8 pairs)", "8"),
        row("lineage_tier_a", "Tier A only (3 pairs)", "3"),
        row("lineage_tier_b", "Tier B only (2 pairs)", "2"),
        row("lineage_full_direction_strict", "Full, direction-strict (EP300->CREBBP relabelled)", "12"),
    ])
    overview.to_csv(OUT / "b2_subset_auroc_overview.csv", index=False)
    print(f"[b] subset overview written ({len(overview)} frames)")

    # ── (c) 8 positive entries composition ────────────────────────
    pos = df[df["is_known_paralog_sl"]].copy()
    comp_rows = []
    for r in pos.itertuples():
        ukey = tuple(sorted([r.driver_gene, r.paralog_gene]))
        # tier keyed by the validated direction in Table S3 where present
        tkey = (r.driver_gene, r.paralog_gene) if (r.driver_gene, r.paralog_gene) in TIER_OF else ukey
        comp_rows.append({
            "driver_gene": r.driver_gene, "paralog_gene": r.paralog_gene,
            "cancer_type": r.cancer_type,
            "tier": TIER_OF.get(tkey, "Comparator"),
            "evidence_source": KEY_REF_OF.get(tkey, KEY_REF_OF.get(ukey)),
            "direct_sl_evidence": DIRECT_SL_OF.get(tkey, DIRECT_SL_OF.get(ukey)),
            "signed_dd": round(float(r.dependency_dd), 4),
            "hedges_g": round(float(r.hedges_g), 3),
            "n_mut": int(r.n_mut), "n_wt": int(r.n_wt),
            "composite_score": round(float(r.composite_score), 3),
        })
    comp = pd.DataFrame(comp_rows).sort_values(
        ["tier", "driver_gene", "cancer_type"])
    comp.to_csv(OUT / "b2_positive_entries_composition.csv", index=False)
    print(f"[c] positive-entry composition written ({len(comp)} entries, "
          f"{comp.groupby(['driver_gene','paralog_gene']).ngroups} pairs)")

    # ── (d) STK11->SIK1 label-identity audit ──────────────────────
    vr = json.loads((ROOT / "output/validation_report.json").read_text())
    audit = {
        "pair": "STK11->SIK1",
        "manuscript_text_identity": {
            "identity": "mechanistic comparator / specificity reference, NOT a sequence "
                        "paralog, NOT a benchmark positive",
            "evidence": [
                "manuscript.tex Results: 'two serve as mechanistic comparators "
                "(BRCA1<->BRCA2, STK11->SIK1)'",
                "manuscript.tex Methods: 'not sequence paralogs and serve as "
                "specificity references only'",
                "supplementary Table S3: Tier='Comparator', Inclusion='Comparator', "
                "Direct_SL='Pathway axis', Key_Ref=Hollstein 2019",
            ],
        },
        "evaluation_code_identity": {
            "identity": "positive (is_known_paralog_sl = True) in every frame driven by "
                        "config.KNOWN_PARALOG_SL",
            "evidence": [
                "config.py: ('STK11','SIK1') is a member of KNOWN_PARALOG_SL",
                "TableS2_FullResults.tsv: STK11->SIK1 Cervical entry has "
                "is_known_paralog_sl=True (1 of the 8 positive entries; signed DD=0.0363)",
                "output/analysis_summary.txt: STK11->SIK1 listed with the star 'Known' flag",
                "supplementary Table S9 (DWS/composite): STK11->SIK1 row carries the star flag",
                f"validation_report.json per-pair frame: {vr['negative_control']['n_known']} "
                "positive pairs of 72 include STK11->SIK1 (labels derive from KNOWN_PARALOG_SL)",
            ],
        },
        "numerical_consequences": [
            "Headline full-frame AUROC 0.629 (110 entries, 8 positives) counts STK11->SIK1 "
            "as a positive although the text designates it a comparator.",
            "comparator-excluded frame removes the ordered tuples (STK11->SIK1 and "
            "BRCA1->BRCA2): 108 entries, 6 positives, AUROC 0.560.",
            "The ordered-tuple exclusion leaves BRCA2->BRCA1 (reverse direction of the "
            "other comparator pair) inside the comparator-excluded frame AS A POSITIVE; "
            "the Tier A u B frame excludes comparators as UNORDERED pairs and is clean.",
        ],
        "options_for_repair": [
            "Option 1: remove STK11->SIK1 (and both BRCA1<->BRCA2 directions) from "
            "KNOWN_PARALOG_SL so evaluation labels match the text (frames recompute).",
            "Option 2: keep labels, state explicitly in text/Table S3 that comparators "
            "are scored as positives in the full and per-pair frames.",
        ],
    }
    (OUT / "b2_stk11_sik1_label_audit.json").write_text(json.dumps(audit, indent=2))
    print("[d] STK11->SIK1 label audit written")

    print("\nDone. Files in output/revision_stage4/: "
          "b2_comparator_excluded_frame.json, b2_subset_auroc_overview.csv, "
          "b2_positive_entries_composition.csv, b2_stk11_sik1_label_audit.json")


if __name__ == "__main__":
    main()
