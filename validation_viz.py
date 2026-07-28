"""
Statistical validation suite for the paralog-SL pipeline (used by main.py).

NOTE (2026-07-26): the legacy standalone figure generators
(fig1_schema, fig2_benchmark, fig3_cross_cancer, fig4_validation,
fig5_survival, fig6_proteomics) were REMOVED. They hardcoded analysis
values — including fabricated TCGA survival hazard ratios, a fabricated
CNV R2 panel (annotated in-source as "estimated R2 values"), a
cross-cancer transfer matrix with no source artifact, and
survival/co-occurrence numbers that contradicted the recomputed pipeline
outputs. Manuscript figures are produced exclusively by the
artifact-driven R scripts (R_fig1.R ... R_figS10.R), which read
single-source-of-truth outputs and fail loudly (stop()) when an artifact
is missing. Simulated, random, or hardcoded data are forbidden in this
repository.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


# ═══════════════════════════════════════════════════════════════
# Validation functions (used by main.py)
# ═══════════════════════════════════════════════════════════════

def run_full_validation(results, n_permutations: int = 10000,
                        n_bootstrap: int = 1000, seed: int = 42):
    """
    Run the validation suite on paralog-SL results and return a metrics dict.

    Frameworks (manuscript "Evaluation frameworks" paragraph):
      * per-pair  — unique driver->paralog pairs in the three gynecological
        lineages (Ovarian/Endometrial/Cervical), scored by max |DD| across
        lineages. This reproduces the 77-pair / 8-positive framework cited in
        the manuscript (observed AUROC 0.6685 on the frozen artifact).
        Bootstrap and permutation analyses use this framework.
      * lineage-level — each driver x paralog x lineage entry separately
        (206 entries, 11 positives; AUROC 0.794).

    The pre-2026-07-25 output/validation_report.json (cited by the manuscript
    and by R_figS8.R) was produced by a no-longer-present script whose
    label-null had mean 0.58 — inconsistent with a true label shuffle (which
    must have mean 0.5). That historical version is preserved under
    output/backup_prerun_20260725/; main.py now regenerates
    output/validation_report.json reproducibly with a seeded, correct
    label-shuffle null, and writes the raw null to
    output/permutation_null_10000.csv for figure scripts.
    """
    rng = np.random.default_rng(seed)

    # ── Per-pair framework: gyn3 lineages, max signed DD across lineages ──
    gyn = results[results["cancer_type"].isin(["Ovarian", "Endometrial", "Cervical"])]
    pp = (gyn.groupby(["driver_gene", "paralog_gene"])
             .agg(score=("dependency_dd", "max"),
                  known=("is_known_paralog_sl", "max"))
             .reset_index())
    yt = pp["known"].astype(int).values
    ys = pp["score"].fillna(0).values
    n_known = int(yt.sum())
    n_total = len(pp)
    dd_auroc = roc_auc_score(yt, ys) if n_known >= 2 else float("nan")

    # Negative control: seeded label shuffle (null mean is ~0.5 by construction)
    null_aurocs = np.array([
        roc_auc_score(rng.permutation(yt), ys) for _ in range(n_permutations)
    ])
    null_mean = float(np.mean(null_aurocs))
    null_std = float(np.std(null_aurocs))
    emp_p = float((np.sum(null_aurocs >= dd_auroc) + 1) / (len(null_aurocs) + 1))

    # Bootstrap CI on the per-pair frame
    bs_aurocs = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n_total, n_total)
        if yt[idx].sum() >= 2:
            bs_aurocs.append(roc_auc_score(yt[idx], ys[idx]))
    bs_mean = float(np.mean(bs_aurocs)) if bs_aurocs else dd_auroc
    bs_ci_low = float(np.percentile(bs_aurocs, 2.5)) if bs_aurocs else 0.0
    bs_ci_high = float(np.percentile(bs_aurocs, 97.5)) if bs_aurocs else 1.0

    # ── Lineage-level frame (gyn3 entries; manuscript's "cancer-type-specific"
    # evaluation: 118 driver x paralog x lineage entries, 11 positives) ──
    yt_lin = gyn["is_known_paralog_sl"].astype(int).values
    ys_lin = gyn["dependency_dd"].fillna(0).values
    nk_lin = int(yt_lin.sum())
    lin_dd_auroc = roc_auc_score(yt_lin, ys_lin) if nk_lin >= 2 else float("nan")
    ys_comp = gyn.get("composite_score", pd.Series(0, index=gyn.index)).fillna(0).values
    lin_comp_auroc = roc_auc_score(yt_lin, ys_comp) if nk_lin >= 2 else float("nan")

    component_metrics = {}
    if "delta_expression" in results.columns:
        # component "expression_only" follows the historical artifact:
        # |delta_expression| AUROC over the full all-lineage frame
        yt_all = results["is_known_paralog_sl"].astype(int).values
        ys_expr = results["delta_expression"].abs().fillna(0).values
        try:
            component_metrics["expression_only"] = (
                roc_auc_score(yt_all, ys_expr) if yt_all.sum() >= 2 else float("nan"))
        except Exception:
            component_metrics["expression_only"] = float("nan")

    return {
        "framework": "per_pair: gyn3 unique pairs, score = max signed DD across lineages",
        "note": ("Reproducible companion to the frozen historical "
                 "output/validation_report.json; label-shuffle null is seeded "
                 "and has mean ~0.5 by construction."),
        "negative_control": {
            "observed_auroc": dd_auroc,
            "null_auroc_mean": null_mean,
            "null_auroc_std": null_std,
            "empirical_p_value": emp_p,
            "n_known": str(n_known),
            "n_total": n_total,
            "n_permutations": n_permutations,
            "seed": seed,
        },
        "component_decomposition": component_metrics,
        "bootstrap": {
            "auroc_mean": bs_mean,
            "auroc_ci_low": bs_ci_low,
            "auroc_ci_high": bs_ci_high,
            "n_bootstrap": n_bootstrap,
        },
        # Raw bootstrap resample AUROCs — popped by main.py and written to
        # output/bootstrap_perpair_1000.csv for the Fig. S8a histogram
        # (replaces the old rnorm-simulated shape with the real draws).
        "bootstrap_distribution": bs_aurocs,
        "lineage_level": {
            "frame": "gyn3 (Ovarian/Endometrial/Cervical) lineage-level entries",
            "dd_auroc": lin_dd_auroc,
            "composite_auroc": lin_comp_auroc,
            "n_entries": int(len(gyn)),
            "n_positives": nk_lin,
        },
        "null_distribution": null_aurocs,
    }


def cross_cancer_validation(all_results):
    """
    Cross-cancer paralog-SL validation.
    all_results: dict of {cancer_type: results_dataframe}
    """
    cancers = list(all_results.keys())
    if len(cancers) < 2:
        return

    print(f"\n  Cross-cancer validation ({len(cancers)} cancer types):")
    print(f"  {'Cancer':15s} {'Pairs':>6s} {'Known':>6s} {'DD AUROC':>9s}")
    print("  " + "-" * 42)

    summary = []
    for ct in cancers:
        r = all_results[ct]
        yt = r["is_known_paralog_sl"].astype(int).values
        ys = r["dependency_dd"].fillna(0).values
        nk = int(yt.sum())
        auc = roc_auc_score(yt, ys) if nk >= 2 else float("nan")
        astr = f"{auc:.3f}" if not np.isnan(auc) else "N/A"
        print(f"  {ct:15s} {len(r):>6d} {nk:>6d} {astr:>9s}")
        summary.append({"cancer": ct, "n_pairs": len(r), "n_known": nk, "dd_auroc": auc})

    if summary:
        pd.DataFrame(summary).to_csv("output/cross_cancer_summary.csv", index=False)


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    raise SystemExit(
        "validation_viz.py no longer generates figures (the legacy figure "
        "functions were removed because they hardcoded fabricated values). "
        "Use the artifact-driven R figure scripts (R_fig1.R ... R_figS10.R)."
    )
