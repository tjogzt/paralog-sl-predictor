"""
DWS Robustness: ranking sensitivity + bootstrap confidence intervals
====================================================================
Round-6 review items (therapeutic-window module), updated for the
round-7 signed-DWS decision:

  The PRIMARY DWS is now the signed formula
      DWS = max(DD, 0) / max(|mean Chronos|, pan-essential fraction, 0.01)
  so that only compensation-direction dependency (DD > 0, WT − MUT as in
  manuscript Eq. 1, pcs.py, and paralogSL::compute_dd) contributes. The
  pre-revision |DD| formula is retained as an explicit sensitivity
  analysis in both parts below.

  Part A — deterministic ranking sensitivity. Recomputes per-context DWS
  under alternative denominator components (|mean Chronos| only,
  pan-essential fraction only, mean of the two), alternative floors
  (0.001 / 0.05), a 1%/99% winsorized denominator, removal of
  floor-dominated pairs (any context whose denominator touches the 0.01
  floor), and the alternative |DD| numerator, then compares the
  resulting pair rankings (mean DWS across contexts) against the
  production signed ranking with Spearman rho and top-5 overlap.
  HIGH_SELECTIVITY counts are also reported for selectivity thresholds
  0.10 / 0.15 / 0.20.

  Part B — bootstrap 95% CIs for mean DWS (signed primary AND |DD|
  sensitivity) and mean selectivity per pair, plus a bootstrap 95% CI
  for each pair's RANK. Cell lines are resampled with replacement within
  the mutant, wild-type, and all-lines groups of every evaluated pair x
  context row (1,000 iterations, seed 42); per-iteration values are
  averaged across a pair's evaluable contexts to form the bootstrap
  distribution of the pair-level means reported in Table S5, and the
  full pair set is re-ranked in every iteration to form the per-pair
  rank distribution.

Inputs:  output/therapeutic_window_all_results.csv (production rows),
         DepMap dependency/mutation/expression data via data_loader
         (same driver-rule mutation sets as therapeutic_window.py).
Output:  output/dws_robustness.json

No simulated data: every number derives from the DepMap 26Q1 dependency
matrix and the production mutation sets.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL, GYN_CANCER_TYPES, MIN_MUT_SAMPLES
from data_loader import (
    load_dependency, load_expression, load_models,
    load_mutations, build_mutation_matrix,
)
from therapeutic_window import DRIVER_PARALOG_PAIRS, THERAPEUTIC_CANCERS, CERES_ESSENTIAL_THRESHOLD

BOOT_ITERS = 1000
BOOT_SEED = 42
FLOOR = 0.01


# ═══════════════════════════════════════════════════════════════
# Part A — deterministic sensitivity from production rows
# ═══════════════════════════════════════════════════════════════
def sensitivity(df: pd.DataFrame) -> dict:
    df = df.copy()
    df["dd_pos"] = df["dd"].clip(lower=0.0)   # signed numerator: max(DD, 0)
    df["abs_dd"] = df["dd"].abs()             # |DD| sensitivity numerator
    df["abs_mu"] = df["paralog_mean_ceres"].abs()
    df["f"] = df["paralog_pan_essential_frac"]

    variants = {
        "base signed max(DD,0) (max of |mu|, f, 0.01)": lambda r: r["dd_pos"] / np.maximum.reduce([r["abs_mu"], r["f"], np.full(len(r), 0.01)]),
        "denominator |mu| only": lambda r: r["dd_pos"] / np.maximum(r["abs_mu"], 0.01),
        "denominator f only": lambda r: r["dd_pos"] / np.maximum(r["f"], 0.01),
        "denominator mean(|mu|, f)": lambda r: r["dd_pos"] / np.maximum((r["abs_mu"] + r["f"]) / 2.0, 0.01),
        "floor 0.001": lambda r: r["dd_pos"] / np.maximum.reduce([r["abs_mu"], r["f"], np.full(len(r), 0.001)]),
        "floor 0.05": lambda r: r["dd_pos"] / np.maximum.reduce([r["abs_mu"], r["f"], np.full(len(r), 0.05)]),
        "numerator |DD| (pre-revision formula, sensitivity)": lambda r: r["abs_dd"] / np.maximum.reduce([r["abs_mu"], r["f"], np.full(len(r), 0.01)]),
    }

    def pair_table(dws: pd.Series) -> pd.DataFrame:
        t = df.assign(dws_v=dws)
        g = t.groupby(["driver", "paralog"], as_index=False).agg(
            mean_dws=("dws_v", "mean"),
            mean_sel=("selectivity", "mean"),
            mean_pan=("f", "mean"))
        g["pair"] = g["driver"] + "->" + g["paralog"]
        return g.sort_values("mean_dws", ascending=False).reset_index(drop=True)

    # validate the signed base variant against the production column
    base_check = variants["base signed max(DD,0) (max of |mu|, f, 0.01)"](df)
    max_diff = float(np.max(np.abs(base_check - df["dws_signed"])))
    # validate the |DD| sensitivity variant against the production column
    abs_check = variants["numerator |DD| (pre-revision formula, sensitivity)"](df)
    max_diff_abs = float(np.max(np.abs(abs_check - df["therapeutic_index"])))

    base_tab = pair_table(df["dws_signed"])
    base_rank = base_tab["pair"].tolist()

    def classify(tab, sel_thr):
        def cls(r):
            if r["mean_pan"] > 0.5:
                return "PAN_ESSENTIAL"
            if r["mean_sel"] > sel_thr and r["mean_dws"] > 1.0:
                return "HIGH_SELECTIVITY"
            if r["mean_sel"] > 0:
                return "MODERATE"
            return "LOW_SELECTIVITY"
        return tab.apply(cls, axis=1)

    base_cls = classify(base_tab, 0.15)

    out = []
    pos_base = pd.Series({p: i for i, p in enumerate(base_rank)})
    for name, fn in variants.items():
        tab = pair_table(fn(df))
        pos_var = pd.Series({p: i for i, p in enumerate(tab["pair"].tolist())})
        rho = float(pos_base.corr(pos_var, method="spearman"))
        top5 = tab["pair"].head(5).tolist()
        entry = {
            "variant": name,
            "spearman_rho_vs_base": round(rho, 4),
            "top5_overlap_with_base": len(set(top5) & set(base_rank[:5])),
            "top5": top5,
        }
        if name.startswith("base"):
            entry["max_abs_diff_vs_production_column"] = max_diff
        if name.startswith("numerator |DD|"):
            entry["max_abs_diff_vs_production_column"] = max_diff_abs
        out.append(entry)

    # (a) Winsorized denominator: clip the baseline-essentiality denominator
    # max(|mu|, f, 0.01) at its 1st/99th percentiles across all pair x
    # context rows, limiting the influence of extreme denominators.
    denom = np.maximum.reduce([df["abs_mu"], df["f"], np.full(len(df), FLOOR)])
    w_lo, w_hi = np.percentile(denom, [1, 99])
    tab_w = pair_table(df["dd_pos"] / np.clip(denom, w_lo, w_hi))
    pos_w = pd.Series({p: i for i, p in enumerate(tab_w["pair"].tolist())})
    out.append({
        "variant": "denominator winsorized 1%/99%",
        "spearman_rho_vs_base": round(float(pos_base.corr(pos_w, method="spearman")), 4),
        "top5_overlap_with_base": len(set(tab_w["pair"].head(5)) & set(base_rank[:5])),
        "top5": tab_w["pair"].head(5).tolist(),
        "winsorize_bounds": [float(w_lo), float(w_hi)],
    })

    # (b) Floor-dominated pairs removed: drop any pair whose denominator
    # touches the 0.01 floor in at least one evaluable context (the floor,
    # not the data, sets that context's denominator), then re-rank the
    # remaining pairs. Spearman rho is computed on the retained pairs.
    floor_active = np.maximum(df["abs_mu"], df["f"]) < FLOOR
    fa_by_pair = df.assign(fa=floor_active).groupby(["driver", "paralog"])["fa"].any()
    excluded_pairs = sorted(fa_by_pair[fa_by_pair].index.tolist())
    keep = np.array([(d, p) not in excluded_pairs
                     for d, p in zip(df["driver"], df["paralog"])])
    sub = df.loc[keep].assign(dws_v=(df["dd_pos"][keep] / denom[keep]))
    tab_nf = (sub.groupby(["driver", "paralog"], as_index=False)
                 .agg(mean_dws=("dws_v", "mean"))
                 .assign(pair=lambda t: t["driver"] + "->" + t["paralog"])
                 .sort_values("mean_dws", ascending=False).reset_index(drop=True))
    excluded_names = {f"{d}->{p}" for d, p in excluded_pairs}
    pos_base_kept = pd.Series(
        {p: i for i, p in enumerate([p for p in base_rank if p not in excluded_names])})
    pos_nf = pd.Series({p: i for i, p in enumerate(tab_nf["pair"].tolist())})
    top5_nf = tab_nf["pair"].head(5).tolist()
    out.append({
        "variant": "floor-dominated pairs removed (denominator at 0.01 floor in >=1 context)",
        "spearman_rho_vs_base": round(float(pos_base_kept.corr(pos_nf, method="spearman")), 4),
        "top5_overlap_with_base": len(set(top5_nf) & set(base_rank[:5])),
        "top5": top5_nf,
        "n_pairs_excluded": len(excluded_pairs),
        "excluded_pairs": sorted(excluded_names),
    })

    for thr in (0.10, 0.15, 0.20):
        cls_v = classify(base_tab, thr)
        out.append({
            "variant": f"HIGH_SELECTIVITY selectivity threshold {thr}",
            "n_high_selectivity": int((cls_v == "HIGH_SELECTIVITY").sum()),
            "high_selectivity_pairs": base_tab.loc[cls_v == "HIGH_SELECTIVITY", "pair"].tolist(),
            "classification_flips_vs_thr_0.15": int((cls_v != base_cls).sum()),
        })
    return {"checks": out, "base_pair_order": base_rank,
            "base_validation_max_abs_diff": max_diff,
            "abs_validation_max_abs_diff": max_diff_abs}


# ═══════════════════════════════════════════════════════════════
# Part B — bootstrap CIs (resample cell lines within groups)
# ═══════════════════════════════════════════════════════════════
def bootstrap(df: pd.DataFrame) -> list:
    print("  loading DepMap data for bootstrap ...")
    dep = load_dependency()
    expr = load_expression()
    models = load_models()
    mutations = load_mutations()

    drivers = sorted({d for d, _ in DRIVER_PARALOG_PAIRS})
    mut_matrix = build_mutation_matrix(mutations, dep.index.tolist(), drivers)
    mut_sets = {g: set(mut_matrix.index[mut_matrix[g] == 1]) for g in mut_matrix.columns}

    id_pos = {c: i for i, c in enumerate(dep.index)}
    rng = np.random.default_rng(BOOT_SEED)

    # Rebuild the exact pair x context sets used by the production table
    rows = []
    for context_name, disease_patterns in THERAPEUTIC_CANCERS.items():
        if disease_patterns is None:
            model_subset = models.copy()
        else:
            pat = "|".join(disease_patterns)
            mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
            model_subset = models[mask].copy()
        valid_ids = [c for c in model_subset["DepMap_ID"].tolist()
                     if c in dep.index and c in expr.index]
        if len(valid_ids) < 20:
            continue
        for driver, paralog in DRIVER_PARALOG_PAIRS:
            if driver not in dep.columns or paralog not in dep.columns:
                continue
            driver_mut_set = mut_sets.get(driver, set())
            mut_ids = [c for c in valid_ids if c in driver_mut_set]
            wt_ids = [c for c in valid_ids if c not in driver_mut_set]
            if len(mut_ids) < MIN_MUT_SAMPLES:
                continue
            hit = df[(df.driver == driver) & (df.paralog == paralog) & (df.context == context_name)]
            if len(hit) != 1:
                continue  # production pipeline dropped this row for another reason
            rows.append(dict(driver=driver, paralog=paralog, context=context_name,
                             mut=np.array([id_pos[c] for c in mut_ids]),
                             wt=np.array([id_pos[c] for c in wt_ids]),
                             all=np.array([id_pos[c] for c in valid_ids]),
                             pvals=dep[paralog].to_numpy()))

    # sanity: coverage of production rows
    print(f"  bootstrap rows: {len(rows)} (production rows: {len(df)})")

    # Per-row bootstrap distributions. The signed (primary) and |DD|
    # (sensitivity) DWS are computed from the SAME resamples; the rng draw
    # order (mut, WT, all) is unchanged from the pre-revision script, so
    # the |DD| sensitivity distribution reproduces the old |DD| run
    # exactly under the same seed.
    pair_boot = {}
    for r in rows:
        pv = r["pvals"]
        n_mut, n_wt, n_all = len(r["mut"]), len(r["wt"]), len(r["all"])
        dws_b = np.empty(BOOT_ITERS)      # signed max(DD,0) — primary
        dws_abs_b = np.empty(BOOT_ITERS)  # |DD| — sensitivity
        sel_b = np.empty(BOOT_ITERS)
        for i in range(BOOT_ITERS):
            mi = rng.integers(0, n_mut, n_mut)
            wi = rng.integers(0, n_wt, n_wt)
            ai = rng.integers(0, n_all, n_all)
            mv, wv, av = pv[r["mut"][mi]], pv[r["wt"][wi]], pv[r["all"][ai]]
            dd = wv.mean() - mv.mean()
            pan_mean = av.mean()
            f_pan = float((av < CERES_ESSENTIAL_THRESHOLD).mean())
            denom = max(abs(pan_mean), f_pan, FLOOR)
            dws_b[i] = max(dd, 0.0) / denom
            dws_abs_b[i] = abs(dd) / denom
            sel_b[i] = float((mv < CERES_ESSENTIAL_THRESHOLD).mean()
                             - (wv < CERES_ESSENTIAL_THRESHOLD).mean())
        key = (r["driver"], r["paralog"])
        pair_boot.setdefault(key, {"dws": [], "dws_abs": [], "sel": []})
        pair_boot[key]["dws"].append(dws_b)
        pair_boot[key]["dws_abs"].append(dws_abs_b)
        pair_boot[key]["sel"].append(sel_b)

    out = []
    pair_mean_dists = {}   # (driver, paralog) -> per-iteration pair mean signed DWS
    for (driver, paralog), bd in sorted(pair_boot.items()):
        dws_mat = np.vstack(bd["dws"])        # contexts x iters
        dws_abs_mat = np.vstack(bd["dws_abs"])
        sel_mat = np.vstack(bd["sel"])
        dws_mean_dist = dws_mat.mean(axis=0)
        dws_abs_mean_dist = dws_abs_mat.mean(axis=0)
        sel_mean_dist = sel_mat.mean(axis=0)
        pair_mean_dists[(driver, paralog)] = dws_mean_dist
        out.append({
            "driver": driver, "paralog": paralog,
            "n_contexts": dws_mat.shape[0],
            "mean_dws": float(dws_mat.mean()),
            "dws_ci95": [float(np.percentile(dws_mean_dist, 2.5)),
                         float(np.percentile(dws_mean_dist, 97.5))],
            "mean_dws_abs": float(dws_abs_mat.mean()),
            "dws_abs_ci95": [float(np.percentile(dws_abs_mean_dist, 2.5)),
                             float(np.percentile(dws_abs_mean_dist, 97.5))],
            "mean_selectivity": float(sel_mat.mean()),
            "selectivity_ci95": [float(np.percentile(sel_mean_dist, 2.5)),
                                 float(np.percentile(sel_mean_dist, 97.5))],
        })

    # Per-pair bootstrap RANK distribution: in every iteration the full
    # pair set is re-ranked by the pair-level mean signed DWS (rank 1 =
    # highest DWS); percentile CI over iterations per pair.
    keys = sorted(pair_mean_dists)
    M = np.vstack([pair_mean_dists[k] for k in keys])   # pairs x iters
    ranks = np.empty_like(M)
    for i in range(BOOT_ITERS):
        ranks[:, i] = pd.Series(M[:, i]).rank(ascending=False, method="average").values
    rank_out = []
    for j, (driver, paralog) in enumerate(keys):
        rank_out.append({
            "driver": driver, "paralog": paralog,
            "rank_median": float(np.median(ranks[j])),
            "rank_ci95": [float(np.percentile(ranks[j], 2.5)),
                          float(np.percentile(ranks[j], 97.5))],
        })
    return out, rank_out


def main():
    print("=" * 70)
    print("  DWS robustness: sensitivity + bootstrap (seed", BOOT_SEED, ")")
    print("=" * 70)
    df = pd.read_csv(OUTPUT_DIR / "therapeutic_window_all_results.csv")

    sens = sensitivity(df)
    print(f"  base-variant validation: max |recomputed signed - production dws_signed| = "
          f"{sens['base_validation_max_abs_diff']:.2e}")
    print(f"  |DD|-variant validation: max |recomputed |DD| - production therapeutic_index| = "
          f"{sens['abs_validation_max_abs_diff']:.2e}")
    for c in sens["checks"]:
        if "spearman_rho_vs_base" in c:
            print(f"  {c['variant']:55s} rho={c['spearman_rho_vs_base']:.3f} "
                  f"top5-overlap={c['top5_overlap_with_base']}")
        else:
            print(f"  {c['variant']:55s} n_HS={c['n_high_selectivity']} "
                  f"flips={c['classification_flips_vs_thr_0.15']}")

    boot, rank_ci = bootstrap(df)
    rank_by_pair = {(r["driver"], r["paralog"]): r for r in rank_ci}
    for b in boot:
        if (b["driver"], b["paralog"]) in [("ARID1A", "ARID1B"), ("SMARCA4", "SMARCA2"), ("NF1", "RASA2")]:
            rk = rank_by_pair[(b["driver"], b["paralog"])]
            print(f"  {b['driver']}->{b['paralog']}: DWS {b['mean_dws']:.2f} "
                  f"[{b['dws_ci95'][0]:.2f}, {b['dws_ci95'][1]:.2f}]  "
                  f"|DD| sens {b['mean_dws_abs']:.2f} "
                  f"[{b['dws_abs_ci95'][0]:.2f}, {b['dws_abs_ci95'][1]:.2f}]  "
                  f"rank med {rk['rank_median']:.0f} "
                  f"[{rk['rank_ci95'][0]:.0f}, {rk['rank_ci95'][1]:.0f}]  "
                  f"sel {b['mean_selectivity']:+.2f} "
                  f"[{b['selectivity_ci95'][0]:+.2f}, {b['selectivity_ci95'][1]:+.2f}]")

    payload = {
        "method": {
            "dws_formula": "PRIMARY signed: max(DD,0) / max(|mean Chronos|, "
                           "pan-essential fraction, 0.01); DD = WT − MUT "
                           "(manuscript Eq. 1, pcs.py, paralogSL::compute_dd). "
                           "|DD| numerator retained as sensitivity only.",
            "sensitivity": "deterministic recomputation of DWS variants from "
                           "output/therapeutic_window_all_results.csv; rankings = mean DWS "
                           "across contexts per pair; comparison by Spearman rho and top-5 overlap",
            "bootstrap": f"{BOOT_ITERS} stratified resamples (mut/WT/all groups) per "
                         f"pair x context, seed {BOOT_SEED}; per-iteration means averaged "
                         "across a pair's evaluable contexts; percentile 95% CI; "
                         "signed and |DD| DWS share the same resamples",
            "rank_bootstrap": "full pair set re-ranked by pair-level mean signed DWS "
                              "in every bootstrap iteration; per-pair percentile 95% CI "
                              "of the rank (rank 1 = highest DWS)",
            "essential_threshold": CERES_ESSENTIAL_THRESHOLD,
            "floor": FLOOR,
        },
        "sensitivity": sens["checks"],
        "bootstrap": boot,
        "bootstrap_rank_ci95": rank_ci,
    }
    out_path = OUTPUT_DIR / "dws_robustness.json"
    out_path.write_text(json.dumps(payload, indent=1))
    print(f"  written: {out_path}")


if __name__ == "__main__":
    main()
