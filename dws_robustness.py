"""
DWS Robustness: ranking sensitivity + bootstrap confidence intervals
====================================================================
Round-6 review items (therapeutic-window module):

  Part A — deterministic ranking sensitivity. Recomputes per-context DWS
  under alternative denominator components (|mean Chronos| only,
  pan-essential fraction only, mean of the two) and alternative floors
  (0.001 / 0.05), then compares the resulting pair rankings (mean DWS
  across contexts) against the production ranking with Spearman rho and
  top-5 overlap. HIGH_SELECTIVITY counts are also reported for
  selectivity thresholds 0.10 / 0.15 / 0.20.

  Part B — bootstrap 95% CIs for mean DWS and mean selectivity per pair.
  Cell lines are resampled with replacement within the mutant,
  wild-type, and all-lines groups of every evaluated pair x context row
  (1,000 iterations, seed 42); per-iteration values are averaged across
  a pair's evaluable contexts to form the bootstrap distribution of the
  pair-level means reported in Table S6.

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
    df["abs_dd"] = df["dd"].abs()
    df["abs_mu"] = df["paralog_mean_ceres"].abs()
    df["f"] = df["paralog_pan_essential_frac"]

    variants = {
        "base (max of |mu|, f, 0.01)": lambda r: r["abs_dd"] / np.maximum.reduce([r["abs_mu"], r["f"], np.full(len(r), 0.01)]),
        "denominator |mu| only": lambda r: r["abs_dd"] / np.maximum(r["abs_mu"], 0.01),
        "denominator f only": lambda r: r["abs_dd"] / np.maximum(r["f"], 0.01),
        "denominator mean(|mu|, f)": lambda r: r["abs_dd"] / np.maximum((r["abs_mu"] + r["f"]) / 2.0, 0.01),
        "floor 0.001": lambda r: r["abs_dd"] / np.maximum.reduce([r["abs_mu"], r["f"], np.full(len(r), 0.001)]),
        "floor 0.05": lambda r: r["abs_dd"] / np.maximum.reduce([r["abs_mu"], r["f"], np.full(len(r), 0.05)]),
    }

    def pair_table(dws: pd.Series) -> pd.DataFrame:
        t = df.assign(dws_v=dws)
        g = t.groupby(["driver", "paralog"], as_index=False).agg(
            mean_dws=("dws_v", "mean"),
            mean_sel=("selectivity", "mean"),
            mean_pan=("f", "mean"))
        g["pair"] = g["driver"] + "->" + g["paralog"]
        return g.sort_values("mean_dws", ascending=False).reset_index(drop=True)

    # validate the base variant against the production column
    base_check = variants["base (max of |mu|, f, 0.01)"](df)
    max_diff = float(np.max(np.abs(base_check - df["therapeutic_index"])))

    base_tab = pair_table(df["therapeutic_index"])
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
        out.append(entry)

    for thr in (0.10, 0.15, 0.20):
        cls_v = classify(base_tab, thr)
        out.append({
            "variant": f"HIGH_SELECTIVITY selectivity threshold {thr}",
            "n_high_selectivity": int((cls_v == "HIGH_SELECTIVITY").sum()),
            "high_selectivity_pairs": base_tab.loc[cls_v == "HIGH_SELECTIVITY", "pair"].tolist(),
            "classification_flips_vs_thr_0.15": int((cls_v != base_cls).sum()),
        })
    return {"checks": out, "base_pair_order": base_rank,
            "base_validation_max_abs_diff": max_diff}


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

    # Per-row bootstrap distributions
    pair_boot = {}
    for r in rows:
        pv = r["pvals"]
        n_mut, n_wt, n_all = len(r["mut"]), len(r["wt"]), len(r["all"])
        dws_b = np.empty(BOOT_ITERS)
        sel_b = np.empty(BOOT_ITERS)
        for i in range(BOOT_ITERS):
            mi = rng.integers(0, n_mut, n_mut)
            wi = rng.integers(0, n_wt, n_wt)
            ai = rng.integers(0, n_all, n_all)
            mv, wv, av = pv[r["mut"][mi]], pv[r["wt"][wi]], pv[r["all"][ai]]
            dd = wv.mean() - mv.mean()
            pan_mean = av.mean()
            f_pan = float((av < CERES_ESSENTIAL_THRESHOLD).mean())
            dws_b[i] = abs(dd) / max(abs(pan_mean), f_pan, FLOOR)
            sel_b[i] = float((mv < CERES_ESSENTIAL_THRESHOLD).mean()
                             - (wv < CERES_ESSENTIAL_THRESHOLD).mean())
        key = (r["driver"], r["paralog"])
        pair_boot.setdefault(key, {"dws": [], "sel": []})
        pair_boot[key]["dws"].append(dws_b)
        pair_boot[key]["sel"].append(sel_b)

    out = []
    for (driver, paralog), bd in sorted(pair_boot.items()):
        dws_mat = np.vstack(bd["dws"])   # contexts x iters
        sel_mat = np.vstack(bd["sel"])
        dws_mean_dist = dws_mat.mean(axis=0)
        sel_mean_dist = sel_mat.mean(axis=0)
        out.append({
            "driver": driver, "paralog": paralog,
            "n_contexts": dws_mat.shape[0],
            "mean_dws": float(dws_mat.mean()),
            "dws_ci95": [float(np.percentile(dws_mean_dist, 2.5)),
                         float(np.percentile(dws_mean_dist, 97.5))],
            "mean_selectivity": float(sel_mat.mean()),
            "selectivity_ci95": [float(np.percentile(sel_mean_dist, 2.5)),
                                 float(np.percentile(sel_mean_dist, 97.5))],
        })
    return out


def main():
    print("=" * 70)
    print("  DWS robustness: sensitivity + bootstrap (seed", BOOT_SEED, ")")
    print("=" * 70)
    df = pd.read_csv(OUTPUT_DIR / "therapeutic_window_all_results.csv")

    sens = sensitivity(df)
    print(f"  base-variant validation: max |recomputed - production| = "
          f"{sens['base_validation_max_abs_diff']:.2e}")
    for c in sens["checks"]:
        if "spearman_rho_vs_base" in c:
            print(f"  {c['variant']:45s} rho={c['spearman_rho_vs_base']:.3f} "
                  f"top5-overlap={c['top5_overlap_with_base']}")
        else:
            print(f"  {c['variant']:45s} n_HS={c['n_high_selectivity']} "
                  f"flips={c['classification_flips_vs_thr_0.15']}")

    boot = bootstrap(df)
    for b in boot:
        if (b["driver"], b["paralog"]) in [("ARID1A", "ARID1B"), ("SMARCA4", "SMARCA2"), ("NF1", "RASA2")]:
            print(f"  {b['driver']}->{b['paralog']}: DWS {b['mean_dws']:.2f} "
                  f"[{b['dws_ci95'][0]:.2f}, {b['dws_ci95'][1]:.2f}]  "
                  f"sel {b['mean_selectivity']:+.2f} "
                  f"[{b['selectivity_ci95'][0]:+.2f}, {b['selectivity_ci95'][1]:+.2f}]")

    payload = {
        "method": {
            "sensitivity": "deterministic recomputation of DWS variants from "
                           "output/therapeutic_window_all_results.csv; rankings = mean DWS "
                           "across contexts per pair; comparison by Spearman rho and top-5 overlap",
            "bootstrap": f"{BOOT_ITERS} stratified resamples (mut/WT/all groups) per "
                         f"pair x context, seed {BOOT_SEED}; per-iteration means averaged "
                         "across a pair's evaluable contexts; percentile 95% CI",
            "essential_threshold": CERES_ESSENTIAL_THRESHOLD,
            "floor": FLOOR,
        },
        "sensitivity": sens["checks"],
        "bootstrap": boot,
    }
    out_path = OUTPUT_DIR / "dws_robustness.json"
    out_path.write_text(json.dumps(payload, indent=1))
    print(f"  written: {out_path}")


if __name__ == "__main__":
    main()
