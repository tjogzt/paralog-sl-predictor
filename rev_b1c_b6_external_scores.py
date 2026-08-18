#!/usr/bin/env python3
"""
rev_b1c_b6_external_scores.py  (Stage-4 revision, items B1c + B6)
================================================================
External-frame (Harle 2025 / Flister 2025 digenic screens) ranking-statistic
sensitivity (B1c) and lineage-matched Layer-1 sensitivity (B6).

Engine: mirrors external_holdout.maxabs_dd_scores (same >=3-mutant / >=3-WT
strata, same max-across-both-orientations selection) but additionally
computes, per (pair, lineage), the best |signed DD|, |Welch t| and
|Hedges' g| (each statistic independently maximized across the two
orientations within the lineage):

  dd  = mean(Chronos paralog | driver-mut) - mean(... | driver-WT)  [external
        sign convention, identical to external_holdout.py / in4mer dd_min3]
  t   = dd / sqrt(var_m/n_m + var_w/n_w)        (Welch, ddof=1 variances)
  g   = (dd / pooled_sd) * (1 - 3/(4N-9)), pooled_sd = sqrt((var_m+var_w)/2)

Fidelity: per-pair max|DD| across ALL lineages must equal the frozen
output/external_holdout_{harle,flister}_layer1.csv dd_min3 for every pair
(tol 1e-8); the run aborts otherwise.

B1c outputs: AUROC(max|DD|) vs AUROC(max|t|) vs AUROC(max|g|) on each
screen's Layer-1 evaluable frame + Spearman rank correlations + delta AUROC
(+ label-permutation p for each statistic, 10,000, seed 42).

B6 outputs: max|DD| recomputed using ONLY the DepMap lineages matching the
screened models:
  Harle  (melanoma/lung-NSCLC/pancreas) -> {"Melanoma",
         "Non-Small Cell Lung Cancer", "Pancreatic Adenocarcinoma"}
  Flister(NCI-H1299 lung + MDA-MB-231 breast) -> {"Non-Small Cell Lung
         Cancer", "Invasive Breast Carcinoma"}; per-model matched analyses
         (lung-only vs H1299 hits; breast-only vs MDAMB231 hits) included.

Outputs (output/revision_stage4/):
  b1c_harle_three_stats.csv / b1c_flister_three_stats.csv   per-pair scores
  b1c_external_three_stats_summary.json
  b6_lineage_matched_layer1.json
  b6_lineage_matched_scores.csv

Usage: python rev_b1c_b6_external_scores.py   (run from repo root)
"""

import json
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

from external_holdout import (  # noqa: E402
    auroc_with_permutation, harle_pair_frame, flister_pair_frame,
    load_harle, load_flister_s3, STATE_PKL, CACHE,
)
from scipy.stats import spearmanr  # noqa: E402

MIN_N = 3
SEED = 42
N_PERM = 10_000

HARLE_MATCHED = {"Melanoma", "Non-Small Cell Lung Cancer", "Pancreatic Adenocarcinoma"}
FLISTER_MATCHED = {"Non-Small Cell Lung Cancer", "Invasive Breast Carcinoma"}


def per_lineage_stat_scores(pairs, dep, mat, lin_map, min_n=MIN_N, tag=""):
    """For each unordered pair: {lineage: (best_dd, best_t, best_g)} with each
    statistic independently maximized (absolute value) across the two
    orientations within that lineage. Mirrors maxabs_dd_scores arm rules."""
    gene2col = {g: j for j, g in enumerate(dep.columns)}
    dep_np = dep.to_numpy(dtype=np.float64)
    lineages = sorted(set(lin_map.dropna()))
    lin_vals = lin_map.to_numpy()
    lin_idx = {lin: np.flatnonzero(lin_vals == lin) for lin in lineages}
    L = np.zeros((len(lineages), dep_np.shape[0]), dtype=np.uint8)
    for i, lin in enumerate(lineages):
        L[i, lin_idx[lin]] = 1

    partners_of = defaultdict(set)
    for a, b in pairs:
        partners_of[a].add(b)
        partners_of[b].add(a)

    best = defaultdict(dict)  # frozenset pair -> {lineage: [dd, t, g]}
    t0 = time.time()
    for gi, (drv, partners) in enumerate(sorted(partners_of.items())):
        if tag and (gi + 1) % 2000 == 0:
            print(f"      [{tag}] driver {gi + 1}/{len(partners_of)} ({time.time() - t0:.0f}s)")
        if drv not in mat.columns:
            continue
        pj = [(p, gene2col[p]) for p in partners if p in gene2col]
        if not pj:
            continue
        m = mat[drv].to_numpy(dtype=np.uint8)
        n_mut = L @ m
        n_wt = L.sum(axis=1) - n_mut
        ok = np.flatnonzero((n_mut >= min_n) & (n_wt >= min_n))
        if len(ok) == 0:
            continue
        cols = [j for _, j in pj]
        for i in ok:
            idx = lin_idx[lineages[i]]
            mm = m[idx].astype(bool)
            sub = dep_np[np.ix_(idx, cols)]
            ms, ws = sub[mm], sub[~mm]
            nvm = np.sum(~np.isnan(ms), axis=0)
            nvw = np.sum(~np.isnan(ws), axis=0)
            valid = (nvm >= min_n) & (nvw >= min_n)
            if not valid.any():
                continue
            sm = np.nansum(ms, axis=0); sw = np.nansum(ws, axis=0)
            s2m = np.nansum(ms * ms, axis=0); s2w = np.nansum(ws * ws, axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                mu_m = sm / nvm; mu_w = sw / nvw
                var_m = (s2m - sm * sm / nvm) / (nvm - 1)
                var_w = (s2w - sw * sw / nvw) / (nvw - 1)
                dd = mu_m - mu_w
                se = np.sqrt(var_m / nvm + var_w / nvw)
                t = np.where(se > 0, dd / se, np.nan)
                pooled = np.sqrt((var_m + var_w) / 2)
                ntot = nvm + nvw
                d = np.where(pooled > 0, dd / pooled, 0.0)
                g = np.where(ntot > 3, d * (1 - 3 / (4 * ntot - 9)), d)
            lin = lineages[i]
            for k, (p, _) in enumerate(pj):
                if not valid[k] or np.isnan(dd[k]):
                    continue
                key = frozenset((drv, p))
                cur = best[key].get(lin)
                if cur is None:
                    best[key][lin] = [float(dd[k]),
                                      float(t[k]) if not np.isnan(t[k]) else None,
                                      float(g[k])]
                else:
                    if abs(dd[k]) > abs(cur[0]):
                        cur[0] = float(dd[k])
                    if not np.isnan(t[k]) and (cur[1] is None or abs(t[k]) > abs(cur[1])):
                        cur[1] = float(t[k])
                    if abs(g[k]) > abs(cur[2]):
                        cur[2] = float(g[k])
    return best


def aggregate(best, pairs, lineages_subset=None):
    """Per-pair frame scores: max over lineages of |dd|, |t|, |g|."""
    out = {}
    for p in pairs:
        rec = best.get(frozenset(p), {})
        if lineages_subset is not None:
            rec = {k: v for k, v in rec.items() if k in lineages_subset}
        if not rec:
            out[p] = (None, None, None)
            continue
        dds = [v[0] for v in rec.values() if v[0] is not None]
        ts = [v[1] for v in rec.values() if v[1] is not None]
        gs = [v[2] for v in rec.values() if v[2] is not None]
        out[p] = (max(dds, key=abs) if dds else None,
                  max(ts, key=abs) if ts else None,
                  max(gs, key=abs) if gs else None)
    return out


def fidelity(frozen_csv, agg):
    fr = pd.read_csv(ROOT / "output" / frozen_csv)
    bad = 0
    for r in fr.itertuples():
        ref = r.dd_min3 if pd.notna(r.dd_min3) else None
        got = agg[(r.gene_a, r.gene_b)][0]
        if (ref is None) != (got is None):
            bad += 1
        elif ref is not None and abs(ref - got) > 1e-8:
            bad += 1
    if bad:
        raise RuntimeError(f"FIDELITY FAIL vs {frozen_csv}: {bad} mismatches")
    print(f"    fidelity vs {frozen_csv}: OK ({len(fr)} pairs)")


def three_stat_eval(labels, agg, pairs, permute=True):
    df = pd.DataFrame({
        "y": labels,
        "abs_dd": [abs(agg[p][0]) if agg[p][0] is not None else np.nan for p in pairs],
        "abs_t": [abs(agg[p][1]) if agg[p][1] is not None else np.nan for p in pairs],
        "abs_g": [abs(agg[p][2]) if agg[p][2] is not None else np.nan for p in pairs],
    })
    res = {"n_pairs": int(len(df))}
    ev = df.dropna(subset=["abs_dd"])
    res["n_scored"] = int(len(ev))
    res["n_hits"] = int(ev["y"].sum())
    for col in ("abs_dd", "abs_t", "abs_g"):
        e2 = ev.dropna(subset=[col])
        r = auroc_with_permutation(e2["y"].to_numpy(), e2[col].to_numpy()) if permute else \
            {"auroc": float(pd.Series(e2[col]).rank().corr())}
        res[col] = {"auroc": r["auroc"], "permutation_p": r.get("permutation_p"),
                    "n": int(len(e2)), "n_hits": int(e2["y"].sum())}
    # spearman among the three score vectors (dd-scored universe)
    res["spearman"] = {
        "abs_dd_vs_abs_t": float(spearmanr(ev["abs_dd"], ev["abs_t"], nan_policy="omit").statistic),
        "abs_dd_vs_abs_g": float(spearmanr(ev["abs_dd"], ev["abs_g"], nan_policy="omit").statistic),
        "abs_t_vs_abs_g": float(spearmanr(ev["abs_t"], ev["abs_g"], nan_policy="omit").statistic),
    }
    res["delta_auroc"] = {
        "t_minus_dd": res["abs_t"]["auroc"] - res["abs_dd"]["auroc"],
        "g_minus_dd": res["abs_g"]["auroc"] - res["abs_dd"]["auroc"],
    }
    return res, df


def main():
    print("=" * 72)
    print("  rev B1c + B6: external-frame three-statistic + lineage-matched Layer 1")
    print("=" * 72)
    st = pd.read_pickle(STATE_PKL)
    dep, mat, lin_map = st["dep"], st["mat"], st["lin_map"]
    print(f"  state loaded: dep {dep.shape}, mat {mat.shape}")

    harle_s5, _ = load_harle()
    fl_s3 = load_flister_s3()
    pairs_h, hit_h = harle_pair_frame(harle_s5)
    pairs_f, h1299, mda, union = flister_pair_frame(fl_s3)

    print("\n  scoring Harle pairs (three statistics, per lineage) ...")
    best_h = per_lineage_stat_scores(pairs_h, dep, mat, lin_map, tag="harle")
    agg_h = aggregate(best_h, pairs_h)
    fidelity("external_holdout_harle_layer1.csv", agg_h)

    print("\n  scoring Flister pairs (three statistics, per lineage) ...")
    best_f = per_lineage_stat_scores(pairs_f, dep, mat, lin_map, tag="flister")
    agg_f = aggregate(best_f, pairs_f)
    fidelity("external_holdout_flister_layer1.csv", agg_f)

    # ── B1c ───────────────────────────────────────────────────────
    print("\n[B1c] three-statistic AUROC on the external frames")
    res_h, dfh = three_stat_eval(hit_h, agg_h, pairs_h)
    res_f, dff = three_stat_eval(union, agg_f, pairs_f)
    # Flister per-model sensitivity
    res_f["per_model"] = {}
    for name, arr in (("hit_h1299", h1299), ("hit_mdamb231", mda)):
        r, _ = three_stat_eval(np.asarray(arr), agg_f, pairs_f)
        res_f["per_model"][name] = r

    dfh_out = pd.DataFrame({
        "pair": ["|".join(sorted(p)) for p in pairs_h],
        "gene_a": [p[0] for p in pairs_h], "gene_b": [p[1] for p in pairs_h],
        "hit_any_line": hit_h,
        "signed_dd_maxabs": [agg_h[p][0] for p in pairs_h],
        "welch_t_maxabs": [agg_h[p][1] for p in pairs_h],
        "hedges_g_maxabs": [agg_h[p][2] for p in pairs_h]})
    dfh_out.to_csv(OUT / "b1c_harle_three_stats.csv", index=False)
    dff_out = pd.DataFrame({
        "pair": [f"{p[0]}_{p[1]}" for p in pairs_f],
        "gene_a": [p[0] for p in pairs_f], "gene_b": [p[1] for p in pairs_f],
        "hit_union": union, "hit_h1299": np.asarray(h1299), "hit_mdamb231": np.asarray(mda),
        "signed_dd_maxabs": [agg_f[p][0] for p in pairs_f],
        "welch_t_maxabs": [agg_f[p][1] for p in pairs_f],
        "hedges_g_maxabs": [agg_f[p][2] for p in pairs_f]})
    dff_out.to_csv(OUT / "b1c_flister_three_stats.csv", index=False)

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "engine": "external_holdout.maxabs_dd_scores mirror; per (pair, lineage) best "
                  "|dd|/|t|/|g| independently maximized across orientations; frame score = "
                  "max across lineages; >=3/>=3 strata; sign convention dd = mut - wt "
                  "(|.| used throughout, as in the frozen Layer 1)",
        "fidelity": "per-pair max|DD| == frozen external_holdout_*_layer1.csv dd_min3 (tol 1e-8)",
        "seed": SEED, "n_perm": N_PERM,
        "harle": res_h, "flister": res_f,
    }
    (OUT / "b1c_external_three_stats_summary.json").write_text(
        json.dumps(summary, indent=2, default=str))
    print(f"  Harle: dd={res_h['abs_dd']['auroc']:.4f} t={res_h['abs_t']['auroc']:.4f} "
          f"g={res_h['abs_g']['auroc']:.4f}")
    print(f"  Flister: dd={res_f['abs_dd']['auroc']:.4f} t={res_f['abs_t']['auroc']:.4f} "
          f"g={res_f['abs_g']['auroc']:.4f}")

    # ── B6 lineage-matched Layer 1 ────────────────────────────────
    print("\n[B6] lineage-matched Layer 1")
    print(f"  DepMap lineages available: {len(set(lin_map.dropna()))}; "
          f"Harle-matched: {sorted(HARLE_MATCHED)}; Flister-matched: {sorted(FLISTER_MATCHED)}")
    for req, nm in ((HARLE_MATCHED, "Harle"), (FLISTER_MATCHED, "Flister")):
        missing = req - set(lin_map.dropna())
        if missing:
            raise RuntimeError(f"{nm} matched lineage(s) absent from DepMap lin_map: {missing}")

    agg_hm = aggregate(best_h, pairs_h, HARLE_MATCHED)
    agg_fm = aggregate(best_f, pairs_f, FLISTER_MATCHED)

    def dd_eval(labels, agg, pairs):
        s = np.array([abs(agg[p][0]) if agg[p][0] is not None else np.nan for p in pairs])
        y = np.asarray(labels, dtype=int)
        ok = ~np.isnan(s)
        r = auroc_with_permutation(y[ok], s[ok])
        r["n_pairs_scored"] = int(ok.sum())
        r["n_hits"] = int(y[ok].sum())
        return r

    b6 = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "design": "max|DD| restricted to DepMap lineages matching the screened models; "
                  "same >=3/>=3 strata and max-across-orientations selection; pairs with "
                  "no evaluable matched-lineage stratum are dropped (n reported)",
        "harle": {
            "matched_lineages": sorted(HARLE_MATCHED),
            "lineage_matched": dd_eval(hit_h, agg_hm, pairs_h),
            "all_lineage_reference": {"auroc": res_h["abs_dd"]["auroc"],
                                      "n_pairs_scored": res_h["n_scored"]},
        },
        "flister": {
            "matched_lineages": sorted(FLISTER_MATCHED),
            "lineage_matched_union": dd_eval(union, agg_fm, pairs_f),
            "all_lineage_reference": {"auroc": res_f["abs_dd"]["auroc"],
                                      "n_pairs_scored": res_f["n_scored"]},
            "per_model_matched": {
                "h1299_lung_only": dd_eval(
                    np.asarray(h1299),
                    aggregate(best_f, pairs_f, {"Non-Small Cell Lung Cancer"}), pairs_f),
                "mdamb231_breast_only": dd_eval(
                    np.asarray(mda),
                    aggregate(best_f, pairs_f, {"Invasive Breast Carcinoma"}), pairs_f),
            },
        },
    }
    (OUT / "b6_lineage_matched_layer1.json").write_text(json.dumps(b6, indent=2))

    b6_rows = pd.DataFrame({
        "pair": ["|".join(sorted(p)) for p in pairs_h],
        "harle_hit": hit_h,
        "abs_dd_all_lineages": [abs(agg_h[p][0]) if agg_h[p][0] is not None else np.nan
                                for p in pairs_h],
        "abs_dd_matched_lineages": [abs(agg_hm[p][0]) if agg_hm[p][0] is not None else np.nan
                                    for p in pairs_h],
    })
    b6_rows.to_csv(OUT / "b6_lineage_matched_scores.csv", index=False)

    hm_, fu_ = b6["harle"]["lineage_matched"], b6["flister"]["lineage_matched_union"]
    print(f"  Harle matched: AUROC={hm_['auroc']:.4f} (p={hm_['permutation_p']:.4f}, "
          f"n={hm_['n_pairs_scored']}) vs all-lineage {res_h['abs_dd']['auroc']:.4f}")
    print(f"  Flister matched union: AUROC={fu_['auroc']:.4f} (p={fu_['permutation_p']:.4f}, "
          f"n={fu_['n_pairs_scored']}) vs all-lineage {res_f['abs_dd']['auroc']:.4f}")
    print(f"\n  wrote b1c_* + b6_* to {OUT}")


if __name__ == "__main__":
    main()
