#!/usr/bin/env python3
"""
rev2_b11_msi_cluster_null.py  (Stage-4 revision, item B11)
===========================================================
Problem: the shipped MSI interaction aggregation (msi_interaction_test.py,
output/msi_interaction_summary.json) combines per-pair interaction p-values
with (i) a one-sided binomial test of the nominal p<0.05 count against the
5% null rate and (ii) a signed Stouffer combination. Both assume the 47
(colorectal) / 24 (endometrial) per-pair p-values are independent. They are
not: pairs share driver genes (identical mutation vector reused across
paralogs), share paralog genes, and are fitted on the SAME cell lines, so
cross-pair p-value correlation is structural.

Fix implemented here: a joint-label permutation null that preserves the
cross-pair dependence structure. In every permutation the MSI-H/MSS labels
of the lineage's cell lines are permuted ONCE and applied jointly to all
pairs of that lineage (mutation vectors and gene-effect vectors untouched).
Every per-pair OLS is refit and the aggregate statistics (nominal count,
signed Stouffer z) recomputed, yielding a null distribution under the true
dependence structure. Empirical p-values replace the binomial test; the
Stouffer z is calibrated by its permutation SD (dependence-corrected
Stouffer). BH per-pair q-values (valid under positive regression
dependence) are retained unchanged from the shipped artifact.

Procedure:
  1. Rebuild every evaluable pair's design exactly as msi_interaction_test.py
     (same estimability rules: >=3 mutants, >=2 per MSI class; driver calls
     via build_mutation_matrix gene-class rules; official MSIsensor2 labels).
  2. Fidelity: refit all pairs with statsmodels OLS and require per-pair
     interaction p-values to match output/msi_interaction_results.csv
     (max abs diff reported; aborts above 1e-8).
  3. B=10,000 joint-label permutations (seed 42), batched numpy OLS.
  4. Report: observed n_nom / binomial p (as shipped), permutation
     empirical p for the count statistic, observed Stouffer z/p,
     permutation SD of z, dependence-corrected Stouffer p, and the
     min BH q (unchanged, from the shipped artifact).

Output: output/revision_stage4/rev2_b11_msi_cluster_null.json
Usage:  python rev2_b11_msi_cluster_null.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

from config import OUTPUT_DIR, KNOWN_PARALOG_SL  # noqa: E402
from data_loader import (  # noqa: E402
    load_dependency, load_models, load_mutations, load_paralogs,
    build_mutation_matrix,
)
from msi_analysis import MSI_CANCERS, MSI_DRIVERS, load_official_msi  # noqa: E402

MIN_MUT_TOTAL = 3
MIN_MUT_PER_MSI_CLASS = 2
B = 10_000
SEED = 42
FID_TOL = 1e-8


def rebuild_pairs():
    """Reconstruct every evaluable pair's row set exactly as
    msi_interaction_test.py, keeping raw vectors for refitting."""
    dep = load_dependency()
    models = load_models()
    mutations = load_mutations()
    paralogs = load_paralogs()
    sig_df = load_official_msi()
    msi_by_id = sig_df.set_index("DepMap_ID")

    pairs = {c: [] for c in MSI_CANCERS}
    for cancer_name, disease_patterns in MSI_CANCERS.items():
        pat = "|".join(disease_patterns)
        mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
        cancer_models = models[mask]
        cell_ids = [c for c in cancer_models["DepMap_ID"]
                    if c in dep.index and c in msi_by_id.index]
        msi_status = (msi_by_id.loc[cell_ids, "msi_official"] == "MSI_H").astype(float).to_numpy()
        id_pos = {c: i for i, c in enumerate(cell_ids)}

        for driver in MSI_DRIVERS.get(cancer_name, []):
            driver_paralogs = paralogs[paralogs["gene_A"] == driver]["gene_B"].unique()
            valid_paralogs = [p for p in driver_paralogs if p in dep.columns]
            if not valid_paralogs:
                continue
            mut_matrix = build_mutation_matrix(mutations, cell_ids, [driver])
            if driver not in mut_matrix.columns:
                continue
            mut_flag = mut_matrix[driver].reindex(cell_ids).fillna(0).astype(float).to_numpy()

            for paralog in valid_paralogs:
                ge = dep.loc[cell_ids, paralog].to_numpy(dtype=float)
                ok = ~np.isnan(ge)
                idx = np.where(ok)[0]
                m = mut_flag[idx]
                s = msi_status[idx]
                n_mut = int(m.sum())
                n_mut_h = int(((m == 1) & (s == 1)).sum())
                n_mut_s = n_mut - n_mut_h
                if (n_mut < MIN_MUT_TOTAL or n_mut_h < MIN_MUT_PER_MSI_CLASS
                        or n_mut_s < MIN_MUT_PER_MSI_CLASS):
                    continue
                pairs[cancer_name].append({
                    "driver": driver, "paralog": paralog,
                    "rows": idx, "mut": m, "ge": ge[idx],
                    "is_known": (driver.upper(), paralog.upper())
                                in {(a.upper(), b.upper()) for a, b in KNOWN_PARALOG_SL},
                })
    return pairs, MSI_CANCERS


def ols_fit_statsmodels(mut, msi, ge):
    X = sm.add_constant(pd.DataFrame({
        "mut": mut, "msi_h": msi, "mut_x_msi": mut * msi,
    }))
    fb = sm.OLS(ge, X).fit()
    return float(fb.pvalues["mut_x_msi"]), float(fb.params["mut_x_msi"])


def batch_ols_pvalues(S, m, y):
    """Vectorized OLS of y ~ 1 + m + s + m*s for many label vectors s.

    S: (B, n) binary label draws; m: (n,) mutation; y: (n,) gene effect.
    Returns interaction p-value and coefficient per draw (B,)."""
    n = len(y)
    X = np.empty((S.shape[0], n, 4))
    X[:, :, 0] = 1.0
    X[:, :, 1] = m[None, :]
    X[:, :, 2] = S
    X[:, :, 3] = S * m[None, :]
    XtX = np.einsum("bni,bnj->bij", X, X)
    Xty = np.einsum("bni,n->bi", X, y)
    rank = np.linalg.matrix_rank(XtX)
    degenerate = rank < 4
    XtX_inv = np.linalg.pinv(XtX)
    beta = np.einsum("bij,bj->bi", XtX_inv, Xty)
    resid = y[None, :] - np.einsum("bni,bi->bn", X, beta)
    rss = (resid ** 2).sum(axis=1)
    df = n - 4
    sigma2 = rss / df
    se = np.sqrt(np.maximum(sigma2 * XtX_inv[:, 3, 3], 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = beta[:, 3] / se
    p = 2.0 * stats.t.sf(np.abs(tstat), df)
    # interaction not estimable under collinear label draws: mark NaN
    p[degenerate] = np.nan
    return p, beta[:, 3]


def main():
    t0 = datetime.now()
    print("=" * 72)
    print("  rev2 B11: MSI interaction aggregation — joint-label permutation null")
    print("=" * 72)

    shipped = pd.read_csv(OUTPUT_DIR / "msi_interaction_results.csv")
    with open(OUTPUT_DIR / "msi_interaction_summary.json") as fh:
        shipped_summary = json.load(fh)

    pairs, _ = rebuild_pairs()

    # ── observed fits + fidelity check ──────────────────────────────
    # Lineage-level official MSI labels, attached per cancer.
    max_diff = 0.0
    dep = load_dependency()
    models = load_models()
    sig_df = load_official_msi()
    msi_by_id = sig_df.set_index("DepMap_ID")
    labels = {}
    for cancer_name, disease_patterns in MSI_CANCERS.items():
        pat = "|".join(disease_patterns)
        mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
        cell_ids = [c for c in models.loc[mask, "DepMap_ID"]
                    if c in dep.index and c in msi_by_id.index]
        labels[cancer_name] = {
            "cell_ids": cell_ids,
            "msi": (msi_by_id.loc[cell_ids, "msi_official"] == "MSI_H").astype(float).to_numpy(),
        }

    results = {}
    rng = np.random.default_rng(SEED)
    for cancer_name, plist in pairs.items():
        msi_obs = labels[cancer_name]["msi"]
        obs_stats = {"n_pairs": len(plist), "n_lines": len(msi_obs),
                     "n_msi_h": int(msi_obs.sum())}
        p_obs, delta_obs = [], []
        sub = shipped[shipped["cancer"] == cancer_name]
        ship_p = {(r.driver, r.paralog): r.interaction_p for r in sub.itertuples()}
        for p in plist:
            s = msi_obs[p["rows"]]
            pv, bint = ols_fit_statsmodels(p["mut"], s, p["ge"])
            ref = ship_p[(p["driver"], p["paralog"])]
            d = abs(pv - ref)
            max_diff = max(max_diff, d)
            if d > FID_TOL:
                raise RuntimeError(
                    f"FIDELITY FAIL {cancer_name} {p['driver']}->{p['paralog']}: "
                    f"recomputed p={pv!r} vs shipped {ref!r}")
            p_obs.append(pv)
            delta_obs.append(-bint)  # delta_dd = -interaction_beta
        p_obs = np.asarray(p_obs)
        delta_obs = np.asarray(delta_obs)
        n_nom = int((p_obs < 0.05).sum())
        z_obs = float((stats.norm.isf(p_obs / 2) * np.sign(delta_obs)).sum()
                      / np.sqrt(len(p_obs)))
        stouffer_p_obs = float(2 * stats.norm.sf(abs(z_obs)))
        binom_p_obs = float(stats.binomtest(n_nom, len(p_obs), 0.05,
                                            alternative="greater").pvalue)

        # ── joint-label permutation null ────────────────────────────
        n_lines = len(msi_obs)
        P = np.empty((B, n_lines))
        for b in range(B):
            P[b] = msi_obs[rng.permutation(n_lines)]
        cnt_null = np.zeros(B, dtype=int)
        z_null = np.zeros(B)
        n_degenerate = 0
        for p in plist:
            S = P[:, p["rows"]]
            pv, bint = batch_ols_pvalues(S, p["mut"], p["ge"])
            n_degenerate += int(np.isnan(pv).sum())
            valid = ~np.isnan(pv)
            cnt_null += ((pv < 0.05) & valid)
            zc = np.where(valid, stats.norm.isf(np.clip(pv, 1e-300, 1) / 2)
                          * np.sign(-bint), 0.0)
            z_null += zc
        z_null /= np.sqrt(len(plist))

        p_count_emp = float((1 + np.sum(cnt_null >= n_nom)) / (1 + B))
        p_z_emp = float((1 + np.sum(np.abs(z_null) >= abs(z_obs))) / (1 + B))
        sd_z_null = float(z_null.std(ddof=1))
        z_corr = z_obs / sd_z_null if sd_z_null > 0 else np.nan
        p_stouffer_corrected = float(2 * stats.norm.sf(abs(z_corr)))

        ship_c = shipped_summary["per_cancer"][cancer_name]
        results[cancer_name] = {
            **obs_stats,
            "shipped_binom_p_vs_5pct": {"recomputed": binom_p_obs,
                                        "artifact": ship_c["binom_p_vs_5pct"]},
            "shipped_stouffer": {"z_recomputed": z_obs,
                                 "z_artifact": ship_c["stouffer_z"],
                                 "p_recomputed": stouffer_p_obs,
                                 "p_artifact": ship_c["stouffer_p"]},
            "n_nominal_p05": n_nom,
            "n_degenerate_draws_total": n_degenerate,
            "permutation_null": {
                "scheme": "MSI-H/MSS labels of the lineage's cell lines permuted "
                          "once per draw and applied jointly to all pairs; mutation "
                          "vectors and gene effects fixed; all per-pair OLS refit; "
                          f"B={B}, seed={SEED}",
                "count_stat": {"null_mean": float(cnt_null.mean()),
                               "null_max": int(cnt_null.max()),
                               "empirical_p_count_ge_observed": p_count_emp},
                "stouffer_z": {"null_sd": sd_z_null,
                               "independence_sd": 1.0,
                               "dependence_inflation_var_ratio": sd_z_null ** 2,
                               "empirical_p_two_sided": p_z_emp,
                               "dependence_corrected_p": p_stouffer_corrected},
            },
            "bh_unchanged": {"min_q": ship_c["min_q"], "n_fdr05": ship_c["n_fdr05"]},
        }
        print(f"  {cancer_name}: {len(plist)} pairs; n_nom={n_nom} "
              f"(binom p={binom_p_obs:.3f} -> perm p={p_count_emp:.4f}); "
              f"Stouffer z={z_obs:+.3f} (p={stouffer_p_obs:.3f}) "
              f"nullSD={sd_z_null:.3f} -> corrected p={p_stouffer_corrected:.4f}; "
              f"min BH q={ship_c['min_q']:.3f}")

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "rev2_b11_msi_cluster_null.py",
        "item": "B11 — replace independence-assuming binomial + Stouffer aggregation "
                "of MSI interaction p-values with a joint-label permutation null "
                "that preserves cross-pair dependence (shared drivers, paralogs, "
                "and cell lines)",
        "fidelity": {"max_abs_diff_recomputed_vs_shipped_interaction_p": max_diff,
                     "tolerance": FID_TOL},
        "per_cancer": results,
    }
    out_path = OUT / "rev2_b11_msi_cluster_null.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n  fidelity: max |p_recomputed - p_shipped| = {max_diff:.2e}")
    print(f"  written: {out_path}  ({(datetime.now()-t0).total_seconds():.0f}s)")


if __name__ == "__main__":
    main()
