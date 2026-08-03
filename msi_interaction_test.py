"""
MSI Interaction Test — Higher-Powered Alternative to Subgroup AUROCs
====================================================================
Motivation: the MSI-H/MSS subgroup AUROC analysis (msi_analysis.py)
recomputes DD inside tiny subgroups (14–45 lines; 3–4 gold-standard
positives per subgroup), so both the subgroup AUROCs and their contrast
are severely underpowered. This script answers the same question —
does a hypermutated (MSI-H) background modulate driver-conditioned
paralog dependency? — using ALL lines of each lineage in a single
per-pair regression:

    gene_effect_P ~ mut_D * MSI_H            (binary, primary)
    gene_effect_P ~ mut_D * log1p(MSIscore)  (continuous, sensitivity)

The interaction coefficient is the difference in DD between MSI-H and
MSS backgrounds (ΔDD). Driver mutation calls use the pipeline's
gene-class rules (build_mutation_matrix, apply_driver_rules=True:
TSG = LikelyLoF, oncogene = Hotspot). A pair is evaluable when the
lineage cohort has >=3 mutant lines with >=2 mutants in EACH MSI class
(interaction estimable), matching the >=3 sensitivity frame.

Outputs
-------
output/msi_interaction_results.csv   per-pair coefficients and p-values,
                                     including the interaction coefficient
                                     (interaction_beta) with its Wald 95% CI
                                     (ci_lo, ci_hi; same SE convention as the
                                     model's p-values)
output/msi_interaction_summary.json  aggregate statistics for the
                                     manuscript + number audit
"""

import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from config import OUTPUT_DIR, KNOWN_PARALOG_SL
from data_loader import (
    load_dependency, load_models, load_mutations, load_paralogs,
    build_mutation_matrix,
)
from msi_analysis import MSI_CANCERS, MSI_DRIVERS, load_official_msi

MIN_MUT_TOTAL = 3   # matches the >=3 sensitivity frame
MIN_MUT_PER_MSI_CLASS = 2   # interaction estimability


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values (monotone step-up)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def fit_pair(df: pd.DataFrame):
    """Fit binary and continuous interaction models for one pair."""
    Xb = sm.add_constant(pd.DataFrame({
        "mut": df["mut"],
        "msi_h": df["msi_h"],
        "mut_x_msi": df["mut"] * df["msi_h"],
    }))
    fb = sm.OLS(df["ge"], Xb).fit()
    Xc = sm.add_constant(pd.DataFrame({
        "mut": df["mut"],
        "lscore": df["lscore"],
        "mut_x_lscore": df["mut"] * df["lscore"],
    }))
    fc = sm.OLS(df["ge"], Xc).fit()
    b_mut = fb.params["mut"]
    b_int = fb.params["mut_x_msi"]
    ci_int = fb.conf_int().loc["mut_x_msi"]  # Wald 95% CI, same (non-robust)
                                             # SE convention as the model's p-values
    return {
        # DD sign convention: positive DD = stronger dependency in mutant
        # (DD = mean ge|WT − mean ge|MUT = −β_mut).
        "dd_mss": -b_mut,
        "dd_msi_h": -(b_mut + b_int),
        "delta_dd": -b_int,          # >0 = stronger compensation in MSI-H
        "interaction_beta": b_int,   # raw OLS coefficient on mut x MSI-H
        "ci_lo": ci_int[0],
        "ci_hi": ci_int[1],
        "interaction_p": fb.pvalues["mut_x_msi"],
        "interaction_p_continuous": fc.pvalues["mut_x_lscore"],
    }


def run_interaction_test():
    print("=" * 65)
    print("  MSI Interaction Test (full-cohort per-pair regression)")
    print("=" * 65)

    dep = load_dependency()
    models = load_models()
    mutations = load_mutations()
    paralogs = load_paralogs()
    sig_df = load_official_msi()

    known_set = set()
    for a, b in KNOWN_PARALOG_SL:
        known_set.add((a.upper(), b.upper()))
        known_set.add((b.upper(), a.upper()))

    msi_by_id = sig_df.set_index("DepMap_ID")
    rows = []

    for cancer_name, disease_patterns in MSI_CANCERS.items():
        pat = "|".join(disease_patterns)
        mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
        cancer_models = models[mask]
        cell_ids = [c for c in cancer_models["DepMap_ID"]
                    if c in dep.index and c in msi_by_id.index]
        if len(cell_ids) < 10:
            print(f"  {cancer_name}: insufficient lines ({len(cell_ids)}) — skipped")
            continue

        msi_status = (msi_by_id.loc[cell_ids, "msi_official"] == "MSI_H").astype(float)
        msi_score = msi_by_id.loc[cell_ids, "msi_score"].astype(float)
        n_h = int((msi_status == 1).sum())
        print(f"\n  {cancer_name}: {len(cell_ids)} lines "
              f"({n_h} MSI-H, {len(cell_ids) - n_h} MSS)")

        n_driver_pairs = 0
        for driver in MSI_DRIVERS.get(cancer_name, []):
            driver_paralogs = paralogs[paralogs["gene_A"] == driver]["gene_B"].unique()
            valid_paralogs = [p for p in driver_paralogs if p in dep.columns]
            if not valid_paralogs:
                continue
            mut_matrix = build_mutation_matrix(mutations, cell_ids, [driver])
            if driver not in mut_matrix.columns:
                continue
            mut_flag = mut_matrix[driver].reindex(cell_ids).fillna(0).astype(float)

            for paralog in valid_paralogs:
                df = pd.DataFrame({
                    "ge": dep.loc[cell_ids, paralog],
                    "mut": mut_flag,
                    "msi_h": msi_status,
                    "lscore": np.log1p(msi_score),
                }).dropna()
                n_mut = int(df["mut"].sum())
                n_mut_h = int(((df["mut"] == 1) & (df["msi_h"] == 1)).sum())
                n_mut_s = n_mut - n_mut_h
                if (n_mut < MIN_MUT_TOTAL or n_mut_h < MIN_MUT_PER_MSI_CLASS
                        or n_mut_s < MIN_MUT_PER_MSI_CLASS):
                    continue
                try:
                    fit = fit_pair(df)
                except Exception:
                    continue
                rows.append({
                    "cancer": cancer_name, "driver": driver, "paralog": paralog,
                    "is_known_paralog_sl": (driver.upper(), paralog.upper()) in known_set,
                    "n_lines": len(df), "n_mut": n_mut,
                    "n_mut_msi_h": n_mut_h, "n_mut_mss": n_mut_s,
                    **fit,
                })
                n_driver_pairs += 1
        print(f"    evaluable pairs: {n_driver_pairs}")

    idf = pd.DataFrame(rows)
    if idf.empty:
        raise SystemExit("No evaluable pairs — check data inputs")

    # BH within cancer on the binary interaction p
    idf["interaction_q"] = np.nan
    for cancer_name, idx in idf.groupby("cancer").groups.items():
        idf.loc[idx, "interaction_q"] = bh_adjust(idf.loc[idx, "interaction_p"].values)
    idf.to_csv(OUTPUT_DIR / "msi_interaction_results.csv", index=False)

    # ── Aggregate per cancer ──
    summary = {}
    for cancer_name, sub in idf.groupby("cancer"):
        p = sub["interaction_p"].values
        n = len(sub)
        n_nom = int((p < 0.05).sum())
        binom_p = float(stats.binomtest(n_nom, n, 0.05, alternative="greater").pvalue)
        # Stouffer combined test, signed by ΔDD direction
        z = stats.norm.isf(p / 2) * np.sign(sub["delta_dd"].values)
        z_comb = float(z.sum() / np.sqrt(n))
        stouffer_p = float(2 * stats.norm.sf(abs(z_comb)))
        n_q = int((sub["interaction_q"] < 0.05).sum())
        spans_zero = (sub["ci_lo"] <= 0) & (sub["ci_hi"] >= 0)
        summary[cancer_name] = {
            "n_lines": int(sub["n_lines"].max()),
            "n_pairs": n,
            "n_nominal_p05": n_nom,
            "binom_p_vs_5pct": binom_p,
            "stouffer_z": z_comb,
            "stouffer_p": stouffer_p,
            "n_fdr05": n_q,
            "min_q": float(sub["interaction_q"].min()),
            "median_abs_delta_dd": float(sub["delta_dd"].abs().median()),
            "median_abs_beta": float(sub["interaction_beta"].abs().median()),
            "ci_span_zero": bool(spans_zero.all()),
            "n_ci_excludes_zero": int((~spans_zero).sum()),
            "known_pairs_evaluable": int(sub["is_known_paralog_sl"].sum()),
        }
        print(f"\n  {cancer_name}: {n} pairs, {n_nom} nominal p<0.05 "
              f"(binomial p={binom_p:.3f}), Stouffer p={stouffer_p:.3f}, "
              f"FDR<0.05: {n_q}, min q={summary[cancer_name]['min_q']:.3f}")

    # Lead-pair rows for manuscript traceability (if evaluable)
    lead = {}
    for drv, par in [("ARID1A", "ARID1B"), ("EP300", "CREBBP"),
                     ("PIK3CA", "PIK3CB"), ("SMARCA4", "SMARCA2"),
                     ("PTEN", "TNS1"), ("KRAS", "NRAS")]:
        hit = idf[(idf["driver"] == drv) & (idf["paralog"] == par)]
        for _, r in hit.iterrows():
            lead[f"{r['cancer']}:{drv}->{par}"] = {
                "dd_mss": round(float(r["dd_mss"]), 4),
                "dd_msi_h": round(float(r["dd_msi_h"]), 4),
                "delta_dd": round(float(r["delta_dd"]), 4),
                "interaction_p": float(r["interaction_p"]),
                "interaction_q": float(r["interaction_q"]),
                "n_mut": int(r["n_mut"]),
                "n_mut_msi_h": int(r["n_mut_msi_h"]),
                "n_mut_mss": int(r["n_mut_mss"]),
            }

    out = {
        "model": "OLS per pair: gene_effect_paralog ~ mut_driver * MSI_H "
                 "(binary, primary) and ~ mut_driver * log1p(MSIsensor2 score) "
                 "(continuous sensitivity); driver calls use gene-class rules",
        "estimability": f">={MIN_MUT_TOTAL} mutant lines with "
                        f">={MIN_MUT_PER_MSI_CLASS} mutants per MSI class",
        "per_cancer": summary,
        "lead_pairs": lead,
    }
    with open(OUTPUT_DIR / "msi_interaction_summary.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nSaved {OUTPUT_DIR}/msi_interaction_results.csv and msi_interaction_summary.json")


if __name__ == "__main__":
    run_interaction_test()
