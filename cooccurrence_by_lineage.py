"""
Task 6 — Paralog mutation co-occurrence by lineage, with TMB adjustment
========================================================================
Reviewer request: break the pan-cancer co-occurrence analysis
(cooccurrence_analysis.py; 2x2 Fisher OR for driver/paralog mutation
co-occurrence across 1,208 DepMap lines) down by lineage and adjust for
total mutation burden (TMB).

Data flow (identical to cooccurrence_analysis.py):
  - universe: cell lines with dependency data (n=1,208)
  - mutation status: build_mutation_matrix(apply_driver_rules=True)
    (TSG: LikelyLoF; oncogene: Hotspot) from data/OmicsSomaticMutations.csv
  - lineage: Model.csv OncotreeLineage
  - TMB covariate: per-line count of nonsilent mutation records in the same
    filtered mutation file (default-entry profiles only), log1p-transformed

Per pair x lineage: Fisher exact OR + Wald 95% CI (0.5 Haldane correction
when a cell is zero), plus a logistic regression
    paralog_mut ~ driver_mut + log1p(TMB)
whose driver coefficient gives the TMB-adjusted OR. Lineages are evaluated
when they have >= 30 lines and >= 3 mutant lines for each gene of the pair.
Pan-cancer rows are included for reference.

Output: output/cooccurrence_by_lineage.json
"""

import json
import warnings

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import statsmodels.api as sm

from config import OUTPUT_DIR
from data_loader import (
    load_mutations, build_mutation_matrix, load_models, load_dependency,
)

PAIRS = [
    ("ARID1A", "ARID1B"),
    ("PIK3CA", "PIK3CB"),
    ("BRCA1", "BRCA2"),
    ("EP300", "CREBBP"),
    ("SMARCA4", "SMARCA2"),
]

MIN_LINES = 30
MIN_MUT_PER_GENE = 3


def fisher_pack(a_mut, b_mut):
    """2x2 Fisher OR + Wald CI (Haldane 0.5 if a zero cell)."""
    both = int(((a_mut == 1) & (b_mut == 1)).sum())
    a_only = int(((a_mut == 1) & (b_mut == 0)).sum())
    b_only = int(((a_mut == 0) & (b_mut == 1)).sum())
    neither = int(((a_mut == 0) & (b_mut == 0)).sum())
    or_val, p_val = fisher_exact([[both, a_only], [b_only, neither]],
                                 alternative="two-sided")
    bb, aa, cc, dd = both, a_only, b_only, neither
    if min(bb, aa, cc, dd) == 0:
        bb, aa, cc, dd = bb + 0.5, aa + 0.5, cc + 0.5, dd + 0.5
    log_or = np.log((bb * dd) / (aa * cc))
    se = np.sqrt(1 / bb + 1 / aa + 1 / cc + 1 / dd)
    return {
        "or": float(or_val),
        "ci": [float(np.exp(log_or - 1.96 * se)), float(np.exp(log_or + 1.96 * se))],
        "p": float(p_val),
        "n_both": both, "n_driver_only": a_only,
        "n_paralog_only": b_only, "n_neither": neither,
    }


def tmb_adjusted(driver_mut, paralog_mut, log_tmb):
    """Logistic regression paralog ~ driver + log_tmb -> adjusted OR for driver."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            X = sm.add_constant(np.column_stack([driver_mut, log_tmb]))
            model = sm.Logit(paralog_mut, X)
            fit = model.fit(disp=0, maxiter=100)
        if not fit.mle_retvals.get("converged", False):
            return None, "not_converged"
        coef = fit.params[1]
        se = fit.bse[1]
        if abs(coef) > 10 or se > 10:
            return None, "separation"
        return {
            "tmb_adjusted_or": float(np.exp(coef)),
            "ci_adj": [float(np.exp(coef - 1.96 * se)), float(np.exp(coef + 1.96 * se))],
            "p_adj": float(fit.pvalues[1]),
            "tmb_coef_p": float(fit.pvalues[2]),
        }, None
    except Exception as e:
        return None, f"fit_error: {type(e).__name__}"


def main():
    print("=" * 70)
    print("  Co-occurrence by lineage with TMB adjustment (DepMap 26Q1)")
    print("=" * 70)

    dep = load_dependency()
    cell_lines = list(dep.index)
    models = load_models().set_index("DepMap_ID")
    mut = load_mutations()

    genes = sorted({g for pair in PAIRS for g in pair})
    mat = build_mutation_matrix(mut, cell_lines, genes, apply_driver_rules=True)

    # TMB: count of nonsilent mutation records per line (same filtered file)
    tmb = mut.groupby("DepMap_ID").size().reindex(cell_lines).fillna(0)
    log_tmb = np.log1p(tmb)

    df = pd.DataFrame({"lineage": models["OncotreeLineage"].reindex(cell_lines)})
    df = df.join(mat).join(log_tmb.rename("log_tmb"))
    df = df.dropna(subset=["lineage"])
    print(f"  Lines with lineage annotation: {len(df)} "
          f"({df['lineage'].nunique()} lineages)")

    records = []
    skipped = []
    scopes = ["PAN-CANCER"] + sorted(df["lineage"].unique())
    for a, b in PAIRS:
        for scope in scopes:
            sub = df if scope == "PAN-CANCER" else df[df["lineage"] == scope]
            n = len(sub)
            n_a = int(sub[a].sum())
            n_b = int(sub[b].sum())
            if n < MIN_LINES or n_a < MIN_MUT_PER_GENE or n_b < MIN_MUT_PER_GENE:
                if scope != "PAN-CANCER":
                    skipped.append({"pair": f"{a}/{b}", "lineage": scope,
                                    "n": n, "n_driver_mut": n_a, "n_paralog_mut": n_b})
                continue
            fp = fisher_pack(sub[a].values, sub[b].values)
            adj, adj_err = tmb_adjusted(sub[a].values.astype(float),
                                        sub[b].values.astype(float),
                                        sub["log_tmb"].values.astype(float))
            rec = {
                "pair": f"{a}/{b}", "driver": a, "paralog": b,
                "lineage": scope, "n_lines": n,
                "n_driver_mut": n_a, "n_paralog_mut": n_b,
                "or": fp["or"], "ci": fp["ci"], "p": fp["p"],
                "n_both": fp["n_both"], "n_driver_only": fp["n_driver_only"],
                "n_paralog_only": fp["n_paralog_only"], "n_neither": fp["n_neither"],
                "tmb_adjusted_or": adj["tmb_adjusted_or"] if adj else None,
                "ci_adj": adj["ci_adj"] if adj else None,
                "p_adj": adj["p_adj"] if adj else None,
                "tmb_coef_p": adj["tmb_coef_p"] if adj else None,
                "adj_note": adj_err,
            }
            records.append(rec)
            if scope != "PAN-CANCER" or True:
                astr = (f"adjOR={adj['tmb_adjusted_or']:.2f} p_adj={adj['p_adj']:.2e}"
                        if adj else f"adj=NA({adj_err})")
                print(f"  {a}/{b:9s} {scope:22s} n={n:4d} OR={fp['or']:6.2f} "
                      f"p={fp['p']:.1e}  {astr}")

    out_obj = {
        "method": {
            "universe": "1,208 DepMap 26Q1 lines with dependency data",
            "mutation_rules": "TSG: LikelyLoF; oncogene: Hotspot "
                              "(build_mutation_matrix, apply_driver_rules=True)",
            "tmb": "per-line nonsilent mutation record count, log1p-transformed",
            "adjusted_model": "Logit(paralog_mut ~ driver_mut + log1p(TMB)) per "
                              "pair x lineage; adjusted OR = exp(driver coefficient)",
            "evaluability": f">= {MIN_LINES} lines and >= {MIN_MUT_PER_GENE} mutant "
                            "lines per gene within lineage",
        },
        "results": records,
        "skipped_lineages": skipped,
    }
    out = OUTPUT_DIR / "cooccurrence_by_lineage.json"
    with open(out, "w") as fh:
        json.dump(out_obj, fh, indent=2)
    print(f"\nSaved: {out} ({len(records)} pair x lineage tests, "
          f"{len(skipped)} skipped)")


if __name__ == "__main__":
    main()
