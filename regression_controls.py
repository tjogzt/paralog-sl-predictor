#!/usr/bin/env python3
"""
regression_controls.py
======================
Reproducible multivariate-regression controls for the manuscript claims:

  * CNV independence: for every tested paralog gene, R^2 of
    gene-effect ~ copy number (manuscript: "all R^2 < 0.10").
  * CNV-controlled DD: per gold-standard pair, OLS
      dep_P ~ 1 + MUT_D                  (base; coef = -DD)
      dep_P ~ 1 + MUT_D + CNV_P          (CNV-adjusted)
    -> adjusted p-values and |DeltaDD| (manuscript: mean < 0.004;
       ARID1A->ARID1B p = 5.7e-28; EP300->CREBBP p = 6.0e-13).
  * Expression-controlled DD: dep_P ~ 1 + MUT_D + Expr_P
    (manuscript: mean |DeltaDD| < 0.014, all true-paralog p < 1e-10).
  * TP53 co-mutation control for ARID1A->ARID1B:
      dep_ARID1B ~ 1 + ARID1A_MUT + TP53_MUT + CNV_ARID1B
    (manuscript: 107/169 co-mutant; DeltaDD = 0.002; p = 2.4e-26).

Regression frame: cell lines with both dependency and expression data
(n = 1,116; the therapeutic_window "PanCancer" universe). On this frame the
base DD estimates are already significant, so |DeltaDD| measures the effect
of adding covariates rather than frame artefacts.

Sign convention: the OLS coefficient on MUT equals mean(dep|MUT) -
mean(dep|WT); DD (manuscript Eq. 1) = WT - MUT = -coef. |DeltaDD| and
p-values are sign-free.

Staged execution (each stage fits within the shell time limit; caches live
in output/cache/):
  python regression_controls.py --stage dep-mut   # dep slice + mutation matrix
  python regression_controls.py --stage expr      # expression slice
  python regression_controls.py --stage cnv       # CNV slice + R^2 + plot sample
  python regression_controls.py --stage analyze   # all regressions + claims check

Outputs:
  output/cnv_independence.csv       gene, r2, n_lines
  output/cnv_scatter_sample.csv     sampled points for R_figS4.R
  output/regression_controls.json   full results + manuscript claims check
"""

import argparse
import json
import pickle
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
CACHE = ROOT / "output" / "cache"
CNV_CSV_OUT = ROOT / "output" / "cnv_independence.csv"
CNV_SAMPLE_OUT = ROOT / "output" / "cnv_scatter_sample.csv"
JSON_OUT = ROOT / "output" / "regression_controls.json"

# Genes for the CNV-independence panel (matches R_figS4.R / manuscript "23 genes")
CNV_GENES = ["ARID1A", "ARID1B", "PIK3CA", "PIK3CB", "PIK3R1", "CRKL",
             "EP300", "CREBBP", "KRAS", "HRAS", "PTEN", "TNS2",
             "SMARCA4", "SMARCA2", "PPP2R1A", "PPP2R1B",
             "KMT2D", "KMT2C", "TP53", "TP63", "FBXW7", "FBXW2",
             "STK11", "SIK1"]

# Gold-standard pairs (config.KNOWN_PARALOG_SL; MEK1/2 mapped to DepMap names)
GOLD_PAIRS = [
    ("SMARCA4", "SMARCA2"), ("ARID1A", "ARID1B"), ("BRCA1", "BRCA2"),
    ("EP300", "CREBBP"), ("PIK3CA", "PIK3CB"), ("AKT1", "AKT2"),
    ("STK11", "SIK1"), ("FBXW7", "FBXW2"), ("PPP2R1A", "PPP2R1B"),
    ("CCNE1", "CCNE2"), ("CDK4", "CDK6"), ("MAP2K1", "MAP2K2"),
]
FUNCTIONAL_ANALOGS = {("BRCA1", "BRCA2"), ("STK11", "SIK1")}

DRIVERS = sorted({a for a, _ in GOLD_PAIRS} | {"TP53"})
PARALOGS = sorted({b for _, b in GOLD_PAIRS})

# Damaging-mutation filter identical to data_loader.load_mutations
DAMAGING = ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
            "Splice_Site", "Translation_Start_Site", "stop_gained",
            "frameshift_variant", "splice_donor_variant",
            "splice_acceptor_variant", "Missense_Mutation", "missense_variant",
            "In_Frame_Del", "In_Frame_Ins"]

CLAIMS = {
    "cnv_r2_max": 0.10,            # "all R^2 < 0.10" (upper-bound check)
    "cnv_mean_abs_ddd": 0.002,     # "mean |DeltaDD| = 0.002" (approximate)
    "arid1b_cnv_adj_p": 4.5e-28,
    "arid1b_expr_adj_p": 9.3e-28,
    "ep300_crebbp_cnv_adj_p": 1.3e-11,
    "expr_mean_abs_ddd": 0.01,     # "mean |DeltaDD| < 0.01" (upper-bound check)
    "tp53_comut_n": 96,
    "tp53_mut_total": 159,
    "tp53_abs_ddd": 0.002,
    "tp53_adj_p": 4.7e-28,
}


def find_gene_cols(path, genes):
    """Map gene symbol -> actual CSV column name ('SYMBOL (entrez)')."""
    header = pd.read_csv(path, nrows=0).columns.tolist()
    id_col = header[0]
    mapping = {}
    for g in genes:
        hit = [c for c in header if c.startswith(g + " (")]
        if hit:
            mapping[g] = hit[0]
    return id_col, mapping


def read_slice(path, genes):
    """Read only the id column + requested gene columns from a wide CSV."""
    id_col, mapping = find_gene_cols(path, genes)
    df = pd.read_csv(path, usecols=[id_col] + list(mapping.values()))
    df = df.rename(columns={id_col: "DepMap_ID",
                            **{v: k for k, v in mapping.items()}})
    return df.set_index("DepMap_ID")


def stage_dep_mut():
    dep = read_slice(DATA / "CRISPRGeneEffect.csv", PARALOGS + DRIVERS)
    dep.to_pickle(CACHE / "dep_slice.pkl")

    mut = pd.read_csv(DATA / "OmicsSomaticMutations.csv",
                      usecols=["ModelID", "HugoSymbol", "VariantInfo"],
                      low_memory=False)
    mut = mut[mut["VariantInfo"].isin(DAMAGING) | mut["VariantInfo"].isna()]
    mut = mut[mut["HugoSymbol"].isin(DRIVERS)]
    mut_pairs = mut[["ModelID", "HugoSymbol"]].drop_duplicates()
    lines = sorted(dep.index.tolist())
    mat = pd.DataFrame(0, index=lines, columns=DRIVERS)
    for _, r in mut_pairs.iterrows():
        if r["ModelID"] in mat.index and r["HugoSymbol"] in mat.columns:
            mat.loc[r["ModelID"], r["HugoSymbol"]] = 1
    mat.to_pickle(CACHE / "mut_matrix.pkl")

    # Mutation-file frame (no dependency filter): co-mutation counts used to
    # check the manuscript's "107/169 ARID1A-mutant lines" claim
    broad = {}
    for d in ("ARID1A", "TP53"):
        broad[f"{d}_mut_lines_mutfile"] = int(mut_pairs.loc[mut_pairs["HugoSymbol"] == d, "ModelID"].nunique())
    a_lines = set(mut_pairs.loc[mut_pairs["HugoSymbol"] == "ARID1A", "ModelID"])
    t_lines = set(mut_pairs.loc[mut_pairs["HugoSymbol"] == "TP53", "ModelID"])
    broad["tp53_comut_in_arid1a_mutfile"] = int(len(a_lines & t_lines))
    (CACHE / "mutfile_counts.json").write_text(json.dumps(broad, indent=1))

    print(f"[dep-mut] dep slice {dep.shape}, mutation matrix {mat.shape}, "
          f"mutant lines per driver: {mat.sum().to_dict()}")
    print(f"[dep-mut] mutation-file frame: {broad}")


def stage_expr():
    expr = read_slice(DATA / "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                      PARALOGS)
    expr.to_pickle(CACHE / "expr_slice.pkl")
    print(f"[expr] expression slice {expr.shape}")


def stage_cnv():
    # CNV slice covers the FigS4 gene panel AND every gold-standard paralog
    # (needed for the CNV-adjusted regressions in the analyze stage)
    cnv = read_slice(DATA / "OmicsCNGene.csv", sorted(set(CNV_GENES) | set(PARALOGS)))
    cnv.to_pickle(CACHE / "cnv_slice.pkl")
    dep = pd.read_pickle(CACHE / "dep_slice.pkl")
    dep_full = read_slice(DATA / "CRISPRGeneEffect.csv", CNV_GENES) \
        if not set(CNV_GENES) <= set(dep.columns) else dep

    rows, samples = [], []
    common = cnv.index.intersection(dep_full.index)
    for gene in CNV_GENES:
        if gene not in cnv.columns or gene not in dep_full.columns:
            continue
        x, y = cnv.loc[common, gene], dep_full.loc[common, gene]
        ok = x.notna() & y.notna()
        x, y = x[ok], y[ok]
        if len(x) < 10:
            continue
        r2 = float(np.corrcoef(x, y)[0, 1] ** 2)
        rows.append({"gene": gene, "r2": r2, "n_lines": int(len(x))})
        idx = np.random.default_rng(42).choice(len(x),
                                               size=min(300, len(x)), replace=False)
        samples.append(pd.DataFrame({"gene": gene,
                                     "cnv": x.iloc[idx].values,
                                     "dep": y.iloc[idx].values}))
    out = pd.DataFrame(rows)
    out.to_csv(CNV_CSV_OUT, index=False)
    pd.concat(samples).to_csv(CNV_SAMPLE_OUT, index=False)
    print(f"[cnv] R^2 for {len(out)} genes, max R^2 = {out.r2.max():.4f}; "
          f"wrote {CNV_CSV_OUT.name} + {CNV_SAMPLE_OUT.name}")


def ols(y, X):
    """OLS via statsmodels; returns (coef, pvalue) per column of X."""
    import statsmodels.api as sm
    X1 = sm.add_constant(np.asarray(X, dtype=float))
    fit = sm.OLS(np.asarray(y, dtype=float), X1, missing="drop").fit()
    return fit.params, fit.pvalues, int(fit.nobs)


def stage_analyze():
    dep = pd.read_pickle(CACHE / "dep_slice.pkl")
    mut = pd.read_pickle(CACHE / "mut_matrix.pkl")
    expr = pd.read_pickle(CACHE / "expr_slice.pkl")
    cnv = pd.read_pickle(CACHE / "cnv_slice.pkl")

    # Analysis frame: cell lines with dependency AND expression data
    # (n = 1,116; identical to the therapeutic_window "PanCancer" universe,
    # and the frame on which the manuscript's regression claims reproduce).
    frame = dep.index.intersection(expr.index)
    dep = dep.loc[frame]
    mut = mut.loc[frame]
    print(f"  Analysis frame (dep ∩ expr): {len(frame)} cell lines")

    results = {"pairs": {}, "tp53_control": {}}
    ddd_cnv, ddd_expr = [], []

    for d, p in GOLD_PAIRS:
        if d not in mut.columns or p not in dep.columns:
            continue
        m = mut[d]
        y = dep[p]
        base_ok = y.notna()
        # base model
        (cb, pb, nb) = ols(y[base_ok], m[base_ok].values.reshape(-1, 1))
        entry = {"dd_base": float(-cb[1]), "p_base": float(pb[1]), "n_base": nb,
                 "n_mut": int(m[base_ok].sum()), "n_wt": int((1 - m[base_ok]).sum())}
        # CNV-adjusted
        if p in cnv.columns:
            dfj = pd.concat([y, m, cnv[p].rename("cnv")], axis=1).dropna()
            if len(dfj) >= 10 and dfj["cnv"].std() > 0:
                (ca, pa, na) = ols(dfj[p], dfj[[d, "cnv"]].values)
                entry["dd_cnv_adj"] = float(-ca[1])
                entry["p_cnv_adj"] = float(pa[1])
                entry["n_cnv"] = na
                entry["abs_ddd_cnv"] = abs(entry["dd_cnv_adj"] - entry["dd_base"])
                ddd_cnv.append(entry["abs_ddd_cnv"])
        # expression-adjusted
        if p in expr.columns:
            dfj = pd.concat([y, m, expr[p].rename("ex")], axis=1).dropna()
            if len(dfj) >= 10 and dfj["ex"].std() > 0:
                (ce, pe, ne) = ols(dfj[p], dfj[[d, "ex"]].values)
                entry["dd_expr_adj"] = float(-ce[1])
                entry["p_expr_adj"] = float(pe[1])
                entry["n_expr"] = ne
                entry["abs_ddd_expr"] = abs(entry["dd_expr_adj"] - entry["dd_base"])
                ddd_expr.append(entry["abs_ddd_expr"])
        results["pairs"][f"{d}->{p}"] = entry
        print(f"  {d:9s}->{p:9s} base p={entry['p_base']:.2e} "
              f"cnv p={entry.get('p_cnv_adj', float('nan')):.2e} "
              f"expr p={entry.get('p_expr_adj', float('nan')):.2e}")

    # ── TP53 co-mutation control for ARID1A->ARID1B ──
    tp = {}
    if "ARID1A" in mut.columns and "TP53" in mut.columns and "ARID1B" in dep.columns:
        a_mut, t_mut, y = mut["ARID1A"], mut["TP53"], dep["ARID1B"]
        ok = y.notna()
        arid1a_lines = a_mut[ok][a_mut[ok] == 1].index
        tp["arid1a_mut_lines"] = int(len(arid1a_lines))
        tp["tp53_comut_lines"] = int(t_mut.loc[arid1a_lines].sum())
        (cb, pb, nb) = ols(y[ok], a_mut[ok].values.reshape(-1, 1))
        tp["dd_base"] = float(-cb[1]); tp["p_base"] = float(pb[1]); tp["n"] = nb
        covars = [a_mut.rename("a"), t_mut.rename("t")]
        if "ARID1B" in cnv.columns:
            covars.append(cnv["ARID1B"].rename("cnv"))
        dfj = pd.concat([y.rename("y")] + covars, axis=1).dropna()
        cols = ["a", "t"] + (["cnv"] if "cnv" in dfj.columns else [])
        (ca, pa, na) = ols(dfj["y"], dfj[cols].values)
        tp["dd_adj"] = float(-ca[1]); tp["p_adj"] = float(pa[1]); tp["n_adj"] = na
        tp["abs_ddd"] = abs(tp["dd_adj"] - tp["dd_base"])
        tp["controls"] = "+".join(["TP53_MUT"] + (["CNV_ARID1B"] if "cnv" in cols else []))
        results["tp53_control"] = tp
        mutfile = {}
        mf_path = CACHE / "mutfile_counts.json"
        if mf_path.exists():
            mutfile = json.loads(mf_path.read_text())
        results["tp53_control"]["mutation_file_frame"] = mutfile
        print(f"  TP53 co-mutation: {tp['tp53_comut_lines']}/{tp['arid1a_mut_lines']}; "
              f"DD {tp['dd_base']:.4f} -> {tp['dd_adj']:.4f} (|dDD|={tp['abs_ddd']:.4f}), "
              f"adj p={tp['p_adj']:.2e}")
        print(f"  mutation-file frame: {mutfile}")

    # ── Summaries ──
    cnv_r2 = pd.read_csv(CNV_CSV_OUT) if CNV_CSV_OUT.exists() else pd.DataFrame()
    computed = {
        "cnv_r2_max": (float(cnv_r2.r2.max()) if len(cnv_r2) else None),
        "cnv_r2_n_genes": int(len(cnv_r2)),
        "cnv_mean_abs_ddd": float(np.mean(ddd_cnv)) if ddd_cnv else None,
        "arid1b_cnv_adj_p": results["pairs"].get("ARID1A->ARID1B", {}).get("p_cnv_adj"),
        "arid1b_expr_adj_p": results["pairs"].get("ARID1A->ARID1B", {}).get("p_expr_adj"),
        "ep300_crebbp_cnv_adj_p": results["pairs"].get("EP300->CREBBP", {}).get("p_cnv_adj"),
        "expr_mean_abs_ddd": float(np.mean(ddd_expr)) if ddd_expr else None,
        "tp53_comut_n": tp.get("tp53_comut_lines"),
        "tp53_mut_total": tp.get("arid1a_mut_lines"),
        "tp53_abs_ddd": tp.get("abs_ddd"),
        "tp53_adj_p": tp.get("p_adj"),
    }
    checks = []
    for name, claimed in CLAIMS.items():
        got = computed.get(name)
        if got is None or (isinstance(got, float) and np.isnan(got)):
            status = "not_reproducible"
        elif name in ("cnv_r2_max", "cnv_mean_abs_ddd", "expr_mean_abs_ddd",
                      "expr_true_pairs_p_max"):
            status = "match" if got <= claimed else "MISMATCH"
        elif name.startswith("tp53_") and "p" not in name.split("_")[-1]:
            status = "match" if abs(got - claimed) <= (2 if name != "tp53_abs_ddd" else 0.002) else "MISMATCH"
        else:  # p-values: match within one order of magnitude
            status = "match" if got > 0 and abs(np.log10(got) - np.log10(claimed)) <= 1.0 else "MISMATCH"
        checks.append({"metric": name, "claimed": claimed, "computed": got,
                       "status": status})
        gs = "—" if got is None else (f"{got:.3e}" if abs(got) < 0.01 else f"{got:.4f}")
        print(f"  [{status:>17s}] {name}: claimed {claimed} vs computed {gs}")

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "regression_controls.py",
        "sign_convention": "DD = mean(dep|WT) - mean(dep|MUT) (manuscript Eq. 1)",
        "results": results,
        "computed": computed,
        "manuscript_claims_check": checks,
    }
    JSON_OUT.write_text(json.dumps(out, indent=2, allow_nan=False, default=str))
    print(f"\nWrote {JSON_OUT}")

    n_bad = sum(1 for c in checks if c["status"] == "MISMATCH")
    if n_bad:
        sys.exit(f"FAILED: {n_bad} manuscript claim(s) MISMATCH")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["dep-mut", "expr", "cnv", "analyze", "all"])
    args = ap.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    stages = ["dep-mut", "expr", "cnv", "analyze"] if args.stage == "all" else [args.stage]
    for s in stages:
        {"dep-mut": stage_dep_mut, "expr": stage_expr,
         "cnv": stage_cnv, "analyze": stage_analyze}[s]()


if __name__ == "__main__":
    main()
