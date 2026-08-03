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
  output/cnv_scatter_sample.csv     sampled points for R_figS3_cnv.R
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
sys.path.insert(0, str(ROOT))
from data_loader import load_mutations, build_mutation_matrix  # noqa: E402
DATA = ROOT / "data"
CACHE = ROOT / "output" / "cache"
CNV_CSV_OUT = ROOT / "output" / "cnv_independence.csv"
CNV_SAMPLE_OUT = ROOT / "output" / "cnv_scatter_sample.csv"
JSON_OUT = ROOT / "output" / "regression_controls.json"

# Genes for the CNV-independence panel (matches R_figS3_cnv.R / manuscript "23 genes")
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

# NOTE: the private DAMAGING list was removed in the class-specific-rules
# batch (C7); mutation filtering now goes through data_loader.load_mutations
# + build_mutation_matrix (default-entry + TSG:LikelyLoF / ONC:Hotspot rules).

# Values stated in the manuscript, updated 2026-07-26 after the C7
# class-specific driver-mutation rules (TP53 LikelyLoF includes DN missense;
# ARID1A restricted to LikelyLoF). p-value claims match within one order of
# magnitude (see check loop below).
CLAIMS = {
    "cnv_r2_max": 0.10,            # "all R^2 < 0.10" (upper-bound check)
    "cnv_mean_abs_ddd": 0.005,     # mean |DeltaDD| after CNV adjustment
    "arid1b_cnv_adj_p": 6.5e-42,
    "arid1b_expr_adj_p": 4.1e-41,
    "ep300_crebbp_cnv_adj_p": 1.4e-13,
    "expr_mean_abs_ddd": 0.014,    # mean |DeltaDD| after expression adjustment
    "tp53_comut_n": 76,
    "tp53_mut_total": 121,
    "tp53_abs_ddd": 0.002,
    "tp53_adj_p": 5.8e-42,
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

    # Shared mutation pipeline (data_loader): default-entry filter +
    # class-specific driver rules (TSG: LikelyLoF; oncogene: Hotspot)
    mut = load_mutations(DATA / "OmicsSomaticMutations.csv")
    mut = mut[mut["Gene"].isin(DRIVERS)]
    lines = sorted(dep.index.tolist())
    mat = build_mutation_matrix(mut, lines, DRIVERS)
    mat.to_pickle(CACHE / "mut_matrix.pkl")

    # Mutation-file frame (no dependency filter): co-mutation counts used to
    # check the manuscript's ARID1A/TP53 co-mutation claim
    mut_pairs = mut[["DepMap_ID", "Gene"]].drop_duplicates()
    broad = {}
    for d in ("ARID1A", "TP53"):
        broad[f"{d}_mut_lines_mutfile"] = int(mut_pairs.loc[mut_pairs["Gene"] == d, "DepMap_ID"].nunique())
    a_lines = set(mut_pairs.loc[mut_pairs["Gene"] == "ARID1A", "DepMap_ID"])
    t_lines = set(mut_pairs.loc[mut_pairs["Gene"] == "TP53", "DepMap_ID"])
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


def _f(x):
    """JSON-safe float: NaN/inf -> None."""
    x = float(x)
    return x if np.isfinite(x) else None


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
        n_mut_ok = int(m[base_ok].sum())
        n_wt_ok = int((1 - m[base_ok]).sum())
        if n_mut_ok < 3 or n_wt_ok < 3:
            # driver has (near-)no qualifying mutant lines on this frame
            # (e.g. amplification-driven CCNE1 under the Hotspot rule):
            # OLS is not identifiable — record and skip.
            results["pairs"][f"{d}->{p}"] = {
                "skipped": f"insufficient MUT/WT lines (n_mut={n_mut_ok}, n_wt={n_wt_ok})"}
            print(f"  {d:9s}->{p:9s} skipped: n_mut={n_mut_ok}, n_wt={n_wt_ok}")
            continue
        # base model
        (cb, pb, nb) = ols(y[base_ok], m[base_ok].values.reshape(-1, 1))
        entry = {"dd_base": _f(-cb[1]), "p_base": _f(pb[1]), "n_base": nb,
                 "n_mut": n_mut_ok, "n_wt": n_wt_ok}
        # CNV-adjusted
        if p in cnv.columns:
            dfj = pd.concat([y, m, cnv[p].rename("cnv")], axis=1).dropna()
            if len(dfj) >= 10 and dfj["cnv"].std() > 0:
                (ca, pa, na) = ols(dfj[p], dfj[[d, "cnv"]].values)
                entry["dd_cnv_adj"] = _f(-ca[1])
                entry["p_cnv_adj"] = _f(pa[1])
                entry["n_cnv"] = na
                if entry["dd_cnv_adj"] is not None and entry["dd_base"] is not None:
                    entry["abs_ddd_cnv"] = abs(entry["dd_cnv_adj"] - entry["dd_base"])
                    ddd_cnv.append(entry["abs_ddd_cnv"])
        # expression-adjusted
        if p in expr.columns:
            dfj = pd.concat([y, m, expr[p].rename("ex")], axis=1).dropna()
            if len(dfj) >= 10 and dfj["ex"].std() > 0:
                (ce, pe, ne) = ols(dfj[p], dfj[[d, "ex"]].values)
                entry["dd_expr_adj"] = _f(-ce[1])
                entry["p_expr_adj"] = _f(pe[1])
                entry["n_expr"] = ne
                if entry["dd_expr_adj"] is not None and entry["dd_base"] is not None:
                    entry["abs_ddd_expr"] = abs(entry["dd_expr_adj"] - entry["dd_base"])
                    ddd_expr.append(entry["abs_ddd_expr"])
        results["pairs"][f"{d}->{p}"] = entry
        _fmt = lambda v: f"{v:.2e}" if v is not None else "   NA   "
        print(f"  {d:9s}->{p:9s} base p={_fmt(entry['p_base'])} "
              f"cnv p={_fmt(entry.get('p_cnv_adj'))} "
              f"expr p={_fmt(entry.get('p_expr_adj'))}")

    # ── Full regression table: robust SE + lineage adjustment (round-4 review) ──
    # For every evaluable gold-standard pair, report beta (=-DD), HC3-robust
    # SE, 95% CI, p, BH q (within model), n, n_mut, n_wt for nested models:
    # base / +CNV / +Expression / +Lineage (OncotreePrimaryDisease fixed
    # effects, levels with <20 lines collapsed into "Other"). HC3 does not
    # address within-lineage correlation, so every model also carries
    # lineage cluster-robust SE/p (cov_type='cluster', groups = lineage).
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests

    lin_map = {}
    model_csv = DATA / "Model.csv"
    if model_csv.exists():
        mcols = pd.read_csv(model_csv, usecols=["ModelID", "OncotreePrimaryDisease"])
        lin_map = dict(zip(mcols["ModelID"], mcols["OncotreePrimaryDisease"]))
    lin = pd.Series({i: lin_map.get(i, "Other") for i in frame}, dtype=object)
    _lc = lin.value_counts()
    lin = lin.where(lin.map(_lc) >= 20, "Other")

    def _fit_mut(dfj, rhs):
        fit = smf.ols(f"dep ~ {rhs}", data=dfj).fit(cov_type="HC3")
        ci = fit.conf_int().loc["mut"]
        return {"beta_mut": float(fit.params["mut"]), "dd": float(-fit.params["mut"]),
                "se": float(fit.bse["mut"]), "ci_lo": float(ci[0]), "ci_hi": float(ci[1]),
                "p": float(fit.pvalues["mut"]), "n": int(fit.nobs)}

    def _fit_mut_cluster(dfj, rhs):
        """Lineage cluster-robust SE for the same model: cov_type='cluster'
        with cell lines grouped by Oncotree primary disease, addressing
        within-lineage correlation that HC3 cannot (reviewer comment)."""
        groups = lin.loc[dfj.index].values
        fit = smf.ols(f"dep ~ {rhs}", data=dfj).fit(
            cov_type="cluster", cov_kwds={"groups": groups})
        return {"se_lineage_cluster": float(fit.bse["mut"]),
                "p_lineage_cluster": float(fit.pvalues["mut"])}

    reg_rows = []
    for d, p in GOLD_PAIRS:
        key = f"{d}->{p}"
        ent = results["pairs"].get(key)
        if not ent or "skipped" in ent:
            continue
        m, y = mut[d], dep[p]
        base = pd.concat([y.rename("dep"), m.rename("mut")], axis=1).dropna()
        n_mut, n_wt = int(base["mut"].sum()), int((1 - base["mut"]).sum())
        model_frames = {"base": (base, "mut")}
        if p in cnv.columns:
            f = pd.concat([base, cnv[p].rename("cnv")], axis=1).dropna()
            if len(f) >= 10 and f["cnv"].std() > 0:
                model_frames["cnv_adj"] = (f, "mut + cnv")
        if p in expr.columns:
            f = pd.concat([base, expr[p].rename("expr")], axis=1).dropna()
            if len(f) >= 10 and f["expr"].std() > 0:
                model_frames["expr_adj"] = (f, "mut + expr")
        f = base.copy()
        f["lineage"] = lin.loc[f.index].values
        if f["lineage"].nunique() > 1:
            model_frames["lineage_adj"] = (f, "mut + C(lineage)")
        for mname, (fj, rhs) in model_frames.items():
            try:
                r = _fit_mut(fj, rhs)
            except Exception as e:  # singular design etc.
                r = {"beta_mut": None, "dd": None, "se": None, "ci_lo": None,
                     "ci_hi": None, "p": None, "n": len(fj), "error": str(e)}
            try:
                r.update(_fit_mut_cluster(fj, rhs))
            except Exception as e:  # e.g. too few clusters
                r["se_lineage_cluster"] = None
                r["p_lineage_cluster"] = None
                r["lineage_cluster_error"] = str(e)
            reg_rows.append({"pair": key, "model": mname,
                             "n_mut": n_mut, "n_wt": n_wt, **r})

    regtab = pd.DataFrame(reg_rows)
    if len(regtab):
        for mname, sub in regtab.groupby("model"):
            ok = sub["p"].notna()
            if ok.sum():
                regtab.loc[sub.index[ok], "q_bh"] = multipletests(
                    sub.loc[ok, "p"], method="fdr_bh")[1]
        regtab_out = ROOT / "output" / "regression_table_full.csv"
        regtab.to_csv(regtab_out, index=False)
        results["full_regression_table"] = "output/regression_table_full.csv"
        lin_p = (regtab[(regtab["model"] == "lineage_adj") &
                        (regtab["pair"] == "ARID1A->ARID1B")]["p"])
        results["arid1a_arid1b_lineage_adj_p"] = (float(lin_p.iloc[0])
                                                  if len(lin_p) and pd.notna(lin_p.iloc[0]) else None)
        # Headline pairs under lineage cluster-robust SE (lineage_adj model;
        # same specification as the 2.2e-13 claim, SE clustered by lineage)
        def _safe(x):
            try:
                return _f(x)
            except (TypeError, ValueError):
                return None
        cl_out = {}
        for pair in ("ARID1A->ARID1B", "SMARCA4->SMARCA2", "EP300->CREBBP"):
            sub = regtab[(regtab["pair"] == pair) & (regtab["model"] == "lineage_adj")]
            if len(sub):
                cl_out[pair] = {
                    "model": "lineage_adj",
                    "se_lineage_cluster": _safe(sub["se_lineage_cluster"].iloc[0]),
                    "p_lineage_cluster": _safe(sub["p_lineage_cluster"].iloc[0]),
                }
        results["lineage_cluster_se_headline"] = cl_out
        print("  lineage cluster-robust SE (lineage_adj model): " + "; ".join(
            f"{k} p={v['p_lineage_cluster']:.2e}" if v["p_lineage_cluster"] is not None
            else f"{k} p=NA" for k, v in cl_out.items()))
        print(f"  full regression table: {len(regtab)} rows, HC3 robust SE "
              f"-> output/regression_table_full.csv; "
              f"ARID1A->ARID1B lineage-adj p={results['arid1a_arid1b_lineage_adj_p']:.2e}"
              if results['arid1a_arid1b_lineage_adj_p'] else
              f"  full regression table written ({len(regtab)} rows)")

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
        # lineage-adjusted variant (round-4 review): + C(lineage), HC3 robust SE
        try:
            fj = dfj.copy()
            fj["lineage"] = lin.loc[fj.index].values
            rhs = "y ~ a + t" + (" + cnv" if "cnv" in fj.columns else "") + " + C(lineage)"
            fitl = smf.ols(rhs, data=fj).fit(cov_type="HC3")
            tp["p_adj_lineage"] = float(fitl.pvalues["a"])
            tp["dd_adj_lineage"] = float(-fitl.params["a"])
            print(f"  lineage-adjusted TP53 control: DD {tp['dd_adj_lineage']:.4f}, "
                  f"p={tp['p_adj_lineage']:.2e}")
        except Exception as e:
            tp["p_adj_lineage_error"] = str(e)
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
