#!/usr/bin/env python3
"""
rev_b1_ranking_stat_sensitivity.py  (Stage-4 revision, item B1 a+b)
=================================================================
Ranking-statistic sensitivity for the manuscript's signed Delta Dependency
(DD) ranking metric. Recomputes three ranking statistics on IDENTICAL
mutant/WT arms and compares their AUROC:

  1. signed DD      = mean(Chronos|WT) - mean(Chronos|MUT)   (manuscript Eq.1)
  2. Welch t        = signed DD / sqrt(s2_mut/n_mut + s2_wt/n_wt)
                      (sign-aligned with DD: positive = stronger dependency
                      in driver-mutant lines; this is minus the t statistic
                      pcs.py passes to scipy.ttest_ind)
  3. Hedges' g      = Cohen's d * J(N), J = 1 - 3/(4N-9)     (already shipped
                      in TableS2 / solid_*_results.csv; recomputed here as a
                      fidelity check and used directly)

Frames:
  (a) PRIMARY lineage-level frame: output/tables/TableS2_FullResults.tsv
      (110 driver x paralog x lineage entries, 8 positives, gyn3 lineages)
      -> overall AUROC x3, Spearman rank correlations, delta AUROC.
  (b) Fig. 1a per-lineage frame: output/solid_<lineage>_results.csv
      (23 solid lineages; 8 evaluable in the primary >=5/>=5 frame)
      -> per-lineage AUROC x3.
  Also writes the Methods-appendix n_mut/n_wt distribution table for the
  primary frame (per driver x lineage stratum).

Engine fidelity: every recomputed DD / Hedges' g / n_mut / n_wt must match
the frozen artifact value within 1e-9 (relative); the run aborts otherwise,
guaranteeing the Welch t statistics are computed on byte-identical arms.

No simulated data. Seed not needed (no resampling in this script except the
paired bootstrap of delta-AUROC, seed 42).

Outputs (output/revision_stage4/):
  b1_primary_frame_three_stats.csv     per-entry dd / welch_t / hedges_g + label
  b1_primary_frame_summary.json        AUROC x3, Spearman, delta AUROC + paired bootstrap CI
  b1_lineage_three_stats.csv           per-lineage AUROC x3 (23 lineages)
  b1_n_mut_n_wt_distribution.csv       per driver x lineage stratum n_mut / n_wt (primary frame)

Usage: python rev_b1_ranking_stat_sensitivity.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

from compute_headline_metrics import auroc  # noqa: E402  (rank-based, manuscript-identical)
from data_loader import (  # noqa: E402
    load_dependency, load_models, load_mutations,
    build_mutation_matrix, filter_gynecological_cell_lines, classify_cancer_type,
)
from pancancer import SOLID_TUMORS, SOLID_DRIVERS  # noqa: E402
from config import DEPMAP_FILES  # noqa: E402

TOL = 1e-9
BOOT_SEED = 42
N_BOOT = 10_000

GYN3 = ["Ovarian", "Endometrial", "Cervical"]


def welch_t_dd(mut_vals: np.ndarray, wt_vals: np.ndarray) -> float:
    """Welch t sign-aligned with DD: (mean_wt - mean_mut)/sqrt(v_m/n_m+v_w/n_w)."""
    n_m, n_w = len(mut_vals), len(wt_vals)
    if n_m < 2 or n_w < 2:
        return float("nan")
    v_m = mut_vals.var(ddof=1)
    v_w = wt_vals.var(ddof=1)
    se = np.sqrt(v_m / n_m + v_w / n_w)
    if se == 0:
        return float("nan")
    return float((wt_vals.mean() - mut_vals.mean()) / se)


def compute_entry(dep_sub: pd.DataFrame, m: pd.Series, paralog: str):
    """Given dep submatrix (cell lines x [paralog]) and mutation flag Series m,
    replicate pcs.py arm construction and return (dd, welch_t, hedges_g, n_mut, n_wt,
    valid_dep_counts)."""
    mut_cl = m[m == 1].index
    wt_cl = m[m == 0].index
    d = dep_sub[paralog]
    dep_mut = d.loc[mut_cl].dropna().to_numpy(dtype=float)
    dep_wt = d.loc[wt_cl].dropna().to_numpy(dtype=float)
    n_mut, n_wt = len(mut_cl), len(wt_cl)
    if len(dep_mut) >= 5 and len(dep_wt) >= 5:
        dd = dep_wt.mean() - dep_mut.mean()
        # Cohen's d: pcs.py uses pandas .var() (ddof=1) on the Series
        pooled_std = np.sqrt((pd.Series(dep_mut).var() + pd.Series(dep_wt).var()) / 2)
        cohens_d = dd / pooled_std if pooled_std > 0 else 0.0
        n_tot = len(dep_mut) + len(dep_wt)
        g = cohens_d * (1 - 3 / (4 * n_tot - 9)) if n_tot > 3 else cohens_d
        t = welch_t_dd(dep_mut, dep_wt)
        return float(dd), t, float(g), n_mut, n_wt, True
    # pcs.py fallback: dd=0, d=0, g=0, p=1 (entry still exported)
    return 0.0, float("nan"), 0.0, n_mut, n_wt, False


def check_close(name, got, ref, tol=TOL):
    if not (np.isnan(got) and (ref is None or (isinstance(ref, float) and np.isnan(ref)))):
        if ref is None or abs(got - ref) > tol * max(1.0, abs(ref)):
            raise RuntimeError(f"FIDELITY FAIL {name}: got {got!r} vs artifact {ref!r}")


def main():
    t0 = datetime.now()
    print("=" * 72)
    print("  rev B1: ranking-statistic sensitivity (primary + Fig.1a lineages)")
    print("=" * 72)

    tables2 = pd.read_csv(ROOT / "output" / "tables" / "TableS2_FullResults.tsv", sep="\t")
    tables2["is_known_paralog_sl"] = tables2["is_known_paralog_sl"].astype(bool)

    # ── load raw data ─────────────────────────────────────────────
    print("  loading DepMap raw data ...")
    dep = load_dependency()
    models = load_models()
    mut = load_mutations()
    expr_ids = set(pd.read_csv(DEPMAP_FILES["expression"], usecols=[0]).iloc[:, 0])
    print(f"  dep {dep.shape}; models {len(models)}; mut records {len(mut)}")

    # ════════════ (a) PRIMARY gyn3 frame ════════════
    print("\n[a] primary gyn3 frame (TableS2, 110 entries)")
    cell_line_types = classify_cancer_type(models)
    gyn_models = filter_gynecological_cell_lines(models, None)
    gyn_ids_all = [c for c in gyn_models["DepMap_ID"].tolist()
                   if c in dep.index and c in expr_ids]

    rows = []
    n_fallback = 0
    for ctype in GYN3:
        cl_subset = [c for c in gyn_ids_all if cell_line_types.get(c, "Other") == ctype]
        sub = tables2[tables2["cancer_type"] == ctype]
        for driver in sorted(sub["driver_gene"].unique()):
            mat = build_mutation_matrix(mut, cl_subset, [driver])
            if driver not in mat.columns:
                raise RuntimeError(f"driver {driver} absent from mutation matrix")
            m = mat[driver]
            entries = sub[sub["driver_gene"] == driver]
            for r in entries.itertuples():
                dd, t, g, n_mut, n_wt, valid = compute_entry(dep.loc[cl_subset], m, r.paralog_gene)
                # fidelity vs frozen artifact
                check_close(f"dd {driver}->{r.paralog_gene} {ctype}", dd, float(r.dependency_dd))
                check_close(f"g {driver}->{r.paralog_gene} {ctype}", g, float(r.hedges_g))
                if n_mut != int(r.n_mut) or n_wt != int(r.n_wt):
                    raise RuntimeError(f"FIDELITY FAIL n {driver} {ctype}: {n_mut}/{n_wt} vs {r.n_mut}/{r.n_wt}")
                if not valid:
                    n_fallback += 1
                rows.append({
                    "driver_gene": driver, "paralog_gene": r.paralog_gene,
                    "cancer_type": ctype, "signed_dd": dd, "welch_t": t,
                    "hedges_g": g, "n_mut": n_mut, "n_wt": n_wt,
                    "sufficient_dep_n": valid,
                    "is_known_paralog_sl": bool(r.is_known_paralog_sl),
                })
    prim = pd.DataFrame(rows)
    prim.to_csv(OUT / "b1_primary_frame_three_stats.csv", index=False)
    print(f"    {len(prim)} entries recomputed; fidelity OK; "
          f"{n_fallback} entries below per-paralog >=5/>=5 dep counts (pcs fallback dd=0)")

    yt = prim["is_known_paralog_sl"].astype(int).to_numpy()
    scores = {
        "signed_dd": prim["signed_dd"].fillna(0).to_numpy(),
        "welch_t": prim["welch_t"].fillna(0).to_numpy(),
        "hedges_g": prim["hedges_g"].fillna(0).to_numpy(),
    }
    aucs = {k: auroc(yt, v) for k, v in scores.items()}
    sp_dd_t = spearmanr(scores["signed_dd"], scores["welch_t"]).statistic
    sp_dd_g = spearmanr(scores["signed_dd"], scores["hedges_g"]).statistic
    sp_t_g = spearmanr(scores["welch_t"], scores["hedges_g"]).statistic

    # paired bootstrap for delta AUROC (t - dd, g - dd), entry resampling
    rng = np.random.default_rng(BOOT_SEED)
    n = len(yt)
    d_t, d_g, b_dd, b_t, b_g = [], [], [], [], []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        yb = yt[idx]
        if yb.sum() == 0 or yb.sum() == n:
            continue
        a_dd = auroc(yb, scores["signed_dd"][idx])
        a_t = auroc(yb, scores["welch_t"][idx])
        a_g = auroc(yb, scores["hedges_g"][idx])
        b_dd.append(a_dd); b_t.append(a_t); b_g.append(a_g)
        d_t.append(a_t - a_dd); d_g.append(a_g - a_dd)
    d_t, d_g = np.array(d_t), np.array(d_g)

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "frame": "PRIMARY lineage-level (TableS2): 110 entries, 8 positives, gyn3",
        "definitions": {
            "signed_dd": "mean(Chronos|WT) - mean(Chronos|MUT) (manuscript Eq.1)",
            "welch_t": "signed_dd / sqrt(s2_mut/n_mut + s2_wt/n_wt), sign-aligned with DD",
            "hedges_g": "Cohen's d * J(N), J = 1-3/(4N-9) (pcs.py definition)",
            "scoring": "AUROC with NaN->0 fill (mirrors compute_headline_metrics.lineage_metrics)",
        },
        "n_entries": int(len(prim)),
        "n_positives": int(yt.sum()),
        "n_entries_below_dep_count_fallback": int(n_fallback),
        "auroc": {k: float(v) for k, v in aucs.items()},
        "delta_auroc": {
            "welch_t_minus_signed_dd": float(aucs["welch_t"] - aucs["signed_dd"]),
            "hedges_g_minus_signed_dd": float(aucs["hedges_g"] - aucs["signed_dd"]),
        },
        "spearman": {
            "signed_dd_vs_welch_t": float(sp_dd_t),
            "signed_dd_vs_hedges_g": float(sp_dd_g),
            "welch_t_vs_hedges_g": float(sp_t_g),
        },
        "paired_bootstrap": {
            "n_boot": N_BOOT, "seed": BOOT_SEED, "scheme": "entry resampling (same draws for all three statistics)",
            "welch_t_minus_signed_dd": {
                "mean": float(d_t.mean()),
                "ci95": [float(np.percentile(d_t, 2.5)), float(np.percentile(d_t, 97.5))],
                "frac_below_0": float((d_t < 0).mean())},
            "hedges_g_minus_signed_dd": {
                "mean": float(d_g.mean()),
                "ci95": [float(np.percentile(d_g, 2.5)), float(np.percentile(d_g, 97.5))],
                "frac_below_0": float((d_g < 0).mean())},
            "auroc_ci95": {
                "signed_dd": [float(np.percentile(b_dd, 2.5)), float(np.percentile(b_dd, 97.5))],
                "welch_t": [float(np.percentile(b_t, 2.5)), float(np.percentile(b_t, 97.5))],
                "hedges_g": [float(np.percentile(b_g, 2.5)), float(np.percentile(b_g, 97.5))],
            },
        },
    }

    # ── Methods-appendix n_mut/n_wt distribution (per driver x lineage stratum) ──
    strata = (prim.groupby(["cancer_type", "driver_gene"], as_index=False)
              .agg(n_mut=("n_mut", "first"), n_wt=("n_wt", "first"),
                   n_paralogs=("paralog_gene", "count")))
    strata = strata.sort_values(["cancer_type", "driver_gene"])
    strata.to_csv(OUT / "b1_n_mut_n_wt_distribution.csv", index=False)
    dist = {}
    for ctype in GYN3:
        s = strata[strata["cancer_type"] == ctype]
        dist[ctype] = {
            "n_strata": int(len(s)),
            "n_mut_min": int(s["n_mut"].min()), "n_mut_median": float(s["n_mut"].median()),
            "n_mut_max": int(s["n_mut"].max()),
            "n_wt_min": int(s["n_wt"].min()), "n_wt_median": float(s["n_wt"].median()),
            "n_wt_max": int(s["n_wt"].max()),
        }
    summary["n_mut_n_wt_distribution_by_lineage"] = dist
    summary["n_mut_n_wt_distribution_table"] = "b1_n_mut_n_wt_distribution.csv"

    # ════════════ (b) Fig. 1a per-lineage frame (23 solid lineages) ════════════
    print("\n[b] per-lineage frames (solid_*_results.csv, 23 lineages)")
    lineage_cells = {}
    for cancer_name, patterns in SOLID_TUMORS.items():
        pat = "|".join(patterns)
        mask = models["OncotreePrimaryDisease"].str.contains(pat, case=False, na=False)
        lin_mod = models[mask]
        if len(lin_mod) < 6:
            continue
        cell_ids = [c for c in lin_mod["DepMap_ID"].tolist()
                    if c in dep.index and c in expr_ids]
        if len(cell_ids) >= 6:
            lineage_cells[cancer_name] = cell_ids

    solid_summary = pd.read_csv(ROOT / "output" / "solid_tumor_summary.csv")
    lin_rows = []
    for cancer_name, cell_ids in lineage_cells.items():
        safe = cancer_name.replace("/", "_").replace(" ", "_")
        f = ROOT / "output" / f"solid_{safe}_results.csv"
        if not f.exists():
            continue
        ref = pd.read_csv(f)
        drivers = sorted(ref["driver_gene"].unique())
        stats_rows = []
        for driver in drivers:
            if driver not in dep.columns:
                continue
            mat = build_mutation_matrix(mut, cell_ids, [driver])
            if driver not in mat.columns:
                continue
            m = mat[driver]
            if m.sum() < 5 or (len(m) - m.sum()) < 5:
                # pancancer.py: compute_pcs_for_driver returns empty -> no entries
                continue
            for r in ref[ref["driver_gene"] == driver].itertuples():
                dd, t, g, n_mut, n_wt, valid = compute_entry(dep.loc[cell_ids], m, r.paralog_gene)
                check_close(f"dd {driver}->{r.paralog_gene} {cancer_name}", dd, float(r.dependency_dd))
                check_close(f"g {driver}->{r.paralog_gene} {cancer_name}", g, float(r.hedges_g))
                stats_rows.append({"signed_dd": dd, "welch_t": t, "hedges_g": g,
                                   "known": bool(r.is_known_paralog_sl)})
        ldf = pd.DataFrame(stats_rows)
        row = {"cancer": cancer_name, "n_lines": len(cell_ids),
               "n_pairs": len(ldf)}
        if ldf.empty:
            row.update({"n_known": 0, "auroc_signed_dd": float("nan"),
                        "auroc_welch_t": float("nan"), "auroc_hedges_g": float("nan"),
                        "delta_t_minus_dd": float("nan"), "delta_g_minus_dd": float("nan")})
            lin_rows.append(row)
            continue
        yl = ldf["known"].astype(int).to_numpy()
        nk = int(yl.sum())
        row["n_known"] = nk
        if nk >= 2:
            row["auroc_signed_dd"] = auroc(yl, ldf["signed_dd"].fillna(0).to_numpy())
            row["auroc_welch_t"] = auroc(yl, ldf["welch_t"].fillna(0).to_numpy())
            row["auroc_hedges_g"] = auroc(yl, ldf["hedges_g"].fillna(0).to_numpy())
            # fidelity vs solid_tumor_summary.csv
            ref_auc = solid_summary.loc[solid_summary["cancer"] == cancer_name, "dd_auroc"]
            if len(ref_auc) and not np.isnan(ref_auc.iloc[0]):
                check_close(f"lineage auroc {cancer_name}", row["auroc_signed_dd"],
                            float(ref_auc.iloc[0]), tol=1e-9)
            row["delta_t_minus_dd"] = row["auroc_welch_t"] - row["auroc_signed_dd"]
            row["delta_g_minus_dd"] = row["auroc_hedges_g"] - row["auroc_signed_dd"]
        else:
            row["auroc_signed_dd"] = row["auroc_welch_t"] = row["auroc_hedges_g"] = float("nan")
            row["delta_t_minus_dd"] = row["delta_g_minus_dd"] = float("nan")
        lin_rows.append(row)
        if nk >= 2:
            print(f"    {cancer_name:20s} nk={nk}  dd={row['auroc_signed_dd']:.4f}  "
                  f"t={row['auroc_welch_t']:.4f}  g={row['auroc_hedges_g']:.4f}")
    lins = pd.DataFrame(lin_rows)
    lins.to_csv(OUT / "b1_lineage_three_stats.csv", index=False)

    ev = lins.dropna(subset=["auroc_signed_dd"])
    summary["lineage_level"] = {
        "n_lineages_total": int(len(lins)),
        "n_lineages_evaluable_primary_ge5": int(len(ev)),
        "evaluable_lineages": ev["cancer"].tolist(),
        "mean_delta_t_minus_dd": float(ev["delta_t_minus_dd"].mean()),
        "mean_delta_g_minus_dd": float(ev["delta_g_minus_dd"].mean()),
        "median_delta_t_minus_dd": float(ev["delta_t_minus_dd"].median()),
        "median_delta_g_minus_dd": float(ev["delta_g_minus_dd"].median()),
        "n_lineages_t_beats_dd": int((ev["delta_t_minus_dd"] > 0).sum()),
        "n_lineages_g_beats_dd": int((ev["delta_g_minus_dd"] > 0).sum()),
        "table": "b1_lineage_three_stats.csv",
    }

    (OUT / "b1_primary_frame_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False, default=str))
    print(f"\n  primary frame: dd={aucs['signed_dd']:.4f} t={aucs['welch_t']:.4f} g={aucs['hedges_g']:.4f}")
    print(f"  spearman dd-t={sp_dd_t:.4f} dd-g={sp_dd_g:.4f} t-g={sp_t_g:.4f}")
    print(f"  wrote {OUT}/b1_primary_frame_summary.json + 3 CSVs "
          f"({(datetime.now() - t0).total_seconds():.0f}s)")


if __name__ == "__main__":
    main()
