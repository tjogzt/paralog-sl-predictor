#!/usr/bin/env python3
"""
external_holdout.py
===================
Locked external hold-out validation of the Delta Dependency (DD) framework
against two independent combinatorial CRISPR screens that became public
after our candidate set was frozen:

1. Harle et al. 2025 (Genome Biology, DOI 10.1186/s13059-025-03737-w)
   - 472 predicted SL pairs x 27 cell lines (melanoma/pancreas/lung).
   - Local file: data/external/harle2025_moesm1.xlsx (pre-cached).
   - Table S5 (header row 4): 0/1 hit-call matrix (pairs x 27 lines).
   - Table S4 (header row 5): long-format pair x line GI table with
     continuous normalized GI score (mean_norm_gi), is_bassik_hit call,
     fdr and depMapID. Their code defines
     is_bassik_hit = mean_norm_gi < -0.5 & fdr < 0.01 & neither single
     target depleted (07_combine_intermediate_datasets.R in the Harle
     GitHub repo mirrored at data/external/harle2025).

2. Flister et al. 2025 (Cell Reports 44:116512, DOI
   10.1016/j.celrep.2025.116512; bioRxiv 10.1101/2024.07.16.603642)
   - Digenic enAsCas12a screen of 36,648 paralog pairs in NCI-H1299 and
     MDA-MB-231 + meta-analysis of 462 pairs across 49 cancer models.
   - cell.com / ScienceDirect / PMC were all inaccessible from this
     machine (Cloudflare 403; paper not yet deposited in PMC as of
     2026-08-01, verified via NCBI idconv). The identical study's
     supplementary tables were therefore obtained from the bioRxiv
     preprint (v3, posted 2025-01-16; v1 2024-07-19), files cached as
     data/external/flister2025_biorxiv_*.{xlsx,pdf,txt}.
   - Table S3 Sheet1: per-dataset "Lethal" calls for 35,108 pairs
     (their two models + six external digenic studies).
   - Table S6 'SL_all_combat': 462 pairs x 49 cell models continuous
     dzL2FC (ComBat-normalized); hit = dzL2FC <= -2 per their methods
     ("dzL2FC > 2 ... considered significant"); verified ~98%
     concordant with their '# SL' column.

Analysis (two layers, exactly as pre-specified):
  Layer 1 (mutation-agnostic): for each screen, pairs overlapping our
    universe are scored by max |DD| (>=3-mutant/>=3-WT frame, max across
    BOTH orientations x all lineages — the exact scoring of
    in4mer_benchmark.py / in4mer_seed_sensitivity.py). AUROC vs the
    screen's experimental hits, 10,000-permutation p (seed 42).
  Layer 2 (genotype-stratified): for our 12 curated directed pairs
    (driver -> paralog) present in a screen, external cell lines are
    classified driver-mutant vs WT using the project's TableS6 mutation
    rules (TSG: LikelyLoF; ONC: Hotspot — via
    build_mutation_matrix(apply_driver_rules=True)). Hit rates compared
    by Fisher exact; continuous GI scores by Mann-Whitney (one-sided in
    the pre-registered direction, two-sided also reported). Per-pair
    BH-adjusted q-values and a pooled analysis are reported.

Engine fidelity: the vectorized |DD| scorer is validated against the
frozen output/in4mer_benchmark.csv (all 413 pairs, dd_min3 column);
the run aborts if any pair disagrees.

No simulated, random, or hard-coded data. Everything derives from the
cached raw files; downloads are skipped when caches are present and
valid (offline re-runnable).

Outputs (all in output/):
  external_holdout.json
  external_holdout_harle_layer1.csv
  external_holdout_harle_layer2.csv
  external_holdout_harle_layer2_lines.csv
  external_holdout_flister_layer1.csv
  external_holdout_flister_layer2.csv
  external_holdout_flister_layer2_models.csv

Usage:
  python external_holdout.py            # full end-to-end run
  python external_holdout.py <stage>    # one stage only:
      inputs | depmap | harle | flister_l1 | flister_l2
"""

import hashlib
import json
import sys
import time
import urllib.request
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, rankdata
from sklearn.metrics import roc_auc_score

from config import OUTPUT_DIR, DATA_DIR
from data_loader import (
    load_dependency,
    load_models,
    load_mutations,
    build_mutation_matrix,
)

# ─────────────────────────── constants ───────────────────────────
ROOT = Path(__file__).resolve().parent
EXT = DATA_DIR / "external"
CACHE = OUTPUT_DIR / "cache"
CACHE.mkdir(parents=True, exist_ok=True)
STATE_PKL = CACHE / "ext_holdout_state.pkl"
PARTIAL_DIR = CACHE / "ext_holdout_partials"
PARTIAL_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_PERM = 10_000
MIN_N = 3                     # >=3-mutant / >=3-WT frame (in4mer precedent)
LINEAGE_COL = "OncotreePrimaryDisease"

HARLE_XLSX = EXT / "harle2025_moesm1.xlsx"

# bioRxiv 10.1101/2024.07.16.603642 v3 (2025-01-16) supplementary files
BIORXIV_BASE = ("https://www.biorxiv.org/content/biorxiv/early/2025/01/16/"
                "2024.07.16.603642")
FLISTER_REMOTE = {
    f"flister2025_biorxiv_suppl_tableS{i}.xlsx":
        f"{BIORXIV_BASE}/DC{i}/embed/media-{i}.xlsx?download=true"
    for i in range(1, 11)
}
FLISTER_REMOTE["flister2025_biorxiv_suppl_figures.pdf"] = (
    f"{BIORXIV_BASE}/DC11/embed/media-11.pdf?download=true")
FLISTER_REMOTE["flister2025_biorxiv_fulltext.txt"] = (
    "https://www.biorxiv.org/content/10.1101/2024.07.16.603642v3.full.txt")
# Only the files this analysis actually parses need to be valid on re-run.
FLISTER_REQUIRED = [
    "flister2025_biorxiv_suppl_tableS3.xlsx",
    "flister2025_biorxiv_suppl_tableS6.xlsx",
]

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
      "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36")

# 12 curated directed pairs (driver -> paralog) with tiers
CURATED = [
    ("AKT1", "AKT2", "A"), ("CDK4", "CDK6", "A"),
    ("MAP2K1", "MAP2K2", "A"),
    ("SMARCA4", "SMARCA2", "B"), ("ARID1A", "ARID1B", "B"),
    ("EP300", "CREBBP", "C"), ("PIK3CA", "PIK3CB", "C"),
    ("CCNE1", "CCNE2", "C"), ("FBXW7", "FBXW2", "C"),
    ("PPP2R1A", "PPP2R1B", "C"),
    ("BRCA1", "BRCA2", "comparator"), ("STK11", "SIK1", "comparator"),
]

# Harle 27 screened lines (Table S5 column order); C092 has no DepMap ID
HARLE_LINES = ["A-375", "A2058", "A549", "AsPC-1", "BxPC-3", "C092",
               "Capan-1", "CFPAC-1", "CHL-1", "COR-L23", "EBC-1",
               "HPAF-II", "KP-1N", "KP4", "LCLC-97TM1", "LK-2", "MeWo",
               "MIA PaCa-2", "NCI-H1299", "NCI-H1568", "NCI-H1975",
               "NCI-H23", "SK-MEL-2", "SK-MEL-28", "SK-MEL-5",
               "SK-MES-1", "SU.86.86"]


# ─────────────────────── download / cache ───────────────────────
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def valid_zip_or_other(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    if path.suffix == ".xlsx":
        try:
            with zipfile.ZipFile(path) as z:
                return z.testzip() is None
        except zipfile.BadZipFile:
            return False
    return True


def download(url: str, dest: Path) -> None:
    print(f"    downloading {dest.name} ...")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)


def ensure_inputs() -> dict:
    """Verify cached raw files; download Flister files if missing/corrupt."""
    prov = {}
    if not HARLE_XLSX.exists():
        raise FileNotFoundError(
            f"Harle et al. supplementary file missing: {HARLE_XLSX}\n"
            "This file was pre-cached by the project (Genome Biology "
            "10.1186/s13059-025-03737-w, Additional file 1). Restore it "
            "into data/external/ before running.")
    prov["harle2025_moesm1.xlsx"] = {
        "path": str(HARLE_XLSX), "sha256": sha256(HARLE_XLSX),
        "size_bytes": HARLE_XLSX.stat().st_size,
        "origin": "pre-cached project download (Genome Biology additional "
                  "file, DOI 10.1186/s13059-025-03737-w)"}
    for name, url in FLISTER_REMOTE.items():
        dest = EXT / name
        if not valid_zip_or_other(dest):
            download(url, dest)
        if name in FLISTER_REQUIRED and not valid_zip_or_other(dest):
            raise RuntimeError(f"could not obtain valid file: {name}")
        if dest.exists():
            prov[name] = {"path": str(dest), "sha256": sha256(dest),
                          "size_bytes": dest.stat().st_size,
                          "origin": "bioRxiv 10.1101/2024.07.16.603642 v3 "
                                    "supplementary (preprint of Cell Reports "
                                    "10.1016/j.celrep.2025.116512)"}
    return prov


# ─────────────────── vectorized |DD| scoring engine ───────────────────
def maxabs_dd_scores(pairs, dep, mat, lin_map, min_n=MIN_N, tag=""):
    """
    Mirror of in4mer_benchmark.pair_score(a, b, min_n), vectorized:
    for each unordered pair, max |DD| across BOTH orientations x all
    lineages; a lineage stratum is evaluable when it has >= min_n mutant
    and >= min_n WT lines for the driver (and >= min_n non-NaN dependency
    values in each arm). DD = mean(dep[paralog] | driver-mut) -
    mean(dep[paralog] | driver-WT)  (signed; |DD| used downstream).
    Returns dict {(a, b): dd or None}.
    """
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

    best = {}
    t0 = time.time()
    for gi, (drv, partners) in enumerate(sorted(partners_of.items())):
        if tag and (gi + 1) % 2000 == 0:
            print(f"      [{tag}] driver {gi + 1}/{len(partners_of)} "
                  f"({time.time() - t0:.0f}s)")
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
            mu_m = np.full(len(cols), np.nan)
            mu_w = np.full(len(cols), np.nan)
            mu_m[valid] = np.nansum(ms[:, valid], axis=0) / nvm[valid]
            mu_w[valid] = np.nansum(ws[:, valid], axis=0) / nvw[valid]
            dd = mu_m - mu_w
            for k, (p, _) in enumerate(pj):
                if not valid[k] or np.isnan(dd[k]):
                    continue
                key = frozenset((drv, p))
                if key not in best or abs(dd[k]) > abs(best[key]):
                    best[key] = float(dd[k])
    return {pair: best.get(frozenset(pair)) for pair in pairs}


def check_engine_fidelity(dep, mat, lin_map) -> dict:
    """Re-score the frozen in4mer benchmark pairs and demand agreement."""
    frozen = pd.read_csv(OUTPUT_DIR / "in4mer_benchmark.csv")
    pairs = [(r.gene_a, r.gene_b) for r in frozen.itertuples()]
    print(f"  engine fidelity check on {len(pairs)} frozen in4mer pairs ...")
    scores = maxabs_dd_scores(pairs, dep, mat, lin_map, tag="fidelity")
    max_diff, n_both_none, n_mismatch = 0.0, 0, 0
    mismatches = []
    for r in frozen.itertuples():
        ref = r.dd_min3 if pd.notna(r.dd_min3) else None
        got = scores[(r.gene_a, r.gene_b)]
        if ref is None and got is None:
            n_both_none += 1
        elif ref is None or got is None:
            n_mismatch += 1
            mismatches.append({"pair": r.pair, "ref": ref, "got": got})
        else:
            d = abs(ref - got)
            max_diff = max(max_diff, d)
            if d > 1e-8:
                n_mismatch += 1
                mismatches.append({"pair": r.pair, "ref": ref, "got": got})
    out = {"n_pairs": len(pairs), "n_both_none": n_both_none,
           "max_abs_diff": max_diff, "n_mismatch": n_mismatch,
           "mismatches": mismatches[:10]}
    print(f"    fidelity: max|diff|={max_diff:.2e}, mismatches={n_mismatch}")
    if n_mismatch:
        raise RuntimeError(
            "DD engine failed fidelity check against frozen "
            "in4mer_benchmark.csv — refusing to continue: "
            f"{mismatches[:5]}")
    return out


# ─────────────────────────── statistics ───────────────────────────
def auroc_with_permutation(y, s, seed=SEED, n_perm=N_PERM):
    """AUROC + label-permutation p (mirrors in4mer_benchmark).

    AUROC uses tie-averaged ranks (Mann-Whitney form), numerically
    identical to sklearn roc_auc_score — asserted on the observed
    statistic (< 1e-12). The permutation null reuses the fixed rank
    vector and sums random n_pos-sized subsets, which is exactly the
    label-permutation distribution but ~1000x faster for large universes.
    """
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    n_pos = int(y.sum())
    n = len(y)
    ranks = rankdata(s)                      # average ranks for ties
    auc_rank = float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2)
                     / (n_pos * (n - n_pos)))
    auc_sk = float(roc_auc_score(y, s))
    assert abs(auc_rank - auc_sk) < 1e-12, (auc_rank, auc_sk)
    rng = np.random.default_rng(seed)
    rank_sums = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.permutation(n)[:n_pos]
        rank_sums[i] = ranks[idx].sum()
    null = (rank_sums - n_pos * (n_pos + 1) / 2) / (n_pos * (n - n_pos))
    return {
        "auroc": auc_sk,
        "permutation_p": float((np.sum(null >= auc_sk) + 1) / (n_perm + 1)),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std()),
        "n_perm": n_perm,
        "seed": seed,
    }


def bh_qvalues(pvals):
    """Benjamini-Hochberg q-values, input order preserved."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order]
    cum = ranked * n / (np.arange(n) + 1)
    cum = np.minimum.accumulate(cum[::-1])[::-1]
    q[order] = np.minimum(cum, 1.0)
    return q


def stratified_test(df, hit_col, gi_col):
    """Fisher (hit rate mut > wt, pre-registered direction) + Mann-Whitney
    (GI more negative in mutants = stronger SL). One-sided primary tests
    plus two-sided versions."""
    mut = df[df["driver_mut"] == 1]
    wt = df[df["driver_mut"] == 0]
    res = {"n_mut": int(len(mut)), "n_wt": int(len(wt)),
           "n_mut_hit": int(mut[hit_col].sum()),
           "n_wt_hit": int(wt[hit_col].sum())}
    if len(mut) == 0 or len(wt) == 0:
        res["evaluable"] = False
        return res
    res["evaluable"] = True
    res["hit_rate_mut"] = float(mut[hit_col].mean())
    res["hit_rate_wt"] = float(wt[hit_col].mean())
    table = [[int(mut[hit_col].sum()), int((1 - mut[hit_col]).sum())],
             [int(wt[hit_col].sum()), int((1 - wt[hit_col]).sum())]]
    res["fisher_p_1sided_greater"] = float(
        fisher_exact(table, alternative="greater")[1])
    res["fisher_p_2sided"] = float(fisher_exact(table)[1])
    gm = mut[gi_col].dropna().to_numpy()
    gw = wt[gi_col].dropna().to_numpy()
    res["gi_mean_mut"] = float(np.mean(gm)) if len(gm) else None
    res["gi_mean_wt"] = float(np.mean(gw)) if len(gw) else None
    res["gi_median_mut"] = float(np.median(gm)) if len(gm) else None
    res["gi_median_wt"] = float(np.median(gw)) if len(gw) else None
    if len(gm) >= 1 and len(gw) >= 1:
        res["mw_p_1sided_mut_more_negative"] = float(
            mannwhitneyu(gm, gw, alternative="less")[1])
        res["mw_p_2sided"] = float(
            mannwhitneyu(gm, gw, alternative="two-sided")[1])
    return res


# ─────────────────────── external data loading ───────────────────────
def load_harle():
    s5_pq = CACHE / "harle_tableS5.pkl"
    s4_pq = CACHE / "harle_tableS4.pkl"
    if s5_pq.exists() and s4_pq.exists():
        return pd.read_pickle(s5_pq), pd.read_pickle(s4_pq)
    print("  parsing Harle xlsx (Table S5 + S4, one-time cache) ...")
    s5 = pd.read_excel(HARLE_XLSX, sheet_name="Table S5", header=3)
    s4 = pd.read_excel(HARLE_XLSX, sheet_name="Table S4", header=4)
    s5.to_pickle(s5_pq)
    s4.to_pickle(s4_pq)
    return s5, s4


def load_flister_s3():
    pq = CACHE / "flister_tableS3_sheet1.pkl"
    if pq.exists():
        return pd.read_pickle(pq)
    df = pd.read_excel(EXT / "flister2025_biorxiv_suppl_tableS3.xlsx",
                       sheet_name="Sheet1")
    df.to_pickle(pq)
    return df


def load_flister_s6():
    pq = CACHE / "flister_tableS6_sl_all_combat.pkl"
    if pq.exists():
        return pd.read_pickle(pq)
    df = pd.read_excel(EXT / "flister2025_biorxiv_suppl_tableS6.xlsx",
                       sheet_name="SL_all_combat")
    df.to_pickle(pq)
    return df


# ─────────────────────── pair-set helpers ───────────────────────
def harle_pair_frame(harle_s5):
    pairs = [tuple(p.split("|")) for p in
             harle_s5["sorted_gene_pair"].astype(str)]
    hit = (harle_s5[HARLE_LINES].sum(axis=1) >= 1).astype(int).to_numpy()
    return pairs, hit


def flister_pair_frame(fl_s3):
    pairs = [tuple(p.split("_")) for p in fl_s3["label"].astype(str)]
    h1299 = (fl_s3["NCIH1299"] == "Lethal").astype(int)
    mda = (fl_s3["MDAMB231"] == "Lethal").astype(int)
    union = (h1299 | mda).astype(int)
    return pairs, h1299.to_numpy(), mda.to_numpy(), union.to_numpy()


# ─────────────────────── stage: inputs ───────────────────────
def stage_inputs():
    print("\n=== stage: inputs ===")
    provenance = ensure_inputs()
    harle_s5, harle_s4 = load_harle()
    fl_s3 = load_flister_s3()
    fl_s6 = load_flister_s6()
    pairs_h, hit_h = harle_pair_frame(harle_s5)
    pairs_f, h1299, mda, union = flister_pair_frame(fl_s3)
    meta = {
        "provenance": provenance,
        "harle": {"n_pairs": len(pairs_h),
                  "n_hit_pairs": int(hit_h.sum()),
                  "n_pairline_hits_S5": int(
                      harle_s5[HARLE_LINES].to_numpy().sum()),
                  "n_pairline_hits_S4": int(harle_s4["is_bassik_hit"].sum()),
                  "n_lines_S4_with_depmapid": int(
                      harle_s4["depMapID"].notna().sum())},
        "flister": {"n_pairs_with_calls": len(pairs_f),
                    "n_lethal_h1299": int(h1299.sum()),
                    "n_lethal_mdamb231": int(mda.sum()),
                    "n_lethal_both": int(((h1299 == 1) & (mda == 1)).sum()),
                    "n_lethal_union": int(union.sum()),
                    "n_s6_rows": int(len(fl_s6))},
    }
    print(f"  Harle: {meta['harle']['n_pairs']} pairs, "
          f"{meta['harle']['n_hit_pairs']} hit pairs (paper: 117), "
          f"{meta['harle']['n_pairline_hits_S5']} pair-line hits (S4: "
          f"{meta['harle']['n_pairline_hits_S4']})")
    print(f"  Flister: {meta['flister']['n_pairs_with_calls']} pairs with "
          f"calls; Lethal H1299={meta['flister']['n_lethal_h1299']}, "
          f"MDAMB231={meta['flister']['n_lethal_mdamb231']}, "
          f"both={meta['flister']['n_lethal_both']}, "
          f"union={meta['flister']['n_lethal_union']}")
    with open(PARTIAL_DIR / "inputs.json", "w") as f:
        json.dump(meta, f, indent=2)


# ─────────────────────── stage: depmap ───────────────────────
def stage_depmap():
    print("\n=== stage: depmap (DepMap load + mutation matrix) ===")
    t0 = time.time()
    dep = load_dependency()
    models = load_models()
    mut = load_mutations()
    assert dep.columns.is_unique, "duplicate gene symbols in dep matrix"
    cell_lines = list(dep.index)
    lin_map = models.set_index("DepMap_ID")[LINEAGE_COL].reindex(cell_lines)

    harle_s5, _ = load_harle()
    fl_s3 = load_flister_s3()
    pairs_h, _ = harle_pair_frame(harle_s5)
    pairs_f, _, _, _ = flister_pair_frame(fl_s3)
    frozen_in4mer = pd.read_csv(OUTPUT_DIR / "in4mer_benchmark.csv")
    all_genes = sorted(
        {g for p in pairs_h + pairs_f for g in p}
        | {d for d, _, _ in CURATED}
        | set(frozen_in4mer["gene_a"]) | set(frozen_in4mer["gene_b"]))
    print(f"  building mutation matrix: {len(all_genes)} genes x "
          f"{len(cell_lines)} lines ...")
    mat = build_mutation_matrix(mut, cell_lines, all_genes,
                                apply_driver_rules=True)
    # driver-only matrix over ALL Model.csv lines (Layer-2 genotyping of
    # external cell lines, including lines absent from the dep matrix)
    driver_genes = sorted({d for d, _, _ in CURATED})
    all_model_ids = list(models["DepMap_ID"])
    mat_drivers = build_mutation_matrix(mut, all_model_ids, driver_genes,
                                        apply_driver_rules=True)
    print(f"  saving state ... ({time.time() - t0:.0f}s so far)")
    with open(STATE_PKL, "wb") as f:
        pd.to_pickle({"dep": dep, "models": models, "lin_map": lin_map,
                      "mat": mat, "mat_drivers": mat_drivers}, f)
    print(f"  state saved to {STATE_PKL} ({time.time() - t0:.0f}s total)")


def load_state():
    return pd.read_pickle(STATE_PKL)


# ─────────────────────── stage: harle ───────────────────────
def stage_harle():
    print("\n=== stage: HARLE (fidelity + Layer 1 + Layer 2) ===")
    st = load_state()
    dep, mat, lin_map = st["dep"], st["mat"], st["lin_map"]
    mat_drivers = st["mat_drivers"]
    harle_s5, harle_s4 = load_harle()
    pairs, hit = harle_pair_frame(harle_s5)
    notes = []

    fidelity = check_engine_fidelity(dep, mat, lin_map)

    print("\n  HARLE Layer 1 (mutation-agnostic AUROC) ...")
    scores = maxabs_dd_scores(pairs, dep, mat, lin_map, tag="harle")
    h_df = pd.DataFrame({
        "pair": ["|".join(sorted(p)) for p in pairs],
        "gene_a": [p[0] for p in pairs],
        "gene_b": [p[1] for p in pairs],
        "hit_any_line": hit,
        "dd_min3": [scores[p] for p in pairs],
    })
    h_df["abs_dd"] = h_df["dd_min3"].abs()
    h_eval = h_df.dropna(subset=["abs_dd"])
    harle_l1 = auroc_with_permutation(
        h_eval["hit_any_line"].to_numpy(), h_eval["abs_dd"].to_numpy())
    harle_l1.update({
        "n_pairs_screen": len(h_df), "n_pairs_scored": len(h_eval),
        "n_hits": int(h_eval["hit_any_line"].sum()),
        "n_nonhits": int((1 - h_eval["hit_any_line"]).sum()),
        "median_abs_dd_hits": float(
            h_eval.loc[h_eval["hit_any_line"] == 1, "abs_dd"].median()),
        "median_abs_dd_nonhits": float(
            h_eval.loc[h_eval["hit_any_line"] == 0, "abs_dd"].median()),
    })
    print(f"    scored {len(h_eval)}/{len(h_df)} pairs; AUROC="
          f"{harle_l1['auroc']:.4f}, perm p={harle_l1['permutation_p']:.4f}")
    h_df.to_csv(OUTPUT_DIR / "external_holdout_harle_layer1.csv", index=False)

    print("\n  HARLE Layer 2 (genotype-stratified) ...")
    s4 = harle_s4[harle_s4["depMapID"].notna()].copy()
    harle_ids = sorted(s4["depMapID"].unique())
    print(f"    {len(harle_ids)} lines with DepMap ID (C092 excluded, "
          "not in DepMap)")
    notes.append("Harle line C092 has no DepMap ID and is excluded from "
                 "genotype stratification (26/27 lines used).")
    harle_l2_rows, harle_l2_line_rows = [], []
    harle_keyed = {frozenset((d, p)): (d, p, t) for d, p, t in CURATED}
    for key, (drv, prl, tier) in harle_keyed.items():
        sp = "|".join(sorted(key))
        sub = s4[s4["sorted_gene_pair"] == sp]
        if sub.empty:
            harle_l2_rows.append({"driver": drv, "paralog": prl,
                                  "tier": tier, "in_screen": False,
                                  "evaluable": False})
            continue
        sub = sub[["sorted_gene_pair", "cell_line_label", "cancer_type",
                   "depMapID", "is_bassik_hit", "mean_norm_gi",
                   "median_norm_gi", "gi_t_score", "fdr"]].copy()
        sub["driver"] = drv
        sub["paralog"] = prl
        geno = sub["depMapID"].map(mat_drivers[drv])
        if geno.isna().any():
            raise RuntimeError(f"Harle line(s) without genotype for {drv}: "
                               f"{sub.loc[geno.isna(), 'depMapID'].unique()}")
        sub["driver_mut"] = geno.astype(int)
        sub = sub.rename(columns={"sorted_gene_pair": "pair",
                                  "is_bassik_hit": "hit"})
        harle_l2_line_rows.append(sub)
        res = stratified_test(sub, "hit", "mean_norm_gi")
        hit_lines = sorted(sub.loc[sub["hit"] == 1, "cell_line_label"])
        mut_hit_lines = sorted(
            sub.loc[(sub["hit"] == 1) & (sub["driver_mut"] == 1),
                    "cell_line_label"])
        mut_lines = sorted(sub.loc[sub["driver_mut"] == 1, "cell_line_label"])
        harle_l2_rows.append({
            "driver": drv, "paralog": prl, "tier": tier, "in_screen": True,
            "driver_mutant_lines": ";".join(mut_lines),
            "hit_lines_all": ";".join(hit_lines),
            "hit_lines_driver_mutant": ";".join(mut_hit_lines),
            **res})
    harle_l2 = pd.DataFrame(harle_l2_rows)
    ev = harle_l2["evaluable"] == True
    if ev.sum() >= 2:
        harle_l2.loc[ev, "fisher_q_BH"] = bh_qvalues(
            harle_l2.loc[ev, "fisher_p_1sided_greater"])
        harle_l2.loc[ev, "mw_q_BH"] = bh_qvalues(
            harle_l2.loc[ev, "mw_p_1sided_mut_more_negative"])
    harle_l2.to_csv(OUTPUT_DIR / "external_holdout_harle_layer2.csv",
                    index=False)
    harle_lines_df = (pd.concat(harle_l2_line_rows, ignore_index=True)
                      if harle_l2_line_rows else pd.DataFrame())
    if not harle_lines_df.empty:
        harle_lines_df.to_csv(
            OUTPUT_DIR / "external_holdout_harle_layer2_lines.csv",
            index=False)
        harle_pooled = stratified_test(harle_lines_df, "hit", "mean_norm_gi")
    else:
        harle_pooled = {"evaluable": False}
    for r in harle_l2_rows:
        if r.get("in_screen") and r.get("evaluable"):
            print(f"    {r['driver']}->{r['paralog']}: "
                  f"mut {r['n_mut_hit']}/{r['n_mut']} hit vs wt "
                  f"{r['n_wt_hit']}/{r['n_wt']} hit, "
                  f"Fisher p1={r['fisher_p_1sided_greater']:.4g}, "
                  f"MW p1={r.get('mw_p_1sided_mut_more_negative', float('nan')):.4g}")
        elif r.get("in_screen"):
            print(f"    {r['driver']}->{r['paralog']}: in screen but NOT "
                  f"evaluable (mut={r.get('n_mut', 0)}, wt={r.get('n_wt', 0)})")
        else:
            print(f"    {r['driver']}->{r['paralog']}: not in Harle screen")
    out = {"engine_fidelity_vs_frozen_in4mer": fidelity,
           "layer1": harle_l1,
           "layer2": {
               "n_lines_genotyped": len(harle_ids),
               "pairs_in_screen": int(harle_l2["in_screen"].sum()),
               "pairs_evaluable": int(ev.sum()),
               "per_pair": harle_l2.to_dict(orient="records"),
               "pooled": harle_pooled},
           "notes": notes}
    with open(PARTIAL_DIR / "harle.json", "w") as f:
        json.dump(out, f, indent=2, default=str)


# ─────────────────────── stage: flister_l1 ───────────────────────
def stage_flister_l1():
    print("\n=== stage: FLISTER Layer 1 (mutation-agnostic AUROC) ===")
    st = load_state()
    dep, mat, lin_map = st["dep"], st["mat"], st["lin_map"]
    fl_s3 = load_flister_s3()
    pairs, h1299, mda, union = flister_pair_frame(fl_s3)

    scores = maxabs_dd_scores(pairs, dep, mat, lin_map, tag="flister")
    f_df = pd.DataFrame({
        "pair": fl_s3["label"].astype(str),
        "gene_a": [p[0] for p in pairs],
        "gene_b": [p[1] for p in pairs],
        "hit_union": union,
        "hit_h1299": h1299,
        "hit_mdamb231": mda,
        "dd_min3": [scores[p] for p in pairs],
    })
    f_df["abs_dd"] = f_df["dd_min3"].abs()
    f_eval = f_df.dropna(subset=["abs_dd"])
    fl_l1 = auroc_with_permutation(
        f_eval["hit_union"].to_numpy(), f_eval["abs_dd"].to_numpy())
    fl_l1.update({
        "n_pairs_screen": len(f_df), "n_pairs_scored": len(f_eval),
        "n_hits": int(f_eval["hit_union"].sum()),
        "n_nonhits": int((1 - f_eval["hit_union"]).sum()),
        "median_abs_dd_hits": float(
            f_eval.loc[f_eval["hit_union"] == 1, "abs_dd"].median()),
        "median_abs_dd_nonhits": float(
            f_eval.loc[f_eval["hit_union"] == 0, "abs_dd"].median()),
    })
    print(f"    scored {len(f_eval)}/{len(f_df)} pairs; AUROC="
          f"{fl_l1['auroc']:.4f}, perm p={fl_l1['permutation_p']:.4f}")
    fl_l1["sensitivity_per_model"] = {}
    for col, name in (("hit_h1299", "NCI-H1299 only"),
                      ("hit_mdamb231", "MDA-MB-231 only")):
        r = auroc_with_permutation(f_eval[col].to_numpy(),
                                   f_eval["abs_dd"].to_numpy())
        r["n_hits"] = int(f_eval[col].sum())
        fl_l1["sensitivity_per_model"][col] = {"label": name, **r}
        print(f"      {name}: AUROC={r['auroc']:.4f}, "
              f"p={r['permutation_p']:.4f} (hits={int(f_eval[col].sum())})")
    f_df.to_csv(OUTPUT_DIR / "external_holdout_flister_layer1.csv",
                index=False)
    with open(PARTIAL_DIR / "flister_l1.json", "w") as f:
        json.dump({"layer1": fl_l1}, f, indent=2)


# ─────────────────────── stage: flister_l2 + final JSON ───────────────────────
def stage_flister_l2():
    print("\n=== stage: FLISTER Layer 2 + final JSON assembly ===")
    st = load_state()
    mat, models = st["mat"], st["models"]
    mat_drivers = st["mat_drivers"]
    fl_s6 = load_flister_s6()
    notes = []

    model_cols = [c for c in fl_s6.columns
                  if c.endswith(("_merged", "_AbbVie", "_Vakoc", "_Sellers"))]
    s6_all = fl_s6.set_index("Unnamed: 0")
    s6_vals = s6_all[model_cols].apply(pd.to_numeric, errors="coerce")
    name2id = dict(zip(models["StrippedCellLineName"].astype(str).str.upper(),
                       models["DepMap_ID"]))
    col2id, col2study = {}, {}
    for c in model_cols:
        base, study = c.rsplit("_", 1)
        col2id[c] = name2id.get(base.upper())
        col2study[c] = study
    unmapped = [c for c, i in col2id.items() if i is None]
    print(f"    {len(model_cols)} model columns; unmapped to DepMap: "
          f"{unmapped if unmapped else 'none'}")
    n_nonnum = int(s6_all[model_cols].notna().sum().sum()
                   - s6_vals.notna().sum().sum())
    if n_nonnum:
        notes.append(f"Flister S6: {n_nonnum} non-numeric cells coerced to "
                     "NaN (blank placeholders).")
    pred_sl = (s6_vals <= -2).sum(axis=1)
    concord = float((pred_sl == s6_all["# SL"]).mean())
    print(f"    dzL2FC<=-2 vs their '# SL' concordance: {concord:.3f}")
    notes.append(f"Flister S6 hit call dzL2FC<=-2 reproduces their '# SL' "
                 f"for {concord:.1%} of pairs (their methods: dzL2FC>2 "
                 "significant).")

    s6_index = {}
    for lbl in s6_vals.index.astype(str):
        parts = lbl.split("_")
        if len(parts) == 2:
            s6_index[frozenset(parts)] = lbl
    fl_l2_rows, fl_l2_model_rows = [], []
    for drv, prl, tier in CURATED:
        key = frozenset((drv, prl))
        if key not in s6_index:
            fl_l2_rows.append({"driver": drv, "paralog": prl, "tier": tier,
                               "in_screen": False, "evaluable": False})
            continue
        lbl = s6_index[key]
        row = s6_vals.loc[lbl]
        sub = pd.DataFrame({
            "pair": lbl, "driver": drv, "paralog": prl,
            "model_col": model_cols,
            "study": [col2study[c] for c in model_cols],
            "depMapID": [col2id[c] for c in model_cols],
            "dzL2FC": [row[c] for c in model_cols]})
        sub = sub[sub["depMapID"].notna() & sub["dzL2FC"].notna()].copy()
        sub["hit"] = (sub["dzL2FC"] <= -2).astype(int)
        geno = sub["depMapID"].map(mat_drivers[drv])
        if geno.isna().any():
            raise RuntimeError(f"Flister model(s) without genotype for "
                               f"{drv}: "
                               f"{sub.loc[geno.isna(), 'depMapID'].unique()}")
        sub["driver_mut"] = geno.astype(int)
        fl_l2_model_rows.append(sub)
        res = stratified_test(sub, "hit", "dzL2FC")
        hit_models = sorted(sub.loc[sub["hit"] == 1, "model_col"])
        mut_hit_models = sorted(
            sub.loc[(sub["hit"] == 1) & (sub["driver_mut"] == 1),
                    "model_col"])
        mut_models = sorted(sub.loc[sub["driver_mut"] == 1, "model_col"])
        fl_l2_rows.append({
            "driver": drv, "paralog": prl, "tier": tier, "in_screen": True,
            "s6_label": lbl,
            "driver_mutant_models": ";".join(mut_models),
            "hit_models_all": ";".join(hit_models),
            "hit_models_driver_mutant": ";".join(mut_hit_models),
            **res})
    fl_l2 = pd.DataFrame(fl_l2_rows)
    evf = fl_l2["evaluable"] == True
    if evf.sum() >= 2:
        fl_l2.loc[evf, "fisher_q_BH"] = bh_qvalues(
            fl_l2.loc[evf, "fisher_p_1sided_greater"])
        fl_l2.loc[evf, "mw_q_BH"] = bh_qvalues(
            fl_l2.loc[evf, "mw_p_1sided_mut_more_negative"])
    fl_l2.to_csv(OUTPUT_DIR / "external_holdout_flister_layer2.csv",
                 index=False)
    fl_models_df = (pd.concat(fl_l2_model_rows, ignore_index=True)
                    if fl_l2_model_rows else pd.DataFrame())
    if not fl_models_df.empty:
        fl_models_df.to_csv(
            OUTPUT_DIR / "external_holdout_flister_layer2_models.csv",
            index=False)
        fl_pooled = stratified_test(fl_models_df, "hit", "dzL2FC")
    else:
        fl_pooled = {"evaluable": False}
    for r in fl_l2_rows:
        if r.get("in_screen") and r.get("evaluable"):
            print(f"    {r['driver']}->{r['paralog']}: "
                  f"mut {r['n_mut_hit']}/{r['n_mut']} hit vs wt "
                  f"{r['n_wt_hit']}/{r['n_wt']} hit, "
                  f"Fisher p1={r['fisher_p_1sided_greater']:.4g}, "
                  f"MW p1={r.get('mw_p_1sided_mut_more_negative', float('nan')):.4g}")
        elif r.get("in_screen"):
            print(f"    {r['driver']}->{r['paralog']}: in S6 but NOT "
                  f"evaluable (mut={r.get('n_mut', 0)}, wt={r.get('n_wt', 0)})")
        else:
            print(f"    {r['driver']}->{r['paralog']}: not in Flister S6")

    # ── assemble final JSON from partials ──
    with open(PARTIAL_DIR / "inputs.json") as f:
        inputs_meta = json.load(f)
    with open(PARTIAL_DIR / "harle.json") as f:
        harle = json.load(f)
    with open(PARTIAL_DIR / "flister_l1.json") as f:
        fl_l1 = json.load(f)["layer1"]

    notes.append(
        "Flister Table S3 (bioRxiv v3) gives 117/51 Lethal calls for "
        "NCI-H1299/MDA-MB-231 (22 shared, matching the text's '22'), "
        "whereas the preprint results text cites 144 and 45 pairs; "
        "we use the table values as the authoritative hit calls.")

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED, "n_perm": N_PERM, "min_n_frame": MIN_N,
        "framework": {
            "dd_score": "max |DD| across both orientations x all lineages, "
                        ">=3 mutant / >=3 WT (mirrors in4mer_benchmark.py "
                        "pair_score min_n=3); DD = mean Chronos dep of "
                        "paralog in driver-mutant lines minus WT lines",
            "mutation_rules": "project TableS6 driver rules via "
                              "build_mutation_matrix(apply_driver_rules=True): "
                              "TSG = LikelyLoF, ONC = Hotspot (DepMap 26Q1)",
            "layer2_tests": "Fisher exact (hit rate, 1-sided 'greater') + "
                            "Mann-Whitney (GI, 1-sided mutant-more-negative); "
                            "2-sided p also in CSVs; BH q within each dataset",
        },
        "engine_fidelity_vs_frozen_in4mer":
            harle["engine_fidelity_vs_frozen_in4mer"],
        "data_provenance": inputs_meta["provenance"],
        "harle2025": {
            "screen": "472 predicted SL pairs x 27 cell lines "
                      "(10 melanoma / 9 lung NSCLC / 8 pancreas), "
                      "dual-sgRNA combinatorial CRISPR (Genome Biology 2025, "
                      "DOI 10.1186/s13059-025-03737-w)",
            "hit_definition": "is_bassik_hit: mean_norm_gi < -0.5 & "
                              "fdr < 0.01 & neither single target depleted "
                              "(their code, 07_combine_intermediate_"
                              "datasets.R); pair hit = >=1 hit line",
            "screen_counts": inputs_meta["harle"],
            "layer1": harle["layer1"],
            "layer2": harle["layer2"],
        },
        "flister2025": {
            "screen": "36,648 paralog pairs digenic enAsCas12a screen in "
                      "NCI-H1299 + MDA-MB-231 (Cell Reports 2025, DOI "
                      "10.1016/j.celrep.2025.116512; data from bioRxiv "
                      "10.1101/2024.07.16.603642 v3 supplementary)",
            "hit_definition_layer1": "'Lethal' call in Table S3 Sheet1 for "
                                     "NCIH1299 and/or MDAMB231 (union)",
            "hit_definition_layer2": "dzL2FC <= -2 in Table S6 SL_all_combat "
                                     "(their methods: dzL2FC>2 significant)",
            "n_pairs_screened": 36648,
            "screen_counts": inputs_meta["flister"],
            "layer1": fl_l1,
            "layer2": {
                "n_models": len(model_cols),
                "pairs_in_S6": int(fl_l2["in_screen"].sum()),
                "pairs_evaluable": int(evf.sum()),
                "s6_concordance": concord,
                "per_pair": fl_l2.to_dict(orient="records"),
                "pooled": fl_pooled,
            },
        },
        "notes": harle.get("notes", []) + notes,
    }
    with open(OUTPUT_DIR / "external_holdout.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSaved {OUTPUT_DIR}/external_holdout.json")


# ───────────────────────────── main ─────────────────────────────
STAGES = {"inputs": stage_inputs, "depmap": stage_depmap,
          "harle": stage_harle, "flister_l1": stage_flister_l1,
          "flister_l2": stage_flister_l2}
ORDER = ["inputs", "depmap", "harle", "flister_l1", "flister_l2"]


def main():
    t_start = time.time()
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    print("=" * 72)
    print("  Locked external hold-out validation (Harle 2025 + Flister 2025)")
    print("=" * 72)
    if stage == "all":
        for s in ORDER:
            STAGES[s]()
    elif stage in STAGES:
        STAGES[stage]()
    else:
        raise SystemExit(f"unknown stage '{stage}'; choose from "
                         f"{ORDER} or 'all'")
    print(f"\nTotal runtime: {time.time() - t_start:.0f}s")


if __name__ == "__main__":
    main()
