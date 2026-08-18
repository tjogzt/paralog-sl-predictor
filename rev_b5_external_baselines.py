#!/usr/bin/env python3
"""
rev_b5_external_baselines.py  (Stage-4 revision, item B5)
==========================================================
Baseline-comparison analysis for the external digenic-screen frames
(Harle 2025 Layer 1, Flister 2025 Layer 1): is mutation-agnostic max|DD|
better than trivial monogenic / sequence baselines?

Baselines computed per pair on each screen's evaluable Layer-1 frame
(the frozen output/external_holdout_*_layer1.csv universes; max|DD| taken
from those frozen, engine-verified artifacts):

  1. max_abs_mean_chronos  = max(|mean Chronos A|, |mean Chronos B|)
     ("max|mean Chronos|", necessity-type baseline; global mean across all
     1,208 DepMap lines, NaN-skipped)
  2. min_abs_mean_chronos  = min(|mean Chronos A|, |mean Chronos B|)
     ("min single-gene effect": the pair's weaker single-gene effect; the
     SL-consistent reading is that LOW values rank digenic hits higher --
     both members individually non-depleted, cf. the Harle hit definition)
  3. mrna_coexpression     = Pearson r of log2(TPM+1) between the two genes
     across DepMap lines (pairwise complete)
  4. kmer_identity         = k-mer Jaccard (k=3) -- NOT COMPUTABLE locally:
     data/uniprot_sequences.rds covers only the 26 main-universe genes;
     data/ensembl_paralogs.csv identity_pct is uniformly NA. Coverage is
     quantified and the item recorded as NOT COMPUTABLE.
  5. De Kegel 2021 public scores -- no local copy exists in data/ or
     output/ (verified); recorded as NOT COMPUTABLE.

Each baseline: AUROC + label-permutation p (10,000, seed 42). Then
delta-AUROC (max|DD| - best baseline) with a paired bootstrap 95% CI
(10,000 pair resamples with replacement, same draws for both scores,
seed 42). All statistics on identical evaluable universes per comparison.

Outputs (output/revision_stage4/):
  b5_harle_baseline_scores.csv, b5_flister_baseline_scores.csv
  b5_external_baselines.json

Usage: python rev_b5_external_baselines.py   (run from repo root)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

from external_holdout import auroc_with_permutation, STATE_PKL  # noqa: E402
from compute_headline_metrics import auroc  # noqa: E402
from config import DEPMAP_FILES  # noqa: E402

SEED = 42
N_PERM = 10_000
N_BOOT = 10_000


def paired_boot_delta(y, s_ref, s_base, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(y)
    d = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yb = y[idx]
        if yb.sum() == 0 or yb.sum() == n:
            continue
        d.append(auroc(yb, s_ref[idx]) - auroc(yb, s_base[idx]))
    d = np.asarray(d)
    return {"mean_delta": float(d.mean()),
            "ci95": [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
            "frac_delta_below_0": float((d < 0).mean()),
            "n_valid": int(len(d)), "seed": seed}


def eval_frame(frozen_csv, mean_dep, coexpr):
    fr = pd.read_csv(ROOT / "output" / frozen_csv)
    fr = fr.dropna(subset=["dd_min3"]).copy()  # evaluable Layer-1 universe
    fr["abs_dd"] = fr["dd_min3"].abs()
    fr["max_abs_mean_chronos"] = [
        max(abs(mean_dep.get(a, np.nan)), abs(mean_dep.get(b, np.nan)))
        if a in mean_dep and b in mean_dep else np.nan
        for a, b in zip(fr["gene_a"], fr["gene_b"])]
    fr["min_abs_mean_chronos"] = [
        min(abs(mean_dep.get(a, np.nan)), abs(mean_dep.get(b, np.nan)))
        if a in mean_dep and b in mean_dep else np.nan
        for a, b in zip(fr["gene_a"], fr["gene_b"])]
    fr["mrna_coexpression"] = [
        coexpr.get(frozenset((a, b)), np.nan)
        for a, b in zip(fr["gene_a"], fr["gene_b"])]
    return fr


def score_all(fr, hit_col):
    y = fr[hit_col].astype(int).to_numpy()
    out = {}
    for col in ("abs_dd", "max_abs_mean_chronos", "min_abs_mean_chronos", "mrna_coexpression"):
        e = fr.dropna(subset=[col])
        r = auroc_with_permutation(e[hit_col].astype(int).to_numpy(), e[col].to_numpy())
        r["n"] = int(len(e))
        r["n_hits"] = int(e[hit_col].sum())
        out[col] = r
    return out, y


def main():
    print("=" * 72)
    print("  rev B5: external-frame baseline comparison")
    print("=" * 72)
    st = pd.read_pickle(STATE_PKL)
    dep = st["dep"]
    mean_dep = dep.mean(axis=0, skipna=True).to_dict()
    print(f"  global mean Chronos computed for {len(mean_dep)} genes")

    # genes needed for coexpression
    fr_h = pd.read_csv(ROOT / "output" / "external_holdout_harle_layer1.csv")
    fr_f = pd.read_csv(ROOT / "output" / "external_holdout_flister_layer1.csv")
    genes = sorted(set(fr_h["gene_a"]) | set(fr_h["gene_b"])
                   | set(fr_f["gene_a"]) | set(fr_f["gene_b"]))
    print(f"  {len(genes)} unique genes across both external frames")

    # expression matrix, needed columns only
    hdr = pd.read_csv(DEPMAP_FILES["expression"], nrows=0)
    id_col = hdr.columns[0]
    col_of = {}
    for c in hdr.columns[1:]:
        g = c.split(" ")[0] if " (" in c else c
        if g in set(genes):
            col_of[g] = c
    usecols = [id_col] + list(col_of.values())
    expr = pd.read_csv(DEPMAP_FILES["expression"], usecols=usecols).set_index(id_col)
    expr.columns = [g for g, c in col_of.items()]
    print(f"  expression slice: {expr.shape}")

    def coexpr_of(a, b):
        if a not in expr.columns or b not in expr.columns:
            return np.nan
        x, z = expr[a], expr[b]
        ok = x.notna() & z.notna()
        if ok.sum() < 10:
            return np.nan
        return float(np.corrcoef(x[ok], z[ok])[0, 1])

    uniq_pairs = {frozenset((a, b)) for a, b in
                  pd.concat([fr_h[["gene_a", "gene_b"]], fr_f[["gene_a", "gene_b"]]]).values}
    print(f"  computing mRNA coexpression for {len(uniq_pairs)} unique pairs ...")
    coexpr = {p: coexpr_of(*tuple(p)) for p in uniq_pairs}

    # ── evaluate ──────────────────────────────────────────────────
    results = {}
    for name, frozen_csv, hit_col in (
            ("harle", "external_holdout_harle_layer1.csv", "hit_any_line"),
            ("flister", "external_holdout_flister_layer1.csv", "hit_union")):
        fr = eval_frame(frozen_csv, mean_dep, coexpr)
        fr.to_csv(OUT / f"b5_{name}_baseline_scores.csv", index=False)
        scored, _ = score_all(fr, hit_col)
        # best baseline = highest AUROC among non-DD baselines
        baselines = {k: v for k, v in scored.items() if k != "abs_dd"}
        best_name = max(baselines, key=lambda k: baselines[k]["auroc"])
        e_ref = fr.dropna(subset=["abs_dd", best_name])
        delta = paired_boot_delta(e_ref[hit_col].astype(int).to_numpy(),
                                  e_ref["abs_dd"].to_numpy(), e_ref[best_name].to_numpy())
        results[name] = {
            "n_evaluable": int(len(fr)), "hit_column": hit_col,
            "auroc": {k: {"auroc": v["auroc"], "permutation_p": v["permutation_p"],
                          "n": v["n"], "n_hits": v["n_hits"]}
                      for k, v in scored.items()},
            "best_baseline": best_name,
            "delta_auroc_maxabsDD_minus_best_baseline": scored["abs_dd"]["auroc"]
                                                        - baselines[best_name]["auroc"],
            "paired_bootstrap_delta": delta,
            "coexpression_coverage": float(fr["mrna_coexpression"].notna().mean()),
        }
        print(f"\n  [{name}] n={len(fr)}")
        for k, v in scored.items():
            print(f"    {k:24s} AUROC={v['auroc']:.4f} (p={v['permutation_p']:.4f})")
        print(f"    best baseline: {best_name}; "
              f"delta max|DD|-best = {results[name]['delta_auroc_maxabsDD_minus_best_baseline']:+.4f} "
              f"CI [{delta['ci95'][0]:+.4f}, {delta['ci95'][1]:+.4f}]")

    # ── k-mer coverage check (expected: ~0) ───────────────────────
    import subprocess
    r = subprocess.run(["/usr/local/bin/Rscript", "-e",
                        'cat(paste(names(readRDS("data/uniprot_sequences.rds")), collapse=","))'],
                       capture_output=True, text=True, cwd=ROOT)
    seq_genes = set(r.stdout.strip().split(",")) if r.returncode == 0 else set()
    cov = {}
    for name, fr in (("harle", fr_h), ("flister", fr_f)):
        both = [frozenset((a, b)) <= seq_genes
                for a, b in zip(fr["gene_a"], fr["gene_b"])]
        cov[name] = {"n_pairs": int(len(fr)), "n_pairs_both_genes_sequenced": int(sum(both)),
                     "fraction": float(np.mean(both))}

    results["not_computable"] = {
        "kmer_identity": {
            "reason": "bundled sequence cache data/uniprot_sequences.rds covers only "
                      f"{len(seq_genes)} main-universe genes; data/ensembl_paralogs.csv "
                      "identity_pct is NA for all 66,595 pairs; computing k-mer identity "
                      "for the external-frame genes would require downloading new sequence "
                      "data (out of scope of the artifact-only revision).",
            "coverage_check": cov},
        "de_kegel_public_scores": {
            "reason": "no local copy of De Kegel et al. 2021 paralog-SL prediction scores "
                      "exists in data/ or output/ (verified 2026-08-09); fetching external "
                      "data is out of scope of the artifact-only revision."},
    }

    results["meta"] = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "seed": SEED, "n_perm": N_PERM, "n_boot": N_BOOT,
        "frames": "frozen external_holdout Layer-1 evaluable universes; max|DD| from the "
                  "engine-verified frozen layer1 CSVs",
        "baseline_definitions": {
            "max_abs_mean_chronos": "max(|mean Chronos A|, |mean Chronos B|), global mean "
                                    "across all DepMap lines (necessity-type)",
            "min_abs_mean_chronos": "min(|mean Chronos A|, |mean Chronos B|) -- the pair's "
                                    "weaker single-gene effect; SL-consistent direction is "
                                    "LOW values ranking hits higher",
            "mrna_coexpression": "Pearson r of log2(TPM+1) across DepMap lines (pairwise "
                                 "complete, >=10 lines)",
        },
    }
    (OUT / "b5_external_baselines.json").write_text(
        json.dumps(results, indent=2, default=str))
    print(f"\n  wrote {OUT}/b5_external_baselines.json + per-frame score CSVs")


if __name__ == "__main__":
    main()
