#!/usr/bin/env python3
"""
in4mer_seed_sensitivity.py
==========================
Control-sampling seed sensitivity for the in4mer external benchmark.

in4mer_benchmark.py samples 400 unlabeled control paralog pairs with
seed 42 and scores every pair (13 in4mer gold standards + 400 controls)
by max |DD| across orientations x lineages. The frozen outputs
(output/in4mer_benchmark.json/.csv) are cited by the audit and are NOT
modified here. This script re-runs ONLY the >=3-mutant/>=3-WT (min3)
sensitivity frame with 20 different control-sampling seeds (42-61),
recomputing the in4mer AUROC each time, to quantify how much the
reported min3 AUROC depends on the control draw.

Fidelity to the frozen benchmark (same code path, copied verbatim from
in4mer_benchmark.py; constants imported from it):
  * same control pool: Ensembl/HGNC paralog pairs with both genes in the
    dependency matrix, excluding KNOWN_PARALOG_SL, MEK1/MEK2 aliases and
    the 13 in4mer pairs;
  * same sampler: np.random.default_rng(seed).choice(len(pa), size=400,
    replace=False) — seed 42 therefore reproduces the frozen control set
    exactly (verified against output/in4mer_benchmark_summary.json and
    output/in4mer_benchmark.csv as a consistency check);
  * same score: max |DD| across BOTH orientations x all lineages,
    requiring >=3 mutant and >=3 WT lines in the lineage stratum;
  * same metric: sklearn roc_auc_score on |DD| with in4mer_gold = 1.

The mutation matrix is built ONCE for the union of genes across all 20
control draws (the build is per-gene independent, so this is numerically
identical to building it per seed).

Output: output/in4mer_seed_sensitivity.json
Usage:  python3 in4mer_seed_sensitivity.py   (run from repo root)
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import OUTPUT_DIR, KNOWN_PARALOG_SL
from data_loader import (
    load_dependency,
    load_models,
    load_mutations,
    build_mutation_matrix,
    load_paralogs,
)
from in4mer_benchmark import IN4MER_PAIRS, N_CONTROLS, LINEAGE_COL

SEEDS = list(range(42, 62))  # 42..61 inclusive -> 20 seeds
MIN_N = 3                    # min3 frame only (>=3 mutant / >=3 WT)
FROZEN_SUMMARY = OUTPUT_DIR / "in4mer_benchmark_summary.json"
FROZEN_CSV = OUTPUT_DIR / "in4mer_benchmark.csv"
JSON_OUT = OUTPUT_DIR / "in4mer_seed_sensitivity.json"


def main():
    print("=" * 70)
    print("  in4mer control-sampling seed sensitivity (min3 frame, seeds 42-61)")
    print("=" * 70)

    dep = load_dependency()
    models = load_models()
    mut = load_mutations()
    paralogs = load_paralogs()

    cell_lines = list(dep.index)
    lin_map = models.set_index("DepMap_ID")[LINEAGE_COL].reindex(cell_lines)

    # ── Control pool (identical construction to in4mer_benchmark.py) ──
    exclude = {frozenset(p) for p in KNOWN_PARALOG_SL}
    exclude |= {frozenset({"MEK1", "MEK2"}), frozenset({"MAP2K1", "MAP2K2"})}
    exclude |= {frozenset(p) for p in IN4MER_PAIRS}
    pa = paralogs[["gene_A", "gene_B"]].dropna()
    ok = pa["gene_A"].isin(dep.columns) & pa["gene_B"].isin(dep.columns)
    pa = pa[ok]
    keep = [frozenset((r.gene_A, r.gene_B)) not in exclude for r in pa.itertuples()]
    pa = pa[keep]
    print(f"  Control pool: {len(pa)} paralog pairs")

    # ── 20 seeded control draws ──
    draws = {}
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        take = rng.choice(len(pa), size=min(N_CONTROLS, len(pa)), replace=False)
        draws[seed] = [tuple(x) for x in pa.iloc[take].values]

    # ── One mutation matrix for the union of genes (per-gene independent) ──
    all_genes = set(g for p in IN4MER_PAIRS for g in p)
    for cps in draws.values():
        all_genes |= set(g for p in cps for g in p)
    genes = sorted(all_genes)
    print(f"  Union genes across 20 draws + in4mer: {len(genes)}")
    mat = build_mutation_matrix(mut, cell_lines, genes, apply_driver_rules=True)

    lineages = sorted(set(lin_map.dropna()))
    lin_idx = {lin: lin_map[lin_map == lin].index for lin in lineages}

    def pair_score(a, b, min_n):
        """max |DD| across orientations x lineages; None if not evaluable.
        Verbatim logic from in4mer_benchmark.py."""
        best = None
        n_eval = 0
        for drv, prl in ((a, b), (b, a)):
            if drv not in mat.columns or prl not in dep.columns:
                continue
            for lin in lineages:
                idx = lin_idx[lin]
                m = mat.loc[idx, drv]
                if m.sum() < min_n or (len(m) - m.sum()) < min_n:
                    continue
                d = dep.loc[idx, prl]
                mut_d = d[m == 1].dropna()
                wt_d = d[m == 0].dropna()
                if len(mut_d) < min_n or len(wt_d) < min_n:
                    continue
                n_eval += 1
                dd = float(mut_d.mean() - wt_d.mean())
                if best is None or abs(dd) > abs(best):
                    best = dd
        return best, n_eval

    # Score the 13 in4mer positives ONCE (seed-independent)
    pos_set = {frozenset(p) for p in IN4MER_PAIRS}
    pos_scores = {}
    for a, b in IN4MER_PAIRS:
        dd, n_ev = pair_score(a, b, MIN_N)
        pos_scores[(a, b)] = {"dd_min3": dd, "n_lineages_eval_3": n_ev}
    n_pos_eval = sum(1 for v in pos_scores.values() if v["dd_min3"] is not None)
    print(f"  in4mer positives evaluable on min3 frame: {n_pos_eval}/13 "
          f"(seed-independent)")

    # ── Per-seed scoring and AUROC ──
    per_seed = []
    for seed in SEEDS:
        rows = []
        for (a, b), v in pos_scores.items():
            rows.append({"pair": f"{a}/{b}", "label": "in4mer_gold",
                         "dd_min3": v["dd_min3"]})
        for a, b in draws[seed]:
            dd, _ = pair_score(a, b, MIN_N)
            rows.append({"pair": f"{a}/{b}", "label": "unlabeled_control",
                         "dd_min3": dd})
        df = pd.DataFrame(rows).dropna(subset=["dd_min3"])
        y = (df["label"] == "in4mer_gold").astype(int).values
        s = df["dd_min3"].abs().values
        auc = float(roc_auc_score(y, s)) if y.sum() >= 2 else float("nan")
        per_seed.append({
            "seed": seed,
            "auroc_min3": auc,
            "n_in4mer_evaluable_of_13": int(y.sum()),
            "n_controls_evaluable_of_400": int(len(y) - y.sum()),
        })
        print(f"  seed {seed}: AUROC={auc:.4f} "
              f"(pos {int(y.sum())}/13, ctrl {int(len(y) - y.sum())}/400)")

    aucs = np.array([r["auroc_min3"] for r in per_seed], dtype=float)

    # ── Consistency check: seed 42 must reproduce the frozen benchmark ──
    check = {"frozen_summary": str(FROZEN_SUMMARY.relative_to(OUTPUT_DIR.parent)),
             "frozen_csv": str(FROZEN_CSV.relative_to(OUTPUT_DIR.parent))}
    if FROZEN_SUMMARY.exists():
        frozen = json.loads(FROZEN_SUMMARY.read_text())
        frozen_auc = frozen["min3"]["auroc"]
        seed42 = next(r for r in per_seed if r["seed"] == 42)
        check.update({
            "frozen_min3_auroc": frozen_auc,
            "rerun_seed42_auroc": seed42["auroc_min3"],
            "abs_diff": abs(frozen_auc - seed42["auroc_min3"]),
            "frozen_min3_n_pos": frozen["min3"]["n_pos"],
            "frozen_min3_n_neg": frozen["min3"]["n_neg"],
            "rerun_seed42_n_pos": seed42["n_in4mer_evaluable_of_13"],
            "rerun_seed42_n_neg": seed42["n_controls_evaluable_of_400"],
            "status": "match" if abs(frozen_auc - seed42["auroc_min3"]) < 1e-9 else "MISMATCH",
        })
    if FROZEN_CSV.exists():
        fz = pd.read_csv(FROZEN_CSV)
        fz_pos = fz[fz["label"] == "in4mer_gold"].set_index("pair")["dd_min3"]
        diffs = []
        for (a, b), v in pos_scores.items():
            key = f"{a}/{b}"
            if key in fz_pos.index:
                old, new = fz_pos[key], v["dd_min3"]
                if pd.isna(old) and new is None:
                    continue
                if pd.isna(old) != (new is None) or (new is not None and abs(old - new) > 1e-9):
                    diffs.append(key)
        check["frozen_csv_positive_dd_min3_mismatches"] = diffs
        check["csv_check_status"] = "match" if not diffs else "MISMATCH"

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "in4mer_seed_sensitivity.py",
        "design": "repeat ONLY the control-pair sampling of in4mer_benchmark.py with "
                  "seeds 42-61; recompute the >=3-mutant/>=3-WT (min3) frame AUROC "
                  "(max |DD| across both orientations x all lineages; sklearn "
                  "roc_auc_score on |DD|, in4mer_gold=1). Frozen in4mer outputs "
                  "untouched; seed 42 reproduces them as a built-in check.",
        "frame": "min3 (>=3 mutant and >=3 WT lines per lineage stratum)",
        "n_controls_per_seed": N_CONTROLS,
        "seeds": SEEDS,
        "per_seed": per_seed,
        "summary": {
            "min": float(np.nanmin(aucs)),
            "median": float(np.nanmedian(aucs)),
            "max": float(np.nanmax(aucs)),
            "mean": float(np.nanmean(aucs)),
            "std": float(np.nanstd(aucs)),
            "n_seeds": len(SEEDS),
            "n_in4mer_evaluable_per_seed": sorted({r["n_in4mer_evaluable_of_13"]
                                                   for r in per_seed}),
        },
        "frozen_benchmark_consistency_check": check,
    }

    JSON_OUT.write_text(json.dumps(out, indent=2, allow_nan=False, default=str))
    print(f"\n  20-seed AUROC: min={aucs.min():.4f} median={np.median(aucs):.4f} "
          f"max={aucs.max():.4f}")
    print(f"  frozen check: {check.get('status')} / csv {check.get('csv_check_status')}")
    print(f"  Saved: {JSON_OUT}")


if __name__ == "__main__":
    main()
