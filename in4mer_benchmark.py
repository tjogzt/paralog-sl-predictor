"""
in4mer Extended External Benchmark
===================================
The in4mer meta-analysis (Aregger et al., bioRxiv 2023.01.03.522655;
"Combined genome-scale fitness and paralog synthetic lethality screens
with just 44k clones") defined 13 cross-platform paralog synthetic-lethal
gold standards (paralog score >= 0.25, hit in >1 of five digenic-KO
screens; their Fig. 1J).

These pairs are MUTATION-AGNOSTIC gold standards (digenic knockout),
a different evidence class from our driver-conditioned Tier A/B pairs.
We therefore do NOT merge them into Table S3; instead we test whether
DD ranks them above unlabeled paralog pairs when either member's
natural mutation status is used as the conditioning event.

Design
------
* Positives: the 13 in4mer pairs. For each pair we try BOTH orientations
  (A driver/B paralog and B driver/A paralog); a pair's score is the
  max |DD| across orientations x lineages, mirroring the per-pair
  framework of the main benchmark.
* Controls: 400 unlabeled paralog pairs (seed 42) from the same
  Ensembl/HGNC pair universe, excluding the known-SL set and the in4mer
  pairs, evaluated identically (positive-unlabeled caveat disclosed in
  the manuscript — contamination biases the AUROC DOWN, i.e. the
  estimate is conservative).
* Frames: >=5 mutant/>=5 WT per lineage (primary, but only 3/13
  positives evaluable) and >=3/>=3 (sensitivity, 9/13 evaluable).
  Pairs with no evaluable lineage are excluded from that frame in BOTH
  arms.
* Statistics: AUROC with seeded bootstrap 95% CI (1,000 resamples) and
  a seeded label-permutation p-value (10,000 shuffles), matching the
  main validation suite.

Outputs (single source of truth for any manuscript numbers):
  output/in4mer_benchmark.csv          — per-pair scores + labels
  output/in4mer_feasibility.csv        — per-pair evaluability detail
  output/in4mer_benchmark_summary.json — AUROC, CI, p, n per frame
"""

import json

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

# ── in4mer 13 cross-platform paralog gold standards (their Fig. 1J) ──
IN4MER_PAIRS = [
    ("CNOT7", "CNOT8"), ("PITPNA", "PITPNB"), ("TIA1", "TIAL1"),
    ("SAR1A", "SAR1B"), ("PTP4A1", "PTP4A2"), ("GSK3A", "GSK3B"),
    ("CSNK2A1", "CSNK2A2"), ("CSNK1D", "CSNK1E"), ("MAPK1", "MAPK3"),
    ("ARFGEF1", "ARFGEF2"), ("HDAC1", "HDAC2"), ("ASF1A", "ASF1B"),
    ("SLC25A28", "SLC25A37"),
]

N_CONTROLS = 400
LINEAGE_COL = "OncotreePrimaryDisease"


def main():
    print("=" * 70)
    print("  in4mer Extended External Benchmark")
    print("=" * 70)

    dep = load_dependency()
    models = load_models()
    mut = load_mutations()
    paralogs = load_paralogs()

    cell_lines = list(dep.index)
    print(f"  Dependency matrix: {len(cell_lines)} lines")
    lin_map = models.set_index("DepMap_ID")[LINEAGE_COL].reindex(cell_lines)

    # ── Control pair sampling (seeded) ──
    exclude = {frozenset(p) for p in KNOWN_PARALOG_SL}
    exclude |= {frozenset({"MEK1", "MEK2"}), frozenset({"MAP2K1", "MAP2K2"})}
    exclude |= {frozenset(p) for p in IN4MER_PAIRS}
    pa = paralogs[["gene_A", "gene_B"]].dropna()
    ok = pa["gene_A"].isin(dep.columns) & pa["gene_B"].isin(dep.columns)
    pa = pa[ok]
    keep = [frozenset((r.gene_A, r.gene_B)) not in exclude for r in pa.itertuples()]
    pa = pa[keep]
    rng = np.random.default_rng(42)
    take = rng.choice(len(pa), size=min(N_CONTROLS, len(pa)), replace=False)
    control_pairs = [tuple(x) for x in pa.iloc[take].values]
    print(f"  Control pairs sampled: {len(control_pairs)} (seed 42)")

    # ── Mutation matrix for all genes involved ──
    all_pairs = IN4MER_PAIRS + control_pairs
    genes = sorted({g for pair in all_pairs for g in pair})
    mat = build_mutation_matrix(mut, cell_lines, genes, apply_driver_rules=True)

    lineages = sorted(set(lin_map.dropna()))
    lin_idx = {lin: lin_map[lin_map == lin].index for lin in lineages}

    def pair_score(a, b, min_n):
        """max |DD| across orientations x lineages; None if not evaluable."""
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

    # ── Score all pairs on both frames ──
    rows = []
    pos_set = {frozenset(p) for p in IN4MER_PAIRS}
    total = len(all_pairs)
    for i, (a, b) in enumerate(all_pairs):
        if (i + 1) % 50 == 0:
            print(f"    scoring {i + 1}/{total}...")
        dd5, n5 = pair_score(a, b, 5)
        dd3, n3 = pair_score(a, b, 3)
        rows.append({
            "pair": f"{a}/{b}", "gene_a": a, "gene_b": b,
            "label": "in4mer_gold" if frozenset((a, b)) in pos_set else "unlabeled_control",
            "dd_min5": dd5, "n_lineages_eval_5": n5,
            "dd_min3": dd3, "n_lineages_eval_3": n3,
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_DIR / "in4mer_benchmark.csv", index=False)

    # ── Metrics per frame ──
    summary = {"n_in4mer_total": len(IN4MER_PAIRS),
               "n_controls_sampled": len(control_pairs), "seed": 42}
    for frame, col, ncol in (("min5", "dd_min5", "n_lineages_eval_5"),
                             ("min3", "dd_min3", "n_lineages_eval_3")):
        sub = df.dropna(subset=[col]).copy()
        # in4mer gold standards are mutation-agnostic and carry no driver
        # direction, so the natural score here is |DD|; signed DD is reported
        # as a sensitivity value (arbitrary driver orientation makes it
        # uninformative by construction).
        sub["score"] = sub[col].abs()
        sub["score_signed"] = sub[col]
        y = (sub["label"] == "in4mer_gold").astype(int).values
        s = sub["score"].values
        s_signed = sub["score_signed"].values
        n_pos = int(y.sum())
        n_neg = int(len(y) - n_pos)
        entry = {"n_pos": n_pos, "n_neg": n_neg}
        if n_pos >= 2 and n_neg >= 2:
            auc = float(roc_auc_score(y, s))
            brng = np.random.default_rng(42)
            bs = []
            for _ in range(1000):
                bi = brng.integers(0, len(y), len(y))
                if y[bi].sum() >= 1 and y[bi].sum() < len(y):
                    bs.append(roc_auc_score(y[bi], s[bi]))
            null = np.array([roc_auc_score(brng.permutation(y), s)
                             for _ in range(10000)])
            entry.update({
                "auroc": auc,
                "auroc_signed_dd_sensitivity": float(roc_auc_score(y, s_signed)),
                "bootstrap_ci_low": float(np.percentile(bs, 2.5)),
                "bootstrap_ci_high": float(np.percentile(bs, 97.5)),
                "permutation_p": float((np.sum(null >= auc) + 1) / 10001),
                "median_dd_pos": float(sub.loc[y == 1, "score"].median()),
                "median_dd_neg": float(sub.loc[y == 0, "score"].median()),
            })
        else:
            entry["note"] = "too few evaluable pairs for AUROC"
        summary[frame] = entry
        print(f"\n  [{frame}] positives={n_pos}, controls={n_neg}")
        for k, v in entry.items():
            if isinstance(v, float):
                print(f"    {k}: {v:.4f}")

    with open(OUTPUT_DIR / "in4mer_benchmark_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Saved: {OUTPUT_DIR}/in4mer_benchmark.csv")
    print(f"  Saved: {OUTPUT_DIR}/in4mer_benchmark_summary.json")


if __name__ == "__main__":
    main()
