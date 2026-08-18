#!/usr/bin/env python3
"""
rev2_b16_kmer_length_correction.py  (Stage-4 revision, item B16)
================================================================
The k-mer Jaccard identity proxy (k=3, compute_sequence_identity.R) is
length-inflated: for unrelated proteins the number of chance-shared 3-mers
grows with chain length, so long proteins (ARID1A/B ~2,300 aa, KMT2C/D
~4,900/5,500 aa, EP300/CREBBP ~2,400 aa) receive inflated "identity"
independent of homology. This script quantifies the inflation within the
26-gene local sequence cache (data/uniprot_sequences.rds; the only genes
with sequences, covering 16 of the 72 TableS2 pairs) and applies a
length correction:

  1. k-mer Jaccard (k=3, unique k-mers, R-identical implementation) for all
     325 unordered pairs of the 26 cached genes; fidelity check against
     output/paralog_identity.csv for the pairs scored there.
  2. Inflation: Pearson/Spearman correlation of Jaccard with pair length
     (min and geometric-mean length) over the 312 non-benchmark pairs
     (13 curated benchmark pairs excluded from the fit population).
  3. Correction: OLS of Jaccard on log(len_A)+log(len_B) fitted on the
     non-benchmark pairs; corrected identity = residual + intercept
     (length-neutralized, same center as raw for the reference population).
  4. Consequence for the manuscript's DD + identity filter
     (headline_metrics.json identity_filter): recompute the filter subsets
     and signed-DD AUROC with corrected identity at rank-equivalent
     thresholds (same number of pairs retained as id_ge_0.2 / id_ge_0.3),
     plus Spearman between raw and corrected identity on the 16 scored
     pairs.

No simulated data. Deterministic (OLS, no resampling).
Output: output/revision_stage4/rev2_b16_kmer_length_correction.{json,csv}
Usage: python rev2_b16_kmer_length_correction.py   (run from repo root)
"""

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

from config import OUTPUT_DIR, KNOWN_PARALOG_SL  # noqa: E402
from compute_headline_metrics import auroc  # noqa: E402

K = 3
FID_TOL = 1e-12


def kmers(seq, k=K):
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def jaccard(a, b):
    u = a | b
    return len(a & b) / len(u) if u else np.nan


def main():
    seqs = {}
    with open(OUT / "rev2_b16_sequences_dump.tsv") as fh:
        next(fh)
        for line in fh:
            gene, seq = line.rstrip("\n").split("\t")
            seqs[gene] = seq
    print(f"{len(seqs)} sequences loaded")
    lengths = {g: len(s) for g, s in seqs.items()}
    km = {g: kmers(s) for g, s in seqs.items()}

    # ── all-pairs Jaccard + fidelity vs paralog_identity.csv ────────
    ident = pd.read_csv(OUTPUT_DIR / "paralog_identity.csv")
    id_map = {(r.driver_gene, r.paralog_gene): r.kmer_jaccard
              for r in ident.itertuples() if pd.notna(r.kmer_jaccard)}
    genes = sorted(seqs)
    rows = []
    max_diff = 0.0
    for a, b in combinations(genes, 2):
        j = jaccard(km[a], km[b])
        rows.append({"gene_A": a, "gene_B": b, "jaccard": j,
                     "len_A": lengths[a], "len_B": lengths[b],
                     "min_len": min(lengths[a], lengths[b]),
                     "geomean_len": float(np.sqrt(lengths[a] * lengths[b]))})
        for key in ((a, b), (b, a)):
            if key in id_map:
                d = abs(j - id_map[key])
                max_diff = max(max_diff, d)
    if max_diff > FID_TOL:
        raise RuntimeError(f"FIDELITY FAIL: k-mer Jaccard differs from artifact "
                           f"by {max_diff}")
    print(f"fidelity: max |Jaccard recomputed - paralog_identity.csv| = {max_diff:.2e}")
    pairs = pd.DataFrame(rows)

    known_unordered = {tuple(sorted((a.upper(), b.upper()))) for a, b in KNOWN_PARALOG_SL}
    pairs["is_benchmark_pair"] = [
        tuple(sorted((r.gene_A.upper(), r.gene_B.upper()))) in known_unordered
        for r in pairs.itertuples()]
    ref = pairs[~pairs["is_benchmark_pair"]].copy()
    print(f"{len(pairs)} pairs total, {pairs['is_benchmark_pair'].sum()} benchmark, "
          f"{len(ref)} reference pairs for the length model")

    # ── inflation quantification ────────────────────────────────────
    infl = {}
    for var in ("min_len", "geomean_len"):
        pr, pp = pearsonr(ref["jaccard"], ref[var])
        sr, sp = spearmanr(ref["jaccard"], ref[var])
        infl[var] = {"pearson_r": pr, "pearson_p": pp,
                     "spearman_r": sr, "spearman_p": sp}
        print(f"inflation ({var}): Pearson r={pr:.3f} (p={pp:.2e}), "
              f"Spearman rho={sr:.3f}")

    # ── length correction: Jaccard ~ log(len_A) + log(len_B) ────────
    X = np.column_stack([np.ones(len(ref)),
                         np.log(ref["len_A"]), np.log(ref["len_B"])])
    beta, *_ = np.linalg.lstsq(X, ref["jaccard"], rcond=None)
    ref_pred = X @ beta
    r2 = 1 - float(((ref["jaccard"] - ref_pred) ** 2).sum()
                   / ((ref["jaccard"] - ref["jaccard"].mean()) ** 2).sum())
    Xall = np.column_stack([np.ones(len(pairs)),
                            np.log(pairs["len_A"]), np.log(pairs["len_B"])])
    pairs["expected_jaccard"] = Xall @ beta
    pairs["jaccard_corrected"] = pairs["jaccard"] - pairs["expected_jaccard"] + beta[0]
    print(f"length model: R^2 = {r2:.3f} on reference pairs; coefficients "
          f"intercept={beta[0]:.4f}, log_len_A={beta[1]:.4f}, log_len_B={beta[2]:.4f}")

    # ── consequence for the DD + identity filter ────────────────────
    tab = pd.read_csv(OUTPUT_DIR / "tables" / "TableS2_FullResults.tsv", sep="\t")
    corr_map = {}
    raw_map = {}
    for r in pairs.itertuples():
        for key in ((r.gene_A, r.gene_B), (r.gene_B, r.gene_A)):
            corr_map[key] = r.jaccard_corrected
            raw_map[key] = r.jaccard
    tab["_raw"] = [raw_map.get((d, p), np.nan)
                   for d, p in zip(tab["driver_gene"], tab["paralog_gene"])]
    tab["_corr"] = [corr_map.get((d, p), np.nan)
                    for d, p in zip(tab["driver_gene"], tab["paralog_gene"])]

    # verify raw join reproduces the artifact's 16 pairs-with-identity
    n_with = int(tab["_raw"].notna().sum())
    with open(OUTPUT_DIR / "headline_metrics.json") as fh:
        hm = json.load(fh)["identity_filter"]

    filt = {}
    scored = tab.dropna(subset=["_raw"])
    srho, _ = spearmanr(scored["_raw"], scored["_corr"])
    filt["spearman_raw_vs_corrected_on_scored_entries"] = srho
    for thr_name, thr in (("id_ge_0.2", 0.2), ("id_ge_0.3", 0.3)):
        art = hm[thr_name]
        # raw-threshold subset (fidelity vs artifact)
        sub_raw = tab[tab["_raw"] >= thr]
        auc_raw = auroc(sub_raw["is_known_paralog_sl"].astype(int).values,
                        sub_raw["dependency_dd"].fillna(0).values)
        # rank-equivalent corrected subset: same number of entries, top corrected
        sub_corr = tab.dropna(subset=["_corr"]).nlargest(len(sub_raw), "_corr")
        auc_corr = auroc(sub_corr["is_known_paralog_sl"].astype(int).values,
                         sub_corr["dependency_dd"].fillna(0).values)
        pairs_raw = sorted({f"{d}->{p}" for d, p in
                            zip(sub_raw["driver_gene"], sub_raw["paralog_gene"])})
        pairs_corr = sorted({f"{d}->{p}" for d, p in
                             zip(sub_corr["driver_gene"], sub_corr["paralog_gene"])})
        filt[thr_name] = {
            "raw_threshold": thr,
            "artifact": art,
            "raw_recomputed": {"auroc": auc_raw, "n_entries": len(sub_raw),
                               "n_positives": int(sub_raw["is_known_paralog_sl"].sum()),
                               "pairs": pairs_raw},
            "corrected_rank_equivalent": {"auroc": auc_corr,
                                          "n_entries": len(sub_corr),
                                          "n_positives": int(sub_corr["is_known_paralog_sl"].sum()),
                                          "pairs": pairs_corr},
        }
        print(f"{thr_name}: raw AUROC {auc_raw:.4f} (artifact {art['auroc']:.4f}) -> "
              f"corrected AUROC {auc_corr:.4f} (rank-equivalent n={len(sub_corr)})")
        print(f"   raw pairs: {pairs_raw}")
        print(f"   corr pairs: {pairs_corr}")

    pairs.to_csv(OUT / "rev2_b16_kmer_length_correction.csv", index=False)
    out = {
        "method": {
            "kmer": "Jaccard over unique 3-mers, R-identical (compute_sequence_identity.R)",
            "inflation_fit_population": f"{len(ref)} non-benchmark pairs of the 26 cached genes",
            "correction": "OLS Jaccard ~ log(len_A)+log(len_B) on non-benchmark pairs; "
                          "corrected = residual + intercept",
            "filter_consequence": "signed-DD AUROC recomputed on length-corrected identity "
                                  "at rank-equivalent subset sizes (same n entries as the "
                                  "raw id>=0.2 / id>=0.3 filters)",
            "scope_note": "sequence cache covers 26 genes = 16 of 72 TableS2 pairs; "
                          "external frames not computable (see B5)",
        },
        "fidelity_max_abs_diff_vs_paralog_identity_csv": max_diff,
        "inflation": infl,
        "length_model": {"intercept": float(beta[0]), "coef_log_len_A": float(beta[1]),
                         "coef_log_len_B": float(beta[2]), "r_squared": r2,
                         "n_reference_pairs": len(ref)},
        "n_tableS2_entries_with_identity": n_with,
        "filter_consequence": filt,
    }
    out_path = OUT / "rev2_b16_kmer_length_correction.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nwritten: {out_path}")


if __name__ == "__main__":
    main()
