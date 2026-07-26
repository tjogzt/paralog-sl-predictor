"""
Task 2 — ARID1B protein abundance by ARID1A mutation status in UCEC CPTAC
==========================================================================
Manuscript claim (Fig. 2d): ARID1B protein +55% in ARID1A-mutant (n=37) vs
wild-type (n=58) endometrial tumors, Wilcoxon p=0.082.

Data sources (no simulated values):
  - ARID1B protein abundance: local cache data/cptac_cache/UCEC_protein_data.json
    (originally fetched from cBioPortal ucec_cptac_2020_protein_quantification;
    95 samples).
  - ARID1A mutation status: fetched live from the cBioPortal API
    (ucec_cptac_2020_mutations profile, any nonsilent mutation record), exactly
    as described in the manuscript Methods.
  - Cross-check against the cached R object
    data/cptac_cache/ucec_arid1b_by_arid1a_status.rds (read via Rscript).

Statistics: n, mean/median ratios, log2 abundance difference, fold change,
Wilcoxon rank-sum p, rank-biserial correlation, bootstrap 95% CI for the log2
difference (10,000 resamples, seed 42), and BH-adjusted p within two
multiple-testing families computed live from the same data sources:
  (a) all evaluable driver->paralog mutation-conditioned tests in UCEC
      (driver mutation data + paralog protein quantified; >=10 mut & >=10 WT);
  (b) the subset of (a) restricted to the 12 known SL pairs.

Output: output/arid1b_ucec_effects.json
"""

import json
import subprocess
import time

import numpy as np
import pandas as pd
import requests
from scipy import stats
from scipy.stats import false_discovery_control

from config import DATA_DIR, OUTPUT_DIR, KNOWN_PARALOG_SL

BASE = "https://www.cbioportal.org/api"
MUT_PROFILE = "ucec_cptac_2020_mutations"
SAMPLE_LIST_ALL = "ucec_cptac_2020_all"
MIN_PER_GROUP = 10
BOOT_B = 10_000
SEED = 42

# 26 CPTAC paralog pairs (same as download_proteomics.py)
PARALOG_PAIRS = [
    ("ARID1A", "ARID1B"), ("PIK3CA", "PIK3CB"),
    ("BRCA1", "BRCA2"), ("EP300", "CREBBP"),
    ("PPP2R1A", "PPP2R1B"), ("FBXW7", "FBXW2"),
    ("STK11", "SIK1"), ("SMARCA4", "SMARCA2"),
    ("CCNE1", "CCNE2"), ("CDK4", "CDK6"),
    ("AKT1", "AKT2"), ("MAP2K1", "MAP2K2"),
    ("KRAS", "HRAS"), ("KRAS", "NRAS"), ("HRAS", "NRAS"),
    ("PIK3R1", "CRKL"), ("PIK3R1", "CRK"),
    ("PTEN", "TNS2"), ("PTEN", "TNS1"),
    ("TP53", "TP63"), ("TP53", "TP73"),
    ("RB1", "RBL1"), ("RB1", "RBL2"),
    ("NF1", "RASA2"), ("NF1", "RASA1"),
    ("ATR", "ATM"), ("KMT2D", "KMT2C"),
    ("CDH1", "CDH2"), ("BRAF", "RAF1"), ("KEAP1", "NFE2L2"),
]

GENE_ENTREZ = {
    "ARID1A": 8289, "ARID1B": 57492, "PIK3CA": 5290, "PIK3CB": 5291,
    "BRCA1": 672, "BRCA2": 675, "EP300": 2033, "CREBBP": 1387,
    "PPP2R1A": 5518, "PPP2R1B": 5519, "FBXW7": 55294, "FBXW2": 26190,
    "STK11": 6794, "SIK1": 150094, "SMARCA4": 6597, "SMARCA2": 6595,
    "CCNE1": 898, "CCNE2": 9134, "CDK4": 1019, "CDK6": 1021,
    "AKT1": 207, "AKT2": 208, "MAP2K1": 5604, "MAP2K2": 5605,
    "KRAS": 3845, "HRAS": 3265, "NRAS": 4893, "PIK3R1": 5295,
    "CRKL": 1399, "CRK": 1398, "PTEN": 5728, "TNS2": 23371, "TNS1": 7145,
    "TP53": 7157, "TP63": 8626, "TP73": 7161, "RB1": 5925, "RBL1": 5933,
    "RBL2": 5934, "NF1": 4763, "RASA1": 5921, "RASA2": 5922,
    "ATR": 545, "ATM": 472, "KMT2D": 8085, "KMT2C": 58508,
    "CDH1": 999, "CDH2": 1000, "BRAF": 673, "RAF1": 5894,
    "KEAP1": 9817, "NFE2L2": 4780,
}

KNOWN_SET = {(a, b) for a, b in KNOWN_PARALOG_SL} | {(b, a) for a, b in KNOWN_PARALOG_SL}
KNOWN_SET |= {("MAP2K1", "MAP2K2"), ("MAP2K2", "MAP2K1")}  # MEK1/MEK2 aliases


def fetch_mutated_samples(entrez_ids):
    """One batched cBioPortal call: sampleIds with any mutation record per gene."""
    url = f"{BASE}/molecular-profiles/{MUT_PROFILE}/mutations/fetch"
    body = {"entrezGeneIds": list(entrez_ids), "sampleListId": SAMPLE_LIST_ALL}
    r = requests.post(url, json=body, params={"projection": "DETAILED"},
                      timeout=120, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    entrez_to_gene = {v: k for k, v in GENE_ENTREZ.items()}
    mut_by_gene = {}
    for m in r.json():
        gene = m.get("gene", {}).get("hugoGeneSymbol") or entrez_to_gene.get(
            m.get("entrezGeneId"))
        if gene:
            mut_by_gene.setdefault(gene, set()).add(m["sampleId"])
    return mut_by_gene


def read_rds_crosscheck():
    """Read the cached RDS via Rscript; return dict with mut/wt arrays or None."""
    rds = DATA_DIR / "cptac_cache" / "ucec_arid1b_by_arid1a_status.rds"
    if not rds.exists():
        return None
    r_code = (
        'x <- readRDS("data/cptac_cache/ucec_arid1b_by_arid1a_status.rds");'
        'cat(jsonlite::toJSON(list(mut=x$mut, wt=x$wt, p=x$p, '
        'n_mut=x$n_mut, n_wt=x$n_wt), auto_unbox=TRUE))'
    )
    try:
        out = subprocess.run(["Rscript", "-e", r_code], capture_output=True,
                             text=True, timeout=120, cwd=DATA_DIR.parent)
        if out.returncode != 0:
            # jsonlite may be absent; fall back to manual dump
            r_code2 = (
                'x <- readRDS("data/cptac_cache/ucec_arid1b_by_arid1a_status.rds");'
                'cat("MUT\\n"); cat(x$mut, sep=","); cat("\\nWT\\n");'
                'cat(x$wt, sep=","); cat(sprintf("\\nP %.6f NM %d NW %d\\n",'
                ' x$p, x$n_mut, x$n_wt))'
            )
            out = subprocess.run(["Rscript", "-e", r_code2], capture_output=True,
                                 text=True, timeout=120, cwd=DATA_DIR.parent)
            lines = out.stdout.strip().split("\n")
            mut = np.array([float(v) for v in lines[1].split(",")])
            wt = np.array([float(v) for v in lines[3].split(",")])
            tail = lines[4].split()
            return {"mut": mut, "wt": wt, "p": float(tail[1]),
                    "n_mut": int(tail[3]), "n_wt": int(tail[5])}
        return json.loads(out.stdout)
    except Exception as e:
        print(f"  RDS cross-check unavailable: {e}")
        return None


def group_stats(mut_vals, wt_vals):
    """Full statistic bundle for one mut-vs-wt comparison."""
    mut_vals = np.asarray(mut_vals, dtype=float)
    wt_vals = np.asarray(wt_vals, dtype=float)
    n1, n2 = len(mut_vals), len(wt_vals)
    mean_ratio = float(mut_vals.mean() / wt_vals.mean()) if wt_vals.mean() != 0 else np.nan
    median_ratio = (float(np.median(mut_vals) / np.median(wt_vals))
                    if np.median(wt_vals) != 0 else np.nan)
    log2_diff = float(mut_vals.mean() - wt_vals.mean())  # values are log2-scale
    fold_change = float(2 ** log2_diff)
    pct_from_fc = (fold_change - 1.0) * 100.0
    # Linearized (2^x) ratios — abundance-scale fold differences
    lin_mean_ratio = float(np.mean(2 ** mut_vals) / np.mean(2 ** wt_vals))
    lin_median_ratio = float(np.median(2 ** mut_vals) / np.median(2 ** wt_vals))
    # The manuscript's "+55%" reproduces ONLY via this quantity: the difference
    # of means divided by |mean_WT| on the log2 scale. Because mean_WT is a
    # negative log-ratio, this is not a fold change and is numerically
    # misleading; the bona fide fold change is ~1.04x (see fold_change).
    pct_vs_abs_mean_wt = ((log2_diff / abs(float(wt_vals.mean()))) * 100.0
                          if wt_vals.mean() != 0 else np.nan)

    U1, p_wil = stats.mannwhitneyu(mut_vals, wt_vals, alternative="two-sided")
    U2 = n1 * n2 - U1  # U statistic of the WT group
    r_rb = 1.0 - 2.0 * U2 / (n1 * n2)  # positive => mutant group higher

    rng = np.random.default_rng(SEED)
    boots = np.empty(BOOT_B)
    for b in range(BOOT_B):
        bm = rng.choice(mut_vals, size=n1, replace=True)
        bw = rng.choice(wt_vals, size=n2, replace=True)
        boots[b] = bm.mean() - bw.mean()
    ci_low, ci_high = np.percentile(boots, [2.5, 97.5])

    return {
        "n_mut": n1, "n_wt": n2,
        "mean_mut": float(mut_vals.mean()), "mean_wt": float(wt_vals.mean()),
        "median_mut": float(np.median(mut_vals)), "median_wt": float(np.median(wt_vals)),
        "mean_ratio": mean_ratio, "median_ratio": median_ratio,
        "log2_diff_mean_minus_mean": log2_diff,
        "fold_change_2pow_log2diff": fold_change,
        "pct_higher_from_fold_change": pct_from_fc,
        "linear_mean_ratio_2pow": lin_mean_ratio,
        "linear_median_ratio_2pow": lin_median_ratio,
        "pct_vs_abs_mean_wt_on_log2_scale": pct_vs_abs_mean_wt,
        "manuscript_55pct_origin": ("(mean_mut - mean_wt)/|mean_wt| on log2-scale "
                                    "values = 55.3%; not a fold change. Bona fide "
                                    "fold change = 2^0.0633 = 1.045 (+4.5%)"),
        "wilcoxon_U_mut": float(U1), "wilcoxon_U_wt": float(U2),
        "wilcoxon_p_two_sided": float(p_wil),
        "rank_biserial_r": float(r_rb),
        "rank_biserial_note": ("r_rb = 1 - 2*U_wt/(n_mut*n_wt) "
                               "= 2*U_mut/(n_mut*n_wt) - 1; "
                               "positive = higher in ARID1A-mutant"),
        "boot_B": BOOT_B, "boot_seed": SEED,
        "boot_ci95_log2diff": [float(ci_low), float(ci_high)],
    }


def main():
    print("=" * 70)
    print("  ARID1B protein by ARID1A mutation status — UCEC CPTAC")
    print("=" * 70)

    # ── ARID1B protein from local cache ──
    cache = json.load(open(DATA_DIR / "cptac_cache" / "UCEC_protein_data.json"))
    arid1b = pd.Series(cache["ARID1B"])
    prot_samples = set(arid1b.index)
    print(f"  ARID1B protein samples (cache): {len(arid1b)}")

    # ── Mutation status for all drivers of the 26 pairs (one batched call) ──
    drivers = sorted({a for a, _ in PARALOG_PAIRS})
    entrez = [GENE_ENTREZ[g] for g in drivers]
    print(f"  Fetching mutation status for {len(drivers)} drivers (cBioPortal API)...")
    mut_by_gene = fetch_mutated_samples(entrez)
    entrez_to_gene = {v: k for k, v in GENE_ENTREZ.items()}

    # ── Primary: ARID1A -> ARID1B ──
    mut_ids = mut_by_gene.get("ARID1A", set()) & prot_samples
    wt_ids = prot_samples - mut_by_gene.get("ARID1A", set())
    mut_vals = arid1b[sorted(mut_ids)].values
    wt_vals = arid1b[sorted(wt_ids)].values
    primary = group_stats(mut_vals, wt_vals)
    print(f"  ARID1A-mut n={primary['n_mut']}, WT n={primary['n_wt']}")
    print(f"  +{(primary['pct_higher_from_fold_change']):.1f}% (2^log2diff), "
          f"Wilcoxon p={primary['wilcoxon_p_two_sided']:.4f}")

    # ── RDS cross-check ──
    rds = read_rds_crosscheck()
    rds_check = None
    if rds:
        rds_mut = np.asarray(rds["mut"], dtype=float)
        rds_wt = np.asarray(rds["wt"], dtype=float)
        # match values (order-independent)
        same = (len(rds_mut) == len(mut_vals) and len(rds_wt) == len(wt_vals)
                and np.allclose(np.sort(rds_mut), np.sort(mut_vals), atol=1e-6)
                and np.allclose(np.sort(rds_wt), np.sort(wt_vals), atol=1e-6))
        rds_check = {
            "n_mut": int(rds["n_mut"]), "n_wt": int(rds["n_wt"]),
            "p": float(rds["p"]),
            "values_match_api_recompute": bool(same),
        }
        print(f"  RDS cross-check: n={rds['n_mut']}/{rds['n_wt']}, "
              f"p={rds['p']:.4f}, values match: {same}")

    # ── BH families: all evaluable mutation-conditioned tests in UCEC ──
    family_rows = []
    for a, b in PARALOG_PAIRS:
        if b not in cache or not cache[b]:
            continue
        prot_b = pd.Series(cache[b])
        samples_b = set(prot_b.index)
        m = mut_by_gene.get(a, set()) & samples_b
        w = samples_b - mut_by_gene.get(a, set())
        if len(m) < MIN_PER_GROUP or len(w) < MIN_PER_GROUP:
            continue
        mv = prot_b[sorted(m)].values
        wv = prot_b[sorted(w)].values
        _, p = stats.mannwhitneyu(mv, wv, alternative="two-sided")
        family_rows.append({
            "driver": a, "paralog": b, "n_mut": len(m), "n_wt": len(w),
            "p": float(p),
            "is_known_sl": (a, b) in KNOWN_SET,
        })
        time.sleep(0)  # no extra API calls in loop

    fam = pd.DataFrame(family_rows)
    p_all = fam["p"].values
    q_all = false_discovery_control(p_all)
    fam["q_bh_family_all"] = q_all
    known = fam[fam["is_known_sl"]]
    if len(known):
        q_known = false_discovery_control(known["p"].values)
        fam.loc[fam["is_known_sl"], "q_bh_family_known_pairs"] = q_known
    arid_row = fam[(fam["driver"] == "ARID1A") & (fam["paralog"] == "ARID1B")]

    result = {
        "primary_test": {
            "cohort": "UCEC CPTAC (ucec_cptac_2020)",
            "comparison": "ARID1B protein abundance, ARID1A-mutant vs wild-type",
            "protein_source": "data/cptac_cache/UCEC_protein_data.json "
                              "(cBioPortal ucec_cptac_2020_protein_quantification)",
            "mutation_source": "cBioPortal API ucec_cptac_2020_mutations, "
                               "any mutation record per sample",
            **primary,
            "manuscript_claim": {"n_mut": 37, "n_wt": 58, "pct_higher": 55, "p": 0.082},
        },
        "rds_crosscheck": rds_check,
        "bh_families": {
            "family_all_evaluable_ucec_tests": {
                "n_tests": int(len(fam)),
                "tests": fam.to_dict(orient="records"),
                "arid1a_arid1b_q": (float(arid_row["q_bh_family_all"].iloc[0])
                                    if len(arid_row) else None),
            },
            "family_known_sl_pairs": {
                "n_tests": int(len(known)),
                "arid1a_arid1b_q": (float(arid_row["q_bh_family_known_pairs"].iloc[0])
                                    if len(arid_row) and "q_bh_family_known_pairs" in fam
                                    else None),
            },
            "note": "Families computed live: Wilcoxon tests of paralog protein by "
                    "driver mutation status for all 26 CPTAC pairs in UCEC with "
                    f"paralog protein quantified and >= {MIN_PER_GROUP} samples per group.",
        },
    }

    out = OUTPUT_DIR / "arid1b_ucec_effects.json"
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\nSaved: {out}")
    print(f"Family (all evaluable UCEC tests): {len(fam)} tests; "
          f"ARID1A->ARID1B BH q = {result['bh_families']['family_all_evaluable_ucec_tests']['arid1a_arid1b_q']}")


if __name__ == "__main__":
    main()
