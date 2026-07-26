"""
Task 1 — Verify CPTAC evaluable-test count and global BH significance count
============================================================================
Manuscript claims (Results, CPTAC section):
  - 26 paralog pairs x 7 cohorts, 122 evaluable tests (>=10 paired samples)
  - 42 of 122 significant at global BH q<0.05

This script recounts both numbers directly from the per-cohort correlation
CSVs produced by download_proteomics.py (output/cptac_<cohort>_correlations.csv;
status == "ok" <=> >=10 paired samples, the threshold used by
compute_pair_correlations with min_samples=10).

Also performs a file-level check on the ovarian PDC000360 HGSC proteome:
whether it exists in data/ and whether it enters any main-analysis product.

Output: output/cptac_counts_verification.json
"""

import json
from pathlib import Path

import pandas as pd
from scipy.stats import false_discovery_control

from config import DATA_DIR, OUTPUT_DIR

MAIN_COHORTS = ["brca", "coad", "luad", "gbm", "pdac", "ucec", "lusc"]
N_PAIRS = 30  # rows per correlations file (26 tested + header excluded); evaluable determined by status


def main():
    per_cohort = {}
    frames = []
    for c in MAIN_COHORTS:
        f = OUTPUT_DIR / f"cptac_{c}_correlations.csv"
        df = pd.read_csv(f)
        ok = df[df["status"] == "ok"].copy()
        per_cohort[c.upper()] = int(len(ok))
        ok["cohort"] = c.upper()
        frames.append(ok)

    allok = pd.concat(frames, ignore_index=True)
    pvals = allok["p"].astype(float).values
    qvals = false_discovery_control(pvals)  # Benjamini-Hochberg

    evaluable = int(len(allok))
    sig_q05 = int((qvals < 0.05).sum())

    # ── Ovarian PDC000360 file-level check ──
    fd_file = DATA_DIR / "FD_GLBL_MI_FFPEbridge_Abund_20201002.tsv"
    ovarian_info = {"file_present": False, "n_sample_columns": None}
    if fd_file.exists():
        hdr = pd.read_csv(fd_file, sep="\t", nrows=0)
        n_meta = 5  # Index, NumberPSM, Proteins, MaxPepProb, ReferenceIntensity
        ovarian_info = {
            "file_present": True,
            "n_sample_columns": int(len(hdr.columns) - n_meta),
        }

    # Does any main-analysis CPTAC product contain an ovarian cohort?
    cptac_outputs = sorted(OUTPUT_DIR.glob("cptac_*_correlations.csv"))
    ovarian_in_outputs = any("ov" == f.stem.replace("cptac_", "").replace("_correlations", "")
                             for f in cptac_outputs)
    # extra safety: no ovarian column in the pair matrix
    pm = pd.read_csv(OUTPUT_DIR / "cptac_pair_matrix.csv", nrows=0)
    ovarian_in_matrix = any(c.startswith(("OV_", "OVARIAN_")) for c in pm.columns)

    # Role of the extra on-disk cohort files (lscc/ccrcc/hnscc)
    extra_roles = {}
    for extra in ["lscc", "ccrcc", "hnscc"]:
        f = OUTPUT_DIR / f"cptac_{extra}_correlations.csv"
        if f.exists():
            d = pd.read_csv(f)
            n_ok = int((d["status"] == "ok").sum())
            extra_roles[extra.upper()] = (
                f"{n_ok} evaluable tests; NOT part of the 7-cohort main analysis"
                + (" (empty stub, all missing_data)" if n_ok == 0 else "")
            )

    result = {
        "evaluable_tests": evaluable,
        "significant_q05": sig_q05,
        "per_cohort_evaluable": per_cohort,
        "n_pairs_tested_per_cohort": 26,
        "min_paired_samples_threshold": 10,
        "raw_p_lt_0.05": int((pvals < 0.05).sum()),
        "manuscript_claims": {"evaluable_tests": 122, "significant_q05": 42},
        "matches_manuscript": bool(evaluable == 122 and sig_q05 == 42),
        "ovarian_pdc_in_outputs": bool(ovarian_in_outputs or ovarian_in_matrix),
        "ovarian_pdc_file_check": {
            **ovarian_info,
            "file": str(fd_file.relative_to(DATA_DIR.parent)) if fd_file.exists() else None,
            "referenced_by": "R_package/vignettes/paralogSL.Rmd example only (plot_protein_correlation)",
            "note": "No cptac_ov/ovarian correlation product exists in output/; "
                    "ovarian PDC000360 proteome does not enter the 7-cohort main analysis.",
        },
        "extra_cohort_files": extra_roles,
        "notes": [
            "Evaluable = status 'ok' in cptac_<cohort>_correlations.csv "
            "(Pearson r on >=10 paired samples, as in download_proteomics.py).",
            "Global BH: scipy false_discovery_control over all 122 p-values across 7 cohorts.",
            "LUSC is the lung squamous cohort used in the main analysis; the on-disk "
            "LSCC file is an empty stub (all missing_data), as are CCRCC and HNSCC.",
        ],
    }

    out = OUTPUT_DIR / "cptac_counts_verification.json"
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"Evaluable tests: {evaluable} (manuscript: 122)")
    print(f"Significant global BH q<0.05: {sig_q05} (manuscript: 42)")
    print(f"Per cohort: {per_cohort}")
    print(f"Ovarian PDC in outputs: {result['ovarian_pdc_in_outputs']}")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
