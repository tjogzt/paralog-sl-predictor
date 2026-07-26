"""
Task 5 — TCGA BRCA survival analysis v2 (continuous, adjusted, PH-tested)
==========================================================================
Upgrade of tcga_survival.py (median-split univariate Cox, ARID1B HR=1.613
[1.163-2.238], p=0.004, n=1,069 / 151 events). lifelines is unavailable in
this runtime, so all Cox fits use statsmodels.duration.hazard_regression.PHReg
and the Schoenfeld-residual proportional-hazards test is implemented by hand.

Upgrades:
  1. Continuous expression (z-score) univariate Cox for all 32 paralog genes.
  2. Multivariable Cox adjusting for diagnosis age and AJCC pathologic stage
     (ordinal I-IV); stage missingness is reported, and an age-only model is
     also fitted so results remain interpretable if stage coverage is poor.
  3. PH assumption: correlation between Schoenfeld residuals and ranked event
     times (Grambsch-Therneau-style approximation, no variance scaling).
  4. BH-FDR across the 32 genes (univariate and multivariable families).

All data are fetched live from the cBioPortal API
(brca_tcga_pan_can_atlas_2018, rna_seq_v2_mrna expression + patient clinical);
the merged per-sample table is saved to output/tcga_survival_v2_samples.csv.

Outputs: output/tcga_survival_v2.json, output/tcga_survival_v2_samples.csv
"""

import json
import time

import numpy as np
import pandas as pd
import requests
from scipy import stats as sstats
from scipy.stats import false_discovery_control
from statsmodels.duration.hazard_regression import PHReg

from config import OUTPUT_DIR

BASE = "https://www.cbioportal.org/api"
STUDY = "brca_tcga_pan_can_atlas_2018"
SAMPLE_LIST = f"{STUDY}_all"

GENE_ENTREZ = {
    "ARID1A": 8289, "ARID1B": 57492, "SMARCA4": 6597, "SMARCA2": 6595,
    "BRCA1": 672, "BRCA2": 675, "EP300": 2033, "CREBBP": 1387,
    "PIK3CA": 5290, "PIK3CB": 5291, "PPP2R1A": 5518, "PPP2R1B": 5519,
    "FBXW7": 55294, "FBXW2": 26190, "STK11": 6794, "SIK1": 150094,
    "KRAS": 3845, "HRAS": 3265, "PIK3R1": 5295, "CRKL": 1399,
    "PTEN": 5728, "TNS2": 23371, "KMT2D": 8085, "KMT2C": 58508,
    "NF1": 4763, "RASA2": 5922, "ATR": 545, "ATM": 472,
    "RB1": 5925, "RBL1": 5933, "BRAF": 673, "RAF1": 5894,
}

STAGE_MAP = {"I": 1, "II": 2, "III": 3, "IV": 4}


def _req(method, url, retries=4, backoff=3.0, **kw):
    last = None
    for att in range(retries):
        try:
            r = requests.request(method, url, **kw)
            if r.status_code < 500:
                return r
            last = RuntimeError(f"HTTP {r.status_code}")
        except Exception as e:
            last = e
        time.sleep(backoff * (2 ** att))
    raise RuntimeError(f"request failed: {last}")


def parse_stage(val):
    """'STAGE IIIA'/'Stage IIIA' -> 3; 'STAGE X'/None/'' -> NaN."""
    if val is None:
        return np.nan
    s = str(val).upper().replace("STAGE", "").strip()
    if not s or s.startswith("X"):
        return np.nan
    for roman in ("IV", "III", "II", "I"):
        if s.startswith(roman):
            return float(STAGE_MAP[roman])
    return np.nan


def fetch_clinical():
    url = f"{BASE}/studies/{STUDY}/clinical-data"
    params = {"clinicalDataType": "PATIENT", "projection": "DETAILED", "pageSize": 100000}
    r = _req("GET", url, params=params, timeout=180)
    r.raise_for_status()
    clin = {}
    for item in r.json():
        pid = item.get("patientId", "")
        attr = item.get("clinicalAttributeId", "")
        val = item.get("value", "")
        if pid and attr and val != "":
            clin.setdefault(pid, {})[attr] = val
    return clin


def fetch_sample_patient():
    r = _req("GET", f"{BASE}/studies/{STUDY}/samples", params={"pageSize": 10000}, timeout=120)
    r.raise_for_status()
    return {s["sampleId"]: s["patientId"] for s in r.json()}


def fetch_expression_all():
    """One batched call for all 32 genes -> DataFrame samples x genes."""
    r = _req("GET", f"{BASE}/studies/{STUDY}/molecular-profiles", timeout=60)
    r.raise_for_status()
    prof = next(p["molecularProfileId"] for p in r.json()
                if "rna_seq_v2_mrna" in p["molecularProfileId"]
                and "zscores" not in p["molecularProfileId"].lower())
    url = f"{BASE}/molecular-profiles/{prof}/molecular-data/fetch"
    body = {"sampleListId": SAMPLE_LIST, "entrezGeneIds": list(GENE_ENTREZ.values())}
    r2 = _req("POST", url, json=body, timeout=300,
              headers={"Content-Type": "application/json"})
    r2.raise_for_status()
    rows = r2.json()
    df = pd.DataFrame(rows)
    entrez_to_gene = {v: k for k, v in GENE_ENTREZ.items()}
    df["gene"] = df["entrezGeneId"].map(entrez_to_gene)
    mat = df.pivot_table(index="sampleId", columns="gene", values="value",
                         aggfunc="first")
    return mat, prof


def cox_fit(times, exog, status):
    """PHReg fit -> params, bse, pvalues. Returns None on failure/NaN."""
    try:
        fit = PHReg(np.asarray(times, float), np.asarray(exog, float),
                    status=np.asarray(status, float)).fit(disp=0)
        if not np.all(np.isfinite(fit.params)) or not np.all(np.isfinite(fit.bse)):
            return None
        return fit.params, fit.bse, fit.pvalues
    except Exception:
        return None


def schoenfeld_ph_test(times, x, status, beta):
    """
    Manual Schoenfeld residuals for a single covariate (Breslow handling of
    ties): at each event time, r_i = x_i - weighted mean of x over the risk
    set. PH test = Pearson correlation of residuals with ranked event times
    (approximation of the Grambsch-Therneau test without variance scaling).
    Returns p-value (two-sided).
    """
    times = np.asarray(times, float)
    x = np.asarray(x, float)
    status = np.asarray(status, float)
    order = np.argsort(times)
    times, x, status = times[order], x[order], status[order]
    w = np.exp(x * beta)
    n = len(times)
    resid, evt_times = [], []
    for i in range(n):
        if status[i] != 1:
            continue
        risk = np.arange(n)[times >= times[i]]
        wr = w[risk]
        xbar = float(np.dot(wr, x[risk]) / wr.sum())
        resid.append(x[i] - xbar)
        evt_times.append(times[i])
    if len(resid) < 10 or np.std(resid) == 0:
        return np.nan
    ranks = sstats.rankdata(evt_times)
    _, p = sstats.pearsonr(resid, ranks)
    return float(p)


def hr_pack(coef, se, p):
    return {
        "hr": float(np.exp(coef)),
        "ci_low": float(np.exp(coef - 1.96 * se)),
        "ci_high": float(np.exp(coef + 1.96 * se)),
        "p": float(p),
    }


def main():
    print("=" * 70)
    print("  TCGA BRCA survival v2 (continuous + adjusted + PH test)")
    print("=" * 70)

    print("  Fetching clinical data...")
    clin = fetch_clinical()
    print(f"    {len(clin)} patients")
    s2p = fetch_sample_patient()
    print(f"    {len(s2p)} samples mapped")
    expr_mat, prof = fetch_expression_all()
    print(f"    expression: {expr_mat.shape[0]} samples x {expr_mat.shape[1]} genes ({prof})")

    # ── Build per-sample master table ──
    rows = []
    for sid, pid in s2p.items():
        c = clin.get(pid)
        if c is None:
            continue
        os_m, os_s = c.get("OS_MONTHS"), c.get("OS_STATUS")
        if os_m is None or os_s is None:
            continue
        try:
            os_m = float(os_m)
        except ValueError:
            continue
        if os_m <= 0:
            continue
        age = c.get("AGE")
        try:
            age = float(age)
        except (TypeError, ValueError):
            age = np.nan
        rows.append({
            "sample_id": sid, "patient_id": pid,
            "os_months": os_m,
            "event": 1 if str(os_s).startswith("1") else 0,
            "age": age,
            "stage_ord": parse_stage(c.get("AJCC_PATHOLOGIC_TUMOR_STAGE")),
        })
    master = pd.DataFrame(rows).set_index("sample_id")
    master = master.join(expr_mat, how="inner")
    master = master.dropna(subset=list(GENE_ENTREZ.keys()), how="all")
    n_stage_missing = int(master["stage_ord"].isna().sum())
    stage_missing_pct = 100.0 * n_stage_missing / len(master)
    print(f"    merged cohort: {len(master)} samples, {int(master['event'].sum())} events")
    print(f"    stage missing: {n_stage_missing} ({stage_missing_pct:.1f}%), "
          f"age missing: {int(master['age'].isna().sum())}")

    master.to_csv(OUTPUT_DIR / "tcga_survival_v2_samples.csv")

    results = []
    for gene in GENE_ENTREZ:
        sub = master.dropna(subset=[gene]).copy()
        if len(sub) < 100 or sub["event"].sum() < 10:
            continue
        # rna_seq_v2_mrna values are raw RSEM (highly right-skewed, e.g.
        # PIK3CA max 27,979): log2(x+1) before z-scoring, standard for
        # continuous-expression Cox models.
        log_expr = np.log2(np.clip(sub[gene].astype(float), 0, None) + 1.0)
        z = (log_expr - log_expr.mean()) / log_expr.std(ddof=0)
        sub["z"] = z

        # 1) univariate continuous
        uni = cox_fit(sub["os_months"], sub[["z"]], sub["event"])
        if uni is None:
            continue
        coef, se, p = uni[0][0], uni[1][0], uni[2][0]
        uni_pack = hr_pack(coef, se, p)

        # PH test on the univariate fit
        ph_p = schoenfeld_ph_test(sub["os_months"].values, z.values,
                                  sub["event"].values, coef)

        # 2) multivariable: z + age + stage (complete cases)
        mv = sub.dropna(subset=["age", "stage_ord"]).copy()
        mv_pack, mv_n, mv_events = None, int(len(mv)), int(mv["event"].sum())
        if len(mv) >= 100 and mv["event"].sum() >= 10:
            fit_mv = cox_fit(mv["os_months"], mv[["z", "age", "stage_ord"]], mv["event"])
            if fit_mv is not None:
                mv_pack = hr_pack(fit_mv[0][0], fit_mv[1][0], fit_mv[2][0])

        # 2b) age-only model (robustness if stage missingness is high)
        ma = sub.dropna(subset=["age"]).copy()
        ma_pack = None
        if len(ma) >= 100 and ma["event"].sum() >= 10:
            fit_ma = cox_fit(ma["os_months"], ma[["z", "age"]], ma["event"])
            if fit_ma is not None:
                ma_pack = hr_pack(fit_ma[0][0], fit_ma[1][0], fit_ma[2][0])

        results.append({
            "gene": gene,
            "n": int(len(sub)), "n_events": int(sub["event"].sum()),
            **{f"hr_continuous": uni_pack["hr"],
               "ci": [uni_pack["ci_low"], uni_pack["ci_high"]],
               "p": uni_pack["p"]},
            "ph_test_p": ph_p,
            "multivar_age_stage": ({"n": mv_n, "n_events": mv_events, **{
                "hr_multivar": mv_pack["hr"], "ci_multivar": [mv_pack["ci_low"], mv_pack["ci_high"]],
                "p_multivar": mv_pack["p"]}} if mv_pack else None),
            "multivar_age_only": ({"n": int(len(ma)), **{
                "hr_multivar": ma_pack["hr"], "ci_multivar": [ma_pack["ci_low"], ma_pack["ci_high"]],
                "p_multivar": ma_pack["p"]}} if ma_pack else None),
        })
        print(f"    {gene:9s} HR/SD={uni_pack['hr']:.3f} p={uni_pack['p']:.4f} "
              f"PH_p={ph_p if not np.isnan(ph_p) else float('nan'):.3f}")

    # ── BH-FDR across genes ──
    pv = np.array([r["p"] for r in results], dtype=float)
    valid = np.isfinite(pv)
    q = np.full(len(pv), np.nan)
    if valid.any():
        q[valid] = false_discovery_control(pv[valid])
    for r, qq in zip(results, q):
        r["q_fdr"] = float(qq) if np.isfinite(qq) else None
    mv_pv = np.array([r["multivar_age_stage"]["p_multivar"] if r["multivar_age_stage"] else np.nan
                      for r in results], dtype=float)
    mv_valid = np.isfinite(mv_pv)
    mv_q = np.full(len(mv_pv), np.nan)
    if mv_valid.any():
        mv_q[mv_valid] = false_discovery_control(mv_pv[mv_valid])
    for r, qq in zip(results, mv_q):
        if r["multivar_age_stage"]:
            r["multivar_age_stage"]["q_fdr_multivar"] = float(qq) if np.isfinite(qq) else None

    arid1b = next(r for r in results if r["gene"] == "ARID1B")

    # ── Sanity: reproduce the published median-split ARID1B result ──
    sub = master.dropna(subset=["ARID1B"]).copy()
    med = sub["ARID1B"].median()
    sub["high"] = (sub["ARID1B"] > med).astype(float)
    fit_ms = cox_fit(sub["os_months"], sub[["high"]], sub["event"])
    median_split_check = None
    if fit_ms is not None:
        ms = hr_pack(fit_ms[0][0], fit_ms[1][0], fit_ms[2][0])
        median_split_check = {
            "hr": ms["hr"], "ci": [ms["ci_low"], ms["ci_high"]], "p": ms["p"],
            "n": int(len(sub)), "n_events": int(sub["event"].sum()),
            "published": {"hr": 1.613, "ci": [1.163, 2.238], "p": 0.004,
                          "n": 1069, "n_events": 151},
        }
        print(f"\n  Median-split sanity check: HR={ms['hr']:.3f} "
              f"[{ms['ci_low']:.3f}-{ms['ci_high']:.3f}] p={ms['p']:.4f} "
              f"n={len(sub)}, events={int(sub['event'].sum())}")

    out_obj = {
        "study": STUDY, "expression_profile": prof,
        "cohort": {"n_samples": int(len(master)), "n_events": int(master["event"].sum()),
                   "stage_missing_pct": round(stage_missing_pct, 1),
                   "stage_note": ("stage coverage adequate; primary multivariable model "
                                  "adjusts for age + AJCC stage (complete cases)"
                                  if stage_missing_pct <= 30 else
                                  "stage missingness high (>30%); rely on age-only model"),
                   "age_missing": int(master["age"].isna().sum())},
        "models": {
            "univariate": "PHReg OS ~ z(log2(expression+1)); HR per +1 SD of "
                          "log2-transformed RSEM",
            "multivariable": "PHReg OS ~ z(log2 expression) + age + AJCC stage "
                             "(ordinal I-IV)",
            "ph_test": "Schoenfeld residuals (manual, Breslow ties) vs ranked event "
                       "times, Pearson correlation p; approximation without variance scaling",
            "fdr": "Benjamini-Hochberg across 32 genes (scipy false_discovery_control)",
        },
        "per_gene": results,
        "arid1b_highlight": arid1b,
        "arid1b_median_split_replication": median_split_check,
    }
    out = OUTPUT_DIR / "tcga_survival_v2.json"
    with open(out, "w") as fh:
        json.dump(out_obj, fh, indent=2)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
