"""
Task 5 extension — TCGA UCEC + OV survival analysis (continuous, adjusted,
PH-tested) with ARID1A-mutation-stratified ARID1B analysis in UCEC
====================================================================
Reviewer-requested extension of tcga_survival_v2.py (BRCA cohort) to two
additional TCGA PanCan Atlas gynecologic cohorts:

  * ucec_tcga_pan_can_atlas_2018  (uterine corpus endometrial carcinoma)
  * ov_tcga_pan_can_atlas_2018    (ovarian serous cystadenocarcinoma)

Methodology is IDENTICAL to tcga_survival_v2.py — the core functions
(_req, parse_stage, cox_fit, schoenfeld_ph_test, hr_pack, GENE_ENTREZ)
are imported from that module so there is a single source of truth:
  1. Continuous expression: log2(RSEM+1) z-scored, univariate PHReg for
     all 32 paralog genes (HR per +1 SD).
  2. Multivariable PHReg adjusting for diagnosis age and AJCC pathologic
     stage (ordinal I-IV, complete cases); age-only robustness model too.
  3. Manual Schoenfeld-residual PH test (Breslow ties, Pearson correlation
     of residuals with ranked event times — Grambsch-Therneau approximation
     without variance scaling).
  4. BH-FDR (scipy false_discovery_control) across the 32 genes, separately
     for the univariate and multivariable families, within each cohort.

Additionally, in UCEC an ARID1A-mutation-stratified analysis is performed:
ARID1A mutations are fetched from the mutations profile; mutant = any
non-silent ARID1A mutation (missense / nonsense / frameshift / splice /
in-frame indel / nonstop, counted by type). The ARID1B continuous
expression Cox model (univariate; multivariable when stratum events allow)
is fitted SEPARATELY within ARID1A-mutant and within wild-type patients.
Strata with <10 events are labelled underpowered/descriptive via a
`power_note` field rather than presented as tests.

All data are fetched live from the cBioPortal API. No simulated data.

Outputs:
  output/tcga_survival_ucec.json,  output/tcga_survival_ucec_samples.csv
  output/tcga_survival_ov.json,    output/tcga_survival_ov_samples.csv
  output/tcga_survival_ucec_ov_summary.json
"""

import json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from config import OUTPUT_DIR
from tcga_survival_v2 import (
    BASE,
    GENE_ENTREZ,
    _req,
    cox_fit,
    hr_pack,
    parse_stage,
    schoenfeld_ph_test,
)

STUDIES = {
    "UCEC": "ucec_tcga_pan_can_atlas_2018",
    "OV": "ov_tcga_pan_can_atlas_2018",
}

ARID1A_ENTREZ = GENE_ENTREZ["ARID1A"]  # 8289

# cBioPortal mutationType values considered protein-altering ("non-silent").
# "Silent" and non-coding-region types (Intron, UTRs, IGR, Flanks, RNA)
# are excluded from the mutant definition.
NON_SILENT_TYPES = {
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Splice_Site", "Translation_Start_Site",
}

MIN_EVENTS_POWERED = 10  # strata below this are underpowered/descriptive


def fetch_clinical(study):
    url = f"{BASE}/studies/{study}/clinical-data"
    params = {"clinicalDataType": "PATIENT", "projection": "DETAILED",
              "pageSize": 100000}
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


def fetch_sample_patient(study):
    r = _req("GET", f"{BASE}/studies/{study}/samples",
             params={"pageSize": 10000}, timeout=120)
    r.raise_for_status()
    return {s["sampleId"]: s["patientId"] for s in r.json()}


def fetch_expression_all(study):
    """One batched call for all 32 genes -> DataFrame samples x genes."""
    r = _req("GET", f"{BASE}/studies/{study}/molecular-profiles", timeout=60)
    r.raise_for_status()
    prof = next(p["molecularProfileId"] for p in r.json()
                if "rna_seq_v2_mrna" in p["molecularProfileId"]
                and "zscores" not in p["molecularProfileId"].lower())
    url = f"{BASE}/molecular-profiles/{prof}/molecular-data/fetch"
    body = {"sampleListId": f"{study}_all",
            "entrezGeneIds": list(GENE_ENTREZ.values())}
    r2 = _req("POST", url, json=body, timeout=300,
              headers={"Content-Type": "application/json"})
    r2.raise_for_status()
    df = pd.DataFrame(r2.json())
    entrez_to_gene = {v: k for k, v in GENE_ENTREZ.items()}
    df["gene"] = df["entrezGeneId"].map(entrez_to_gene)
    mat = df.pivot_table(index="sampleId", columns="gene", values="value",
                         aggfunc="first")
    return mat, prof


def fetch_arid1a_mutations(study):
    """All ARID1A mutation records for the study -> list of dicts."""
    profile = f"{study}_mutations"
    url = f"{BASE}/molecular-profiles/{profile}/mutations/fetch"
    body = {"sampleListId": f"{study}_all", "entrezGeneIds": [ARID1A_ENTREZ]}
    r = _req("POST", url, json=body, timeout=300,
             headers={"Content-Type": "application/json"})
    r.raise_for_status()
    return r.json()


def build_master(study, clin, s2p, expr_mat):
    """Per-sample master table, identical schema to tcga_survival_v2."""
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
    return master


def run_gene_panel(master):
    """32-gene continuous Cox panel — exact tcga_survival_v2 methodology."""
    results = []
    genes_missing = [g for g in GENE_ENTREZ if g not in master.columns]
    for gene in GENE_ENTREZ:
        if gene not in master.columns:
            continue
        sub = master.dropna(subset=[gene]).copy()
        if len(sub) < 100 or sub["event"].sum() < 10:
            continue
        # rna_seq_v2_mrna values are raw RSEM (highly right-skewed):
        # log2(x+1) before z-scoring, standard for continuous-expression Cox.
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
            fit_mv = cox_fit(mv["os_months"], mv[["z", "age", "stage_ord"]],
                             mv["event"])
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
                "hr_multivar": mv_pack["hr"],
                "ci_multivar": [mv_pack["ci_low"], mv_pack["ci_high"]],
                "p_multivar": mv_pack["p"]}} if mv_pack else None),
            "multivar_age_only": ({"n": int(len(ma)), **{
                "hr_multivar": ma_pack["hr"],
                "ci_multivar": [ma_pack["ci_low"], ma_pack["ci_high"]],
                "p_multivar": ma_pack["p"]}} if ma_pack else None),
        })
        print(f"    {gene:9s} HR/SD={uni_pack['hr']:.3f} p={uni_pack['p']:.4f} "
              f"PH_p={ph_p if not np.isnan(ph_p) else float('nan'):.3f}")

    # ── BH-FDR across genes (same as v2: per family, within cohort) ──
    pv = np.array([r["p"] for r in results], dtype=float)
    valid = np.isfinite(pv)
    q = np.full(len(pv), np.nan)
    if valid.any():
        q[valid] = false_discovery_control(pv[valid])
    for r, qq in zip(results, q):
        r["q_fdr"] = float(qq) if np.isfinite(qq) else None
    mv_pv = np.array([r["multivar_age_stage"]["p_multivar"]
                      if r["multivar_age_stage"] else np.nan
                      for r in results], dtype=float)
    mv_valid = np.isfinite(mv_pv)
    mv_q = np.full(len(mv_pv), np.nan)
    if mv_valid.any():
        mv_q[mv_valid] = false_discovery_control(mv_pv[mv_valid])
    for r, qq in zip(results, mv_q):
        if r["multivar_age_stage"]:
            r["multivar_age_stage"]["q_fdr_multivar"] = \
                float(qq) if np.isfinite(qq) else None
    # age-only family FDR (becomes the primary adjusted family when stage
    # is unavailable in a cohort, per the v2 fallback rule)
    ma_pv = np.array([r["multivar_age_only"]["p_multivar"]
                      if r["multivar_age_only"] else np.nan
                      for r in results], dtype=float)
    ma_valid = np.isfinite(ma_pv)
    ma_q = np.full(len(ma_pv), np.nan)
    if ma_valid.any():
        ma_q[ma_valid] = false_discovery_control(ma_pv[ma_valid])
    for r, qq in zip(results, ma_q):
        if r["multivar_age_only"]:
            r["multivar_age_only"]["q_fdr_multivar"] = \
                float(qq) if np.isfinite(qq) else None
    return results, genes_missing


def cox_stratum(sub, label):
    """ARID1B z-expression Cox within one ARID1A stratum.

    Univariate always attempted (if fit numerically succeeds); multivariable
    (age+stage) attempted when the stratum has >= MIN_EVENTS_POWERED events.
    Strata with < MIN_EVENTS_POWERED events get a `power_note` and their HRs
    must be read as underpowered/descriptive, not as tests.
    """
    sub = sub.dropna(subset=["ARID1B"]).copy()
    n, events = int(len(sub)), int(sub["event"].sum())
    out = {"stratum": label, "n": n, "n_events": events}
    if events < MIN_EVENTS_POWERED:
        out["power_note"] = (
            f"UNDERPOWERED/DESCRIPTIVE: only {events} events "
            f"(< {MIN_EVENTS_POWERED}); HR estimate reported for "
            "descriptive purposes only, not as a hypothesis test.")
    if n < 15 or events < 3:
        out["fit_note"] = "stratum too small for a stable Cox fit; not fitted"
        return out

    log_expr = np.log2(np.clip(sub["ARID1B"].astype(float), 0, None) + 1.0)
    z = (log_expr - log_expr.mean()) / log_expr.std(ddof=0)
    sub["z"] = z

    uni = cox_fit(sub["os_months"], sub[["z"]], sub["event"])
    if uni is not None:
        pk = hr_pack(uni[0][0], uni[1][0], uni[2][0])
        out["univariate"] = {"hr": pk["hr"],
                             "ci": [pk["ci_low"], pk["ci_high"]],
                             "p": pk["p"]}
        out["ph_test_p"] = schoenfeld_ph_test(
            sub["os_months"].values, z.values, sub["event"].values, uni[0][0])
    else:
        out["univariate"] = None
        out["fit_note"] = "univariate PHReg fit failed to converge"

    # multivariable only if events allow (>= MIN_EVENTS_POWERED events).
    # Covariates mirror the v2 primary model (age + AJCC stage); when stage
    # is entirely unavailable in the cohort (as for UCEC/OV in cBioPortal)
    # fall back to the v2 age-only robustness model.
    if events >= MIN_EVENTS_POWERED:
        stage_available = sub["stage_ord"].notna().sum() > 0
        covars = ["z", "age", "stage_ord"] if stage_available else ["z", "age"]
        mv = sub.dropna(subset=covars[1:]).copy()
        if len(mv) >= 20 and mv["event"].sum() >= MIN_EVENTS_POWERED:
            fit_mv = cox_fit(mv["os_months"], mv[covars], mv["event"])
            if fit_mv is not None:
                pk = hr_pack(fit_mv[0][0], fit_mv[1][0], fit_mv[2][0])
                out["multivariable"] = {
                    "covariates": covars,
                    "n": int(len(mv)), "n_events": int(mv["event"].sum()),
                    "hr": pk["hr"], "ci": [pk["ci_low"], pk["ci_high"]],
                    "p": pk["p"]}
                if int(mv["event"].sum()) < 10 * len(covars):
                    out["multivariable"]["power_note"] = (
                        f"only {int(mv['event'].sum())} events for a "
                        f"{len(covars)}-covariate model (<10 events/variable); "
                        "interpret with caution")
            else:
                out["multivariable"] = None
        else:
            out["multivariable"] = None
            out["mv_note"] = "multivariable not fitted: insufficient complete cases"
    else:
        out["multivariable"] = None
        out["mv_note"] = (f"multivariable not fitted: {events} events "
                          f"< {MIN_EVENTS_POWERED}")
    return out


def stratified_arid1a(master, study):
    """ARID1A-mutation-stratified ARID1B Cox analysis (UCEC)."""
    muts = fetch_arid1a_mutations(study)
    type_counts = {}
    mut_samples = set()
    for m in muts:
        mtype = m.get("mutationType", "Unknown")
        type_counts[mtype] = type_counts.get(mtype, 0) + 1
        if mtype in NON_SILENT_TYPES:
            mut_samples.add(m.get("sampleId"))
    # patient-level mutation status within the analysis cohort
    cohort = master.copy()
    cohort["arid1a_mutant"] = [
        sid in mut_samples for sid in cohort.index
    ]
    n_mut = int(cohort["arid1a_mutant"].sum())
    print(f"    ARID1A mutation records fetched: {len(muts)} "
          f"(types: {type_counts})")
    print(f"    cohort patients with non-silent ARID1A mutation: {n_mut} "
          f"/ {len(cohort)}")

    res = {
        "mutation_definition": ("any non-silent ARID1A mutation in profile "
                                f"{study}_mutations; types counted: "
                                f"{sorted(NON_SILENT_TYPES)}; Silent and "
                                "non-coding-region records excluded"),
        "mutation_type_counts_all_records": type_counts,
        "n_cohort": int(len(cohort)),
        "arid1b_model": ("PHReg OS ~ z(log2(ARID1B RSEM+1)) within each "
                         "ARID1A stratum; multivariable adds age + AJCC "
                         "stage when stratum events allow (age-only fallback "
                         "when stage is unpopulated, as in this study)"),
        "strata": [
            cox_stratum(cohort[cohort["arid1a_mutant"]],
                        "ARID1A_mutant"),
            cox_stratum(cohort[~cohort["arid1a_mutant"]],
                        "ARID1A_wildtype"),
        ],
    }
    for s in res["strata"]:
        uni = s.get("univariate")
        line = (f"    {s['stratum']:17s} n={s['n']:4d} events={s['n_events']:3d}")
        if uni:
            line += (f"  ARID1B HR/SD={uni['hr']:.3f} "
                     f"[{uni['ci'][0]:.3f}-{uni['ci'][1]:.3f}] p={uni['p']:.4f}")
        else:
            line += "  (no univariate fit)"
        if s.get("power_note"):
            line += "  [UNDERPOWERED]"
        print(line)
    return res


def run_cohort(label, study):
    print("=" * 70)
    print(f"  TCGA {label} survival ({study})")
    print("=" * 70)
    print("  Fetching clinical data...")
    clin = fetch_clinical(study)
    print(f"    {len(clin)} patients")
    s2p = fetch_sample_patient(study)
    print(f"    {len(s2p)} samples mapped")
    expr_mat, prof = fetch_expression_all(study)
    print(f"    expression: {expr_mat.shape[0]} samples x "
          f"{expr_mat.shape[1]} genes ({prof})")

    master = build_master(study, clin, s2p, expr_mat)
    n_stage_missing = int(master["stage_ord"].isna().sum())
    stage_missing_pct = 100.0 * n_stage_missing / len(master)
    print(f"    merged cohort: {len(master)} samples, "
          f"{int(master['event'].sum())} events")
    print(f"    stage missing: {n_stage_missing} ({stage_missing_pct:.1f}%), "
          f"age missing: {int(master['age'].isna().sum())}")

    master.to_csv(OUTPUT_DIR / f"tcga_survival_{label.lower()}_samples.csv")

    results, genes_missing = run_gene_panel(master)
    arid1b = next((r for r in results if r["gene"] == "ARID1B"), None)

    out_obj = {
        "study": study, "expression_profile": prof,
        "methodology": "identical to tcga_survival_v2.py (BRCA cohort); "
                       "core functions imported from that module",
        "cohort": {"n_samples": int(len(master)),
                   "n_events": int(master["event"].sum()),
                   "stage_missing_pct": round(stage_missing_pct, 1),
                   "stage_note": ("stage coverage adequate; primary "
                                  "multivariable model adjusts for age + "
                                  "AJCC stage (complete cases)"
                                  if stage_missing_pct <= 30 else
                                  "stage missingness high (>30%); rely on "
                                  "age-only model"),
                   "age_missing": int(master["age"].isna().sum())},
        "genes_missing_from_profile": genes_missing,
        "models": {
            "univariate": "PHReg OS ~ z(log2(expression+1)); HR per +1 SD of "
                          "log2-transformed RSEM",
            "multivariable": "PHReg OS ~ z(log2 expression) + age + AJCC stage "
                             "(ordinal I-IV)",
            "ph_test": "Schoenfeld residuals (manual, Breslow ties) vs ranked "
                       "event times, Pearson correlation p; approximation "
                       "without variance scaling",
            "fdr": "Benjamini-Hochberg across 32 genes "
                   "(scipy false_discovery_control), separately per family "
                   "and per cohort",
        },
        "per_gene": results,
        "arid1b_highlight": arid1b,
    }

    if label == "UCEC":
        print("  ARID1A-mutation-stratified ARID1B analysis...")
        out_obj["arid1a_stratified"] = stratified_arid1a(master, study)

    out = OUTPUT_DIR / f"tcga_survival_{label.lower()}.json"
    with open(out, "w") as fh:
        json.dump(out_obj, fh, indent=2)
    print(f"\nSaved: {out}")
    return out_obj


def n_pass_fdr(results, key):
    return sum(1 for r in results
               if r.get(key) is not None and r[key] < 0.05)


def main():
    np.random.seed(42)
    print("=" * 70)
    print("  TCGA UCEC + OV survival (continuous + adjusted + PH test)")
    print("=" * 70)

    cohort_objs = {}
    for label, study in STUDIES.items():
        cohort_objs[label] = run_cohort(label, study)

    # ── machine-readable summary ──
    summary = {
        "generated_by": "tcga_survival_ucec_ov.py",
        "methodology": "identical to tcga_survival_v2.py (BRCA): continuous "
                       "z(log2(RSEM+1)) expression, statsmodels PHReg, "
                       "univariable + age/AJCC-stage-adjusted models, BH-FDR "
                       "per family per cohort",
        "ph_test_method": ("manual Schoenfeld-type test: Schoenfeld residuals "
                           "(Breslow ties) vs ranked event times, Pearson "
                           "correlation p-value; Grambsch-Therneau-style "
                           "approximation without variance scaling "
                           "(lifelines unavailable)"),
        "cohorts": {},
    }
    for label, obj in cohort_objs.items():
        hb = obj["arid1b_highlight"] or {}
        mv = hb.get("multivar_age_stage") or {}
        entry = {
            "study": obj["study"],
            "n": obj["cohort"]["n_samples"],
            "events": obj["cohort"]["n_events"],
            "stage_missing_pct": obj["cohort"]["stage_missing_pct"],
            "age_missing": obj["cohort"]["age_missing"],
            "genes_missing_from_profile": obj["genes_missing_from_profile"],
            "arid1b": {
                "univariable": {
                    "hr": hb.get("hr_continuous"),
                    "ci": hb.get("ci"),
                    "p": hb.get("p"),
                    "q_fdr": hb.get("q_fdr"),
                },
                "multivariable_age_stage": {
                    "hr": mv.get("hr_multivar"),
                    "ci": mv.get("ci_multivar"),
                    "p": mv.get("p_multivar"),
                    "q_fdr": mv.get("q_fdr_multivar"),
                    "n": mv.get("n"),
                    "n_events": mv.get("n_events"),
                } if mv else None,
                "multivariable_age_only": ({
                    "hr": hb["multivar_age_only"].get("hr_multivar"),
                    "ci": hb["multivar_age_only"].get("ci_multivar"),
                    "p": hb["multivar_age_only"].get("p_multivar"),
                    "q_fdr": hb["multivar_age_only"].get("q_fdr_multivar"),
                    "n": hb["multivar_age_only"].get("n"),
                } if hb.get("multivar_age_only") else None),
                "adjusted_model_note": (
                    "AJCC stage entirely unpopulated for this study in "
                    "cBioPortal (verified: attribute defined, zero values at "
                    "patient and sample level); the age+stage model has zero "
                    "complete cases, so the v2 age-only fallback model is the "
                    "adjusted estimate of record"
                    if obj["cohort"]["stage_missing_pct"] > 30 else
                    "adjusted estimate adjusts for age + AJCC stage "
                    "(complete cases)"),
                "ph_test_p": hb.get("ph_test_p"),
            },
            "n_genes_fdr_pass_univariable": n_pass_fdr(obj["per_gene"], "q_fdr"),
            "n_genes_fdr_pass_multivariable": sum(
                1 for r in obj["per_gene"]
                if r.get("multivar_age_stage")
                and r["multivar_age_stage"].get("q_fdr_multivar") is not None
                and r["multivar_age_stage"]["q_fdr_multivar"] < 0.05),
            "n_genes_fdr_pass_age_only": sum(
                1 for r in obj["per_gene"]
                if r.get("multivar_age_only")
                and r["multivar_age_only"].get("q_fdr_multivar") is not None
                and r["multivar_age_only"]["q_fdr_multivar"] < 0.05),
        }
        if "arid1a_stratified" in obj:
            entry["arid1a_stratified"] = obj["arid1a_stratified"]
        summary["cohorts"][label] = entry

    sout = OUTPUT_DIR / "tcga_survival_ucec_ov_summary.json"
    with open(sout, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nSaved: {sout}")

    # ── console digest ──
    print("\n" + "=" * 70)
    print("  DIGEST")
    print("=" * 70)
    for label, e in summary["cohorts"].items():
        a = e["arid1b"]
        print(f"  {label}: n={e['n']}, events={e['events']}, "
              f"stage missing {e['stage_missing_pct']}%")
        if a["univariable"]["hr"] is not None:
            u = a["univariable"]
            print(f"    ARID1B uni : HR={u['hr']:.3f} "
                  f"[{u['ci'][0]:.3f}-{u['ci'][1]:.3f}] p={u['p']:.4f} "
                  f"q={u['q_fdr']:.4f}" if u["q_fdr"] is not None else
                  f"    ARID1B uni : HR={u['hr']:.3f} p={u['p']:.4f}")
        if a["multivariable_age_stage"]:
            m = a["multivariable_age_stage"]
            qstr = (f" q={m['q_fdr']:.4f}" if m["q_fdr"] is not None else "")
            print(f"    ARID1B mult(age+stage): HR={m['hr']:.3f} "
                  f"[{m['ci'][0]:.3f}-{m['ci'][1]:.3f}] p={m['p']:.4f}{qstr} "
                  f"(n={m['n']}, events={m['n_events']})")
        if a["multivariable_age_only"]:
            m = a["multivariable_age_only"]
            qstr = (f" q={m['q_fdr']:.4f}" if m["q_fdr"] is not None else "")
            print(f"    ARID1B mult(age-only) : HR={m['hr']:.3f} "
                  f"[{m['ci'][0]:.3f}-{m['ci'][1]:.3f}] p={m['p']:.4f}{qstr} "
                  f"(n={m['n']})")
        print(f"    FDR<0.05 genes: univariable="
              f"{e['n_genes_fdr_pass_univariable']}, multivariable(age+stage)="
              f"{e['n_genes_fdr_pass_multivariable']}, age-only="
              f"{e['n_genes_fdr_pass_age_only']}")


if __name__ == "__main__":
    main()
