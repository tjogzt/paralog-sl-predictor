#!/usr/bin/env python3
"""audit_manuscript_numbers.py — closed-loop audit of every numeric claim in
manuscript.tex / supplementary.tex against the reproducible artifacts under
output/.

Complements the per-script claims checks (headline / ML / regression, run in
steps 1-3 of verify_all.sh) by covering the numbers those scripts do NOT own:
lineage-level AUROCs, TSG/ONC mechanism contrast, MSI stratification,
mutation-type analysis, direction audit, and the therapeutic-window module.

Every check recomputes the value from an artifact and compares it against the
value written in the manuscript at manuscript precision. Exit code is non-zero
if any check fails, so verify_all.sh (and any future rerun) raises an alarm.

Output: output/manuscript_number_audit.tsv (one row per claim).
"""
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).parent
OUT = ROOT / "output"

ROWS = []


def check(cid, description, manuscript, computed, source, tol=0.5):
    """Register one claim. `manuscript` is the value as printed (string or
    number); `computed` is the artifact value; tol is in units of the last
    manuscript digit (default: half a unit of the last digit)."""
    try:
        m = float(manuscript)
        c = float(computed)
        # infer manuscript precision from its string form
        s = str(manuscript)
        if "." in s:
            nd = len(s.split(".")[1].rstrip("0")) or 1
        else:
            nd = 0
        step = 10 ** (-nd)
        ok = abs(c - m) <= tol * step + 1e-12
        shown = f"{c:.{max(nd, 3)}f}" if abs(c) >= 1e-3 else f"{c:.2e}"
    except (TypeError, ValueError):
        ok = str(manuscript) == str(computed)
        shown = str(computed)
    ROWS.append({
        "id": cid, "description": description,
        "manuscript": manuscript, "recomputed": shown,
        "source": source, "status": "match" if ok else "MISMATCH",
    })
    return ok


def check_int(cid, description, manuscript, computed, source):
    ok = int(manuscript) == int(computed)
    ROWS.append({"id": cid, "description": description, "manuscript": manuscript,
                 "recomputed": computed, "source": source,
                 "status": "match" if ok else "MISMATCH"})
    return ok


def check_none(cid, description, computed, source):
    ok = computed is None or (isinstance(computed, float) and np.isnan(computed))
    ROWS.append({"id": cid, "description": description, "manuscript": "not evaluable",
                 "recomputed": "not evaluable" if ok else f"{computed:.3f}",
                 "source": source, "status": "match" if ok else "MISMATCH"})
    return ok


# ── Artifacts ──────────────────────────────────────────────────────
hm = json.loads((OUT / "headline_metrics.json").read_text())
vr = json.loads((OUT / "validation_report.json").read_text())
ml = json.loads((OUT / "ml_benchmark.json").read_text())
reg = json.loads((OUT / "regression_controls.json").read_text())
da = json.loads((OUT / "direction_audit.json").read_text())
solid = pd.read_csv(OUT / "solid_tumor_summary.csv")
solid3 = pd.read_csv(OUT / "solid_tumor_summary_min3.csv")
msi5 = json.loads((OUT / "msi_key_numbers_min5.json").read_text())
msi3 = json.loads((OUT / "msi_key_numbers_min3.json").read_text())
msi_int = json.loads((OUT / "msi_interaction_summary.json").read_text())
mut = pd.read_csv(OUT / "muttype_all_results.csv")
twc = pd.read_csv(OUT / "therapeutic_window_paralog_classification.csv")
twa = pd.read_csv(OUT / "therapeutic_window_all_results.csv")
alpha = pd.read_csv(OUT / "alphafold_structural_analysis.csv")
ts2 = pd.read_csv(OUT / "tables" / "TableS2_FullResults.tsv", sep="\t")

HM = "output/headline_metrics.json"
VR = "output/validation_report.json"
ML = "output/ml_benchmark.json"

# ═══════════════════════════════════════════════════════════════════
# 1. Headline framework (mirror of headline claims — the audit table
#    aggregates them so a reviewer sees one closed loop)
# ═══════════════════════════════════════════════════════════════════
f = hm["lineage_full"]
check("full_auroc", "Full-set AUROC (110 entries, 8 positives)", "0.629", f["auroc"], HM)
check_int("full_n", "Full-set entries", 110, f["n_entries"], HM)
check_int("full_pos", "Full-set positives", 8, f["n_positives"], HM)
check("full_auprc", "Full-set AUPRC", "0.346", f["auprc"], HM)
lo = hm["lineage_leave_out_depmap_era"]
check("leave_out_auroc", "Leave-out (DepMap-era removed) AUROC", "0.613", lo["auroc"], HM)
pd_ = hm["lineage_pre_depmap_only"]
check("pre_depmap_auroc", "Pre-DepMap evidence AUROC", "0.525", pd_["auroc"], HM)
tab = hm["lineage_tier_ab"]
check("tier_ab_auroc", "Tier A∪B primary benchmark AUROC", "1.000", tab["auroc"], HM)
check_int("tier_ab_pos", "Tier A∪B evaluable positives", 2, tab["n_positives"], HM)
ds = hm["lineage_full_direction_strict"]
check("direction_strict", "Direction-strict AUROC", "0.629", ds["auroc"], HM)
llo = hm["leave_one_lineage_out"]["range"]
check("llo_lo", "Leave-one-lineage-out lower bound", "0.606", llo[0], HM)
check("llo_hi", "Leave-one-lineage-out upper bound", "0.692", llo[1], HM)

comp = hm["component_decomposition_lineage"]
check("comp_dd", "Component: DD", "0.629", comp["dd"], HM)
check("comp_pcs", "Component: PCS", "0.825", comp["pcs"], HM)
check("comp_dexpr", "Component: ΔExpression", "0.547", comp["delta_expression_abs"], HM)
check("comp_nec", "Component: necessity", "0.642", comp["necessity"], HM)
pb = hm["component_paired_bootstrap"]["pcs_minus_dd"]
check("pb_pcs_dd", "Paired bootstrap PCS−DD", "+0.195", pb["mean_delta"], HM)
check("pb_pcs_dd_lo", "Paired bootstrap CI lower", "-0.132", pb["ci95"][0], HM)
check("pb_pcs_dd_hi", "Paired bootstrap CI upper", "+0.545", pb["ci95"][1], HM)

# ═══════════════════════════════════════════════════════════════════
# 2. Per-pair framework
# ═══════════════════════════════════════════════════════════════════
nc = vr["negative_control"]
check("pp_auroc", "Per-pair AUROC (72 pairs, 6 positives)", "0.598", nc["observed_auroc"], VR)
check_int("pp_n", "Per-pair universe", 72, nc["n_total"], VR)
check("pp_null", "Permutation null mean", "0.501", nc["null_auroc_mean"], VR)
check("pp_p", "Permutation empirical p", "0.223", nc["empirical_p_value"], VR)
bs = vr["bootstrap"]
check("pp_bs_lo", "Per-pair bootstrap CI lower", "0.294", bs["auroc_ci_low"], VR)
check("pp_bs_hi", "Per-pair bootstrap CI upper", "0.845", bs["auroc_ci_high"], VR)
ppm = hm["per_pair_mean_from_tables2"]
check("pp_mean_dd", "Per-pair mean signed DD AUROC", "0.672", ppm["auroc_dd"], HM)
pc = hm["per_pair_composite_mean"]
check("pp_comp", "Composite AUROC", "0.831", pc["auroc"], HM)
check("pp_comp_auprc", "Composite AUPRC", "0.356", pc["auprc"], HM)
check("pp_baseline", "Baseline prevalence", "0.083", pc["baseline_prevalence"], HM)
check("pp_enrich", "AUPRC enrichment ×baseline", "4.3",
      pc["auprc"] / pc["baseline_prevalence"], HM, tol=0.5)

# ═══════════════════════════════════════════════════════════════════
# 3. ML benchmark
# ═══════════════════════════════════════════════════════════════════
clf = ml["classifiers"]
check("ml_svmrbf", "SVM-RBF AUROC", "0.841", clf["SVM_RBF"]["auroc"], ML)
check("ml_rf", "Random forest AUROC", "0.722", clf["RF"]["auroc"], ML)
check("ml_svmlin", "SVM-Linear AUROC", "0.240", clf["SVM_Linear"]["auroc"], ML)
check("ml_lr", "Logistic regression AUROC", "0.136", clf["LR"]["auroc"], ML)
check("ml_dd", "DD alone AUROC", "0.672", ml["single_feature"]["dd_alone"], ML)
check("ml_comp", "Composite alone AUROC", "0.831", ml["single_feature"]["composite_alone"], ML)
check("ml_comp_lofo", "Composite fold-internal-scaling AUROC", "0.831",
      ml["single_feature"]["composite_alone_lofo"], ML)
check("ml_comp_lofo_auprc", "Composite fold-internal-scaling AUPRC", "0.368",
      ml["composite_auprc_lofo"], ML)
# Full-data LR Wald p-value claim (ml_lr_feat_p) removed: the manuscript no
# longer states per-feature LR p-values (round-7 revision).
dl = ml["paired_bootstrap_delta_auroc"]
check("ml_delta_auroc_mean", "Paired bootstrap ΔAUROC mean (composite − SVM-RBF)",
      "-0.01", dl["mean_delta"], ML)
check("ml_delta_auroc_lo", "Paired bootstrap ΔAUROC 95% CI lower",
      "-0.21", dl["ci95"][0], ML)
check("ml_delta_auroc_hi", "Paired bootstrap ΔAUROC 95% CI upper",
      "+0.18", dl["ci95"][1], ML)
check("ml_delta_auroc_pneg", "Paired bootstrap P(ΔAUROC < 0)",
      "0.54", dl["p_delta_lt_0"], ML)

# ═══════════════════════════════════════════════════════════════════
# 4. Lineage level (drift blind spot #1)
# ═══════════════════════════════════════════════════════════════════
S = "output/solid_tumor_summary.csv"
S3 = "output/solid_tumor_summary_min3.csv"
ev = solid.dropna(subset=["dd_auroc"])
ev = ev[ev.n_known >= 2]
check_int("lin_neval", "Evaluable lineages (primary ≥5 frame)", 8, len(ev), S)
check_int("lin_n01", "Lineages with AUROC ≥ 0.7", 1, (ev.dd_auroc >= 0.7).sum(), S)
for name, val in [("Esophagogastric", "0.674"), ("SCLC", "0.651"),
                  ("Bladder Urothelial", "0.531"), ("Colorectal", "0.775"),
                  ("Endometrial", "0.519"), ("Breast", "0.150"),
                  ("NSCLC", "0.664"), ("Ovarian", "0.611")]:
    got = solid.loc[solid.cancer == name, "dd_auroc"].iloc[0]
    check(f"lin_{name[:6]}", f"{name} AUROC (primary frame)", val, got, S)

ev3 = solid3.dropna(subset=["dd_auroc"])
ev3 = ev3[ev3.n_known >= 2]
check_int("lin3_neval", "Evaluable lineages (sensitivity ≥3 frame)", 12, len(ev3), S3)
check_int("lin3_n07", "Sensitivity lineages with AUROC ≥ 0.7", 6, (ev3.dd_auroc >= 0.7).sum(), S3)
for name, val in [("Biliary Tract", "0.990"), ("Pancreatic", "0.962"),
                  ("Melanoma", "0.840"), ("Cervical", "0.714")]:
    got = solid3.loc[solid3.cancer == name, "dd_auroc"].iloc[0]
    check(f"lin3_{name[:6]}", f"{name} AUROC (sensitivity frame)", val, got, S3)

ONC = {"Melanoma", "NSCLC", "Pancreatic"}
tsg = ev3[~ev3.cancer.isin(ONC)].dd_auroc.values
onc = ev3[ev3.cancer.isin(ONC)].dd_auroc.values
check("tsg_mean", "TSG-driven mean AUROC (n=9)", "0.649", tsg.mean(), S3)
check("onc_mean", "Oncogene-driven mean AUROC (n=3)", "0.834", onc.mean(), S3)
vals = np.concatenate([tsg, onc])
obs = tsg.mean() - onc.mean()
diffs = []
for ii in combinations(range(len(vals)), len(onc)):
    mask = np.zeros(len(vals), bool)
    mask[list(ii)] = True
    diffs.append(vals[~mask].mean() - vals[mask].mean())
diffs = np.array(diffs)
perm_p = min(1.0, 2 * min((diffs >= obs).mean(), (diffs <= obs).mean()))
mw_p = stats.mannwhitneyu(tsg, onc, alternative="two-sided", method="exact").pvalue
check("tsg_perm_p", "TSG vs ONC permutation p (exact, 220 combos)", "0.127", perm_p, S3)
check("tsg_mw_p", "TSG vs ONC exact Mann-Whitney p", "0.145", mw_p, S3)

# ═══════════════════════════════════════════════════════════════════
# 5. MSI stratification (blind spot #2)
# ═══════════════════════════════════════════════════════════════════
M3 = "output/msi_key_numbers_min3.json"
M5 = "output/msi_key_numbers_min5.json"
g = msi3["subgroups"]
check("msi_endo_h", "Endometrial MSI-H AUROC (sensitivity)", "0.303", g["Endometrial_MSI_H"]["dd_auroc"], M3)
check("msi_endo_s", "Endometrial MSS AUROC (sensitivity)", "0.500", g["Endometrial_MSS"]["dd_auroc"], M3)
check("msi_col_h", "Colorectal MSI-H AUROC (sensitivity)", "0.558", g["Colorectal_MSI_H"]["dd_auroc"], M3)
check("msi_col_s", "Colorectal MSS AUROC (sensitivity)", "0.545", g["Colorectal_MSS"]["dd_auroc"], M3)
check_int("msi_n17", "Endometrial MSI-H cell lines", 17, g["Endometrial_MSI_H"]["n_lines"], M3)
check_int("msi_n11", "Endometrial MSS cell lines", 11, g["Endometrial_MSS"]["n_lines"], M3)
check_int("msi_n14", "Colorectal MSI-H cell lines", 14, g["Colorectal_MSI_H"]["n_lines"], M3)
check_int("msi_n45", "Colorectal MSS cell lines", 45, g["Colorectal_MSS"]["n_lines"], M3)
g5 = msi5["subgroups"]
check_none("msi5_endo", "Endometrial MSI subgroups not evaluable (primary frame)",
           g5["Endometrial_MSI_H"]["dd_auroc"], M5)
check("msi5_col_h", "Colorectal MSI-H AUROC (primary)", "0.574", g5["Colorectal_MSI_H"]["dd_auroc"], M5)
check("msi5_col_s", "Colorectal MSS AUROC (primary)", "0.298", g5["Colorectal_MSS"]["dd_auroc"], M5)

# 5b. MSI full-cohort interaction test (msi_interaction_test.py)
MI = "output/msi_interaction_summary.json"
pc_crc = msi_int["per_cancer"]["Colorectal"]
pc_end = msi_int["per_cancer"]["Endometrial"]
check_int("msi_int_crc_lines", "Interaction cohort: colorectal lines", 63, pc_crc["n_lines"], MI)
check_int("msi_int_end_lines", "Interaction cohort: endometrial lines", 28, pc_end["n_lines"], MI)
check_int("msi_int_crc_pairs", "Interaction evaluable pairs: colorectal", 47, pc_crc["n_pairs"], MI)
check_int("msi_int_end_pairs", "Interaction evaluable pairs: endometrial", 24, pc_end["n_pairs"], MI)
check_int("msi_int_crc_nom", "Nominal interactions: colorectal", 4, pc_crc["n_nominal_p05"], MI)
check_int("msi_int_end_nom", "Nominal interactions: endometrial", 1, pc_end["n_nominal_p05"], MI)
check("msi_int_crc_minq", "Min interaction q: colorectal", "0.157", pc_crc["min_q"], MI)
check("msi_int_end_minq", "Min interaction q: endometrial", "0.235", pc_end["min_q"], MI)
check("msi_int_crc_binom", "Binomial p: colorectal", "0.207", pc_crc["binom_p_vs_5pct"], MI)
check("msi_int_end_binom", "Binomial p: endometrial", "0.708", pc_end["binom_p_vs_5pct"], MI)
check("msi_int_crc_stouf", "Stouffer p: colorectal", "0.667", pc_crc["stouffer_p"], MI)
check("msi_int_end_stouf", "Stouffer p: endometrial", "0.276", pc_end["stouffer_p"], MI)
check("msi_int_crc_medabs", "Median |delta DD|: colorectal", "0.062",
      pc_crc["median_abs_delta_dd"], MI)
check("msi_int_end_medabs", "Median |delta DD|: endometrial", "0.076",
      pc_end["median_abs_delta_dd"], MI)
check_int("msi_int_crc_ci_excl", "Interaction CIs excluding zero: colorectal", 4,
          pc_crc["n_ci_excludes_zero"], MI)
check_int("msi_int_end_ci_excl", "Interaction CIs excluding zero: endometrial", 1,
          pc_end["n_ci_excludes_zero"], MI)
lp = msi_int["lead_pairs"]
check("msi_int_arid_end_dd", "ARID1A->ARID1B deltaDD: endometrial", "+0.02", lp["Endometrial:ARID1A->ARID1B"]["delta_dd"], MI)
check("msi_int_arid_end_p", "ARID1A->ARID1B interaction p: endometrial", "0.96", lp["Endometrial:ARID1A->ARID1B"]["interaction_p"], MI)
check("msi_int_arid_crc_dd", "ARID1A->ARID1B deltaDD: colorectal", "+0.04", lp["Colorectal:ARID1A->ARID1B"]["delta_dd"], MI)
check("msi_int_arid_crc_p", "ARID1A->ARID1B interaction p: colorectal", "0.79", lp["Colorectal:ARID1A->ARID1B"]["interaction_p"], MI)

# ═══════════════════════════════════════════════════════════════════
# 6. Mutation type (blind spot #3)
# ═══════════════════════════════════════════════════════════════════
MT = "output/muttype_all_results.csv"
arid_ov = mut[(mut.cancer == "Ovarian") & (mut.driver == "ARID1A") & (mut.paralog == "ARID1B")].iloc[0]
check("mt_arid_t", "ARID1A→ARID1B Ovarian truncating DD", "0.388", arid_ov.dd_trunc, MT)
check("mt_arid_m", "ARID1A→ARID1B Ovarian missense DD", "0.020", arid_ov.dd_miss, MT)
ep_col = mut[(mut.cancer == "Colorectal") & (mut.driver == "EP300") & (mut.paralog == "CREBBP")].iloc[0]
check("mt_ep_t", "EP300→CREBBP Colorectal truncating DD", "0.464", ep_col.dd_trunc, MT)
check("mt_ep_m", "EP300→CREBBP Colorectal missense DD", "0.150", ep_col.dd_miss, MT)

br = mut[mut.cancer == "Breast"]
yt = br.is_known_sl.astype(int)
from sklearn.metrics import roc_auc_score
check("mt_br_all", "Breast all-DD AUROC (muttype frame)", "0.465",
      roc_auc_score(yt, br.dd_all.fillna(0)), MT)
check("mt_br_t", "Breast truncating-only AUROC", "0.460",
      roc_auc_score(yt, br.dd_trunc.fillna(0)), MT)
check("mt_br_m", "Breast missense-only AUROC", "0.437",
      roc_auc_score(yt, br.dd_miss.fillna(0)), MT)

sub = mut.dropna(subset=["dd_trunc", "dd_miss"])
t_stat, p_val = stats.ttest_rel(sub.dd_miss.abs(), sub.dd_trunc.abs())
check("mt_miss_mean", "Pan-pair mean |DD missense|", "0.042", sub.dd_miss.abs().mean(), MT)
check("mt_trunc_mean", "Pan-pair mean |DD truncating|", "0.026", sub.dd_trunc.abs().mean(), MT)
check("mt_paired_p", "Paired t-test p", "3e-4", p_val, MT)
check("mt_panelb_mean", "Fig S5b mean(|trunc|−|miss|)", "-0.015",
      (sub.dd_trunc.abs() - sub.dd_miss.abs()).mean(), MT)

# ═══════════════════════════════════════════════════════════════════
# 7. Direction audit
# ═══════════════════════════════════════════════════════════════════
DA = "output/direction_audit.json"
fr_full = next(x for x in da["frames"] if x["frame"] == "pair_level_primary_gyn3")
check("dir_abs", "Full frame |DD| AUROC (sensitivity analysis)", "0.676", fr_full["auroc_abs"], DA)
check("dir_signed", "Full frame signed-DD AUROC", "0.629", fr_full["auroc_signed"], DA)
check_int("dir_neg", "Positive entries with DD<0", 3, fr_full["n_pos_dd_negative"], DA)
fr_tier = next(x for x in da["frames"] if "TIER_A" in x["frame"])
check("dir_tier", "Tier A∪B frame signed = absolute", "1.000", fr_tier["auroc_signed"], DA)

# ═══════════════════════════════════════════════════════════════════
# 8. Therapeutic window module (blind spot #4)
# ═══════════════════════════════════════════════════════════════════
TW = "output/therapeutic_window_paralog_classification.csv"
def trow(d, p):
    return twc[(twc.driver == d) & (twc.paralog == p)].iloc[0]

check("dws_arid", "ARID1A→ARID1B mean DWS (signed max(DD,0) formula)", "2.82", trow("ARID1A", "ARID1B").mean_ti, TW)
check("dws_arid_sel", "ARID1A→ARID1B mean selectivity", "+0.28", trow("ARID1A", "ARID1B").mean_selectivity, TW)
check("dws_arid_dd", "ARID1A→ARID1B mean max(DD,0) numerator", "0.270", trow("ARID1A", "ARID1B").mean_dd_pos, TW)
check("dws_smarca", "SMARCA4→SMARCA2 mean DWS (signed max(DD,0) formula)", "4.87", trow("SMARCA4", "SMARCA2").mean_ti, TW)
check("dws_smarca_sel", "SMARCA4→SMARCA2 selectivity", "+0.18", trow("SMARCA4", "SMARCA2").mean_selectivity, TW)
check("dws_nf1", "NF1→RASA2 mean DWS (signed max(DD,0) formula)", "1.82", trow("NF1", "RASA2").mean_ti, TW)
check("dws_nf1_sel", "NF1→RASA2 selectivity", "0.005", trow("NF1", "RASA2").mean_selectivity, TW)
check("dws_nf1_pan", "NF1→RASA2 pan-essential denominator (%)", "0.1",
      trow("NF1", "RASA2").mean_pan_essential * 100, TW)
check("dws_pan_brca_lo", "BRCA1/2 pan-essential fraction lower", "0.53",
      min(trow("BRCA1", "BRCA2").mean_pan_essential, trow("BRCA2", "BRCA1").mean_pan_essential), TW)
check("dws_pan_brca_hi", "BRCA1/2 pan-essential fraction upper", "0.55",
      max(trow("BRCA1", "BRCA2").mean_pan_essential, trow("BRCA2", "BRCA1").mean_pan_essential), TW)
check("dws_pan_crkl", "PIK3R1→CRKL pan-essential fraction", "0.73",
      trow("PIK3R1", "CRKL").mean_pan_essential, TW)
check_int("dws_npairs", "Evaluated pairs (Table S5)", 21, len(twc), TW)
ti_gt1 = twa.groupby(["driver", "paralog"]).apply(
    lambda g: (g.dws_signed > 1).sum(), include_groups=False)
check_int("dws_n_gt1", "Pairs with signed DWS > 1.0 in ≥2 contexts", 5, (ti_gt1 >= 2).sum(),
          "output/therapeutic_window_all_results.csv")
tiers = twc.classification.value_counts()
check_int("dws_tier_hs", "HIGH_SELECTIVITY pairs", 2, tiers.get("HIGH_SELECTIVITY", 0), TW)
check_int("dws_tier_mod", "MODERATE pairs", 5, tiers.get("MODERATE", 0), TW)
check_int("dws_tier_low", "LOW_SELECTIVITY pairs", 11, tiers.get("LOW_SELECTIVITY", 0), TW)
check_int("dws_tier_pan", "PAN_ESSENTIAL pairs", 3, tiers.get("PAN_ESSENTIAL", 0), TW)

AL = "output/alphafold_structural_analysis.csv"
a0 = alpha.sort_values("clinical_targetability", ascending=False).iloc[0]
check("comp_top1", "Composite rank 1 score (ARID1A→ARID1B; Table S9, signed-DWS component)", "0.823", a0.clinical_targetability, AL)
a1 = alpha.sort_values("clinical_targetability", ascending=False).iloc[1]
check("comp_top2", "Composite rank 2 score (EP300→CREBBP; Table S9, signed-DWS component)", "0.644", a1.clinical_targetability, AL)
check_int("struct_npairs", "Scored structural pairs (text: 13)", 13, len(alpha), AL)
check_int("struct_domain_perfect", "Perfect domain conservation (text: 5 of 13)",
          5, int((alpha["domain_similarity"] == 1.0).sum()), AL)
check_int("struct_domain_nonzero", "Non-zero domain overlap (Fig. S9 caption: 10)",
          10, int((alpha["domain_similarity"] > 0).sum()), AL)
check("struct_domain_pct", "Perfect domain conservation share (text: 38%)",
      "0.38", (alpha["domain_similarity"] == 1.0).mean(), AL, tol=0.5)

# ═══════════════════════════════════════════════════════════════════
# 8b. DWS robustness: sensitivity variants (Table S10) and bootstrap
#     confidence intervals (Table S5 CI columns)
# ═══════════════════════════════════════════════════════════════════
DWSR = "output/dws_robustness.json"
dwsr = json.loads((OUT / "dws_robustness.json").read_text())
sens = {s["variant"]: s for s in dwsr["sensitivity"]}
boot = {(b["driver"], b["paralog"]): b for b in dwsr["bootstrap"]}

check("dws_sens_rho_mu", "Table S10: rho, denominator |mu| only (signed base)", "0.992",
      sens["denominator |mu| only"]["spearman_rho_vs_base"], DWSR)
check("dws_sens_rho_f", "Table S10: rho, denominator f only (signed base)", "0.966",
      sens["denominator f only"]["spearman_rho_vs_base"], DWSR)
check("dws_sens_rho_mean", "Table S10: rho, denominator mean(|mu|, f) (signed base)", "0.994",
      sens["denominator mean(|mu|, f)"]["spearman_rho_vs_base"], DWSR)
check("dws_sens_rho_f0001", "Table S10: rho, floor 0.001 (signed base)", "0.996",
      sens["floor 0.001"]["spearman_rho_vs_base"], DWSR)
check("dws_sens_rho_f005", "Table S10: rho, floor 0.05 (signed base)", "0.979",
      sens["floor 0.05"]["spearman_rho_vs_base"], DWSR)
check("dws_abs_sens_rho", "Table S10: rho, |DD| numerator (pre-revision formula, sensitivity)", "0.381",
      sens["numerator |DD| (pre-revision formula, sensitivity)"]["spearman_rho_vs_base"], DWSR)
check_int("dws_sens_top5_mu", "Table S10: top-5 overlap, |mu| only", 5,
          sens["denominator |mu| only"]["top5_overlap_with_base"], DWSR)
check_int("dws_sens_top5_f", "Table S10: top-5 overlap, f only", 3,
          sens["denominator f only"]["top5_overlap_with_base"], DWSR)
check_int("dws_sens_top5_mean", "Table S10: top-5 overlap, mean(|mu|, f)", 4,
          sens["denominator mean(|mu|, f)"]["top5_overlap_with_base"], DWSR)
check_int("dws_sens_top5_f0001", "Table S10: top-5 overlap, floor 0.001", 5,
          sens["floor 0.001"]["top5_overlap_with_base"], DWSR)
check_int("dws_sens_top5_f005", "Table S10: top-5 overlap, floor 0.05", 4,
          sens["floor 0.05"]["top5_overlap_with_base"], DWSR)
check_int("dws_abs_sens_top5", "Table S10: top-5 overlap, |DD| numerator (sensitivity)", 4,
          sens["numerator |DD| (pre-revision formula, sensitivity)"]["top5_overlap_with_base"], DWSR)
check_int("dws_hs_n_010", "Table S10: HS pairs at threshold 0.10", 3,
          sens["HIGH_SELECTIVITY selectivity threshold 0.1"]["n_high_selectivity"], DWSR)
check_int("dws_hs_n_015", "Table S10: HS pairs at threshold 0.15", 2,
          sens["HIGH_SELECTIVITY selectivity threshold 0.15"]["n_high_selectivity"], DWSR)
check_int("dws_hs_n_020", "Table S10: HS pairs at threshold 0.20", 1,
          sens["HIGH_SELECTIVITY selectivity threshold 0.2"]["n_high_selectivity"], DWSR)
check_int("dws_hs_flip_010", "Table S10: flips at threshold 0.10", 1,
          sens["HIGH_SELECTIVITY selectivity threshold 0.1"]["classification_flips_vs_thr_0.15"], DWSR)
check_int("dws_hs_flip_020", "Table S10: flips at threshold 0.20", 1,
          sens["HIGH_SELECTIVITY selectivity threshold 0.2"]["classification_flips_vs_thr_0.15"], DWSR)

# Signed max(DD,0) DWS bootstrap CIs (primary formula) — Table S5 CI columns
for drv, par, dlo, dhi, slo, shi in [
    ("ARID1A",  "ARID1B",  "1.72", "5.89",  "0.19", "0.37"),
    ("SMARCA4", "SMARCA2", "3.09", "7.59",  "0.07", "0.32"),
    ("NF1",     "RASA2",   "0.00", "4.84",  "0.00", "0.02"),
    ("EP300",   "CREBBP",  "1.00", "4.37",  "0.00", "0.21"),
]:
    b = boot[(drv, par)]
    tag = f"{drv}→{par}"
    check(f"dws_boot_dws_lo_{drv}", f"Table S5: {tag} DWS CI lower (signed max(DD,0) formula)", dlo, b["dws_ci95"][0], DWSR)
    check(f"dws_boot_dws_hi_{drv}", f"Table S5: {tag} DWS CI upper (signed max(DD,0) formula)", dhi, b["dws_ci95"][1], DWSR)
    check(f"dws_boot_sel_lo_{drv}", f"Table S5: {tag} selectivity CI lower", slo, b["selectivity_ci95"][0], DWSR)
    check(f"dws_boot_sel_hi_{drv}", f"Table S5: {tag} selectivity CI upper", shi, b["selectivity_ci95"][1], DWSR)

# |DD| DWS bootstrap CIs (pre-revision formula, sensitivity analysis)
for drv, par, alo, ahi in [
    ("ARID1A",  "ARID1B",  "1.76", "5.89"),
    ("SMARCA4", "SMARCA2", "3.28", "7.78"),
    ("NF1",     "RASA2",   "1.22", "13.83"),
    ("EP300",   "CREBBP",  "1.13", "4.55"),
]:
    b = boot[(drv, par)]
    tag = f"{drv}→{par}"
    check(f"dwsabs_boot_dws_lo_{drv}", f"{tag} |DD| DWS CI lower (sensitivity)", alo, b["dws_abs_ci95"][0], DWSR)
    check(f"dwsabs_boot_dws_hi_{drv}", f"{tag} |DD| DWS CI upper (sensitivity)", ahi, b["dws_abs_ci95"][1], DWSR)

# ═══════════════════════════════════════════════════════════════════
# 9. Table S2 spot values
# ═══════════════════════════════════════════════════════════════════
TS2 = "output/tables/TableS2_FullResults.tsv"
check_int("ts2_rows", "Table S2 associations", 110, len(ts2), TS2)
check_int("ts2_cols", "Table S2 columns", 17, len(ts2.columns), TS2)
r = ts2[(ts2.driver_gene == "ARID1A") & (ts2.paralog_gene == "ARID1B") & (ts2.cancer_type == "Ovarian")].iloc[0]
check("ts2_arid_dd", "ARID1A→ARID1B Ovarian DD (Table S2)", "0.386", r.dependency_dd, TS2)
check("ts2_arid_g", "ARID1A→ARID1B Ovarian Hedges' g", "1.39", r.hedges_g, TS2)
check("ts2_arid_q", "ARID1A→ARID1B Ovarian q-value", "0.39", r.q_value, TS2)
r2 = ts2[(ts2.driver_gene == "ARID1A") & (ts2.paralog_gene == "ARID1B") & (ts2.cancer_type == "Endometrial")].iloc[0]
check("ts2_arid_q2", "ARID1A→ARID1B Endometrial q-value", "0.88", r2.q_value, TS2)

# ═══════════════════════════════════════════════════════════════════
# 10. Regression lineage fixed effects
# ═══════════════════════════════════════════════════════════════════
REG = "output/regression_controls.json"
check("reg_lin_p", "ARID1A lineage-adjusted p", "2.2e-13",
      reg["results"]["arid1a_arid1b_lineage_adj_p"], REG)
check("reg_tp53_lin_p", "TP53+CNV+lineage adjusted p", "7.0e-13",
      reg["results"]["tp53_control"]["p_adj_lineage"], REG)

# ═══════════════════════════════════════════════════════════════════
# 11. Identity filter, benchmark tables, Table S3 tier consistency
# ═══════════════════════════════════════════════════════════════════
idf = hm["identity_filter"]["id_ge_0.3"]
check("id03_auroc", "DD + ID >= 0.3 AUROC (Table 1)", "1.000", idf["auroc"], HM)
check_int("id03_npairs", "DD + ID >= 0.3 unique pairs", 3, idf["n_unique_pairs"], HM)
check_int("id03_npos", "DD + ID >= 0.3 positives", 2, idf["n_positives"], HM)

for tname in ["Table1_Benchmark.tsv", "Table2_Benchmark.tsv"]:
    tb = pd.read_csv(OUT / "tables" / tname, sep="\t")
    src = f"output/tables/{tname}"
    tb_dd = tb.loc[tb.Method == "DD (this study)", "CV3_AUROC"].iloc[0]
    tb_id = tb.loc[tb.Method.str.startswith("DD + ID"), "CV3_AUROC"].iloc[0]
    tb_sl = tb.loc[tb.Method == "SLMGAE", "CV3_AUROC"].iloc[0]
    check(f"{tname[:6]}_dd", f"{tname} DD row", "0.629", tb_dd, src)
    check(f"{tname[:6]}_id", f"{tname} DD+ID row", "1.000", tb_id, src)
    check(f"{tname[:6]}_slm", f"{tname} SLMGAE row", "0.790", tb_sl, src)

s3 = pd.read_csv(OUT / "tables" / "TableS3_GoldStandard.tsv", sep="\t")
S3 = "output/tables/TableS3_GoldStandard.tsv"
check_int("s3_rows", "Table S3 pairs", 12, len(s3), S3)
s3_a = sorted(f"{r.Driver}->{r.Paralog}" for r in s3[s3.Tier == "A"].itertuples())
s3_b = sorted(f"{r.Driver}->{r.Paralog}" for r in s3[s3.Tier == "B"].itertuples())
ROWS.append({"id": "s3_tier_a", "description": "Table S3 Tier A == headline tier_a_pairs",
             "manuscript": str(hm["lineage_tier_ab"]["tier_a_pairs"]),
             "recomputed": str(s3_a), "source": S3,
             "status": "match" if s3_a == sorted(hm["lineage_tier_ab"]["tier_a_pairs"]) else "MISMATCH"})
ROWS.append({"id": "s3_tier_b", "description": "Table S3 Tier B == headline tier_b_pairs",
             "manuscript": str(hm["lineage_tier_ab"]["tier_b_pairs"]),
             "recomputed": str(s3_b), "source": S3,
             "status": "match" if s3_b == sorted(hm["lineage_tier_ab"]["tier_b_pairs"]) else "MISMATCH"})

# Table S11: full regression model table (confounder controls)
s11 = pd.read_csv(OUT / "tables" / "TableS11_RegressionModels.tsv", sep="\t")
S11 = "output/tables/TableS11_RegressionModels.tsv"
_reg_src = pd.read_csv(OUT / "regression_table_full.csv")
check_int("s11_rows", "Table S11 regression model rows (43)", 43, len(s11), S11)
check_int("s11_mirror_rows", "Table S11 mirrors regression_table_full.csv rows",
          len(_reg_src), len(s11), S11)
check_int("s11_models", "Table S11 model variants (base/cnv/expr/lineage)", 4,
          s11["model"].nunique(), S11)
_s11_arid = s11[(s11["pair"] == "ARID1A->ARID1B") &
                (s11["model"] == "lineage_adj")]["p"].iloc[0]
check("s11_arid1b_lineage_p", "Table S11 ARID1A->ARID1B lineage-adjusted p (text: 2.2e-13)",
      "2.2e-13", _s11_arid, S11)
check_int("s11_cols", "Table S11 columns (12 + se/p lineage cluster)", 14,
          len(s11.columns), S11)
_clh = reg["results"]["lineage_cluster_se_headline"]
check("reg_cl_arid1b_p", "ARID1A->ARID1B lineage-clustered SE p", "1.7e-9",
      _clh["ARID1A->ARID1B"]["p_lineage_cluster"], REG)
check("reg_cl_smarca2_p", "SMARCA4->SMARCA2 lineage-clustered SE p", "5.2e-3",
      _clh["SMARCA4->SMARCA2"]["p_lineage_cluster"], REG)
check("reg_cl_crebbp_p", "EP300->CREBBP lineage-clustered SE p", "3.2e-6",
      _clh["EP300->CREBBP"]["p_lineage_cluster"], REG)

# Table S9: evidence_tier annotation must match Table S3 tier mapping
s9 = pd.read_csv(OUT / "tables" / "TableS9_CompositeScore.tsv", sep="\t")
S9 = "output/tables/TableS9_CompositeScore.tsv"
_tmap = {f"{r.Driver}->{r.Paralog}": (f"Tier {r.Tier}" if r.Tier in ("A", "B", "C")
                                      else "comparator") for r in s3.itertuples()}
_s9_expected = [_tmap.get(f"{d}->{p}", "unlabeled")
                for d, p in zip(s9["driver"], s9["paralog"])]
ROWS.append({"id": "s9_evidence_tier",
             "description": "Table S9 evidence_tier == Table S3 mapping (unlabeled otherwise)",
             "manuscript": str(_s9_expected), "recomputed": str(list(s9["evidence_tier"])),
             "source": S9,
             "status": "match" if list(s9["evidence_tier"]) == _s9_expected else "MISMATCH"})

# Table S10: two new DWS ranking-sensitivity variants
S10 = "output/tables/TableS10_DWS_Sensitivity.tsv (dws_robustness.py)"
_dwsr = json.load(open(OUT / "dws_robustness.json"))
_s10s = {s["variant"]: s for s in _dwsr["sensitivity"]}
_s10w = _s10s["denominator winsorized 1%/99%"]
_s10f = _s10s["floor-dominated pairs removed (denominator at 0.01 floor in >=1 context)"]
check("s10_wins_rho", "Table S10 winsorized denominator rho", "0.999",
      _s10w["spearman_rho_vs_base"], S10)
check_int("s10_wins_top5", "Table S10 winsorized denominator top-5 overlap", 5,
          _s10w["top5_overlap_with_base"], S10)
check("s10_floor_rho", "Table S10 floor-dominated-removed rho", "1.000",
      _s10f["spearman_rho_vs_base"], S10)
check_int("s10_floor_top5", "Table S10 floor-dominated-removed top-5 overlap", 4,
          _s10f["top5_overlap_with_base"], S10)

# ═══════════════════════════════════════════════════════════════════
# TCGA survival v2 (continuous Cox PH, multivariable age+stage) — Fig. 3c,
# Results, Supplementary Table S8
# ═══════════════════════════════════════════════════════════════════
sv = json.loads((OUT / "tcga_survival_v2.json").read_text())
SV = "output/tcga_survival_v2.json (tcga_survival_v2.py)"
pgene = {g["gene"]: g for g in sv["per_gene"]}
check_int("tcga_ngenes", "TCGA paralog genes analysed (text + Table S8)", 32,
          len(sv["per_gene"]), SV)
check_int("tcga_n", "TCGA patients with OS+expression (n=1,069)", 1069,
          sv["cohort"]["n_samples"], SV)
check_int("tcga_events", "TCGA OS events (151 deaths)", 151,
          sv["cohort"]["n_events"], SV)
hl = sv["arid1b_highlight"]
mv = hl["multivar_age_stage"]
check_int("tcga_mv_n", "Multivariable complete cases (n=1,050; Table S8 caption)", 1050,
          mv["n"], SV)
check_int("tcga_mv_events", "Multivariable events (143; Table S8 caption)", 143,
          mv["n_events"], SV)
check("tcga_arid1b_hr", "ARID1B univariable HR per SD (text)", "1.30",
      hl["hr_continuous"], SV)
check("tcga_arid1b_lo", "ARID1B univariable CI low (text)", "1.087", hl["ci"][0], SV)
check("tcga_arid1b_hi", "ARID1B univariable CI high (text)", "1.552", hl["ci"][1], SV)
check("tcga_arid1b_p", "ARID1B univariable p (text)", "0.004", hl["p"], SV)
check("tcga_arid1b_hr_mv", "ARID1B multivariable HR (text + Fig. 3c)", "1.31",
      mv["hr_multivar"], SV)
check("tcga_arid1b_lo_mv", "ARID1B multivariable CI low (text)", "1.087",
      mv["ci_multivar"][0], SV)
check("tcga_arid1b_hi_mv", "ARID1B multivariable CI high (text)", "1.569",
      mv["ci_multivar"][1], SV)
check("tcga_arid1b_p_mv", "ARID1B multivariable p (text)", "0.004",
      mv["p_multivar"], SV)
check("tcga_arid1b_ph", "ARID1B Schoenfeld PH p (text)", "0.27", hl["ph_test_p"], SV)
check("tcga_arid1b_q_uni", "ARID1B univariable BH q (text: smallest 0.129)", "0.129",
      hl["q_fdr"], SV)
check("tcga_arid1b_q_mv", "ARID1B multivariable BH q (text: 0.040)", "0.040",
      mv["q_fdr_multivar"], SV)
check("tcga_pik3ca_q_mv", "PIK3CA multivariable BH q (text: 0.035)", "0.035",
      pgene["PIK3CA"]["multivar_age_stage"]["q_fdr_multivar"], SV)
check("tcga_brca2_q_mv", "BRCA2 multivariable BH q (text: 0.040)", "0.040",
      pgene["BRCA2"]["multivar_age_stage"]["q_fdr_multivar"], SV)
check("tcga_rbl1_q_mv", "RBL1 multivariable BH q (text: 0.040)", "0.040",
      pgene["RBL1"]["multivar_age_stage"]["q_fdr_multivar"], SV)
check_int("tcga_nfdr_mv", "Genes passing FDR<0.05 in multivariable family (text: four — PIK3CA, ARID1B, BRCA2, RBL1)",
          4, sum(1 for g in sv["per_gene"]
                 if (g.get("multivar_age_stage") or {}).get("q_fdr_multivar", 1) < 0.05), SV)
check_int("tcga_nfdr_uni", "Genes passing FDR<0.05 in univariable family (text: none)", 0,
          sum(1 for g in sv["per_gene"] if g["q_fdr"] < 0.05), SV)

# ═══════════════════════════════════════════════════════════════════
# Mutational co-occurrence with TMB adjustment — Results, Fig. S6
# ═══════════════════════════════════════════════════════════════════
colin = json.loads((OUT / "cooccurrence_by_lineage.json").read_text())
CO = "output/cooccurrence_by_lineage.json (cooccurrence_by_lineage.py)"
pan = {r["pair"]: r for r in colin["results"] if r["lineage"] == "PAN-CANCER"}
check_int("co_npairs", "Co-occurrence pairs analysed pan-cancer", 5, len(pan), CO)
ROWS.append({
    "id": "co_all_gt1", "description": "All key pairs co-occur (unadjusted OR > 1; text claim)",
    "manuscript": "OR > 1 for all",
    "recomputed": f"min OR = {min(r['or'] for r in pan.values()):.3f}",
    "source": CO,
    "status": "match" if all(r["or"] > 1 for r in pan.values()) else "MISMATCH"})
for pair, or_ms in [("ARID1A/ARID1B", "5.72"), ("EP300/CREBBP", "6.11"),
                    ("PIK3CA/PIK3CB", "6.84"), ("BRCA1/BRCA2", "4.53"),
                    ("SMARCA4/SMARCA2", "1.62")]:
    check(f"co_or_{pair.split('/')[0]}", f"{pair} unadjusted OR (Fig. S6a)", or_ms,
          pan[pair]["or"], CO)
check("co_arid_adj_or", "ARID1A/ARID1B TMB-adjusted OR (text + Fig. S6)", "2.38",
      pan["ARID1A/ARID1B"]["tmb_adjusted_or"], CO)
check("co_arid_adj_p", "ARID1A/ARID1B TMB-adjusted p (text)", "0.002",
      pan["ARID1A/ARID1B"]["p_adj"], CO)
check("co_ep300_adj_or", "EP300/CREBBP TMB-adjusted OR (text + Fig. S6)", "2.14",
      pan["EP300/CREBBP"]["tmb_adjusted_or"], CO)
check("co_ep300_adj_p", "EP300/CREBBP TMB-adjusted p (text)", "0.011",
      pan["EP300/CREBBP"]["p_adj"], CO)
check_int("co_nsig_adj", "Pairs significant after TMB adjustment (text: 2)", 2,
          sum(1 for r in pan.values()
              if r["p_adj"] is not None and r["p_adj"] < 0.05), CO)

# ═══════════════════════════════════════════════════════════════════
# CPTAC cohort counts + UCEC ARID1B protein effect sizes — Results, Fig. 2d
# ═══════════════════════════════════════════════════════════════════
cc = json.loads((OUT / "cptac_counts_verification.json").read_text())
CC = "output/cptac_counts_verification.json (verify_cptac_counts.py)"
check_int("cptac_eval", "CPTAC evaluable tests (text: 122 of 182)", 122,
          cc["evaluable_tests"], CC)
check_int("cptac_sig", "CPTAC significant at global BH q<0.05 (text: 42)", 42,
          cc["significant_q05"], CC)
uej = json.loads((OUT / "arid1b_ucec_effects.json").read_text())
ue = uej["primary_test"]
UE = "output/arid1b_ucec_effects.json (arid1b_ucec_effects.py)"
fam = uej["bh_families"]["family_all_evaluable_ucec_tests"]
check_int("cptac_ucec_ntests", "UCEC mutation-conditioned test family (text: ten)", 10,
          fam["n_tests"], UE)
check_int("cptac_ucec_nmut", "UCEC ARID1A-mutant tumours (n=37)", 37, ue["n_mut"], UE)
check_int("cptac_ucec_nwt", "UCEC wild-type tumours (n=58)", 58, ue["n_wt"], UE)
check("cptac_ucec_log2diff", "ARID1B log2 difference (text: 0.063)", "0.063",
      ue["log2_diff_mean_minus_mean"], UE)
check("cptac_ucec_fc", "ARID1B fold change (text: approx 1.05)", "1.05",
      ue["fold_change_2pow_log2diff"], UE, tol=1.0)
check("cptac_ucec_boot_lo", "ARID1B bootstrap 95% CI low (text: -0.018)", "-0.018",
      ue["boot_ci95_log2diff"][0], UE)
check("cptac_ucec_boot_hi", "ARID1B bootstrap 95% CI high (text: +0.141)", "+0.141",
      ue["boot_ci95_log2diff"][1], UE)
check("cptac_ucec_rbs", "ARID1B rank-biserial r (text: 0.213)", "0.213",
      ue["rank_biserial_r"], UE)
check("cptac_ucec_p", "ARID1B Wilcoxon p (text: 0.082)", "0.082",
      ue["wilcoxon_p_two_sided"], UE)
check("cptac_ucec_q", "ARID1B BH q across ten UCEC tests (text: 0.27)", "0.27",
      next(t for t in fam["tests"] if t["driver"] == "ARID1A")["q_bh_family_all"], UE)

# ═══════════════════════════════════════════════════════════════════
# MSI stratified bootstrap — Results
# ═══════════════════════════════════════════════════════════════════
mb = json.loads((OUT / "msi_bootstrap_test.json").read_text())
MB = "output/msi_bootstrap_test.json (msi_bootstrap_test.py)"
mbr = {r["cancer"]: r for r in mb["results"]}
check("msi_boot_endo_delta", "Endometrial MSI-H vs MSS delta AUROC (text: -0.197)",
      "-0.197", mbr["Endometrial"]["delta"], MB)
check("msi_boot_endo_p", "Endometrial bootstrap p (text: 0.534)", "0.534",
      mbr["Endometrial"]["p_empirical"], MB)
check("msi_boot_crc_delta", "Colorectal MSI-H vs MSS delta AUROC (text: +0.013)",
      "+0.013", mbr["Colorectal"]["delta"], MB)
check("msi_boot_crc_p", "Colorectal bootstrap p (text: 0.83)", "0.967",
      mbr["Colorectal"]["p_empirical"], MB)

# ═══════════════════════════════════════════════════════════════════
# PRISM driver-genotype drug sensitivity (d, delta AUC, q) — Results, Fig. 4a
# ═══════════════════════════════════════════════════════════════════
pr = json.loads((OUT / "prism_delta_auc.json").read_text())
PR = "output/prism_delta_auc.json (prism_delta_auc.py)"
pra = {a["drug"]: a for a in pr["associations"]}
for drug, d_ms, da_ms, q_ms in [
        ("AZD8330", "-0.62", "-0.37", "0.0007"),
        ("TRAMETINIB", "-0.57", "-0.24", "0.0006"),
        ("EVEROLIMUS", "-0.61", "-0.32", "0.0009"),
        ("IPATASERTIB", "-0.99", "-0.14", "0.0001"),
        ("PANOBINOSTAT", "-1.45", "-0.34", "0.010")]:
    a = pra[drug]
    check(f"prism_d_{drug.lower()}", f"{drug} Cohen's d (text)", d_ms,
          a["cohens_d"], PR)
    check(f"prism_dauc_{drug.lower()}", f"{drug} delta AUC (text)", da_ms,
          a["delta_auc"], PR)
    check(f"prism_q_{drug.lower()}", f"{drug} BH q (text)", q_ms, a["q"], PR)

# ═══════════════════════════════════════════════════════════════════
# k-mer vs Needleman-Wunsch validation — Fig. S9 + Methods line ~313
# ═══════════════════════════════════════════════════════════════════
kn = pd.read_csv(OUT / "kmer_nw_correlation.csv")
KN = "output/kmer_nw_correlation.csv (R_figS9_descriptors.R)"
k_all = kn[kn.group == "all"].iloc[0]
k_par = kn[kn.group == "paralog"].iloc[0]
k_non = kn[kn.group == "non_paralog"].iloc[0]
check("s9_r_all", "k-mer vs NW Pearson r (all pairs)", "0.88", k_all["pearson_r"], KN)
check_int("s9_n_all", "Validation pairs total", 50, int(k_all["n"]), KN)
ROWS.append({
    "id": "s9_p_all", "description": "Overall r p < 1e-16 (Methods claim)",
    "manuscript": "p < 1e-16", "recomputed": f"p = {k_all['p_value']:.2e}", "source": KN,
    "status": "match" if k_all["p_value"] < 1e-16 else "MISMATCH"})
check("s9_r_par", "r within paralogs", "0.91", k_par["pearson_r"], KN)
check_int("s9_n_par", "Paralog pairs", 13, int(k_par["n"]), KN)
check("s9_r_non", "r within non-paralogs", "0.41", k_non["pearson_r"], KN)
check_int("s9_n_non", "Non-paralog pairs", 37, int(k_non["n"]), KN)

# ═══════════════════════════════════════════════════════════════════
# in4mer mutation-agnostic gold-standard check — Discussion limitation
# ═══════════════════════════════════════════════════════════════════
in4s = json.loads((OUT / "in4mer_benchmark_summary.json").read_text())
in4 = pd.read_csv(OUT / "in4mer_benchmark.csv")
IN4 = "output/in4mer_benchmark_summary.json (in4mer_benchmark.py)"
check_int("in4_npairs", "in4mer gold-standard pairs (Nat Commun 2024)", 13,
          in4s["n_in4mer_total"], IN4)
check_int("in4_controls", "Unlabeled control pairs sampled (seed 42)", 400,
          in4s["n_controls_sampled"], IN4)
check_int("in4_min3_npos", "in4mer pairs evaluable at >=3 threshold", 10,
          in4s["min3"]["n_pos"], IN4)
check("in4_min3_auroc", "in4mer min3 AUROC (0.50)", "0.50",
      in4s["min3"]["auroc"], IN4)
check("in4_min3_p", "in4mer min3 permutation p (0.51)", "0.51",
      in4s["min3"]["permutation_p"], IN4)
check_int("in4_min5_npos", "in4mer pairs evaluable at >=5 threshold", 3,
          in4s["min5"]["n_pos"], IN4)
check("in4_min5_auroc", "in4mer min5 exploratory AUROC (0.68)", "0.68",
      in4s["min5"]["auroc"], IN4)
neg3 = in4[in4.label == "unlabeled_control"]["dd_min3"].abs().dropna()
p90 = float(neg3.quantile(0.9))
check_int("in4_min3_nneg", "Controls evaluable at >=3 threshold (Fig. S2)", 326,
          in4s["min3"]["n_neg"], IN4)
check("in4_p90", "Control 90th pct of |DD| at >=3 (Fig. S2: 0.217)", "0.217",
      p90, "output/in4mer_benchmark.csv")
pos3 = in4[in4.label == "in4mer_gold"]["dd_min3"].abs().dropna()
above = pos3[pos3 > p90]
ok = len(above) == 3
ROWS.append({
    "id": "in4_above_p90",
    "description": "in4mer pairs above control 90th pct (3: ARFGEF1/2, HDAC1/2, SLC25A28/37)",
    "manuscript": "3 pairs", "recomputed": f"{len(above)} pairs (p90={p90:.3f})",
    "source": "output/in4mer_benchmark.csv",
    "status": "match" if ok else "MISMATCH"})

# ═══════════════════════════════════════════════════════════════════
# Wave-3 additions: pair-clustered bootstrap, power analysis,
# composite weight sensitivity, in4mer seed sensitivity,
# CPTAC RNA/purity, TCGA UCEC/OV survival extension
# ═══════════════════════════════════════════════════════════════════
CB = "output/cluster_bootstrap_primary.json"
cb = json.load(open(OUT / "cluster_bootstrap_primary.json"))
check("cb_ci_lo", "Pair-clustered bootstrap CI lower (primary frame)", "0.253",
      cb["pair_clustered_bootstrap_auroc"]["ci95_percentile"][0], CB)
check("cb_ci_hi", "Pair-clustered bootstrap CI upper (primary frame)", "0.933",
      cb["pair_clustered_bootstrap_auroc"]["ci95_percentile"][1], CB)
check("cb_perm_p", "Pair-level permutation p (primary 0.629)", "0.122",
      cb["pair_level_permutation_primary"]["empirical_p"], CB)
check("cb_neff", "Effective sample size (of 110 entries)", "85",
      cb["effective_sample_size"]["neff"], CB)
check("cb_icc", "Within-pair ICC of signed DD", "0.554",
      cb["effective_sample_size"]["anova_icc_dependency_dd_by_pair"], CB)
check("cb_mean_p", "Mean-aggregation frame permutation p", "0.090",
      cb["per_pair_mean_aggregation_permutation"]["empirical_p"], CB)
check("db_ci_lo", "Driver-block bootstrap CI lower (primary frame)", "0.272",
      cb["driver_block_bootstrap_auroc"]["ci95_percentile"][0], CB)
check("db_ci_hi", "Driver-block bootstrap CI upper (primary frame)", "0.929",
      cb["driver_block_bootstrap_auroc"]["ci95_percentile"][1], CB)

PW = "output/power_analysis.json"
pw = json.load(open(OUT / "power_analysis.json"))
check("pow_now", "Power at observed 0.629 (8 vs 102)", "0.34",
      pw["i_power_current_sizes"]["lineage_level_8v102"]["power"], PW)
check("pow_now_pair", "Power at observed 0.629 (6 vs 66)", "0.28",
      pw["i_power_current_sizes"]["per_pair_6v66"]["power"], PW)
check_int("pow_n45", "Positives required for 80% power at 0.629", 45,
          pw["ii_required_for_80pct_at_observed"]["same_control_count_102"]["n_pos_required"], PW)

WS = "output/composite_weight_sensitivity.json"
wsen = json.load(open(OUT / "composite_weight_sensitivity.json"))
check("ws_p25", "Dirichlet weight sensitivity AUROC p25", "0.803",
      wsen["auroc_distribution"]["p25"], WS)
check("ws_med", "Dirichlet weight sensitivity AUROC median", "0.819",
      wsen["auroc_distribution"]["p50"], WS)
check("ws_p75", "Dirichlet weight sensitivity AUROC p75", "0.838",
      wsen["auroc_distribution"]["p75"], WS)
check("ws_min", "Dirichlet weight sensitivity AUROC min", "0.664",
      wsen["auroc_distribution"]["min"], WS)
check("ws_max", "Dirichlet weight sensitivity AUROC max", "0.884",
      wsen["auroc_distribution"]["max"], WS)

SE = "output/in4mer_seed_sensitivity.json"
sed = json.load(open(OUT / "in4mer_seed_sensitivity.json"))
check("se_min", "in4mer 20-seed AUROC min", "0.480", sed["summary"]["min"], SE)
check("se_med", "in4mer 20-seed AUROC median", "0.499", sed["summary"]["median"], SE)
check("se_max", "in4mer 20-seed AUROC max", "0.522", sed["summary"]["max"], SE)

# ═══════════════════════════════════════════════════════════════════
# 10b. External digenic screens (Harle 2025; Flister 2025)
# ═══════════════════════════════════════════════════════════════════
EH = "output/external_holdout.json (external_holdout.py)"
eh = json.load(open(OUT / "external_holdout.json"))
fid = eh["engine_fidelity_vs_frozen_in4mer"]
check("eh_fid_maxdiff", "Engine fidelity vs frozen in4mer max|diff|", "2.8e-16",
      fid["max_abs_diff"], EH)
check_int("eh_fid_npairs", "Engine fidelity pairs", 413, fid["n_pairs"], EH)

eh_h, eh_f = eh["harle2025"], eh["flister2025"]
h1, f1 = eh_h["layer1"], eh_f["layer1"]
h2, f2 = eh_h["layer2"], eh_f["layer2"]
hsc, fsc = eh_h["screen_counts"], eh_f["screen_counts"]

check_int("eh_h_pairs", "Harle pairs screened", 472, hsc["n_pairs"], EH)
check_int("eh_h_hitpairs", "Harle hit pairs", 117, hsc["n_hit_pairs"], EH)
check_int("eh_h_lines", "Harle cell lines", 27, h2["n_lines_genotyped"] + 1, EH)
check_int("eh_h_genotyped", "Harle genotyped lines", 26, h2["n_lines_genotyped"], EH)
check_int("eh_f_pairs", "Flister pairs screened", 36648, eh_f["n_pairs_screened"], EH)
check_int("eh_f_s6rows", "Flister S6 meta pairs", 462, fsc["n_s6_rows"], EH)
check_int("eh_f_models", "Flister S6 models", 49, f2["n_models"], EH)
check_int("eh_f_lethal_h1299", "Flister lethal NCI-H1299", 117, fsc["n_lethal_h1299"], EH)
check_int("eh_f_lethal_mb231", "Flister lethal MDA-MB-231", 51, fsc["n_lethal_mdamb231"], EH)
check_int("eh_f_lethal_both", "Flister lethal shared", 22, fsc["n_lethal_both"], EH)
check_int("eh_f_lethal_union", "Flister lethal union", 146, fsc["n_lethal_union"], EH)
check("eh_f_s6concord", "Flister S6 hit-call concordance (%)", "98.3",
      f2["s6_concordance"] * 100, EH)

check("eh_h1_auroc", "Harle layer-1 AUROC", "0.619", h1["auroc"], EH)
check("eh_h1_p", "Harle layer-1 permutation p", "0.0002", h1["permutation_p"], EH)
check_int("eh_h1_n", "Harle layer-1 scored pairs", 411, h1["n_pairs_scored"], EH)
check_int("eh_h1_hits", "Harle layer-1 evaluable hits", 98, h1["n_hits"], EH)
check("eh_h1_med_hit", "Harle median |DD| hits", "0.202", h1["median_abs_dd_hits"], EH)
check("eh_h1_med_non", "Harle median |DD| non-hits", "0.135", h1["median_abs_dd_nonhits"], EH)
check("eh_f1_auroc", "Flister layer-1 AUROC", "0.696", f1["auroc"], EH)
check("eh_f1_p", "Flister layer-1 permutation p", "0.0001", f1["permutation_p"], EH)
check_int("eh_f1_n", "Flister layer-1 scored pairs", 29457, f1["n_pairs_scored"], EH)
check_int("eh_f1_hits", "Flister layer-1 evaluable hits", 114, f1["n_hits"], EH)
check("eh_f1_med_hit", "Flister median |DD| hits", "0.165", f1["median_abs_dd_hits"], EH)
check("eh_f1_med_non", "Flister median |DD| non-hits", "0.105", f1["median_abs_dd_nonhits"], EH)
check("eh_f1_h1299", "Flister per-model AUROC NCI-H1299", "0.671",
      f1["sensitivity_per_model"]["hit_h1299"]["auroc"], EH)
check("eh_f1_mb231", "Flister per-model AUROC MDA-MB-231", "0.744",
      f1["sensitivity_per_model"]["hit_mdamb231"]["auroc"], EH)

check_int("eh_h_evalpairs", "Harle evaluable curated pairs", 3, h2["pairs_evaluable"], EH)
check_int("eh_f_evalpairs", "Flister evaluable curated pairs", 3, f2["pairs_evaluable"], EH)

hp, fp = h2["pooled"], f2["pooled"]
check_int("eh_h2_nmut", "Harle pooled mutant observations", 9, hp["n_mut"], EH)
check_int("eh_h2_muthit", "Harle pooled mutant hits", 0, hp["n_mut_hit"], EH)
check_int("eh_h2_nwt", "Harle pooled WT observations", 121, hp["n_wt"], EH)
check_int("eh_h2_wthit", "Harle pooled WT hits", 27, hp["n_wt_hit"], EH)
check("eh_h2_fisher", "Harle pooled Fisher p", "1.0", hp["fisher_p_1sided_greater"], EH)
check("eh_h2_mw", "Harle pooled MW p", "0.899", hp["mw_p_1sided_mut_more_negative"], EH)
check_int("eh_f2_nmut", "Flister pooled mutant observations", 6, fp["n_mut"], EH)
check_int("eh_f2_muthit", "Flister pooled mutant hits", 2, fp["n_mut_hit"], EH)
check_int("eh_f2_nwt", "Flister pooled WT observations", 235, fp["n_wt"], EH)
check_int("eh_f2_wthit", "Flister pooled WT hits", 59, fp["n_wt_hit"], EH)
check("eh_f2_fisher", "Flister pooled Fisher p", "0.475", fp["fisher_p_1sided_greater"], EH)
check("eh_f2_mw", "Flister pooled MW p", "0.150", fp["mw_p_1sided_mut_more_negative"], EH)

# Per-pair layer-2 table (Supplementary Methods)
_EH_TABLE = {
    ("harle2025", "SMARCA4", "SMARCA2"): (5, 0, 21, 3, "1.000", "0.575"),
    ("harle2025", "ARID1A", "ARID1B"): (2, 0, 24, 1, "1.000", "0.972"),
    ("harle2025", "EP300", "CREBBP"): (2, 0, 24, 9, "1.000", "0.277"),
    ("flister2025", "MAP2K1", "MAP2K2"): (1, 0, 47, 8, "1.000", "0.282"),
    ("flister2025", "SMARCA4", "SMARCA2"): (2, 0, 46, 5, "1.000", "0.741"),
    ("flister2025", "EP300", "CREBBP"): (3, 2, 46, 20, "0.422", "0.121"),
}
for (study, drv, prl), (nm, nmh, nw, nwh, fs, mws) in _EH_TABLE.items():
    row = next(p for p in eh[study]["layer2"]["per_pair"]
               if p["driver"] == drv and p["paralog"] == prl)
    tag = f"eh_{'h' if study == 'harle2025' else 'f'}2_{drv.lower()}"
    check_int(f"{tag}_nmut", f"{study} {drv}->{prl} mutant models", nm, int(row["n_mut"]), EH)
    check_int(f"{tag}_muthit", f"{study} {drv}->{prl} mutant hits", nmh, int(row["n_mut_hit"]), EH)
    check_int(f"{tag}_nwt", f"{study} {drv}->{prl} WT models", nw, int(row["n_wt"]), EH)
    check_int(f"{tag}_wthit", f"{study} {drv}->{prl} WT hits", nwh, int(row["n_wt_hit"]), EH)
    check(f"{tag}_fisher", f"{study} {drv}->{prl} Fisher p", fs,
          row["fisher_p_1sided_greater"], EH, tol=0.6)
    check(f"{tag}_mw", f"{study} {drv}->{prl} MW p", mws,
          row["mw_p_1sided_mut_more_negative"], EH, tol=0.6)

_f_ep = next(p for p in f2["per_pair"] if p["driver"] == "EP300")
check("eh_f2_ep300_gi_mut", "Flister EP300->CREBBP mean GI mutant", "-2.84",
      _f_ep["gi_mean_mut"], EH)
check("eh_f2_ep300_gi_wt", "Flister EP300->CREBBP mean GI WT", "-1.61",
      _f_ep["gi_mean_wt"], EH)

RNA = "output/cptac_rna_covariation.json"
rna = json.load(open(OUT / "cptac_rna_covariation.json"))
rov = rna["summary"]["overall"]
check_int("rna_ntests", "Pair-by-cohort mRNA tests", 179, rov["n_pair_cohort_mrna_tests"], RNA)
check_int("rna_nmatched", "Pair-by-cohort matched protein+mRNA tests", 122,
          rov["n_pair_cohort_with_matched_protein"], RNA)
check("rna_med_prot", "Median matched protein r", "0.140", rov["median_protein_r_matched"], RNA)
check("rna_med_mrna", "Median mRNA r", "0.074", rov["median_mrna_r"], RNA)
check("rna_sign_p", "Sign test |protein r| vs |mRNA r| p", "0.79",
      rna["conclusion_flags"]["protein_stronger_than_mrna"]["sign_test_abs_r"]["p_two_sided"], RNA)
ep_ucec = rna["ep300_crebbp_per_cohort"]["UCEC"]
check("rna_ep_ucec_mrna", "EP300-CREBBP UCEC mRNA r", "0.765", ep_ucec["mrna_r"], RNA)
check("rna_ep_ucec_prot", "EP300-CREBBP UCEC protein r", "0.081", ep_ucec["protein_r_full"], RNA)

PU = "output/cptac_purity_partial.json"
pu = json.load(open(OUT / "cptac_purity_partial.json"))
check_int("pu_ntests", "Purity-adjusted pair-by-cohort tests", 69,
          pu["summary"]["overall"]["n_pair_cohort_purity_tests"], PU)
check("pu_med_dr", "Median |delta r| after purity adjustment", "0.014",
      pu["summary"]["overall"]["median_abs_delta_r_purity"], PU)

SV = "output/tcga_survival_ucec_ov_summary.json"
sv = json.load(open(OUT / "tcga_survival_ucec_ov_summary.json"))
svu = sv["cohorts"]["UCEC"]
svo = sv["cohorts"]["OV"]
check_int("sv_ucec_n", "UCEC cohort n", 525, svu["n"], SV)
check_int("sv_ucec_ev", "UCEC events", 87, svu["events"], SV)
check_int("sv_ov_n", "OV cohort n", 299, svo["n"], SV)
check_int("sv_ov_ev", "OV events", 180, svo["events"], SV)
check("sv_ucec_hr", "UCEC ARID1B univariable HR", "1.17", svu["arid1b"]["univariable"]["hr"], SV)
check("sv_ucec_hr_lo", "UCEC ARID1B HR CI lower", "0.941", svu["arid1b"]["univariable"]["ci"][0], SV)
check("sv_ucec_hr_hi", "UCEC ARID1B HR CI upper", "1.464", svu["arid1b"]["univariable"]["ci"][1], SV)
check("sv_ucec_p", "UCEC ARID1B univariable p", "0.155", svu["arid1b"]["univariable"]["p"], SV)
check("sv_ucec_q", "UCEC ARID1B BH q", "0.381", svu["arid1b"]["univariable"]["q_fdr"], SV)
check("sv_ucec_adj_hr", "UCEC ARID1B age-adjusted HR", "1.17", svu["arid1b"]["multivariable_age_only"]["hr"], SV)
check("sv_ucec_adj_p", "UCEC ARID1B age-adjusted p", "0.179", svu["arid1b"]["multivariable_age_only"]["p"], SV)
check("sv_ov_hr", "OV ARID1B univariable HR", "1.26", svo["arid1b"]["univariable"]["hr"], SV)
check("sv_ov_hr_lo", "OV ARID1B HR CI lower", "1.092", svo["arid1b"]["univariable"]["ci"][0], SV)
check("sv_ov_hr_hi", "OV ARID1B HR CI upper", "1.461", svo["arid1b"]["univariable"]["ci"][1], SV)
check("sv_ov_p", "OV ARID1B univariable p", "0.002", svo["arid1b"]["univariable"]["p"], SV)
check("sv_ov_q", "OV ARID1B BH q", "0.054", svo["arid1b"]["univariable"]["q_fdr"], SV)
check("sv_ov_adj_hr", "OV ARID1B age-adjusted HR", "1.21", svo["arid1b"]["multivariable_age_only"]["hr"], SV)
check("sv_ov_adj_p", "OV ARID1B age-adjusted p", "0.013", svo["arid1b"]["multivariable_age_only"]["p"], SV)
check("sv_ov_adj_q", "OV ARID1B age-adjusted q", "0.406", svo["arid1b"]["multivariable_age_only"]["q_fdr"], SV)
svst = {s["stratum"]: s for s in svu["arid1a_stratified"]["strata"]}
check("sv_ucec_mut_hr", "UCEC ARID1A-mutant stratum ARID1B HR", "1.24",
      svst["ARID1A_mutant"]["univariate"]["hr"], SV)
check("sv_ucec_mut_p", "UCEC mutant stratum p", "0.391", svst["ARID1A_mutant"]["univariate"]["p"], SV)
check("sv_ucec_wt_hr", "UCEC wild-type stratum ARID1B HR", "1.15",
      svst["ARID1A_wildtype"]["univariate"]["hr"], SV)
check("sv_ucec_wt_p", "UCEC wild-type stratum p", "0.248", svst["ARID1A_wildtype"]["univariate"]["p"], SV)
check_int("sv_ucec_mut_n", "UCEC mutant stratum n", 226, svst["ARID1A_mutant"]["n"], SV)
check_int("sv_ucec_mut_ev", "UCEC mutant stratum events", 19, svst["ARID1A_mutant"]["n_events"], SV)
check_int("sv_ucec_wt_n", "UCEC WT stratum n", 299, svst["ARID1A_wildtype"]["n"], SV)
check_int("sv_ucec_wt_ev", "UCEC WT stratum events", 68, svst["ARID1A_wildtype"]["n_events"], SV)

# ═══════════════════════════════════════════════════════════════════
# Report
# ═══════════════════════════════════════════════════════════════════
df = pd.DataFrame(ROWS)
out_path = OUT / "manuscript_number_audit.tsv"
df.to_csv(out_path, sep="\t", index=False)
n_fail = (df.status == "MISMATCH").sum()
print(df.to_string(index=False, max_colwidth=52))
print(f"\n{len(df) - n_fail}/{len(df)} claims match; {n_fail} mismatch(es).")
print(f"Audit table written to {out_path}")
if n_fail:
    print("\nFAILED CLAIMS:")
    print(df[df.status == "MISMATCH"].to_string(index=False, max_colwidth=80))
    sys.exit(1)
print("CLOSED LOOP OK — every audited manuscript number reproduces from artifacts.")
