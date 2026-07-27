"""Generate the machine-readable Additional-file TSV mirrors (Tables S5, S6,
S8, S9, S10) from pipeline outputs. All values are copied or reshaped from the
named output artifacts — nothing is recomputed, simulated, or hard-coded.

Inputs (all produced upstream by the analysis pipeline):
  output/therapeutic_window_paralog_classification.csv  (therapeutic_window.py)
  output/dws_robustness.json                            (dws_robustness.py)
  output/driver_mutation_rules.csv                      (make_driver_mutation_rules_table.py)
  output/tcga_survival_v2.json                          (tcga_survival_v2.py)
  output/alphafold_structural_analysis.csv              (alphafold_analysis.py)

Outputs (referenced in manuscript "Additional files" and supplementary.tex):
  output/tables/TableS5_DWS.tsv
  output/tables/TableS6_MutationRules.tsv
  output/tables/TableS8_BRCA_Survival.tsv
  output/tables/TableS9_CompositeScore.tsv
  output/tables/TableS10_DWS_Sensitivity.tsv
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)


def table_s5_dws() -> None:
    tw = pd.read_csv(OUT / "therapeutic_window_paralog_classification.csv")
    boot = json.loads((OUT / "dws_robustness.json").read_text())["bootstrap"]
    ci = pd.DataFrame(boot)[["driver", "paralog", "dws_ci95", "selectivity_ci95"]]
    df = tw.merge(ci, on=["driver", "paralog"], how="left", validate="one_to_one")
    df["dws_ci95_lo"] = df["dws_ci95"].map(lambda x: x[0])
    df["dws_ci95_hi"] = df["dws_ci95"].map(lambda x: x[1])
    df["selectivity_ci95_lo"] = df["selectivity_ci95"].map(lambda x: x[0])
    df["selectivity_ci95_hi"] = df["selectivity_ci95"].map(lambda x: x[1])

    def pair_type(row) -> str:
        pair = {row["driver"], row["paralog"]}
        if pair == {"BRCA1", "BRCA2"}:
            return "functional analog"
        if pair == {"STK11", "SIK1"}:
            return "partial homolog"
        return "sequence paralog"

    df["pair_type"] = df.apply(pair_type, axis=1)
    df = df.rename(columns={
        "mean_ti": "dws",
        "mean_dd": "abs_dd",
        "mean_selectivity": "selectivity",
        "mean_pan_essential": "pan_essential_fraction",
    })
    cols = ["driver", "paralog", "abs_dd", "dws", "dws_ci95_lo", "dws_ci95_hi",
            "selectivity", "selectivity_ci95_lo", "selectivity_ci95_hi",
            "pan_essential_fraction", "n_contexts", "classification", "pair_type"]
    df = df[cols].sort_values("dws", ascending=False)
    df.to_csv(TAB / "TableS5_DWS.tsv", sep="\t", index=False)
    print(f"  TableS5_DWS.tsv: {len(df)} rows")
    assert len(df) == 21, "expected 21 DWS pairs"


def table_s6_mutation_rules() -> None:
    df = pd.read_csv(OUT / "driver_mutation_rules.csv")
    df.to_csv(TAB / "TableS6_MutationRules.tsv", sep="\t", index=False)
    print(f"  TableS6_MutationRules.tsv: {len(df)} rows")


def table_s8_brca_survival() -> None:
    d = json.loads((OUT / "tcga_survival_v2.json").read_text())
    rows = []
    for g in d["per_gene"]:
        mv = g.get("multivar_age_stage", {})
        rows.append({
            "gene": g["gene"],
            "n": g["n"],
            "n_events": g["n_events"],
            "hr_univariable": g["hr_continuous"],
            "ci95_lo_univariable": g["ci"][0],
            "ci95_hi_univariable": g["ci"][1],
            "p_univariable": g["p"],
            "q_fdr_univariable": g["q_fdr"],
            "n_multivariable": mv.get("n"),
            "n_events_multivariable": mv.get("n_events"),
            "hr_multivariable": mv.get("hr_multivar"),
            "ci95_lo_multivariable": mv.get("ci_multivar", [None, None])[0],
            "ci95_hi_multivariable": mv.get("ci_multivar", [None, None])[1],
            "p_multivariable": mv.get("p_multivar"),
            "q_fdr_multivariable": mv.get("q_fdr_multivar"),
            "ph_test_p": g.get("ph_test_p"),
        })
    df = pd.DataFrame(rows).sort_values("p_multivariable")
    df.to_csv(TAB / "TableS8_BRCA_Survival.tsv", sep="\t", index=False)
    print(f"  TableS8_BRCA_Survival.tsv: {len(df)} rows")
    assert len(df) == 32, "expected 32 paralog genes"


def table_s9_composite_score() -> None:
    df = pd.read_csv(OUT / "alphafold_structural_analysis.csv")
    df = df.rename(columns={"mean_ti": "dws", "mean_dd": "abs_dd",
                            "mean_selectivity": "selectivity"})
    cols = ["driver", "paralog", "is_known_sl", "dws", "abs_dd", "selectivity",
            "structural_similarity", "domain_similarity", "druggability",
            "protac_score", "targetability", "clinical_targetability"]
    df = df[cols].sort_values("clinical_targetability", ascending=False)
    df.to_csv(TAB / "TableS9_CompositeScore.tsv", sep="\t", index=False)
    print(f"  TableS9_CompositeScore.tsv: {len(df)} rows")
    assert len(df) == 13, "expected 13 scored pairs"


def table_s10_dws_sensitivity() -> None:
    sens = json.loads((OUT / "dws_robustness.json").read_text())["sensitivity"]
    rows = []
    for s in sens:
        rows.append({
            "variant": s["variant"],
            "spearman_rho_vs_base": s.get("spearman_rho_vs_base"),
            "top5_overlap_with_base": s.get("top5_overlap_with_base"),
            "n_high_selectivity": s.get("n_high_selectivity"),
            "classification_flips_vs_thr_0.15": s.get("classification_flips_vs_thr_0.15"),
        })
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "TableS10_DWS_Sensitivity.tsv", sep="\t", index=False)
    print(f"  TableS10_DWS_Sensitivity.tsv: {len(df)} rows")
    assert len(df) == 9, "expected 9 sensitivity variants"


if __name__ == "__main__":
    table_s5_dws()
    table_s6_mutation_rules()
    table_s8_brca_survival()
    table_s9_composite_score()
    table_s10_dws_sensitivity()
    print(f"Additional-file TSVs written to {TAB}")
