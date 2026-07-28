"""Generate manuscript tables for the pan-cancer revision."""
import pandas as pd; import numpy as np
from pathlib import Path
OUT = Path(__file__).resolve().parent / "output" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

# ── Table S2: Full gynecological-cancer results (118 entries) ──
# Rebuilt from the pipeline output so q_values/DD always track the current run.
_cand = pd.read_csv(OUT.parent / "paralog_sl_candidates.csv")
_gyn3 = ["Ovarian", "Endometrial", "Cervical"]
_ts2 = _cand[_cand["cancer_type"].isin(_gyn3)].copy()
_ts2["novelty"] = _ts2["is_known_paralog_sl"].map({True: "Known", False: "Novel"})
_ts2_cols = ["driver_gene","paralog_gene","cancer_type","pcs","delta_expression",
             "necessity","dependency_dd","cohens_d","hedges_g","dd_p_value",
             "composite_score","novelty","q_value",
             "mutation_frequency","n_mut","n_wt","is_known_paralog_sl"]
_ts2 = (_ts2.sort_values(["is_known_paralog_sl","composite_score"],
                         ascending=[False,False])[_ts2_cols])
_ts2.to_csv(OUT / "TableS2_FullResults.tsv", sep="\t", index=False)
print(f"  TableS2_FullResults.tsv: {len(_ts2)} rows")

summary = pd.read_csv(OUT.parent / "solid_tumor_summary.csv")
summary = summary.dropna(subset=["dd_auroc"]).sort_values("dd_auroc", ascending=False)

# ── Table 1 (updated): Top de novo candidates ──
all_candidates = []
for _, r in summary.iterrows():
    safe = r["cancer"].replace("/","_").replace(" ","_")
    f = OUT.parent / f"solid_{safe}_results.csv"
    if not f.exists(): continue
    df = pd.read_csv(f)
    if "is_known_paralog_sl" not in df.columns: continue
    novel = df[~df["is_known_paralog_sl"]]
    if novel.empty: continue
    top = novel.nlargest(3, "composite_score")
    top["cancer_type"] = r["cancer"]
    all_candidates.append(top[["driver_gene","paralog_gene","cancer_type","pcs","dependency_dd","composite_score","is_known_paralog_sl"]])

if all_candidates:
    cand = pd.concat(all_candidates, ignore_index=True)
    cand = cand.sort_values("composite_score", ascending=False).head(12)
    cand.columns = ["Driver","Paralog","Cancer","PCS","DD","Composite","Known"]
    cand["Status"] = cand["Known"].map({True:"Known",False:"Novel"})
    cand.to_csv(OUT / "Table1_DeNovoCandidates.tsv", sep="\t", index=False)

# ── Table 2 (benchmark) ──
# Headline AUROCs are read from output/headline_metrics.json (written by
# compute_headline_metrics.py) — never hard-code this-study values here.
# Published values are literature constants (Feng et al. 2024, Suppl. Data 1,
# CV3 gene-pair isolation, NSMRand 1:1, complete dataset).
import json
_metrics_path = OUT.parent / "headline_metrics.json"
if not _metrics_path.exists():
    raise SystemExit("output/headline_metrics.json not found — "
                     "run compute_headline_metrics.py first")
_hm = json.loads(_metrics_path.read_text())
_pub = _hm["published_benchmarks"]["values"]
_dd_val = _hm["lineage_full"]["auroc"]
_dd_id_val = _hm.get("identity_filter", {}).get("id_ge_0.3", {}).get("auroc")
if _dd_id_val is None:
    raise SystemExit("identity-filter metric missing — run compute_sequence_identity.R "
                     "then compute_headline_metrics.py")
_THIS = "This study (recomputed)"
_LIT = "Feng et al. 2024 (CV3, SD1)"

bench = pd.DataFrame({
    "Method": ["SLMGAE","NSF4SL","GCATSL","GRSMF","PiLSL","KG4SL","SLGNN","PTGNN","DD (this study)","DD + ID≥0.3"],
    "CV3_AUROC": [f"{_pub['SLMGAE']:.3f}", f"{_pub['NSF4SL']:.3f}", f"{_pub['GCATSL']:.3f}",
                  f"{_pub['GRSMF']:.3f}", f"{_pub['PiLSL']:.3f}", f"{_pub['KG4SL']:.3f}",
                  f"{_pub['SLGNN']:.3f}", f"{_pub['PTGNN']:.3f}",
                  f"{_dd_val:.3f}", f"{_dd_id_val:.3f}"],
    "Architecture": ["Multi-view GAE","Contrastive","Graph Attention","Matrix Factor.","Pairwise GNN","KG Embed.","KG GNN","Pre-trained GNN","Univariate","Univariate+Filter"],
    "Interpretability": ["Low","Low","Low","Low","Low","Low","Low","Low","High","High"],
    "Source": [_LIT]*8 + [_THIS, _THIS],
})
bench.to_csv(OUT / "Table2_Benchmark.tsv", sep="\t", index=False)

# Manuscript Table 1 layout (Method / AUROC / Interpretability) — same data,
# fewer columns. Kept script-owned so it can never drift from Table 2 again.
bench[["Method", "CV3_AUROC", "Interpretability"]].to_csv(
    OUT / "Table1_Benchmark.tsv", sep="\t", index=False)

# ── Table S3: evidence-tiered gold standard ──
# Machine-readable mirror of supplementary.tex Table S3. Tier membership MUST
# match the TIER_A/TIER_B/TIER_C/FUNCTIONAL_ANALOGS constants in
# compute_headline_metrics.py — audit_manuscript_numbers.py enforces this.
_s3_rows = [
    # tier, driver, paralog, assay, model, direction, dual, indep, direct_sl, inclusion, ref
    ("A", "AKT1", "AKT2", "Combinatorial CRISPR digenic KO", "Cancer cell lines",
     "AKT1->AKT2", "Yes", "Yes", "Yes", "Primary", "Najm 2018"),
    ("A", "CDK4", "CDK6", "Digenic KO (pgPEN library)", "Cancer cell lines",
     "CDK4->CDK6", "Yes", "Yes", "Yes", "Primary", "Parrish 2021"),
    ("A", "MAP2K1", "MAP2K2", "Digenic KO (pgPEN library)", "Cancer cell lines",
     "MAP2K1->MAP2K2", "Yes", "Yes", "Yes", "Primary", "Parrish 2021"),
    ("B", "SMARCA4", "SMARCA2", "CRISPR KO conditioned on natural SMARCA4 mutation",
     "SMARCA4-mutant cancer lines", "SMARCA4->SMARCA2", "No", "Yes", "Conditional",
     "Primary", "Hoffman 2014"),
    ("B", "ARID1A", "ARID1B", "shRNA knockdown conditioned on natural ARID1A mutation",
     "ARID1A-mutant cancer lines", "ARID1A->ARID1B", "No", "Yes", "Conditional",
     "Primary", "Helming 2014"),
    ("C", "EP300", "CREBBP", "p300 degradation / CRISPR in CREBBP-deficient lines",
     "CREBBP-mutant cancer lines", "Reciprocal only (CREBBP->EP300)", "No", "Yes",
     "Reciprocal", "Secondary", "Ogiwara 2016; Nie 2021"),
    ("C", "PIK3CA", "PIK3CB", "shRNA / PI3K inhibitor in PTEN-deficient lines",
     "PTEN-deficient cancer lines", "PTEN->PIK3CB only", "No", "Yes",
     "No (other driver)", "Secondary", "Wee 2008"),
    ("C", "CCNE1", "CCNE2", "Mouse developmental double knockout", "Mouse embryo",
     "CCNE1<->CCNE2 redundancy", "Yes (mouse)", "Yes", "No (redundancy)",
     "Secondary", "Geng 2003"),
    ("C", "FBXW7", "FBXW2", "DepMap computational analysis", "700+ cell lines",
     "FBXW7->FBXW2", "No", "No", "No (computational)", "Secondary", "DepMap 26Q1 release"),
    ("C", "PPP2R1A", "PPP2R1B", "DepMap computational analysis", "700+ cell lines",
     "PPP2R1A->PPP2R1B", "No", "No", "No (computational)", "Secondary", "DepMap 26Q1 release"),
    ("Comparator", "BRCA1", "BRCA2", "PARP inhibition / genetic screen (functional analogs)",
     "BRCA-mutant cancer lines", "BRCA1/2->PARP axis", "No", "Yes",
     "Functional analog", "Comparator", "Bryant 2005"),
    ("Comparator", "STK11", "SIK1", "Mouse genetics, LKB1-SIK axis (partial homolog)",
     "Mouse NSCLC models", "STK11->SIK1/3 axis", "No", "Yes",
     "Pathway axis", "Comparator", "Hollstein 2019"),
]
s3 = pd.DataFrame(_s3_rows, columns=[
    "Tier", "Driver", "Paralog", "Assay", "Model", "Validated_Direction",
    "Dual_Perturbation", "DepMap_Independent", "Direct_SL", "Inclusion", "Key_Ref"])
s3.to_csv(OUT / "TableS3_GoldStandard.tsv", sep="\t", index=False)

# ── Table S1: Cell line counts ──
cl_counts = summary[["cancer","n_lines"]].copy()
cl_counts.columns = ["Cancer Type","Cell Lines (with data)"]
cl_counts.to_csv(OUT / "TableS1_CellLineCounts.tsv", sep="\t", index=False)

# ── Table S11: Full regression model table (confounder controls) ──
# Machine-readable mirror of output/regression_table_full.csv written by
# regression_controls.py (base / CNV- / expression- / lineage-adjusted models
# with beta, HC3 robust SE, 95% CI, nominal p and BH q per pair).
_reg_path = OUT.parent / "regression_table_full.csv"
if not _reg_path.exists():
    raise SystemExit("output/regression_table_full.csv not found — "
                     "run regression_controls.py first")
_regtab = pd.read_csv(_reg_path)
_regtab.to_csv(OUT / "TableS11_RegressionModels.tsv", sep="\t", index=False)
print(f"  TableS11_RegressionModels.tsv: {len(_regtab)} rows")

print(f"Tables saved to {OUT}")
print(f"  Table1_DeNovoCandidates.tsv: {len(cand) if all_candidates else 0} rows")
print(f"  Table1_Benchmark.tsv: {len(bench)} rows")
print(f"  Table2_Benchmark.tsv: {len(bench)} rows")
print(f"  TableS1_CellLineCounts.tsv: {len(cl_counts)} rows")
print(f"  TableS3_GoldStandard.tsv: {len(s3)} rows")
