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
             "necessity","dependency_dd","composite_score","novelty","q_value",
             "mutation_frequency","n_mut","n_wt","is_known_paralog_sl"]
_ts2 = (_ts2.sort_values(["is_known_paralog_sl","composite_score"],
                         ascending=[False,False])[_ts2_cols])
_ts2.to_csv(OUT / "TableS2_FullResults.tsv", sep="\t", index=False)
print(f"  TableS2_FullResults.tsv: {len(_ts2)} rows")

summary = pd.read_csv(OUT.parent / "solid_tumor_summary.csv")
summary = summary.dropna(subset=["dd_auroc"]).sort_values("dd_auroc", ascending=False)

# ── Table S4 (NEW): Pan-cancer DD AUROC table ──
t = summary[["cancer","n_lines","n_pairs","n_known","dd_auroc"]].copy()
t.columns = ["Cancer Type","Cell Lines","SL Pairs","Known SL","DD AUROC"]
t["DD AUROC"] = t["DD AUROC"].apply(lambda x: f"{x:.3f}")
t.to_csv(OUT / "TableS4_PancancerSummary.tsv", sep="\t", index=False)

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

# ── Table S1: Cell line counts ──
cl_counts = summary[["cancer","n_lines"]].copy()
cl_counts.columns = ["Cancer Type","Cell Lines (with data)"]
cl_counts.to_csv(OUT / "TableS1_CellLineCounts.tsv", sep="\t", index=False)

print(f"Tables saved to {OUT}")
print(f"  Table1_DeNovoCandidates.tsv: {len(cand) if all_candidates else 0} rows")
print(f"  Table2_Benchmark.tsv: {len(bench)} rows")
print(f"  TableS1_CellLineCounts.tsv: {len(cl_counts)} rows")
print(f"  TableS4_PancancerSummary.tsv: {len(t)} rows")
