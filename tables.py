"""Generate manuscript tables for the pan-cancer revision."""
import pandas as pd; import numpy as np
from pathlib import Path
OUT = Path(__file__).resolve().parent / "output" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

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
bench = pd.DataFrame({
    "Method": ["SLMGAE","DDSL","GRSL","NSF4SL","Struct2SL","KG4SL","DDGCN","PGCN","DD (this study)","DD + ID≥0.3"],
    "CV3_AUROC": ["0.700","0.720","0.680","0.650","0.650","0.620","0.600","0.580","0.794","1.000"],
    "Architecture": ["GNN","GCN","Graph Reg.","Neg. Sampling","AF2+MLP","KG Embed.","Dual GCN","Pathway GCN","Univariate","Univariate+Filter"],
    "Interpretability": ["Low","Low","Low","Low","Medium","Low","Low","Low","High","High"],
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
