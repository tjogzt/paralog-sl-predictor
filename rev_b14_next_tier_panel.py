#!/usr/bin/env python3
"""
rev_b14_next_tier_panel.py  (Stage-4 revision, item B14)
=========================================================
Next-tier candidate panel: the top-10 NON-BENCHMARK driver x paralog x
lineage entries of the primary frame by composite score, with the evidence
columns requested by the reviewers (DD / DWS / selectivity / composite /
library membership), so that KMT2D->KMT2C is presented as one panel member
among several rather than a singled-out nomination.

Universe: output/tables/TableS2_FullResults.tsv (110 gyn3 entries) minus the
12 curated benchmark pairs (unordered match; comparators included in the
exclusion). Ranked by composite_score (shipped, pcs.py).

Columns:
  driver, paralog, lineage, signed_dd, hedges_g, pcs, necessity, composite,
  dws, selectivity, dws_classification (TableS5_DWS; NaN = not computed for
  that pair), in_harle_library (472 pairs), in_flister_library (35,108
  Table S3 pairs), in_synlethdb (local data/synlethdb_sl_pairs.csv),
  is_previous_nomination (KMT2D->KMT2C).

Outputs (output/revision_stage4/):
  b14_next_tier_panel.csv
  b14_next_tier_panel.md   (human-readable table)

Usage: python rev_b14_next_tier_panel.py   (run from repo root)
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output" / "revision_stage4"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT))
from config import KNOWN_PARALOG_SL  # noqa: E402

BENCH = {frozenset(p) for p in KNOWN_PARALOG_SL}
PREV_NOMINATION = frozenset({"KMT2D", "KMT2C"})


def main():
    print("=" * 72)
    print("  rev B14: next-tier candidate panel")
    print("=" * 72)

    df = pd.read_csv(ROOT / "output" / "tables" / "TableS2_FullResults.tsv", sep="\t")
    df["is_known_paralog_sl"] = df["is_known_paralog_sl"].astype(bool)
    non_bench = df[~df.apply(lambda r: frozenset((r["driver_gene"], r["paralog_gene"])) in BENCH,
                             axis=1)].copy()
    top = non_bench.sort_values("composite_score", ascending=False).head(10).copy()

    # DWS / selectivity (pair level)
    dws = pd.read_csv(ROOT / "output" / "tables" / "TableS5_DWS.tsv", sep="\t")
    dmap = {frozenset((r.driver, r.paralog)): r for r in dws.itertuples()}

    # Harle / Flister library membership
    h5 = pd.read_pickle(ROOT / "output" / "cache" / "harle_tableS5.pkl")
    harle_lib = {frozenset(str(p).split("|")) for p in h5["sorted_gene_pair"]}
    fl = pd.read_pickle(ROOT / "output" / "cache" / "flister_tableS3_sheet1.pkl")
    flister_lib = {frozenset(str(p).split("_")) for p in fl["label"]}

    # SynLethDB (local)
    sldb = pd.read_csv(ROOT / "data" / "synlethdb_sl_pairs.csv")
    sl_pairs = {frozenset((r.gene_A, r.gene_B)) for r in sldb.itertuples()}

    rows = []
    for r in top.itertuples():
        key = frozenset((r.driver_gene, r.paralog_gene))
        d = dmap.get(key)
        rows.append({
            "rank": len(rows) + 1,
            "driver": r.driver_gene, "paralog": r.paralog_gene,
            "lineage": r.cancer_type,
            "signed_dd": round(float(r.dependency_dd), 4),
            "hedges_g": round(float(r.hedges_g), 3),
            "pcs": round(float(r.pcs), 4),
            "necessity": round(float(r.necessity), 4),
            "composite_score": round(float(r.composite_score), 4),
            "dws": round(float(d.dws), 3) if d is not None else None,
            "selectivity": round(float(d.selectivity), 4) if d is not None else None,
            "dws_classification": d.classification if d is not None else "not computed",
            "in_harle_library": key in harle_lib,
            "in_flister_library": key in flister_lib,
            "in_synlethdb": key in sl_pairs,
            "is_previous_nomination": key == PREV_NOMINATION,
        })
    panel = pd.DataFrame(rows)
    panel.to_csv(OUT / "b14_next_tier_panel.csv", index=False)

    # human-readable markdown
    md = ["# B14 next-tier candidate panel (top-10 non-benchmark, primary frame)",
          "",
          "Universe: TableS2 gyn3 entries minus the 12 curated benchmark pairs; ranked by composite score.",
          "DWS/selectivity from TableS5 (pair level; 'not computed' = pair absent from the DWS table).",
          "",
          panel.to_markdown(index=False),
          "",
          f"Harle library = 472 screened pairs (Harle et al. 2025); Flister library = "
          f"35,108 pairs with lethality calls (Flister et al. 2025 Table S3); "
          f"SynLethDB membership checked against the LOCAL data/synlethdb_sl_pairs.csv, "
          f"which is the pipeline's minimal literature-curated fallback set "
          f"({len(sl_pairs)} unordered pairs, evidence='literature') -- NOT the full "
          f"SynLethDB database; membership in the full database is NOT COMPUTABLE "
          f"offline and is therefore annotated as 'not checked (full DB)'.",
          ""]
    (OUT / "b14_next_tier_panel.md").write_text("\n".join(md))

    print(panel[["rank", "driver", "paralog", "lineage", "signed_dd",
                "composite_score", "dws", "in_harle_library",
                "in_flister_library", "in_synlethdb"]].to_string(index=False))
    print(f"\n  wrote {OUT}/b14_next_tier_panel.csv + .md")


if __name__ == "__main__":
    main()
