#!/usr/bin/env python3
"""
make_driver_mutation_rules_table.py
====================================
Generate output/driver_mutation_rules.csv — the complete variant
classification table required by the Methods section:

For every gene in config.GENE_DRIVER_CLASS, report
  gene, driver_class (TSG/ONC), rule (LikelyLoF/Hotspot),
  n_lines_old   — cell lines carrying any non-silent variant
                  (pre-C7 permissive definition),
  n_lines_new   — cell lines carrying a rule-qualifying driver variant
                  (class-specific definition used throughout the pipeline),
  top VariantInfo terms among the non-silent variants of that gene.

The table is the audit trail for the C7 change: TSG driver status now
requires LikelyLoF, oncogene driver status requires Hotspot, both on the
default omics profile per model (DepMap 26Q1).
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from config import DEPMAP_FILES, GENE_DRIVER_CLASS, driver_mutation_rule  # noqa: E402
from data_loader import load_mutations  # noqa: E402

OUT = ROOT / "output" / "driver_mutation_rules.csv"


def main():
    mut = load_mutations(DEPMAP_FILES["mutations"])
    rows = []
    for gene in sorted(GENE_DRIVER_CLASS):
        cls = driver_mutation_rule(gene)
        rule = "LikelyLoF" if cls == "TSG" else "Hotspot"
        sub = mut[mut["Gene"] == gene]
        old_lines = sub["DepMap_ID"].nunique()
        new_lines = sub.loc[sub[rule] == True, "DepMap_ID"].nunique()
        terms = (
            sub["VariantInfo"].value_counts().head(5)
            .apply(lambda n: str(n)).to_dict()
        )
        top = "; ".join(f"{k} ({v})" for k, v in terms.items())
        rows.append({
            "gene": gene,
            "driver_class": cls,
            "rule": rule,
            "n_lines_old_nonsilent": old_lines,
            "n_lines_new_rule": new_lines,
            "pct_retained": round(100 * new_lines / old_lines, 1) if old_lines else 0.0,
            "top_variant_info": top,
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(df.to_string(index=False))
    print(f"\nWrote {OUT} ({len(df)} genes)")


if __name__ == "__main__":
    main()
