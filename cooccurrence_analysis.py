"""
Mutational Co-occurrence Analysis for Key Paralog Pairs
========================================================
For each key driver–paralog pair, tests whether mutations in the two
genes co-occur or are mutually exclusive across the 1,208 DepMap 26Q1
cell lines. Mutation status uses the same gene-class-specific driver
rules as the main pipeline (TSG: LikelyLoF; oncogene: Hotspot; see
data_loader.build_mutation_matrix and manuscript Methods).

For each pair a 2x2 contingency table (both / A-only / B-only / neither)
is built and Fisher's exact odds ratio + p-value computed. Results are
written to output/cooccurrence_analysis.csv — the single source of truth
for Fig. 3d (R_fig3.R panel_d) and the manuscript co-occurrence statement.

No simulated, random, or hardcoded values: every number in the output is
computed from data/OmicsSomaticMutations.csv.
"""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from config import OUTPUT_DIR
from data_loader import (
    load_mutations,
    build_mutation_matrix,
    load_models,
    load_dependency,
)

# Key pairs shown in Fig. 3d (driver, paralog)
PAIRS = [
    ("ARID1A", "ARID1B"),
    ("PIK3CA", "PIK3CB"),
    ("BRCA1", "BRCA2"),
    ("EP300", "CREBBP"),
    ("SMARCA4", "SMARCA2"),
]


def main():
    print("=" * 70)
    print("  Mutational Co-occurrence Analysis (DepMap 26Q1)")
    print("=" * 70)

    # Cell-line universe: models with dependency data (same 1,208 as the
    # main pipeline)
    dep = load_dependency()
    cell_lines = list(dep.index)
    print(f"  Cell lines with dependency data: {len(cell_lines)}")

    genes = sorted({g for pair in PAIRS for g in pair})
    mut = load_mutations()
    mat = build_mutation_matrix(mut, cell_lines, genes,
                                apply_driver_rules=True)

    rows = []
    for a, b in PAIRS:
        both = int(((mat[a] == 1) & (mat[b] == 1)).sum())
        a_only = int(((mat[a] == 1) & (mat[b] == 0)).sum())
        b_only = int(((mat[a] == 0) & (mat[b] == 1)).sum())
        neither = int(((mat[a] == 0) & (mat[b] == 0)).sum())
        table = [[both, a_only], [b_only, neither]]
        or_val, p_val = fisher_exact(table, alternative="two-sided")

        # Wald 95% CI on log(OR) with 0.5 Haldane correction if any cell
        # is zero (avoids division by zero)
        bb, aa, cc, dd = both, a_only, b_only, neither
        if min(bb, aa, cc, dd) == 0:
            bb, aa, cc, dd = bb + 0.5, aa + 0.5, cc + 0.5, dd + 0.5
        log_or = np.log((bb * dd) / (aa * cc))
        se = np.sqrt(1 / bb + 1 / aa + 1 / cc + 1 / dd)
        ci_low = float(np.exp(log_or - 1.96 * se))
        ci_high = float(np.exp(log_or + 1.96 * se))

        rows.append({
            "driver": a,
            "paralog": b,
            "pair": f"{a}/{b}",
            "odds_ratio": float(or_val),
            "ci_low": ci_low,
            "ci_high": ci_high,
            "p_value": float(p_val),
            "n_both": both,
            "n_driver_only": a_only,
            "n_paralog_only": b_only,
            "n_neither": neither,
            "n_total": len(cell_lines),
        })
        print(f"  {a:9s}/{b:9s}  OR={or_val:6.3f} [{ci_low:.3f}-{ci_high:.3f}]  "
              f"p={p_val:.2e}  (both={both}, {a}-only={a_only}, "
              f"{b}-only={b_only}, neither={neither})")

    out = pd.DataFrame(rows)
    out_path = OUTPUT_DIR / "cooccurrence_analysis.csv"
    out.to_csv(out_path, index=False)
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()
