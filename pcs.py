"""
Paralog SL Predictor — PCS (Paralog Compensation Score) Core
=============================================================
Computes the Paralog Compensation Score for gene pairs across
gynecological cancer cell lines.

PCS_AB = ΔExpression(B, MUT vs WT) × Necessity(B)

Where:
  ΔExpression = mean expression of paralog B in A-mutant cell lines
               minus mean expression of B in A-wildtype cell lines
  Necessity(B) = -mean(CERES/Chronos dependency of B)
                 (negative CERES = essential, so -CERES = positive necessity)

DD (Delta Dependency, manuscript Eq. 1):
  DD = mean(Chronos of paralog B | A-WT) − mean(Chronos of B | A-MUT)
  Positive DD = stronger paralog dependency in driver-mutant lines,
  consistent with paralog compensation. All AUROC/scoring uses |DD|
  (see manuscript), so the sign convention is presentation-only here;
  it matches manuscript Eq. 1 and the paralogSL R package (compute_dd).
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import false_discovery_control
from typing import Dict, List, Tuple, Optional
from tqdm import tqdm

from config import (
    DRIVER_GENES, GYN_CANCER_TYPES, KNOWN_PARALOG_SL,
    MIN_MUT_SAMPLES, MIN_WT_SAMPLES,
    MIN_DELTA_EXPR, PCS_THRESHOLD, SIGNIFICANCE_ALPHA,
)
from data_loader import (
    build_mutation_matrix, filter_gynecological_cell_lines,
    classify_cancer_type,
)


def _normalize_series(s: pd.Series) -> pd.Series:
    """Min-max normalize to [0, 1]. Returns 0 if all values equal."""
    if s.max() == s.min():
        return pd.Series(0.0, index=s.index)
    return (s - s.min()) / (s.max() - s.min())


def _bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted q-values (step-up with ranked monotonicity).

    Thin wrapper around scipy.stats.false_discovery_control. NaN p-values are
    treated as 1.0 (never significant) so they cannot corrupt the correction.
    Replaces a previous hand-rolled implementation that enforced monotonicity
    along the dataframe row order instead of the p-value rank order.
    """
    p = np.asarray(p_values, dtype=float)
    p = np.where(np.isnan(p), 1.0, p).clip(0.0, 1.0)
    if p.size <= 1:
        return p.copy()
    return false_discovery_control(p, method="bh")


class ParalogCompensationScore:
    """
    Compute PCS for all driver-paralog pairs in gynecological cancers.

    Parameters
    ----------
    dependency : pd.DataFrame
        CRISPR dependency matrix (rows=cell lines, cols=genes).
    expression : pd.DataFrame
        Gene expression matrix (rows=cell lines, cols=genes).
    models : pd.DataFrame
        Cell line annotations with OncotreePrimaryDisease.
    mutations : pd.DataFrame
        Damaging mutation calls.
    paralogs : pd.DataFrame
        Paralog pairs with columns [gene_A, gene_B].
    """

    def __init__(self,
                 dependency: pd.DataFrame,
                 expression: pd.DataFrame,
                 models: pd.DataFrame,
                 mutations: pd.DataFrame,
                 paralogs: pd.DataFrame):
        self.dep = dependency
        self.expr = expression
        self.models = models
        self.mutations = mutations
        self.paralogs = paralogs

        # Pre-compute: global mean necessity for each gene
        self._necessity = {}
        self._compute_necessity()

    def _compute_necessity(self):
        """Compute global mean necessity (-mean CERES) for all genes."""
        for gene in self.dep.columns:
            ceres = self.dep[gene].dropna()
            if len(ceres) > 0:
                # CERES: negative = essential. necessity = -CERES (positive = more essential)
                self._necessity[gene] = -ceres.mean()
            else:
                self._necessity[gene] = 0.0

    def get_necessity(self, gene: str) -> float:
        """Get global necessity for a gene. Returns 0 if unknown."""
        return self._necessity.get(gene, 0.0)

    def _get_cell_lines(self, cancer_type: Optional[str] = None) -> List[str]:
        """Get gynecological cell line IDs, optionally filtered by cancer type."""
        gyn = filter_gynecological_cell_lines(self.models, cancer_type)
        ids = gyn["DepMap_ID"].tolist()
        # Intersect with available data
        return [cl for cl in ids if cl in self.dep.index and cl in self.expr.index]

    def compute_pcs_for_driver(self,
                                driver_gene: str,
                                cell_lines: List[str],
                                cancer_label: str = "Gyn",
                                mut_matrix: Optional[pd.DataFrame] = None
                                ) -> pd.DataFrame:
        """
        Compute PCS for all paralogs of a given driver gene.

        Parameters
        ----------
        driver_gene : str
            The potentially mutated driver gene.
        cell_lines : list
            List of cell line IDs to analyze.
        cancer_label : str
            Label for the cancer type (Ovarian/Endometrial/Cervical/Gyn).
        mut_matrix : pd.DataFrame, optional
            Precomputed binary mutation matrix (rows = cell lines, must
            include a column named `driver_gene`). When provided, the internal
            build_mutation_matrix call is skipped and rows are sliced to
            `cell_lines` — used by pancancer.py to avoid rebuilding the
            matrix once per lineage x driver combination.

        Returns
        -------
        pd.DataFrame with columns:
            driver_gene, paralog_gene, pcs, delta_expression,
            necessity, p_value, cohens_d,
            n_mut, n_wt, cancer_type,
            mutation_frequency, pass_threshold
        """
        # Get paralogs for this driver
        driver_paralogs = self.paralogs[self.paralogs["gene_A"] == driver_gene]["gene_B"].unique()
        if len(driver_paralogs) == 0:
            return pd.DataFrame()

        # Filter paralogs that exist in both dependency and expression matrices
        valid_paralogs = [p for p in driver_paralogs
                          if p in self.dep.columns and p in self.expr.columns]
        if not valid_paralogs:
            return pd.DataFrame()

        # Build (or reuse) the mutation matrix for this driver only
        if mut_matrix is None:
            mut_matrix = build_mutation_matrix(self.mutations, cell_lines, [driver_gene])
        else:
            mut_matrix = mut_matrix.reindex(index=cell_lines, fill_value=0)
        if driver_gene not in mut_matrix.columns:
            return pd.DataFrame()

        mut_cl = mut_matrix[mut_matrix[driver_gene] == 1].index.tolist()
        wt_cl = mut_matrix[mut_matrix[driver_gene] == 0].index.tolist()

        # Check sample sizes
        mut_cl = [c for c in mut_cl if c in self.dep.index and c in self.expr.index]
        wt_cl = [c for c in wt_cl if c in self.dep.index and c in self.expr.index]

        if len(mut_cl) < MIN_MUT_SAMPLES or len(wt_cl) < MIN_WT_SAMPLES:
            return pd.DataFrame()

        # Mutation frequency
        mut_freq = len(mut_cl) / (len(mut_cl) + len(wt_cl))

        results = []
        for paralog in valid_paralogs:
            # ── 1. Expression change ──
            expr_mut = self.expr.loc[mut_cl, paralog].mean()
            expr_wt = self.expr.loc[wt_cl, paralog].mean()
            delta_expr = expr_mut - expr_wt

            # Welch's t-test on expression change
            expr_mut_vals = self.expr.loc[mut_cl, paralog].dropna()
            expr_wt_vals = self.expr.loc[wt_cl, paralog].dropna()

            if len(expr_mut_vals) < MIN_MUT_SAMPLES or len(expr_wt_vals) < MIN_WT_SAMPLES:
                t_stat, p_val = 0.0, 1.0
            else:
                t_stat, p_val = stats.ttest_ind(expr_mut_vals, expr_wt_vals, equal_var=False)

            # ── 2. Dependency change: Delta Dependency (manuscript Eq. 1) ──
            dep_mut_vals = self.dep.loc[mut_cl, paralog].dropna()
            dep_wt_vals = self.dep.loc[wt_cl, paralog].dropna()

            if len(dep_mut_vals) >= MIN_MUT_SAMPLES and len(dep_wt_vals) >= MIN_WT_SAMPLES:
                # DD = mean(Chronos | WT) − mean(Chronos | MUT).
                # Positive DD = stronger paralog dependency in driver-mutant
                # lines (paralog compensation); matches manuscript Eq. 1 and
                # paralogSL::compute_dd. Validation/AUROC ranks by |DD|.
                dd = dep_wt_vals.mean() - dep_mut_vals.mean()
                # Cohen's d on dependency (positive = stronger in mutant)
                pooled_std = np.sqrt((dep_mut_vals.var() + dep_wt_vals.var()) / 2)
                cohens_d = dd / pooled_std if pooled_std > 0 else 0.0
                # Welch's t-test on dependency scores; mirrors the p_value
                # returned by paralogSL::compute_dd (kept distinct from the
                # expression-based expr_p_value above).
                _, dd_p_val = stats.ttest_ind(dep_mut_vals, dep_wt_vals, equal_var=False)
            else:
                dd = 0.0
                cohens_d = 0.0
                dd_p_val = 1.0

            # ── 3. Paralog compensation score ──
            necessity = self.get_necessity(paralog)
            # PCS = expression upregulation × functional necessity
            # Only positive delta_expr contributes (upregulation in mutant context)
            pcs = max(delta_expr, 0.0) * max(necessity, 0.0)

            results.append({
                "driver_gene": driver_gene,
                "paralog_gene": paralog,
                "pcs": pcs,
                "delta_expression": delta_expr,
                "necessity": necessity,
                "dependency_dd": dd,
                "cohens_d": cohens_d,
                "dd_p_value": dd_p_val,
                "expr_p_value": p_val,
                "expr_t_stat": t_stat,
                "n_mut": len(mut_cl),
                "n_wt": len(wt_cl),
                "cancer_type": cancer_label,
                "mutation_frequency": mut_freq,
            })

        if not results:
            return pd.DataFrame()

        df = pd.DataFrame(results)

        # ── 4. Benjamini-Hochberg correction ──
        df["q_value"] = _bh_adjust(df["expr_p_value"].values)

        # ── 5. Threshold filtering ──
        df["pass_threshold"] = (
            (df["pcs"] >= PCS_THRESHOLD) &
            (df["delta_expression"] >= MIN_DELTA_EXPR) &
            (df["q_value"] < SIGNIFICANCE_ALPHA)
        )

        # ── 6. Normalized ranks ──
        df["pcs_rank"] = df["pcs"].rank(ascending=False)
        df["pcs_norm"] = _normalize_series(df["pcs"])
        df["dd_norm"] = _normalize_series(df["dependency_dd"].abs())
        df["composite_score"] = (
            0.50 * df["pcs_norm"] +
            0.20 * df["dd_norm"] +
            0.15 * _normalize_series(-np.log10(df["q_value"] + 1e-10)) +
            0.15 * _normalize_series(df["mutation_frequency"])
        )

        return df.sort_values("composite_score", ascending=False)


def run_full_analysis(dependency: pd.DataFrame,
                      expression: pd.DataFrame,
                      models: pd.DataFrame,
                      mutations: pd.DataFrame,
                      paralogs: pd.DataFrame,
                      cancer_types: Optional[List[str]] = None
                      ) -> pd.DataFrame:
    """
    Run the complete PCS analysis across all specified cancer types.

    Parameters
    ----------
    cancer_types : list of str, optional
        List of cancer type keys: ['Ovarian', 'Endometrial', 'Cervical'].
        Default: all three.

    Returns
    -------
    pd.DataFrame with all PCS results, sorted by composite_score.
    """
    if cancer_types is None:
        cancer_types = list(GYN_CANCER_TYPES.keys())

    pcs = ParalogCompensationScore(
        dependency, expression, models, mutations, paralogs
    )

    all_results = []
    cell_lines_all = pcs._get_cell_lines(cancer_type=None)
    cell_line_types = classify_cancer_type(models)

    for ctype in cancer_types:
        print(f"\n{'='*50}")
        print(f"Analyzing: {ctype} Cancer")
        print(f"{'='*50}")

        # Get cell lines for this cancer type
        cl_subset = [c for c in cell_lines_all
                     if cell_line_types.get(c, "Other") == ctype]
        print(f"  Cell lines: {len(cl_subset)}")

        if len(cl_subset) < (MIN_MUT_SAMPLES + MIN_WT_SAMPLES):
            print(f"  Skipping: insufficient cell lines")
            continue

        driver_genes = DRIVER_GENES.get(ctype, [])
        print(f"  Driver genes: {len(driver_genes)}")

        for driver in tqdm(driver_genes, desc=f"  {ctype} drivers"):
            if driver not in pcs.dep.columns:
                continue

            df = pcs.compute_pcs_for_driver(driver, cl_subset, cancer_label=ctype)
            if not df.empty:
                all_results.append(df)

    if not all_results:
        print("\nNo results generated. Check data availability.")
        return pd.DataFrame()

    final = pd.concat(all_results, ignore_index=True)

    # Global normalization for composite_score (fixes per-cancer collapse)
    if len(final) > 1:
        final["pcs_norm"] = _normalize_series(final["pcs"])
        final["dd_norm"] = _normalize_series(final["dependency_dd"].abs())
        final["q_penalty"] = _normalize_series(-np.log10(final["q_value"].clip(lower=1e-10)))
        final["freq_norm"] = _normalize_series(final["mutation_frequency"])
        final["composite_score"] = (
            0.50 * final["pcs_norm"] +
            0.20 * final["dd_norm"] +
            0.15 * final["q_penalty"] +
            0.15 * final["freq_norm"]
        )

    # Add known-SL annotation
    final["is_known_paralog_sl"] = final.apply(
        lambda r: (r["driver_gene"], r["paralog_gene"]) in KNOWN_PARALOG_SL or
                  (r["paralog_gene"], r["driver_gene"]) in KNOWN_PARALOG_SL,
        axis=1
    )

    # Sort: known-positive pairs first, then by composite score
    final = final.sort_values(
        ["is_known_paralog_sl", "composite_score"],
        ascending=[False, False]
    )

    print(f"\n{'='*50}")
    print(f"Analysis complete.")
    print(f"  Total pairs analyzed: {len(final)}")
    print(f"  Passing threshold:    {final['pass_threshold'].sum()}")
    print(f"  Known paralog-SL:     {final['is_known_paralog_sl'].sum()}")
    print(f"{'='*50}")

    return final
