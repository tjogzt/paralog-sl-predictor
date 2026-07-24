"""Unit tests for pcs.py — ParalogCompensationScore and run_full_analysis."""

import numpy as np
import pandas as pd
import pytest
from pathlib import Path
import sys

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pcs import ParalogCompensationScore, run_full_analysis, _normalize_series


@pytest.fixture
def small_dependency():
    """6 cell lines x 5 genes, with BRCA1 dependency pattern."""
    np.random.seed(42)
    cl = [f"CL{i}" for i in range(1, 7)]
    df = pd.DataFrame(np.random.default_rng(42).normal(-0.5, 0.3, (6, 5)),
                      index=cl,
                      columns=["BRCA1", "BRCA2", "TP53", "ARID1A", "ARID1B"])
    # Make BRCA2 more essential in BRCA1-mutant lines (CL1-3)
    df.loc[["CL1", "CL2", "CL3"], "BRCA2"] = np.array([-0.9, -0.85, -0.95])
    df.loc[["CL4", "CL5", "CL6"], "BRCA2"] = np.array([-0.2, -0.15, -0.1])
    return df


@pytest.fixture
def small_expression():
    """6 cell lines x 5 genes."""
    np.random.seed(42)
    cl = [f"CL{i}" for i in range(1, 7)]
    return pd.DataFrame(np.random.default_rng(43).uniform(0, 5, (6, 5)),
                        index=cl,
                        columns=["BRCA1", "BRCA2", "TP53", "ARID1A", "ARID1B"])


@pytest.fixture
def small_models():
    """Cell line annotations."""
    return pd.DataFrame({
        "DepMap_ID": [f"CL{i}" for i in range(1, 7)],
        "OncotreePrimaryDisease": ["Ovarian Cancer", "Ovarian Cancer",
                                     "Ovarian Cancer", "Ovarian Cancer",
                                     "Ovarian Cancer", "Ovarian Cancer"],
        "StrippedCellLineName": [f"CL{i}_name" for i in range(1, 7)],
    })


@pytest.fixture
def small_mutations():
    """BRCA1-mutant in CL1-3."""
    return pd.DataFrame({
        "DepMap_ID": ["CL1", "CL2", "CL3"],
        "Gene": ["BRCA1", "BRCA1", "BRCA1"],
        "VariantInfo": ["frameshift_variant", "stop_gained", "missense_variant"],
    })


@pytest.fixture
def small_paralogs():
    """Paralog pairs: BRCA1-BRCA2, ARID1A-ARID1B."""
    return pd.DataFrame({
        "gene_A": ["BRCA1", "BRCA1", "ARID1A", "ARID1A"],
        "gene_B": ["BRCA2", "TP53", "ARID1B", "TP53"],
        "homology_type": ["paralog_known"] * 4,
        "identity_pct": [45.0, 10.0, 40.0, 8.0],
    })


class TestNormalizeSeries:
    def test_normalize_zero_range(self):
        s = pd.Series([0.5, 0.5, 0.5])
        result = _normalize_series(s)
        assert (result == 0.0).all()

    def test_normalize_typical(self):
        s = pd.Series([0, 5, 10])
        result = _normalize_series(s)
        assert result.iloc[0] == 0.0
        assert result.iloc[2] == 1.0


class TestParalogCompensationScore:
    def test_init_and_necessity(self, small_dependency, small_expression,
                                 small_models, small_mutations, small_paralogs):
        pcs = ParalogCompensationScore(
            small_dependency, small_expression, small_models,
            small_mutations, small_paralogs
        )
        assert "BRCA1" in pcs._necessity
        assert isinstance(pcs._necessity["BRCA1"], float)

    def test_compute_pcs_for_driver(self, small_dependency, small_expression,
                                     small_models, small_mutations, small_paralogs):
        pcs = ParalogCompensationScore(
            small_dependency, small_expression, small_models,
            small_mutations, small_paralogs
        )
        cell_lines = [f"CL{i}" for i in range(1, 7)]
        result = pcs.compute_pcs_for_driver("BRCA1", cell_lines, "Ovarian")

        assert isinstance(result, pd.DataFrame)
        assert not result.empty
        assert "driver_gene" in result.columns
        assert "paralog_gene" in result.columns
        assert "pcs" in result.columns
        assert "dependency_dd" in result.columns
        # BRCA2 is a paralog of BRCA1
        assert "BRCA2" in result["paralog_gene"].values
        # DD should be positive (BRCA2 more essential in BRCA1-mutant lines)
        brca2_row = result[result["paralog_gene"] == "BRCA2"]
        assert len(brca2_row) == 1

    def test_compute_pcs_for_driver_insufficient_samples(
            self, small_dependency, small_expression,
            small_models, small_mutations, small_paralogs):
        pcs = ParalogCompensationScore(
            small_dependency, small_expression, small_models,
            small_mutations, small_paralogs
        )
        # Only 1 cell line, below MIN_MUT_SAMPLES=3
        result = pcs.compute_pcs_for_driver("BRCA1", ["CL1"], "Ovarian")
        assert result.empty

    def test_compute_pcs_missing_driver(self, small_dependency, small_expression,
                                         small_models, small_mutations, small_paralogs):
        pcs = ParalogCompensationScore(
            small_dependency, small_expression, small_models,
            small_mutations, small_paralogs
        )
        cell_lines = [f"CL{i}" for i in range(1, 7)]
        result = pcs.compute_pcs_for_driver("NONEXISTENT", cell_lines, "Ovarian")
        assert isinstance(result, pd.DataFrame)
        assert result.empty


class TestRunFullAnalysis:
    def test_run_full_analysis_return_type(self, small_dependency,
                                            small_expression,
                                            small_models,
                                            small_mutations,
                                            small_paralogs):
        result = run_full_analysis(
            small_dependency, small_expression, small_models,
            small_mutations, small_paralogs,
            cancer_types=["Ovarian"]
        )
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert "composite_score" in result.columns
            assert "is_known_paralog_sl" in result.columns
            assert "q_value" in result.columns

    def test_run_full_analysis_with_no_results(self, small_dependency,
                                                small_expression,
                                                small_models,
                                                small_paralogs):
        """Empty mutations should produce empty results."""
        empty_mut = pd.DataFrame(columns=["DepMap_ID", "Gene"])
        result = run_full_analysis(
            small_dependency, small_expression, small_models,
            empty_mut, small_paralogs,
            cancer_types=["Ovarian"]
        )
        assert result.empty


class TestDDSignConventionAndBH:
    """Numerical-correctness tests for the DD sign convention (manuscript
    Eq. 1: DD = mean WT - mean MUT) and the BH q-value implementation."""

    def test_dd_sign_matches_manuscript_eq1(self, small_dependency, small_expression,
                                            small_models, small_mutations, small_paralogs):
        """Mutant lines MORE dependent on the paralog (lower Chronos) must
        yield POSITIVE DD (compensation), per manuscript Eq. 1 and the
        paralogSL R package."""
        pcs = ParalogCompensationScore(
            small_dependency, small_expression, small_models,
            small_mutations, small_paralogs
        )
        result = pcs.compute_pcs_for_driver("BRCA1", [f"CL{i}" for i in range(1, 7)], "Ovarian")
        brca2 = result[result["paralog_gene"] == "BRCA2"].iloc[0]
        # Fixture: BRCA2 Chronos ~ -0.9 in MUT vs ~ -0.15 in WT (stronger
        # dependency in mutant) -> DD = mean(WT) - mean(MUT) ~ +0.75
        assert brca2["dependency_dd"] > 0
        assert brca2["dependency_dd"] == pytest.approx(0.75, abs=0.05)
        assert brca2["cohens_d"] > 0

    def test_dd_p_value_on_dependency(self, small_dependency, small_expression,
                                      small_models, small_mutations, small_paralogs):
        """dd_p_value is the Welch t-test on dependency scores (distinct
        statistic from expr_p_value, mirroring R compute_dd()$p_value)."""
        pcs = ParalogCompensationScore(
            small_dependency, small_expression, small_models,
            small_mutations, small_paralogs
        )
        result = pcs.compute_pcs_for_driver("BRCA1", [f"CL{i}" for i in range(1, 7)], "Ovarian")
        assert "dd_p_value" in result.columns
        brca2 = result[result["paralog_gene"] == "BRCA2"].iloc[0]
        assert brca2["dd_p_value"] < 0.01  # groups are strongly separated in fixture

    def test_bh_adjust_matches_r_reference(self):
        """_bh_adjust must equal R's p.adjust(method='BH') reference values
        (verified against R 4.5: see compute_headline_metrics.py docs)."""
        from pcs import _bh_adjust
        v = np.array([0.9, 0.001, 0.5, 0.02, 0.2, 0.01, 0.03, 0.05])
        r_ref = np.array([0.9, 0.008, 0.5714285714286, 0.0533333333333,
                          0.2666666666667, 0.04, 0.06, 0.08])
        assert np.allclose(_bh_adjust(v), r_ref, atol=1e-10)

    def test_bh_adjust_monotone_and_bounded(self):
        """For sorted p-values, q-values must be non-decreasing and within [0, 1]
        (the property the previous row-order implementation could violate)."""
        from pcs import _bh_adjust
        rng = np.random.default_rng(7)
        for _ in range(50):
            p = np.sort(rng.uniform(0, 1, size=rng.integers(2, 40)))
            q = _bh_adjust(p)
            assert np.all(np.diff(q) >= -1e-12)
            assert (q >= 0).all() and (q <= 1).all()

    def test_bh_adjust_nan_safe(self):
        """NaN p-values are treated as 1.0 (never significant) and cannot
        poison the correction of the remaining p-values."""
        from pcs import _bh_adjust
        q = _bh_adjust(np.array([0.01, np.nan, 0.04]))
        # Reference: R p.adjust(c(0.01, 1.0, 0.04), 'BH') = c(0.03, 1.0, 0.06)
        assert np.allclose(q, [0.03, 1.0, 0.06], atol=1e-12)
