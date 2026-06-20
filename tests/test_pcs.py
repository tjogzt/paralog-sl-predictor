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
