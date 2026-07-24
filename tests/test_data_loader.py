"""Unit tests for data_loader.py — DepMap data loading and preprocessing."""

import numpy as np
import pandas as pd
import pytest
from pathlib import Path
import sys
import tempfile
import os

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from data_loader import (
    _parse_gene_column,
    build_mutation_matrix,
    filter_gynecological_cell_lines,
    classify_cancer_type,
)


class TestParseGeneColumn:
    def test_parse_with_gene_id(self):
        assert _parse_gene_column("BRCA1 (672)") == "BRCA1"

    def test_parse_plain_gene(self):
        assert _parse_gene_column("TP53") == "TP53"

    def test_parse_empty(self):
        assert _parse_gene_column("") == ""

    def test_parse_gene_with_spaces(self):
        assert _parse_gene_column("ARID1A (12345)") == "ARID1A"

    def test_parse_depmap_id_column(self):
        assert _parse_gene_column("DepMap_ID") == "DepMap_ID"


class TestBuildMutationMatrix:
    def test_basic_matrix(self):
        mut_df = pd.DataFrame({
            "DepMap_ID": ["CL1", "CL1", "CL2", "CL4"],
            "Gene": ["TP53", "BRCA1", "TP53", "BRCA1"],
        })
        mat = build_mutation_matrix(
            mut_df,
            cell_lines=["CL1", "CL2", "CL3"],
            genes=["TP53", "BRCA1"]
        )
        assert mat.shape == (3, 2)
        assert mat.loc["CL1", "TP53"] == 1
        assert mat.loc["CL1", "BRCA1"] == 1
        assert mat.loc["CL2", "TP53"] == 1
        assert mat.loc["CL2", "BRCA1"] == 0
        assert mat.loc["CL3", "TP53"] == 0
        assert mat.loc["CL3", "BRCA1"] == 0

    def test_empty_mutations(self):
        mut_df = pd.DataFrame(columns=["DepMap_ID", "Gene"])
        mat = build_mutation_matrix(
            mut_df,
            cell_lines=["CL1", "CL2"],
            genes=["TP53"]
        )
        assert mat.shape == (2, 1)
        assert (mat == 0).all().all()

    def test_no_matching_cell_lines(self):
        mut_df = pd.DataFrame({
            "DepMap_ID": ["CL5", "CL6"],
            "Gene": ["TP53", "BRCA1"],
        })
        mat = build_mutation_matrix(
            mut_df,
            cell_lines=["CL1", "CL2"],
            genes=["TP53", "BRCA1"]
        )
        assert mat.shape == (2, 2)
        assert (mat == 0).all().all()

    def test_binary_output(self):
        """Multiple mutations in same gene/cell should still be 1."""
        mut_df = pd.DataFrame({
            "DepMap_ID": ["CL1", "CL1", "CL1"],
            "Gene": ["TP53", "TP53", "TP53"],
        })
        mat = build_mutation_matrix(
            mut_df,
            cell_lines=["CL1", "CL2"],
            genes=["TP53"]
        )
        assert mat.loc["CL1", "TP53"] == 1


class TestFilterGynecologicalCellLines:
    def test_filter_ovarian(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1", "CL2", "CL3", "CL4"],
            "OncotreePrimaryDisease": [
                "Ovarian Cancer", "Lung Adenocarcinoma",
                "Ovarian Serous Cystadenocarcinoma", "Breast Cancer"
            ],
        })
        result = filter_gynecological_cell_lines(models, "Ovarian")
        assert len(result) == 2
        assert "CL1" in result["DepMap_ID"].values
        assert "CL3" in result["DepMap_ID"].values

    def test_filter_endometrial(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1", "CL2"],
            "OncotreePrimaryDisease": [
                "Endometrial Carcinoma", "Uterine Carcinosarcoma"
            ],
        })
        result = filter_gynecological_cell_lines(models, "Endometrial")
        assert len(result) == 2

    def test_filter_no_cancer_type_all_gyn(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1", "CL2", "CL3"],
            "OncotreePrimaryDisease": [
                "Ovarian Cancer", "Endometrial Carcinoma", "Cervical Cancer"
            ],
        })
        result = filter_gynecological_cell_lines(models, None)
        assert len(result) == 3

    def test_filter_no_match(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1"],
            "OncotreePrimaryDisease": ["Lung Adenocarcinoma"],
        })
        result = filter_gynecological_cell_lines(models, "Ovarian")
        assert len(result) == 0


class TestClassifyCancerType:
    def test_ovarian_match(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1"],
            "OncotreePrimaryDisease": ["Ovarian Serous Cystadenocarcinoma"],
        })
        result = classify_cancer_type(models)
        assert result.get("CL1") == "Ovarian"

    def test_endometrial_match(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1"],
            "OncotreePrimaryDisease": ["Endometrial Endometrioid Adenocarcinoma"],
        })
        result = classify_cancer_type(models)
        assert result.get("CL1") == "Endometrial"

    def test_cervical_match(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1"],
            "OncotreePrimaryDisease": ["Cervical Squamous Cell Carcinoma"],
        })
        result = classify_cancer_type(models)
        assert result.get("CL1") == "Cervical"

    def test_multiple_cell_lines(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1", "CL2", "CL3"],
            "OncotreePrimaryDisease": [
                "Ovarian Cancer", "Lung Adenocarcinoma",
                "Endometrial Carcinoma"
            ],
        })
        result = classify_cancer_type(models)
        assert result.get("CL1") == "Ovarian"
        assert result.get("CL2") in ("Other", "Lung")
        assert result.get("CL3") == "Endometrial"

    def test_case_insensitive(self):
        models = pd.DataFrame({
            "DepMap_ID": ["CL1"],
            "OncotreePrimaryDisease": ["OVARIAN CANCER"],
        })
        result = classify_cancer_type(models)
        assert result.get("CL1") == "Ovarian"
