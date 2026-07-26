#!/usr/bin/env python3
"""Update test_pcs.py fixtures to the min>=5 framework (5 MUT / 5 WT lines)."""
from pathlib import Path

p = Path("tests/test_pcs.py")
t = p.read_text()

REPS = [
    # dependency fixture: 10 lines, BRCA2 strongly essential in CL1-5
    ('''    np.random.seed(42)
    cl = [f"CL{i}" for i in range(1, 7)]
    df = pd.DataFrame(np.random.default_rng(42).normal(-0.5, 0.3, (6, 5)),
                      index=cl,
                      columns=["BRCA1", "BRCA2", "TP53", "ARID1A", "ARID1B"])
    # Make BRCA2 more essential in BRCA1-mutant lines (CL1-3)
    df.loc[["CL1", "CL2", "CL3"], "BRCA2"] = np.array([-0.9, -0.85, -0.95])
    df.loc[["CL4", "CL5", "CL6"], "BRCA2"] = np.array([-0.2, -0.15, -0.1])
    return df''',
     '''    np.random.seed(42)
    cl = [f"CL{i}" for i in range(1, 11)]
    df = pd.DataFrame(np.random.default_rng(42).normal(-0.5, 0.3, (10, 5)),
                      index=cl,
                      columns=["BRCA1", "BRCA2", "TP53", "ARID1A", "ARID1B"])
    # Make BRCA2 more essential in BRCA1-mutant lines (CL1-5)
    df.loc[["CL1", "CL2", "CL3", "CL4", "CL5"], "BRCA2"] = np.array(
        [-0.9, -0.85, -0.95, -0.88, -0.92])
    df.loc[["CL6", "CL7", "CL8", "CL9", "CL10"], "BRCA2"] = np.array(
        [-0.2, -0.15, -0.1, -0.18, -0.12])
    return df'''),
    # expression fixture: 10 lines
    ('''    np.random.seed(42)
    cl = [f"CL{i}" for i in range(1, 7)]
    return pd.DataFrame(np.random.default_rng(43).uniform(0, 5, (6, 5)),
                        index=cl,
                        columns=["BRCA1", "BRCA2", "TP53", "ARID1A", "ARID1B"])''',
     '''    np.random.seed(42)
    cl = [f"CL{i}" for i in range(1, 11)]
    return pd.DataFrame(np.random.default_rng(43).uniform(0, 5, (10, 5)),
                        index=cl,
                        columns=["BRCA1", "BRCA2", "TP53", "ARID1A", "ARID1B"])'''),
    # models fixture: 10 rows
    ('''        "DepMap_ID": [f"CL{i}" for i in range(1, 7)],
        "OncotreePrimaryDisease": ["Ovarian Cancer", "Ovarian Cancer",
                                     "Ovarian Cancer", "Ovarian Cancer",
                                     "Ovarian Cancer", "Ovarian Cancer"],
        "StrippedCellLineName": [f"CL{i}_name" for i in range(1, 7)],''',
     '''        "DepMap_ID": [f"CL{i}" for i in range(1, 11)],
        "OncotreePrimaryDisease": ["Ovarian Cancer"] * 10,
        "StrippedCellLineName": [f"CL{i}_name" for i in range(1, 11)],'''),
    # mutations fixture: CL1-5 mutant (LikelyLoF-type variants for the TSG rule)
    ('''    """BRCA1-mutant in CL1-3."""
    return pd.DataFrame({
        "DepMap_ID": ["CL1", "CL2", "CL3"],
        "Gene": ["BRCA1", "BRCA1", "BRCA1"],
        "VariantInfo": ["frameshift_variant", "stop_gained", "missense_variant"],
    })''',
     '''    """BRCA1-mutant in CL1-5."""
    return pd.DataFrame({
        "DepMap_ID": ["CL1", "CL2", "CL3", "CL4", "CL5"],
        "Gene": ["BRCA1"] * 5,
        "VariantInfo": ["frameshift_variant", "stop_gained", "frameshift_variant",
                        "stop_gained", "frameshift_variant"],
    })'''),
    # compute_pcs_for_driver tests: 10 lines
    ('''        cell_lines = [f"CL{i}" for i in range(1, 7)]
        result = pcs.compute_pcs_for_driver("BRCA1", cell_lines, "Ovarian")

        assert isinstance(result, pd.DataFrame)''',
     '''        cell_lines = [f"CL{i}" for i in range(1, 11)]
        result = pcs.compute_pcs_for_driver("BRCA1", cell_lines, "Ovarian")

        assert isinstance(result, pd.DataFrame)'''),
    ('''        # Only 1 cell line, below MIN_MUT_SAMPLES=3
        result = pcs.compute_pcs_for_driver("BRCA1", ["CL1"], "Ovarian")''',
     '''        # Only 1 cell line, below MIN_MUT_SAMPLES=5
        result = pcs.compute_pcs_for_driver("BRCA1", ["CL1"], "Ovarian")'''),
    ('''        cell_lines = [f"CL{i}" for i in range(1, 7)]
        result = pcs.compute_pcs_for_driver("NONEXISTENT", cell_lines, "Ovarian")''',
     '''        cell_lines = [f"CL{i}" for i in range(1, 11)]
        result = pcs.compute_pcs_for_driver("NONEXISTENT", cell_lines, "Ovarian")'''),
    ('''        result = pcs.compute_pcs_for_driver("BRCA1", [f"CL{i}" for i in range(1, 7)], "Ovarian")
        brca2 = result[result["paralog_gene"] == "BRCA2"].iloc[0]
        # Fixture: BRCA2 Chronos ~ -0.9 in MUT vs ~ -0.15 in WT (stronger''',
     '''        result = pcs.compute_pcs_for_driver("BRCA1", [f"CL{i}" for i in range(1, 11)], "Ovarian")
        brca2 = result[result["paralog_gene"] == "BRCA2"].iloc[0]
        # Fixture: BRCA2 Chronos ~ -0.9 in MUT vs ~ -0.15 in WT (stronger'''),
    ('''        result = pcs.compute_pcs_for_driver("BRCA1", [f"CL{i}" for i in range(1, 7)], "Ovarian")
        assert "dd_p_value" in result.columns''',
     '''        result = pcs.compute_pcs_for_driver("BRCA1", [f"CL{i}" for i in range(1, 11)], "Ovarian")
        assert "dd_p_value" in result.columns'''),
]

missing = [(t.count(old), old[:80]) for old, _ in REPS if t.count(old) != 1]
if missing:
    for n, s in missing:
        print(f"MATCH {n}: {s}")
    raise SystemExit("ABORT")
for old, new in REPS:
    t = t.replace(old, new)
p.write_text(t)
print(f"OK: {len(REPS)} fixture replacements")
