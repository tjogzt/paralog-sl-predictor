# Tests for paralogSL core functions (v2 - fixed matrix dimensions)
source("../../R/paralogSL.R", chdir = TRUE)
test_that("compute_dd returns correct structure", {
  vals <- c(-0.5, -0.3, -0.1, -0.8, -0.4, -0.2,
            -0.3, -0.1, 0.1, -0.5, -0.2, -0.3)
  dep <- matrix(vals, nrow = 6, ncol = 2,
                dimnames = list(paste0("CL", 1:6), c("BRCA1", "BRCA2")))
  res <- compute_dd(dep, "BRCA1", "BRCA2",
                    mut_lines = c("CL1", "CL2", "CL3"),
                    wt_lines  = c("CL4", "CL5", "CL6"),
                    min_samples = 2)
  expect_type(res, "list")
  expect_named(res, c("DD", "p_value", "cohens_d", "n_mut", "n_wt"))
  expect_true(is.numeric(res$DD))
  expect_equal(res$n_mut, 3L)
  expect_equal(res$n_wt, 3L)
})
test_that("compute_dd returns NA for insufficient samples", {
  vals <- c(-0.5, -0.3, -0.8, -0.4,
            -0.1, -0.2, -0.3, -0.6)
  dep <- matrix(vals, nrow = 4, ncol = 2,
                dimnames = list(paste0("CL", 1:4), c("BRCA1", "BRCA2")))
  res <- compute_dd(dep, "BRCA1", "BRCA2",
                    mut_lines = c("CL1"),
                    wt_lines  = c("CL2", "CL3", "CL4"),
                    min_samples = 3)
  expect_true(is.na(res$DD))
  expect_equal(res$n_mut, 1L)
})
test_that("build_mutation_matrix produces binary matrix", {
  mut_df <- data.frame(
    DepMap_ID = c("CL1", "CL1", "CL2"),
    Gene = c("TP53", "BRCA1", "TP53"),
    stringsAsFactors = FALSE
  )
  mat <- build_mutation_matrix(mut_df,
                                cell_lines = c("CL1", "CL2", "CL3"),
                                genes = c("TP53", "BRCA1"))
  expect_equal(dim(mat), c(3L, 2L))
  expect_equal(mat["CL1", "TP53"], 1)
  expect_equal(mat["CL1", "BRCA1"], 1)
  expect_equal(mat["CL2", "TP53"], 1)
  expect_equal(mat["CL3", "TP53"], 0)
  expect_equal(mat["CL3", "BRCA1"], 0)
  expect_true(all(mat == 0 | mat == 1))
})
test_that("compute_pcs handles edge cases", {
  dep_vals <- c(-0.5, -0.3, -0.1, -0.8, -0.4, -0.2,
                -0.3, -0.1, 0.1,  -0.5, -0.2, -0.3)
  dep <- matrix(dep_vals, nrow = 6, ncol = 2,
                dimnames = list(paste0("CL", 1:6), c("BRCA1", "BRCA2")))
  expr_vals <- c(2.0, 2.5, 1.8, 1.5, 1.2, 1.1,
                 1.8, 2.1, 1.6, 1.3, 1.4, 1.0)
  expr <- matrix(expr_vals, nrow = 6, ncol = 2,
                 dimnames = list(paste0("CL", 1:6), c("BRCA1", "BRCA2")))
  res <- compute_pcs(dep, expr, "BRCA1", "BRCA2",
                     mut_lines = c("CL1", "CL2", "CL3"),
                     wt_lines  = c("CL4", "CL5", "CL6"))
  expect_type(res, "list")
  expect_named(res, c("PCS", "delta_expression", "necessity"))
  expect_true(is.numeric(res$PCS))
  expect_true(is.numeric(res$delta_expression))
})
test_that("classify_msi_status returns correct labels", {
  mut_df <- data.frame(
    DepMap_ID = c("CL1", "CL3"),
    Gene = c("MLH1", "MSH6"),
    stringsAsFactors = FALSE
  )
  status <- classify_msi_status(mut_df, c("CL1", "CL2", "CL3", "CL4"))
  expect_named(status, c("CL1", "CL2", "CL3", "CL4"))
  expect_equal(status[["CL1"]], "MSI_H")
  expect_equal(status[["CL2"]], "MSS")
  expect_equal(status[["CL3"]], "MSI_H")
  expect_equal(status[["CL4"]], "MSS")
})
test_that("compute_therapeutic_window classifies correctly", {
  set.seed(42)
  vals <- rnorm(200, mean = -0.3, sd = 0.3)
  dep <- matrix(vals, nrow = 100, ncol = 2,
                dimnames = list(paste0("CL", 1:100), c("ARID1B", "OTHER")))
  dep[1:30, "ARID1B"] <- rnorm(30, mean = -0.9, sd = 0.2)
  res <- compute_therapeutic_window(dep, "ARID1A", "ARID1B",
                                     mut_lines = paste0("CL", 1:30),
                                     wt_lines  = paste0("CL", 31:100),
                                     all_cell_lines = paste0("CL", 1:100))
  expect_type(res, "list")
  expect_named(res, c("TI", "selectivity", "pan_essential_frac",
                      "DD", "classification"))
  expect_true(is.numeric(res$TI))
  expect_true(res$classification %in% c("PAN_ESSENTIAL", "HIGH_SELECTIVITY",
                                         "MODERATE", "LOW_SELECTIVITY",
                                         "INSUFFICIENT_DATA"))
})
test_that("stratify_by_mutation_type splits correctly", {
  mut_df <- data.frame(
    DepMap_ID = c("CL1", "CL2", "CL3", "CL4"),
    Gene = rep("BRCA1", 4),
    VariantInfo = c("frameshift_variant", "missense_variant",
                    "stop_gained", "missense_variant"),
    stringsAsFactors = FALSE
  )
  res <- stratify_by_mutation_type(mut_df,
                                    cell_lines = paste0("CL", 1:4),
                                    driver_gene = "BRCA1")
  expect_type(res, "list")
  expect_equal(sort(res$trunc_lines), c("CL1", "CL3"))
  expect_equal(sort(res$miss_lines), c("CL2", "CL4"))
})
