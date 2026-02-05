library(testthat)

test_that("DESeq2 results CSV exists and is well-formed", {
  expect_true(file.exists("results/tables/deseq2_results.csv"))

  res <- read.csv("results/tables/deseq2_results.csv", stringsAsFactors = FALSE)
  required_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  expect_true(all(required_cols %in% names(res)))

  expect_gt(nrow(res), 0)
  expect_true(any(!is.na(res$padj)))
  expect_true(all(nchar(res$gene) > 0))
})

test_that("key figures exist", {
  expect_true(file.exists("results/figures/volcano_plot.png"))
  expect_true(file.exists("results/figures/pca_plot.png"))
})
