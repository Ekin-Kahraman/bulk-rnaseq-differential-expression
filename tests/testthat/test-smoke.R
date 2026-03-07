library(testthat)

find_project_root <- function() {
  dir <- getwd()
  for (i in 0:10) {
    if (file.exists(file.path(dir, "renv.lock")) && file.exists(file.path(dir, "run_all.R"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) break
    dir <- parent
  }
  stop("Could not locate project root (expected renv.lock + run_all.R).")
}

root <- find_project_root()
source(file.path(root, "scripts/config.R"), local = TRUE)

test_that("DESeq2 results CSV exists and is well-formed", {
  path <- file.path(root, "results/tables/deseq2_results.csv")
  expect_true(file.exists(path))
  skip_if_not(file.exists(path))

  res <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  expect_true(all(required_cols %in% names(res)))

  expect_gt(nrow(res), 1000)
  expect_true(any(!is.na(res$padj)))
  expect_true(all(nchar(res$gene) > 0))
  expect_gt(sum(res$padj < analysis_config$de_padj_cutoff, na.rm = TRUE), 0)

  pvals <- res$pvalue[!is.na(res$pvalue)]
  expect_true(length(pvals) > 0)
  expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that("key figures exist", {
  volcano <- file.path(root, "results/figures/volcano_plot.png")
  pca <- file.path(root, "results/figures/pca_plot.png")

  expect_true(file.exists(volcano))
  expect_true(file.exists(pca))
})

test_that("key derived tables exist", {
  top_genes <- file.path(root, "results/tables/top_genes.csv")
  go_terms <- file.path(root, "results/tables/go_biological_process.csv")
  kegg_terms <- file.path(root, "results/tables/kegg_pathways.csv")

  expect_true(file.exists(top_genes))
  expect_true(file.exists(go_terms))
  expect_true(file.exists(kegg_terms))

  top_df <- read.csv(top_genes, stringsAsFactors = FALSE)
  go_df <- read.csv(go_terms, stringsAsFactors = FALSE)
  kegg_df <- read.csv(kegg_terms, stringsAsFactors = FALSE)

  expect_true(all(c("gene", "log2FoldChange", "padj", "baseMean") %in% names(top_df)))
  expect_gt(nrow(top_df), 0)

  expect_true(all(c("Description", "p.adjust", "Count") %in% names(go_df)))
  expect_gt(nrow(go_df), 0)

  expect_true(all(c("Description", "p.adjust", "Count") %in% names(kegg_df)))
  expect_gt(nrow(kegg_df), 0)
})

test_that("legacy tracked outputs are removed", {
  legacy <- file.path(root, "results/tables/top_volcano_genes.csv")
  expect_false(file.exists(legacy))
})
