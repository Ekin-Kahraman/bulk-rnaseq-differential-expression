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
  shrunken_path <- file.path(root, "results/tables/deseq2_results_shrunken.csv")
  expect_true(file.exists(path))
  expect_true(file.exists(shrunken_path))
  skip_if_not(file.exists(path))

  res <- read.csv(path, stringsAsFactors = FALSE)
  res_shrunk <- read.csv(shrunken_path, stringsAsFactors = FALSE)
  required_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  expect_true(all(required_cols %in% names(res)))
  expect_true(all(required_cols %in% names(res_shrunk)))
  expect_equal(sort(res$gene), sort(res_shrunk$gene))

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
  sensitivity <- file.path(root, "results/figures/sensitivity_lfc_scatter.png")

  expect_true(file.exists(volcano))
  expect_true(file.exists(pca))
  expect_true(file.exists(sensitivity))
})

test_that("key derived tables exist", {
  top_genes <- file.path(root, "results/tables/top_genes.csv")
  go_terms <- file.path(root, "results/tables/go_biological_process.csv")
  kegg_terms <- file.path(root, "results/tables/kegg_pathways.csv")
  full_cohort <- file.path(root, "results/tables/full_cohort_deseq2_results.csv")
  summary_path <- file.path(root, "results/tables/analysis_summary.csv")

  expect_true(file.exists(top_genes))
  expect_true(file.exists(go_terms))
  expect_true(file.exists(kegg_terms))
  expect_true(file.exists(full_cohort))
  expect_true(file.exists(summary_path))

  top_df <- read.csv(top_genes, stringsAsFactors = FALSE)
  go_df <- read.csv(go_terms, stringsAsFactors = FALSE)
  kegg_df <- read.csv(kegg_terms, stringsAsFactors = FALSE)
  full_df <- read.csv(full_cohort, stringsAsFactors = FALSE)
  summary_df <- read.csv(summary_path, stringsAsFactors = FALSE)

  expect_true(all(c("gene", "log2FoldChange", "padj", "baseMean") %in% names(top_df)))
  expect_gt(nrow(top_df), 0)

  expect_true(all(c("Description", "p.adjust", "Count") %in% names(go_df)))
  expect_gt(nrow(go_df), 0)

  expect_true(all(c("Description", "p.adjust", "Count") %in% names(kegg_df)))
  expect_gt(nrow(kegg_df), 0)

  expect_true(all(c("gene", "baseMean", "log2FoldChange", "padj", "shrunkenLog2FoldChange") %in% names(full_df)))
  expect_gt(nrow(full_df), 1000)

  expect_equal(nrow(summary_df), 1)
  expect_true(all(c(
    "qc_samples_total",
    "balanced_sig_genes",
    "full_cohort_sig_genes",
    "shared_threshold_genes",
    "shared_sign_concordance",
    "shrunken_lfc_spearman"
  ) %in% names(summary_df)))
  expect_gt(summary_df$qc_samples_total[[1]], 100)
  expect_gt(summary_df$shared_threshold_genes[[1]], 100)
  expect_gt(summary_df$shared_sign_concordance[[1]], 0.8)
  expect_gt(summary_df$shrunken_lfc_spearman[[1]], 0.8)
})

test_that("output manifest tracks committed result artefacts", {
  manifest_path <- file.path(root, "results/tables/output_manifest.csv")
  expect_true(file.exists(manifest_path))

  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  expect_true(all(c("path", "bytes", "md5") %in% names(manifest)))
  expect_gt(nrow(manifest), 20)
  expect_true(all(file.exists(file.path(root, manifest$path))))
  expect_true(all(manifest$bytes > 0))
  expect_true(all(grepl("^[0-9a-f]{32}$", manifest$md5)))

  actual_sizes <- file.info(file.path(root, manifest$path))$size
  expect_equal(as.numeric(manifest$bytes), as.numeric(actual_sizes))
  expect_true("results/tables/analysis_summary.csv" %in% manifest$path)
  expect_true("results/figures/volcano_plot.png" %in% manifest$path)
  expect_true("results/figures/sensitivity_lfc_scatter.png" %in% manifest$path)
})

test_that("KEGG enrichment uses pinned human pathway references", {
  links_path <- file.path(root, "data/reference/kegg_hsa_pathway_links.tsv")
  names_path <- file.path(root, "data/reference/kegg_hsa_pathway_names.tsv")

  expect_true(file.exists(links_path))
  expect_true(file.exists(names_path))

  links <- read.delim(links_path, header = FALSE, stringsAsFactors = FALSE)
  names <- read.delim(names_path, header = FALSE, stringsAsFactors = FALSE)

  expect_gt(nrow(links), 9000)
  expect_gt(nrow(names), 300)
  expect_true(all(grepl("^path:hsa", links[[1]])))
  expect_true(all(grepl("^hsa:", links[[2]])))
  expect_true(all(grepl("^hsa", names[[1]])))
})

test_that("legacy tracked outputs are removed", {
  legacy <- file.path(root, "results/tables/top_volcano_genes.csv")
  expect_false(file.exists(legacy))
})
