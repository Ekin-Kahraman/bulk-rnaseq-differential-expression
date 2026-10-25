source(testthat::test_path("..", "..", "scripts", "enrichment_input.R"), local = TRUE)

testthat::test_that("the background contains exactly genes eligible for DE selection", {
  results <- data.frame(
    gene = c("up", "down", "nonsignificant", "independent_filtered", "unestimated"),
    padj = c(0.01, 0.02, 0.8, NA, 0.03),
    log2FoldChange = c(2, -2, 0.1, 2, NA)
  )
  genes <- select_enrichment_genes(results, 0.05, 1)
  testthat::expect_identical(genes$universe, c("up", "down", "nonsignificant"))
  testthat::expect_identical(genes$selected, c("up", "down"))
  testthat::expect_true(all(genes$selected %in% genes$universe))

  # One-to-many mappings are applied once to both sets; off-background database
  # genes and unmapped IDs must never enlarge the comparison population.
  mapping <- data.frame(
    ENSEMBL = c("up", "up", "down", "nonsignificant", "database_only", "up"),
    ENTREZID = c("1", "2", "3", "4", "5", NA)
  )
  mapped <- map_enrichment_genes(genes, mapping, "ENSEMBL")
  testthat::expect_identical(mapped$selected, c("1", "2", "3"))
  testthat::expect_identical(mapped$universe, c("1", "2", "3", "4"))
  testthat::expect_false("5" %in% mapped$universe)
})

testthat::test_that("published background and both enrichment denominators agree", {
  root <- testthat::test_path("..", "..")
  results <- read.csv(file.path(root, "results/tables/deseq2_results.csv"))
  background <- read.csv(file.path(root, "results/tables/enrichment_background.csv"),
                         colClasses = c("character", "character", "logical"))
  summary <- read.csv(file.path(root, "results/tables/enrichment_summary.csv"))
  analysis <- read.csv(file.path(root, "results/tables/analysis_summary.csv"))
  genes <- select_enrichment_genes(results, 0.05, 1)
  testthat::expect_true(all(background$gene %in% genes$universe))
  testthat::expect_identical(background$selected, background$gene %in% genes$selected)
  testthat::expect_true(all(summary$eligible_genes == length(genes$universe)))
  testthat::expect_true(all(summary$mapped_background == length(unique(background$ENTREZID))))
  for (i in seq_len(nrow(summary))) {
    name <- if (summary$database[i] == "GO_BP") "go_biological_process.csv" else "kegg_pathways.csv"
    enriched <- read.csv(file.path(root, "results/tables", name))
    denominator <- as.integer(sub(".*/", "", enriched$BgRatio))
    testthat::expect_true(all(denominator == summary$annotated_background[i]))
    testthat::expect_equal(nrow(enriched), summary$enriched_terms[i])
    aggregate_count <- if (summary$database[i] == "GO_BP") {
      analysis$go_terms_enriched
    } else {
      analysis$kegg_pathways_enriched
    }
    testthat::expect_equal(aggregate_count, nrow(enriched))
  }
})
