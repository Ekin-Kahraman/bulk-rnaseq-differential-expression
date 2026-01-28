# ============================================================
# Script: 00_get_data.R
# Purpose: Download and prepare raw RNA-seq count data
# Project: Bulk RNA-seq Differential Expression Analysis
# Author: Ekin Kahraman
# ============================================================

# ---- 1. Load required libraries ----
# GEOquery: access GEO datasets
# tidyverse: data manipulation
library(GEOquery)
library(readr)
library(dplyr)
library(tibble)

# ---- 2. Define GEO accession ----
geo_accession <- "GSE152075"

# ---- 3. Create directory structure ----
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

# ---- 4. Download supplementary files (raw counts) ----
# GEO supplementary files contain the author-provided count matrix

getGEOSuppFiles(
  geo_accession,
  makeDirectory = TRUE,
  baseDir = "data/raw"
)

# ---- 5. Locate downloaded count file ----
# We expect a single space-delimited count matrix

count_files <- list.files(
  path = file.path("data/raw", geo_accession),
  pattern = "txt|tsv|counts",
  full.names = TRUE
)

stopifnot(length(count_files) == 1)
count_files

# ---- 6. Read raw count matrix (robust parsing) ----
# Rationale:
# GEO count files may contain mixed types or inconsistent delimiters.
# We explicitly parse gene IDs and coerce all count columns to numeric.

counts_raw <- read.table(
  count_files,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)

# First column = gene identifiers
gene_ids <- counts_raw[[1]]

# Remaining columns = counts
counts_mat <- counts_raw[, -1]

# Force numeric conversion (introduces NA if parsing failed)
counts_mat[] <- lapply(counts_mat, function(x) as.numeric(x))

# Convert to matrix
counts_mat <- as.matrix(counts_mat)
rownames(counts_mat) <- gene_ids

# ---- 7. Remove genes with parsing failures ----
# Genes with any NA values are removed entirely
na_genes <- rowSums(is.na(counts_mat)) > 0
sum(na_genes)  # diagnostic

counts_mat <- counts_mat[!na_genes, ]

# ---- 8. Enforce raw integer counts ----
mode(counts_mat) <- "integer"

# ---- 9. Sanity checks ----
stopifnot(
  all(counts_mat >= 0),
  all(colSums(counts_mat) > 0)
)

# ---- 10. Save clean raw count matrix ----
saveRDS(counts_mat, "data/counts.rds")


# Save frozen counts
saveRDS(count_matrix, "data/counts.rds")

# ---- 11. Create minimal metadata table ----
# Only variables required for downstream modelling are retained

metadata <- tibble(
  sample_id = colnames(count_matrix),
  condition = if_else(
    sample_id %in% pos_keep,
    "SARS_CoV_2_positive",
    "SARS_CoV_2_negative"
  )
)

write_csv(metadata, "data/metadata.csv")

# ---- 12. Final sanity checks ----
dim(count_matrix)
table(metadata$condition)

# ---- Final dataset summary (FROZEN) ----
# 35,784 genes × 60 samples
# 30 SARS-CoV-2 positive
# 30 SARS-CoV-2 negative
#
# This balanced subset was constructed from a larger public GEO dataset
# to maximise interpretability, statistical comparability, and downstream
# QC/DE clarity for a portfolio-grade bulk RNA-seq analysis.
