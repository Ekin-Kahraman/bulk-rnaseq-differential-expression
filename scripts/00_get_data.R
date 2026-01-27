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

# ---- 6. Read raw count matrix ----
# Expectation:
# - first column = gene identifiers
# - remaining columns = sample counts

counts_raw <- read_delim(
  count_files,
  delim = " ",
  col_names = TRUE,
  trim_ws = TRUE,
  show_col_types = FALSE
)

# Inspect structure after parsing
dim(counts_raw)
head(counts_raw[, 1:6])

# ---- 7. Fetch GEO sample metadata ----
# Purpose: identify experimental conditions for subsetting

gse <- getGEO(geo_accession, GSEMatrix = TRUE)
stopifnot(length(gse) == 1)

meta <- pData(gse[[1]])

# ---- Decision log: dataset subset for CV-grade analysis ----
# GEO metadata shows strong class imbalance (pos=430, neg=54; all Homo sapiens).
# Analysing all 484 samples would introduce unnecessary heterogeneity and batch
# structure without increasing the signal a reviewer cares about.
# We therefore construct a balanced 2-condition subset (n=30/30) to maximise
# interpretability, statistical comparability, and interview explainability.

# ---- 8. Select balanced subset of samples (POS IDs) ----
set.seed(42)

neg_samples <- meta$title[meta$`sars-cov-2 positivity:ch1` == "neg"]
pos_samples <- meta$title[meta$`sars-cov-2 positivity:ch1` == "pos"]

neg_keep <- sample(neg_samples, 30)
pos_keep <- sample(pos_samples, 30)

selected_samples <- c(neg_keep, pos_keep)
stopifnot(length(selected_samples) == 60)

# ---- 9. Subset count matrix to selected samples ----

common_samples <- intersect(selected_samples, colnames(counts_raw))
stopifnot(length(common_samples) == 60)

counts_subset <- counts_raw %>%
  select(1, all_of(common_samples))

dim(counts_subset)

# ---- 10. Finalise count matrix ----
# Convert to numeric gene × sample matrix

gene_ids <- counts_subset[[1]]
count_matrix <- as.matrix(counts_subset[, -1])

rownames(count_matrix) <- gene_ids
mode(count_matrix) <- "numeric"

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
