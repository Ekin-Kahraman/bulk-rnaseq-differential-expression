# Script: 00_get_data.R
# Purpose: Download raw RNA-seq counts and create a balanced dataset
# Project: Bulk RNA-seq Differential Expression Analysis
# Author: Ekin Kahraman

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(tibble)
})

# GEO accession
geo_accession <- "GSE152075"

# Create data directories
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

# Download supplementary GEO files
getGEOSuppFiles(
  geo_accession,
  makeDirectory = TRUE,
  baseDir = "data/raw"
)

# Locate raw count matrix
count_files <- list.files(
  path = file.path("data/raw", geo_accession),
  pattern = "raw_counts.*\\.txt(\\.gz)?$",
  full.names = TRUE
)

stopifnot(length(count_files) >= 1)
count_file <- count_files[1]

message("Using count file: ", basename(count_file))

# Read count matrix
counts_df <- read.table(
  count_file,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)

# First column contains gene identifiers
gene_ids <- counts_df[[1]]

counts_raw <- counts_df[, -1]
counts_raw[] <- lapply(counts_raw, as.numeric)

counts_raw <- as.matrix(counts_raw)
rownames(counts_raw) <- gene_ids

# Remove genes with parsing issues
bad_genes <- rowSums(is.na(counts_raw)) > 0
counts_raw <- counts_raw[!bad_genes, , drop = FALSE]

# Confirm data are raw integer counts
if (any(abs(counts_raw - round(counts_raw)) > 1e-6)) {
  stop("Non-integer values detected: input is not raw count data")
}
storage.mode(counts_raw) <- "integer"

# Basic integrity checks
stopifnot(
  nrow(counts_raw) > 0,
  ncol(counts_raw) > 0,
  max(counts_raw) > 0,
  all(colSums(counts_raw) > 0)
)

# Assign infection status from sample names
sample_names <- colnames(counts_raw)

condition <- ifelse(
  grepl("covid|sars|positive|pos", sample_names, ignore.case = TRUE),
  "positive",
  ifelse(
    grepl("control|negative|neg|healthy", sample_names, ignore.case = TRUE),
    "negative",
    NA
  )
)

# Keep only labelled samples
keep <- !is.na(condition)

counts_labeled <- counts_raw[, keep, drop = FALSE]
condition <- condition[keep]

table(condition)

# Create balanced subset (30 positive / 30 negative)
set.seed(42)

pos_idx <- which(condition == "positive")
neg_idx <- which(condition == "negative")

stopifnot(
  length(pos_idx) >= 30,
  length(neg_idx) >= 30
)

keep_idx <- c(sample(pos_idx, 30), sample(neg_idx, 30))
counts_bal <- counts_labeled[, keep_idx, drop = FALSE]

metadata <- data.frame(
  condition = condition[keep_idx],
  row.names = colnames(counts_bal),
  stringsAsFactors = FALSE
)

# Final consistency checks
stopifnot(
  ncol(counts_bal) == 60,
  nrow(metadata) == 60,
  all(colnames(counts_bal) == rownames(metadata))
)

# Save processed inputs
saveRDS(counts_bal, "data/counts_raw.rds")
saveRDS(metadata,    "data/metadata.rds")

message("00_get_data.R complete")
message("Genes: ", nrow(counts_bal), " | Samples: ", ncol(counts_bal))

