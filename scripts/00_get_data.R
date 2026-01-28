# ============================================================
# Script: 00_get_data.R
# Purpose: Download and prepare raw RNA-seq count data
# Project: Bulk RNA-seq Differential Expression Analysis
# Author: Ekin Kahraman
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(tibble)
})

# ------------------------------------------------------------
# 1. GEO accession
# ------------------------------------------------------------

geo_accession <- "GSE152075"

# ------------------------------------------------------------
# 2. Directory structure
# ------------------------------------------------------------

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Download supplementary files (raw counts)
# ------------------------------------------------------------

getGEOSuppFiles(
  geo_accession,
  makeDirectory = TRUE,
  baseDir = "data/raw"
)

# ------------------------------------------------------------
# 4. Locate count matrix file
# ------------------------------------------------------------

count_files <- list.files(
  path = file.path("data/raw", geo_accession),
  pattern = "txt|tsv|counts",
  full.names = TRUE
)

stopifnot(length(count_files) == 1)

count_file <- count_files[1]

# ------------------------------------------------------------
# 5. Read raw count matrix
# ------------------------------------------------------------

counts_df <- read.table(
  count_file,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)

# First column = gene IDs
gene_ids <- counts_df[[1]]

# Remaining columns = counts
counts_raw <- counts_df[, -1]
counts_raw[] <- lapply(counts_raw, as.numeric)

counts_raw <- as.matrix(counts_raw)
rownames(counts_raw) <- gene_ids

# ------------------------------------------------------------
# 6. Remove genes with parsing failures
# ------------------------------------------------------------

bad_genes <- rowSums(is.na(counts_raw)) > 0
counts_raw <- counts_raw[!bad_genes, , drop = FALSE]

# ------------------------------------------------------------
# 7. Enforce integer raw counts (safely)
# ------------------------------------------------------------

if (any(abs(counts_raw - round(counts_raw)) > 1e-6)) {
  stop("Non-integer values detected: input is NOT raw counts")
}
storage.mode(counts_raw) <- "integer"

# ------------------------------------------------------------
# 8. Basic sanity checks
# ------------------------------------------------------------

stopifnot(
  nrow(counts_raw) > 0,
  ncol(counts_raw) > 0,
  max(counts_raw) > 0,
  all(colSums(counts_raw) > 0)
)

# ------------------------------------------------------------
# 9. Infer condition directly from column names (ROBUST)
# ------------------------------------------------------------

sample_names <- colnames(counts_raw)

# Inspect once (debug)
# head(sample_names)

condition <- ifelse(
  grepl("covid|sars|positive|pos", sample_names, ignore.case = TRUE),
  "positive",
  ifelse(
    grepl("control|negative|neg|healthy", sample_names, ignore.case = TRUE),
    "negative",
    NA
  )
)

# Drop samples with unknown status
keep <- !is.na(condition)

counts_labeled <- counts_raw[, keep, drop = FALSE]
condition <- condition[keep]
sample_names <- sample_names[keep]

table(condition)  # SHOULD SHOW ~430 positive / ~54 negative

# ------------------------------------------------------------
# 10. Construct balanced subset (30 / 30)
# ------------------------------------------------------------

set.seed(42)

pos_idx <- which(condition == "positive")
neg_idx <- which(condition == "negative")

stopifnot(
  length(pos_idx) >= 30,
  length(neg_idx) >= 30
)

pos_keep <- sample(pos_idx, 30)
neg_keep <- sample(neg_idx, 30)

keep_idx <- c(pos_keep, neg_keep)

counts_bal <- counts_labeled[, keep_idx, drop = FALSE]

metadata <- data.frame(
  condition = condition[keep_idx],
  row.names = colnames(counts_bal),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 11. Final sanity checks
# ------------------------------------------------------------

stopifnot(
  ncol(counts_bal) == 60,
  nrow(metadata) == 60,
  all(colnames(counts_bal) == rownames(metadata)),
  max(counts_bal) > 0
)

# ------------------------------------------------------------
# 12. Save frozen raw data
# ------------------------------------------------------------

saveRDS(counts_bal, "data/counts_raw.rds")
saveRDS(metadata,    "data/metadata.rds")

message("00_get_data.R complete")
message("Genes: ", nrow(counts_bal), " | Samples: ", ncol(counts_bal))

