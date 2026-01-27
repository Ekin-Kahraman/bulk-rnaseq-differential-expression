# ============================================================
# Script: 01_qc.R
# Purpose: Quality control and initial exploration of bulk RNA-seq data
# Project: Bulk RNA-seq Differential Expression Analysis
# Author: Ekin Kahraman
# Date: 2026-01-25
# ============================================================

# ---- 1. Load required libraries ----
library(DESeq2)
library(edgeR)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)

# ---- 2. Load frozen inputs ----
# These objects are produced by 00_get_data.R and are treated as immutable

counts   <- readRDS("data/counts.rds")
metadata <- read_csv("data/metadata.csv", show_col_types = FALSE)

# Sanity check: initial alignment
stopifnot(all(colnames(counts) == metadata$sample_id))

# ---- 3. Library size QC (pre-cleaning) ----
# Rationale:
# Large differences in sequencing depth can indicate technical artefacts.
# We inspect library sizes before any filtering or modelling decisions.

library_sizes <- colSums(counts, na.rm = TRUE)

qc_libsize <- tibble(
  sample_id    = names(library_sizes),
  library_size = library_sizes,
  condition    = metadata$condition
)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

p_libsize <- ggplot(qc_libsize, aes(x = condition, y = library_size)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_y_log10() +
  labs(
    title = "Library size distribution by condition",
    y = "Total counts (log10)",
    x = NULL
  ) +
  theme_minimal()

ggsave(
  filename = "results/figures/qc_library_size.png",
  plot     = p_libsize,
  width    = 6,
  height   = 4
)

# ---- 4. Remove invalid sample columns (FAILED LIBRARIES) ----
# Rationale:
# CPM, PCA, and DESeq2 require numeric count matrices with no NA values.
# Any sample that is entirely NA or has zero total counts is a failed library
# and must be removed (not imputed).

counts <- as.matrix(counts)

# Identify samples containing ANY NA values
na_cols <- colSums(is.na(counts)) > 0

if (any(na_cols)) {
  bad_samples <- colnames(counts)[na_cols]
  counts   <- counts[, !na_cols, drop = FALSE]
  metadata <- metadata %>% filter(!sample_id %in% bad_samples)
}

# Recompute library sizes after NA-column removal
library_sizes <- colSums(counts)

# Identify zero-library samples
zero_cols <- library_sizes == 0

if (any(zero_cols)) {
  bad_samples <- colnames(counts)[zero_cols]
  counts   <- counts[, !zero_cols, drop = FALSE]
  metadata <- metadata %>% filter(!sample_id %in% bad_samples)
}

# Re-align metadata order to count columns (CRITICAL)
metadata <- metadata %>% slice(match(colnames(counts), sample_id))

# Hard invariants (fail fast if violated)
stopifnot(!any(is.na(counts)))
stopifnot(all(colnames(counts) == metadata$sample_id))
stopifnot(all(colSums(counts) > 0))

# ---- 5. Low-expression gene filtering (CPM-based) ----
# Rationale:
# In deeply sequenced datasets, raw-count thresholds are insufficient.
# We retain genes with CPM ≥ 1 in at least 10 samples, a standard bulk RNA-seq filter.

cpm_values <- cpm(counts)

keep_genes <- rowSums(cpm_values >= 1) >= 10
counts_filtered <- counts[keep_genes, , drop = FALSE]

n_before <- nrow(counts)
n_after  <- nrow(counts_filtered)

n_before
n_after


# ---- QC summary ----
# One failed library was removed due to NA / zero counts.
# CPM-based filtering retained ~14,744 expressed genes across 59 samples.
# Dataset passes QC and is suitable for PCA and DE analysis.


# ---- 6. Ensure correct data types for DESeq2 ----
# Rationale:
# DESeq2 requires integer counts and categorical design variables.

metadata$condition <- factor(
  metadata$condition,
  levels = c("SARS_CoV_2_negative", "SARS_CoV_2_positive")
)

mode(counts_filtered) <- "integer"

# Final sanity checks before modelling
stopifnot(nrow(counts_filtered) > 0)
stopifnot(all(colSums(counts_filtered) > 0))

# ---- 7. Variance stabilisation (for PCA / visualisation) ----
# Rationale:
# PCA assumes approximately homoscedastic data.
# Variance stabilising transformation corrects mean–variance dependence
# while preserving sample relationships.

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = metadata,
  design    = ~ condition
)

vsd <- vst(dds, blind = TRUE)

# Save transformed object for reuse
saveRDS(vsd, "data/vst_object.rds")
