# Script: 01_qc.R
# Purpose: Quality control and filtering of bulk RNA-seq data
# Project: Bulk RNA-seq Differential Expression Analysis
# Author: Ekin Kahraman

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
})

# Load processed raw inputs
counts   <- readRDS("data/counts_raw.rds")   # raw integer counts
metadata <- readRDS("data/metadata.rds")     # sample metadata

# Basic integrity checks
stopifnot(
  is.matrix(counts),
  storage.mode(counts) == "integer",
  nrow(counts) > 0,
  ncol(counts) > 0,
  max(counts) > 0,
  !any(is.na(counts)),
  nrow(metadata) == ncol(counts),
  all(colnames(counts) == rownames(metadata))
)

# Collapse duplicated gene identifiers by summing counts
# (required for DESeq2 compatibility)
if (any(duplicated(rownames(counts)))) {
  counts <- rowsum(counts, group = rownames(counts))
}

stopifnot(
  !any(duplicated(rownames(counts))),
  nrow(counts) > 0
)

# Inspect library sizes
library_sizes <- colSums(counts)

qc_df <- data.frame(
  sample    = colnames(counts),
  library   = library_sizes,
  condition = metadata$condition,
  stringsAsFactors = FALSE
)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

p_libsize <- ggplot(qc_df, aes(x = condition, y = library)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  scale_y_log10() +
  labs(
    title = "Library size distribution by condition",
    y = "Total counts (log10)",
    x = NULL
  ) +
  theme_minimal()

ggsave(
  "results/figures/qc_library_size.png",
  p_libsize,
  width = 6,
  height = 4
)

# Remove samples with zero total counts
keep_samples <- library_sizes > 0

counts   <- counts[, keep_samples, drop = FALSE]
metadata <- metadata[keep_samples, , drop = FALSE]

stopifnot(
  ncol(counts) > 0,
  all(colSums(counts) > 0),
  all(colnames(counts) == rownames(metadata))
)

# Filter lowly expressed genes using CPM
# (raw counts retained for downstream modelling)
cpm_mat <- edgeR::cpm(counts)

keep_genes <- rowSums(cpm_mat >= 1) >= 10

if (sum(keep_genes) == 0) {
  stop("CPM filtering removed all genes — check thresholds.")
}

counts_filt <- counts[keep_genes, , drop = FALSE]

stopifnot(
  nrow(counts_filt) > 0,
  ncol(counts_filt) > 0,
  max(counts_filt) > 0
)

message(
  "QC complete: retained ",
  nrow(counts_filt), " genes across ",
  ncol(counts_filt), " samples"
)

# Save filtered raw counts
saveRDS(counts_filt, "data/counts_clean.rds")
saveRDS(metadata,    "data/metadata_clean.rds")

# Variance stabilisation for exploratory analysis only
metadata$condition <- factor(metadata$condition)

dds_tmp <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = metadata,
  design    = ~ condition
)

dds_tmp <- dds_tmp[rowSums(counts(dds_tmp)) > 0, ]

vsd <- varianceStabilizingTransformation(dds_tmp, blind = TRUE)

saveRDS(vsd, "data/vst_object.rds")

message("01_qc.R complete")
