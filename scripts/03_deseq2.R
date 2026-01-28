# ============================================================
# Script: 03_deseq2.R
# Purpose: Differential expression analysis with DESeq2
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

# ------------------------------------------------------------
# 1. Load inputs (RAW COUNTS ONLY — REQUIRED)
# ------------------------------------------------------------

counts <- readRDS("data/counts_clean.rds")
metadata <- readRDS("data/metadata.rds")

# ------------------------------------------------------------
# 2. Hard invariants (FAIL FAST)
# ------------------------------------------------------------

stopifnot(
  is.matrix(counts),
  storage.mode(counts) == "integer",
  !any(is.na(counts)),
  !any(is.na(metadata)),
  nrow(counts) > 0,
  ncol(counts) > 0,
  max(counts) > 0                     # CRITICAL: prevents all-zero matrices
)

# Ensure metadata contains required samples
if (!all(colnames(counts) %in% rownames(metadata))) {
  stop("Metadata missing samples present in count matrix")
}

# Align metadata to counts (authoritative source)
metadata <- metadata[colnames(counts), , drop = FALSE]

# Final alignment check
stopifnot(
  ncol(counts) == nrow(metadata),
  all(colnames(counts) == rownames(metadata))
)

# Ensure condition exists and is valid
stopifnot(
  "condition" %in% colnames(metadata),
  all(metadata$condition %in% c("positive", "negative"))
)

metadata$condition <- factor(metadata$condition)

# ------------------------------------------------------------
# 3. Construct DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Remove genes with zero counts across all samples
dds <- dds[rowSums(counts(dds)) > 0, ]

# Sanity check post-filtering
stopifnot(
  nrow(dds) > 0,
  max(counts(dds)) > 0
)

# Explicit reference level (interpretability)
dds$condition <- relevel(dds$condition, ref = "negative")

# ------------------------------------------------------------
# 4. Run DESeq2 model
# ------------------------------------------------------------

dds <- DESeq(dds)

# ------------------------------------------------------------
# 5. Extract results
# ------------------------------------------------------------

res <- results(
  dds,
  contrast = c("condition", "positive", "negative"),
  alpha = 0.05
)

res <- res[order(res$padj), ]

# ------------------------------------------------------------
# 6. Export outputs
# ------------------------------------------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

write.csv(
  as.data.frame(res),
  file = "results/tables/deseq2_results.csv",
  row.names = TRUE
)

saveRDS(dds, "data/dds_object.rds")

message("DESeq2 analysis complete: results written to results/tables/")
