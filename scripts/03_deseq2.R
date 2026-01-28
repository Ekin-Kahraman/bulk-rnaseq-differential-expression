# Script: 03_deseq2.R
# Purpose: Identify genes transcriptionally altered by SARS-CoV-2 infection

# This analysis tests whether SARS-CoV-2 infection is associated with
# systematic changes in host gene expression.
#
# The goal is not simply to find “different genes”, but to quantify
# which biological programs are activated or suppressed in infected cells.
# These gene-level changes form the basis for downstream pathway
# and immune-response interpretation.

suppressPackageStartupMessages({
  library(DESeq2)
})

# Load raw, filtered count data
# DESeq2 models integer RNA-seq counts using a negative binomial framework,
# which captures both biological variability and sequencing noise.
counts   <- readRDS("data/counts_clean.rds")
metadata <- readRDS("data/metadata.rds")

# Validate inputs before modelling
# Errors at this stage typically indicate upstream problems that would
# invalidate biological conclusions if ignored.
stopifnot(
  is.matrix(counts),
  storage.mode(counts) == "integer",
  !any(is.na(counts)),
  !any(is.na(metadata)),
  nrow(counts) > 0,
  ncol(counts) > 0,
  max(counts) > 0
)

# Align metadata to expression data
# Correct sample matching is critical: misalignment silently produces
# false biological signals.
if (!all(colnames(counts) %in% rownames(metadata))) {
  stop("Metadata missing samples present in count matrix")
}
metadata <- metadata[colnames(counts), , drop = FALSE]

stopifnot(
  ncol(counts) == nrow(metadata),
  all(colnames(counts) == rownames(metadata))
)

# Ensure experimental condition is correctly defined
# The condition variable represents infection status and is the sole
# biological factor tested in this model.
stopifnot(
  "condition" %in% colnames(metadata),
  all(metadata$condition %in% c("positive", "negative"))
)
metadata$condition <- factor(metadata$condition)

# Construct DESeq2 dataset
# Each gene is modelled independently, testing whether infection status
# explains observed expression differences beyond random variation.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Genes with zero counts across all samples contain no information
# and cannot contribute to differential expression inference.
dds <- dds[rowSums(counts(dds)) > 0, ]

stopifnot(
  nrow(dds) > 0,
  max(counts(dds)) > 0
)

# Set uninfected samples as the reference state
# This makes log2 fold changes biologically interpretable as
# infection-induced up- or down-regulation.
dds$condition <- relevel(dds$condition, ref = "negative")

# Fit the DESeq2 model
# This step estimates size factors, dispersion parameters,
# and gene-wise statistical tests.
dds <- DESeq(dds)

# Extract differential expression results
# Each gene is tested for evidence that its expression differs
# between infected and uninfected samples.
res <- results(
  dds,
  contrast = c("condition", "positive", "negative"),
  alpha = 0.05
)

# Rank genes by statistical significance
# Strongly significant genes often reflect host antiviral responses,
# immune signalling, or virus-driven transcriptional reprogramming.
res <- res[order(res$padj), ]

# Save outputs for downstream interpretation
# These results form the input for pathway enrichment and biological synthesis.
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

write.csv(
  as.data.frame(res),
  file = "results/tables/deseq2_results.csv",
  row.names = TRUE
)

saveRDS(dds, "data/dds_object.rds")

message("DESeq2 analysis complete")
