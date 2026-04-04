#!/usr/bin/env Rscript
# Quality control: filter samples and genes, normalize for visualization

source("scripts/config.R", local = TRUE)

library(DESeq2)
library(edgeR)
library(ggplot2)
library(scales)

if (!file.exists("data/counts_raw.rds")) {
  stop("File not found: data/counts_raw.rds. Run scripts/00_get_data.R first.")
}
if (!file.exists("data/metadata.rds")) {
  stop("File not found: data/metadata.rds. Run scripts/00_get_data.R first.")
}

counts <- readRDS("data/counts_raw.rds")
metadata <- readRDS("data/metadata.rds")

metadata <- metadata[colnames(counts), , drop = FALSE]
if (!all(colnames(counts) == rownames(metadata))) {
  stop("Counts and metadata are not aligned by sample ID.")
}
if (any(is.na(metadata$condition))) {
  stop("Metadata contains missing condition labels. Re-run scripts/00_get_data.R.")
}

message("Starting QC: ", nrow(counts), " genes, ", ncol(counts), " samples")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

qc_data <- data.frame(
  library_size = colSums(counts),
  condition = metadata$condition,
  stringsAsFactors = FALSE
)

ggplot(qc_data, aes(condition, library_size)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85") +
  geom_jitter(width = 0.15, alpha = 0.6) +
  scale_y_log10(labels = comma) +
  labs(title = "Library Size Distribution", y = "Total Reads", x = NULL) +
  theme_classic(base_size = 12)

ggsave("results/figures/qc_library_size.png", width = 6, height = 5, dpi = 300)

keep <- qc_data$library_size > analysis_config$min_library_size
counts <- counts[, keep, drop = FALSE]
metadata <- metadata[keep, , drop = FALSE]

message(sum(keep), " samples passed QC")

rownames(counts) <- sub("\\..*", "", rownames(counts))

if (any(duplicated(rownames(counts)))) {
  n_dup <- sum(duplicated(rownames(counts)))
  counts <- rowsum(counts, rownames(counts))
  message("Collapsed ", n_dup, " duplicates")
}

saveRDS(counts, "data/counts_qc.rds")
saveRDS(metadata, "data/metadata_qc.rds")

# Filter to samples with non-NA gender for covariate-adjusted model
has_gender <- !is.na(metadata$gender)
if (sum(!has_gender) > 0) {
  message("Excluding ", sum(!has_gender), " samples with missing gender annotation")
  metadata <- metadata[has_gender, , drop = FALSE]
  counts <- counts[, rownames(metadata), drop = FALSE]
}

set.seed(analysis_config$sample_seed)
n <- analysis_config$samples_per_group
pos_samples <- rownames(metadata)[metadata$condition == "positive"]
neg_samples <- rownames(metadata)[metadata$condition == "negative"]
target_n <- min(n, length(pos_samples), length(neg_samples))

if (target_n == 0) {
  stop("At least one condition has zero samples after QC filtering.")
}
if (target_n < n) {
  message("Using ", target_n, " samples per condition (limited by available data).")
}

balanced <- c(
  sample(pos_samples, target_n),
  sample(neg_samples, target_n)
)
counts <- counts[, balanced, drop = FALSE]
metadata <- metadata[balanced, , drop = FALSE]

keep_genes <- rowSums(cpm(counts) >= analysis_config$gene_cpm_cutoff) >= analysis_config$gene_min_samples
counts <- counts[keep_genes, , drop = FALSE]

message(nrow(counts), " genes retained after filtering")
message("Final: ", sum(metadata$condition == "negative"), " neg, ",
        sum(metadata$condition == "positive"), " pos")

saveRDS(counts, "data/counts_clean.rds")
saveRDS(metadata, "data/metadata_clean.rds")

dds <- DESeqDataSetFromMatrix(counts, metadata, design = ~ condition + gender)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
saveRDS(vsd, "data/vst_data.rds")

message("QC complete")
