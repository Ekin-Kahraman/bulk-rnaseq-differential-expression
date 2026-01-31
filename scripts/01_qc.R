#!/usr/bin/env Rscript
# Quality control: filter samples and genes, normalize for visualization

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

keep <- qc_data$library_size > 1e5
counts <- counts[, keep, drop = FALSE]
metadata <- metadata[keep, , drop = FALSE]

message(sum(keep), " samples passed QC")

set.seed(123)
n <- 30
pos_samples <- rownames(metadata)[metadata$condition == "positive"]
neg_samples <- rownames(metadata)[metadata$condition == "negative"]
balanced <- c(
  sample(pos_samples, min(n, length(pos_samples))),
  sample(neg_samples, min(n, length(neg_samples)))
)
counts <- counts[, balanced, drop = FALSE]
metadata <- metadata[balanced, , drop = FALSE]

rownames(counts) <- sub("\\..*", "", rownames(counts))

if (any(duplicated(rownames(counts)))) {
  n_dup <- sum(duplicated(rownames(counts)))
  counts <- rowsum(counts, rownames(counts))
  message("Collapsed ", n_dup, " duplicates")
}

keep_genes <- rowSums(cpm(counts) >= 1) >= 10
counts <- counts[keep_genes, , drop = FALSE]

message(nrow(counts), " genes retained after filtering")
message("Final: ", sum(metadata$condition == "negative"), " neg, ",
        sum(metadata$condition == "positive"), " pos")

saveRDS(counts, "data/counts_clean.rds")
saveRDS(metadata, "data/metadata_clean.rds")

dds <- DESeqDataSetFromMatrix(counts, metadata, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
saveRDS(vsd, "data/vst_data.rds")

message("QC complete")