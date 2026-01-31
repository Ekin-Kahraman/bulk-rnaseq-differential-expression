#!/usr/bin/env Rscript
# Differential expression analysis with DESeq2

library(DESeq2)

if (!file.exists("data/counts_clean.rds")) {
  stop("File not found: data/counts_clean.rds. Run scripts/01_qc.R first.")
}
if (!file.exists("data/metadata_clean.rds")) {
  stop("File not found: data/metadata_clean.rds. Run scripts/01_qc.R first.")
}

counts <- readRDS("data/counts_clean.rds")
metadata <- readRDS("data/metadata_clean.rds")

message("Testing ", nrow(counts), " genes across ", ncol(counts), " samples")

dds <- DESeqDataSetFromMatrix(counts, metadata, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "negative")

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "positive", "negative"), alpha = 0.05)
res <- res[order(res$padj), ]

n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
n_down <- sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)

message(n_sig, " significant genes (padj < 0.05)")
message("  ", n_up, " upregulated (log2FC > 1)")
message("  ", n_down, " downregulated (log2FC < -1)")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res_df, "results/tables/deseq2_results.csv", row.names = FALSE)
saveRDS(dds, "data/dds_object.rds")

message("Results saved")