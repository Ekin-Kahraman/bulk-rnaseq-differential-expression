#!/usr/bin/env Rscript
# Differential expression analysis with DESeq2

source("scripts/config.R", local = TRUE)

library(DESeq2)

if (!requireNamespace("apeglm", quietly = TRUE)) {
  stop("Package 'apeglm' is required for log2 fold-change shrinkage. Run Rscript 000_install_dependencies.R.")
}

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

coef_name <- grep("^condition_.*_vs_.*$", resultsNames(dds), value = TRUE)
if (length(coef_name) != 1) {
  stop("Expected exactly one condition coefficient, found: ", paste(resultsNames(dds), collapse = ", "))
}

res_raw <- results(
  dds,
  contrast = c("condition", "positive", "negative"),
  alpha = analysis_config$de_padj_cutoff
)
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

res_raw <- res_raw[order(res_raw$padj), ]
res_shrunk <- res_shrunk[order(res_shrunk$padj), ]

n_sig <- sum(res_raw$padj < analysis_config$de_padj_cutoff, na.rm = TRUE)
n_up <- sum(
  res_raw$padj < analysis_config$de_padj_cutoff & res_raw$log2FoldChange > analysis_config$de_lfc_cutoff,
  na.rm = TRUE
)
n_down <- sum(
  res_raw$padj < analysis_config$de_padj_cutoff & res_raw$log2FoldChange < -analysis_config$de_lfc_cutoff,
  na.rm = TRUE
)

message(n_sig, " significant genes (padj < 0.05)")
message("  ", n_up, " upregulated (log2FC > 1)")
message("  ", n_down, " downregulated (log2FC < -1)")

res_raw_df <- as.data.frame(res_raw)
res_raw_df$gene <- rownames(res_raw_df)
res_raw_df <- res_raw_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

res_shrunk_df <- data.frame(
  gene = rownames(res_shrunk),
  baseMean = res_shrunk$baseMean,
  log2FoldChange = res_shrunk$log2FoldChange,
  lfcSE = res_shrunk$lfcSE,
  stat = res_raw[rownames(res_shrunk), "stat"],
  pvalue = res_raw[rownames(res_shrunk), "pvalue"],
  padj = res_raw[rownames(res_shrunk), "padj"],
  stringsAsFactors = FALSE
)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res_raw_df, "results/tables/deseq2_results.csv", row.names = FALSE)
write.csv(res_shrunk_df, "results/tables/deseq2_results_shrunken.csv", row.names = FALSE)
saveRDS(dds, "data/dds_object.rds")
saveRDS(res_raw, "data/deseq_results.rds")
saveRDS(res_shrunk, "data/deseq_results_shrunken.rds")

message("Results saved")
