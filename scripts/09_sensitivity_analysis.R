#!/usr/bin/env Rscript
# Compare the balanced subset analysis against the full QC-passed cohort.

source("scripts/config.R", local = TRUE)

library(DESeq2)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)

required_files <- c(
  "data/counts_qc.rds",
  "data/metadata_qc.rds",
  "results/tables/deseq2_results.csv",
  "results/tables/deseq2_results_shrunken.csv",
  "results/tables/go_biological_process.csv",
  "results/tables/kegg_pathways.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

if (!requireNamespace("apeglm", quietly = TRUE)) {
  stop("Package 'apeglm' is required for sensitivity analysis. Run Rscript 000_install_dependencies.R.")
}

counts_qc <- readRDS("data/counts_qc.rds")
metadata_qc <- readRDS("data/metadata_qc.rds")

# Filter to samples with non-NA gender (matching the balanced subset model)
has_gender <- !is.na(metadata_qc$gender)
if (sum(!has_gender) > 0) {
  message("Excluding ", sum(!has_gender), " full-cohort samples with missing gender annotation")
  metadata_qc <- metadata_qc[has_gender, , drop = FALSE]
  counts_qc <- counts_qc[, rownames(metadata_qc), drop = FALSE]
}

balanced_raw <- read.csv("results/tables/deseq2_results.csv", stringsAsFactors = FALSE)
balanced_shrunk <- read.csv("results/tables/deseq2_results_shrunken.csv", stringsAsFactors = FALSE)
go_df <- read.csv("results/tables/go_biological_process.csv", stringsAsFactors = FALSE)
kegg_df <- read.csv("results/tables/kegg_pathways.csv", stringsAsFactors = FALSE)
counts_balanced <- readRDS("data/counts_clean.rds")

keep_genes_full <- rowSums(cpm(counts_qc) >= analysis_config$gene_cpm_cutoff) >= analysis_config$gene_min_samples
counts_full <- counts_qc[keep_genes_full, , drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(counts_full, metadata_qc, design = ~ condition + gender)
dds_full$condition <- relevel(dds_full$condition, ref = "negative")
dds_full <- DESeq(dds_full)

coef_name <- grep("^condition_.*_vs_.*$", resultsNames(dds_full), value = TRUE)
if (length(coef_name) != 1) {
  stop("Expected exactly one condition coefficient, found: ", paste(resultsNames(dds_full), collapse = ", "))
}

full_raw <- results(
  dds_full,
  contrast = c("condition", "positive", "negative"),
  alpha = analysis_config$de_padj_cutoff
)
full_shrunk <- lfcShrink(dds_full, coef = coef_name, type = "apeglm")

full_raw_df <- as.data.frame(full_raw)
full_raw_df$gene <- rownames(full_raw_df)
full_raw_df <- full_raw_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

full_shrunk_df <- as.data.frame(full_shrunk)
full_shrunk_df$gene <- rownames(full_shrunk_df)
full_shrunk_df <- full_shrunk_df[, c("gene", "log2FoldChange")]
names(full_shrunk_df)[2] <- "shrunkenLog2FoldChange"

full_results <- merge(full_raw_df, full_shrunk_df, by = "gene", sort = FALSE)
full_results <- full_results[order(full_results$padj), ]

write.csv(full_results, "results/tables/full_cohort_deseq2_results.csv", row.names = FALSE)
saveRDS(full_raw, "data/full_cohort_deseq_results.rds")
saveRDS(full_shrunk, "data/full_cohort_deseq_results_shrunken.rds")

balanced_threshold <- balanced_raw %>%
  dplyr::filter(
    padj < analysis_config$de_padj_cutoff,
    abs(log2FoldChange) > analysis_config$de_lfc_cutoff
  ) %>%
  pull(gene)

full_threshold <- full_results %>%
  dplyr::filter(
    padj < analysis_config$de_padj_cutoff,
    abs(log2FoldChange) > analysis_config$de_lfc_cutoff
  ) %>%
  pull(gene)

comparison <- balanced_shrunk %>%
  dplyr::select(gene, balancedLog2FoldChange = log2FoldChange, balancedPadj = padj) %>%
  dplyr::inner_join(
    full_results %>%
      dplyr::select(
        gene,
        fullLog2FoldChange = shrunkenLog2FoldChange,
        fullPadj = padj
      ),
    by = "gene"
  ) %>%
  dplyr::mutate(
    balancedSelected = gene %in% balanced_threshold,
    fullSelected = gene %in% full_threshold,
    selectionClass = dplyr::case_when(
      balancedSelected & fullSelected ~ "Shared DE",
      balancedSelected ~ "Balanced only",
      fullSelected ~ "Full cohort only",
      TRUE ~ "Not selected"
    )
  )

shared_de <- comparison %>%
  dplyr::filter(selectionClass == "Shared DE")

sign_concordance <- if (nrow(shared_de) > 0) {
  mean(sign(shared_de$balancedLog2FoldChange) == sign(shared_de$fullLog2FoldChange))
} else {
  NA_real_
}

rho_all <- suppressWarnings(
  cor(
    comparison$balancedLog2FoldChange,
    comparison$fullLog2FoldChange,
    method = "spearman",
    use = "complete.obs"
  )
)

summary_df <- data.frame(
  qc_samples_total = nrow(metadata_qc),
  balanced_samples_total = ncol(counts_balanced),
  balanced_genes_tested = nrow(balanced_raw),
  full_cohort_genes_tested = nrow(full_results),
  balanced_sig_genes = sum(balanced_raw$padj < analysis_config$de_padj_cutoff, na.rm = TRUE),
  balanced_threshold_genes = length(unique(balanced_threshold)),
  full_cohort_sig_genes = sum(full_results$padj < analysis_config$de_padj_cutoff, na.rm = TRUE),
  full_cohort_threshold_genes = length(unique(full_threshold)),
  shared_threshold_genes = sum(comparison$selectionClass == "Shared DE"),
  balanced_only_threshold_genes = sum(comparison$selectionClass == "Balanced only"),
  full_only_threshold_genes = sum(comparison$selectionClass == "Full cohort only"),
  shared_sign_concordance = round(sign_concordance, 4),
  shrunken_lfc_spearman = round(rho_all, 4),
  go_terms_enriched = nrow(go_df),
  kegg_pathways_enriched = nrow(kegg_df),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "results/tables/analysis_summary.csv", row.names = FALSE)

plot_df <- comparison %>%
  dplyr::filter(selectionClass != "Not selected")

label_df <- plot_df %>%
  dplyr::mutate(rank_score = pmin(balancedPadj, fullPadj, na.rm = TRUE)) %>%
  dplyr::arrange(rank_score, desc(abs(balancedLog2FoldChange) + abs(fullLog2FoldChange))) %>%
  head(analysis_config$sensitivity_label_count)

p <- ggplot(plot_df, aes(balancedLog2FoldChange, fullLog2FoldChange, color = selectionClass)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_point(alpha = 0.65, size = 1.9) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3,
    seed = analysis_config$sample_seed,
    max.overlaps = 20,
    box.padding = 0.4,
    segment.color = "grey60",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Shared DE" = "#1f78b4",
      "Balanced only" = "#e67e22",
      "Full cohort only" = "#27ae60"
    )
  ) +
  labs(
    title = "Effect-size robustness: balanced subset vs full cohort",
    subtitle = "Shrunken log2 fold changes for genes selected in at least one analysis",
    x = "Balanced subset shrunken log2 fold change",
    y = "Full cohort shrunken log2 fold change",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave("results/figures/sensitivity_lfc_scatter.png", plot = p, width = 8.5, height = 6.5, dpi = 300)

message(
  "Full cohort DE genes (thresholded): ", length(unique(full_threshold)),
  "; shared with balanced subset: ", sum(comparison$selectionClass == "Shared DE"),
  "; sign concordance: ", round(sign_concordance, 3)
)
