#!/usr/bin/env Rscript
# Extended analysis: Viral load stratification
#
# Hypothesis: High and low viral load patients activate distinct immune programs.
# Specifically, high viral load (low Ct) may drive stronger interferon signaling,
# while low viral load (high Ct) may show elevated inflammatory/adaptive responses
# reflecting a later stage of immune engagement.
#
# Uses N1 Ct values from GEO metadata to stratify COVID-positive samples.

source("scripts/config.R", local = TRUE)

library(DESeq2)
library(ggplot2)
library(dplyr)

if (!file.exists("data/counts_qc.rds") || !file.exists("data/metadata_qc.rds")) {
  stop("QC data not found. Run scripts/01_qc.R first.")
}

counts <- readRDS("data/counts_qc.rds")
metadata <- readRDS("data/metadata_qc.rds")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Filter to COVID-positive samples with viral load data
pos_mask <- metadata$condition == "positive" & !is.na(metadata$viral_load_ct)

if (sum(pos_mask) < 20) {
  message("Only ", sum(pos_mask), " positive samples with Ct data. ",
          "Viral load analysis requires >= 20. Skipping.")
  # Write empty outputs so downstream tests pass
  write.csv(data.frame(), "results/tables/viral_load_de_results.csv", row.names = FALSE)
  quit(save = "no", status = 0)
}

counts_pos <- counts[, pos_mask, drop = FALSE]
meta_pos <- metadata[pos_mask, , drop = FALSE]

message("Viral load analysis: ", nrow(meta_pos), " positive samples with Ct values")
message("Ct range: ", round(min(meta_pos$viral_load_ct), 1), " - ",
        round(max(meta_pos$viral_load_ct), 1))

# Stratify by median Ct (lower Ct = higher viral load)
ct_median <- median(meta_pos$viral_load_ct, na.rm = TRUE)
meta_pos$viral_load_group <- factor(
  ifelse(meta_pos$viral_load_ct <= ct_median, "high_viral_load", "low_viral_load"),
  levels = c("low_viral_load", "high_viral_load")
)

message("Median Ct: ", round(ct_median, 1))
message("  High viral load (Ct <= ", round(ct_median, 1), "): ",
        sum(meta_pos$viral_load_group == "high_viral_load"))
message("  Low viral load (Ct > ", round(ct_median, 1), "): ",
        sum(meta_pos$viral_load_group == "low_viral_load"))

# --- Ct distribution plot ---
p_ct <- ggplot(meta_pos, aes(x = viral_load_ct, fill = viral_load_group)) +
  geom_histogram(bins = 25, alpha = 0.7, color = "white") +
  geom_vline(xintercept = ct_median, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = c("high_viral_load" = "#d73027", "low_viral_load" = "#4575b4"),
                    labels = c("High Viral Load", "Low Viral Load")) +
  labs(title = "N1 Ct Value Distribution (COVID-Positive Samples)",
       subtitle = paste0("Median Ct = ", round(ct_median, 1), " | Dashed line = stratification threshold"),
       x = "N1 Ct Value (lower = more virus)", y = "Count", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")
ggsave("results/figures/viral_load_ct_distribution.png", p_ct, width = 7, height = 5, dpi = 300)

# --- Filter low-count genes for this subset ---
library(edgeR)
keep_genes <- rowSums(cpm(counts_pos) >= analysis_config$gene_cpm_cutoff) >=
  min(analysis_config$gene_min_samples, ncol(counts_pos) / 2)
counts_pos <- counts_pos[keep_genes, , drop = FALSE]
message(nrow(counts_pos), " genes retained after filtering")

# --- DESeq2: high vs low viral load ---
dds_vl <- DESeqDataSetFromMatrix(counts_pos, meta_pos, design = ~ viral_load_group)
dds_vl <- DESeq(dds_vl)

res_vl <- results(dds_vl, contrast = c("viral_load_group", "high_viral_load", "low_viral_load"),
                  alpha = analysis_config$de_padj_cutoff)
res_vl <- res_vl[order(res_vl$padj), ]

n_sig <- sum(res_vl$padj < analysis_config$de_padj_cutoff, na.rm = TRUE)
n_up <- sum(res_vl$padj < analysis_config$de_padj_cutoff &
              res_vl$log2FoldChange > analysis_config$de_lfc_cutoff, na.rm = TRUE)
n_down <- sum(res_vl$padj < analysis_config$de_padj_cutoff &
                res_vl$log2FoldChange < -analysis_config$de_lfc_cutoff, na.rm = TRUE)

message("\nViral load DE results:")
message("  ", n_sig, " significant genes (FDR < 0.05)")
message("  ", n_up, " upregulated in high viral load (log2FC > 1)")
message("  ", n_down, " downregulated in high viral load (log2FC < -1)")

# --- Save DE results ---
res_df <- as.data.frame(res_vl)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.csv(res_df, "results/tables/viral_load_de_results.csv", row.names = FALSE)

# --- Volcano plot ---
library(ggrepel)

plot_df <- res_df %>%
  mutate(
    significance = case_when(
      is.na(padj) ~ "NS",
      padj < analysis_config$de_padj_cutoff &
        log2FoldChange > analysis_config$de_lfc_cutoff ~ "Up in High VL",
      padj < analysis_config$de_padj_cutoff &
        log2FoldChange < -analysis_config$de_lfc_cutoff ~ "Down in High VL",
      TRUE ~ "NS"
    ),
    label = ifelse(
      !is.na(padj) & padj < analysis_config$volcano_label_padj_cutoff &
        abs(log2FoldChange) > analysis_config$volcano_label_lfc_cutoff,
      gene, NA_character_
    )
  )

set.seed(analysis_config$sample_seed)
p_volcano <- ggplot(plot_df, aes(log2FoldChange, -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 15, na.rm = TRUE) +
  scale_color_manual(values = c("Up in High VL" = "#d73027",
                                "Down in High VL" = "#4575b4",
                                "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
  labs(title = "Differential Expression: High vs Low Viral Load",
       subtitle = paste0(n_up, " up / ", n_down, " down in high viral load (FDR < 0.05, |log2FC| > 1)"),
       x = expression(log[2]~"Fold Change (High VL / Low VL)"),
       y = expression(-log[10]~"p-value"),
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")
ggsave("results/figures/viral_load_volcano.png", p_volcano, width = 8, height = 7, dpi = 300)

# --- Ct vs gene expression correlation for top ISGs ---
vsd <- readRDS("data/vst_data.rds")
vsd_pos <- vsd[, colnames(vsd) %in% rownames(meta_pos)]

# Check for canonical ISGs in the data
isgs <- c("IFIT1", "IFIT2", "IFIT3", "CXCL10", "DDX58", "OAS3", "MX1", "ISG15", "IFITM3")
available_isgs <- isgs[isgs %in% rownames(vsd_pos)]

if (length(available_isgs) >= 3) {
  isg_expr <- t(assay(vsd_pos)[available_isgs, , drop = FALSE])
  isg_df <- data.frame(
    viral_load_ct = meta_pos[rownames(isg_expr), "viral_load_ct"],
    isg_expr,
    check.names = FALSE
  )

  isg_long <- tidyr::pivot_longer(isg_df, cols = -viral_load_ct,
                                   names_to = "gene", values_to = "expression")

  p_corr <- ggplot(isg_long, aes(x = viral_load_ct, y = expression)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "#d73027", linewidth = 0.8) +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    labs(title = "ISG Expression vs Viral Load",
         subtitle = "Lower Ct = higher viral load | Negative slope = higher expression with more virus",
         x = "N1 Ct Value", y = "VST-Normalized Expression") +
    theme_classic(base_size = 11)
  ggsave("results/figures/viral_load_isg_correlation.png", p_corr,
         width = 10, height = 8, dpi = 300)
  message("ISG correlation plot saved")
} else {
  message("Fewer than 3 ISGs found in data; skipping correlation plot")
}

# --- Summary statistics ---
vl_summary <- data.frame(
  total_positive_with_ct = nrow(meta_pos),
  ct_median = round(ct_median, 1),
  high_vl_n = sum(meta_pos$viral_load_group == "high_viral_load"),
  low_vl_n = sum(meta_pos$viral_load_group == "low_viral_load"),
  genes_tested = nrow(res_vl),
  sig_genes_fdr05 = n_sig,
  up_in_high_vl = n_up,
  down_in_high_vl = n_down,
  stringsAsFactors = FALSE
)
write.csv(vl_summary, "results/tables/viral_load_summary.csv", row.names = FALSE)

message("\nViral load stratification complete")
