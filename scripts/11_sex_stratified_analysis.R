#!/usr/bin/env Rscript
# Novel analysis: Sex-stratified host response
#
# Hypothesis: Males and females mount distinct transcriptional responses to
# SARS-CoV-2, which may contribute to the observed sex disparity in COVID-19
# severity and mortality (males ~1.7x higher mortality, Peckham et al. 2020).
#
# Approach: Interaction model (condition * gender) to identify genes where the
# infection response differs between sexes.

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

# Filter to samples with gender annotation
has_gender <- !is.na(metadata$gender) & metadata$gender %in% c("M", "F")

if (sum(has_gender) < 20) {
  message("Only ", sum(has_gender), " samples with gender data. Skipping sex analysis.")
  write.csv(data.frame(), "results/tables/sex_interaction_de_results.csv", row.names = FALSE)
  quit(save = "no", status = 0)
}

counts_sex <- counts[, has_gender, drop = FALSE]
meta_sex <- metadata[has_gender, , drop = FALSE]
meta_sex$gender <- droplevels(factor(meta_sex$gender, levels = c("F", "M")))

message("Sex-stratified analysis: ", nrow(meta_sex), " samples")
message("  Female: ", sum(meta_sex$gender == "F"),
        " (", sum(meta_sex$gender == "F" & meta_sex$condition == "positive"), " pos, ",
        sum(meta_sex$gender == "F" & meta_sex$condition == "negative"), " neg)")
message("  Male: ", sum(meta_sex$gender == "M"),
        " (", sum(meta_sex$gender == "M" & meta_sex$condition == "positive"), " pos, ",
        sum(meta_sex$gender == "M" & meta_sex$condition == "negative"), " neg)")

# Check minimum group sizes for interaction model
group_counts <- table(meta_sex$condition, meta_sex$gender)
min_group <- min(group_counts)
if (min_group < 5) {
  message("Minimum group size is ", min_group, " (< 5). Interaction model may be underpowered.")
}

# --- Filter low-count genes ---
library(edgeR)
keep_genes <- rowSums(cpm(counts_sex) >= analysis_config$gene_cpm_cutoff) >=
  min(analysis_config$gene_min_samples, ncol(counts_sex) / 2)
counts_sex <- counts_sex[keep_genes, , drop = FALSE]
message(nrow(counts_sex), " genes retained after filtering")

# --- Interaction model: condition * gender ---
dds_sex <- DESeqDataSetFromMatrix(counts_sex, meta_sex, design = ~ condition * gender)
dds_sex <- DESeq(dds_sex)

message("Model terms: ", paste(resultsNames(dds_sex), collapse = ", "))

# Extract interaction term (genes where M vs F response to infection differs)
interaction_coef <- grep("condition.*gender|gender.*condition", resultsNames(dds_sex), value = TRUE)

if (length(interaction_coef) == 0) {
  message("No interaction coefficient found. Falling back to sex-only comparison.")
  res_sex <- results(dds_sex, name = "gender_M_vs_F", alpha = analysis_config$de_padj_cutoff)
  comparison_label <- "Male vs Female (main effect)"
} else {
  res_sex <- results(dds_sex, name = interaction_coef[1], alpha = analysis_config$de_padj_cutoff)
  comparison_label <- "Condition x Sex Interaction"
}

res_sex <- res_sex[order(res_sex$padj), ]

n_sig <- sum(res_sex$padj < analysis_config$de_padj_cutoff, na.rm = TRUE)
n_up <- sum(res_sex$padj < analysis_config$de_padj_cutoff &
              res_sex$log2FoldChange > analysis_config$de_lfc_cutoff, na.rm = TRUE)
n_down <- sum(res_sex$padj < analysis_config$de_padj_cutoff &
                res_sex$log2FoldChange < -analysis_config$de_lfc_cutoff, na.rm = TRUE)

message("\n", comparison_label, " results:")
message("  ", n_sig, " significant genes (FDR < 0.05)")
message("  ", n_up, " with |log2FC| > 1 favoring male-specific response")
message("  ", n_down, " with |log2FC| > 1 favoring female-specific response")

# --- Save results ---
res_df <- as.data.frame(res_sex)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.csv(res_df, "results/tables/sex_interaction_de_results.csv", row.names = FALSE)

# --- Volcano plot ---
library(ggrepel)

plot_df <- res_df %>%
  mutate(
    significance = case_when(
      is.na(padj) ~ "NS",
      padj < analysis_config$de_padj_cutoff &
        log2FoldChange > analysis_config$de_lfc_cutoff ~ "Male-biased",
      padj < analysis_config$de_padj_cutoff &
        log2FoldChange < -analysis_config$de_lfc_cutoff ~ "Female-biased",
      TRUE ~ "NS"
    ),
    label = ifelse(
      !is.na(padj) & padj < 0.01 & abs(log2FoldChange) > 1,
      gene, NA_character_
    )
  )

set.seed(analysis_config$sample_seed)
p_volcano <- ggplot(plot_df, aes(log2FoldChange, -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 15, na.rm = TRUE) +
  scale_color_manual(values = c("Male-biased" = "#2166ac",
                                "Female-biased" = "#b2182b",
                                "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
  labs(title = paste0("Sex-Stratified DE: ", comparison_label),
       subtitle = paste0(n_sig, " genes with sex-differential response (FDR < 0.05)"),
       x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~"p-value"),
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")
ggsave("results/figures/sex_interaction_volcano.png", p_volcano, width = 8, height = 7, dpi = 300)

# --- PCA colored by sex and condition ---
vsd <- readRDS("data/vst_data.rds")
vsd_sex <- vsd[, colnames(vsd) %in% rownames(meta_sex)]

pca_data <- prcomp(t(assay(vsd_sex)), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  condition = meta_sex[rownames(pca_data$x), "condition"],
  gender = meta_sex[rownames(pca_data$x), "gender"]
)
var_explained <- round(100 * summary(pca_data)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = gender)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("negative" = "#4daf4a", "positive" = "#e41a1c")) +
  labs(title = "PCA by Infection Status and Sex",
       x = paste0("PC1 (", var_explained[1], "% variance)"),
       y = paste0("PC2 (", var_explained[2], "% variance)"),
       color = "Condition", shape = "Sex") +
  theme_classic(base_size = 12)
ggsave("results/figures/pca_sex_stratified.png", p_pca, width = 8, height = 6, dpi = 300)

# --- Summary ---
sex_summary <- data.frame(
  total_samples = nrow(meta_sex),
  female_n = sum(meta_sex$gender == "F"),
  male_n = sum(meta_sex$gender == "M"),
  comparison = comparison_label,
  genes_tested = nrow(res_sex),
  sig_genes_fdr05 = n_sig,
  male_biased = n_up,
  female_biased = n_down,
  stringsAsFactors = FALSE
)
write.csv(sex_summary, "results/tables/sex_analysis_summary.csv", row.names = FALSE)

message("\nSex-stratified analysis complete")
