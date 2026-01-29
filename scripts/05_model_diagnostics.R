#!/usr/bin/env Rscript
# Model diagnostics and quality control plots

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

dds <- readRDS("data/dds_object.rds")
vsd <- readRDS("data/vst_data.rds")
res_df <- read.csv("results/tables/deseq2_results.csv")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

message("Generating diagnostic plots...")

# MA plot
res <- results(dds, contrast = c("condition", "positive", "negative"), alpha = 0.05)

png("results/figures/ma_plot.png", width = 7, height = 6, units = "in", res = 300)
plotMA(res, ylim = c(-8, 8), main = "MA Plot: Mean Expression vs Log2 Fold Change")
dev.off()

# Dispersion estimates
png("results/figures/dispersion_plot.png", width = 7, height = 6, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# Sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)

condition_labels <- as.character(colData(vsd)$condition)
sample_labels <- paste0(ifelse(condition_labels == "positive", "P", "N"), "_", 1:ncol(vsd))
rownames(sample_dist_mat) <- sample_labels
colnames(sample_dist_mat) <- sample_labels

annot_df <- data.frame(
  Condition = factor(condition_labels, levels = c("negative", "positive")),
  row.names = sample_labels
)

annot_colors <- list(Condition = c(negative = "#3498db", positive = "#e74c3c"))

pheatmap(
  sample_dist_mat,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_col = annot_df,
  annotation_row = annot_df,
  annotation_colors = annot_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
  main = "Sample Distance Matrix",
  filename = "results/figures/sample_distances.png",
  width = 8, height = 7
)

# Top 50 DE genes heatmap
top50 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

heatmap_mat <- assay(vsd)[top50, ]
heatmap_mat <- t(scale(t(heatmap_mat)))

pheatmap(
  heatmap_mat,
  annotation_col = annot_df,
  annotation_colors = annot_colors,
  show_colnames = FALSE,
  fontsize_row = 6,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  main = "Top 50 DE Genes (Z-score)",
  filename = "results/figures/top50_heatmap.png",
  width = 8, height = 12
)

# PCA scree plot
pca_res <- prcomp(t(assay(vsd)))
var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

scree_df <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  Variance = var_pct[1:10]
)

ggplot(scree_df, aes(PC, Variance)) +
  geom_col(fill = "#2c3e50", alpha = 0.85) +
  geom_text(aes(label = paste0(Variance, "%")), vjust = -0.5, size = 3.5) +
  geom_line(aes(x = as.numeric(PC), y = Variance), color = "#e74c3c", linewidth = 1) +
  geom_point(aes(x = as.numeric(PC), y = Variance), color = "#e74c3c", size = 2) +
  labs(title = "PCA Variance Explained", x = "Principal Component", y = "Variance (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("results/figures/pca_scree.png", width = 7, height = 5, dpi = 300)

message("Diagnostic plots saved")