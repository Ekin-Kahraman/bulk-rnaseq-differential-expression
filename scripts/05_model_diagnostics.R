#!/usr/bin/env Rscript
# Model diagnostics and quality control plots

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(dplyr)

dds <- readRDS("data/dds_object.rds")
vsd <- readRDS("data/vst_data.rds")
res_df <- read.csv("results/tables/deseq2_results.csv")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

message("Generating diagnostic plots...")

# 1. MA plot
png("results/figures/ma_plot.png", width = 7, height = 6, units = "in", res = 300)
plotMA(dds, alpha = 0.05, ylim = c(-8, 8), 
       main = "MA Plot: Mean Expression vs Log2 Fold Change")
dev.off()

# 2. Dispersion plot
png("results/figures/dispersion_plot.png", width = 7, height = 6, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# 3. Sample distances heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)

# Simplified labels
sample_labels <- paste0(as.character(colData(vsd)$condition), "_", 1:ncol(vsd))
rownames(sample_dist_mat) <- sample_labels
colnames(sample_dist_mat) <- sample_labels

# Annotation
annot_df <- data.frame(
  Condition = colData(vsd)$condition,
  row.names = colnames(vsd)
)
rownames(annot_df) <- sample_labels

pheatmap(
  sample_dist_mat,
  annotation_col = annot_df,
  annotation_row = annot_df,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Sample-to-Sample Distance Matrix",
  filename = "results/figures/sample_distances.png",
  width = 8,
  height = 7
)

# 4. Top 50 genes heatmap
top50 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

heatmap_mat <- assay(vsd)[top50, ]
heatmap_mat <- t(scale(t(heatmap_mat)))  # Z-score by row

pheatmap(
  heatmap_mat,
  annotation_col = annot_df,
  show_colnames = FALSE,
  fontsize_row = 7,
  cluster_cols = TRUE,
  main = "Top 50 Differentially Expressed Genes (Z-score)",
  filename = "results/figures/top50_heatmap.png",
  width = 8,
  height = 11
)

# 5. PCA variance scree plot
pca_res <- prcomp(t(assay(vsd)))
var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

scree_df <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  Variance = var_pct[1:10]
)

ggplot(scree_df, aes(PC, Variance)) +
  geom_col(fill = "#2c3e50", alpha = 0.8) +
  geom_text(aes(label = paste0(Variance, "%")), vjust = -0.5, size = 3) +
  labs(
    title = "PCA Scree Plot",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("results/figures/pca_scree.png", width = 7, height = 5, dpi = 300)

message("Diagnostic plots complete:")
message("  - MA plot")
message("  - Dispersion plot")
message("  - Sample distance heatmap")
message("  - Top 50 genes heatmap")
message("  - PCA scree plot")