# Script: 02_pca.R
# Purpose: PCA on variance-stabilised bulk RNA-seq data

# PCA is used as an exploratory step to check:
# - whether samples separate by biological condition
# - whether any obvious outliers or structure remain after QC
#
# Variance-stabilised data are used so that highly expressed genes
# do not dominate distance-based analyses.

library(DESeq2)
library(ggplot2)

# Load variance-stabilised dataset (generated in 01_qc.R)
vst <- readRDS("data/vst_object.rds")

# Run PCA using the most variable genes
pca_data <- plotPCA(
  vst,
  intgroup = "condition",
  returnData = TRUE
)

# Percentage variance explained by each component
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Plot PCA
p <- ggplot(
  pca_data,
  aes(PC1, PC2, color = condition)
) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  labs(
    title = "PCA of variance-stabilised RNA-seq data",
    color = "Condition"
  ) +
  theme_minimal()

# Save figure
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "results/figures/pca_plot.png",
  plot = p,
  width = 6,
  height = 5
)
