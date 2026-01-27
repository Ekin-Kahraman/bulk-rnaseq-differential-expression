# ============================================================
# Script: 02_pca.R
# Purpose: PCA of variance-stabilised bulk RNA-seq data
# ============================================================

# Rationale:
# Principal Component Analysis (PCA) is used here as an *unsupervised*
# exploratory step to assess dominant sources of variance in the dataset.
#
# PCA allows us to:
#   - verify that biological condition (SARS-CoV-2 status) contributes
#     meaningfully to global expression variance
#   - detect potential outliers or batch-driven structure
#   - confirm that QC and filtering steps have produced a coherent dataset
#
# PCA is performed on variance-stabilised data because raw counts violate
# the homoscedasticity assumption underlying Euclidean distance-based methods.

library(DESeq2)
library(ggplot2)

# ---- Load inputs ----
# vst_object.rds is generated in 01_qc.R and contains a DESeq2 VST object
vst <- readRDS("data/vst_object.rds")

# ---- PCA computation ----
# plotPCA computes PCA on the top variable genes and returns sample scores
# when returnData = TRUE
pca_data <- plotPCA(
  vst,
  intgroup = "condition",
  returnData = TRUE
)

# Extract proportion of variance explained for axis labelling
percent_var <- round(100 * attr(pca_data, "percentVar"))

# ---- PCA visualisation ----
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

# ---- Save output ----
# Figures are not version-controlled and are regenerated from scripts
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "results/figures/pca_plot.png",
  plot = p,
  width = 6,
  height = 5
)

