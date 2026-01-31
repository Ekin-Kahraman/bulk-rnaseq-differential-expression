#!/usr/bin/env Rscript
# PCA of variance-stabilized expression data

library(DESeq2)
library(ggplot2)

if (!file.exists("data/vst_data.rds")) {
  stop("File not found: data/vst_data.rds. Run scripts/01_qc.R first.")
}

vsd <- readRDS("data/vst_data.rds")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(
    values = c("negative" = "#3498db", "positive" = "#e74c3c"),
    labels = c("Negative", "SARS-CoV-2")
  ) +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1 (", pct_var[1], "% variance)"),
    y = paste0("PC2 (", pct_var[2], "% variance)"),
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/pca_plot.png", width = 7, height = 6, dpi = 300)

message("PC1: ", pct_var[1], "%, PC2: ", pct_var[2], "%")