# 04_volcano.R
# Volcano plot of DESeq2 results highlighting key antiviral genes

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# ---- Load DE results ----
res <- read.csv(
  "results/tables/deseq2_results.csv",
  row.names = 1
)

res <- res %>%
  mutate(
    gene = rownames(res),
    sig = padj < 0.05 & abs(log2FoldChange) > 1
  )

# ---- Select genes to label (top biological signals) ----
label_genes <- res %>%
  filter(padj < 0.01, abs(log2FoldChange) > 2) %>%
  arrange(padj) %>%
  slice_head(n = 8)

# ---- Volcano plot ----
p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +
  geom_point(
    aes(color = sig),
    alpha = 0.6,
    size = 1.5
  ) +
  scale_color_manual(
    values = c("grey70", "orange"),
    guide = "none"
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  labs(
    title = "Differential expression in SARS-CoV-2 infection",
    x = "Log2 fold change (Positive vs Negative)",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal()

# ---- Save output ----
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "results/figures/volcano_plot.png",
  plot = p,
  width = 7,
  height = 6
)
