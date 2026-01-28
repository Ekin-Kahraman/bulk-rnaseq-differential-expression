# Script: 04_volcano.R
# Purpose: Visualise differential gene expression in SARS-CoV-2 infection
#
# Rationale:
# A volcano plot summarises two biologically meaningful quantities:
#   1) magnitude of transcriptional change
#   2) confidence in that change after multiple-testing correction
#
# Genes with large fold changes and low adjusted p-values are the most
# plausible drivers of host response to infection.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# Load differential expression results
res <- read.csv(
  "results/tables/deseq2_results.csv",
  row.names = 1
)

res <- res %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

# Classify genes by biological relevance
res <- res %>%
  mutate(
    group = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE                              ~ "Not significant"
    )
  )

# Volcano plot
p <- ggplot(
  res,
  aes(log2FoldChange, -log10(padj), color = group)
) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(
    values = c(
      "Upregulated"     = "#D55E00",
      "Downregulated"   = "#0072B2",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    title = "Differential gene expression in SARS-CoV-2 infection",
    x = "Log2 fold change (Positive vs Negative)",
    y = "-Log10 adjusted p-value",
    color = NULL
  ) +
  theme_minimal(base_size = 12)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "results/figures/volcano_plot.png",
  p,
  width = 7,
  height = 5
)

message("Volcano plot saved to results/figures/")
