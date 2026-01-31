#!/usr/bin/env Rscript
# Volcano plot visualisation of differential expression results

library(ggplot2)
library(ggrepel)
library(dplyr)

if (!file.exists("results/tables/deseq2_results.csv")) {
  stop("File not found: results/tables/deseq2_results.csv. Run scripts/03_deseq2.R first.")
}

res <- read.csv("results/tables/deseq2_results.csv", stringsAsFactors = FALSE)

res <- res %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    sig = case_when(
      padj >= 0.05 ~ "NS",
      log2FoldChange > 1 ~ "Up",
      log2FoldChange < -1 ~ "Down",
      TRUE ~ "Marginal"
    )
  )

top_genes <- res %>%
  filter(padj < 0.001, abs(log2FoldChange) > 2) %>%
  arrange(padj) %>%
  head(10)

message("Labelling ", nrow(top_genes), " genes: ", paste(top_genes$gene, collapse = ", "))

ggplot(res, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.8) +
  scale_color_manual(
    values = c("Up" = "#e74c3c", "Down" = "#3498db", "Marginal" = "#95a5a6", "NS" = "grey80"),
    breaks = c("Up", "Down"),
    labels = c("Upregulated", "Downregulated")
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.5,
    fontface = "italic",
    max.overlaps = 20,
    box.padding = 0.5,
    segment.color = "grey50"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(
    title = "Differential Expression in SARS-CoV-2 Infection",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/volcano_plot.png", width = 8, height = 7, dpi = 300)

write.csv(
  top_genes %>% dplyr::select(gene, log2FoldChange, padj, baseMean),
  "results/tables/top_genes.csv",
  row.names = FALSE
)

message("Volcano plot saved")