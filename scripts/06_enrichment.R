#!/usr/bin/env Rscript
# Functional enrichment analysis: GO and KEGG pathways

source("scripts/config.R", local = TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(enrichplot)

if (!file.exists("results/tables/deseq2_results.csv")) {
  stop("File not found: results/tables/deseq2_results.csv. Run scripts/03_deseq2.R first.")
}

res_df <- read.csv("results/tables/deseq2_results.csv", stringsAsFactors = FALSE)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("Running enrichment analysis...")

allow_kegg_failure <- is_truthy_env("ALLOW_KEGG_FAILURE")

sig_genes <- res_df %>%
  filter(
    padj < analysis_config$de_padj_cutoff,
    abs(log2FoldChange) > analysis_config$de_lfc_cutoff
  ) %>%
  pull(gene)
sig_genes <- unique(sig_genes)

n_up <- sum(
  res_df$padj < analysis_config$de_padj_cutoff & res_df$log2FoldChange > analysis_config$de_lfc_cutoff,
  na.rm = TRUE
)
n_down <- sum(
  res_df$padj < analysis_config$de_padj_cutoff & res_df$log2FoldChange < -analysis_config$de_lfc_cutoff,
  na.rm = TRUE
)

message(length(sig_genes), " DE genes (", n_up, " up, ", n_down, " down)")
if (length(sig_genes) == 0) {
  stop("No DE genes passed padj/log2FC thresholds; enrichment cannot proceed.")
}

gene_map <- tryCatch({
  bitr(sig_genes, fromType = "ENSEMBL",
       toType = c("ENTREZID", "SYMBOL"),
       OrgDb = org.Hs.eg.db)
}, error = function(e) {
  message("ENSEMBL IDs not recognized, trying SYMBOL...")
  bitr(sig_genes, fromType = "SYMBOL",
       toType = c("ENTREZID"),
       OrgDb = org.Hs.eg.db)
})

if (is.null(gene_map) || nrow(gene_map) == 0) {
  stop("Gene ID conversion failed. Check input gene names match ENSEMBL or SYMBOL format.")
}

message("Mapped ", nrow(gene_map), "/", length(sig_genes), " genes")

go_bp <- enrichGO(
  gene = gene_map$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

go_df <- as.data.frame(go_bp)
write.csv(go_df, "results/tables/go_biological_process.csv", row.names = FALSE)
if (nrow(go_df) > 0) {
  message("  ", nrow(go_df), " GO terms enriched")
  
  p1 <- dotplot(go_bp, showCategory = 20, font.size = 9) +
    ggtitle("GO Biological Process") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave("results/figures/go_dotplot.png", p1, width = 10, height = 9, dpi = 300)
}

kegg_error <- NULL
kegg <- tryCatch(
  enrichKEGG(
    gene = gene_map$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  ),
  error = function(e) {
    kegg_error <<- conditionMessage(e)
    NULL
  }
)

kegg_df <- data.frame()
if (!is.null(kegg)) {
  kegg_df <- as.data.frame(kegg)
  if (nrow(kegg_df) > 0) {
    kegg_read <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    kegg_df <- as.data.frame(kegg_read)
    message("  ", nrow(kegg_df), " KEGG pathways enriched")

    p2 <- dotplot(kegg_read, showCategory = 20, font.size = 9) +
      ggtitle("KEGG Pathway Enrichment") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave("results/figures/kegg_dotplot.png", p2, width = 10, height = 9, dpi = 300)
  }
}
write.csv(kegg_df, "results/tables/kegg_pathways.csv", row.names = FALSE)

if (is.null(kegg) && !is.null(kegg_error)) {
  if (allow_kegg_failure) {
    warning("KEGG enrichment unavailable; continuing without KEGG outputs: ", kegg_error)
  } else {
    stop(
      "KEGG enrichment unavailable: ", kegg_error,
      "\nSet ALLOW_KEGG_FAILURE=true to continue without KEGG outputs."
    )
  }
}

if (nrow(go_df) > 0) {
  message("\nTop GO terms:")
  print(head(go_df[, c("Description", "p.adjust", "Count")], 5), row.names = FALSE)
}

if (nrow(kegg_df) > 0) {
  message("\nTop KEGG pathways:")
  print(head(kegg_df[, c("Description", "p.adjust", "Count")], 5), row.names = FALSE)
}

message("\nEnrichment analysis complete")
