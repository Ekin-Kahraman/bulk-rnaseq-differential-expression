#!/usr/bin/env Rscript
# Functional enrichment analysis: GO and KEGG pathways

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(enrichplot)

res_df <- read.csv("results/tables/deseq2_results.csv", stringsAsFactors = FALSE)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("Running enrichment analysis...")

# Get significant genes
sig_genes <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(gene)

n_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)
n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)

message(length(sig_genes), " DE genes (", n_up, " up, ", n_down, " down)")

# Convert to Entrez IDs
# Try ENSEMBL first, then SYMBOL if that fails
gene_map <- tryCatch({
  bitr(sig_genes, fromType = "ENSEMBL", 
       toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
}, error = function(e) {
  message("ENSEMBL IDs not recognized, trying SYMBOL...")
  bitr(sig_genes, fromType = "SYMBOL", 
       toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
})

message("Mapped ", nrow(gene_map), "/", length(sig_genes), " genes")

# GO biological process
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

if (nrow(go_df) > 0) {
  write.csv(go_df, "results/tables/go_biological_process.csv", row.names = FALSE)
  message("  ", nrow(go_df), " GO terms enriched")
  
  p1 <- dotplot(go_bp, showCategory = 20, font.size = 9) +
    ggtitle("GO Biological Process") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave("results/figures/go_dotplot.png", p1, width = 10, height = 9, dpi = 300)
}

# KEGG pathways
kegg <- enrichKEGG(
  gene = gene_map$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

kegg_df <- as.data.frame(kegg)

if (nrow(kegg_df) > 0) {
  kegg_read <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(kegg_read)
  write.csv(kegg_df, "results/tables/kegg_pathways.csv", row.names = FALSE)
  message("  ", nrow(kegg_df), " KEGG pathways enriched")
  
  p2 <- dotplot(kegg_read, showCategory = 20, font.size = 9) +
    ggtitle("KEGG Pathway Enrichment") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave("results/figures/kegg_dotplot.png", p2, width = 10, height = 9, dpi = 300)
}

# Print top results
if (nrow(go_df) > 0) {
  message("\nTop GO terms:")
  print(head(go_df[, c("Description", "p.adjust", "Count")], 5), row.names = FALSE)
}

if (nrow(kegg_df) > 0) {
  message("\nTop KEGG pathways:")
  print(head(kegg_df[, c("Description", "p.adjust", "Count")], 5), row.names = FALSE)
}

message("\nEnrichment analysis complete")