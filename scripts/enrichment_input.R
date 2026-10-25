# Over-representation asks which pathways are over-represented among selected
# genes relative to genes that could have passed that same selection.
select_enrichment_genes <- function(results, padj_cutoff, lfc_cutoff) {
  eligible <- is.finite(results$padj) & is.finite(results$log2FoldChange)
  selected <- eligible & results$padj < padj_cutoff &
    abs(results$log2FoldChange) > lfc_cutoff
  list(
    universe = unique(results$gene[eligible]),
    selected = unique(results$gene[selected])
  )
}

map_enrichment_genes <- function(genes, mapping, id_column) {
  mapping <- mapping[mapping[[id_column]] %in% genes$universe &
    !is.na(mapping$ENTREZID), , drop = FALSE]
  selected <- unique(mapping$ENTREZID[mapping[[id_column]] %in% genes$selected])
  universe <- unique(mapping$ENTREZID)
  if (!length(selected) || !length(universe)) {
    stop("Enrichment requires mapped selected genes and a mapped eligible background.")
  }
  list(selected = selected, universe = universe, mapping = mapping)
}
