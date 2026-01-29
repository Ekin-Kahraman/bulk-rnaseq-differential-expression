#!/usr/bin/env Rscript
# Record session information for reproducibility

dir.create("results", showWarnings = FALSE)

message("Recording session info...")

sink("results/session_info.txt")

cat("SARS-CoV-2 Host Transcriptional Response Analysis\n")
cat("Session Information\n")
cat(rep("=", 50), "\n\n", sep = "")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

si <- sessionInfo()

cat("R:", si$R.version$version.string, "\n")
cat("Platform:", si$R.version$platform, "\n")
cat("OS:", si$running, "\n\n")

cat("Attached packages:\n")
for (pkg in names(si$otherPkgs)) {
  cat(sprintf("  %-20s %s\n", pkg, si$otherPkgs[[pkg]]$Version))
}

cat("\nKey packages:\n")
key_pkgs <- c("DESeq2", "edgeR", "GEOquery", "clusterProfiler", 
              "org.Hs.eg.db", "ggplot2", "pheatmap", "dplyr")

for (pkg in key_pkgs) {
  if (pkg %in% rownames(installed.packages())) {
    cat(sprintf("  %-20s %s\n", pkg, as.character(packageVersion(pkg))))
  }
}

if (requireNamespace("BiocManager", quietly = TRUE)) {
  cat("\nBioconductor:", as.character(BiocManager::version()), "\n")
}

sink()

message("Saved to results/session_info.txt")