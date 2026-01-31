#!/usr/bin/env Rscript
# Install all required packages for the analysis pipeline

message("Checking and installing required packages...\n")

# CRAN packages
cran_pkgs <- c("ggplot2", "dplyr", "pheatmap", "RColorBrewer", "ggrepel", "scales")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c("DESeq2", "edgeR", "GEOquery", "clusterProfiler", "org.Hs.eg.db", "enrichplot")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

message("\nAll packages installed. Run: source('run_all.R')")