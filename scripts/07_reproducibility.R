#!/usr/bin/env Rscript
# Record session information for reproducibility

source("scripts/config.R", local = TRUE)

dir.create("results", showWarnings = FALSE)

message("Recording session info...")

git_value <- function(args) {
  out <- tryCatch(
    system2("git", args, stdout = TRUE, stderr = FALSE),
    warning = function(w) character(0),
    error = function(e) character(0)
  )
  if (length(out) == 0) {
    return("unavailable")
  }
  trimws(paste(out, collapse = "\n"))
}

sink("results/session_info.txt")

cat("SARS-CoV-2 Host Transcriptional Response Analysis\n")
cat("Session Information\n")
cat(rep("=", 50), "\n\n", sep = "")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Git commit:", git_value(c("rev-parse", "HEAD")), "\n")
cat("Git branch:", git_value(c("rev-parse", "--abbrev-ref", "HEAD")), "\n")

git_status <- git_value(c("status", "--short"))
cat("Git working tree clean:", if (identical(git_status, "")) "yes" else "no", "\n")
if (!identical(git_status, "")) {
  cat("Git status:\n", git_status, "\n", sep = "")
}
cat("\n")

si <- sessionInfo()
cat("R:", si$R.version$version.string, "\n")
cat("Platform:", si$R.version$platform, "\n")
cat("OS:", si$running, "\n\n")

cat("Analysis configuration:\n")
for (name in names(analysis_config)) {
  cat(sprintf("  %-24s %s\n", name, as.character(analysis_config[[name]])))
}
cat("\n")

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
