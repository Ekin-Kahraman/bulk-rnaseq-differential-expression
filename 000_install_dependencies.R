#!/usr/bin/env Rscript
# Restore all required packages for the analysis pipeline using renv.

repos <- "https://packagemanager.posit.co/cran/latest"

if (!file.exists("renv.lock")) {
  stop("File not found: renv.lock\nThis project uses renv for reproducible dependencies.")
}

options(renv.consent = TRUE)

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv...")
  install.packages("renv", repos = repos)
}

message("Restoring packages from renv.lock...")
renv::restore(prompt = FALSE)

message("\nDone. Run: source('run_all.R')")
