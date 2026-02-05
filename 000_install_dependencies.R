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

deps <- renv::dependencies()
pkgs <- sort(unique(deps$Package))
pkgs <- setdiff(pkgs, "renv")

missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  renv::install(missing)
}

message("Updating renv.lock...")
renv::snapshot(prompt = FALSE)

message("\nDone. Run: source('run_all.R')")
