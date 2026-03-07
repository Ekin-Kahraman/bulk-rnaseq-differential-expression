#!/usr/bin/env Rscript
# Maintainer helper: refresh renv.lock after intentionally changing dependencies.

if (!file.exists("renv.lock")) {
  stop("File not found: renv.lock")
}

options(renv.consent = TRUE)

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("Package 'renv' is required. Run: Rscript 000_install_dependencies.R")
}

message("Updating renv.lock from the current project library...")
renv::snapshot(prompt = FALSE)

status_output <- capture.output(status <- renv::status(), type = "output")
if (!isTRUE(status$synchronized)) {
  if (length(status_output) > 0) {
    message(paste(status_output, collapse = "\n"))
  }
  stop("renv.lock was updated, but the project is still not synchronized.")
}

message("renv.lock updated successfully.")
