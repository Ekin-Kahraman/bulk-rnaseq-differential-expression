#!/usr/bin/env Rscript
# Restore the pinned project library from renv.lock.

is_macos <- identical(Sys.info()[["sysname"]], "Darwin")

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(renv.consent = TRUE)

if (!file.exists("renv.lock")) {
  stop("File not found: renv.lock\nThis project uses renv for reproducible dependencies.")
}

if (is_macos) {
  binary_repo <- utils::contrib.url(getOption("repos"), type = "binary")
  source_repo <- utils::contrib.url(getOption("repos"), type = "source")
  has_cran_binaries <- !identical(binary_repo, source_repo) && grepl("/bin/", binary_repo)
  options(pkgType = if (has_cran_binaries) "binary" else "source")
}

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv...")
  install.packages("renv")
}

message("Restoring packages from renv.lock...")
renv::restore(prompt = FALSE)

status_output <- capture.output(status <- renv::status(), type = "output")
if (!isTRUE(status$synchronized)) {
  if (length(status_output) > 0) {
    message(paste(status_output, collapse = "\n"))
  }
  stop("renv restore completed, but the project library is still not synchronized.")
}

message("\nDone. Run: Rscript run_all.R")
