#!/usr/bin/env Rscript
# Restore all required packages for the analysis pipeline using renv.

is_macos <- identical(Sys.info()[["sysname"]], "Darwin")

# On macOS, prefer CRAN binaries to avoid compilation toolchain issues.
# (Package Manager / RSPM often serves source builds for macOS.)
if (is_macos) {
  options(pkgType = "binary")
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  Sys.setenv(RENV_CONFIG_PPM_ENABLED = "FALSE")
}

if (!file.exists("renv.lock")) {
  stop("File not found: renv.lock\nThis project uses renv for reproducible dependencies.")
}

options(renv.consent = TRUE)

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv...")
  install.packages("renv")
}

message("Restoring packages from renv.lock...")
renv::restore(prompt = FALSE)

deps <- renv::dependencies()
pkgs <- sort(unique(deps$Package))
pkgs <- setdiff(pkgs, "renv")

missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  for (pkg in missing) {
    message("\n==> Installing: ", pkg)
    renv::install(
      pkg,
      prompt = FALSE,
      dependencies = c("Depends", "Imports", "LinkingTo")
    )
  }
}

message("Updating renv.lock...")
renv::snapshot(prompt = FALSE)

message("\nDone. Run: source('run_all.R')")
