#!/usr/bin/env Rscript
# Restore all required packages for the analysis pipeline using renv.

is_macos <- identical(Sys.info()[["sysname"]], "Darwin")
deps_fields <- c("Depends", "Imports", "LinkingTo")

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
  install.packages("renv", dependencies = deps_fields)
}

message("Restoring packages from renv.lock...")
renv::restore(prompt = FALSE)

if (is_macos && !requireNamespace("ggiraph", quietly = TRUE)) {
  message("\nPre-installing ggiraph (CRAN binary) to avoid local compilation...")
  tryCatch(
    install.packages("ggiraph", type = "binary", dependencies = deps_fields),
    error = function(e) {
      message("ggiraph binary install failed: ", conditionMessage(e))
      message("If ggiraph is required and no binary is available, install Xcode CLT with:")
      message("  xcode-select --install")
    }
  )
}

deps <- renv::dependencies()
pkgs <- sort(unique(deps$Package))
pkgs <- setdiff(pkgs, "renv")

missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  for (pkg in missing) {
    message("\n==> Installing: ", pkg)
    tryCatch(
      renv::install(
        pkg,
        prompt = FALSE,
        type = if (is_macos) "binary" else NULL,
        dependencies = deps_fields
      ),
      error = function(e) {
        message("\nFailed installing ", pkg, ": ", conditionMessage(e))
        if (is_macos) {
          message("\nmacOS troubleshooting:")
          message("  - Ensure Xcode Command Line Tools are installed: xcode-select --install")
          message("  - Then retry: Rscript 000_install_dependencies.R")
        }
        stop(e)
      }
    )
  }
}

message("Updating renv.lock...")
renv::snapshot(prompt = FALSE)

message("\nDone. Run: source('run_all.R')")
