#!/usr/bin/env Rscript
# Restore all required packages for the analysis pipeline using renv.

is_macos <- identical(Sys.info()[["sysname"]], "Darwin")
deps_fields <- c("Depends", "Imports", "LinkingTo")

# Use a single, stable CRAN repo for both renv + non-renv installs.
options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(RENV_CONFIG_PPM_ENABLED = "FALSE")

# Not every R distribution on macOS supports CRAN macOS binaries. Detect whether
# a "binary" repo path exists; if not, fall back to source installs.
binary_repo <- utils::contrib.url(getOption("repos"), type = "binary")
source_repo <- utils::contrib.url(getOption("repos"), type = "source")
has_cran_binaries <- is_macos && !identical(binary_repo, source_repo) && grepl("/bin/", binary_repo)
install_type <- if (has_cran_binaries) "binary" else "source"
options(pkgType = install_type)

if (!file.exists("renv.lock")) {
  stop("File not found: renv.lock\nThis project uses renv for reproducible dependencies.")
}

options(renv.consent = TRUE)

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv...")
  install.packages("renv", dependencies = deps_fields)
}

if (requireNamespace("BiocManager", quietly = TRUE)) {
  options(repos = BiocManager::repositories())
}

message("Restoring packages from renv.lock...")
renv::restore(prompt = FALSE)

if (is_macos && !requireNamespace("ggiraph", quietly = TRUE)) {
  message("\nPre-installing ggiraph (CRAN ", install_type, ")...")
  tryCatch(
    install.packages("ggiraph", type = install_type, dependencies = deps_fields),
    error = function(e) {
      message("ggiraph install failed: ", conditionMessage(e))
      message("If compilation is required, ensure Xcode CLT are installed:")
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
        type = install_type,
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

if (!identical(Sys.getenv("CI"), "true")) {
  message("Updating renv.lock...")
  renv::snapshot(prompt = FALSE)
} else {
  message("CI detected; skipping renv.lock snapshot.")
}

message("\nDone. Run: Rscript run_all.R")
