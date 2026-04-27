escape_regex <- function(value) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", value)
}

is_binary_build_revision_only <- function(status_output) {
  lines <- trimws(status_output)
  mismatch_lines <- grep("^-\\s+\\S+\\s+\\[repo:", lines, value = TRUE)
  if (length(mismatch_lines) == 0) {
    return(FALSE)
  }

  all(vapply(mismatch_lines, function(line) {
    match <- regexec(
      "^-\\s+(\\S+)\\s+\\[repo: CRAN != RSPM; ver: (\\S+) != (\\S+)\\]$",
      line
    )
    parts <- regmatches(line, match)[[1]]
    if (length(parts) != 4) {
      return(FALSE)
    }

    lock_version <- parts[[3]]
    library_version <- parts[[4]]
    grepl(
      paste0("^", escape_regex(lock_version), "-[0-9]+$"),
      library_version
    )
  }, logical(1)))
}

assert_renv_synchronized <- function() {
  status_output <- capture.output(status <- renv::status(), type = "output")
  if (is.list(status) && isTRUE(status$synchronized)) {
    message("renv library is synchronized with renv.lock.")
    return(invisible(TRUE))
  }

  if (is_binary_build_revision_only(status_output)) {
    message(
      "renv library differs only by Posit Package Manager binary build ",
      "revision suffixes; lockfile package versions are otherwise synchronized."
    )
    return(invisible(TRUE))
  }

  if (length(status_output) > 0) {
    message(paste(status_output, collapse = "\n"))
  }
  stop("renv restore completed, but the project library is still not synchronized.")
}
