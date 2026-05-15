#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript dev/compare_tables.R <reference_dir> <generated_dir>")
}

reference_dir <- args[[1]]
generated_dir <- args[[2]]
relative_tolerance <- 1e-3  # 0.1% — cross-platform DESeq2 drift (macOS vs Linux CI)
zero_tolerance <- 1e-12
effect_size_absolute_tolerance <- 1e-6

if (!dir.exists(reference_dir)) {
  stop("Reference directory not found: ", reference_dir)
}
if (!dir.exists(generated_dir)) {
  stop("Generated directory not found: ", generated_dir)
}

list_csv_files <- function(path) {
  files <- sort(list.files(path, pattern = "\\.csv$", full.names = FALSE))
  setdiff(files, "output_manifest.csv")
}

reference_files <- list_csv_files(reference_dir)
generated_files <- list_csv_files(generated_dir)

if (!identical(reference_files, generated_files)) {
  missing_files <- setdiff(reference_files, generated_files)
  extra_files <- setdiff(generated_files, reference_files)
  details <- c(
    if (length(missing_files) > 0) paste("missing:", paste(missing_files, collapse = ", ")),
    if (length(extra_files) > 0) paste("extra:", paste(extra_files, collapse = ", "))
  )
  stop("CSV file set mismatch between reference and generated outputs. ", paste(details, collapse = "; "))
}

first_mismatch <- function(reference_values, generated_values) {
  same <- (reference_values == generated_values) | (is.na(reference_values) & is.na(generated_values))
  same[is.na(same)] <- FALSE
  mismatch <- which(!same)
  if (length(mismatch) == 0) {
    return(NA_integer_)
  }
  mismatch[[1]]
}

format_value <- function(value) {
  if (is.na(value)) {
    return("NA")
  }
  as.character(value)
}

compare_numeric_column <- function(reference_col, generated_col, file_label, col_name) {
  is_probability_column <- grepl("^(pvalue|padj|qvalue|p.adjust)$", col_name, ignore.case = TRUE)
  absolute_tolerance <- if (is_probability_column) 0 else effect_size_absolute_tolerance

  same <- (reference_col == generated_col) | (is.na(reference_col) & is.na(generated_col))
  same[is.na(same)] <- FALSE

  non_missing <- !(is.na(reference_col) | is.na(generated_col))
  needs_check <- non_missing & !same
  if (any(needs_check)) {
    reference_vals <- reference_col[needs_check]
    generated_vals <- generated_col[needs_check]
    scale <- pmax(abs(reference_vals), abs(generated_vals))
    allowed_diff <- ifelse(scale == 0, zero_tolerance, pmax(relative_tolerance * scale, absolute_tolerance))
    same[needs_check] <- abs(reference_vals - generated_vals) <= allowed_diff
  }

  mismatch <- which(!same)
  if (length(mismatch) > 0) {
    mismatch <- mismatch[[1]]
    diff <- abs(reference_col[[mismatch]] - generated_col[[mismatch]])
    scale <- max(abs(reference_col[[mismatch]]), abs(generated_col[[mismatch]]))
    allowed_diff <- if (scale == 0) {
      zero_tolerance
    } else {
      max(relative_tolerance * scale, absolute_tolerance)
    }
    stop(
      sprintf(
        "%s: numeric mismatch in column '%s' at row %d (reference=%s, generated=%s; abs_diff=%s, allowed=%s)",
        file_label,
        col_name,
        mismatch,
        format_value(reference_col[[mismatch]]),
        format_value(generated_col[[mismatch]]),
        format_value(diff),
        format_value(allowed_diff)
      )
    )
  }
}

compare_text_column <- function(reference_col, generated_col, file_label, col_name) {
  reference_chr <- ifelse(is.na(reference_col), NA_character_, as.character(reference_col))
  generated_chr <- ifelse(is.na(generated_col), NA_character_, as.character(generated_col))
  mismatch <- first_mismatch(reference_chr, generated_chr)
  if (!is.na(mismatch)) {
    stop(
      sprintf(
        "%s: value mismatch in column '%s' at row %d (reference=%s, generated=%s)",
        file_label,
        col_name,
        mismatch,
        format_value(reference_chr[[mismatch]]),
        format_value(generated_chr[[mismatch]])
      )
    )
  }
}

for (filename in reference_files) {
  file_label <- paste0("results/tables/", filename)
  reference_df <- read.csv(
    file.path(reference_dir, filename),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  generated_df <- read.csv(
    file.path(generated_dir, filename),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!identical(names(reference_df), names(generated_df))) {
    stop(file_label, ": column mismatch between reference and generated outputs")
  }
  if (nrow(reference_df) != nrow(generated_df) || ncol(reference_df) != ncol(generated_df)) {
    stop(
      file_label,
      ": dimension mismatch (reference=",
      nrow(reference_df),
      "x",
      ncol(reference_df),
      ", generated=",
      nrow(generated_df),
      "x",
      ncol(generated_df),
      ")"
    )
  }

  for (col_name in names(reference_df)) {
    reference_col <- reference_df[[col_name]]
    generated_col <- generated_df[[col_name]]
    if (is.numeric(reference_col) && is.numeric(generated_col)) {
      compare_numeric_column(reference_col, generated_col, file_label, col_name)
    } else {
      compare_text_column(reference_col, generated_col, file_label, col_name)
    }
  }
}

message(
  "Tracked tables match regenerated outputs across ",
  length(reference_files),
  " files within relative tolerance ",
  format(relative_tolerance, scientific = TRUE),
  "."
)
