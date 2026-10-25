# Validate raw counts before integer conversion and align metadata explicitly.

validate_sample_ids <- function(ids, label) {
  if (is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(label, " must contain non-empty sample IDs.")
  }
  if (anyDuplicated(ids)) {
    stop(label, " contains duplicate sample IDs.")
  }
  invisible(TRUE)
}

validate_count_matrix <- function(counts) {
  if (!is.matrix(counts) || any(dim(counts) == 0L)) {
    stop("Counts must be a non-empty gene-by-sample matrix.")
  }
  if (!is.numeric(counts) || any(!is.finite(counts))) {
    stop("Counts must contain finite numeric values without missing entries.")
  }
  if (any(counts < 0) || any(counts != floor(counts))) {
    stop("Counts must be non-negative whole numbers; do not supply normalised expression.")
  }
  if (any(counts > .Machine$integer.max)) {
    stop("Counts exceed the supported integer range.")
  }
  validate_sample_ids(colnames(counts), "Count columns")
  invisible(TRUE)
}

align_sample_metadata <- function(counts, metadata) {
  validate_count_matrix(counts)
  if (!is.data.frame(metadata)) {
    stop("Sample metadata must be a data frame.")
  }
  validate_sample_ids(rownames(metadata), "Metadata rows")
  missing <- setdiff(colnames(counts), rownames(metadata))
  if (length(missing) > 0L) {
    stop("Missing metadata for count samples: ", paste(missing, collapse = ", "))
  }
  metadata[colnames(counts), , drop = FALSE]
}
