library(testthat)

project_dir <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(project_dir, "scripts", "input_validation.R"), local = TRUE)

example_counts <- function() {
  matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(c("GENE1", "GENE2"), c("A", "B")))
}

test_that("metadata follows sample IDs rather than its input row order", {
  counts <- example_counts()
  metadata <- data.frame(condition = c("positive", "negative"), row.names = c("B", "A"))
  aligned <- align_sample_metadata(counts, metadata)
  expect_identical(rownames(aligned), c("A", "B"))
  expect_identical(aligned$condition, c("negative", "positive"))
  expect_identical(counts, example_counts())
})

test_that("unmatched count samples fail instead of silently disappearing", {
  metadata <- data.frame(condition = "negative", row.names = "A")
  expect_error(align_sample_metadata(example_counts(), metadata), "Missing metadata for count samples: B")
})

test_that("single-sample alignment keeps a data frame", {
  counts <- example_counts()[, "A", drop = FALSE]
  metadata <- data.frame(condition = c("positive", "negative"), row.names = c("B", "A"))
  aligned <- align_sample_metadata(counts, metadata)
  expect_s3_class(aligned, "data.frame")
  expect_identical(rownames(aligned), "A")
  expect_identical(aligned$condition, "negative")
})

test_that("malformed counts cannot be silently truncated or dropped", {
  for (bad_value in c(NA_real_, NaN, Inf, -1, 0.25, .Machine$integer.max + 1)) {
    counts <- example_counts()
    counts[1, 1] <- bad_value
    expect_error(validate_count_matrix(counts), "Counts")
  }
  expect_error(validate_count_matrix(matrix("x", 1, 1)), "finite numeric")
  expect_error(validate_count_matrix(matrix(numeric(), 0, 0)), "non-empty")
  expect_silent(validate_count_matrix(example_counts()))
})

test_that("ambiguous sample names are rejected before alignment", {
  counts <- example_counts()
  colnames(counts) <- c("A", "A")
  expect_error(validate_count_matrix(counts), "duplicate sample IDs")
  for (bad_id in c("", "  ", NA_character_)) {
    colnames(counts) <- c("A", bad_id)
    expect_error(validate_count_matrix(counts), "non-empty sample IDs")
  }
})
