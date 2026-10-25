library(testthat)

project_dir <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(project_dir, "scripts", "input_validation.R"), local = TRUE)
source(file.path(project_dir, "scripts", "geo_import.R"), local = TRUE)

geo_id <- "GSE152075"
count_name <- paste0(geo_id, "_raw_counts.txt")
matrix_name <- paste0(geo_id, "_series_matrix.txt.gz")

geo_counts <- function(value = 1L) {
  matrix(value + 0:3, nrow = 2, dimnames = list(c("GENE1", "GENE2"), c("POS_1", "NEG_1")))
}

geo_pheno <- function(age = 40L) {
  data.frame(
    title = c("Patient POS_1", "Patient NEG_1"),
    characteristics_ch1 = c("positivity: positive", "positivity: negative"),
    characteristics_ch1.1 = c("n1_ct: 20", "n1_ct: 0"),
    characteristics_ch1.2 = paste0("age: ", age),
    characteristics_ch1.3 = c("gender: M", "gender: F"),
    characteristics_ch1.4 = "sequencing_batch: 1",
    row.names = c("GSM1", "GSM2")
  )
}

write_geo_cache <- function(data_dir, counts = geo_counts(), pheno = geo_pheno()) {
  cache <- file.path(data_dir, "raw", geo_id)
  dir.create(cache, recursive = TRUE)
  write.table(counts, file.path(cache, count_name), quote = FALSE, col.names = NA)
  saveRDS(pheno, file.path(cache, matrix_name))
  saveRDS(counts, file.path(data_dir, "counts_raw.rds"))
  saveRDS(parse_geo_metadata(pheno), file.path(data_dir, "metadata.rds"))
  cache
}

fake_geo_downloaders <- function(counts = geo_counts(10L), pheno = geo_pheno(50L),
                                 on_counts = function(...) NULL, on_metadata = function(...) NULL) {
  state <- new.env(parent = emptyenv())
  state$counts_dirs <- state$metadata_dirs <- character()
  state$counts_cached <- state$metadata_cached <- logical()
  list(
    state = state,
    get_supp_files = function(GEO, makeDirectory, baseDir) {
      stopifnot(identical(GEO, geo_id), !makeDirectory)
      state$counts_dirs <- c(state$counts_dirs, baseDir)
      path <- file.path(baseDir, count_name)
      state$counts_cached <- c(state$counts_cached, file.exists(path))
      on_counts(length(state$counts_dirs), baseDir)
      if (!file.exists(path)) write.table(counts, path, quote = FALSE, col.names = NA)
    },
    get_geo = function(GEO, GSEMatrix, getGPL, destdir) {
      stopifnot(identical(GEO, geo_id), GSEMatrix, !getGPL)
      state$metadata_dirs <- c(state$metadata_dirs, destdir)
      path <- file.path(destdir, matrix_name)
      state$metadata_cached <- c(state$metadata_cached, file.exists(path))
      on_metadata(length(state$metadata_dirs), destdir)
      if (!file.exists(path)) saveRDS(pheno, path)
      list(Biobase::AnnotatedDataFrame(readRDS(path)))
    }
  )
}

run_fake_import <- function(data_dir, fake, force = FALSE) {
  import_geo_data(geo_id, data_dir, force_download = force,
                  get_supp_files = fake$get_supp_files, get_geo = fake$get_geo,
                  sleep = function(...) NULL)
}

test_that("forced refresh bypasses both GEOquery caches and publishes the validated pair", {
  data_dir <- withr::local_tempdir()
  cache <- write_geo_cache(data_dir)
  fake <- fake_geo_downloaders()
  result <- run_fake_import(data_dir, fake, force = TRUE)
  expect_identical(fake$state$counts_cached, FALSE)
  expect_identical(fake$state$metadata_cached, FALSE)
  expect_identical(result$counts, geo_counts(10L))
  expect_identical(result$metadata$age, c(50L, 50L))
  expect_identical(readRDS(file.path(data_dir, "counts_raw.rds")), result$counts)
  expect_identical(readRDS(file.path(data_dir, "metadata.rds")), result$metadata)
  expect_identical(read_geo_counts(cache), result$counts)
  expect_identical(readRDS(file.path(cache, matrix_name)), geo_pheno(50L))
  expect_length(list.files(data_dir, pattern = "^\\.geo-", all.files = TRUE), 0L)
})

test_that("ordinary rebuilds can reuse validated raw caches", {
  data_dir <- withr::local_tempdir()
  write_geo_cache(data_dir)
  fake <- fake_geo_downloaders()
  result <- run_fake_import(data_dir, fake)
  expect_identical(fake$state$counts_cached, TRUE)
  expect_identical(fake$state$metadata_cached, TRUE)
  expect_identical(result$counts, geo_counts())
  expect_identical(result$metadata$age, c(40L, 40L))
})

test_that("a corrupt metadata cache is retried fresh without downloading counts again", {
  data_dir <- withr::local_tempdir()
  cache <- write_geo_cache(data_dir)
  writeLines("invalid matrix", file.path(cache, matrix_name))
  fake <- fake_geo_downloaders(on_metadata = function(attempt, directory) {
    expect_identical(readLines(file.path(cache, matrix_name)), "invalid matrix")
  })
  result <- run_fake_import(data_dir, fake)
  expect_identical(fake$state$metadata_cached, c(TRUE, FALSE))
  expect_length(unique(fake$state$metadata_dirs), 2L)
  expect_length(fake$state$counts_dirs, 1L)
  expect_identical(result$counts, geo_counts())
  expect_identical(result$metadata$age, c(50L, 50L))
})

test_that("stale metadata missing count samples is refreshed before publication", {
  data_dir <- withr::local_tempdir()
  cache <- write_geo_cache(data_dir)
  saveRDS(geo_pheno()[1, , drop = FALSE], file.path(cache, matrix_name))
  fake <- fake_geo_downloaders()
  result <- run_fake_import(data_dir, fake)
  expect_identical(fake$state$metadata_cached, c(TRUE, FALSE))
  expect_identical(rownames(result$metadata), colnames(result$counts))
})

test_that("an interrupted count download cannot poison the next attempt", {
  data_dir <- withr::local_tempdir()
  cache <- write_geo_cache(data_dir)
  fake <- fake_geo_downloaders(on_counts = function(attempt, directory) {
    expect_identical(read_geo_counts(cache), geo_counts())
    if (attempt == 1L) {
      writeLines("partial download", file.path(directory, count_name))
      stop("connection interrupted")
    }
  })
  result <- run_fake_import(data_dir, fake, force = TRUE)
  expect_identical(fake$state$counts_cached, c(FALSE, FALSE))
  expect_length(unique(fake$state$counts_dirs), 2L)
  expect_identical(result$counts, geo_counts(10L))
})

test_that("failed count validation preserves the previous cache and both RDS files", {
  data_dir <- withr::local_tempdir()
  cache <- write_geo_cache(data_dir)
  before <- tools::md5sum(list.files(data_dir, recursive = TRUE, full.names = TRUE))
  fake <- fake_geo_downloaders(counts = geo_counts() + 0.25)
  expect_error(run_fake_import(data_dir, fake, force = TRUE), "GEO counts failed after 3 attempts")
  expect_length(fake$state$counts_dirs, 3L)
  expect_length(fake$state$metadata_dirs, 0L)
  expect_identical(tools::md5sum(names(before)), before)
  expect_identical(read_geo_counts(cache), geo_counts())
})

test_that("exhausted metadata retries preserve all previously usable inputs", {
  data_dir <- withr::local_tempdir()
  write_geo_cache(data_dir)
  before <- tools::md5sum(list.files(data_dir, recursive = TRUE, full.names = TRUE))
  fake <- fake_geo_downloaders(on_metadata = function(attempt, directory) {
    writeLines("partial response", file.path(directory, matrix_name))
    stop("metadata unavailable")
  })
  expect_error(run_fake_import(data_dir, fake, force = TRUE), "GEO metadata failed after 3 attempts")
  expect_identical(fake$state$metadata_cached, c(FALSE, FALSE, FALSE))
  expect_length(fake$state$counts_dirs, 1L)
  expect_identical(tools::md5sum(names(before)), before)
  expect_length(list.files(data_dir, pattern = "^\\.geo-", all.files = TRUE), 0L)
})

test_that("a publish failure restores an already replaced input pair", {
  directory <- withr::local_tempdir()
  staged <- file.path(directory, c("new-counts", "new-metadata"))
  targets <- file.path(directory, c("counts", "metadata"))
  for (path in staged) writeLines("new", path)
  for (path in targets) writeLines("old", path)
  fail_metadata_move <- function(from, to) {
    if (identical(from, staged[2])) return(FALSE)
    file.rename(from, to)
  }
  expect_error(publish_geo_import(staged, targets, rename = fail_metadata_move), "Could not publish")
  expect_identical(readLines(targets[1]), "old")
  expect_identical(readLines(targets[2]), "old")
})
