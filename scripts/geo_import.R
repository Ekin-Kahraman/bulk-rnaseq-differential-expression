# Stage GEO downloads and validate them before replacing a usable local import.
# The entry point and tests source input_validation.R into this environment.

copy_geo_files <- function(files, destination) {
  if (length(files) > 0L && !all(file.copy(files, destination, recursive = TRUE))) {
    stop("Could not copy staged GEO files.")
  }
}

retry_geo_download <- function(label, download, sleep = Sys.sleep) {
  for (attempt in seq_len(3L)) {
    result <- tryCatch(
      list(value = download(attempt)),
      error = function(error) {
        if (attempt == 3L) {
          stop(label, " failed after 3 attempts: ", conditionMessage(error))
        }
        message(label, " attempt ", attempt, " failed: ", conditionMessage(error), "; retrying fresh.")
        sleep(attempt)
        NULL
      }
    )
    if (!is.null(result)) return(result$value)
  }
}

read_geo_counts <- function(directory) {
  tar_files <- list.files(directory, pattern = "\\.tar$", full.names = TRUE)
  invisible(lapply(tar_files, untar, exdir = directory))
  files <- list.files(directory, recursive = TRUE, full.names = TRUE)
  count_files <- files[grepl("count.*\\.txt(\\.gz)?$", basename(files), ignore.case = TRUE)]
  if (length(count_files) != 1L) stop("Expected exactly one GEO count matrix.")

  message("Reading: ", basename(count_files))
  counts <- as.matrix(read.table(count_files, header = TRUE, row.names = 1,
                                  check.names = FALSE, stringsAsFactors = FALSE))
  validate_count_matrix(counts) # nolint: object_usage_linter.
  storage.mode(counts) <- "integer"
  counts
}

publish_geo_import <- function(staged, destinations, rename = file.rename) {
  # Each rename stays on the same filesystem. Restore the previous import if a
  # later move fails, so ordinary I/O errors cannot leave a mixed input pair.
  backups <- rep(NA_character_, length(destinations))
  published <- rep(FALSE, length(destinations))
  complete <- FALSE
  on.exit({
    if (!complete) {
      for (i in rev(seq_along(destinations))) {
        if (published[i]) unlink(destinations[i], recursive = TRUE)
        if (!is.na(backups[i]) && file.exists(backups[i])) {
          if (!file.rename(backups[i], destinations[i])) {
            warning("Could not restore previous input; backup retained at ", backups[i])
          }
        }
      }
    } else {
      unlink(backups[!is.na(backups)], recursive = TRUE)
    }
  }, add = TRUE)

  for (i in seq_along(destinations)) {
    if (file.exists(destinations[i])) {
      backups[i] <- tempfile(".geo-previous-", tmpdir = dirname(destinations[i]))
      if (!rename(destinations[i], backups[i])) stop("Could not preserve previous GEO import.")
    }
    if (!rename(staged[i], destinations[i])) stop("Could not publish validated GEO import.")
    published[i] <- TRUE
  }
  complete <- TRUE
  invisible(NULL)
}

import_geo_data <- function(geo_id, data_dir = "data", force_download = FALSE,
                            get_supp_files = GEOquery::getGEOSuppFiles,
                            get_geo = GEOquery::getGEO, sleep = Sys.sleep) {
  raw_dir <- file.path(data_dir, "raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(raw_dir, geo_id)
  staging <- tempfile(".geo-stage-", tmpdir = data_dir)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  cached_files <- list.files(cache_dir, full.names = TRUE)
  matrix_file <- grepl("_series_matrix.*\\.txt(\\.gz)?$", basename(cached_files))
  counts_result <- retry_geo_download("GEO counts", function(attempt) {
    directory <- file.path(staging, paste0("counts-", attempt))
    dir.create(directory)
    if (!force_download && attempt == 1L) {
      copy_geo_files(cached_files[!matrix_file], directory)
    }
    # GEOquery skips files that already exist; a forced refresh and every retry
    # therefore need a fresh directory, not just another getGEOSuppFiles call.
    get_supp_files(geo_id, makeDirectory = FALSE, baseDir = directory)
    list(counts = read_geo_counts(directory), directory = directory)
  }, sleep = sleep)
  counts <- counts_result$counts

  metadata_result <- retry_geo_download("GEO metadata", function(attempt) {
    directory <- file.path(staging, paste0("metadata-", attempt))
    dir.create(directory)
    if (!force_download && attempt == 1L) {
      copy_geo_files(cached_files[matrix_file], directory)
    }
    gse <- get_geo(geo_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = directory)
    metadata <- parse_geo_metadata(Biobase::pData(gse[[1]]))
    # A parseable but stale metadata cache must not silently drop count samples.
    metadata <- align_sample_metadata(counts, metadata) # nolint: object_usage_linter.
    list(metadata = metadata, directory = directory)
  }, sleep = sleep)
  metadata <- metadata_result$metadata

  missing_condition <- is.na(metadata$condition)
  if (any(missing_condition)) {
    warning(sprintf("Dropping %d samples with unknown condition labels after alignment.", sum(missing_condition)))
    metadata <- metadata[!missing_condition, , drop = FALSE]
    counts <- counts[, rownames(metadata), drop = FALSE]
  }
  if (ncol(counts) == 0L) stop("No samples left after metadata alignment/condition filtering.")

  validated_cache <- file.path(staging, "validated")
  dir.create(validated_cache)
  copy_geo_files(list.files(counts_result$directory, full.names = TRUE), validated_cache)
  copy_geo_files(list.files(metadata_result$directory, full.names = TRUE), validated_cache)
  staged_rds <- file.path(staging, c("counts_raw.rds", "metadata.rds"))
  saveRDS(counts, staged_rds[1])
  saveRDS(metadata, staged_rds[2])
  publish_geo_import(
    c(validated_cache, staged_rds),
    c(cache_dir, file.path(data_dir, c("counts_raw.rds", "metadata.rds")))
  )

  message(nrow(counts), " genes × ", ncol(counts), " samples")
  message(sum(metadata$condition == "negative"), " negative, ",
          sum(metadata$condition == "positive"), " positive")
  message("Covariates extracted: viral_load_ct (", sum(!is.na(metadata$viral_load_ct)),
          " non-NA), age (", sum(!is.na(metadata$age)),
          " non-NA), gender (", sum(!is.na(metadata$gender)), " non-NA)")
  invisible(list(counts = counts, metadata = metadata))
}

parse_geo_metadata <- function(pheno) {
  parse_characteristic <- function(pheno_data, pattern) {
    # Search all characteristics_ch1.* columns for the pattern
    char_cols <- grep("^characteristics_ch1", colnames(pheno_data), value = TRUE)
    values <- rep(NA_character_, nrow(pheno_data))
    for (col in char_cols) {
      matches <- grepl(pattern, pheno_data[[col]], ignore.case = TRUE)
      extracted <- sub(paste0(".*", pattern, "\\s*"), "", pheno_data[[col]][matches],
                       ignore.case = TRUE)
      values[matches] <- trimws(extracted)
    }
    values
  }

  positivity <- parse_characteristic(pheno, "positivity:")
  condition <- ifelse(
    grepl("pos", positivity, ignore.case = TRUE), "positive",
    ifelse(grepl("neg", positivity, ignore.case = TRUE), "negative", NA)
  )

  viral_load_ct <- as.numeric(parse_characteristic(pheno, "n1_ct:"))
  age <- as.integer(parse_characteristic(pheno, "age:"))
  gender <- parse_characteristic(pheno, "gender:")
  batch <- parse_characteristic(pheno, "sequencing_batch:")

  n_na <- sum(is.na(condition))
  if (n_na > 0) {
    warning(sprintf("%d samples have unknown condition; check GEO metadata format", n_na))
  }

  sample_ids <- sub(".*\\b(POS_\\d+|NEG_\\d+)\\b.*", "\\1", pheno$title)
  sample_ids[!grepl("^(POS|NEG)_\\d+$", sample_ids)] <- NA_character_

  # Fall back to GEO accessions if title parsing fails for a subset of records.
  if ("geo_accession" %in% colnames(pheno)) {
    sample_ids[is.na(sample_ids)] <- pheno$geo_accession[is.na(sample_ids)]
  }

  if (any(is.na(sample_ids))) {
    stop("Could not derive sample IDs for all records from GEO metadata.")
  }
  if (any(duplicated(sample_ids))) {
    stop("Duplicate sample IDs detected in GEO metadata.")
  }

  metadata <- data.frame(
    sample_id = sample_ids,
    condition = condition,
    viral_load_ct = viral_load_ct,
    age = age,
    gender = gender,
    batch = batch,
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )

  metadata$condition <- factor(metadata$condition, levels = c("negative", "positive"))
  metadata$gender <- factor(metadata$gender)
  metadata$batch <- factor(metadata$batch)
  metadata
}
