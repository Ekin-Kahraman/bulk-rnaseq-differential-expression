#!/usr/bin/env Rscript
# Download and prepare SARS-CoV-2 RNA-seq data from GEO
# Dataset: GSE152075 (nasopharyngeal swabs, pos/neg)

GEO_ID <- "GSE152075"

counts_out <- "data/counts_raw.rds"
metadata_out <- "data/metadata.rds"

force_download <- tolower(Sys.getenv("FORCE_DOWNLOAD", unset = "false")) %in% c(
  "1", "true", "t", "yes", "y"
)

if (!force_download && file.exists(counts_out) && file.exists(metadata_out)) {
  message(
    "Found existing data outputs (", counts_out, ", ", metadata_out, "). Skipping download.\n",
    "Set FORCE_DOWNLOAD=true to re-download from GEO."
  )
} else {
  library(GEOquery)

  dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
  dir.create("data", recursive = TRUE, showWarnings = FALSE)

  message("Downloading ", GEO_ID)

  getGEOSuppFiles(GEO_ID, makeDirectory = TRUE, baseDir = "data/raw")

  supp_dir <- file.path("data/raw", GEO_ID)
  tar_files <- list.files(supp_dir, pattern = "\\.tar$", full.names = TRUE)
  invisible(lapply(tar_files, untar, exdir = supp_dir))

  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE)
  count_file <- files[grepl("count.*\\.txt", basename(files), ignore.case = TRUE)][1]

  if (is.na(count_file)) stop("Count matrix not found")

  message("Reading: ", basename(count_file))

  raw_counts <- read.table(count_file, header = TRUE, row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE)
  raw_counts <- as.matrix(raw_counts)
  storage.mode(raw_counts) <- "integer"
  raw_counts <- raw_counts[complete.cases(raw_counts), ]

  message(nrow(raw_counts), " genes × ", ncol(raw_counts), " samples")

  gse <- getGEO(GEO_ID, GSEMatrix = TRUE)
  pheno <- pData(gse[[1]])

  positivity <- pheno$characteristics_ch1
  condition <- ifelse(
    grepl("positivity:\\s*pos", positivity, ignore.case = TRUE), "positive",
    ifelse(grepl("positivity:\\s*neg", positivity, ignore.case = TRUE), "negative", NA)
  )

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
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )

  metadata$condition <- factor(metadata$condition, levels = c("negative", "positive"))

  common <- intersect(colnames(raw_counts), rownames(metadata))
  if (length(common) == 0) {
    stop("No overlapping sample IDs between count matrix and metadata.")
  }
  raw_counts <- raw_counts[, common, drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]

  missing_condition <- is.na(metadata$condition)
  if (any(missing_condition)) {
    warning(
      sprintf(
        "Dropping %d samples with unknown condition labels after alignment.",
        sum(missing_condition)
      )
    )
    metadata <- metadata[!missing_condition, , drop = FALSE]
    raw_counts <- raw_counts[, rownames(metadata), drop = FALSE]
  }

  if (ncol(raw_counts) == 0) {
    stop("No samples left after metadata alignment/condition filtering.")
  }

  stopifnot(is.data.frame(metadata))
  stopifnot("condition" %in% colnames(metadata))

  message(sum(metadata$condition == "negative"), " negative, ",
          sum(metadata$condition == "positive"), " positive")

  saveRDS(raw_counts, counts_out)
  saveRDS(metadata, metadata_out)

  message("Data saved")
}
