#!/usr/bin/env Rscript
# Download and prepare SARS-CoV-2 RNA-seq data from GEO
# Dataset: GSE152075 (nasopharyngeal swabs, pos/neg)

GEO_ID <- "GSE152075"
source("scripts/input_validation.R", local = TRUE)

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
  source("scripts/geo_import.R", local = TRUE)
  message("Downloading ", GEO_ID)
  import_geo_data(GEO_ID, force_download = force_download)
  message("Data saved")
}
