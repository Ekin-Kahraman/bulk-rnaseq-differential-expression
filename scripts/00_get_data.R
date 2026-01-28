#!/usr/bin/env Rscript
# Download and prepare SARS-CoV-2 RNA-seq data from GEO
# Dataset: GSE152075 (nasopharyngeal swabs, pos/neg)

library(GEOquery)

GEO_ID <- "GSE152075"

# Setup directories
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

message("Downloading ", GEO_ID)

# Download supplementary files
getGEOSuppFiles(GEO_ID, makeDirectory = TRUE, baseDir = "data/raw")

# Extract archives
supp_dir <- file.path("data/raw", GEO_ID)
tar_files <- list.files(supp_dir, pattern = "\\.tar$", full.names = TRUE)
invisible(lapply(tar_files, untar, exdir = supp_dir))

# Find count matrix
files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE)
count_file <- files[grepl("count.*\\.txt", basename(files), ignore.case = TRUE)][1]

if (is.na(count_file)) stop("Count matrix not found")

message("Reading: ", basename(count_file))

# Load counts with gene IDs as rownames
raw_counts <- read.table(count_file, header = TRUE, row.names = 1, 
                         check.names = FALSE, stringsAsFactors = FALSE)
raw_counts <- as.matrix(raw_counts)
storage.mode(raw_counts) <- "integer"
raw_counts <- raw_counts[complete.cases(raw_counts), ]

message(nrow(raw_counts), " genes × ", ncol(raw_counts), " samples")

# Get sample metadata
gse <- getGEO(GEO_ID, GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

# Parse SARS-CoV-2 status from characteristics
positivity <- pheno$characteristics_ch1
condition <- ifelse(
  grepl("positivity:\\s*pos", positivity, ignore.case = TRUE), "positive",
  ifelse(grepl("positivity:\\s*neg", positivity, ignore.case = TRUE), "negative", NA)
)

# Extract sample IDs (POS_### or NEG_###) from titles
sample_ids <- sub(".*\\b(POS_\\d+|NEG_\\d+)\\b.*", "\\1", pheno$title)

# Build metadata - ensure it stays a dataframe
metadata <- data.frame(
  sample_id = sample_ids,
  condition = condition,
  row.names = sample_ids,
  stringsAsFactors = FALSE
)

# Convert condition to factor
metadata$condition <- factor(metadata$condition, levels = c("negative", "positive"))

# Align counts with metadata
common <- intersect(colnames(raw_counts), rownames(metadata))
raw_counts <- raw_counts[, common, drop = FALSE]
metadata <- metadata[common, , drop = FALSE]

# Verify structure
stopifnot(is.data.frame(metadata))
stopifnot("condition" %in% colnames(metadata))

message(sum(metadata$condition == "negative"), " negative, ",
        sum(metadata$condition == "positive"), " positive")

# Save
saveRDS(raw_counts, "data/counts_raw.rds")
saveRDS(metadata, "data/metadata.rds")

message("Data saved")