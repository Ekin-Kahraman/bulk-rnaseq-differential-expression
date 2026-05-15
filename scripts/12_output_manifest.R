#!/usr/bin/env Rscript
# Write a checksum manifest for committed result artefacts.

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

figure_files <- list.files(
  "results/figures",
  pattern = "\\.(png|pdf)$",
  full.names = TRUE
)
table_files <- list.files(
  "results/tables",
  pattern = "\\.csv$",
  full.names = TRUE
)
table_files <- table_files[basename(table_files) != "output_manifest.csv"]

artefacts <- sort(c(figure_files, table_files))
if (length(artefacts) == 0) {
  stop("No result artefacts found under results/figures or results/tables")
}

info <- file.info(artefacts)
manifest <- data.frame(
  path = gsub("\\\\", "/", artefacts),
  bytes = as.numeric(info$size),
  md5 = unname(tools::md5sum(artefacts)),
  stringsAsFactors = FALSE
)

manifest <- manifest[order(manifest$path), ]
write.csv(manifest, "results/tables/output_manifest.csv", row.names = FALSE)
message("Wrote output manifest for ", nrow(manifest), " artefacts")
