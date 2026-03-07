#!/usr/bin/env Rscript
# Master script: Run complete analysis pipeline

message("\n=== SARS-CoV-2 Host Response Analysis ===\n")

# Check dependencies
required <- c("DESeq2", "edgeR", "GEOquery", "clusterProfiler",
              "org.Hs.eg.db", "enrichplot", "ggplot2", "dplyr", 
              "pheatmap", "RColorBrewer", "ggrepel", "scales", "apeglm")

missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "), 
       "\nRun: Rscript 000_install_dependencies.R")
}

start_time <- Sys.time()

run_script <- function(script_path, step) {
  message("[", step, "] Running ", basename(script_path), "...")
  tryCatch({
    source(script_path, local = new.env())
    message("    Done\n")
  }, error = function(e) {
    stop("Failed at ", script_path, ": ", e$message)
  })
}

run_script("scripts/00_get_data.R", "1/9")
run_script("scripts/01_qc.R", "2/9")
run_script("scripts/02_pca.R", "3/9")
run_script("scripts/03_deseq2.R", "4/9")
run_script("scripts/04_visualisation_volcano.R", "5/9")
run_script("scripts/05_model_diagnostics.R", "6/9")
run_script("scripts/06_enrichment.R", "7/9")
run_script("scripts/09_sensitivity_analysis.R", "8/9")
run_script("scripts/08_pathway_diagram.R", "9/9")


run_script("scripts/07_reproducibility.R", "Session")

runtime <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
message("=== Complete (", runtime, " min) ===")
