analysis_config <- list(
  sample_seed = 123L,
  samples_per_group = 30L,
  min_library_size = 1e5,
  gene_cpm_cutoff = 1,
  gene_min_samples = 10L,
  de_padj_cutoff = 0.05,
  de_lfc_cutoff = 1,
  volcano_label_padj_cutoff = 0.001,
  volcano_label_lfc_cutoff = 2
)

is_truthy_env <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}
