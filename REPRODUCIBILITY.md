# Reproducibility Guide

This repository is set up so a reviewer can reproduce the analysis with a small number of commands.

## What's Included
- **Code**: analysis scripts in `scripts/` and an orchestrator in `run_all.R`.
- **Version pinning**: `renv.lock` pins CRAN + Bioconductor package versions.
- **Pre-computed outputs**: key figures and tables are committed under `results/` for convenience and quick verification.
- **Analysis summary**: `results/tables/analysis_summary.csv` captures the main counts used in the narrative.
- **Pinned pathway snapshot**: `data/reference/kegg_hsa_pathway_*.tsv` freezes the KEGG human pathway universe used by enrichment.

## From a clean checkout (recommended)
Run these commands from the repository root (i.e., a fresh clone):

```sh
# Restore/install pinned dependencies from renv.lock
Rscript 000_install_dependencies.R

# Run the full pipeline (downloads data if needed, then regenerates results)
Rscript run_all.R
```

Maintainers who intentionally change dependencies should refresh the lockfile explicitly:

```sh
Rscript dev/snapshot_lockfile.R
```

## Data Download Behaviour
The data download step (`scripts/00_get_data.R`) is **idempotent**:
- If `data/counts_raw.rds` and `data/metadata.rds` already exist, it will **skip** re-downloading.
- To force a fresh download from GEO:

```sh
FORCE_DOWNLOAD=true Rscript scripts/00_get_data.R
```

## Network dependencies
The GEO download step (`scripts/00_get_data.R`) requires network access on first run. KEGG enrichment does not query live KEGG during routine analysis; it reads the pinned human pathway snapshot in `data/reference/` so exact table comparisons remain meaningful when KEGG changes upstream.

## Determinism
The balanced subset selection uses a fixed seed (`set.seed(123)` in `scripts/01_qc.R`) so repeated runs should yield the same subset and downstream results, given the same package versions. Figure label placement for `ggrepel`-based figures is also seeded, and `results/session_info.txt` now records the active git commit, branch, and analysis configuration.

## Expected Outputs
After a successful run, you should see (among others):
- `results/tables/deseq2_results.csv`
- `results/tables/deseq2_results_shrunken.csv`
- `results/tables/full_cohort_deseq2_results.csv`
- `results/tables/analysis_summary.csv`
- `results/figures/volcano_plot.png`
- `results/figures/sensitivity_lfc_scatter.png`
- `results/figures/pca_plot.png`
- `results/session_info.txt` (records R, package versions, git commit, and config for the run)

## Verification Commands
```sh
# Check environment consistency against renv.lock
Rscript -e 'renv::status()'

# Run output validation tests
Rscript -e 'testthat::test_dir("tests/testthat")'

# Lint the analysis scripts
Rscript dev/lint.R
```

GitHub Actions also performs a clean rebuild of the tracked analysis outputs, compares regenerated tables against the committed versions, and checks that key figures were regenerated successfully.

For a workflow-level comparison against DESeq2, nf-core/rnaseq, `targets`, and `workflowr`, see `WORKFLOW_BENCHMARK.md`.
