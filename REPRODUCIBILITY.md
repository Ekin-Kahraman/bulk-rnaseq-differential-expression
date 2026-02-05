# Reproducibility Guide

This repository is set up so a reviewer can reproduce the analysis with a small number of commands.

## What’s Included
- **Code**: analysis scripts in `scripts/` and an orchestrator in `run_all.R`.
- **Version pinning**: `renv.lock` pins CRAN + Bioconductor package versions.
- **Precomputed outputs**: key figures and tables are committed under `results/` for convenience and fast verification.

## From a Clean Checkout (Recommended)
Run these commands from the repository root:

```sh
# Restore/install pinned dependencies (CRAN + Bioconductor)
Rscript 000_install_dependencies.R

# Run the full pipeline (downloads data if needed, then regenerates results)
Rscript run_all.R
```

## Data Download Behavior
The data download step (`scripts/00_get_data.R`) is **idempotent**:
- If `data/counts_raw.rds` and `data/metadata.rds` already exist, it will **skip** re-downloading.
- To force a fresh download from GEO:

```sh
FORCE_DOWNLOAD=true Rscript scripts/00_get_data.R
```

## Network Dependencies
Some steps require network access:
- GEO download (via `GEOquery`) in `scripts/00_get_data.R`
- KEGG pathway annotation (via KEGG REST) in `scripts/06_enrichment.R`

If you are running in a restricted environment, these steps may fail until network access is available.

## Determinism
The balanced subset selection uses a fixed seed (`set.seed(123)` in `scripts/01_qc.R`) so repeated runs should yield the same subset and downstream results, given the same package versions.

## Expected Outputs
After a successful run, you should see (among others):
- `results/tables/deseq2_results.csv`
- `results/figures/volcano_plot.png`
- `results/figures/pca_plot.png`
- `results/session_info.txt` (records R + package versions for the run)

## Verification Commands
```sh
# Check environment consistency against renv.lock
Rscript -e 'renv::status()'

# Run smoke tests (validates presence/shape of committed artifacts)
Rscript -e 'testthat::test_dir("tests/testthat")'

# Lint the analysis scripts
Rscript dev/lint.R
```
