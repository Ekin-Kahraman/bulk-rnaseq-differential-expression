# Reproducibility Guide

This repository is set up so a reviewer can reproduce the analysis with a small number of commands.

## What's Included
- **Code**: analysis scripts in `scripts/` and an orchestrator in `run_all.R`.
- **Version pinning**: `renv.lock` pins CRAN + Bioconductor package versions.
- **Pre-computed outputs**: key figures and tables are committed under `results/` for convenience and quick verification.
- **Analysis summary**: `results/tables/analysis_summary.csv` captures the main counts used in the narrative.
- **Output manifest**: `results/tables/output_manifest.csv` records file sizes and MD5 checksums for committed figures and tables.
- **Pinned pathway snapshot**: `data/reference/kegg_hsa_pathway_*.tsv` freezes the KEGG human annotations used by enrichment.
- **Explicit enrichment background**: `results/tables/enrichment_background.csv` records the eligible genes, Entrez mappings and selected-gene membership; `enrichment_summary.csv` records each database's effective denominator.

## Enrichment background

The over-representation analysis compares the 1,773 selected genes with all
13,279 genes eligible for the same selection: finite DESeq2 adjusted p-values
and finite log2 fold changes. This includes non-significant genes, but excludes
genes removed by independent filtering or without estimable statistics. It
conditions the pathway question on DESeq2's selection eligibility; it does not
use the entire annotation database or only the significant genes as background.

The selected set and background share one gene-ID mapping. Entrez IDs are
deduplicated and passed explicitly as `gene` and `universe` to both GO and KEGG
enrichment. The databases then retain their annotated genes: 11,520 for GO
Biological Process and 5,699 for KEGG. See the
[clusterProfiler reference manual](https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf)
for the `universe` argument and database intersection.

The 25 October 2026 correction replaces the earlier database-wide background.
It changes the reported enrichment from 397 to 204 GO terms and from 23 to 7
KEGG pathways. The top COVID-19 pathway remains strongly enriched
(adjusted p = 4.13e-32); the 1,773-gene result and antiviral interpretation are
unchanged. Tables, dotplots and checksums were regenerated together. Tests
check selection eligibility, shared mapping and the published denominators.

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

GEO series metadata is cached beside the downloaded counts. `FORCE_DOWNLOAD=true`
bypasses both caches. Downloads are staged separately, validated and then published;
a failed refresh preserves the previous cache and RDS inputs. Counts and metadata
each retry up to three times in fresh directories, so a partial or corrupt cached
file cannot poison later attempts. Metadata retries do not repeat count processing.
Unused platform annotations are not downloaded.

Raw counts must be finite non-negative integers, with unique sample IDs and matching
metadata for every count column. Invalid inputs fail before modelling rather than
silently losing samples or truncating values. Tests exercise these rules and staged
download failures with fake downloaders, without contacting GEO.

## Determinism
The balanced subset selection uses a fixed seed (`set.seed(123)` in `scripts/01_qc.R`) so repeated runs should yield the same subset and downstream results, given the same package versions. Figure label placement for `ggrepel`-based figures is also seeded, and `results/session_info.txt` now records the active git commit, branch, and analysis configuration.

## Expected Outputs
After a successful run, you should see (among others):
- `results/tables/deseq2_results.csv`
- `results/tables/deseq2_results_shrunken.csv`
- `results/tables/full_cohort_deseq2_results.csv`
- `results/tables/analysis_summary.csv`
- `results/tables/output_manifest.csv`
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
