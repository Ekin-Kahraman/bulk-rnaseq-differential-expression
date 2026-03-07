# Workflow Benchmark

This repository starts from a published GEO count matrix, so the benchmark below focuses on **differential-expression and reproducibility workflow quality**, not on FASTQ-level alignment or quantification.

## Reference workflows

- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html): recommends effect-size shrinkage for ranking and visualization.
- [nf-core/rnaseq](https://nf-co.re/rnaseq/latest/): exemplifies top-tier end-to-end RNA-seq workflow engineering, especially standardized QC and reproducible execution from raw reads.
- [`targets` user manual](https://books.ropensci.org/targets/): exemplifies dependency-aware skipping and pipeline orchestration.
- [workflowr](https://jdblischak.github.io/workflowr/): exemplifies research provenance via git-aware reporting and session/environment capture.

## Current alignment

- **DE effect estimation**: raw DESeq2 inference plus `apeglm`-shrunken log2 fold changes for ranking, MA plotting, and volcano visualization.
- **Robustness**: balanced-subset analysis is benchmarked against the full QC-passed cohort, with overlap/concordance written to `results/tables/analysis_summary.csv`.
- **Reproducibility**: `renv`-pinned environment, deterministic seeds, GitHub Actions rebuilds, and explicit git/session provenance in `results/session_info.txt`.
- **Artifact validation**: tracked tables are rebuilt in CI and compared against committed results; key figures are also diff-checked.

## Still narrower than top-tier end-to-end workflows

- **Upstream RNA-seq processing**: unlike `nf-core/rnaseq`, this repo does not perform raw-read QC, alignment, quantification, or MultiQC because the starting point is the GEO count matrix.
- **Pipeline engine**: the workflow is still a sequential R-script orchestrator rather than a declarative DAG like `targets`.
- **Model complexity**: the primary DE model remains `~ condition`; covariates such as age, sex, viral load, or inferred cell composition are not yet modeled explicitly.

## Why this is still a reasonable design

- The codebase is intentionally small and reviewable for a focused secondary analysis.
- The new robustness layer addresses the highest-risk biological weakness without forcing a full re-architecture.
- The remaining gaps are now explicit and documented, which makes future extensions decision-ready rather than hidden assumptions.
