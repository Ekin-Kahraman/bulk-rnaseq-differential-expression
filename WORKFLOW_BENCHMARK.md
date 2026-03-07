# Workflow Benchmark

This repository starts from a published GEO count matrix, so the benchmark below focuses on differential expression and reproducibility workflow quality, not on FASTQ-level alignment or quantification.

## Reference workflows

- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html): recommends effect-size shrinkage for ranking and visualization.
- [nf-core/rnaseq](https://nf-co.re/rnaseq/latest/): exemplifies a mature end-to-end RNA-seq workflow, especially standardized QC and reproducible execution from raw reads.
- [`targets` user manual](https://books.ropensci.org/targets/): exemplifies dependency-aware skipping and pipeline orchestration.
- [workflowr](https://jdblischak.github.io/workflowr/): exemplifies research provenance via git-aware reporting and session/environment capture.

## Current alignment

- DE effect estimation: raw DESeq2 inference plus `apeglm`-shrunken log2 fold changes for ranking, MA plotting, and volcano visualization.
- Robustness: balanced-subset analysis is benchmarked against the full QC-passed cohort, with overlap and concordance written to `results/tables/analysis_summary.csv`.
- Reproducibility: `renv`-pinned environment, deterministic seeds, GitHub Actions rebuilds, and explicit git/session provenance in `results/session_info.txt`.
- Artifact validation: tracked tables are rebuilt in CI and compared against committed results; key figures are also diff-checked.

## Remaining gaps relative to broader workflows

- Upstream RNA-seq processing: unlike `nf-core/rnaseq`, this repo does not perform raw-read QC, alignment, quantification, or MultiQC because the starting point is the GEO count matrix.
- Pipeline engine: the workflow remains a sequential R-script orchestrator rather than a declarative DAG like `targets`.
- Model complexity: the primary DE model remains `~ condition`; covariates such as age, sex, viral load, or inferred cell composition are not yet modeled explicitly.

## Scope and rationale

- The codebase is intentionally small and reviewable for a focused secondary analysis.
- The robustness layer addresses the main biological risk introduced by the balanced-subset design without requiring a full re-architecture.
- The remaining gaps are explicit, documented, and suitable targets for future extension.
