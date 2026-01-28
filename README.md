# Bulk RNA-seq Differential Expression Analysis (R)

A reproducible, modular bulk RNA-seq analysis pipeline in R demonstrating
data acquisition, quality control (QC), exploratory analysis, and differential
expression modelling.

Designed as a portfolio-grade workflow with emphasis on biological
interpretability, statistical correctness, and reproducibility.

---

## Dataset

- **Source:** NCBI GEO — [GSE152075](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075)
- **Organism:** *Homo sapiens*
- **Data type:** Bulk RNA-seq (raw counts)
- **Biological context:** SARS-CoV-2 infection status

**Original dataset:**
- 484 samples (430 positive, 54 negative)

**Subset used in this analysis:**
- 30 SARS-CoV-2 positive samples  
- 30 SARS-CoV-2 negative samples  

A balanced subset was used to improve interpretability of PCA and
differential expression results while reducing heterogeneity unrelated
to the biological contrast of interest.

---

## Analysis workflow

`scripts/`

- **00_get_data.R**  
  Download GEO data and construct a balanced case–control subset

- **01_qc.R**  
  Library size QC, CPM-based gene filtering, variance stabilisation

- **02_pca.R**  
  PCA of variance-stabilised expression data

- **03_deseq2.R**  
  Differential expression analysis using DESeq2

- **04_pathways.R**  
  Pathway enrichment analysis (KEGG / GO)

---

## Quality control summary

- Failed libraries removed
- Library sizes inspected across conditions
- Genes retained based on CPM filtering
- Variance stabilising transformation used for exploratory analysis only

The final dataset used for differential expression consists of
high-quality samples and genes suitable for count-based modelling.

---

## Key results

### PCA of samples
<img src="results/figures/pca_plot.png" width="600">

Principal component analysis shows partial separation between
SARS-CoV-2–positive and negative samples, indicating a condition-associated
transcriptional signal beyond technical variation.

---

### Library size QC
<img src="results/figures/qc_library_size.png" width="550">

Library sizes are broadly comparable across conditions, supporting the
use of standard normalisation and downstream modelling assumptions.

---

### Differential expression (volcano)
<img src="results/figures/volcano_plot.png" width="650">

A small number of genes show large effect sizes with strong statistical
support, consistent with a focused host transcriptional response.

---

## Gene-level interpretation

The most strongly up-regulated genes in SARS-CoV-2–positive samples are
consistent with a canonical interferon-mediated antiviral response.

Several top-ranked signals correspond to interferon-stimulated genes
(ISGs), including members of the **ISG**, **OAS**, and **MX** gene families,
which are known to restrict viral replication and amplify innate immune
signalling.

This pattern indicates that the dominant signal in the data reflects
host immune activation rather than batch or technical artefacts.

---

## Requirements

- R ≥ 4.2
- GEOquery
- DESeq2
- edgeR
- tidyverse
- ggplot2
- pheatmap

---

## Reproducibility

Run the analysis end-to-end from the project root:

```r
source("scripts/00_get_data.R")
source("scripts/01_qc.R")
source("scripts/02_pca.R")
source("scripts/03_deseq2.R")
source("scripts/04_volcano.R")

