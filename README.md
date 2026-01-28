# Bulk RNA-seq Differential Expression Analysis (R)

A reproducible, modular bulk RNA-seq differential expression pipeline in R.

The workflow demonstrates data acquisition, quality control, exploratory
analysis, and statistical modelling using DESeq2, with an emphasis on
biological interpretability, statistical validity, and reproducibility.

Designed as a portfolio-grade analysis suitable for academic or industry review.

---

## Dataset

- **Source:** NCBI GEO — [GSE152075](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075)
- **Organism:** *Homo sapiens*
- **Data type:** Bulk RNA-seq (raw counts)
- **Biological context:** SARS-CoV-2 infection status

**Original dataset:**
- 484 samples (430 positive, 54 negative)

**Subset used in this analysis:**

To improve interpretability and reduce confounding heterogeneity,
a balanced subset of samples was selected:

- 30 SARS-CoV-2 positive
- 30 SARS-CoV-2 negative


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

- **04_volcano.R**  
  Volcano plot highlighting key antiviral response genes

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

Principal component analysis was performed on variance-stabilised
expression values.

- **PC1 (31% variance):** primary source of expression variability
- **PC2 (22% variance):** secondary source of variability
- Each point represents a sample, coloured by infection status

Samples show partial separation by SARS-CoV-2 infection status along
PC1, indicating a condition-associated transcriptional signal that
exceeds technical variation but does not dominate all sources of
variance.


---

### Library size QC
<img src="results/figures/qc_library_size.png" width="550">

Each point represents a sample, with total sequencing depth shown on a
log10 scale.

- **y-axis:** total read counts per sample (library size)
- **x-axis:** experimental condition

Library sizes are broadly comparable between SARS-CoV-2 positive and
negative samples, supporting the use of standard DESeq2 normalisation
without additional depth-based correction.


---

### Differential expression (volcano)
<img src="results/figures/volcano_plot.png" width="650">

Each point represents a single gene.

- **x-axis:** log2 fold change (SARS-CoV-2 positive vs negative)  
  Positive values indicate higher expression in infected samples.
- **y-axis:** −log10 adjusted p-value (Benjamini–Hochberg)  
  Higher values indicate stronger statistical evidence for differential expression.

Genes passing the significance threshold (adjusted p-value < 0.05 and
|log2 fold change| > 1) are highlighted.

A small number of genes exhibit both large effect sizes and strong statistical
support, indicating a focused transcriptional response rather than
global transcriptome-wide disruption.


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

