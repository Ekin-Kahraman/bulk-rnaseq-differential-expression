# Bulk RNA-seq Differential Expression Analysis (R)

A reproducible, modular bulk RNA-seq analysis pipeline in R, demonstrating
end-to-end handling of public transcriptomic data: from acquisition and quality
control to exploratory analysis and differential expression modelling.

The workflow prioritises **biological interpretability**, **statistical rigour**,
and **reproducibility**, and is designed as a portfolio-grade example of standard
bulk RNA-seq analysis practice.

---

## Dataset

- **Source:** NCBI GEO (GSE152075)  
- **Organism:** *Homo sapiens*  
- **Data type:** Bulk RNA-seq (raw counts)  
- **Biological context:** SARS-CoV-2 infection status  

**Original dataset:**  
- 484 samples (430 positive, 54 negative)

**Subset used for analysis:**  
- 30 SARS-CoV-2 positive  
- 30 SARS-CoV-2 negative  

A balanced subset was selected to improve interpretability of exploratory analyses
(PCA), ensure comparability between conditions, and reduce unnecessary
heterogeneity unrelated to the biological contrast of interest.

---

## Analysis workflow

scripts/
├── 00_get_data.R # GEO download and balanced subset construction
├── 01_qc.R # Library QC, gene filtering, variance stabilisation
├── 02_pca.R # PCA of variance-stabilised expression data
├── 03_deseq2.R # Differential expression analysis (DESeq2)
├── 04_pathways.R # Pathway enrichment (KEGG / GO)



---

## Quality control summary

- Failed libraries removed
- Library sizes inspected across conditions
- CPM-based gene filtering applied
- Genes retained if ≥ 1 CPM in ≥ 10 samples
- Variance stabilising transformation (VST) used for exploratory analyses

**After QC:**  
- 14,744 genes  
- 60 samples (30 positive / 30 negative)

---

## Biological interpretation (summary)

Principal component analysis of variance-stabilised expression data shows partial
separation between SARS-CoV-2–positive and negative samples, indicating that
infection status contributes meaningfully to global transcriptional variance
beyond technical noise.

Differential expression analysis using DESeq2 identified approximately 136 genes
with significant expression changes (FDR < 0.05) between conditions, consistent
with a host transcriptional response to viral infection.

Pathway enrichment analysis (KEGG / GO) contextualises these gene-level changes at
the systems level, highlighting immune- and antiviral-response pathways
characteristic of SARS-CoV-2 infection.

---

## Key results

### PCA of samples
<img src="results/figures/pca_plot.png" width="600">

### Library size QC
<img src="results/figures/qc_library_size.png" width="550">

### Differential expression (volcano)
<img src="results/figures/volcano_plot.png" width="650">

---

## Key results (gene-level)

The volcano plot highlights a small subset of genes with large effect sizes and
strong statistical support. The most strongly up-regulated genes in
SARS-CoV-2–positive samples are consistent with a canonical
interferon-mediated antiviral response.

Several top-ranked signals correspond to well-characterised
interferon-stimulated genes (ISGs), including members of the **ISG15**, **OAS**,
and **MX** gene families, which are known to restrict viral replication and amplify
innate immune signalling.

Overall, the dominant transcriptional signal reflects biologically coherent host
immune activation rather than batch- or quality-driven artefacts.

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

source("scripts/00_get_data.R")
source("scripts/01_qc.R")
source("scripts/02_pca.R")
source("scripts/03_deseq2.R")
source("scripts/04_pathways.R")


### PROJECT INTENT:

This project demonstrates correct handling of public RNA-seq data, sound
experimental design, robust quality control, and familiarity with standard
bulk RNA-seq analytical workflows commonly used in academic and applied
bioinformatics.
