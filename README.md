# Bulk RNA-seq Differential Expression Analysis (R)

Reproducible, modular bulk RNA-seq analysis pipeline in R demonstrating RNA-seq
data acquisition, quality control (QC), exploratory analysis, and differential
expression modelling.

Designed as a portfolio-grade workflow prioritising biological interpretability,
statistical correctness, reproducibility, and explicit analytical decision-making.



### DATASET:

Source: NCBI GEO (GSE152075)  
Organism: Homo sapiens  
Data type: Bulk RNA-seq (raw counts)  
Biological context: SARS-CoV-2 infection status  

Original dataset:

484 samples (430 positive, 54 negative)

Balanced subset used: 

30 SARS-CoV-2 positive  
30 SARS-CoV-2 negative  

Balancing improves QC clarity, PCA interpretability, and statistical comparability,
while reducing unnecessary heterogeneity.



### ANALYSIS WORKFLOW:

scripts/  
00_get_data.R – download GEO data and construct balanced subset  
01_qc.R – library QC, gene filtering, variance stabilisation  
02_pca.R – PCA of variance-stabilised expression data  
03_deseq2.R – differential expression analysis (DESeq2)  
04_pathways.R – pathway enrichment (KEGG / GO)



### QUALITY CONTROL SUMMARY: 

Failed libraries removed; library sizes inspected; CPM-based filtering applied.
Genes retained if ≥1 CPM in ≥10 samples. VST used for exploratory analysis.

After QC:  

14,744 genes  
60 samples (30 positive / 30 negative)



### BIOLOGICAL INTERPRETATION (SUMMARY):

Principal component analysis of variance-stabilised expression data shows
partial separation between SARS-CoV-2–positive and negative samples, indicating
condition-associated transcriptional structure beyond technical variation.

Differential expression analysis using DESeq2 identified ~136 genes with
significant expression changes (FDR < 0.05) between conditions, consistent with
a host transcriptional response to viral infection.

Pathway enrichment analysis (KEGG / GO) is used to contextualise these gene-level
changes at the systems level, highlighting immune- and antiviral-response
processes characteristic of SARS-CoV-2 infection.


![PCA Plot](results/figures/pca_plot.png)

*Principal component analysis of variance-stabilised expression data shows partial
separation between SARS-CoV-2–positive and negative samples, indicating
condition-associated transcriptional structure beyond technical variation.*

![Volcano Plot](results/figures/volcano_plot.png)

*Volcano plot of differential expression results (DESeq2). Orange points indicate
significantly differentially expressed genes (FDR < 0.05), highlighting a strong
host antiviral transcriptional response to SARS-CoV-2 infection.*



### KEY RESULTS (GENE-LEVEL):

A volcano plot of differential expression highlights a small number of genes with
large effect sizes and strong statistical support. The most strongly up-regulated
genes in SARS-CoV-2–positive samples are consistent with a canonical interferon-
mediated antiviral response.

Several of the top-ranked signals correspond to interferon-stimulated genes (ISGs),
such as ISG15, OAS family members, and MX genes, which are known to restrict viral
replication and amplify innate immune signalling.

This pattern indicates that the dominant transcriptional signal in the dataset
reflects host innate immune activation rather than technical or batch-driven effects.




### REQUIREMENTS:

R ≥ 4.2; GEOquery; DESeq2; edgeR; tidyverse; ggplot2; pheatmap



### REPRODUCIBILITY:

source("scripts/00_get_data.R")  
source("scripts/01_qc.R")  
source("scripts/02_pca.R")  
source("scripts/03_deseq2.R")  
source("scripts/04_pathways.R")



### PROJECT INTENT:

Demonstrates correct handling of public RNA-seq data, sound experimental design,
robust QC strategy, and familiarity with standard bulk RNA-seq workflows.
