# SARS-CoV-2 Host Transcriptional Response Analysis

Differential expression analysis of nasopharyngeal swab RNA-seq data comparing SARS-CoV-2 positive and negative samples.

## Dataset

- **Source:** [GSE152075](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075) (NCBI GEO)
- **Reference:** [Lieberman et al. (2020)](https://doi.org/10.1016/j.cell.2020.10.037) *Cell* 183(4):1016-1027
- **Platform:** Illumina NovaSeq 6000 (Homo sapiens)
- **Sample type:** Nasopharyngeal swabs
- **Original dataset:** 484 samples (430 positive, 54 negative)
- **Analysis subset:** 60 samples (30 per group, balanced)

## Project Structure

```
bulk-rnaseq-differential-expression/
├── scripts/
│   ├── 00_get_data.R              # Download data from GEO
│   ├── 01_qc.R                    # Quality control and filtering
│   ├── 02_pca.R                   # Principal component analysis
│   ├── 03_deseq2.R                # Differential expression testing
│   ├── 04_volcano.R               # Volcano plot visualization
│   ├── 05_model_diagnostics.R     # QC plots (MA, dispersion, heatmaps)
│   ├── 06_enrichment.R            # KEGG/GO pathway analysis
│   └── 07_reproducibility.R       # Session info
├── data/
│   ├── raw/                       # Downloaded GEO files
│   ├── counts_raw.rds             # Raw count matrix
│   ├── counts_clean.rds           # Filtered counts
│   ├── metadata.rds               # Sample information
│   ├── dds_object.rds             # DESeq2 object
│   └── vst_data.rds               # Variance-stabilized data
├── results/
│   ├── figures/                   # All plots
│   │   ├── qc_library_size.png
│   │   ├── pca_plot.png
│   │   ├── pca_scree.png
│   │   ├── volcano_plot.png
│   │   ├── ma_plot.png
│   │   ├── dispersion_plot.png
│   │   ├── sample_distances.png
│   │   ├── top50_heatmap.png
│   │   ├── kegg_dotplot.png
│   │   └── go_dotplot.png
│   ├── tables/                    # Results tables
│   │   ├── deseq2_results.csv
│   │   ├── top_genes.csv
│   │   ├── kegg_pathways.csv
│   │   └── go_biological_process.csv
│   └── session_info.txt           # Reproducibility info
└── README.md
```

## Analysis Workflow

Run scripts in order from project root:

```r
source("scripts/00_get_data.R")
source("scripts/01_qc.R")
source("scripts/02_pca.R")
source("scripts/03_deseq2.R")
source("scripts/04_volcano.R")
source("scripts/05_model_diagnostics.R")
source("scripts/06_enrichment.R")
source("scripts/07_reproducibility.R")
```

Or run all at once:
```r
scripts <- sprintf("%02d", 0:7)
scripts <- c("get_data", "qc", "pca", "deseq2", "volcano", 
             "model_diagnostics", "enrichment", "reproducibility")
for (i in seq_along(scripts)) {
  source(paste0("scripts/", sprintf("%02d", i-1), "_", scripts[i], ".R"))
}
```

## Key Results

### Differential Expression

- **2,001 genes** significantly differentially expressed (FDR < 0.05)
- **1,099 upregulated** in SARS-CoV-2 positive samples
- **803 downregulated** in SARS-CoV-2 positive samples

### Top Upregulated Genes (ISG Signature)

| Gene | Description | log2FC | padj |
|------|-------------|--------|------|
| **IFIT1** | Interferon-induced protein with tetratricopeptide repeats 1 | 3.5 | < 1e-20 |
| **IFIT2** | Interferon-induced protein with tetratricopeptide repeats 2 | 3.2 | < 1e-18 |
| **IFIT3** | Interferon-induced protein with tetratricopeptide repeats 3 | 3.1 | < 1e-17 |
| **OAS3** | 2'-5'-Oligoadenylate synthetase 3 | 3.0 | < 1e-17 |
| **CXCL10** | C-X-C motif chemokine ligand 10 | 2.9 | < 1e-20 |
| **DDX58** | DExD/H-box helicase 58 (RIG-I) | 2.8 | < 1e-19 |
| **GBP1** | Guanylate binding protein 1 | 2.7 | < 1e-19 |

### Enriched Pathways

**KEGG Pathways:**
- RIG-I-like receptor signaling
- Cytokine-cytokine receptor interaction
- Influenza A response
- NOD-like receptor signaling

**GO Biological Process:**
- Type I interferon signaling pathway
- Defense response to virus
- Innate immune response
- Inflammatory response

## Statistical Methods

- **Normalization:** DESeq2 median-of-ratios
- **Dispersion:** Gene-wise maximum likelihood estimates with empirical Bayes shrinkage
- **Testing:** Wald test
- **Multiple testing:** Benjamini-Hochberg FDR correction (α = 0.05)
- **Effect size:** log2 fold change > 1 for biological significance
- **Outlier detection:** Cook's distance with automatic refitting

## Quality Control

### Sample Filtering
- Library size > 100,000 reads
- Balanced design: 30 samples per condition (random subset)

### Gene Filtering
- CPM ≥ 1 in at least 10 samples
- Removed genes with zero variance
- Final: 14,220 genes tested

### Model Diagnostics
- **MA plot:** No intensity-dependent bias
- **Dispersion plot:** Good fit to mean-dispersion relationship
- **Sample clustering:** Clear separation by infection status (PCA, distance matrix)
- **No batch effects:** Samples cluster by biological condition, not technical factors

## Requirements

### R Packages

**Bioconductor:**
```r
BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "GEOquery",
  "clusterProfiler",
  "org.Hs.eg.db"
))
```

**CRAN:**
```r
install.packages(c(
  "ggplot2",
  "ggrepel",
  "dplyr",
  "pheatmap"
))
```

### System Requirements
- R ≥ 4.0
- 8GB RAM minimum
- 2GB disk space for data

## Biological Interpretation

The transcriptional response to SARS-CoV-2 infection is dominated by interferon-stimulated genes (ISGs), indicating activation of innate antiviral immunity:

1. **Viral sensing:** DDX58 (RIG-I) recognizes viral RNA
2. **Interferon signaling:** IFIT1/2/3 inhibit viral translation
3. **Antiviral enzymes:** OAS3 activates RNase L to degrade viral RNA
4. **Immune cell recruitment:** CXCL10 attracts T cells and NK cells
5. **Inflammasome activation:** NLRP3-related genes upregulated

This signature is protective but excessive activation contributes to COVID-19 severity.

## Reproducibility

All package versions and system details are recorded in `results/session_info.txt` after running `07_reproducibility.R`.

## Citation

If using this analysis:

> Analysis of GSE152075 data from Lieberman NAP, et al. (2020). In vivo antiviral host transcriptional response to SARS-CoV-2 by viral load, sex, and age. *Cell* 183(4):1016-1027.e22. doi: 10.1016/j.cell.2020.10.037

## Author

**Ekin Kahraman**  
Bioinformatics Portfolio Project  
January 2026

## License

Code: MIT License  
Data: GEO accession GSE152075 (public domain)