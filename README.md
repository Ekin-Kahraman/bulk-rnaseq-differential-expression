# SARS-CoV-2 Differential Expression Analysis

Bulk RNA-seq analysis pipeline identifying host transcriptional responses to SARS-CoV-2 infection using nasopharyngeal swab samples.

## Dataset

**GEO Accession:** GSE152075  
**Organism:** *Homo sapiens*  
**Platform:** Illumina RNA-seq  
**Samples:** 60 (30 positive, 30 negative)  
**Tissue:** Nasopharyngeal swabs

## Analysis Pipeline

| Script | Purpose |
|--------|---------|
| `00_get_data.R` | Download raw counts and metadata from GEO |
| `01_qc.R` | Library size QC, gene filtering, variance stabilization |
| `02_pca.R` | Principal component analysis |
| `03_deseq2.R` | Differential expression testing (negative binomial GLM) |
| `04_volcano.R` | Volcano plot with top gene labels |
| `05_model_diagnostics.R` | MA plot, dispersion plot, sample correlation, heatmaps |
| `06_reproducibility.R` | Record package versions and session info |

## Key Results

**2,001 differentially expressed genes** (padj < 0.05):
- 1,099 upregulated in SARS-CoV-2
- 803 downregulated

Top upregulated genes show strong interferon-stimulated gene (ISG) signature:

- **IFIT1/2/3** - Interferon-induced proteins with tetratricopeptide repeats
- **OAS3** - 2'-5'-Oligoadenylate synthetase 3 (antiviral enzyme)
- **CXCL10** - C-X-C motif chemokine ligand 10 (inflammatory)
- **DDX58** (RIG-I) - Viral RNA sensor
- **GBP1** - Guanylate-binding protein 1
- **XAF1** - XIAP-associated factor 1
- **SIGLEC1** - Sialic acid-binding Ig-like lectin 1

These results are consistent with canonical type I/II interferon responses to viral infection.

## Quality Control

All diagnostic plots confirm robust analysis:

1. **MA Plot** - Log fold changes are independent of mean expression
2. **Dispersion Plot** - Gene-wise dispersion estimates fit the mean-dispersion relationship
3. **Sample Correlation** - Samples cluster by infection status, not technical artifacts
4. **PCA** - PC1 captures 33% variance, separating positive/negative samples
5. **Top Genes Heatmap** - Top 50 DE genes show clear condition-specific expression

## Usage

```bash
# Run complete pipeline
Rscript scripts/00_get_data.R
Rscript scripts/01_qc.R
Rscript scripts/02_pca.R
Rscript scripts/03_deseq2.R
Rscript scripts/04_volcano.R
Rscript scripts/05_model_diagnostics.R
Rscript scripts/06_reproducibility.R
```

Or in R:
```r
scripts <- c(
  "00_get_data.R",
  "01_qc.R",
  "02_pca.R",
  "03_deseq2.R",
  "04_volcano.R",
  "05_model_diagnostics.R",
  "06_reproducibility.R"
)

for (script in scripts) {
  message("Running: ", script)
  source(file.path("scripts", script))
}
```

## Output Structure

```
results/
├── figures/
│   ├── qc_library_size.png
│   ├── pca_plot.png
│   ├── pca_variance.png
│   ├── volcano_plot.png
│   ├── ma_plot.png
│   ├── dispersion_plot.png
│   ├── sample_correlation.png
│   └── top_genes_heatmap.png
├── tables/
│   ├── deseq2_results.csv
│   └── top_genes.csv
└── session_info.txt
```

## Statistical Methods

- **Normalization:** DESeq2 median-of-ratios
- **Dispersion estimation:** Gene-wise estimates with empirical Bayes shrinkage
- **Testing:** Wald test for differential expression
- **Multiple testing correction:** Benjamini-Hochberg FDR (padj < 0.05)
- **Variance stabilization:** `varianceStabilizingTransformation()` for exploratory analysis
- **Outlier handling:** Cook's distance filtering with automatic refitting

## Filtering Criteria

- **Sample QC:** Library size > 100,000 reads
- **Gene filtering:** CPM ≥ 1 in at least 10 samples
- **Balanced design:** 30 samples per condition (random sampling from available)

## Requirements

- R ≥ 4.0
- Bioconductor: DESeq2, edgeR, GEOquery
- CRAN: ggplot2, ggrepel, dplyr, pheatmap

Install with:
```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "GEOquery"))
install.packages(c("ggplot2", "ggrepel", "dplyr", "pheatmap"))
```

## Reproducibility

Package versions and system information are recorded in `results/session_info.txt` after running `06_reproducibility.R`.

## Biological Interpretation

The dominant transcriptional signature is upregulation of interferon-stimulated genes (ISGs), consistent with innate antiviral immunity. Key pathways include:

1. **Pattern recognition:** DDX58 (RIG-I) detects viral RNA
2. **Interferon signaling:** IFIT family proteins inhibit translation
3. **Antiviral enzymes:** OAS3 activates RNase L to degrade viral RNA
4. **Immune recruitment:** CXCL10 attracts T cells and NK cells
5. **Monocyte activation:** SIGLEC1 (CD169) marks activated monocytes

This response is protective but can contribute to inflammation in severe COVID-19.

## Author

Ekin Kahraman

## License

Data from GEO accession GSE152075. Analysis code available for educational use.