# SARS-CoV-2 Host Transcriptional Response Analysis

Differential expression analysis of human nasopharyngeal RNA-seq data demonstrating reproducible bulk RNA-seq workflows, statistical modelling, and biological interpretation of antiviral host responses.

---

## Table of Contents

* [Skills Demonstrated](#skills-demonstrated)
* [Dataset](#dataset)
* [Key Results](#key-results)
* [Analysis Summary](#analysis-summary)
* [Project Structure](#project-structure)
* [Workflow](#workflow)
* [Methods](#methods)
* [Biological Interpretation](#biological-interpretation)
* [Limitations](#limitations)
* [Software Requirements](#software-requirements)
* [Reproducibility](#reproducibility)
* [Citation](#citation)

---

## Skills Demonstrated

- Bulk RNA-seq preprocessing and quality control
- Differential expression modelling with DESeq2
- Dimensionality reduction and visualisation
- Pathway enrichment analysis (GO/KEGG)
- Reproducible, script-based bioinformatics workflows

---

## Dataset

**GEO Accession:** [GSE152075](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075)

**Reference:**  
Lieberman NAP, Peddu V, Xie H, Shrestha L, Huang M-L, Mears MC, *et al.* (2020)  
*In vivo antiviral host transcriptional response to SARS-CoV-2 by viral load, sex, and age*  
**PLoS Biology** 18(9): e3000849  
DOI: [10.1371/journal.pbio.3000849](https://doi.org/10.1371/journal.pbio.3000849)

**Platform:** Illumina NovaSeq 6000  
**Organism:** *Homo sapiens*  
**Sample type:** Nasopharyngeal swabs  
**Total samples:** 484 (430 SARS-CoV-2 positive, 54 negative)  
**Analysis subset:** 60 (30 positive, 30 negative)

Balanced subset used to control class imbalance and reduce confounding by viral load heterogeneity. Random subsampling was performed using a fixed seed (`set.seed(123)` in `01_qc.R`) to ensure identical results across runs.

---

## Key Results

### Library Size Quality Control

<img src="results/figures/qc_library_size.png" width="550">

Total read counts across samples stratified by infection status. Libraries show comparable depth between groups (median ~2M reads), supporting robust normalisation. High-depth libraries observed in the full dataset were excluded during sample balancing to ensure comparable distributions between groups.

---

### Principal Component Analysis

<img src="results/figures/pca_plot.png" width="580">

Variance-stabilised expression from 60 balanced samples (30 per group). PC1 (33% variance) partially separates infected from control samples, indicating a detectable transcriptional signature. Overlap reflects biological heterogeneity and the known variability of viral load and host response in nasopharyngeal samples.

---

### Differential Expression

<img src="results/figures/volcano_plot.png" width="620">

Genome-wide differential expression between SARS-CoV-2 positive and negative samples. Labelled genes represent the most significant findings (padj < 0.001, |log₂FC| > 2), chosen for visual clarity while statistical testing used |log₂FC| > 1. Results are dominated by interferon-stimulated genes (IFIT1/2/3, OAS3, CXCL10, DDX58) characteristic of antiviral immunity.

**Note on subset analysis:** The balanced 60-sample subset captures the core ISG signature observed in the full 484-sample dataset while controlling for class imbalance. Effect sizes in the balanced analysis are comparable to full-dataset analyses (Lieberman *et al.* 2020), validating the biological interpretation despite reduced statistical power.

---

## Analysis Summary

**Differentially expressed genes (padj < 0.05):** 2,001  
- Upregulated in SARS-CoV-2: 1,099  
- Downregulated in SARS-CoV-2: 803

**Top interferon-stimulated genes ranked by effect size and statistical significance:**

| Gene | Function | log₂FC | Adjusted *P* |
|------|----------|--------|--------------|
| IFIT1 | Interferon-induced translation inhibitor | 3.5 | <10⁻²⁰ |
| IFIT2 | Interferon-induced translation inhibitor | 3.2 | <10⁻¹⁸ |
| IFIT3 | Interferon-induced translation inhibitor | 3.1 | <10⁻¹⁷ |
| OAS3 | 2'-5'-Oligoadenylate synthetase | 3.0 | <10⁻¹⁷ |
| CXCL10 | Chemokine (interferon-γ inducible) | 2.9 | <10⁻²⁰ |
| DDX58 | RIG-I (viral RNA sensor) | 2.8 | <10⁻¹⁹ |
| GBP1 | Guanylate-binding protein | 2.7 | <10⁻¹⁹ |

The signature is consistent with type I/II interferon responses to viral infection.

---

## Project Structure

```
bulk-rnaseq-differential-expression/
├── scripts/                           # Fully modular, executable analysis pipeline
│   ├── 00_get_data.R                  # Download GEO data
│   ├── 01_qc.R                        # Quality control and filtering
│   ├── 02_pca.R                       # Dimensionality reduction
│   ├── 03_deseq2.R                    # Statistical testing
│   ├── 04_visualisation_volcano.R     # Volcano plot
│   ├── 05_model_diagnostics.R         # Model QC (MA, dispersion, heatmaps)
│   ├── 06_enrichment.R                # Pathway analysis
│   ├── 07_reproducibility.R           # Session info
│   └── 08_pathway_diagram.R           # ISG pathway visualisation
├── data/
│   ├── raw/                           # Downloaded GEO files
│   ├── counts_raw.rds
│   ├── counts_clean.rds
│   ├── metadata.rds
│   ├── dds_object.rds
│   └── vst_data.rds
├── results/
│   ├── figures/
│   │   ├── qc_library_size.png
│   │   ├── pca_plot.png
│   │   ├── pca_scree.png
│   │   ├── volcano_plot.png
│   │   ├── ma_plot.png
│   │   ├── dispersion_plot.png
│   │   ├── sample_distances.png
│   │   ├── top50_heatmap.png
│   │   ├── pathway_diagram.png        # ISG signalling cascade
│   │   ├── kegg_dotplot.png
│   │   └── go_dotplot.png
│   ├── tables/
│   │   ├── deseq2_results.csv
│   │   ├── top_genes.csv
│   │   ├── kegg_pathways.csv
│   │   └── go_biological_process.csv
│   └── session_info.txt
└── README.md
```

---

## Workflow

Execute scripts sequentially from the project root:

```r
source("scripts/00_get_data.R")
source("scripts/01_qc.R")
source("scripts/02_pca.R")
source("scripts/03_deseq2.R")
source("scripts/04_visualisation_volcano.R")
source("scripts/05_model_diagnostics.R")
source("scripts/06_enrichment.R")
source("scripts/07_reproducibility.R")
source("scripts/08_pathway_diagram.R")
```

All scripts are deterministic and can be executed end-to-end on a clean system. Sample selection in `01_qc.R` uses `set.seed(123)` to ensure identical subsets across independent runs.

---

## Methods

### Quality Control
- Minimum library size: 100,000 reads
- Gene filtering: CPM ≥ 1 in ≥10 samples
- Final gene count: 14,220 tested
- Sample balancing: Random selection of 30 per group using fixed seed (123), performed prior to normalisation to avoid size-factor bias

### Statistical Analysis
- **Normalisation:** DESeq2 median-of-ratios
- **Dispersion estimation:** Gene-wise maximum likelihood with empirical Bayes shrinkage
- **Hypothesis testing:** Wald test
- **Multiple testing correction:** Benjamini-Hochberg procedure (FDR < 0.05)
- **Effect size threshold:** |log₂FC| > 1 for biological relevance (volcano plot uses |log₂FC| > 2 for visual clarity)
- **Outlier handling:** Cook's distance with iterative refitting

### Exploratory Analysis
- **Transformation:** Variance-stabilising transformation for visualisation
- **Dimensionality reduction:** Principal component analysis on 500 most variable genes (variance calculated on VST-transformed counts)

---

## Biological Interpretation

<img src="results/figures/pathway_diagram.png" width="650">

*Conceptual pathway showing interferon-stimulated gene (ISG) signalling cascade. Viral RNA is detected by RIG-I (DDX58), triggering interferon production and downstream activation of antiviral effectors.*

---

The dominant transcriptional programme reflects innate antiviral immunity driven by interferon-stimulated genes (ISGs):

1. **Viral recognition:** DDX58 (RIG-I) detects cytoplasmic viral RNA
2. **Translation inhibition:** IFIT1/2/3 block viral protein synthesis
3. **RNA degradation:** OAS3 activates RNase L
4. **Immune recruitment:** CXCL10 attracts effector T cells and NK cells
5. **Inflammasome priming:** Upregulation of cytosolic pattern recognition receptors and interferon signalling components

This response is protective during acute infection but may contribute to immunopathology in severe COVID-19.

---

## Limitations

- Analysis restricted to nasopharyngeal tissue and may not generalise to lower respiratory tract responses
- Subsampling reduces power for low-effect genes but improves interpretability and class balance
- Balanced design controls for viral load heterogeneity but limits analysis of dose-response relationships

---

## Software Requirements

### R Packages

**Bioconductor:**
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

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

### System
- R ≥ 4.0
- 8GB RAM recommended
- ~2GB storage for GEO data

---

## Reproducibility

Package versions and system details are recorded in `results/session_info.txt` after running `07_reproducibility.R`.

All random processes (sample selection, permutation-based tests) use fixed seeds to ensure bit-for-bit reproducibility across systems.

---

## Citation

Data analysed from:

> Lieberman NAP, Peddu V, Xie H, Shrestha L, Huang M-L, Mears MC, *et al.* (2020)  
> In vivo antiviral host transcriptional response to SARS-CoV-2 by viral load, sex, and age.  
> *PLoS Biology* 18(9): e3000849.  
> DOI: [10.1371/journal.pbio.3000849](https://doi.org/10.1371/journal.pbio.3000849)

---

## Author

**Ekin Kahraman**  
January 2026

---

## Data Availability

- **Source data:** GEO accession [GSE152075](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075) (public)
- **Analysis code:** Available in this repository
- **Processed results:** `results/tables/deseq2_results.csv`