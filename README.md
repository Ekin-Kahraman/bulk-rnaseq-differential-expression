#' 🧬 Bulk RNA-seq Differential Expression Analysis (R)
#'
#' This repository implements a reproducible, modular bulk RNA-seq analysis
#' pipeline in R, demonstrating practical competence in RNA-seq data
#' acquisition, quality control, exploratory analysis, and differential
#' expression modelling.
#'
#' The project is designed as a portfolio-grade bioinformatics workflow,
#' prioritising:
#'   - biological interpretability
#'   - statistical correctness
#'   - reproducibility
#'   - explicit analytical decision-making
#'
#' over unnecessary technical complexity.
#'
#' -------------------------------------------------------------------------
#' 📊 Dataset
#' -------------------------------------------------------------------------
#'
#' Source:
#'   - NCBI GEO — GSE152075
#'
#' Organism:
#'   - Homo sapiens
#'
#' Data type:
#'   - Bulk RNA-seq (raw count matrix)
#'
#' Biological context:
#'   - SARS-CoV-2 infection status
#'
#' Original dataset composition:
#'   - 484 total samples
#'   - 430 SARS-CoV-2 positive
#'   - 54 SARS-CoV-2 negative
#'
#' Due to strong class imbalance, a balanced subset was deliberately constructed
#' to maximise interpretability and reduce unnecessary heterogeneity:
#'
#'   - 30 SARS-CoV-2 positive samples
#'   - 30 SARS-CoV-2 negative samples
#'
#' This design choice improves:
#'   - QC clarity
#'   - PCA interpretability
#'   - statistical comparability
#'   - interview-level explainability
#'
#' All data acquisition, filtering, and subsetting steps are fully reproducible
#' via script.
#'
#' -------------------------------------------------------------------------
#' 🧪 Analysis Workflow
#' -------------------------------------------------------------------------
#'
#' The pipeline is organised into modular scripts, each responsible for a
#' single analytical stage:
#'
#' scripts/
#' ├── 00_get_data.R    # Download GEO data and construct balanced sample subset
#' ├── 01_qc.R          # Library QC, gene filtering, variance stabilisation
#' ├── 02_pca.R         # PCA of variance-stabilised expression data
#' ├── 03_deseq2.R      # Differential expression analysis (DESeq2)
#' ├── 04_pathways.R   # Pathway enrichment analysis (KEGG / GO)
#'
#' Scripts are executed sequentially, with each step consuming the outputs of
#' the previous stage.
#'
#' Derived data objects and figures are not version-controlled and are
#' regenerated from scripts to ensure full reproducibility.
#'
#' -------------------------------------------------------------------------
#' 🔍 Quality Control Summary
#' -------------------------------------------------------------------------
#'
#' Key QC steps include:
#'   - Removal of failed libraries (zero or invalid counts)
#'   - Library size inspection across experimental conditions
#'   - CPM-based low-count gene filtering
#'
#' Gene filtering rule:
#'   - Genes retained if expressed at ≥ 1 CPM in ≥ 10 samples
#'
#' Downstream transformation:
#'   - Variance stabilising transformation (VST) for PCA and exploratory analysis
#'
#' After QC and filtering:
#'   - 14,744 genes retained
#'   - 60 samples (30 positive / 30 negative)
#'
#' The dataset passes QC and is suitable for exploratory analysis and
#' differential expression modelling.
#'
#' -------------------------------------------------------------------------
#' 🛠 Requirements
#' -------------------------------------------------------------------------
#'
#' Key R packages used:
#'   - GEOquery
#'   - DESeq2
#'   - edgeR
#'   - tidyverse
#'   - ggplot2
#'   - pheatmap
#'
#' R version:
#'   - R ≥ 4.2
#'
#' -------------------------------------------------------------------------
#' 🎯 Project Intent
#' -------------------------------------------------------------------------
#'
#' This repository is intended to demonstrate:
#'   - Correct handling of public RNA-seq count data
#'   - Sound experimental design decisions
#'   - Robust QC and filtering strategy
#'   - Familiarity with standard bulk RNA-seq workflows
#'   - Clear communication of analytical rationale
#'
#' It is designed for academic review, bioinformatics internship applications,
#' and technical interviews.
#'
#' -------------------------------------------------------------------------
#' 📌 Reproducibility
#' -------------------------------------------------------------------------
#'
#' To reproduce the analysis from scratch, run the scripts sequentially:
#'
#' source("scripts/00_get_data.R")
#' source("scripts/01_qc.R")
#' source("scripts/02_pca.R")
#' source("scripts/03_deseq2.R")
#' source("scripts/04_pathways.R")
#'
#' -------------------------------------------------------------------------

