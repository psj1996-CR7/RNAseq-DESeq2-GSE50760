# RNAseq-DESeq2-GSE50760
Reproducible RNA-seq analysis of GEO dataset GSE50760 using DESeq2 for size-factor normalization and differential expression, followed by GO and KEGG pathway enrichment with clusterProfiler and publication-quality visualizations.

## Overview

This repository provides a fully reproducible RNA-sequencing (RNA-seq) data analysis workflow for the public Gene Expression Omnibus (GEO) dataset **GSE50760**, which examines transcriptional differences between colorectal cancer (CRC) and normal tissue samples. The pipeline is scripted in R which performs the diffrential gene expression analysis along with downstream analysis.

The workflow begins with raw gene-level count files downloaded directly from GEO and performs preprocessing, normalization, differential expression analysis, visualization, and functional enrichment analysis. 


## Dataset
- **GEO Accession**: GSE50760  
- **Organism**: *Homo sapiens*  
- **Study Design**: Colorectal cancer vs normal tissue  
- **Data Type**: RNA-seq raw gene count files  
- **Source**: NCBI Gene Expression Omnibus (GEO)

## Analysis Workflow

1. Download and extraction of raw GEO RNA-seq count files  
2. Merging of individual sample count tables into a single large count matrix  
3. Integration of Metadata 
4. Normalization and differential expression analysis using DESeq2  
5. Exploratory and result visualization  
6. Functional enrichment analysis (GO and KEGG)  
7. Gene–pathway network visualization  

## Normalization and Differential Expression

Normalization is performed using **DESeq2’s median-of-ratios method**, which corrects for differences in library size and sequencing depth across samples. Differential expression analysis is conducted using a generalized linear model framework with the following design:

Statistical significance is assessed using the Wald test, and multiple hypothesis testing is controlled using the Benjamini–Hochberg false discovery rate (FDR). For visualization and dimensionality reduction, variance stabilizing transformation (VST) is applied to normalized counts.

## Visualization Outputs

The pipeline generates publication-quality visualizations, including:

- Principal Component Analysis (PCA) plot  
- Volcano plot of differentially expressed genes  
- Heatmap of the top 20 most significant genes  
- Boxplot of the most significant gene  
- GO Biological Process enrichment dot plot  
- KEGG pathway enrichment dot plot  
- Gene–pathway network visualization  

## Functional Enrichment Analysis

Genes with an adjusted p-value < 0.05 and an absolute log2 fold change > 1 are used for enrichment analysis. Gene symbols are mapped to Entrez IDs using **org.Hs.eg.db**.

Functional interpretation is performed using **clusterProfiler**, including:

- Gene Ontology (GO) Biological Process enrichment  
- KEGG pathway enrichment  

Results are exported as structured CSV files and visualized using dot plots and gene–pathway network plots.


## Software and Dependencies

### R Environment
- R version ≥ 4.1

### Required Packages

**CRAN**  
- ggplot2  
- dplyr  
- pheatmap  
- ggrepel  
- readr  
- R.utils  

**Bioconductor**  
- DESeq2  
- clusterProfiler  
- org.Hs.eg.db  
- enrichplot  
- DOSE  

Installation commands are provided in the analysis script.


## Reproducibility

The analysis is fully reproducible. Raw data are downloaded programmatically, all parameters are explicitly defined, and results are generated deterministically. Users may optionally provide a custom `metadata.csv` file to override inferred sample annotations.

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/psj1996-CR7/RNAseq-DESeq2-GSE50760.git
   cd RNAseq-DESeq2-GSE50760

