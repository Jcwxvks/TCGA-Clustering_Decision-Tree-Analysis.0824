# TCGA cohort analysis: gene-based clustering and classification

This repository demonstrates the workflow of distinguishing different sample types in the provided TCGA cohort using expression patterns for clustering and classification.

## Project Overview

TCGA utilizes genomic, epigenomic, transcriptomic and proteomic data to improve diagnosis, treatment and prevention of cancer. This project aims to screen informative genetic samples from the TCGA cohort, explore the sample structure and cluster them, use decision trees for classification and evaluation.

## Repository Contents

- **`processed_data_sampled.xlsx`**  
  A dataset from the online TCGA database containing 232,000 cases with:
  - Genetic information (`gene_id`, `gene_name`, `gene_type`)
  - TPM (Transcripts per Million)
  - AJCC stage information (`ajcc_pathologic_m`, `ajcc_pathologic_n`, `ajcc_pathologic_t`, `ajcc_staging_system_edition`)
  - Patients case metadata (`age_at_index`, `days_to_death`, etc.)

- **`Script.R`**  
  An R script that performs the complete downstream analysis:
  - Screening for information-rich genes (high variant, high expression)
  - Downscaling to explore sample structure (PCA) and visualize by sample type
  - Visualize clustering of gene expression profiles (heatmaps)
  - Classification modeling of samples using decision trees and confusion matrix evaluation of model performance on test sets
 
## How to Use

1. **Clone the repository:**
   ```bash
   git clone https://github.com/<your-username>/<your-repo>.git
2. **Install required R packages and set random seed**, for example,
   ```R
   library(caret)
   library(ComplexHeatmap)
   library(circlize)
   library(data.table)
   library(ggrepel)
   library(ggplot2)
   # Add other packages as needed
   set.seed(281101) # any random number
3. **Run the analysis:**
   ```R
   source("Script.R")
All intermediate and final results will be saved in the working directory.

## References

- **Tools: [caret](https://cran.r-project.org/web/packages/caret/index.html), [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)**, etc (see `Script.R` requirements for details).

## License

This project is open for academic and educational purposes.
