# Differential Gene Expression Analysis for *Candida albicans*

This repository contains the code, data, and results for the differential gene expression (DGE) analysis of *Candida albicans*.

## Files and Directories

- `dge_analysis.R`: R script for running the DGE analysis.
- `data/`: Input files for the analysis.
  - `data_for_dge_csv.csv`: Raw count data.
  - `DGE_FILE_DATA.txt`: Metadata file.
- `results/`: Output files generated from the analysis.
  - `cts_integer.csv`: Normalized count data.
  - `top_diff_expressed_genes.csv`: Top differentially expressed genes.
- `visualizations/`: Visual plots generated during the analysis.

## Prerequisites

- R (version X.X.X or higher)
- Required R libraries: `DESeq2`, `ggplot2`, `pheatmap`, `writexl`, `reshape2`, `plotly`

## Running the Analysis

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/DGE_Candida_Albicans.git
   cd DGE_Candida_Albicans
   ```

2. Install required R packages:
   ```R
   install.packages(c("DESeq2", "ggplot2", "pheatmap", "writexl", "reshape2", "plotly"))
   ```

3. Run the R script:
   ```R
   source("dge_analysis.R")
   ```

## Contact

For questions or feedback, contact Dhruva Tirumalasetty at tdhruva970@gmail.com .
