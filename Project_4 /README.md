# Project 4: TCGA-LUAD NKX2-1 Expression Analysis

## Overview
**Project 4** is a bioinformatics pipeline focused on analyzing **NKX2-1** gene expression in Lung Adenocarcinoma (TCGA-LUAD) datasets. It includes a custom R script to automate the extraction of NKX2-1 gene expression values (TPM) from RNA-Seq count files, performs statistical testing (Welch's two-sample t-test), and generates visualization plots.

## Project Structure
- **`run_analysis.R`**: The main R script that performs data extraction, statistical testing, and visualization.
- **`README.md`**: Project documentation and instructions.
- **Input Data**: GDC RNA-seq gene expression count files and sample metadata (gdc_sample_sheet.tsv).
- **Output**: Expression comparison results, statistical analysis, and visualization plots.

---

## Prerequisites (Linux/macOS)
- **R (v4.0+)**: Must be installed and available in PATH (`Rscript` binary accessible).
- **Required R packages**: `readr`, `dplyr`, `ggplot2`, `tibble`.
- **Data**: Sample sheet (`gdc_sample_sheet.tsv`) and compressed expression data archive (`tcga_data.tar.gz` or `.zip`).

To verify R is available:
```bash
Rscript --version
```

To check R package installation:
```bash
R -e "library(readr); library(dplyr); library(ggplot2); library(tibble)"
```

---

## Running the Pipeline on Linux

### Step 1: Navigate to Project Directory
```bash
cd /path/to/Project\ 4
```

### Step 2: Ensure Archive Extraction (If applicable)
If your raw GDC data is compressed, extract it first:
```bash
tar -xzf tcga_data.tar.gz
```

### Step 3: Execute the R Script
Run the analysis script via `Rscript`:
```bash
Rscript run_analysis.R
```

---

## Outputs Generated

Upon execution, the script creates a `results/` folder containing:

- **`NKX2-1_PT_vs_STN_boxplot.png`**: The expression comparison plot.
- **`NKX2-1_PT_vs_STN_boxplot.pdf`**: Vector format of the plot.
- **`NKX2-1_expression_PT_STN_sample_values.csv`**: Detailed sample TPM values and log2-transformed values.
- **`NKX2-1_expression_summary_PT_STN.csv`**: Summary statistics per tissue condition.

---

## Additional OS Notes

<details>
<summary>Running on Windows (PowerShell)</summary>


</details>

---

## Troubleshooting

### Error: "No rna_seq.augmented_star_gene_counts.tsv files were found"
The script expects these files to be in the project directory, either already extracted or newly extracted from `tcga_data.tar.gz`. If the archive hasn't been extracted:

1. Verify `tcga_data.tar.gz` exists.
2. Manually extract it:
   ```bash
   tar -xzf tcga_data.tar.gz
   ```
3. Re-run the script.

### Error: "gene_name column missing"
Some RNA-Seq count files have an extra header comment line (`# gene-model: GENCODE v36`). The script includes `comment = "#"` in the file reading step to skip it, but if the column is still missing:

1. Verify the TSV file has `gene_name` and `tpm_unstranded` columns.
2. If the column positions differ, the script may fall back to column 7 for TPM, but the `gene_name` column should still be present.

---

## Linux Notes

- The script uses standard Unix tools for archive extraction (`untar`, `unzip`) and file globbing.
- All file paths are handled in a platform-independent way via R functions.
- The script is fully self-contained; no external package installation is required if dependencies are already present.

This script is now designed for **robust, automated execution** on Linux systems with a comprehensive data processing pipeline, complete error handling, and professional output visualization.