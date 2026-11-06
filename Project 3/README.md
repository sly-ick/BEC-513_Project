# Project P3t: V-Plot Generation

## Overview
This suite generates a 2D V-plot frequency matrix (`vplot_m1.txt`) from mapped genomic BED files (`mapped.bed`), calculating fragment length versus relative distance to reference landmark centers.

**Key Feature**: The script now **automatically generates both** the frequency matrix file **and** the V-plot heatmap image (`vplot_heatmap.png`) from a **single execution**.

---

## Linux Execution Guide

To run this V-plot script on a Linux environment to generate the V-plot suite, follow these steps:

### 1. Ensure Execute Permissions (Optional)
Make sure Python 3 is installed (`python3 --version`). You can make the script executable:
```bash
chmod +x "Project_3 _code.py"
```

### 2. Run the Command
Pipe `mapped.bed` into the Python script to produce the frequency matrix and image simultaneously:
```bash
python3 "Project_3 _code.py" < mapped.bed
```

*(If your input BED file is compressed)*:
```bash
gunzip -c mapped.bed.gz | python3 "Project_3 _code.py"
```

### 3. What This Does
- **Reads `mapped.bed`**: Processes BED format to extract landmark centers (Column 5) and fragment coordinates (Columns 9 & 10)
- **Generates `vplot_m1.txt`**: 2D frequency matrix (fragment length vs. position) and saves it to file and stdout
- **Creates `vplot_heatmap.png`**: Automatically renders a professional V-plot heatmap visualization using matplotlib
- **Logs Progress**: All processing metrics and completion messages to stderr (terminal)

---

## Output Files

After execution, you will have:

1. **`vplot_m1.txt`**: Tab-separated frequency matrix with rows representing fragment lengths (50-400 bp) and columns representing relative positions (-500 to +499 bp)

2. **`vplot_heatmap.png`**: High-resolution (300 DPI) V-plot heatmap visualization of the frequency matrix

3. **Console Output**: Terminal will show processing statistics and completion status

---

## Dependencies

The script requires:
- **Python 3.x**
- **matplotlib** (for generating the V-plot image)
- **numpy** (for matrix operations, usually included with matplotlib)

**Install dependencies** (if not already installed):
```bash
pip install matplotlib numpy
```

*Note: If matplotlib/numpy are unavailable, the matrix file will still be generated. The image generation will be skipped with a clear warning.*

---

## Quick Start

```bash
# Navigate to your project directory
cd /path/to/P3t

# Generate V-plot suite from BED file
python3 "Project_3 _code.py" < mapped.bed

# Or from compressed file
python3 "Project_3 _code.py" < <(gunzip -c mapped.bed.gz)
```

---

## Notes

- The script uses efficient regex matching to extract genomic landmark centers from Column 5
- Implements caching of consecutive identical landmark sites to reduce processing overhead
- Progress monitoring is provided via stderr for large datasets
- Automatic generation of both data file and visualization ensures a complete workflow

**This single script handles the complete V-plot generation pipeline — from raw BED data to frequency matrix and visual heatmap representation.**