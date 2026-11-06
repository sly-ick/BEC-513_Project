# Project 2 — Rescale to Reference

## 1. Input files

Place the following input files in the project directory:

| File | Purpose |
|---|---|
| `query.bed.gz` | Query BED fragments to be rescaled. |
| `reference.hist` | Reference fragment-length distribution. |

The following files may also be present from Project 1, but they are **not used** by this Project 2 pipeline:

- `selection.tsv`
- `shuf.a.bed.gz`
- `shuf.b.bed.gz`

## 2. Code files

Keep these five files in the same directory as the input files:

- `query_hist.py` — calculates the normalized fragment-length distribution.
- `rescale_query.py` — resamples query BED records using length-specific retention probabilities.
- `combine_distributions.py` — combines the original, rescaled, and reference distributions.
- `plot_three_way.gp` — Gnuplot instructions for the comparison plot.
- `run_project2.sh` — runs the complete pipeline.

## 3. Requirements

Install Python 3 and Gnuplot. On Ubuntu/WSL, install Gnuplot with:

```bash
sudo apt update
sudo apt install gnuplot
```

## 4. Execution steps

1. Open the Project 2 directory in Ubuntu/WSL.
2. Make the pipeline script executable:

   ```bash
   chmod +x run_project2.sh
   ```

3. Run the complete pipeline:

   ```bash
   ./run_project2.sh
   ```

The script uses `query.bed.gz` and `reference.hist` and runs all analysis and plotting steps automatically. The random seed is set to `123` in `run_project2.sh` to make the sampling reproducible.

## 5. Output files

After successful execution, the following files are generated:

| File | Description |
|---|---|
| `query_distribution.tsv` | Original query fragment-length distribution. |
| `rescaled_query.bed.gz` | BED records retained after rejection sampling. |
| `rescaled_distribution.tsv` | Fragment-length distribution of the rescaled query. |
| `three_way_comparison.tsv` | Combined original-query, rescaled-query, and reference distributions. |
| `plot_counts.gp` | Gnuplot labels containing fragment counts. |
| `three_way_comparison_distribution.png` | Final plot comparing the three distributions. |

## 6. Project directory structure

```text
Project_2/
├── query.bed.gz
├── reference.hist
├── selection.tsv                 # Project 1 file; not used here
├── shuf.a.bed.gz                 # Project 1 file; not used here
├── shuf.b.bed.gz                 # Project 1 file; not used here
├── query_hist.py
├── rescale_query.py
├── combine_distributions.py
├── plot_three_way.gp
├── run_project2.sh
├── query_distribution.tsv        # Generated
├── rescaled_query.bed.gz         # Generated
├── rescaled_distribution.tsv     # Generated
├── three_way_comparison.tsv      # Generated
├── plot_counts.gp                # Generated
└── three_way_comparison_distribution.png  # Generated
```
