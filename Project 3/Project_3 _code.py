#!/usr/bin/env python3
"""
Project 3: All-in-One V-Plot Generator & Plotter Suite
This script processes `mapped.bed` from stdin, automatically outputs the 2D frequency
matrix (`vplot_m1.txt`), and simultaneously renders a high-resolution V-plot heatmap image
(`vplot_heatmap.png`) using matplotlib/seaborn.
"""

import sys
import re
import os

def main():
    # Axis boundaries specification
    PM, PX = -500, 499  # Relative position range (X-axis)
    LM, LX = 50, 400    # Fragment length range (Y-axis)

    positions = list(range(PM, PX + 1))
    lengths = list(range(LM, LX + 1))

    pos_count = len(positions)
    matrix = [[0] * pos_count for _ in lengths]

    # Regex pattern to extract genomic coordinates from column 5 (e.g., "chr1:1000-2000")
    region_pattern = re.compile(r"^([^:|]+):(\d+)-(\d+)")

    # Caching identical consecutive motif sites to eliminate redundant regex matches
    last_site = ""
    site_center = 0

    total_processed = 0
    valid_counted = 0

    for line in sys.stdin:
        line_stripped = line.strip()
        if not line_stripped or line_stripped.startswith("#"):
            continue

        columns = line_stripped.split()
        total_processed += 1

        # Need at least 10 columns to access column 5 (index 4) and columns 9/10 (indices 8 and 9)
        if len(columns) < 10:
            continue

        # 1. Landmark motif center logic from column 5
        raw_site = columns[4]

        if raw_site != last_site:
            match = region_pattern.match(raw_site)
            if not match:
                continue

            last_site = raw_site
            ref_start = int(match.group(2))
            ref_end = int(match.group(3))

            # Calculate reference landmark midpoint center
            site_center = (ref_start + ref_end) // 2

        # 2. Fragment start and end from columns 9 and 10 (indices 8 and 9)
        try:
            start = int(columns[8])
            end = int(columns[9])
        except ValueError:
            continue

        # Calculate fragment length and relative midpoint position
        length = end - start
        position = (start + end) // 2 - site_center

        # 3. Populate matrix if within defined boundaries
        if PM <= position <= PX and LM <= length <= LX:
            matrix[length - LM][position - PM] += 1
            valid_counted += 1

    # Log metrics to stderr
    sys.stderr.write(f"[P3t V-Plot] Total lines read: {total_processed} | Valid fragments mapped: {valid_counted}\n")

    # Output file paths
    matrix_filename = "vplot_m1.txt"
    plot_filename = "vplot_heatmap.png"

    # Write frequency matrix to file and stdout simultaneously
    with open(matrix_filename, "w") as f_out:
        header = "length\t" + "\t".join(map(str, positions)) + "\n"
        f_out.write(header)
        sys.stdout.write(header)

        for i, length in enumerate(lengths):
            row_data = "\t".join(str(c) for c in matrix[i])
            row_str = f"{length}\t{row_data}\n"
            f_out.write(row_str)
            sys.stdout.write(row_str)

    sys.stderr.write(f"[P3t V-Plot] Matrix successfully saved to '{matrix_filename}'\n")

    # Automatically generate the V-plot heatmap image simultaneously
    try:
        import matplotlib
        matplotlib.use('Agg') # Non-interactive backend for headless Linux/servers
        import matplotlib.pyplot as plt
        import numpy as np

        sys.stderr.write("[P3t V-Plot] Generating heatmap image...\n")

        # Convert matrix to numpy array for plotting
        data_arr = np.array(matrix)

        plt.figure(figsize=(10, 7), dpi=300)

        # Plot using imshow for fast, clean heatmap rendering
        plt.imshow(
            data_arr,
            extent=[PM, PX, LX, LM],
            aspect='auto',
            cmap='viridis',
            origin='lower'
        )

        plt.colorbar(label="Fragment Count")
        plt.title("V-Plot: Fragment Length vs. Distance to Landmark Center", fontsize=14, fontweight='bold')
        plt.xlabel("Distance from Landmark Center (bp)", fontsize=12)
        plt.ylabel("Fragment Length (bp)", fontsize=12)

        plt.tight_layout()
        plt.savefig(plot_filename, dpi=300)
        plt.close()

        sys.stderr.write(f"[P3t V-Plot] Plot successfully saved to '{plot_filename}'\n")

    except ImportError:
        sys.stderr.write("[P3t Warning] 'matplotlib' or 'numpy' not found. Matrix generated, but image plotting was skipped.\n")
        sys.stderr.write("[P3t Warning] Install them via: pip install matplotlib numpy\n")

if __name__ == "__main__":
    main()
