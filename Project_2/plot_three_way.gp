load "plot_counts.gp"

set terminal pngcairo size 1100,700 enhanced font "Arial,11"

set output "three_way_comparison_distribution.png"

set title "Fragment Length Distribution Before and After Rescaling"

set xlabel "Fragment length (bp)"

set ylabel "Normalized frequency"

set xrange [0:700]

set yrange [0:*]

set key top right

set border linewidth 1

set tics out

set grid ytics

plot "three_way_comparison.tsv" using 1:2 \
    with lines linewidth 2 \
    linecolor rgb "#BA55D3" \
    title query_label, \
    "three_way_comparison.tsv" using 1:4 \
    with lines linewidth 2 \
    linecolor rgb "#36B99A" \
    title "Reference", \
    "three_way_comparison.tsv" using 1:3 \
    with linespoints linewidth 1.5 \
    pointtype 7 pointsize 0.35 \
    linecolor rgb "#35A7DF" \
    title rescaled_label
