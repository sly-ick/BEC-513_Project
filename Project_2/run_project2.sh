#!/usr/bin/env bash

set -e

# Step 1: Calculate original query distribution

python3 query_hist.py \
    query.bed.gz \
    reference.hist \
    > query_distribution.tsv


# Step 2: Rescale query BED file

python3 rescale_query.py \
    query.bed.gz \
    query_distribution.tsv \
    --seed 123 \
    | gzip -c > rescaled_query.bed.gz


# Step 3: Calculate rescaled distribution

python3 query_hist.py \
    rescaled_query.bed.gz \
    reference.hist \
    > rescaled_distribution.tsv


# Step 4: Combine distributions

python3 combine_distributions.py \
    query_distribution.tsv \
    rescaled_distribution.tsv \
    > three_way_comparison.tsv


# Step 5: Generate final graph

gnuplot plot_three_way.gp


echo "Project 2 completed successfully!"
