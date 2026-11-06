#!/usr/bin/env python3

import argparse
import gzip
import sys
from collections import Counter


def open_text(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, "r")


def read_reference(filename):
    reference = Counter()

    with open_text(filename) as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            columns = line.split()

            if len(columns) < 2:
                continue

            try:
                length = int(float(columns[0]))
                probability = float(columns[1])
            except ValueError:
                continue

            if length >= 0 and probability >= 0:
                reference[length] += probability

    total = sum(reference.values())

    if total <= 0:
        raise ValueError("Invalid reference distribution.")

    for length in reference:
        reference[length] /= total

    return reference


def count_bed_lengths(filename):
    counts = Counter()
    total = 0
    skipped = 0

    with open_text(filename) as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            columns = line.split()

            if len(columns) < 3:
                skipped += 1
                continue

            try:
                start = int(columns[1])
                end = int(columns[2])
            except ValueError:
                skipped += 1
                continue

            length = end - start

            if length <= 0:
                skipped += 1
                continue

            counts[length] += 1
            total += 1

    if total == 0:
        raise ValueError("No valid BED fragments found.")

    return counts, total, skipped


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("bed_file")
    parser.add_argument("reference_file")

    args = parser.parse_args()

    try:
        reference = read_reference(args.reference_file)

        counts, total, skipped = count_bed_lengths(
            args.bed_file
        )

    except (OSError, ValueError) as error:
        sys.exit("Error: " + str(error))

    query = {
        length: count / total
        for length, count in counts.items()
    }

    all_lengths = sorted(
        set(query) | set(reference)
    )

    print("# fragments = {}".format(total))
    print("# skipped BED lines = {}".format(skipped))

    print(
        "length\tquery_frequency\treference_frequency"
    )

    for length in all_lengths:

        query_frequency = query.get(length, 0.0)

        reference_frequency = reference.get(
            length, 0.0
        )

        print(
            "{}\t{:.10f}\t{:.10f}".format(
                length,
                query_frequency,
                reference_frequency
            )
        )

    print(
        "[query_hist] Valid fragments: {} | Skipped lines: {}".format(
            total, skipped
        ),
        file=sys.stderr
    )


if __name__ == "__main__":
    main()
    
