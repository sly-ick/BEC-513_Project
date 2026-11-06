#!/usr/bin/env python3

import argparse
import gzip
import sys


def open_text(filename):

    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")

    return open(filename, "r")


def load_distribution(filename):

    frequencies = {}

    fragment_total = None

    with open_text(filename) as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):

                if "fragments" in line and "=" in line:

                    try:

                        fragment_total = int(
                            line.split("=", 1)[1].strip()
                        )

                    except ValueError:
                        pass

                continue

            columns = line.split()

            if len(columns) < 3:
                continue

            try:

                length = int(columns[0])

                query_value = float(columns[1])

                reference_value = float(columns[2])

            except ValueError:
                continue

            frequencies[length] = (
                query_value,
                reference_value
            )

    if not frequencies:
        raise ValueError(
            "No distribution values found in " + filename
        )

    return frequencies, fragment_total


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("original_table")

    parser.add_argument("rescaled_table")

    args = parser.parse_args()

    try:

        original, original_total = load_distribution(
            args.original_table
        )

        rescaled, rescaled_total = load_distribution(
            args.rescaled_table
        )

    except (OSError, ValueError) as error:

        sys.exit("Error: " + str(error))

    all_lengths = sorted(
        set(original) | set(rescaled)
    )

    print(
        "# original_fragments = {}".format(
            original_total
            if original_total is not None
            else "unknown"
        )
    )

    print(
        "# rescaled_fragments = {}".format(
            rescaled_total
            if rescaled_total is not None
            else "unknown"
        )
    )

    print(
        "length\toriginal_query\trescaled_query\treference"
    )

    for length in all_lengths:

        original_query, reference = original.get(
            length, (0.0, 0.0)
        )

        rescaled_query, _ = rescaled.get(
            length, (0.0, 0.0)
        )

        print(
            "{}\t{:.10f}\t{:.10f}\t{:.10f}".format(
                length,
                original_query,
                rescaled_query,
                reference
            )
        )

    original_label = "Query (n = {})".format(
        format(original_total, ",")
        if original_total is not None
        else "unknown"
    )

    rescaled_label = "Rep 0 (n = {})".format(
        format(rescaled_total, ",")
        if rescaled_total is not None
        else "unknown"
    )

    with open("plot_counts.gp", "w") as label_file:

        label_file.write(
            'query_label = "{}"\n'.format(original_label)
        )

        label_file.write(
            'rescaled_label = "{}"\n'.format(rescaled_label)
        )


if __name__ == "__main__":
    main()
