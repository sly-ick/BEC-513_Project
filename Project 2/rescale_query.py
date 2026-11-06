#!/usr/bin/env python3

import argparse
import gzip
import random
import sys


def open_text(filename):

    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")

    return open(filename, "r")


def read_frequencies(filename):

    query_frequency = {}
    reference_frequency = {}

    with open_text(filename) as handle:

        for line in handle:

            line = line.strip()

            if not line or line.startswith("#"):
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

            query_frequency[length] = query_value

            reference_frequency[length] = reference_value

    if not query_frequency:
        raise ValueError("No frequency data found.")

    return query_frequency, reference_frequency


def make_keep_probabilities(
    query_frequency,
    reference_frequency
):

    ratios = {}

    for length, query_value in query_frequency.items():

        reference_value = reference_frequency.get(
            length, 0.0
        )

        if query_value > 0:

            ratios[length] = (
                reference_value / query_value
            )

        else:

            ratios[length] = 0.0

    largest_ratio = max(ratios.values())

    if largest_ratio <= 0:
        raise ValueError("No matching frequencies found.")

    keep_probability = {}

    for length, ratio in ratios.items():

        keep_probability[length] = (
            ratio / largest_ratio
        )

    return keep_probability


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("bed_file")

    parser.add_argument("frequency_file")

    parser.add_argument(
        "--seed",
        type=int,
        default=123
    )

    args = parser.parse_args()

    try:

        query_frequency, reference_frequency = (
            read_frequencies(args.frequency_file)
        )

        keep_probability = make_keep_probabilities(
            query_frequency,
            reference_frequency
        )

    except (OSError, ValueError) as error:

        sys.exit("Error: " + str(error))

    missing_lengths = [
        length
        for length, value in reference_frequency.items()
        if value > 0
        and query_frequency.get(length, 0.0) == 0
    ]

    if missing_lengths:

        print(
            "[rescale_query] Some reference lengths are absent "
            "from the query.",
            file=sys.stderr
        )

    rng = random.Random(args.seed)

    input_total = 0
    kept_total = 0
    invalid_total = 0

    try:

        with open_text(args.bed_file) as handle:

            for line in handle:

                stripped = line.strip()

                if not stripped or stripped.startswith("#"):
                    continue

                columns = stripped.split()

                if len(columns) < 3:
                    invalid_total += 1
                    continue

                try:

                    start = int(columns[1])

                    end = int(columns[2])

                except ValueError:

                    invalid_total += 1
                    continue

                length = end - start

                if length <= 0:
                    invalid_total += 1
                    continue

                input_total += 1

                probability = keep_probability.get(
                    length, 0.0
                )

                if rng.random() < probability:

                    sys.stdout.write(line)

                    kept_total += 1

    except OSError as error:

        sys.exit("Error: " + str(error))

    percentage = (
        100.0 * kept_total / max(input_total, 1)
    )

    print(
        "[rescale_query] Input: {} | Kept: {} ({:.2f}%) "
        "| Invalid: {} | Seed: {}".format(
            input_total,
            kept_total,
            percentage,
            invalid_total,
            args.seed
        ),
        file=sys.stderr
    )


if __name__ == "__main__":
    main()
