
#!/usr/bin/env python3

import sys
from operator import itemgetter


def load_chromosome_ranks(filename):

    priority = {}

    with open(filename, "r") as file:

        for chromosome in file:

            name = chromosome.strip()

            if name and name not in priority:
                priority[name] = len(priority)

    return priority


def collect_bed_entries(rank_table):

    selected_records = []
    add_entry = selected_records.append
    find_rank = rank_table.get

    for raw_line in sys.stdin:

        bed_line = raw_line.strip()

        if not bed_line:
            continue

        chromosome, separator, remaining = bed_line.partition("\t")

        if not separator:
            continue

        position = find_rank(chromosome)

        if position is not None:
            add_entry((position, bed_line))

    return selected_records


def arrange_records(entries):

    entries.sort(key=itemgetter(0, 1))

    return entries


def print_bed_file(entries):

    output = sys.stdout.write

    for chromosome_position, bed_line in entries:
        output(bed_line)
        output("\n")


def main():

    chromosome_ranks = load_chromosome_ranks("selection.tsv")

    bed_entries = collect_bed_entries(chromosome_ranks)

    ordered_entries = arrange_records(bed_entries)

    print_bed_file(ordered_entries)


if __name__ == "__main__":
    main()

