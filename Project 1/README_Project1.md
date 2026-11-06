# Project 1 – Quick BED Sort

## Files required

Keep these files together in the `Project1` directory:

```text
Project1/
├── bed_reorder.py
├── selection.tsv
├── shuf.a.bed.gz
├── shuf.b.bed.gz
└── sorted_bed_file_per_sample/
```

- `bed_reorder.py` – Python program.
- `selection.tsv` – chromosome priority/order list.
- `shuf.a.bed.gz` and `shuf.b.bed.gz` – compressed BED input files.
- `sorted_bed_file_per_sample/` – directory for the output file.

## Requirements

Run the project in Ubuntu/WSL with Python 3 and `gzip` available.

Check the installations:

```bash
python3 --version
gzip --version
```

## Execute the program

Open Ubuntu/WSL and move to the project directory:

```bash
cd ~/Project1
```

Create the output directory (only needed the first time):

```bash
mkdir -p sorted_bed_file_per_sample
```

Run the pipeline:

```bash
gzip -cd shuf.a.bed.gz shuf.b.bed.gz | python3 bed_reorder.py > sorted_bed_file_per_sample/X.standard.bed
```

Run this command from the `Project1` directory because the Python script reads `selection.tsv` from the current directory.

## Check the output

List the output file:

```bash
ls -lh sorted_bed_file_per_sample/
```

Display the first 10 records:

```bash
head -n 10 sorted_bed_file_per_sample/X.standard.bed
```

Count the output records:

```bash
wc -l sorted_bed_file_per_sample/X.standard.bed
```

Display the chromosome order present in the output:

```bash
cut -f1 sorted_bed_file_per_sample/X.standard.bed | uniq
```

## Important

The command uses `>` to write the result. Running it again replaces the existing `X.standard.bed` file.

The output filename is `X.standard.bed`. Change `X` in the command only when a different sample name is required.
