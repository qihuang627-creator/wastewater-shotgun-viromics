#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict, OrderedDict
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Build sample-by-vOTU abundance and coverage matrices."
)
parser.add_argument("--input", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--samples", nargs="+", required=True)
args = parser.parse_args()

input_file = Path(args.input)
outdir = Path(args.outdir)
samples_order = args.samples

outdir.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------
# Read CoverM sparse output
# -------------------------------------------------------

rows = []

with input_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    print("CoverM columns:")
    for x in reader.fieldnames:
        print("  ", x)

    for row in reader:

        sample_raw = row["Sample"]

        m = re.search(r"(ERR\d+)", sample_raw)

        if not m:
            raise RuntimeError(
                f"Cannot infer ERR accession from Sample field: "
                f"{sample_raw}"
            )

        sample = m.group(1)

        contig_col = (
            "Contig"
            if "Contig" in row
            else "Genome"
            if "Genome" in row
            else None
        )

        if contig_col is None:
            raise RuntimeError(
                "Could not find Contig column in CoverM output"
            )

        votu = row[contig_col]

        rows.append((sample, votu, row))

# -------------------------------------------------------
# Identify metric columns dynamically
# -------------------------------------------------------

fieldnames = reader.fieldnames

metric_map = {}

for col in fieldnames:

    low = col.lower()

    if low in ("sample", "contig", "genome"):
        continue

    if "covered fraction" in low:
        metric_map["covered_fraction"] = col
    elif low.endswith("count") or low == "count":
        metric_map["count"] = col
    elif low.endswith("mean") or low == "mean":
        metric_map["mean"] = col
    elif "rpkm" in low:
        metric_map["rpkm"] = col
    elif "tpm" in low:
        metric_map["tpm"] = col

print("Detected metrics:")
for k, v in metric_map.items():
    print(f"  {k}: {v}")

required = {
    "count",
    "mean",
    "covered_fraction",
    "rpkm",
    "tpm",
}

missing = required - set(metric_map)

if missing:
    raise RuntimeError(
        f"Missing CoverM metrics: {sorted(missing)}"
    )

# -------------------------------------------------------
# Store values
# -------------------------------------------------------

data = {
    metric: defaultdict(dict)
    for metric in required
}

all_votus = set()

for sample, votu, row in rows:

    all_votus.add(votu)

    for metric, col in metric_map.items():

        value = row[col]

        try:
            value = float(value)
        except ValueError:
            value = 0.0

        data[metric][votu][sample] = value

# vOTU IDs are zero padded, lexical sorting is fine
all_votus = sorted(all_votus)

# -------------------------------------------------------
# Write one matrix per metric
# -------------------------------------------------------

for metric in [
    "count",
    "mean",
    "covered_fraction",
    "rpkm",
    "tpm",
]:

    outfile = outdir / f"votu_{metric}_matrix.tsv"

    with outfile.open("w", newline="") as f:

        writer = csv.writer(f, delimiter="\t")

        writer.writerow(
            ["votu_id"] + samples_order
        )

        for votu in all_votus:

            values = [
                data[metric][votu].get(sample, 0.0)
                for sample in samples_order
            ]

            writer.writerow(
                [votu] + values
            )

# -------------------------------------------------------
# Preliminary detection summary only
#
# Do NOT hard-code presence yet.
# Report how many vOTUs have any mapped reads and several
# breadth thresholds so we can choose sensibly afterward.
# -------------------------------------------------------

summary_file = outdir / "sample_mapping_summary.tsv"

with summary_file.open("w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "sample",
        "votus_count_gt0",
        "votus_breadth_ge_0.10",
        "votus_breadth_ge_0.25",
        "votus_breadth_ge_0.50",
        "votus_breadth_ge_0.75",
    ])

    for sample in samples_order:

        count_gt0 = 0
        thresholds = {
            0.10: 0,
            0.25: 0,
            0.50: 0,
            0.75: 0,
        }

        for votu in all_votus:

            count = data["count"][votu].get(sample, 0.0)
            breadth = data["covered_fraction"][votu].get(
                sample, 0.0
            )

            if count > 0:
                count_gt0 += 1

            for threshold in thresholds:
                if breadth >= threshold:
                    thresholds[threshold] += 1

        writer.writerow([
            sample,
            count_gt0,
            thresholds[0.10],
            thresholds[0.25],
            thresholds[0.50],
            thresholds[0.75],
        ])

print()
print(f"vOTUs: {len(all_votus)}")
print("Matrices written to:")
print(outdir)
