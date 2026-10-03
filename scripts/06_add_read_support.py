#!/usr/bin/env python3

import csv
from pathlib import Path


# ============================================================
# Input / output files
# ============================================================

EVIDENCE = Path(
    "results_summary/viral_contig_evidence.tsv"
)

SUPPORT = Path(
    "results_summary/viral_contig_read_support.tsv"
)

OUT = Path(
    "results_summary/viral_contig_evidence_with_reads.tsv"
)


# ============================================================
# Check input files
# ============================================================

for path in [EVIDENCE, SUPPORT]:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")


# ============================================================
# Load read-mapping support
# ============================================================

support = {}

with SUPPORT.open("r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    required_columns = {
        "contig",
        "breadth1x",
        "breadth5x",
        "breadth10x",
        "mean_depth",
        "length",
    }

    missing = required_columns - set(reader.fieldnames or [])

    if missing:
        raise ValueError(
            f"Missing required read-support columns: "
            f"{sorted(missing)}"
        )

    for row in reader:
        contig = row["contig"].strip()

        support[contig] = {
            "breadth1x": float(row["breadth1x"]),
            "breadth5x": float(row["breadth5x"]),
            "breadth10x": float(row["breadth10x"]),
            "mean_depth": float(row["mean_depth"]),
            "length": int(row["length"]),
        }


# ============================================================
# Merge master evidence table with read support
# ============================================================

output_rows = []

with EVIDENCE.open("r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    if not reader.fieldnames:
        raise ValueError(
            f"No header found in {EVIDENCE}"
        )

    original_fields = reader.fieldnames

    for row in reader:
        contig = row["contig"].strip()

        s = support.get(contig)

        if s is None:
            row["read_breadth1x"] = ""
            row["read_breadth5x"] = ""
            row["read_breadth10x"] = ""
            row["read_mean_depth"] = ""
            row["read_support_class"] = (
                "No_read_support_data"
            )

        else:
            b1 = s["breadth1x"]
            b5 = s["breadth5x"]
            b10 = s["breadth10x"]
            depth = s["mean_depth"]

            # Operational prioritization thresholds.
            # These labels describe read support only.
            # They are NOT species-confirmation criteria.
            if (
                b1 >= 95
                and b10 >= 80
                and depth >= 10
            ):
                support_class = "High_read_support"

            elif (
                b1 >= 80
                and b5 >= 50
                and depth >= 3
            ):
                support_class = "Moderate_read_support"

            else:
                support_class = "Low_read_support"

            row["read_breadth1x"] = f"{b1:.2f}"
            row["read_breadth5x"] = f"{b5:.2f}"
            row["read_breadth10x"] = f"{b10:.2f}"
            row["read_mean_depth"] = f"{depth:.2f}"
            row["read_support_class"] = support_class

        output_rows.append(row)


# ============================================================
# Write Linux-style TSV
# ============================================================

new_fields = original_fields + [
    "read_breadth1x",
    "read_breadth5x",
    "read_breadth10x",
    "read_mean_depth",
    "read_support_class",
]

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=new_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(output_rows)


# ============================================================
# Summary
# ============================================================

counts = {}

for row in output_rows:
    label = row["read_support_class"]
    counts[label] = counts.get(label, 0) + 1

print(f"Wrote merged evidence table:")
print(OUT)

print("\nRead-support classes:")
for label, count in sorted(
    counts.items(),
    key=lambda x: x[1],
    reverse=True,
):
    print(f"  {label}: {count}")
