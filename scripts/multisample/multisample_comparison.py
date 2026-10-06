#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from itertools import combinations

parser = argparse.ArgumentParser(
    description="Compare vOTU prevalence and community composition."
)
parser.add_argument("--input", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--threshold", type=float, required=True)
args = parser.parse_args()

INPUT = Path(args.input)
OUTDIR = Path(args.outdir)
THRESHOLD = args.threshold

OUTDIR.mkdir(parents=True, exist_ok=True)

PRESENCE_TAG = int(round(THRESHOLD * 100))

# -------------------------------------------------------
# Read breadth matrix
# -------------------------------------------------------

with INPUT.open() as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    samples = header[1:]

    breadth = {}

    for row in reader:
        votu = row[0]
        breadth[votu] = {
            sample: float(value)
            for sample, value in zip(samples, row[1:])
        }

votus = sorted(breadth)

# -------------------------------------------------------
# High-confidence presence / absence
# -------------------------------------------------------

presence = {}

for votu in votus:
    presence[votu] = {
        sample: int(breadth[votu][sample] >= THRESHOLD)
        for sample in samples
    }

presence_file = OUTDIR / f"votu_presence_breadth{PRESENCE_TAG}.tsv"

with presence_file.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["votu_id"] + samples)

    for votu in votus:
        writer.writerow(
            [votu] +
            [presence[votu][s] for s in samples]
        )

# -------------------------------------------------------
# vOTU prevalence
# -------------------------------------------------------

prevalence_file = OUTDIR / "votu_prevalence.tsv"

prevalence_counts = {}

with prevalence_file.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "votu_id",
        "n_samples_present",
        "prevalence_fraction",
        "samples_present",
    ])

    for votu in votus:

        positive = [
            s for s in samples
            if presence[votu][s] == 1
        ]

        n = len(positive)

        prevalence_counts[n] = (
            prevalence_counts.get(n, 0) + 1
        )

        writer.writerow([
            votu,
            n,
            n / len(samples),
            ",".join(positive),
        ])

# -------------------------------------------------------
# Prevalence distribution
# -------------------------------------------------------

dist_file = OUTDIR / "prevalence_distribution.tsv"

with dist_file.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "n_samples_present",
        "n_votus"
    ])

    for n in range(len(samples) + 1):
        writer.writerow([
            n,
            prevalence_counts.get(n, 0)
        ])

# -------------------------------------------------------
# Sample richness
# -------------------------------------------------------

richness_file = OUTDIR / "sample_votu_richness.tsv"

with richness_file.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "sample",
        "high_confidence_votus"
    ])

    for sample in samples:

        n = sum(
            presence[v][sample]
            for v in votus
        )

        writer.writerow([
            sample,
            n
        ])

# -------------------------------------------------------
# Pairwise Jaccard similarity
# -------------------------------------------------------

jaccard_file = OUTDIR / "pairwise_jaccard.tsv"

with jaccard_file.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "sample_1",
        "sample_2",
        "shared_votus",
        "union_votus",
        "jaccard"
    ])

    for s1, s2 in combinations(samples, 2):

        set1 = {
            v for v in votus
            if presence[v][s1]
        }

        set2 = {
            v for v in votus
            if presence[v][s2]
        }

        shared = len(set1 & set2)
        union = len(set1 | set2)

        j = shared / union if union else 0.0

        writer.writerow([
            s1,
            s2,
            shared,
            union,
            j
        ])

# -------------------------------------------------------
# Summary categories
# -------------------------------------------------------

sample_specific = []
shared = []
core = []

for votu in votus:

    n = sum(
        presence[votu][s]
        for s in samples
    )

    if n == 1:
        sample_specific.append(votu)

    if n >= 2:
        shared.append(votu)

    if n == len(samples):
        core.append(votu)

summary_file = OUTDIR / "multisample_summary.tsv"

with summary_file.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["metric", "value"])

    writer.writerow([
        "total_catalog_votus",
        len(votus)
    ])

    writer.writerow([
        "breadth_threshold",
        THRESHOLD
    ])

    writer.writerow([
        "sample_specific_votus",
        len(sample_specific)
    ])

    writer.writerow([
        "shared_ge2_samples",
        len(shared)
    ])

    writer.writerow([
        "core_all_7_samples",
        len(core)
    ])

# -------------------------------------------------------
# Lists
# -------------------------------------------------------

for name, values in [
    ("sample_specific_votus.txt", sample_specific),
    ("shared_votus.txt", shared),
    ("core_votus.txt", core),
]:

    with (OUTDIR / name).open("w") as f:
        for v in values:
            f.write(v + "\n")

print("===================================")
print("Multisample comparison complete")
print("===================================")
print(f"Catalog vOTUs       : {len(votus)}")
print(f"Presence threshold  : breadth >= {THRESHOLD}")
print(f"Sample-specific     : {len(sample_specific)}")
print(f"Shared >=2 samples  : {len(shared)}")
print(f"Core all 7 samples  : {len(core)}")
