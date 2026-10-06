#!/usr/bin/env python3

import csv
from itertools import combinations

FILE = (
    "results/multisample/abundance/"
    "votu_covered_fraction_matrix.tsv"
)

THRESHOLDS = [0.50, 0.75]

with open(FILE) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    samples = header[1:]

    data = {}

    for row in reader:
        votu = row[0]

        data[votu] = {
            sample: float(value)
            for sample, value in zip(samples, row[1:])
        }

votus = sorted(data)

print("=== GLOBAL SUMMARY ===")
print(
    "threshold\t"
    "sample_specific\t"
    "shared_ge2\t"
    "present_ge5\t"
    "present_ge6\t"
    "core_7"
)

for threshold in THRESHOLDS:

    prevalence = {
        votu: sum(
            data[votu][sample] >= threshold
            for sample in samples
        )
        for votu in votus
    }

    print(
        f"{threshold}\t"
        f"{sum(n == 1 for n in prevalence.values())}\t"
        f"{sum(n >= 2 for n in prevalence.values())}\t"
        f"{sum(n >= 5 for n in prevalence.values())}\t"
        f"{sum(n >= 6 for n in prevalence.values())}\t"
        f"{sum(n == 7 for n in prevalence.values())}"
    )

print()

for threshold in THRESHOLDS:

    print(f"=== SAMPLE RICHNESS threshold={threshold} ===")

    for sample in samples:

        n = sum(
            data[votu][sample] >= threshold
            for votu in votus
        )

        print(f"{sample}\t{n}")

    print()

    print(f"=== TOP JACCARD threshold={threshold} ===")

    results = []

    for s1, s2 in combinations(samples, 2):

        A = {
            votu for votu in votus
            if data[votu][s1] >= threshold
        }

        B = {
            votu for votu in votus
            if data[votu][s2] >= threshold
        }

        shared = len(A & B)
        union = len(A | B)

        jaccard = shared / union if union else 0.0

        results.append(
            (jaccard, s1, s2, shared, union)
        )

    for jaccard, s1, s2, shared, union in sorted(
        results,
        reverse=True
    )[:10]:

        print(
            f"{s1}\t{s2}\t"
            f"{shared}\t{union}\t"
            f"{jaccard:.4f}"
        )

    print()
