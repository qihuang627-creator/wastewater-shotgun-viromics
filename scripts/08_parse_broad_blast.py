#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path


BLAST = Path(
    "results/broad_blast/"
    "ERR14789190_high_priority_core_nt_megablast.csv"
)

MASTER = Path(
    "results_summary/"
    "viral_contig_evidence_with_reads.tsv"
)

QUERY_LIST = Path(
    "results_summary/"
    "high_priority_panel_external_contigs.txt"
)

OUT = Path(
    "results_summary/"
    "high_priority_broad_nt_summary.tsv"
)


# ------------------------------------------------------------
# Load query lengths from master evidence table
# ------------------------------------------------------------

query_lengths = {}

with MASTER.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        query_lengths[row["contig"]] = int(row["length"])


# ------------------------------------------------------------
# Load the 17 selected query IDs
# ------------------------------------------------------------

queries = []

with QUERY_LIST.open() as f:
    for line in f:
        q = line.strip()
        if q:
            queries.append(q)


# ------------------------------------------------------------
# Parse BLAST CSV
#
# Columns:
# qseqid
# sseqid
# pident
# length
# mismatch
# gapopen
# qstart
# qend
# sstart
# send
# evalue
# bitscore
# ------------------------------------------------------------

hits = defaultdict(lambda: defaultdict(list))

with BLAST.open() as f:
    reader = csv.reader(f)

    for row in reader:
        if len(row) < 12:
            continue

        qseqid = row[0]
        sseqid = row[1]

        hits[qseqid][sseqid].append(
            {
                "pident": float(row[2]),
                "aln_len": int(row[3]),
                "qstart": int(row[6]),
                "qend": int(row[7]),
                "evalue": float(row[10]),
                "bitscore": float(row[11]),
            }
        )


# ------------------------------------------------------------
# Calculate union coverage across HSPs
# ------------------------------------------------------------

def union_length(intervals):
    if not intervals:
        return 0

    intervals = sorted(intervals)

    merged = []

    for start, end in intervals:
        if start > end:
            start, end = end, start

        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(
                merged[-1][1],
                end
            )

    return sum(
        end - start + 1
        for start, end in merged
    )


# ------------------------------------------------------------
# Summarize each query
# ------------------------------------------------------------

results = []

for query in queries:

    qlen = query_lengths.get(query)

    if query not in hits:

        results.append(
            {
                "contig": query,
                "length": qlen,
                "best_subject": "",
                "n_hsps": 0,
                "query_coverage_pct": 0,
                "weighted_identity_pct": "",
                "total_bitscore": "",
                "best_evalue": "",
                "nt_match_class": "No_megablast_hit",
            }
        )

        continue

    subject_summaries = []

    for subject, hsps in hits[query].items():

        intervals = [
            (h["qstart"], h["qend"])
            for h in hsps
        ]

        covered_bp = union_length(intervals)

        qcov = (
            100 * covered_bp / qlen
            if qlen else 0
        )

        total_aln = sum(
            h["aln_len"]
            for h in hsps
        )

        weighted_identity = (
            sum(
                h["pident"] * h["aln_len"]
                for h in hsps
            )
            / total_aln
        )

        total_bitscore = sum(
            h["bitscore"]
            for h in hsps
        )

        best_evalue = min(
            h["evalue"]
            for h in hsps
        )

        subject_summaries.append(
            {
                "subject": subject,
                "n_hsps": len(hsps),
                "qcov": qcov,
                "identity": weighted_identity,
                "bitscore": total_bitscore,
                "evalue": best_evalue,
            }
        )

    # Choose subject with greatest total bitscore
    best = max(
        subject_summaries,
        key=lambda x: x["bitscore"]
    )

    qcov = best["qcov"]
    identity = best["identity"]

    # Operational categories only
    if qcov >= 90 and identity >= 90:
        match_class = "Strong_near_known_nt_match"

    elif qcov >= 50 and identity >= 75:
        match_class = "Moderate_divergent_nt_match"

    elif qcov >= 10:
        match_class = "Partial_or_distant_nt_match"

    else:
        match_class = "Weak_nt_match"

    results.append(
        {
            "contig": query,
            "length": qlen,
            "best_subject": best["subject"],
            "n_hsps": best["n_hsps"],
            "query_coverage_pct": f"{qcov:.2f}",
            "weighted_identity_pct": f"{identity:.2f}",
            "total_bitscore": f"{best['bitscore']:.1f}",
            "best_evalue": best["evalue"],
            "nt_match_class": match_class,
        }
    )


# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

fields = [
    "contig",
    "length",
    "best_subject",
    "n_hsps",
    "query_coverage_pct",
    "weighted_identity_pct",
    "total_bitscore",
    "best_evalue",
    "nt_match_class",
]

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(results)


print(f"Wrote {len(results)} contigs to:")
print(OUT)

print("\nClassification counts:")

counts = defaultdict(int)

for row in results:
    counts[row["nt_match_class"]] += 1

for key, value in sorted(
    counts.items(),
    key=lambda x: x[1],
    reverse=True,
):
    print(f"  {key}: {value}")
