#!/usr/bin/env python3

import csv
from pathlib import Path


# ============================================================
# Input / output files
# ============================================================

GENOMAD = Path(
    "results/genomad/ERR14789190/"
    "ERR14789190_contigs_1kb_summary/"
    "ERR14789190_contigs_1kb_virus_summary.tsv"
)

BLAST = Path(
    "results/assembly_validation/"
    "ERR14789190_contigs_vs_candidates.tsv"
)

ACCESSIONS = Path(
    "refs/candidate_viruses/accessions.tsv"
)

OUT = Path(
    "results_summary/viral_contig_evidence.tsv"
)


# ============================================================
# Check input files
# ============================================================

for path in [GENOMAD, BLAST, ACCESSIONS]:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")


# ============================================================
# Load accession -> candidate name mapping
#
# accessions.tsv format:
# taxid    name    accession
# ============================================================

accession_to_name = {}

with ACCESSIONS.open("r", newline="") as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        if len(row) < 3:
            continue

        taxid, name, accession = row[:3]

        accession_to_name[accession.strip()] = name.strip()


# ============================================================
# Load BLAST results
#
# BLAST columns:
# 1  qseqid
# 2  qlen
# 3  sseqid
# 4  slen
# 5  pident
# 6  alignment length
# 7  qcovs
# 8  evalue
# 9  bitscore
#
# Retain the highest-bitscore hit for each contig.
# ============================================================

best_blast = {}

with BLAST.open("r", newline="") as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        if len(row) < 9:
            continue

        qseqid = row[0].strip()
        blast_ref = row[2].strip()

        hit = {
            "blast_ref": blast_ref,
            "blast_name": accession_to_name.get(blast_ref, ""),
            "blast_identity": float(row[4]),
            "blast_alignment_length": int(row[5]),
            "blast_query_coverage": float(row[6]),
            "blast_evalue": row[7],
            "blast_bitscore": float(row[8]),
        }

        if (
            qseqid not in best_blast
            or hit["blast_bitscore"]
            > best_blast[qseqid]["blast_bitscore"]
        ):
            best_blast[qseqid] = hit


# ============================================================
# Merge geNomad results with candidate-panel BLAST evidence
# ============================================================

output_rows = []

with GENOMAD.open("r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    required_columns = {
        "seq_name",
        "length",
        "topology",
        "n_genes",
        "virus_score",
        "n_hallmarks",
        "marker_enrichment",
        "taxonomy",
    }

    missing = required_columns - set(reader.fieldnames or [])

    if missing:
        raise ValueError(
            f"Missing required geNomad columns: {sorted(missing)}"
        )

    for row in reader:
        contig = row["seq_name"].strip()

        hit = best_blast.get(contig)

        if hit is None:
            blast_ref = ""
            blast_name = ""
            blast_identity = ""
            blast_alignment_length = ""
            blast_query_coverage = ""
            blast_evalue = ""
            blast_bitscore = ""

            evidence_class = "No_match_to_candidate_panel"

        else:
            blast_ref = hit["blast_ref"]
            blast_name = hit["blast_name"]
            blast_identity = hit["blast_identity"]
            blast_alignment_length = hit["blast_alignment_length"]
            blast_query_coverage = hit["blast_query_coverage"]
            blast_evalue = hit["blast_evalue"]
            blast_bitscore = hit["blast_bitscore"]

            # Operational prioritization thresholds.
            # These are not species-confirmation criteria.
            if (
                blast_identity >= 95
                and blast_query_coverage >= 90
            ):
                evidence_class = "Strong_panel_match"

            elif (
                blast_identity >= 90
                and blast_query_coverage >= 50
            ):
                evidence_class = "Moderate_panel_match"

            else:
                evidence_class = "Weak_panel_match"

        output_rows.append(
            {
                "contig": contig,
                "length": row["length"].strip(),
                "topology": row["topology"].strip(),
                "n_genes": row["n_genes"].strip(),
                "virus_score": row["virus_score"].strip(),
                "n_hallmarks": row["n_hallmarks"].strip(),
                "marker_enrichment": row[
                    "marker_enrichment"
                ].strip(),
                "genomad_taxonomy": row["taxonomy"].strip(),
                "blast_ref": blast_ref,
                "blast_name": blast_name,
                "blast_identity": blast_identity,
                "blast_alignment_length": blast_alignment_length,
                "blast_query_coverage": blast_query_coverage,
                "blast_evalue": blast_evalue,
                "blast_bitscore": blast_bitscore,
                "evidence_class": evidence_class,
            }
        )


# Sort by contig length, largest first
output_rows.sort(
    key=lambda row: int(row["length"]),
    reverse=True,
)


# ============================================================
# Write Linux-style TSV
# ============================================================

OUT.parent.mkdir(parents=True, exist_ok=True)

fields = [
    "contig",
    "length",
    "topology",
    "n_genes",
    "virus_score",
    "n_hallmarks",
    "marker_enrichment",
    "genomad_taxonomy",
    "blast_ref",
    "blast_name",
    "blast_identity",
    "blast_alignment_length",
    "blast_query_coverage",
    "blast_evalue",
    "blast_bitscore",
    "evidence_class",
]

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
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
    label = row["evidence_class"]
    counts[label] = counts.get(label, 0) + 1

print(f"Wrote {len(output_rows)} viral contigs to:")
print(OUT)

print("\nEvidence classes:")
for label, count in sorted(
    counts.items(),
    key=lambda x: x[1],
    reverse=True,
):
    print(f"  {label}: {count}")
