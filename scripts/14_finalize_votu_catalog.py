#!/usr/bin/env python3

import sys
from collections import OrderedDict

if len(sys.argv) != 6:
    sys.exit(
        "Usage: 14_finalize_votu_catalog.py "
        "<pooled_fasta> <clusters.tsv> <catalog.fna> "
        "<membership.tsv> <representatives.tsv>"
    )

fasta_file, cluster_file, catalog_file, membership_file, repmap_file = sys.argv[1:]

# ---------------------------------------------------------
# Read pooled FASTA
# ---------------------------------------------------------

seqs = OrderedDict()

current = None
chunks = []

with open(fasta_file) as fh:
    for line in fh:
        line = line.rstrip()

        if line.startswith(">"):
            if current is not None:
                seqs[current] = "".join(chunks)

            current = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)

    if current is not None:
        seqs[current] = "".join(chunks)

print(f"Pooled FASTA sequences : {len(seqs)}")

# ---------------------------------------------------------
# Read aniclust output
#
# Format:
# column 1 = representative
# column 2 = comma-separated cluster members
# ---------------------------------------------------------

clusters = OrderedDict()
assigned = {}

with open(cluster_file) as fh:

    for lineno, line in enumerate(fh, 1):

        line = line.rstrip("\n")

        if not line:
            continue

        fields = line.split("\t")

        if len(fields) < 2:
            raise RuntimeError(
                f"Unexpected aniclust format at line {lineno}: {line}"
            )

        rep = fields[0].strip()

        if rep not in seqs:
            raise RuntimeError(
                f"Representative not found in FASTA at line {lineno}: {rep}"
            )

        members = [
            x.strip()
            for x in fields[1].split(",")
            if x.strip()
        ]

        # Include representative itself
        all_members = [rep] + members

        # Remove duplicates but preserve order
        unique_members = []
        seen_here = set()

        for member in all_members:

            if member in seen_here:
                continue

            seen_here.add(member)

            if member not in seqs:
                raise RuntimeError(
                    f"Cluster member not found in FASTA at line "
                    f"{lineno}: {member}"
                )

            if member in assigned and assigned[member] != rep:
                raise RuntimeError(
                    f"Sequence assigned to multiple clusters: "
                    f"{member} -> {assigned[member]} and {rep}"
                )

            assigned[member] = rep
            unique_members.append(member)

        clusters[rep] = unique_members

# ---------------------------------------------------------
# Add any sequences not reported by aniclust as singleton
# clusters. This makes the catalog exhaustive.
# ---------------------------------------------------------

n_singletons_added = 0

for seq_id in seqs:

    if seq_id not in assigned:
        clusters[seq_id] = [seq_id]
        assigned[seq_id] = seq_id
        n_singletons_added += 1

# ---------------------------------------------------------
# Validation
# ---------------------------------------------------------

if len(assigned) != len(seqs):
    raise RuntimeError(
        f"Assignment mismatch: {len(assigned)} assigned vs "
        f"{len(seqs)} FASTA sequences"
    )

# Deterministic ordering:
# longest representative first, then sequence ID
representatives = sorted(
    clusters.keys(),
    key=lambda rep: (-len(seqs[rep]), rep)
)

rep_to_votu = {
    rep: f"vOTU_{i:05d}"
    for i, rep in enumerate(representatives, 1)
}

# ---------------------------------------------------------
# Write representative catalog FASTA
# ---------------------------------------------------------

with open(catalog_file, "w") as out:

    for rep in representatives:

        votu = rep_to_votu[rep]
        sequence = seqs[rep]

        out.write(f">{votu} representative={rep}\n")

        for i in range(0, len(sequence), 80):
            out.write(sequence[i:i+80] + "\n")

# ---------------------------------------------------------
# Write representative table
# ---------------------------------------------------------

with open(repmap_file, "w") as out:

    out.write(
        "votu_id\t"
        "representative_id\t"
        "representative_sample\t"
        "representative_contig\t"
        "length_bp\t"
        "cluster_size\t"
        "n_source_samples\t"
        "source_samples\n"
    )

    for rep in representatives:

        votu = rep_to_votu[rep]
        members = clusters[rep]

        if "|" in rep:
            rep_sample, rep_contig = rep.split("|", 1)
        else:
            rep_sample, rep_contig = "", rep

        source_samples = sorted({
            m.split("|", 1)[0]
            for m in members
            if "|" in m
        })

        out.write(
            f"{votu}\t"
            f"{rep}\t"
            f"{rep_sample}\t"
            f"{rep_contig}\t"
            f"{len(seqs[rep])}\t"
            f"{len(members)}\t"
            f"{len(source_samples)}\t"
            f"{','.join(source_samples)}\n"
        )

# ---------------------------------------------------------
# Write complete membership table
# ---------------------------------------------------------

with open(membership_file, "w") as out:

    out.write(
        "votu_id\t"
        "representative_id\t"
        "member_id\t"
        "member_sample\t"
        "member_contig\t"
        "cluster_size\n"
    )

    for rep in representatives:

        votu = rep_to_votu[rep]
        members = clusters[rep]
        size = len(members)

        for member in members:

            if "|" in member:
                sample, contig = member.split("|", 1)
            else:
                sample, contig = "", member

            out.write(
                f"{votu}\t"
                f"{rep}\t"
                f"{member}\t"
                f"{sample}\t"
                f"{contig}\t"
                f"{size}\n"
            )

print(f"Clusters from aniclust  : {len(clusters) - n_singletons_added}")
print(f"Singletons added        : {n_singletons_added}")
print(f"Assigned sequences      : {len(assigned)}")
print(f"Final vOTUs             : {len(representatives)}")
print("Catalog complete.")
