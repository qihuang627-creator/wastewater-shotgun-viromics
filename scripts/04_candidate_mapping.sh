#!/usr/bin/env bash
set -euo pipefail

SAMPLE="ERR14789190"
THREADS=8

REFDIR="refs/candidate_viruses"
PANEL="${REFDIR}/candidate_viruses.fa"
INDEX="${REFDIR}/candidate_viruses"

mkdir -p "${REFDIR}" results/viral_validation results_summary logs

# Download candidate viral RefSeq genomes listed in accessions.tsv
> "${PANEL}"

while IFS=$'\t' read -r taxid name accession; do

    fasta="${REFDIR}/${name}_${accession}.fa"

    if [[ ! -s "${fasta}" ]]; then
        curl -fL --retry 3 \
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text" \
          -o "${fasta}"
        sleep 0.4
    fi

    cat "${fasta}" >> "${PANEL}"

done < "${REFDIR}/accessions.tsv"

# Build Bowtie2 index
bowtie2-build "${PANEL}" "${INDEX}"

# Competitive alignment against all candidate viral references
bowtie2 \
  --very-sensitive \
  -x "${INDEX}" \
  -1 results/host_removal/${SAMPLE}_R1.nonhuman.fastq.gz \
  -2 results/host_removal/${SAMPLE}_R2.nonhuman.fastq.gz \
  -p ${THREADS} \
  2> logs/${SAMPLE}_candidate_mapping.log \
| samtools sort \
  -@ ${THREADS} \
  -o results/viral_validation/${SAMPLE}_candidate_panel.bam

samtools index \
  results/viral_validation/${SAMPLE}_candidate_panel.bam

# Coverage summary using MAPQ >=20 and base quality >=20
samtools coverage \
  -q 20 \
  -Q 20 \
  results/viral_validation/${SAMPLE}_candidate_panel.bam \
  > results/viral_validation/${SAMPLE}_candidate_coverage.tsv

cp \
  results/viral_validation/${SAMPLE}_candidate_coverage.tsv \
  results_summary/candidate_coverage.tsv

# Per-base depth
samtools depth \
  -aa \
  -q 20 \
  -Q 20 \
  results/viral_validation/${SAMPLE}_candidate_panel.bam \
  > results/viral_validation/${SAMPLE}_candidate_depth.tsv

# Genome breadth at >=1x, >=5x and >=10x
awk '
BEGIN {
    OFS="\t";
    print "Reference","Breadth1x","Breadth5x","Breadth10x","MeanDepth","Length"
}
{
    n[$1]++;
    if($3>=1)  c1[$1]++;
    if($3>=5)  c5[$1]++;
    if($3>=10) c10[$1]++;
    sum[$1]+=$3
}
END {
    for (r in n)
        printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n",
        r,
        100*c1[r]/n[r],
        100*c5[r]/n[r],
        100*c10[r]/n[r],
        sum[r]/n[r],
        n[r];
}' results/viral_validation/${SAMPLE}_candidate_depth.tsv \
> results_summary/candidate_breadth.tsv
