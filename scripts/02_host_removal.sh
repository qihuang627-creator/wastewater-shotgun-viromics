#!/usr/bin/env bash
set -euo pipefail

SAMPLE="ERR14789190"
THREADS=8
DB="databases/human_GRCh38/GRCh38"

mkdir -p results/host_removal logs

# Map cleaned reads against GRCh38
bowtie2 \
  --very-sensitive \
  -x "${DB}" \
  -1 qc/fastp/${SAMPLE}_R1.clean.fastq.gz \
  -2 qc/fastp/${SAMPLE}_R2.clean.fastq.gz \
  -p ${THREADS} \
  2> logs/${SAMPLE}_bowtie2_GRCh38.log \
| samtools view -b \
  -o results/host_removal/${SAMPLE}_vs_GRCh38.bam -

# Alignment summary
samtools flagstat \
  results/host_removal/${SAMPLE}_vs_GRCh38.bam \
  > results/host_removal/${SAMPLE}_GRCh38.flagstat.txt

# Retain pairs for which both mates are unmapped to the human genome.
# -f 12: read unmapped + mate unmapped
# -F 2304: exclude secondary and supplementary alignments
samtools view \
  -b \
  -f 12 \
  -F 2304 \
  results/host_removal/${SAMPLE}_vs_GRCh38.bam \
| samtools sort \
  -n \
  -@ ${THREADS} \
  -o results/host_removal/${SAMPLE}_both_unmapped.namesort.bam -

# Convert non-human pairs back to FASTQ
samtools fastq \
  -@ ${THREADS} \
  -n \
  -1 results/host_removal/${SAMPLE}_R1.nonhuman.fastq.gz \
  -2 results/host_removal/${SAMPLE}_R2.nonhuman.fastq.gz \
  -0 /dev/null \
  -s /dev/null \
  results/host_removal/${SAMPLE}_both_unmapped.namesort.bam
