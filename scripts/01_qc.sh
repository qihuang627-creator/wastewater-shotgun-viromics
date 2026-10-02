#!/usr/bin/env bash
set -euo pipefail

SAMPLE="ERR14789190"
THREADS=8

mkdir -p qc/fastqc_raw qc/fastp

# Raw-read QC
fastqc \
  raw/${SAMPLE}_1.fastq.gz \
  raw/${SAMPLE}_2.fastq.gz \
  -o qc/fastqc_raw \
  -t ${THREADS}

# Adapter trimming and quality filtering
fastp \
  -i raw/${SAMPLE}_1.fastq.gz \
  -I raw/${SAMPLE}_2.fastq.gz \
  -o qc/fastp/${SAMPLE}_R1.clean.fastq.gz \
  -O qc/fastp/${SAMPLE}_R2.clean.fastq.gz \
  --html qc/fastp/${SAMPLE}.fastp.html \
  --json qc/fastp/${SAMPLE}.fastp.json \
  --thread ${THREADS}
