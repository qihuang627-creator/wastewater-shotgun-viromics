#!/usr/bin/env bash
set -euo pipefail

SAMPLE="ERR14789190"
THREADS=8
DB="databases/kraken2_standard16_20260626"

mkdir -p results/kraken2 logs

kraken2 \
  --db "${DB}" \
  --paired \
  --gzip-compressed \
  --threads ${THREADS} \
  --use-names \
  --report results/kraken2/${SAMPLE}.kreport \
  --output results/kraken2/${SAMPLE}.kraken \
  results/host_removal/${SAMPLE}_R1.nonhuman.fastq.gz \
  results/host_removal/${SAMPLE}_R2.nonhuman.fastq.gz \
  2> logs/${SAMPLE}_kraken2.log
