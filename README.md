# Wastewater Shotgun Viromics

A reproducible workflow for detecting and validating viral signals from untargeted wastewater shotgun metagenomic sequencing.

## Project objective

This project evaluates a reference-based and reference-independent workflow for wastewater viral metagenomics. The current analysis uses a publicly available paired-end wastewater shotgun metagenome (ERR14789190) as a test dataset.

The workflow is designed to distinguish preliminary taxonomic assignments from viral detections supported by genome-wide sequence coverage.

## Workflow

Raw paired-end reads
→ quality control
→ human-read depletion
→ Kraken2 taxonomic screening
→ candidate virus selection
→ competitive reference mapping
→ genome breadth/depth validation
→ de novo assembly and viral contig analysis

## Dataset

- Accession: ERR14789190
- Sequencing: paired-end shotgun metagenomics
- Library strategy: WGS
- Library source: METAGENOMIC
- Selection: RANDOM

After quality filtering, 408,438 read pairs were retained.

Human-reference screening against GRCh38 removed only 9 read pairs, leaving 408,429 non-human read pairs for downstream analysis.

## Reference-based viral screening

Kraken2 Standard-16 classified:

- 28,871 / 408,429 read pairs (7.07%)
- 379,558 read pairs (92.93%) remained unclassified
- 14,628 read pairs (3.58% of all pairs) were assigned within the viral taxonomic clade

Major candidates included several tobamoviruses, Picalivirus A, Shahe picorna-like virus 8, and a crAss-like bacteriophage signal.

Because Kraken2 assignments alone do not establish species-level detection, candidate viruses were subsequently evaluated using competitive reference mapping.

## Competitive mapping validation

Reads were competitively aligned against selected RefSeq viral genomes. Alignments with mapping quality <20 or base quality <20 were excluded from coverage calculations.

Several candidates showed near-complete genome coverage, including:

| Candidate | RefSeq | Breadth ≥1× | Breadth ≥10× | Mean depth |
|---|---|---:|---:|---:|
| Pepper mild mottle virus | NC_003630.1 | 99.98% | 99.18% | 142.84× |
| Tomato mottle mosaic virus | NC_022230.1 | 99.86% | 98.05% | 73.06× |
| Picalivirus A | NC_040594.1 | 99.67% | 99.33% | 50.44× |
| Tobacco mosaic virus | NC_001367.1 | 99.61% | 94.73% | 43.04× |
| Tomato mosaic virus | NC_002692.1 | 99.92% | 94.02% | 45.37× |
| Tobacco mild green mosaic virus | NC_001556.1 | 98.71% | 91.82% | 22.52× |

In contrast, the crAss001 reference showed only 30.88% breadth at ≥1× and 0.07% breadth at ≥10×, illustrating why Kraken2 assignments require independent validation.

## Current interpretation

The reference-based analysis identifies several strong whole-genome viral signals. However, closely related tobamoviruses require additional validation because sequence similarity among related references can complicate species-level interpretation.

The next stage uses de novo assembly to evaluate whether reference-independent contigs independently support these detections.

## Status

Reference-based screening and coverage validation: completed.

De novo assembly and viral contig analysis: in progress.
