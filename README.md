# Wastewater Shotgun Viromics

A reproducible, evidence-aware workflow for **non-target shotgun wastewater metagenomics**, integrating reference-based viral screening with de novo viral discovery and hierarchical sequence characterization.

The current repository contains a completed **single-sample case study** used to develop and validate the workflow. A multi-sample implementation is the next stage of the project.

## Project goals

This project demonstrates how shotgun wastewater metagenomic reads can be processed from raw FASTQ files to interpretable viral signals while avoiding over-reliance on any single classification method.

The workflow combines:

- read quality control and host depletion;
- taxonomic screening;
- reference-based candidate validation;
- de novo metagenomic assembly;
- viral contig detection;
- read-support assessment;
- hierarchical nucleotide and protein homology searches;
- explicit separation of high-confidence, divergent, and unresolved candidates.

The workflow is intended for **non-target metagenomic sequencing**, not hybrid-capture data.

## Workflow

```mermaid
flowchart TD
    A[Paired-end shotgun reads] --> B[fastp QC]
    B --> C[Bowtie2 human depletion]
    C --> D[Kraken2 taxonomic screening]

    C --> E[MEGAHIT de novo assembly]
    E --> F[Contigs >= 1 kb]
    F --> G[geNomad viral detection]
    G --> H[Reads mapped back to viral contigs]
    H --> I[Coverage breadth and depth]
    I --> J[Viral contig evidence table]

    D --> K[Candidate reference panel]
    K --> L[Competitive reference mapping]
    L --> M[Genome breadth/depth]
    E --> N[Contig-to-candidate BLAST]
    M --> J
    N --> J

    J --> O[High-priority contigs]
    O --> P[core_nt megablast]
    P --> Q[Standard BLASTn for unresolved candidates]
    Q --> R[Conditional BLASTx protein follow-up]
The nucleotide/protein characterization branch is optional follow-up, not a step that should be repeated exhaustively for every contig or every sample.

Single-sample case study

Development sample:

Field	Value
BioProject	PRJEB87273
Run	ERR14789190
Library strategy	WGS / METAGENOMIC
Library selection	RANDOM
Layout	PAIRED
Capture	Non-target / non-capture

The sample was used as a workflow-development and smoke-test dataset. No specific virus was assumed to be present a priori.

Key results
Stage	Result
Input read pairs	408,438
Non-human read pairs retained	408,429
Kraken2 classified pairs	7.07%
Kraken2 virus-classified pairs	3.58%
MEGAHIT contigs	8,545
Assembly size	6.36 Mb
Contigs >=1 kb	889
Contigs >=3 kb	63
geNomad viral contigs	270
High read-support viral contigs	39
High-priority >=3 kb panel-external contigs	17

Kraken2 provided an initial screening layer, but taxonomic calls were not treated as final evidence.

Reference-based mapping and de novo assembly provided converging support for several candidate viral genomes. Near-full-length assembled contigs showed strong similarity to reference sequences related to PMMoV, ToMMV, ToMV, TMGMV, and TMV.

geNomad independently detected viral signal on the de novo assembled contigs without using the manually selected candidate panel. Because geNomad itself uses trained models and reference databases, this is described here as candidate-panel-independent detection, rather than fully reference-independent detection.

Viral contig prioritization

The 270 geNomad viral contigs were integrated with:

geNomad score and taxonomy;
candidate-panel BLAST similarity;
read-mapping breadth;
read depth.

Read support was divided into operational prioritization classes:

High_read_support: >=95% breadth at 1x, >=80% breadth at 10x, and mean depth >=10x;
Moderate_read_support: >=80% breadth at 1x, >=50% breadth at 5x, and mean depth >=3x;
Low_read_support: below those thresholds.

These thresholds are used for prioritization only and are not virus-species confirmation criteria.

This produced:

39 high-read-support contigs;
109 moderate-read-support contigs;
122 low-read-support contigs.

Seventeen contigs were selected as high-priority panel-external candidates because they were >=3 kb, had high read support, and did not match the original manually defined candidate panel.

Hierarchical sequence characterization

Rather than performing the most computationally expensive search on every viral contig, high-priority candidates were characterized hierarchically:

core_nt megablast for near-known nucleotide matches;
standard BLASTn for candidates unresolved or weakly resolved by megablast;
BLASTx against ClusteredNR (nr_cluster_seq) only for candidates with absent or weak nucleotide-level evidence.

Among the 17 high-priority panel-external contigs, megablast identified strong near-known matches for two contigs and more distant nucleotide relationships for several others.

Nine candidates were selected for protein-level follow-up. Eight produced detectable BLASTx protein homology, while one (k141_869) remained unresolved at both the sensitive nucleotide and protein-search levels.

An unresolved result is not interpreted as evidence of a novel virus.

Repository structure
.
├── README.md
├── docs/
│   ├── 01_reference_based_detection.md
│   └── 02_de_novo_viral_discovery.md
├── refs/
│   └── candidate_viruses/
├── results_summary/
│   ├── candidate_breadth.tsv
│   ├── candidate_coverage.tsv
│   ├── viral_contig_evidence.tsv
│   ├── viral_contig_read_support.tsv
│   ├── viral_contig_evidence_with_reads.tsv
│   ├── high_priority_broad_nt_summary.tsv
│   └── ERR14789190_best_blastx_hits.csv
├── scripts/
│   ├── 01_qc.sh
│   ├── 02_host_removal.sh
│   ├── 03_kraken2.sh
│   ├── 04_candidate_mapping.sh
│   ├── 05_integrate_viral_evidence.py
│   ├── 06_add_read_support.py
│   ├── 07_broad_nt_blast.sh
│   ├── 08_parse_broad_blast.py
│   ├── 09_sensitive_core_nt_blast.sh
│   ├── 10_blastx_clusterednr.sh
│   └── 10b_retry_single_blastx.sh
└── software_versions.txt

Large raw data, databases, BAM files, complete assembly outputs, complete BLAST outputs, logs, and temporary batch files are excluded from version control.

Interpretation principles

Several distinctions are important throughout this workflow:

Unclassified reads are not equivalent to non-viral reads. Environmental wastewater contains sequences absent from reduced taxonomic databases.

Kraken2 assignments are screening signals, not species confirmation. Closely related viruses may share sequence and produce ambiguous assignments.

Read remapping is not independent biological validation. The same reads were used for assembly; remapping is used to evaluate coverage consistency and identify poorly supported assemblies.

A viral contig is not equivalent to a virus species. Multiple contigs may originate from the same genome, and a single genome may be fragmented.

No database hit does not demonstrate novelty. Candidates without significant nucleotide or protein similarity are reported conservatively as unresolved or divergent.

Documentation

Detailed reference-based validation:

Reference-based viral detection

De novo viral discovery, evidence integration, and hierarchical characterization:

De novo viral discovery
Next stage: multi-sample analysis

The current single-sample analysis was intentionally used to test individual analytical modules.

The multi-sample workflow will be simplified into a reproducible core pipeline:

samples
  -> QC
  -> host depletion
  -> assembly
  -> viral detection
  -> pooled viral catalog
  -> dereplication / clustering
  -> all-sample mapping to a common catalog
  -> sample x viral-contig abundance matrix
  -> prevalence and taxonomic analysis

Deep megablast -> sensitive BLASTn -> BLASTx characterization will remain an optional branch for selected representative viral sequences rather than being repeated independently for every sample.
