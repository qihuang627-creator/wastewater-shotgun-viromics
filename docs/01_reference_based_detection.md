# Reference-based viral detection

## 1. Quality control

Tools:
- FastQC 0.11.9
- fastp

ERR14789190 contained 408,438 paired reads.

After fastp filtering:
- 408,438 pairs retained
- no reads removed for low quality
- 16,788 reads contained adapter trimming
- R1 Q30: 94.70%
- R2 Q30: 96.49%

## 2. Human-read depletion

Human reference:
GRCh38

Aligner:
Bowtie2, --very-sensitive

Only read pairs for which both mates were unmapped to GRCh38 were retained.

Result:
408,429 / 408,438 read pairs remained.

## 3. Kraken2 screening

Database:
Kraken2 Standard-16, 2026-06-26

Results:
- Classified: 28,871 pairs (7.07%)
- Unclassified: 379,558 pairs (92.93%)
- Viral clade: 14,628 pairs (3.58%)

## 4. Candidate validation

Selected Kraken2 candidates were competitively mapped against RefSeq genomes with Bowtie2.

Coverage was calculated using:

samtools coverage -q 20 -Q 20

Genome breadth was additionally evaluated at ≥1×, ≥5×, and ≥10× depth.

## Interpretation

Kraken2 was treated as a candidate-generation step rather than a final detection method.

Near-complete genome breadth at multiple depth thresholds provided substantially stronger evidence than taxonomic assignment alone.

Closely related tobamoviruses remain subject to additional validation through de novo assembly and sequence-level comparison.
