# De novo viral discovery and hierarchical candidate characterization

## 1. Purpose

This analysis extends the reference-based screening workflow by identifying viral sequences from de novo assembled wastewater metagenomic contigs.

The goal is not to force every sequence into a species-level assignment. Instead, complementary evidence is used to distinguish:

1. sequences strongly related to known viruses;
2. divergent sequences with detectable nucleotide or protein homology;
3. unresolved viral candidates.

The development sample is **ERR14789190** from **PRJEB87273**.

---

## 2. De novo assembly

Human-depleted paired reads were assembled with MEGAHIT v1.2.9.

Input files:

```text
results/host_removal/ERR14789190_R1.nonhuman.fastq.gz
results/host_removal/ERR14789190_R2.nonhuman.fastq.gz
```

Representative command:

```bash
megahit \
  -1 results/host_removal/ERR14789190_R1.nonhuman.fastq.gz \
  -2 results/host_removal/ERR14789190_R2.nonhuman.fastq.gz \
  -t 8 \
  --min-contig-len 500 \
  -o results/assembly/ERR14789190_megahit
```

Assembly summary:

| Metric | Result |
|---|---:|
| Total contigs | 8,545 |
| Total length | 6,357,323 bp |
| Contigs >= 1 kb | 889 |
| Contigs >= 3 kb | 63 |
| Contigs >= 5 kb | 19 |

Contigs >=1 kb were retained for viral discovery.

---

## 3. Viral detection with geNomad

geNomad v1.12.0 with database v1.9 was used to analyze the 889 contigs >=1 kb.

The analysis identified:

**270 viral contigs**

This number represents contigs classified as viral, not 270 distinct virus species.

The detected sequences included members or higher-level assignments associated with groups such as Caudoviricetes, Microviridae, Nodaviridae, Tombusviridae, Picornavirales, Crassvirales, and several RNA-virus-associated lineages.

geNomad is considered **candidate-panel-independent** in this workflow because it does not use the manually selected ten-virus reference panel. It is not completely reference-independent because its predictions rely on trained models and reference information.

---

## 4. Concordance with the reference-based branch

Several assembled contigs related to viruses detected during reference-based screening were also classified as viral by geNomad.

Examples include:

| Contig | Reference-related signal | Contig length | Approx. nucleotide identity |
|---|---|---:|---:|
| k141_3066 | ToMMV-related | 6,421 bp | 99.6% |
| k141_9183 | ToMV-related | 6,424 bp | 99.5% |
| k141_599 | PMMoV-related | 6,761 bp | 99.5% over most of the contig |
| k141_4630 | TMGMV-related | 6,414 bp | 98.8% |
| k141_398 | TMV-related | 6,382 bp | 98.6% |

This provides converging sequence evidence from the reference-based and assembly-based branches.

These results are not treated as unconditional species confirmation.

---

## 5. Read support for viral contigs

Non-human reads were competitively mapped against the 270 geNomad viral contigs with Bowtie2.

The overall alignment rate to the viral-contig set was approximately **14.34%**.

Per-contig coverage was evaluated using samtools.

### High read support

```text
breadth >=1x  >= 95%
breadth >=10x >= 80%
mean depth    >= 10x
```

### Moderate read support

```text
breadth >=1x >= 80%
breadth >=5x >= 50%
mean depth   >= 3x
```

Contigs below these thresholds were categorized as low read support.

| Category | Number of contigs |
|---|---:|
| High | 39 |
| Moderate | 109 |
| Low | 122 |
| Total | 270 |

These categories are operational prioritization criteria rather than formal viral-genome quality standards.

Read-back mapping is also not independent validation because the same reads were used to generate the assembly.

---

## 6. Evidence integration

The script:

```text
scripts/05_integrate_viral_evidence.py
```

combines geNomad output with contig-to-candidate-panel BLAST evidence.

The primary output is:

```text
results_summary/viral_contig_evidence.tsv
```

Candidate-panel nucleotide matches were operationally categorized as follows.

### Strong panel match

```text
identity >= 95%
query coverage >= 90%
```

### Moderate panel match

```text
identity >= 90%
query coverage >= 50%
```

Matches below those thresholds were categorized as weak.

Contigs without a hit against the manually selected candidate-virus panel were labeled:

```text
No_match_to_candidate_panel
```

This label means only that the contig did not match the small candidate panel. It does not indicate that the sequence is novel.

Read-support information was then added using:

```text
scripts/06_add_read_support.py
```

to produce:

```text
results_summary/viral_contig_evidence_with_reads.tsv
```

---

## 7. Selection of high-priority panel-external contigs

A subset was selected for deeper sequence characterization using:

```text
length >= 3 kb
AND
No_match_to_candidate_panel
AND
High_read_support
```

Seventeen contigs met these criteria.

This selection was designed for efficient candidate prioritization in the case study and should not be interpreted as a universal viral-discovery threshold.

---

## 8. Broad nucleotide search

The 17 selected contigs were first queried against NCBI `core_nt` using megablast.

The results were summarized with:

```text
scripts/08_parse_broad_blast.py
```

and stored in:

```text
results_summary/high_priority_broad_nt_summary.tsv
```

The resulting categories were:

| Category | Number |
|---|---:|
| Strong near-known nt match | 2 |
| Moderate divergent nt match | 2 |
| Partial or distant nt match | 2 |
| Weak nt match | 1 |
| No megablast hit | 10 |

Two strong near-known examples were:

```text
k141_2232
k141_6149
```

Because megablast is optimized for relatively similar nucleotide sequences, absence of a megablast hit was not considered evidence of novelty.

---

## 9. Sensitive nucleotide follow-up

Candidates not strongly resolved by megablast were subsequently analyzed with standard BLASTn against `core_nt`.

The stable implementation is:

```text
scripts/09_sensitive_core_nt_blast.sh
```

Fifteen contigs entered this stage.

Results:

- 11 produced nucleotide hits;
- 4 produced no significant nucleotide hit.

The four contigs without a standard BLASTn hit were:

```text
k141_10055
k141_4996
k141_869
k141_899
```

Other contigs produced only partial or distant nucleotide matches and were evaluated according to the strength and coverage of those alignments.

---

## 10. Conditional protein-level characterization

Protein-level analysis was reserved for candidates with absent or weak nucleotide evidence.

Nine contigs were selected for BLASTx against NCBI ClusteredNR (`nr_cluster_seq`).

The corresponding scripts are:

```text
scripts/10_blastx_clusterednr.sh
scripts/10b_retry_single_blastx.sh
```

A three-contig BLASTx batch exceeded the NCBI remote-service CPU limit, so those sequences were retried individually.

Final outcome:

| Protein follow-up | Number |
|---|---:|
| Candidates evaluated | 9 |
| Detectable BLASTx homology | 8 |
| No significant BLASTx hit | 1 |

The unresolved candidate was:

```text
k141_869
```

Examples of informative protein-level evidence included:

- `k141_9867`: approximately 1,365 aa aligned at ~40% amino-acid identity, with E-value 0;
- `k141_14034`, `k141_14035`, and `k141_6353`: approximately 400-aa alignments at ~46-48% amino-acid identity;
- `k141_899`: weaker but detectable protein homology despite the lack of a significant standard BLASTn hit.

Best protein hits are summarized in:

```text
results_summary/ERR14789190_best_blastx_hits.csv
```

---

## 11. Evidence hierarchy

The workflow deliberately uses multiple levels of evidence.

### Near-known sequence

A contig with strong nucleotide similarity, broad query coverage, and consistent read support can be interpreted relatively confidently as closely related to a known reference sequence.

### Divergent sequence

A contig with weak or absent nucleotide similarity but significant protein-level homology may represent a more divergent sequence relative to available nucleotide references.

### Unresolved candidate

A contig with viral-prediction evidence and strong read support but no significant nucleotide or protein similarity remains unresolved.

For example:

```text
k141_869
```

is conservatively described as an **unresolved/divergent viral candidate**.

No claim of a novel virus is made from database no-hit alone.

---

## 12. Limitations

### Viral contigs are not equivalent to viral species

One viral genome may be fragmented across several contigs, while multiple viral contigs can belong to the same species or strain.

### Read remapping is not independent biological validation

The same sequencing reads were used for both assembly and remapping.

### Classifiers may disagree

Kraken2, BLAST, and geNomad use different databases and algorithms. Closely related sequences, incomplete references, and divergent genomes can therefore produce different assignments.

### Database no-hit does not establish novelty

Failure to obtain a significant nucleotide or protein hit can result from database incompleteness, divergence, assembly artifacts, limited informative sequence, or search-method limitations.

### Evidence thresholds are operational

The thresholds used here were designed to prioritize candidates in this workflow and are not universal virological standards.

---

## 13. Transition to multi-sample analysis

The ERR14789190 analysis served as a development case for testing each analytical component.

The multi-sample workflow will use a different structure:

```text
per-sample QC
        |
host depletion
        |
per-sample assembly
        |
viral detection
        |
pool viral contigs across samples
        |
dereplication / clustering
        |
nonredundant viral catalog
        |
map every sample to the common catalog
        |
sample-by-virus abundance matrix
        |
prevalence, taxonomy, and comparative analyses
```

The deeper sequence-characterization branch:

```text
megablast
   -> sensitive BLASTn
   -> conditional BLASTx
```

will be applied only to selected representative sequences from the final catalog rather than repeated for every sequence in every sample.
