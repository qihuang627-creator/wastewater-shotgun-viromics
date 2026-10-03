# De novo viral discovery and hierarchical candidate characterization

## 1. Purpose

This analysis extends the reference-based screening workflow by identifying viral sequences directly from de novo assembled wastewater metagenomic contigs.

The objective is not to classify every sequence to species level. Instead, the workflow combines multiple complementary sources of evidence to prioritize viral contigs and distinguish:

1. contigs strongly related to known viruses;
2. divergent but detectable viral sequences;
3. unresolved viral candidates.

The development sample is ERR14789190 from PRJEB87273.

---

## 2. Assembly

Human-depleted paired reads were assembled with MEGAHIT v1.2.9.

Input:

```text
results/host_removal/ERR14789190_R1.nonhuman.fastq.gz
results/host_removal/ERR14789190_R2.nonhuman.fastq.gz

Assembly command:

megahit \
  -1 results/host_removal/ERR14789190_R1.nonhuman.fastq.gz \
  -2 results/host_removal/ERR14789190_R2.nonhuman.fastq.gz \
  -t 8 \
  --min-contig-len 500 \
  -o results/assembly/ERR14789190_megahit

Assembly summary:

Metric	Result
Total contigs	8,545
Total assembly length	6,357,323 bp
Contigs >=1 kb	889
Contigs >=3 kb	63
Contigs >=5 kb	19

Contigs >=1 kb were retained for viral discovery:

results/viral_discovery/ERR14789190_contigs_1kb.fa
3. geNomad viral detection

geNomad v1.12.0 was run in Docker using database v1.9.

Container:

community.wave.seqera.io/library/genomad:1.12.0--27836e6e665e84b5

Representative command:

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$PWD":/work \
  -w /work \
  community.wave.seqera.io/library/genomad:1.12.0--27836e6e665e84b5 \
  genomad end-to-end \
  --threads 8 \
  --splits 8 \
  results/viral_discovery/ERR14789190_contigs_1kb.fa \
  results/genomad/ERR14789190 \
  databases/genomad/genomad_db

geNomad identified:

270 viral contigs

from the 889 contigs >=1 kb.

This number represents viral contigs, not viral species.

The output included diverse predicted viral groups, including Caudoviricetes, Microviridae, Nodaviridae, Tombusviridae, Picornavirales, Crassvirales, and multiple RNA-virus-associated groups.

geNomad detection is described as candidate-panel-independent, because it does not depend on the manually selected ten-virus panel used in the reference-based branch. It is not fully reference-independent because geNomad uses trained models, marker information, and taxonomic databases.

4. Concordance with reference-based candidates

Several contigs previously linked to candidate viruses by reference-based analysis were independently classified as viral by geNomad.

Examples included:

k141_599
k141_9183
k141_3066
k141_4630
k141_398

These contigs were assigned within Riboviria / Kitrinoviricota / Martellivirales by geNomad.

Independent contig-to-reference BLAST analysis had previously shown near-full-length sequence similarity for several of these assemblies:

Contig	Reference-related signal	Contig length	Approx. nucleotide identity
k141_3066	ToMMV-related	6,421 bp	99.6%
k141_9183	ToMV-related	6,424 bp	99.5%
k141_599	PMMoV-related	6,761 bp	99.5% over most of the contig
k141_4630	TMGMV-related	6,414 bp	98.8%
k141_398	TMV-related	6,382 bp	98.6%

These assignments are treated as converging sequence evidence rather than unconditional species confirmation.

5. Read support for viral contigs

Non-human reads were mapped competitively against the 270 geNomad viral contigs.

Bowtie2 overall alignment rate to this viral-contig set was:

14.34%

Per-contig coverage was calculated with samtools.

For prioritization, the following operational classes were used.

High read support
breadth >=1x : >=95%
breadth >=10x: >=80%
mean depth   : >=10x
Moderate read support
breadth >=1x: >=80%
breadth >=5x: >=50%
mean depth  : >=3x
Low read support

Anything below the above thresholds.

Result:

Read-support category	Viral contigs
High	39
Moderate	109
Low	122
Total	270

These thresholds are operational prioritization rules and are not formal viral-genome quality standards.

Read remapping is also not independent biological confirmation because the mapped reads are the same reads from which the assembly was generated. Instead, the mapping step evaluates consistency of read support across each assembled contig.

6. Evidence integration

scripts/05_integrate_viral_evidence.py integrates geNomad results with the candidate-reference BLAST results.

The resulting table is:

results_summary/viral_contig_evidence.tsv

Candidate-panel BLAST matches were operationally divided into:

Strong panel match
nucleotide identity >=95%
query coverage      >=90%
Moderate panel match
nucleotide identity >=90%
query coverage      >=50%
Weak panel match

Hits below those thresholds.

No hit to the ten manually selected candidate references is labeled:

No_match_to_candidate_panel

This label does not mean that the contig is novel or absent from public sequence databases.

scripts/06_add_read_support.py then combines the sequence-level and read-mapping information into:

results_summary/viral_contig_evidence_with_reads.tsv
7. High-priority panel-external contigs

A focused candidate set was generated using:

length >=3 kb
AND
No_match_to_candidate_panel
AND
High_read_support

Seventeen contigs met these criteria.

This subset was used only for deeper characterization. The threshold is a prioritization strategy for this case study rather than a universal viral-discovery rule.

8. Broad nucleotide characterization

The 17 high-priority contigs were first searched using megablast against NCBI core_nt.

Results were summarized with:

scripts/08_parse_broad_blast.py

and written to:

results_summary/high_priority_broad_nt_summary.tsv

Initial classification:

Megablast category	Number
Strong near-known nt match	2
Moderate divergent nt match	2
Partial/distant nt match	2
Weak nt match	1
No megablast hit	10

Two particularly strong near-known examples were:

k141_2232
k141_6149

Both had near-complete query coverage against known nucleotide sequences.

Megablast is optimized for relatively similar sequences, so absence of a megablast hit was not interpreted as evidence of novelty.

9. Sensitive nucleotide follow-up

Candidates not already strongly resolved by megablast were followed with standard BLASTn.

A full nt remote search proved computationally expensive for the remote service, so the stable workflow uses:

standard BLASTn
+
core_nt
+
small sequential batches

The implementation is:

scripts/09_sensitive_core_nt_blast.sh

Among 15 candidates entering this step:

11 produced nucleotide hits
4 produced no significant nucleotide hit

The four candidates without standard BLASTn hits were:

k141_10055
k141_4996
k141_869
k141_899

Some additional candidates had only short or low-similarity nucleotide matches and were therefore retained for protein-level follow-up.

10. Conditional protein-level characterization

Protein-level searches were not applied to every viral contig.

Nine candidates with absent or weak nucleotide-level evidence were selected for BLASTx against NCBI ClusteredNR:

nr_cluster_seq

Scripts:

scripts/10_blastx_clusterednr.sh
scripts/10b_retry_single_blastx.sh

Long contigs were submitted in small batches. Candidates from a batch that exceeded the remote NCBI CPU limit were retried individually.

Final result:

9 candidates investigated
8 with detectable protein-level homology
1 with no significant ClusteredNR BLASTx hit

The unresolved candidate was:

k141_869

Examples of informative protein-level matches included:

k141_9867: a 1,365-aa alignment at approximately 40% amino-acid identity with E-value 0;
k141_14034, k141_14035, and k141_6353: approximately 400-aa protein alignments at ~46-48% identity;
k141_899: weaker but detectable protein-level homology despite the absence of a standard nucleotide BLAST hit.

Best protein hits are summarized in:

results_summary/ERR14789190_best_blastx_hits.csv
11. Interpretation framework

The workflow intentionally uses a hierarchy of evidence.

A contig with:

strong nucleotide similarity
+ near-full query coverage
+ consistent read support

is considered much more readily interpretable than a contig supported only by a taxonomic classifier.

A contig with:

no nucleotide hit
+ detectable protein homology

is treated as potentially divergent relative to currently represented nucleotide sequences.

A contig with:

viral-prediction evidence
+ strong read support
+ no significant nucleotide hit
+ no significant protein hit

is reported conservatively as:

unresolved/divergent viral candidate

It is not described as a novel virus without substantially more evidence.

12. Important limitations
Viral contigs are not viral species

Multiple contigs can derive from a single viral genome, while fragmented assemblies can split one genome into several contigs.

Read remapping is not independent validation

The reads used for read-support analysis are the same reads used for assembly.

Taxonomic classifiers may disagree

Closely related viruses, incomplete references, divergent sequences, and shared conserved regions can produce disagreement between Kraken2, BLAST, and geNomad.

Database absence does not establish novelty

Failure to identify significant nucleotide or protein similarity can result from incomplete databases, highly divergent sequence, short informative regions, assembly artifacts, or limitations of the search strategy.

Operational thresholds are prioritization rules

The evidence classes used in this case study were developed to organize candidates and should not be interpreted as universal virological standards.

13. Role of this single-sample analysis

ERR14789190 served as a development sample for understanding and validating individual analysis modules.

The final multi-sample workflow will not repeat the complete deep-characterization branch independently for every sample.

Instead, the planned multi-sample structure is:

per-sample QC
-> host depletion
-> assembly
-> viral detection
-> pooled viral sequences
-> dereplication / clustering
-> nonredundant viral catalog
-> all samples mapped to the same catalog
-> sample x viral-contig abundance matrix
-> prevalence / taxonomy / ecological comparisons

Detailed megablast, sensitive BLASTn, BLASTx, genome architecture, or phylogenetic analysis will be reserved for selected representative sequences from the final catalog.
