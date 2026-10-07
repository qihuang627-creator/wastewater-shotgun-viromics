configfile: "config/config.yaml"

# ============================================================
# Wastewater shotgun virome workflow
#
# FASTQ
#   -> QC
#   -> human depletion
#   -> assembly
#   -> viral discovery
#   -> multisample vOTU catalog
#   -> common-catalog mapping
#   -> abundance matrices
#   -> comparative analysis
#   -> ordination
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

SAMPLES = config["samples"]
THREADS = int(config.get("threads", 8))

P = config["paths"]
TOOLS = config["tools"]

RAW = P["raw"]
FASTP_DIR = P["fastp"]
HOST_DIR = P["host_removal"]
ASSEMBLY_DIR = P["assembly"]
VIRAL_DIR = P["viral_discovery"]
GENOMAD_DIR = P["genomad"]
MULTI_DIR = P["multisample"]

FASTQC_DIR = "qc/fastqc_raw"

HUMAN_INDEX = config["human_bowtie2_index"]
GENOMAD_DB = config["genomad_db"]

GENOMAD_IMAGE = config["genomad"]["docker_image"]
GENOMAD_SPLITS = int(config["genomad"]["splits"])

MIN_CONTIG = int(config["votu"]["min_contig_length"])
MIN_ANI = float(config["votu"]["min_ani"])
MIN_AF = float(config["votu"]["min_alignment_fraction"])

MIN_READ_ID = float(config["mapping"]["min_read_identity"])
MIN_READ_ALIGNED = float(
    config["mapping"]["min_read_aligned_percent"]
)

PRESENCE_BREADTH = float(
    config["mapping"]["presence_breadth"]
)
PRESENCE_TAG = int(round(PRESENCE_BREADTH * 100))


# ------------------------------------------------------------
# Important paths
# ------------------------------------------------------------

GENOMAD_VIRUS = (
    GENOMAD_DIR
    + "/{sample}/{sample}_contigs_1kb_summary/"
      "{sample}_contigs_1kb_virus.fna"
)

POOLED = MULTI_DIR + "/pooled/all_genomad_viruses.prefixed.fna"
IDMAP = MULTI_DIR + "/pooled/all_genomad_viruses.id_map.tsv"

CAT_DIR = MULTI_DIR + "/catalog"

BLAST_DB = CAT_DIR + "/viral_pool_db"
BLAST_DB_DONE = CAT_DIR + "/viral_pool_db.done"

BLAST_OUT = CAT_DIR + "/all_vs_all_blast.tsv"
ANI_OUT = CAT_DIR + "/all_vs_all_ani.tsv"
CLUSTERS = CAT_DIR + "/votu_clusters.tsv"

VOTU_CATALOG = CAT_DIR + "/votu_catalog.fna"
VOTU_MEMBERSHIP = CAT_DIR + "/votu_membership.tsv"
VOTU_REPS = CAT_DIR + "/votu_representatives.tsv"

ABUND_DIR = MULTI_DIR + "/abundance"
COMP_DIR = MULTI_DIR + "/comparative"
ORD_DIR = MULTI_DIR + "/ordination"


# ------------------------------------------------------------
# Input collections
# ------------------------------------------------------------

HUMAN_INDEX_FILES = [
    HUMAN_INDEX + ".1.bt2",
    HUMAN_INDEX + ".2.bt2",
    HUMAN_INDEX + ".3.bt2",
    HUMAN_INDEX + ".4.bt2",
    HUMAN_INDEX + ".rev.1.bt2",
    HUMAN_INDEX + ".rev.2.bt2",
]

NONHUMAN_READS = []

for sample in SAMPLES:
    NONHUMAN_READS.extend([
        HOST_DIR + f"/{sample}_R1.nonhuman.fastq.gz",
        HOST_DIR + f"/{sample}_R2.nonhuman.fastq.gz",
    ])


# ============================================================
# Final targets
# ============================================================

ALL_TARGETS = (
    expand(
        [
            FASTQC_DIR + "/{sample}_1_fastqc.html",
            FASTQC_DIR + "/{sample}_2_fastqc.html",
        ],
        sample=SAMPLES,
    )
    + [
        # Main vOTU catalog
        VOTU_CATALOG,
        VOTU_MEMBERSHIP,
        VOTU_REPS,

        # Abundance matrices
        ABUND_DIR + "/votu_count_matrix.tsv",
        ABUND_DIR + "/votu_mean_matrix.tsv",
        ABUND_DIR + "/votu_covered_fraction_matrix.tsv",
        ABUND_DIR + "/votu_rpkm_matrix.tsv",
        ABUND_DIR + "/votu_tpm_matrix.tsv",
        ABUND_DIR + "/sample_mapping_summary.tsv",

        # Comparative outputs
        COMP_DIR + "/multisample_summary.tsv",
        COMP_DIR + "/votu_prevalence.tsv",
        COMP_DIR + "/prevalence_distribution.tsv",
        COMP_DIR + "/sample_votu_richness.tsv",
        COMP_DIR + "/pairwise_jaccard.tsv",

        # Ordination outputs
        ORD_DIR + "/jaccard_pcoa.tsv",
        ORD_DIR + "/braycurtis_sqrt_tpm_pcoa.tsv",
        ORD_DIR + "/ordination_summary.tsv",
    ]
)


rule all:
    input:
        ALL_TARGETS


# ============================================================
# 1. Raw-read FastQC
# ============================================================

rule fastqc_raw:
    input:
        r1=RAW + "/{sample}_1.fastq.gz",
        r2=RAW + "/{sample}_2.fastq.gz"
    output:
        r1_html=FASTQC_DIR + "/{sample}_1_fastqc.html",
        r1_zip=FASTQC_DIR + "/{sample}_1_fastqc.zip",
        r2_html=FASTQC_DIR + "/{sample}_2_fastqc.html",
        r2_zip=FASTQC_DIR + "/{sample}_2_fastqc.zip"
    threads:
        THREADS
    params:
        fastqc=TOOLS["fastqc"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p {FASTQC_DIR}

        "{params.fastqc}" \
            {input.r1} \
            {input.r2} \
            -o {FASTQC_DIR} \
            -t {threads}
        """


# ============================================================
# 2. fastp
# ============================================================

rule fastp:
    input:
        r1=RAW + "/{sample}_1.fastq.gz",
        r2=RAW + "/{sample}_2.fastq.gz"
    output:
        r1=FASTP_DIR + "/{sample}_R1.clean.fastq.gz",
        r2=FASTP_DIR + "/{sample}_R2.clean.fastq.gz",
        html=FASTP_DIR + "/{sample}.fastp.html",
        json=FASTP_DIR + "/{sample}.fastp.json"
    threads:
        THREADS
    params:
        fastp=TOOLS["fastp"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p {FASTP_DIR}

        "{params.fastp}" \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads}
        """


# ============================================================
# 3. Human host removal
# ============================================================

rule host_removal:
    input:
        r1=FASTP_DIR + "/{sample}_R1.clean.fastq.gz",
        r2=FASTP_DIR + "/{sample}_R2.clean.fastq.gz",
        index=HUMAN_INDEX_FILES
    output:
        bam=HOST_DIR + "/{sample}_vs_GRCh38.bam",
        flagstat=HOST_DIR + "/{sample}_GRCh38.flagstat.txt",
        r1=HOST_DIR + "/{sample}_R1.nonhuman.fastq.gz",
        r2=HOST_DIR + "/{sample}_R2.nonhuman.fastq.gz"
    threads:
        THREADS
    params:
        bowtie2=TOOLS["bowtie2"],
        samtools=TOOLS["samtools"],
        index=HUMAN_INDEX,
        namesort=lambda wc: (
            HOST_DIR
            + f"/{wc.sample}_both_unmapped.namesort.bam"
        )
    shell:
        r"""
        set -euo pipefail

        mkdir -p {HOST_DIR}

        "{params.bowtie2}" \
            --very-sensitive \
            -x "{params.index}" \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
        | "{params.samtools}" view -b \
            -o {output.bam} -

        "{params.samtools}" flagstat \
            {output.bam} \
            > {output.flagstat}

        "{params.samtools}" view \
            -b \
            -f 12 \
            -F 2304 \
            {output.bam} \
        | "{params.samtools}" sort \
            -n \
            -@ {threads} \
            -o {params.namesort} -

        "{params.samtools}" fastq \
            -@ {threads} \
            -n \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null \
            -s /dev/null \
            {params.namesort}

        rm -f {params.namesort}
        """


# ============================================================
# 4. MEGAHIT assembly
# ============================================================

rule megahit:
    input:
        r1=HOST_DIR + "/{sample}_R1.nonhuman.fastq.gz",
        r2=HOST_DIR + "/{sample}_R2.nonhuman.fastq.gz"
    output:
        contigs=(
            ASSEMBLY_DIR
            + "/{sample}_megahit/final.contigs.fa"
        )
    threads:
        THREADS
    params:
        megahit=TOOLS["megahit"],
        outdir=lambda wc: (
            ASSEMBLY_DIR
            + f"/{wc.sample}_megahit"
        )
    shell:
        r"""
        set -euo pipefail

        mkdir -p {ASSEMBLY_DIR}

        rm -rf "{params.outdir}"

        "{params.megahit}" \
            -1 {input.r1} \
            -2 {input.r2} \
            -o "{params.outdir}" \
            --min-contig-len 500 \
            -t {threads}
        """


# ============================================================
# 5. Retain contigs >= 1 kb
# ============================================================

rule filter_contigs_1kb:
    input:
        contigs=(
            ASSEMBLY_DIR
            + "/{sample}_megahit/final.contigs.fa"
        )
    output:
        filtered=(
            VIRAL_DIR
            + "/{sample}_contigs_1kb.fa"
        )
    run:

        import os

        os.makedirs(VIRAL_DIR, exist_ok=True)

        min_len = MIN_CONTIG

        with open(input.contigs) as src, \
             open(output.filtered, "w") as dst:

            header = None
            seq = []

            def write_record(h, s):
                if h is None:
                    return

                sequence = "".join(s)

                if len(sequence) >= min_len:
                    dst.write(h + "\n")
                    dst.write(sequence + "\n")

            for line in src:
                line = line.rstrip()

                if line.startswith(">"):
                    write_record(header, seq)
                    header = line
                    seq = []
                else:
                    seq.append(line)

            write_record(header, seq)


# ============================================================
# 6. geNomad de novo viral discovery
# ============================================================

rule genomad:
    input:
        contigs=(
            VIRAL_DIR
            + "/{sample}_contigs_1kb.fa"
        )
    output:
        virus=GENOMAD_VIRUS
    threads:
        THREADS
    params:
        docker=TOOLS["docker"],
        image=GENOMAD_IMAGE,
        database=GENOMAD_DB,
        splits=GENOMAD_SPLITS,
        outdir=lambda wc: (
            GENOMAD_DIR
            + f"/{wc.sample}"
        )
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}"

        DB_ABS="$(cd "{params.database}" && pwd -P)"

        "{params.docker}" run --rm \
            -u "$(id -u):$(id -g)" \
            -v "$PWD:$PWD" \
            -v "${{DB_ABS}}:${{DB_ABS}}:ro" \
            -w "$PWD" \
            "{params.image}" \
            genomad end-to-end \
            {input.contigs} \
            "{params.outdir}" \
            "${{DB_ABS}}" \
            --threads {threads} \
            --splits {params.splits}
        """


# ============================================================
# 7. Pool all geNomad viral contigs
# ============================================================

rule pool_viral_contigs:
    input:
        expand(
            GENOMAD_VIRUS,
            sample=SAMPLES
        )
    output:
        fasta=POOLED,
        idmap=IDMAP
    run:

        import os

        os.makedirs(
            os.path.dirname(output.fasta),
            exist_ok=True
        )

        with open(output.fasta, "w") as outfa, \
             open(output.idmap, "w") as outmap:

            outmap.write(
                "prefixed_id\t"
                "sample\t"
                "original_contig\n"
            )

            for sample, fasta in zip(
                SAMPLES,
                list(input)
            ):

                with open(fasta) as fh:

                    for line in fh:

                        if line.startswith(">"):

                            original = (
                                line[1:]
                                .strip()
                                .split()[0]
                            )

                            prefixed = (
                                sample
                                + "|"
                                + original
                            )

                            outfa.write(
                                ">"
                                + prefixed
                                + "\n"
                            )

                            outmap.write(
                                prefixed
                                + "\t"
                                + sample
                                + "\t"
                                + original
                                + "\n"
                            )

                        else:
                            outfa.write(line)


# ============================================================
# 8. BLAST database
# ============================================================

rule build_votu_blast_db:
    input:
        POOLED
    output:
        BLAST_DB_DONE
    params:
        makeblastdb=TOOLS["makeblastdb"],
        db=BLAST_DB
    shell:
        r"""
        set -euo pipefail

        mkdir -p {CAT_DIR}

        "{params.makeblastdb}" \
            -in {input} \
            -dbtype nucl \
            -out {params.db}

        touch {output}
        """


# ============================================================
# 9. All-vs-all BLASTn
# ============================================================

rule all_vs_all_blast:
    input:
        fasta=POOLED,
        dbdone=BLAST_DB_DONE
    output:
        BLAST_OUT
    threads:
        THREADS
    params:
        blastn=TOOLS["blastn"],
        db=BLAST_DB
    shell:
        r"""
        set -euo pipefail

        "{params.blastn}" \
            -query {input.fasta} \
            -db {params.db} \
            -outfmt '6 std qlen slen' \
            -max_target_seqs 10000 \
            -num_threads {threads} \
            -out {output}
        """


# ============================================================
# 10. Pairwise ANI calculation
# ============================================================

rule anicalc:
    input:
        BLAST_OUT
    output:
        ANI_OUT
    params:
        command=TOOLS["anicalc"]
    shell:
        r"""
        set -euo pipefail

        {params.command} \
            -i {input} \
            -o {output}
        """


# ============================================================
# 11. Species-level vOTU clustering
#
# >=95% ANI
# >=85% alignment fraction
# ============================================================

rule aniclust:
    input:
        fasta=POOLED,
        ani=ANI_OUT
    output:
        CLUSTERS
    params:
        command=TOOLS["aniclust"]
    shell:
        r"""
        set -euo pipefail

        {params.command} \
            --fna {input.fasta} \
            --ani {input.ani} \
            --out {output} \
            --min_ani {MIN_ANI} \
            --min_tcov {MIN_AF} \
            --min_qcov 0
        """


# ============================================================
# 12. Finalize non-redundant vOTU catalog
# ============================================================

rule finalize_votu_catalog:
    input:
        fasta=POOLED,
        clusters=CLUSTERS
    output:
        catalog=VOTU_CATALOG,
        membership=VOTU_MEMBERSHIP,
        representatives=VOTU_REPS
    params:
        python=TOOLS.get("python", "python3")
    shell:
        r"""
        set -euo pipefail

        "{params.python}" \
            scripts/multisample/finalize_votu_catalog.py \
            {input.fasta} \
            {input.clusters} \
            {output.catalog} \
            {output.membership} \
            {output.representatives}
        """


# ============================================================
# 13. Map all samples against common catalog
# ============================================================

rule coverm_common_catalog:
    input:
        catalog=VOTU_CATALOG,
        reads=NONHUMAN_READS
    output:
        ABUND_DIR + "/votu_coverm_long.tsv"
    threads:
        THREADS
    params:
        coverm=TOOLS["coverm"],
        coverm_env=config["coverm_env"],
        reads=" ".join(NONHUMAN_READS)
    shell:
        r"""
        set -euo pipefail

        mkdir -p {ABUND_DIR}

        conda run -n {params.coverm_env} coverm contig \
            -r {input.catalog} \
            -c {params.reads} \
            -p minimap2-sr \
            --min-read-percent-identity {MIN_READ_ID} \
            --min-read-aligned-percent {MIN_READ_ALIGNED} \
            --exclude-supplementary \
            --contig-end-exclusion 0 \
            --methods \
                count \
                mean \
                covered_fraction \
                rpkm \
                tpm \
            --output-format sparse \
            -t {threads} \
            -o {output}
        """


# ============================================================
# 14. Build abundance matrices
# ============================================================

rule abundance_matrices:
    input:
        ABUND_DIR + "/votu_coverm_long.tsv"
    output:
        count_matrix=ABUND_DIR + "/votu_count_matrix.tsv",
        mean=ABUND_DIR + "/votu_mean_matrix.tsv",
        breadth=ABUND_DIR + "/votu_covered_fraction_matrix.tsv",
        rpkm=ABUND_DIR + "/votu_rpkm_matrix.tsv",
        tpm=ABUND_DIR + "/votu_tpm_matrix.tsv",
        summary=ABUND_DIR + "/sample_mapping_summary.tsv"
    params:
        samples=" ".join(SAMPLES)
    shell:
        r"""
        set -euo pipefail

        python3 scripts/multisample/build_votu_matrices.py \
            --input "{input}" \
            --outdir "{ABUND_DIR}" \
            --samples {params.samples}
        """


# ============================================================
# 15. Multisample comparison
# ============================================================

rule comparative_analysis:
    input:
        ABUND_DIR + "/votu_covered_fraction_matrix.tsv"
    output:
        presence=(COMP_DIR + f"/votu_presence_breadth{PRESENCE_TAG}.tsv"),
        prevalence=COMP_DIR + "/votu_prevalence.tsv",
        distribution=COMP_DIR + "/prevalence_distribution.tsv",
        richness=COMP_DIR + "/sample_votu_richness.tsv",
        jaccard=COMP_DIR + "/pairwise_jaccard.tsv",
        summary=COMP_DIR + "/multisample_summary.tsv",
        specific=COMP_DIR + "/sample_specific_votus.txt",
        shared=COMP_DIR + "/shared_votus.txt",
        core=COMP_DIR + "/core_votus.txt"
    shell:
        r"""
        set -euo pipefail

        python3 scripts/multisample/multisample_comparison.py \
            --input "{input}" \
            --outdir "{COMP_DIR}" \
            --threshold {PRESENCE_BREADTH}

        # Normalize line endings for Linux downstream tools
        sed -i 's/\r$//' {COMP_DIR}/*.tsv
        """


# ============================================================
# 16. Jaccard + Bray-Curtis ordination / clustering
# ============================================================

rule ordination:
    input:
        breadth=(
            ABUND_DIR
            + "/votu_covered_fraction_matrix.tsv"
        ),
        tpm=(
            ABUND_DIR
            + "/votu_tpm_matrix.tsv"
        ),
        comparison=(
            COMP_DIR
            + "/multisample_summary.tsv"
        )
    output:
        jaccard_distance=(
            ORD_DIR
            + "/jaccard_distance.tsv"
        ),
        bray_distance=(
            ORD_DIR
            + "/braycurtis_sqrt_tpm_distance.tsv"
        ),
        jaccard_pcoa=(
            ORD_DIR
            + "/jaccard_pcoa.tsv"
        ),
        bray_pcoa=(
            ORD_DIR
            + "/braycurtis_sqrt_tpm_pcoa.tsv"
        ),
        jaccard_tree=(
            ORD_DIR
            + "/jaccard_upgma.newick"
        ),
        bray_tree=(
            ORD_DIR
            + "/braycurtis_sqrt_tpm_upgma.newick"
        ),
        jaccard_merges=(
            ORD_DIR
            + "/jaccard_upgma_merges.tsv"
        ),
        bray_merges=(
            ORD_DIR
            + "/braycurtis_sqrt_tpm_upgma_merges.tsv"
        ),
        summary=(
            ORD_DIR
            + "/ordination_summary.tsv"
        )
    shell:
        r"""
        set -euo pipefail

        python3 scripts/multisample/multisample_ordination.py \
            --breadth "{input.breadth}" \
            --tpm "{input.tpm}" \
            --outdir "{ORD_DIR}" \
            --threshold {PRESENCE_BREADTH}
        """
