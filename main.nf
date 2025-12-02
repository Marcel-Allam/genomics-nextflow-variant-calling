// main.nf
//
// Genomics – Variant Calling Pipeline (Nextflow DSL2)
//
// Overview:
//   For each sample (paired-end FASTQ), this pipeline:
//     1) Performs raw read QC with FastQC
//     2) Optionally trims reads using fastp (quality-based trimming)
//     3) Aligns reads to a reference genome with bwa-mem2
//     4) Sorts and indexes BAM files using samtools
//     5) Calls variants using bcftools mpileup + call
//     6) Annotates variants using VEP (offline cache, GRCh38)
//
// Notes:
//   - This is a teaching / portfolio pipeline, not a clinical pipeline.
//   - It demonstrates a realistic structure similar to pipelines used in
//     cancer / WES workflows.
//
// Requirements:
//   - Input FASTQ(s) under: data/fastq/*_{R1,R2}.fastq.gz
//   - Reference FASTA at:   data/reference/genome.fa
//   - Conda env in:         envs/genomics_env.yml
//   - VEP offline cache installed in: ~/.vep (e.g. homo_sapiens/GRCh38)

nextflow.enable.dsl = 2

// -----------------------------
// PARAMETERS (with fallbacks)
// -----------------------------

params.reads     = params.reads     ?: "data/fastq/*_{R1,R2}.fastq.gz"
params.reference = params.reference ?: "data/reference/genome.fa"
params.outdir    = params.outdir    ?: "results"

// Ensure main output directory exists
new File(params.outdir).mkdirs()

// -----------------------------
// CHANNELS
// -----------------------------
//
// 1) Paired-end FASTQ files
//
// This creates tuples of:
//   ( sample_id, [R1.fastq.gz, R2.fastq.gz] )
// based on the pattern defined in params.reads.

Channel
    .fromFilePairs( params.reads, flat: false )
    .set { read_pairs_ch }

// -----------------------------
// PROCESSES
// -----------------------------

// 1) RAW FASTQ QC (FastQC)
//
// Always run this first on the untrimmed reads to understand
// library quality, adapter content, and base quality profiles.

process FASTQC {

    tag "${sample_id}"

    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: true

    input:
    // tuple of (sample_id, [R1, R2])
    tuple val(sample_id), path(reads)

    output:
    // All FastQC output files for this sample
    path "*_fastqc.*"

    """
    fastqc \
        --outdir . \
        ${reads.join(' ')}
    """
}

// 1b) OPTIONAL TRIMMING (fastp)
//
// This step is only used when params.trim_reads == true.
//
// It uses fastp to:
//   - detect adapters for paired-end reads
//   - trim low-quality bases from 5' and 3' ends
//   - enforce a Phred quality threshold controlled by params.trim_qual
//
// params.trim_qual:
//   - 20 = standard (Q20 ~ 1% error rate)
//   - 30 = strict (Q30 ~ 0.1% error rate)
//
// You can enable trimming and choose the threshold using:
//   --trim_reads true --trim_qual 20
//   --trim_reads true --trim_qual 30

process TRIM_READS {

    tag "${sample_id}"

    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true

    input:
    // (sample_id, [R1, R2])
    tuple val(sample_id), path(reads)

    output:
    // (sample_id, [R1.trimmed, R2.trimmed])
    tuple val(sample_id), path("${sample_id}_R{1,2}.trimmed.fastq.gz")

    """
    # fastp trimming based on Phred quality score.
    # params.trim_qual controls the threshold:
    #   - 20 = standard (1% error rate)
    #   - 30 = strict (0.1% error rate)
    # This can be set on the Nextflow command line, e.g.:
    #   nextflow run main.nf --trim_reads true --trim_qual 30

    fastp \
      --in1 ${reads[0]} \
      --in2 ${reads[1]} \
      --out1 ${sample_id}_R1.trimmed.fastq.gz \
      --out2 ${sample_id}_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      --cut_front \
      --cut_front_window_size 4 \
      --cut_front_mean_quality ${params.trim_qual} \
      --cut_tail \
      --cut_tail_window_size 4 \
      --cut_tail_mean_quality ${params.trim_qual} \
      --qualified_quality_phred ${params.trim_qual} \
      --thread ${params.cpus}
    """
}

// 2) ALIGN + SORT BAM
//
// Aligns reads to the reference genome using bwa-mem2,
// converts to BAM, sorts, and indexes the BAM.
//
// IMPORTANT:
// We assume the reference FASTA (params.reference) and all
// its index files already exist in data/reference/.
// We do NOT attempt to build a new index inside Nextflow,
// because indexing a whole GRCh38 genome requires >50GB RAM.

process ALIGN_AND_SORT {

    tag "${sample_id}"

    publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true

    input:
    // Paired reads (raw or trimmed, depending on workflow)
    tuple val(sample_id), path(reads)

    output:
    // (sample_id, sorted BAM)
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    """
    # Align reads with bwa-mem2 using the pre-indexed reference.
    # params.reference is an absolute path to data/reference/genome.fa

    bwa-mem2 mem -t ${params.cpus} ${params.reference} ${reads[0]} ${reads[1]} \
        | samtools view -bS - \
        | samtools sort -o ${sample_id}.sorted.bam -

    # Index the BAM for downstream tools.
    samtools index ${sample_id}.sorted.bam
    """
}

// 3) CALL VARIANTS (bcftools)
//
// Uses bcftools mpileup + call to produce a raw VCF of variants.
// The reference genome is taken directly from params.reference,
// so we do not need to stage it as a separate channel.

process CALL_VARIANTS {

    tag "${sample_id}"

    publishDir "${params.outdir}/vcf", mode: 'copy', overwrite: true

    input:
    // (sample_id, sorted BAM)
    tuple val(sample_id), path(sorted_bam)

    output:
    // Emit (sample_id, raw VCF)
    tuple val(sample_id), path("${sample_id}.raw.vcf.gz")

    """
    # Generate pileup and call variants using bcftools.
    # The reference is supplied via params.reference (absolute path).
    bcftools mpileup -f ${params.reference} ${sorted_bam} \
        | bcftools call -mv -Oz -o ${sample_id}.raw.vcf.gz

    # Index the raw VCF for downstream use.
    bcftools index ${sample_id}.raw.vcf.gz
    """
}

// 4) ANNOTATE VARIANTS WITH VEP
//
// Annotates variants using Ensembl VEP in OFFLINE mode, assuming that
// a cache is installed under ~/.vep (e.g. homo_sapiens / GRCh38).
//
// To avoid VEP stdin issues, we:
//   - decompress the bgzipped VCF to a plain text file: input.vcf
//   - tell VEP explicitly that the input format is VCF
//   - stream annotated VCF back out and re-compress with bgzip

process ANNOTATE_VEP {

    tag "${sample_id}"

    publishDir "${params.outdir}/vcf", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf.gz")

    """
    # Decompress the bgzipped VCF to a regular text file so VEP can read it.
    zcat ${vcf} > input.vcf

    # Run VEP in offline mode using the local cache.
    # We explicitly specify --format vcf so VEP doesn't have to guess.
    vep \
        --format vcf \
        --input_file input.vcf \
        --output_file STDOUT \
        --vcf \
        --cache \
        --offline \
        --species homo_sapiens \
        --assembly GRCh38 \
        --everything \
        --force_overwrite \
        --quiet \
        | bgzip -c > ${sample_id}.annotated.vcf.gz

    # Index the annotated VCF for fast querying.
    bcftools index ${sample_id}.annotated.vcf.gz
    """
}

// -----------------------------
// WORKFLOW
// -----------------------------
//
// This block wires all processes together in DSL2 style.
// Processes are called like functions, taking channels as input
// and returning new channels as output.

workflow {

    // 1) Run FastQC on RAW reads (always)
    //
    // This allows you to inspect initial quality and decide whether
    // to re-run the pipeline with trimming enabled.
    fastqc_raw_ch = FASTQC(read_pairs_ch)

    // 2) Optionally trim reads based on a parameter
    //
    // If params.trim_reads == true:
    //   - TRIM_READS(read_pairs_ch) will generate trimmed reads
    //   - These are used for alignment
    //
    // Else:
    //   - The original read_pairs_ch is used directly.
    reads_for_alignment_ch = params.trim_reads ? TRIM_READS(read_pairs_ch)
                                               : read_pairs_ch

    // 3) Align reads and produce sorted BAMs
    sorted_bam_ch = ALIGN_AND_SORT(reads_for_alignment_ch)

    // 4) Call variants
    raw_vcf_ch = CALL_VARIANTS(sorted_bam_ch)

    // 5) Annotate variants (VEP)
    annotated_vcf_ch = ANNOTATE_VEP(raw_vcf_ch)

    // Final outputs are written to:
    //   - results/fastqc   (FastQC HTML + zip)
    //   - results/trimmed  (trimmed FASTQs, if trimming is enabled)
    //   - results/bam      (sorted BAM + BAM index)
    //   - results/vcf      (raw + annotated VCFs)
}
