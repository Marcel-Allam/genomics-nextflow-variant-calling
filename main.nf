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

// RAW FASTQ QC (FastQC)
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

// OPTIONAL TRIMMING (fastp)
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

// ALIGN + SORT BAM
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

// -----------------------------------------------------------------------------
// Process: ALIGN_QC
// Purpose:
//   Generate alignment-level QC metrics from the sorted BAM using samtools.
//
// Outputs per sample:
//   - <sample>.flagstat.txt       : High-level mapping summary
//   - <sample>.samtools_stats.txt : Detailed alignment statistics
//   - <sample>.idxstats.txt       : Per-contig read counts
//
// Notes:
//   - This mimics standard QC in clinical and research labs,
//     where alignment quality must be documented.
// -----------------------------------------------------------------------------
process ALIGN_QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/alignment_qc", mode: 'copy', overwrite: true

    input:
    // From ALIGN_AND_SORT: (sample_id, sorted BAM)
    tuple val(sample_id), path(sorted_bam)

    output:
    // Emit all three QC files as a tuple
    tuple val(sample_id),
          path("${sample_id}.flagstat.txt"),
          path("${sample_id}.samtools_stats.txt"),
          path("${sample_id}.idxstats.txt")

    script:
    """
    # 1) High-level mapping summary: mapped %, properly paired, duplicates, etc.
    samtools flagstat ${sorted_bam} > ${sample_id}.flagstat.txt

    # 2) Detailed alignment statistics (read lengths, error rates, etc.).
    samtools stats ${sorted_bam} > ${sample_id}.samtools_stats.txt

    # 3) Per-contig read counts (useful for off-target or contamination checks).
    samtools idxstats ${sorted_bam} > ${sample_id}.idxstats.txt
    """
}

// -----------------------------------------------------------------------------
// Process: COVERAGE_QC
// Purpose:
//   Generate coverage metrics from the sorted BAM using samtools.
//
// Outputs per sample:
//   - <sample>.coverage.txt : Summary coverage metrics (mean depth, % bases ≥ X)
//   - <sample>.depth.txt    : Per-base depth (may be large)
// -----------------------------------------------------------------------------
process COVERAGE_QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/coverage_qc", mode: 'copy', overwrite: true

    input:
    // From ALIGN_AND_SORT: (sample_id, sorted BAM)
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id),
          path("${sample_id}.coverage.txt"),
          path("${sample_id}.depth.txt")

    script:
    """
    # 1) Summary coverage table: mean depth, % bases with adequate coverage
    samtools coverage ${sorted_bam} > ${sample_id}.coverage.txt

    # 2) Raw per-base depth
    # NOTE: This file can be large for WGS. For WES or targeted panels it's fine.
    samtools depth ${sorted_bam} > ${sample_id}.depth.txt
    """
}





// CALL VARIANTS (bcftools)
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

// ANNOTATE VARIANTS WITH VEP
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

// -----------------------------------------------------------------------------
// Process: EXPORT_VARIANTS_CSV
// Purpose:
//   Convert the VEP-annotated VCF into a flat CSV table per sample.
//   This CSV is easier to inspect in Excel / R / Python.
//
// Input:
//   - (sample_id, annotated_vcf_gz) from ANNOTATE_VEP
//
// Output:
//   - <sample>.variants.csv in results/variants_csv/
// -----------------------------------------------------------------------------
process EXPORT_VARIANTS_CSV {

    tag "${sample_id}"

    publishDir "${params.outdir}/variants_csv", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(annotated_vcf_gz)

    output:
    tuple val(sample_id), path("${sample_id}.variants.csv")

    script:
    """
    ${projectDir}/bin/vep_to_csv.py \\
        --vep_vcf ${annotated_vcf_gz} \\
        --out ${sample_id}.variants.csv
    """
}


// -----------------------------------------------------------------------------
// Process: VARIANT_QC
// Purpose:
//   Generate variant-level QC statistics from the raw VCF using bcftools stats.
//   This mimics typical clinical / research lab variant QC and is run on the
//   caller output *before* annotation (best practice).
//
// Inputs:
//   - sample_id : String identifier for the sample
//   - vcf_gz    : Compressed VCF (.vcf.gz) from CALL_VARIANTS
//
// Outputs:
//   - <sample>.vcfstats.txt : Text summary of variant metrics
// -----------------------------------------------------------------------------
process VARIANT_QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/variant_qc", mode: 'copy', overwrite: true

    input:
    // We expect:
    //   - sample_id (a string value)
    //   - vcf_gz    (path to .vcf.gz from CALL_VARIANTS)
    tuple val(sample_id), path(vcf_gz)

    output:
    // We emit:
    //   - the same sample_id
    //   - the bcftools stats file only
    tuple val(sample_id),
          path("${sample_id}.vcfstats.txt")

    script:
    """
    # 1) Generate summary statistics from the VCF.
    #    This includes total SNPs, INDELs, Ts/Tv ratio, per-chromosome counts, etc.
    bcftools stats ${vcf_gz} > ${sample_id}.vcfstats.txt

    # 2) (Optional) Create visual QC plots from the stats file.
    #    This step is disabled by default to avoid LaTeX/TeX dependencies.
    #    To enable plots, install a TeX engine and uncomment the block below.
    #
    # mkdir -p ${sample_id}.vcfstats_plots
    # plot-vcfstats \\
    #   -p ${sample_id}.vcfstats_plots \\
    #   ${sample_id}.vcfstats.txt
    """
}


// -----------------------------------------------------------------------------
// Process: CLINICAL_SUMMARY
// Purpose:
//   Generate a human-readable Markdown report that looks like a simplified
//   clinical variant summary for each sample, including:
//
//   - Alignment QC (flagstat / stats)
//   - Coverage QC (samtools coverage)
//   - Variant QC (bcftools stats on RAW VCF)
//   - VEP annotation summary (high-impact variants)
//
// Inputs:
//   - From ANNOTATE_VEP : (sample_id, annotated_vcf_gz)
//   - From VARIANT_QC   : (sample_id, vcfstats_txt)
//   - From ALIGN_QC     : (sample_id, flagstat, samtools_stats, idxstats)
//   - From COVERAGE_QC  : (sample_id, coverage_txt, depth_txt)
//
// Output:
//   - <sample>_clinical_summary.md : Markdown report summarising key metrics.
//
// Notes:
//   - This is for training / portfolio use only, not for real diagnostics.
// -----------------------------------------------------------------------------
process CLINICAL_SUMMARY {

    tag "${sample_id}"

    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

    input:
    // Annotated VCF from ANNOTATE_VEP
    tuple val(sample_id), path(annotated_vcf_gz)

    // Variant QC from VARIANT_QC (bcftools stats on RAW VCF)
    tuple val(sample_id_qc), path(vcfstats_txt)

    // Alignment QC from ALIGN_QC
    tuple val(sample_id_align),
          path(flagstat_txt),
          path(align_stats_txt),
          path(idxstats_txt)

    // Coverage QC from COVERAGE_QC
    tuple val(sample_id_cov),
          path(coverage_txt),
          path(depth_txt)

    output:
    path("${sample_id}_clinical_summary.*")

    when:
    // Safety check: only run when all the sample IDs match
    sample_id == sample_id_qc &&
    sample_id == sample_id_align &&
    sample_id == sample_id_cov

       script:
    """
    # 1) Generate the Markdown clinical summary
    ${projectDir}/bin/make_clinical_summary.py \
        --sample ${sample_id} \
        --vep_vcf ${annotated_vcf_gz} \
        --vcf_stats ${vcfstats_txt} \
        --flagstat ${flagstat_txt} \
        --align_stats ${align_stats_txt} \
        --idxstats ${idxstats_txt} \
        --coverage ${coverage_txt} \
        --out ${sample_id}_clinical_summary.md

    # 2) Convert Markdown → standalone HTML using pandoc
    pandoc ${sample_id}_clinical_summary.md \
        -f markdown -t html -s \
        -o ${sample_id}_clinical_summary.html
    """
}

// -----------------------------------------------------------------------------
// Process: MULTIQC
// Purpose:
//   Aggregate all QC metrics into a single HTML MultiQC report.
//
// Inputs:
//   - One or more clinical summary reports, just to make sure this process
//     runs only after the rest of the pipeline has finished.
//   - MultiQC will scan the global results/ directory, not the work dir.
//
// Output:
//   - multiqc_report.html in results/multiqc/
// -----------------------------------------------------------------------------
process MULTIQC {

    tag "multiqc"

    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    // We don't actually use these files inside the script; they just ensure
    // that MULTIQC runs after CLINICAL_SUMMARY has completed.
    path(summary_reports)

    output:
    path("multiqc_report.html")

    script:
    """
    # Scan the top-level results directory in the project, not the work dir
    multiqc ${projectDir}/${params.outdir} --outdir .
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

    // 3a) Align reads and produce sorted BAMs
    sorted_bam_ch = ALIGN_AND_SORT(reads_for_alignment_ch)

    // 3c) Alignment QC (samtools flagstat / stats / idxstats)
    align_qc_ch = ALIGN_QC(sorted_bam_ch)

    // 3b) Coverage QC (samtools coverage / depth)
    coverage_qc_ch = COVERAGE_QC(sorted_bam_ch)

    // 4) Call variants
    raw_vcf_ch = CALL_VARIANTS(sorted_bam_ch)

    // 5) Variant QC on RAW VCF (best practice)
    variant_qc_ch = VARIANT_QC(raw_vcf_ch)

    // 5) Annotate variants (VEP)
    annotated_vcf_ch = ANNOTATE_VEP(raw_vcf_ch)

    // 6) Export annotated variants to CSV for manual review
    variants_csv_ch = EXPORT_VARIANTS_CSV(annotated_vcf_ch)

    // 7) Clinical-style Markdown summary per sample
     clinical_summary_ch = CLINICAL_SUMMARY(annotated_vcf_ch, variant_qc_ch, align_qc_ch, coverage_qc_ch)
     
    // 8) MultiQC summary across all QC outputs
    multiqc_ch = MULTIQC(clinical_summary_ch)


    // Final outputs are written to:
    //   - results/fastqc      (FastQC HTML + zip)
    //   - results/trimmed     (trimmed FASTQs, if trimming is enabled)
    //   - results/bam         (sorted BAM + BAM index)
    //   - results/alignment_qc   (flagstat, stats, idxstats)
    //   - results/coverage_qc    (coverage metrics)
    //   - results/duplicate_qc   (duplicate rate metrics)
    //   - results/vcf         (raw + annotated VCFs)
    //   - results/variant_qc  (variant QC stats)
    //   - results/reports     (clinical-style Markdown summaries)
}
